"""Build Tables 1–6 and link each experimental component to the source records."""
import argparse
import csv
import json
import os
from collections import defaultdict
from functools import lru_cache
from pathlib import Path
from urllib.parse import quote

import yaml

from main_table_pages import render_page


PROJECT = Path(__file__).resolve().parents[1]
PURPOSE = "Experimental family and purpose"
FACTORS = "Experimental factors (number of conditions)"
MEASUREMENTS = "Replication and measurement details"
HEADINGS = {
    1: ["Experiment ID", PURPOSE, "Experimental factors (number of values or conditions)",
        "Absorbance measurements", "Plate type and configuration", MEASUREMENTS],
    2: ["Experiment ID", PURPOSE, "Experimental factors (number of values or conditions)",
        "Absorbance measurements", "Plate type and configuration", MEASUREMENTS],
    3: ["Experiment ID", PURPOSE, "Experimental factors (number of values or conditions)", MEASUREMENTS],
    4: ["Experiment ID", "Organism", "Medium", "Temperature (°C)", FACTORS, MEASUREMENTS],
    5: ["Experiment ID", "Organism", "Temperature (°C)", FACTORS, MEASUREMENTS],
    6: ["Experiment ID", "Organism", "Medium", "Temperature (°C)", FACTORS],
}


def yaml_document(path):
    return yaml.safe_load(path.read_bytes())


@lru_cache(maxsize=None)
def yaml_tree(path):
    return yaml.compose(path.read_text(encoding="utf-8"))


def yaml_lines(path, *keys):
    """Locate a mapping or list entry in the tracked YAML without hard-coded line numbers."""
    node = yaml_tree(path)
    first = 1
    for key in keys:
        if isinstance(key, int):
            node = node.value[key]
            first = node.start_mark.line + 1
        else:
            field, node = next((k, v) for k, v in node.value if k.value == key)
            first = field.start_mark.line + 1
    end = node.end_mark.line + int(node.end_mark.column > node.start_mark.column)
    return first, max(first, end)


def url(data, path, first=None, last=None):
    if not path.exists():
        raise FileNotFoundError(path)
    result = quote(path.relative_to(data).as_posix())
    if first is not None:
        result += f"#L{first}"
        if last is not None and last != first:
            result += f"-L{last}"
    return result


def yaml_url(data, path, *keys):
    return url(data, path, *yaml_lines(path, *keys))


def csv_rows(path):
    """Number physical CSV lines, including the header, for the source links."""
    with path.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            yield reader.line_num, row


def line_urls(data, path, numbers):
    """Link exactly the consecutive runs of matching records, across interleaved cohorts."""
    numbers = sorted(set(numbers))
    if not numbers:
        raise ValueError(f"No matching records in {path}")
    ranges = [[numbers[0], numbers[0]]]
    for n in numbers[1:]:
        if n == ranges[-1][1] + 1:
            ranges[-1][1] = n
        else:
            ranges.append([n, n])
    return [url(data, path, first, last) for first, last in ranges]


def unique_lines(items):
    return "\n".join(dict.fromkeys(items))


def optical_memberships(filename, acq):
    experiment = acq["experiment_id"]
    if filename == "lid_comparison.yaml":
        return ["pH-02" if acq["family"] == "acid" else "FeIII-03"]
    ph = filename == "ratiometric_pH_indicator.yaml"
    code = "pH" if ph else "FeIII"
    step = acq.get("incubation_step_id")
    storage = experiment in {"final_indicator_acid", "final_indicator_ferric"}
    if storage:
        groups = []
        if acq.get("read_index") == 0 or acq["reading_label"] == "spectrum":
            groups.append(f"{code}-01")
        if acq["reading_label"] == "endpoint":
            groups.append("pH-03" if ph else "FeIII-04")
        return groups
    if experiment.startswith("final_indicator_ferrous") or experiment == "pH_UV_20uM_stock_initial_and_spectra":
        return [f"{code}-01"]
    groups = []
    if step == "initial_23_C" and ("_UV" in experiment):
        groups.append("pH-01" if ph else "FeIII-02")
    if step in ("initial_23_C", "incubation_4_C", "incubation_30_C", "incubation_75_C"):
        groups.append("pH-04" if ph else "FeIII-05")
    if step in ("overnight_23_C", "membrane_23_C"):
        groups.append("pH-05" if ph else "FeIII-06")
    return groups


def optical_components(data):
    root = data / "01_indicator_calibration"
    references = root / "reference_measurements.csv"
    reference_rows = list(csv_rows(references))
    for path in sorted(root.rglob("*.yaml")):
        doc = yaml_document(path)
        raw = {p.name: p for p in path.parent.rglob("*.xls*")}
        by_group = defaultdict(list)
        for acq in doc["acquisitions"]:
            groups = optical_memberships(path.name, acq)
            if not groups:
                raise ValueError(f"Unassigned acquisition {acq['source_filename']}")
            for group in groups:
                # Separate sequential incubations, endpoints, spectra and workbook sheets.
                detail = acq.get("incubation_step_id") or acq["reading_label"]
                by_group[group, acq["experiment_id"], detail].append(acq)
        for (group, experiment, detail), records in by_group.items():
            refs = sorted({r["reference_set_id"] for r in records})
            reference_links = [yaml_url(data, path, "reference_sets", ref) for ref in refs]
            chemical_series = {"final_indicator_ferric":"ferric",
                               "final_indicator_ferrous_1M":"ferrous_1M",
                               "final_indicator_ferrous_0.001M":"ferrous_0.001M"}
            series = chemical_series.get(experiment)
            if group == "FeIII-03":
                series = "ferric"
            if series:
                reference_links.extend(line_urls(data, references,
                     [n for n,r in reference_rows if r["series_ID"] == series]))
            sheets = unique_lines(r["source_sheet"] for r in records if "source_sheet" in r)
            yield group, {"Component":experiment, "Acquisition or condition":detail,
                "Experimental design URL":yaml_url(data, path, "experiments", experiment),
                "Standard and reference URLs":unique_lines(reference_links),
                "Physical plate IDs":unique_lines(r["physical_plate_id"] for r in records),
                "Workbook sheets":sheets,
                "Measurement URLs":unique_lines(url(data, raw[r["source_filename"]]) for r in records)}


def media_components(data):
    root = data / "02_media_validation"
    path = root / "media_validation.yaml"
    doc = yaml_document(path)
    files = {p.name:p for p in root.rglob("*.xls")}
    reference_csv = root / "reference_measurements.csv"
    references = list(csv_rows(reference_csv))
    group_ids = {"media_F2S":"F2S-01", "media_AFM1_AFM3":"AFM-01",
                 "media_AFM4_AFM5":"9K-01", "media_Starkey":"Starkey-01"}
    for name, group in doc["groups"].items():
        experiment = group["experiment_id"]
        primary = group.get("reference_set_id") or "media_ferric_stock"
        reference = [yaml_url(data, path, "reference_sets", primary)]
        reference.extend(line_urls(data, reference_csv, [n for n,r in references if r["reference_set_ID"] == primary]))
        acqs = [a for a in doc["acquisitions"] if a["experiment_id"] == experiment]
        if not acqs:
            raise ValueError(experiment)
        yield group_ids[experiment], {"Component": name, "Acquisition or condition":group["condition_role"],
            "Experimental design URL":yaml_url(data, path, "experiments", experiment),
            "Well assignments URL":yaml_url(data, path, "groups", name),
            "Standard and reference URLs":unique_lines(reference),
            "Measured plate rows":", ".join(group["plate_rows"]),
            "Technical replicate columns":", ".join(map(str,group["technical_replicate_columns"])),
            "Physical plate IDs":unique_lines(a["physical_plate_id"] for a in acqs),
            "Measurement URLs":unique_lines(url(data,files[a["source_filename"]]) for a in acqs)}


def microplate_components(data):
    root = data / "03_biological_microplates"
    definitions = root / "biological_definitions.yaml"
    index = root / "chronology" / "condition_timecourse_index.csv"
    acquisitions = root / "chronology" / "acquisitions.csv"
    conditions = defaultdict(list)
    by_plate = defaultdict(list)
    for lineno, r in csv_rows(index):
        conditions[r["experiment_id"]].append((lineno, r["physical_plate_id"]))
    for lineno, r in csv_rows(acquisitions):
        by_plate[r["physical_plate_id"]].append((lineno, r["path"]))
    categories = {("AF","AFM1"):"AF-01", ("LF","HH"):"LF-01", ("AF","AFM4"):"AF-02",
                  ("AF","AFM1_BS"):"AF-03", ("AF","SM4_BS"):"AF-04", ("AT","SM4_BS"):"AT-01"}
    for path in sorted(root.rglob("*_experiments.yaml")):
        doc = yaml_document(path)
        for i, exp in enumerate(doc["experiments"]):
            key = exp["experiment_id"]
            rows = conditions[key]
            plates = sorted({p for _,p in rows})
            files = sorted(item for plate in plates for item in by_plate[plate])
            if not rows or not files:
                raise ValueError(key)
            category = categories[doc["species_id"], exp["medium_id"]]
            dirs = sorted({(data.parent / f).parent for _,f in files})
            yield category, {"Component":key, "Acquisition or condition":exp["name"],
                "Experimental design URL":yaml_url(data,path,"experiments",i),
                "Readout definition URL":yaml_url(data,definitions,"readouts",exp["readout_id"]),
                "Preparation definition URL":yaml_url(data,definitions,"preparations",exp["preparation_id"]),
                "Condition time-course records":unique_lines(line_urls(data,index,[n for n,_ in rows])),
                "Acquisition records":unique_lines(line_urls(data,acquisitions,[n for n,_ in files])),
                "Raw measurement directories":unique_lines(url(data,p) for p in dirs),
                "Example raw acquisition URL":url(data,data.parent/files[0][1])}


def vessel_components(data, family):
    flask = family == "flask"
    folder = data / ("04_flask_validation" if flask else "05_test_tube_validation")
    descriptor = folder / ("flask_validation.yaml" if flask else "test_tube_validation.yaml")
    measurements = folder / ("flask_measurements.csv" if flask else "test_tube_measurements.csv")
    doc = yaml_document(descriptor)
    records = defaultdict(list)
    for n, r in csv_rows(measurements):
        records[r["condition_ID"]].append(n)
    if not flask:
        yield "TB-01", {"Component":doc["lf_starter_followup"]["observation_id"],
            "Acquisition or condition":"14-day visual observation of the LF microplate starters",
            "Recorded observation URL":yaml_url(data, descriptor, "lf_starter_followup")}
    flask_groups = {"2025-07-11":"FL-01", "2025-12-17":"FL-01", "2025-07-16":"FL-02",
                    "2025-07-22":"FL-03", "2025-08-20":"FL-04", "2025-08-29":"FL-05"}
    tube_groups = {"HH":"TB-02","SM4_S0":"TB-03","SM4_tetrathionate":"TB-04"}
    formulations = data / "00_materials_and_strains" / "media_formulations.yaml"
    for i, condition in enumerate(doc["conditions"]):
        key = condition["condition_ID"]
        if flask:
            group = flask_groups[condition["start_date"]]
            component = condition["description"]
            medium = "AFM1"
        else:
            medium = condition["medium_ID"]
            group = tube_groups[medium]
            if group in ("TB-03", "TB-04"):
                group += ":" + str(condition["incubation_temperature_C"])
            component = " / ".join(condition[k] for k in ("organism_code","medium_ID","treatment_ID"))
        yield group, {"Component":key, "Acquisition or condition":component,
           "Experimental design URL":yaml_url(data, descriptor, "conditions",i),
           "Measurement definitions URL":yaml_url(data, descriptor, "columns"),
           "Measurement methods URL":yaml_url(data, descriptor, "measurement_methods"),
           "Medium recipe URL":yaml_url(data, formulations, "media",medium),
           "Measurement URLs":unique_lines(line_urls(data, measurements, records[key]))}



def write_guide(descriptions, destination):
    lines = [
        "# Main-text experimental overview tables", "",
        "Each Markdown file contains the manuscript overview and its caption. Select an "
        "experiment ID to open the component records for that group. Each record identifies "
        "the acquisition or condition, the source field, and its values or links.", "",
        "Links are relative to the table file. For text sources, the link label gives the "
        "physical line numbers, counting the header as line 1. A viewer that opens only the "
        "file can use the displayed lines or the component ID to locate the record. "
        "For Table 4, join the condition and acquisition registries on `physical_plate_id`. "
        "For Tables 5 and 6, use `condition_ID` to locate the vessel measurements.", "",
        "| Table | Experimental overview |", "|---|---|",
    ]
    for number, spec in descriptions.items():
        lines.append(f"| [Table {number}](Table_{number}.md) | {spec['title']} |")
    (destination / "README.md").write_text("\n".join(lines)+"\n", encoding="utf-8")


def build(data_dir, destination, *, link_directory=None):
    """Write captioned Markdown tables with paths relative to their final location."""
    data, destination = Path(data_dir).resolve(), Path(destination)
    link_directory = destination if link_directory is None else Path(link_directory)
    data_href = quote(Path(os.path.relpath(data, link_directory)).as_posix())
    descriptions = json.loads(Path(__file__).with_name("main_table_descriptions.json").read_text())
    groups = defaultdict(list)
    for key, entry in optical_components(data):
        groups[key].append(entry)
    for series in (media_components(data), microplate_components(data),
                   vessel_components(data, "flask"), vessel_components(data, "tube")):
        for key, entry in series:
            groups[key].append(entry)
    expected = {k for d in descriptions.values() for k in d["rows"]}
    if set(groups) != expected:
        raise ValueError(f"Experimental groups differ: {set(groups) ^ expected}")
    destination.mkdir(parents=True, exist_ok=True)
    for number, definition in descriptions.items():
        render_page(number, definition, HEADINGS[int(number)], groups, destination, data_href)
    write_guide(descriptions, destination)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=PROJECT / "data")
    parser.add_argument("--output", type=Path, default=Path(__file__).parent/"main_tables")
    args = parser.parse_args()
    build(args.data, args.output)
