"""Optical."""
from __future__ import annotations
from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path
from typing import Any
import xml.etree.ElementTree as ET
from zipfile import ZipFile
import numpy as np
import pandas as pd
from .metadata import load_metadata


DESCRIPTORS = (
    "data/01_indicator_calibration/01_ratiometric_pH_indicator/ratiometric_pH_indicator.yaml",
    "data/01_indicator_calibration/02_iron_III_indicator/iron_III_indicator.yaml",
    "data/01_indicator_calibration/03_lid_comparison/lid_comparison.yaml",
    "data/02_media_validation/media_validation.yaml",
)


IRON_REFERENCES = "data/01_indicator_calibration/reference_measurements.csv"


MEDIA_REFERENCES = "data/02_media_validation/reference_measurements.csv"


_XML_NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


def _json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, default=str, allow_nan=False)


def _safe_path(root: Path, relative: str) -> Path:
    """Resolve a source path within the project root."""
    candidate = Path(relative)
    if candidate.is_absolute() or ".." in candidate.parts:
        raise ValueError(f"Unsafe source path: {relative!r}")
    result = (root / candidate).resolve()
    if not result.is_relative_to(root):
        raise ValueError(f"Source path escapes the project root: {relative!r}")
    if not result.is_file():
        raise FileNotFoundError(f"Missing source: {relative}")
    return result


def _load_descriptor(root: Path, relative: str) -> tuple[Path, dict]:
    path = _safe_path(root, relative)
    document = load_metadata(path)
    if not isinstance(document, dict) or document.get("metadata_format") not in {"indicator_calibration", "media_validation"}:
        raise ValueError(f"Unsupported optical metadata format: {relative}")
    # Descriptors can omit format_version. Validate it when present.
    if "format_version" in document and document["format_version"] != 1:
        raise ValueError(f"Unsupported optical format_version: {relative}")
    required = {"experiments": dict, "reference_sets": dict, "acquisitions": list}
    required.update({"physical_plates": dict, "readouts": dict}
                    if document["metadata_format"] == "indicator_calibration"
                    else {"groups": dict, "reader_settings": dict})
    for key, kind in required.items():
        if not isinstance(document.get(key), kind) or not document[key]:
            raise ValueError(f"Missing or invalid optical {key}: {relative}")
    return path, document


def _workbook_sheets(path: Path) -> dict[str, bool]:
    """Return sheet names and whether each contains stored cell content.

    Empty template sheets are not acquisitions. Nonempty unmapped sheets fail
    registry validation rather than being silently dropped.
    """
    with ZipFile(path) as archive:
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        rels = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        targets = {rel.attrib["Id"]: rel.attrib["Target"] for rel in rels}
        sheets = {}
        for sheet in workbook.findall("m:sheets/m:sheet", _XML_NS):
            rid = sheet.attrib["{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"]
            target = targets[rid]
            member = target.lstrip("/") if target.startswith("/") else f"xl/{target}"
            tree = ET.fromstring(archive.read(member))
            populated = any(
                cell.find("m:v", _XML_NS) is not None
                or cell.find("m:is", _XML_NS) is not None
                or cell.find("m:f", _XML_NS) is not None
                for cell in tree.findall(".//m:sheetData/m:row/m:c", _XML_NS)
            )
            sheets[sheet.attrib["name"]] = populated
        return sheets


def _descriptor_rows(root: Path, relative: str) -> list[dict]:
    path, document = _load_descriptor(root, relative)
    family = document.get("family", "media_validation")
    source_index: dict[str, list[Path]] = defaultdict(list)
    for source in path.parent.rglob("*"):
        if source.is_file() and source.suffix.lower() in {".xls", ".xlsx"}:
            resolved = _safe_path(root, source.relative_to(root).as_posix())
            source_index[source.name].append(resolved)
    seen: set[tuple[str, str]] = set()
    source_sheets: dict[Path, set[str]] = defaultdict(set)
    rows = []
    for acquisition in document["acquisitions"]:
        filename = acquisition["source_filename"]
        if Path(filename).name != filename or "\\" in filename:
            raise ValueError(f"source_filename must be an exact basename: {filename!r}")
        candidates = source_index.get(filename, [])
        if not candidates:
            raise FileNotFoundError(f"{relative}: missing export {filename!r}")
        if len(candidates) != 1:
            raise ValueError(f"{relative}: ambiguous export {filename!r}")
        source = candidates[0]
        sheet = acquisition.get("source_sheet", "")
        source_key = (filename, sheet)
        if source_key in seen:
            raise ValueError(f"Duplicate acquisition for {filename!r}, sheet {sheet!r}")
        seen.add(source_key)
        source_sheets[source].add(sheet)
        if (source.suffix.lower() == ".xlsx") != bool(sheet):
            raise ValueError(f"Workbook acquisition requires source_sheet: {filename!r}")
        if acquisition.get("sha256"):
            actual = hashlib.sha256(source.read_bytes()).hexdigest()
            if actual != acquisition["sha256"]:
                raise ValueError(f"Checksum mismatch: {filename}")
        if "size_bytes" in acquisition and source.stat().st_size != acquisition["size_bytes"]:
            raise ValueError(f"File-size mismatch: {filename}")
        experiment_id = acquisition["experiment_id"]
        if experiment_id not in document["experiments"]:
            raise ValueError(f"Unknown experiment {experiment_id!r}: {filename}")
        experiment = document["experiments"][experiment_id]
        plate_id = acquisition["physical_plate_id"]
        if plate_id not in experiment.get("physical_plate_ids", []):
            raise ValueError(f"Plate {plate_id!r} not assigned to {experiment_id!r}")
        reference_id = acquisition.get("reference_set_id", "")
        readout_id = acquisition.get("readout_id", "")
        if document["metadata_format"] == "indicator_calibration":
            if plate_id not in document["physical_plates"]:
                raise ValueError(f"Unknown physical plate: {plate_id}")
            plate = document["physical_plates"][plate_id]
            if reference_id not in document["reference_sets"]:
                raise ValueError(f"Unknown reference set: {reference_id}")
            if plate.get("reference_set_id") != reference_id:
                raise ValueError(f"Plate/reference mismatch: {plate_id}, {reference_id}")
            if readout_id not in document["readouts"]:
                raise ValueError(f"Unknown optical readout: {readout_id}")
            readout = document["readouts"][readout_id]
            references = {reference_id: document["reference_sets"][reference_id]}
            groups = {}
        else:
            plate = {"stock_set_id": experiment["stock_set_id"]}
            readout = document["reader_settings"]
            groups = {gid: document["groups"][gid] for gid in experiment["group_ids"]}
            references = {}
            for group_id, group in groups.items():
                if group["experiment_id"] != experiment_id:
                    raise ValueError(f"Group/experiment mismatch: {group_id}")
                if group.get("reference_set_id"):
                    ref_id = group["reference_set_id"]
                    references[ref_id] = document["reference_sets"][ref_id]
        row = dict(acquisition)
        # Preserve the descriptor family separately from lid acquisition acid/iron.
        row.update(
            family=family,
            acquisition_family=acquisition.get("family", family),
            source_path=source.relative_to(root).as_posix(),
            source_sheet=sheet,
            metadata_path=relative,
            acquisition_id=acquisition.get("acquisition_id", f"{family}:{filename}:{sheet}"),
            reference_set_id=reference_id,
            readout_id=readout_id,
            physical_plate_metadata_json=_json(plate),
            experiment_metadata_json=_json(experiment),
            reference_metadata_json=_json(references),
            readout_metadata_json=_json(readout),
            group_metadata_json=_json(groups),
            preparation_metadata_json=_json(document.get("preparation", document.get("study_design", {}))),
            acquisition_metadata_json=_json(acquisition),
        )
        # Nested source fields have a CSV-stable representation; the full entry
        # is also preserved above so no schema-specific metadata is discarded.
        for key, value in list(row.items()):
            if isinstance(value, (dict, list, tuple)):
                row[key] = _json(value)
        rows.append(row)
    unmapped = set(source_index) - {name for name, _ in seen}
    if unmapped:
        raise ValueError(f"Unmapped optical exports in {relative}: {sorted(unmapped)}")
    for source, mapped_sheets in source_sheets.items():
        if source.suffix.lower() == ".xlsx":
            available = _workbook_sheets(source)
            missing = mapped_sheets - set(available)
            unmapped_sheets = {name for name, populated in available.items() if populated} - mapped_sheets
            if missing or unmapped_sheets:
                raise ValueError(f"Workbook sheet coverage: missing={sorted(missing)}, unmapped={sorted(unmapped_sheets)}")
    return rows


def load_optical_registry(root: str | Path) -> pd.DataFrame:
    """Load all four required descriptors, failing on missing or unmapped sources.

    ``root`` is the parent directory of ``data/``. Shared physical plates have
    the same ID across descriptors. Each acquisition has its own record.
    """
    root = Path(root).resolve()
    rows = [row for relative in DESCRIPTORS for row in _descriptor_rows(root, relative)]
    frame = pd.DataFrame(rows)
    if frame["acquisition_id"].duplicated().any():
        raise ValueError("Optical acquisition IDs are not globally unique")
    if frame.duplicated(["source_path", "source_sheet"]).any():
        raise ValueError("A source acquisition is registered in multiple descriptors")
    return frame


def final_well_concentration(
    value: float,
    volume: float,
    *,
    reference_volume_uL: float = 200,
    concentration_basis: str = "before_final_addition",
    quantity: str = "concentration",
) -> float:
    """Convert a standard concentration using the actual well preparation volumes.

    ``volume`` is the final well volume in uL.  Thus undyed 200-uL standards are
    unchanged, dyed 210-uL wells receive 200/210, and a 100-uL media stock mixed to
    200 uL uses ``reference_volume_uL=100``. Already-final concentrations are
    returned unchanged with ``concentration_basis='final_well'``. No analytical
    dilution or salt-to-metal stoichiometric correction occurs here. pH is not a
    concentration and is deliberately rejected.
    """
    if quantity.lower() == "ph":
        raise ValueError("pH must not receive a linear dilution correction")
    if concentration_basis not in {"before_final_addition", "200_uL_standard", "final_well"}:
        raise ValueError(f"Unknown concentration basis: {concentration_basis}")
    if not math.isfinite(volume) or not math.isfinite(reference_volume_uL) or volume <= 0 or reference_volume_uL <= 0:
        raise ValueError("Well and reference volumes must be finite and positive")
    if concentration_basis == "final_well" or volume == reference_volume_uL:
        return value
    if volume < reference_volume_uL:
        raise ValueError("Final volume is less than the dispensed reference volume")
    return value * reference_volume_uL / volume


def load_iron_references(root: str | Path) -> pd.DataFrame:
    """Read the long iron table without reapplying analytical or well dilution.

    IDs and literal NA in the inapplicable reference_replicate field are retained;
    reference replicates must not be paired to optical plates by their numbers.
    """
    root = Path(root).resolve()
    _, document = _load_descriptor(root, DESCRIPTORS[1])
    reference_metadata = document["reference_measurements"]
    table = reference_metadata["table"]
    frame = pd.read_csv(_safe_path(root, IRON_REFERENCES), encoding="utf-8-sig", keep_default_na=False)
    required = set(reference_metadata["csv_format"]["column_order"])
    if not required.issubset(frame.columns):
        raise ValueError(f"Missing iron reference columns: {sorted(required - set(frame.columns))}")
    if frame.empty or frame["measurement_ID"].eq("").any() or frame["measurement_ID"].duplicated().any():
        raise ValueError("Iron measurement_ID must be nonempty and unique")
    for column in ("measured_concentration_M", "reference_volume_uL", "optical_final_volume_uL", "plate_column", "nominal_H2SO4_concentration_M"):
        frame[column] = pd.to_numeric(frame[column], errors="raise")
        if not np.isfinite(frame[column]).all():
            raise ValueError(f"Nonfinite iron reference values: {column}")
    if not frame["concentration_basis"].eq(table["concentration_basis"]).all():
        raise ValueError("Unexpected iron reference concentration basis")
    if not frame["analyte"].isin(["total_Fe", "Fe_II"]).all():
        raise ValueError("Unknown iron reference analyte")
    if (frame["measured_concentration_M"] < 0).any():
        raise ValueError("Negative iron reference concentration")
    if len(frame) != table["expected_rows"]:
        raise ValueError("Iron reference row coverage differs from metadata")
    if frame["reference_sample_ID"].nunique() != table["reference_samples"]:
        raise ValueError("Iron reference sample coverage differs from metadata")
    for analyte in ("total_Fe", "Fe_II"):
        if frame["analyte"].eq(analyte).sum() != table[f"{analyte}_rows"]:
            raise ValueError(f"Iron reference analyte coverage differs: {analyte}")
    if not frame["reference_volume_uL"].eq(table["reference_solution_volume_uL"]).all():
        raise ValueError("Iron reference volume differs from metadata")
    dyed_volume = table["optical_final_well_volume_uL_by_condition"]["thymol_blue"]
    if not frame["optical_final_volume_uL"].eq(dyed_volume).all():
        raise ValueError("Iron dyed optical volume differs from metadata")
    series_metadata = {series["series_ID"]: series for series in reference_metadata["series"]}
    if set(frame["series_ID"]) != set(series_metadata):
        raise ValueError("Iron reference series coverage differs from metadata")
    # The final calibration reference sets identify the exact sample/column
    # associations. Configuration arrays are independent references and must not
    # be joined to this CSV merely because a numerical concentration coincides.
    reference_sets = {}
    for reference in document["reference_sets"].values():
        axis = reference.get("axis_details", {})
        if axis.get("reference_filename") == reference_metadata["data_filename"]:
            series_id = axis["series_ID"]
            if series_id in reference_sets:
                raise ValueError(f"Ambiguous iron reference set for series: {series_id}")
            reference_sets[series_id] = reference
    if set(reference_sets) != set(series_metadata):
        raise ValueError("Iron reference series are not fully mapped by YAML reference_sets")
    for series_id, metadata in series_metadata.items():
        series = frame[frame["series_ID"].eq(series_id)]
        if series["reference_sample_ID"].nunique() != metadata["reference_samples"]:
            raise ValueError(f"Iron sample coverage differs for series: {series_id}")
        for analyte in ("total_Fe", "Fe_II"):
            if series["analyte"].eq(analyte).sum() != metadata[f"{analyte}_rows"]:
                raise ValueError(f"Iron analyte coverage differs for series: {series_id}")
        for field in ("preparation_salt", "nominal_H2SO4_concentration_M"):
            if not series[field].eq(metadata[field]).all():
                raise ValueError(f"Iron series {field} differs from metadata: {series_id}")
        reference = reference_sets[series_id]
        axis = reference["axis_details"]
        columns = reference["row_groups"][0]["plate_columns"]
        total_values = axis["measured_total_Fe_before_final_addition_M_by_column"]
        if len(columns) != len(total_values) or len(columns) != metadata["reference_samples"]:
            raise ValueError(f"Invalid iron reference column mapping: {series_id}")
        if set(series["plate_column"]) != set(columns):
            raise ValueError(f"Iron reference plate-column coverage differs: {series_id}")
        replicate_count, remainder = divmod(metadata["Fe_II_rows"], metadata["reference_samples"])
        if remainder or replicate_count < 1 or metadata["total_Fe_rows"] != len(columns):
            raise ValueError(f"Invalid analyte replication metadata: {series_id}")
        expected_keys = {("total_Fe", "NA")} | {("Fe_II", str(rep)) for rep in range(1, replicate_count + 1)}
        pattern = axis["reference_sample_ID_pattern"]
        if not pattern.endswith("cNN"):
            raise ValueError(f"Unsupported iron sample-ID pattern: {pattern}")
        for column, concentration in zip(columns, total_values):
            sample_id = pattern[:-2] + f"{column:02d}"
            sample = series[series["plate_column"].eq(column)]
            if not sample["reference_sample_ID"].eq(sample_id).all():
                raise ValueError(f"Iron reference sample/column mapping differs: {sample_id}")
            keys = list(zip(sample["analyte"], sample["reference_replicate"].astype(str)))
            if len(keys) != len(expected_keys) or set(keys) != expected_keys:
                raise ValueError(f"Iron sample analyte/replicate coverage differs: {sample_id}")
            total = sample.loc[sample["analyte"].eq("total_Fe"), "measured_concentration_M"].iloc[0]
            if not math.isclose(total, concentration, rel_tol=1e-12, abs_tol=0):
                raise ValueError(f"Iron total-Fe reference disagrees with YAML axis: {sample_id}")
    return frame


def load_media_references(root: str | Path) -> pd.DataFrame:
    """Read all independently identified media references and validate YAML joins."""
    root = Path(root).resolve()
    _, document = _load_descriptor(root, DESCRIPTORS[-1])
    frame = pd.read_csv(_safe_path(root, MEDIA_REFERENCES), encoding="utf-8-sig", keep_default_na=False)
    required = {"reference_ID", "reference_set_ID", "level_ID", "plate_row", "quantity", "value", "unit", "measurement_basis", "method_ID"}
    if not required.issubset(frame.columns):
        raise ValueError(f"Missing media reference columns: {sorted(required - set(frame.columns))}")
    if frame.empty or frame["reference_ID"].eq("").any() or frame["reference_ID"].duplicated().any():
        raise ValueError("Media reference_ID must be nonempty and unique")
    frame["value"] = pd.to_numeric(frame["value"], errors="raise")
    if not np.isfinite(frame["value"]).all():
        raise ValueError("Media reference values must be finite")
    expected_count = document.get("reference_table", {}).get("record_count")
    if expected_count is not None and len(frame) != expected_count:
        raise ValueError("Media reference record count differs from metadata")
    if not frame["method_ID"].isin(document["reference_methods"]).all():
        raise ValueError("Unknown media reference method_ID")
    lookup = frame.set_index("reference_ID", verify_integrity=True)
    referenced = set()
    for reference_set_id, reference in document["reference_sets"].items():
        identifiers = reference.get("reference_ids_by_row", reference.get("reference_ids", []))
        rows = reference.get("plate_rows")
        if rows is not None and len(rows) != len(identifiers):
            raise ValueError(f"Reference row/ID length mismatch: {reference_set_id}")
        for position, identifier in enumerate(identifiers):
            if identifier not in lookup.index:
                raise ValueError(f"Missing media reference: {identifier}")
            record = lookup.loc[identifier]
            if record["reference_set_ID"] != reference_set_id:
                raise ValueError(f"Reference set mismatch: {identifier}")
            if rows is not None and record["plate_row"] != rows[position]:
                raise ValueError(f"Reference row mismatch: {identifier}")
            referenced.add(identifier)
    if referenced != set(lookup.index):
        raise ValueError("Media references contain records not mapped by reference_sets")
    return frame
