"""Biological."""
from __future__ import annotations
from functools import lru_cache
import json
from pathlib import Path
from typing import Iterable
import numpy as np
import pandas as pd
from .metadata import load_metadata
from .biological_chronology import load_chronology
from .plate_reader import PlateData, read_plate


BIOLOGICAL_DIR = Path("data/03_biological_microplates")


WAVELENGTHS = (380, 750, 545, 485, 600, 428, 430)


MOLAR_TO_MM = {"uM": 0.001, "mM": 1.0, "M": 1000.0}


NUMERIC_OVER_RANGE = 4.0


def _json(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), default=str)


def _selection(value: str | Iterable[str] | None) -> set[str] | None:
    return None if value is None else {value} if isinstance(value, str) else set(value)


def _paths(root: str | Path) -> tuple[Path, Path]:
    root = Path(root).resolve()
    # Accept repository root or its biological data directory.
    if (root / "biological_definitions.yaml").is_file():
        return root.parent.parent, root
    directory = root / BIOLOGICAL_DIR
    if not (directory / "biological_definitions.yaml").is_file():
        raise FileNotFoundError(f"Biological definitions not found under {root}")
    return root, directory


@lru_cache(maxsize=8192)
def _read_cached(path: str, size: int, mtime_ns: int) -> PlateData:
    """Cache by physical path and file state, never by experiment membership."""
    return read_plate(Path(path))


def _read(path: Path) -> PlateData:
    state = path.stat()
    return _read_cached(str(path.resolve()), state.st_size, state.st_mtime_ns)


def _documents(root: str | Path):
    repository, directory = _paths(root)
    definitions = load_metadata(directory / "biological_definitions.yaml")
    replication = definitions.get("culture_preparation", {})
    layout = definitions.get("plate_layout", {})
    required_replication = {
        "replicate_microplates_use_independent_culture_sets": True,
        "biological_replicate_unit": "independent_microplate",
        "biotic_rows_counted_as_separate_biological_replicates": False,
        "separately_grown_starter_cultures_per_microplate": len(layout.get("default_biotic_rows", [])),
    }
    if any(replication.get(field) != value for field, value in required_replication.items()):
        raise ValueError("Unsupported biological replication metadata: independent microplates required")
    if layout.get("biotic_culture_number_scope") != "within_microplate":
        raise ValueError("Unsupported culture numbering scope: within_microplate required")
    documents = []
    for path in sorted(directory.glob("*/*_experiments.yaml")):
        document = load_metadata(path)
        species = document["species_id"]
        if document.get("definitions_id") != definitions["document_id"]:
            raise ValueError(f"{path}: unresolved definitions_id")
        if species not in definitions["species"]:
            raise ValueError(f"{path}: unknown species {species}")
        experiment_ids, group_ids = set(), set()
        for experiment in document["experiments"]:
            identity = experiment["experiment_id"]
            if identity in experiment_ids:
                raise ValueError(f"{path}: duplicate experiment_id {identity}")
            experiment_ids.add(identity)
            for field, registry in (("medium_id", "media"), ("preparation_id", "preparations"),
                                    ("readout_id", "readouts")):
                if experiment[field] not in definitions[registry]:
                    raise ValueError(f"{identity}: unresolved {field}={experiment[field]}")
            preparation = definitions["preparations"][experiment["preparation_id"]]
            if preparation["readout_id"] != experiment["readout_id"]:
                raise ValueError(f"{identity}: preparation/readout mismatch")
            for group in experiment["plate_groups"]:
                if group["group_id"] in group_ids:
                    raise ValueError(f"{path}: duplicate group_id {group['group_id']}")
                group_ids.add(group["group_id"])
                if group["concentration_unit"] not in MOLAR_TO_MM:
                    raise ValueError(f"{identity}: unsupported concentration unit {group['concentration_unit']}")
                for plate_id in group["plate_ids"]:
                    if plate_id not in document["plates"]:
                        raise ValueError(f"{identity}: unresolved physical plate {plate_id}")
                    blanks = document["plates"][plate_id].get("blank_condition_columns", [])
                    if set(group["column_conditions"]).intersection(blanks):
                        raise ValueError(f"{identity}: condition assigned to a declared blank column")
                for column, condition in group["column_conditions"].items():
                    if column not in definitions["plate_layout"]["default_condition_columns"]:
                        raise ValueError(f"{identity}: invalid condition column {column}")
                    unknown = set(condition) - set(experiment["compounds"])
                    if unknown:
                        raise ValueError(f"{identity}: unresolved compounds {sorted(unknown)}")
        documents.append((path, document))
    if not documents:
        raise ValueError(f"No biological experiment documents under {directory}")
    supplied_species = [document["species_id"] for _, document in documents]
    missing_species = set(definitions["species"]) - set(supplied_species)
    if missing_species:
        raise ValueError(f"Missing experiment documents for declared species: {sorted(missing_species)}")
    duplicated_species = {species for species in supplied_species if supplied_species.count(species) > 1}
    if duplicated_species:
        raise ValueError(f"Duplicate experiment documents for species: {sorted(duplicated_species)}")
    return repository, definitions, documents


def initial_ph_metadata(experiment: dict) -> dict:
    """Transcribe initial pH and preserve its stated experimental scope.

    ``initial_pH_measured`` contains the recorded value.
    its explicit status/scope determine whether it is a measurement or summary.
    A plus/minus value is an SD only when the deposit explicitly identifies it
    as such. A measurement count describes the experiment-level summary, not
    a measurement assigned to every condition or acquisition.
    """
    reported = experiment.get("initial_pH") or {}
    status = experiment.get("initial_pH_status", reported.get("status"))
    if status is None:
        status = ("measured_experiment_summary" if reported.get("measurement_count") is not None
                  else "reported_value_without_explicit_status" if reported else "not_recorded")
    return {
        "initial_pH_measured": reported.get("value", np.nan),
        "initial_pH_reported_plus_minus": reported.get("reported_plus_minus", np.nan),
        "initial_pH_status": status,
        "initial_pH_scope": reported.get("scope", experiment.get("initial_pH_scope", "not_recorded")),
        "initial_pH_uncertainty_type": reported.get("uncertainty_type", "not_specified"),
        "initial_pH_measurement_count": reported.get("measurement_count", np.nan),
    }


def _condition_metadata(experiment: dict, group: dict, condition: dict,
                        definitions: dict) -> dict:
    preparation = definitions["preparations"][experiment["preparation_id"]]
    final_volume = preparation["final_well_volume_uL"]
    values, units, basis, prep_mm, final_mm, components = {}, {}, {}, {}, {}, {}
    for name, concentration in condition.items():
        components[name] = {**experiment["compounds"][name], "concentration": concentration,
                            "unit": group["concentration_unit"],
                            "volume_basis_uL": group["volume_basis_uL"], "role": "variable"}
    for fixed in experiment.get("fixed_components", []):
        name = fixed["compound_id"]
        if name in components:
            raise ValueError(f"{experiment['experiment_id']}: duplicate fixed/variable component {name}")
        components[name] = {**experiment["compounds"].get(name, {}), **fixed, "role": "fixed"}
    default_basis = definitions["concentration_conventions"]["additive_quantity"]
    for name, component in components.items():
        unit = component["unit"]
        if unit not in MOLAR_TO_MM:
            raise ValueError(f"Unsupported component concentration unit {unit}")
        concentration = float(component["concentration"])
        if not np.isfinite(concentration) or concentration < 0:
            raise ValueError(f"Invalid concentration for {name}: {concentration}")
        if float(component["volume_basis_uL"]) <= 0 or final_volume <= 0:
            raise ValueError("Concentration conversion requires positive volumes")
        values[name] = component["concentration"]
        units[name] = unit
        basis[name] = component.get("concentration_basis", default_basis)
        prep_mm[name] = concentration * MOLAR_TO_MM[unit]
        final_mm[name] = prep_mm[name] * component["volume_basis_uL"] / final_volume
    od = experiment["initial_OD600_200_uL"]
    return {
        "condition_json": _json(values), "component_units_json": _json(units),
        "component_basis_json": _json(basis), "components_json": _json(components),
        "component_preparation_mM_json": _json(prep_mm), "component_final_mM_json": _json(final_mm),
        "is_zero_variable_addition": all(value == 0 for value in condition.values()),
        "initial_OD600_200_uL": od,
        "initial_OD600_final": od * preparation["preparation_volume_uL"] / final_volume,
        "final_well_volume_uL": final_volume, "dye_volume_uL": preparation["dye_volume_uL"],
        **initial_ph_metadata(experiment),
        "nominal_initial_pH_json": _json(experiment.get("nominal_initial_pH")),
        "additive_stock_pH": experiment.get("additive_stock_pH", np.nan),
    }


def _well_values(reading: PlateData, well: str) -> dict:
    result = {}
    for wavelength in WAVELENGTHS:
        name = f"A{wavelength}"
        try:
            value = reading.value(well, wavelength)
            parse_error = reading.is_parse_error(well, wavelength)
            censored = reading.is_censored_high(well, wavelength)
        except KeyError:
            value, parse_error, censored = np.nan, False, False
        result[name] = value  # Preserve source values even when the channel is ineligible.
        result[f"{name}_parse_error"] = parse_error
        result[f"{name}_censored_high"] = censored
        result[f"{name}_overrange"] = bool(censored or (np.isfinite(value) and value > NUMERIC_OVER_RANGE))
        result[f"{name}_eligible"] = bool(np.isfinite(value) and not parse_error
                                           and not result[f"{name}_overrange"])
    result["A380_corrected"] = result["A380"] - result["A750"]
    result["A545_corrected"] = result["A545"] - result["A750"]
    result["A485_corrected"] = result["A485"] - result["A750"]
    denominator = result["A485_corrected"]
    result["R545_485_corrected"] = (result["A545_corrected"] / denominator
                                    if np.isfinite(denominator) and denominator != 0 else np.nan)
    result["R545_485"] = (result["A545"] / result["A485"]
                          if np.isfinite(result["A485"]) and result["A485"] != 0 else np.nan)
    result["eligible_A380_corrected"] = result["A380_eligible"] and result["A750_eligible"]
    result["eligible_R545_485_corrected"] = bool(
        all(result[f"A{wave}_eligible"] for wave in (545, 485, 750)) and denominator > 0)
    return result


def load_biological_observations(root: str | Path, species: str | Iterable[str] | None = None,
                                 experiment_ids: str | Iterable[str] | None = None) -> pd.DataFrame:
    """Return condition-well observations, preserving physical and culture IDs.

    ``condition_json`` includes fixed and variable *additions*, not the basal
    medium ingredients. Chemical concentration basis is retained separately; no
    formula-unit-to-elemental conversion is silently applied. OD columns record
    the experimental biotic inoculum target; abiotic rows have no ``culture_id``.
    Culture IDs identify separately grown inocula within each microplate. The
    ``biological_replicate_id`` identifies the independent microplate; its rows
    are technical replicates to be averaged before independent inference.
    ``eligible`` refers to the scatter-corrected primary readout. All raw values
    remain present, with wavelength-specific flags for sensitivity analyses.
    """
    repository, definitions, documents = _documents(root)
    selected_species, selected_experiments = _selection(species), _selection(experiment_ids)
    known_species = {document["species_id"] for _, document in documents}
    known_experiments = {experiment["experiment_id"] for _, document in documents
                         for experiment in document["experiments"]}
    if selected_species is not None and selected_species - known_species:
        raise ValueError(f"Unknown species: {sorted(selected_species - known_species)}")
    if selected_experiments is not None and selected_experiments - known_experiments:
        raise ValueError(f"Unknown experiments: {sorted(selected_experiments - known_experiments)}")
    documents = [(path, document) for path, document in documents
                 if selected_species is None or document["species_id"] in selected_species]
    # Registry timing always includes every acquisition of a selected physical
    # plate, including alternate labels, before any experiment filtering.
    registry, condition_index = load_chronology(repository, definitions, documents, _read)
    timecourse_ids = {
        (row.species_id, row.experiment_id, row.group_id, row.plate_id, row.column): row.timecourse_id
        for row in condition_index.itertuples(index=False)
    }
    plate_records = {key: group.to_dict("records")
                     for key, group in registry.groupby(["species", "plate_id"], sort=False)}
    rows, value_cache = [], {}
    layout = definitions["plate_layout"]
    row_kinds = {row: kind for kind in ("abiotic", "biotic")
                 for row in layout[f"default_{kind}_rows"]}
    for _, document in documents:
        organism = document["species_id"]
        for experiment in document["experiments"]:
            if selected_experiments is not None and experiment["experiment_id"] not in selected_experiments:
                continue
            for group in experiment["plate_groups"]:
                for column, condition in group["column_conditions"].items():
                    condition_fields = _condition_metadata(experiment, group, condition, definitions)
                    common = {"experiment_id": experiment["experiment_id"], "group_id": group["group_id"],
                              "medium_id": experiment["medium_id"], "preparation_id": experiment["preparation_id"],
                              "readout": experiment["readout_id"], "readout_id": experiment["readout_id"],
                              "column": int(column), **condition_fields}
                    for plate_id in group["plate_ids"]:
                        timecourse_id = timecourse_ids[(organism, experiment["experiment_id"],
                                                      group["group_id"], plate_id, int(column))]
                        for acquisition in plate_records[(organism, plate_id)]:
                            reading = _read(repository / acquisition["source"])
                            for row, kind in row_kinds.items():
                                well = f"{row}{column}"
                                physical_key = f"{acquisition['acquisition_id']}#{well}"
                                if physical_key not in value_cache:
                                    value_cache[physical_key] = _well_values(reading, well)
                                values = value_cache[physical_key]
                                primary_waves = (380, 750) if common["readout"] == "iron" else (545, 485, 750)
                                eligibility = "eligible_A380_corrected" if common["readout"] == "iron" else "eligible_R545_485_corrected"
                                rows.append({**acquisition, **common, **values, "timecourse_id": timecourse_id,
                                             "well": well, "row": row,
                                             "kind": kind, "physical_observation_key": physical_key,
                                             "culture_id": f"{acquisition['physical_plate_id']}_{row}" if kind == "biotic" else None,
                                             "culture_number": layout["biotic_culture_numbers_by_row"].get(row),
                                             "parse_error": any(values[f"A{wave}_parse_error"] for wave in primary_waves),
                                             "censored_high": any(values[f"A{wave}_censored_high"] for wave in primary_waves),
                                             "overrange": any(values[f"A{wave}_overrange"] for wave in primary_waves),
                                             "eligible": values[eligibility]})
    frame = pd.DataFrame(rows)
    if frame.empty:
        raise ValueError("No biological observations match the supplied species/experiment selection")
    identity = ["experiment_id", "group_id", "source", "well"]
    if frame.duplicated(identity).any():
        raise ValueError("Duplicate experiment-membership observations")
    return frame.sort_values(["species", "experiment_id", "group_id", "plate_id", "timestamp", "source", "well"]).reset_index(drop=True)
