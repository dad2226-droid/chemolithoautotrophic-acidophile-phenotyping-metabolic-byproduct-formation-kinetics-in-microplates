"""Vessels."""
from __future__ import annotations
from decimal import Decimal
from pathlib import Path
from typing import Any
import numpy as np
import pandas as pd
import yaml
from .metadata import load_metadata


class VesselDataError(ValueError):
    """A measurement file is inconsistent with its declared design or statuses."""


_FLASK_DIR = "data/04_flask_validation"


_TUBE_DIR = "data/05_test_tube_validation"


_NOMINAL_COLUMNS = (
    "nominal_Ga_mM", "nominal_In_mM", "nominal_Cu_mM", "nominal_Al_mM",
    "nominal_V_mM", "nominal_NaCl_mM", "nominal_citrate_mM",
)


def _metadata(root: str | Path, directory: str, filename: str) -> dict[str, Any]:
    path = Path(root) / directory / filename
    result = load_metadata(path)
    if (not isinstance(result, dict)
            or ("schema_version" in result and result["schema_version"] != "1.0")):
        raise VesselDataError(f"Unsupported vessel descriptor: {path}")
    for key, kind in {"files": dict, "columns": list, "data_format": dict,
                      "conditions": list}.items():
        if not isinstance(result.get(key), kind) or not result[key]:
            raise VesselDataError(f"Missing or invalid vessel {key}: {path}")
    if not isinstance(result["files"].get("measurements"), str):
        raise VesselDataError(f"Missing vessel measurements filename: {path}")
    if not isinstance(result["data_format"].get("unique_key"), list):
        raise VesselDataError(f"Missing vessel measurement key: {path}")
    return result


def get_flask_metadata(root: str | Path) -> dict[str, Any]:
    """Return the flask descriptor without modifying or normalizing its values."""
    return _metadata(root, _FLASK_DIR, "flask_validation.yaml")


def get_tube_metadata(root: str | Path) -> dict[str, Any]:
    """Return tube metadata, including the separate qualitative LF follow-up."""
    return _metadata(root, _TUBE_DIR, "test_tube_validation.yaml")


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise VesselDataError(message)


def _equal(left: Any, right: Any) -> bool:
    """Compare metadata literal NA with pandas missing values, never with zero."""
    right_missing = right is None or right == "NA"
    if pd.isna(left):
        return bool(right_missing)
    return False if right_missing else bool(left == right)


def _read_measurements(root: str | Path, directory: str, meta: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    path = Path(root) / directory / meta["files"]["measurements"]
    raw = pd.read_csv(path, dtype=str, keep_default_na=False, encoding="utf-8")
    columns = meta["columns"]
    names = [item["name"] for item in columns]
    _require(len(names) == len(set(names)), "Duplicate descriptor column names")
    _require(raw.columns.tolist() == names, f"CSV columns disagree with descriptor: {path}")
    _require(not raw.eq("").any().any(), f"Empty CSV fields must be explicit NA: {path}")
    frame = raw.mask(raw.eq("NA"))
    for definition in columns:
        name = definition["name"]
        if definition["type"] in ("number", "integer"):
            try:
                numeric = pd.to_numeric(frame[name], errors="raise")
            except (ValueError, TypeError) as exc:
                raise VesselDataError(f"Invalid numeric entry in {name}") from exc
            _require(np.isfinite(numeric.dropna()).all(), f"Nonfinite numeric entry in {name}")
            if definition["type"] == "integer":
                _require(numeric.notna().all() and numeric.eq(np.floor(numeric)).all(),
                         f"Missing or noninteger identifier in {name}")
                numeric = numeric.astype("int64")
            frame[name] = numeric
        if "allowed_values" in definition:
            allowed = definition["allowed_values"]
            allowed_mask = frame[name].isin(allowed)
            if "NA" in allowed:
                allowed_mask |= frame[name].isna()
            _require(allowed_mask.all(),
                     f"Unsupported status/value in {name}")
    key = meta["data_format"]["unique_key"]
    _require(frame[key].notna().all().all(), "Missing vessel record key")
    _require(not frame.duplicated(key).any(), f"Duplicate record key: {key}")
    _require(frame["condition_ID"].notna().all(), "Missing condition_ID")
    return frame, raw


def _index_conditions(meta: dict) -> dict[str, dict]:
    conditions = meta["conditions"]
    result = {item["condition_ID"]: item for item in conditions}
    _require(len(result) == len(conditions), "Duplicate condition_ID in descriptor")
    return result


def _validate_count(meta: dict, key: str, observed: int) -> None:
    if key in meta.get("counts", {}):
        _require(meta["counts"][key] == observed,
                 f"Declared count {key}={meta['counts'][key]}, observed {observed}")


def _status_matches(frame: pd.DataFrame, value: str, status: str, numeric_statuses: set[str]) -> None:
    _require(frame[value].notna().eq(frame[status].isin(numeric_statuses)).all(),
             f"Numerical availability disagrees with {status}")


def _flask_additives(condition: dict) -> dict[str, float]:
    result = dict.fromkeys(_NOMINAL_COLUMNS, 0.0)
    aliases = {"sodium chloride": "nominal_NaCl_mM", "citrate": "nominal_citrate_mM",
               "citric acid": "nominal_citrate_mM"}
    seen = set()
    for item in condition["additives"]:
        key = aliases.get(item["name"])
        if key is None and "measured_element" in item:
            key = f"nominal_{item['measured_element']}_mM"
        _require(key in result, f"Unsupported additive definition: {item}")
        _require(key not in seen, f"Duplicate additive component: {key}")
        concentration = float(item["nominal_concentration_mM"])
        _require(np.isfinite(concentration) and concentration >= 0,
                 f"Invalid nominal concentration: {item}")
        result[key] = concentration
        seen.add(key)
    if "nominal_citric_acid_to_Ga_molar_ratio" in condition:
        _require(result["nominal_Ga_mM"] > 0, "Citrate/Ga ratio is undefined without Ga")
        _require(np.isclose(condition["nominal_citric_acid_to_Ga_molar_ratio"],
                            result["nominal_citrate_mM"] / result["nominal_Ga_mM"]),
                 "Nominal citrate/Ga ratio disagrees with additions")
    return result


def _od600_values(frame: pd.DataFrame, rule: dict) -> tuple[pd.Series, pd.Series]:
    """Apply the documented proxy without changing source RFU or missingness."""
    input_status = frame["planktonic_fluorescence_status"]
    rfu = frame[rule["input_column"]]
    expected = pd.Series(np.nan, index=frame.index, dtype=float)
    eligible = rfu.notna() & rfu.ge(0)
    expected.loc[eligible] = ((rfu.loc[eligible] / rule["coefficient"]) **
                              (1 / rule["exponent"]) - rule["baseline_OD600"])
    expected_status = input_status.copy()
    expected_status.loc[rfu.lt(0)] = "negative_rfu"
    expected_status.loc[eligible & expected.lt(0)] = "below_baseline"
    expected_status.loc[eligible & expected.ge(0)] = "calculated"
    return expected, expected_status


def _verify_od600(frame: pd.DataFrame, raw: pd.DataFrame, meta: dict) -> None:
    """Check stored apparent OD600 values against the fluorescence calibration."""
    rule = meta["planktonic_OD600_calculation"]
    value, output, status = rule["input_column"], rule["output_column"], rule["status_column"]
    expected, expected_status = _od600_values(frame, rule)
    eligible = frame[value].notna() & frame[value].ge(0)
    matching_status = frame[status].eq(expected_status) | (frame[status].isna() & expected_status.isna())
    _require(matching_status.all(), "OD600 status disagrees with calculation rule")
    _require(frame[output].notna().eq(eligible).all(), "OD600 numerical availability disagrees with RFU")
    # Half one unit of the final *stored* decimal place, including scientific
    # notation. Machine-roundoff allowance does not imply measurement precision.
    for index in frame.index[eligible]:
        text = raw.at[index, output]
        unit = float(Decimal(10) ** Decimal(text).as_tuple().exponent)
        tolerance = unit / 2 + 8 * np.finfo(float).eps * max(1.0, abs(expected.at[index]))
        _require(abs(frame.at[index, output] - expected.at[index]) <= tolerance,
                 f"Stored OD600 proxy disagrees with formula at CSV row {index + 2}")


def load_flasks(root: str | Path) -> pd.DataFrame:
    """Validate every flask record and add analysis quantities in memory.

    ``eligible_feii_mM`` excludes every vanadyl condition and every non-reported
    result. ``derived_feiii_mM`` is total soluble Fe minus eligible Fe(II) only
    on the *same row*. Negative differences remain numerical and are flagged;
    no interpolation, initial-total-Fe substitution or clipping is performed.
    OD600 equivalents use the fluorescence calibration specified in metadata.
    When the input contains calculated OD600 values, the loader checks them
    against the calibration specified by its descriptor.
    """
    meta = get_flask_metadata(root)
    frame, raw = _read_measurements(root, _FLASK_DIR, meta)
    conditions = _index_conditions(meta)
    _require(set(frame.condition_ID) == set(conditions), "Flask condition coverage disagrees with descriptor")
    common = meta["common_experimental_conditions"]
    assignments, metadata_rows = {}, {}
    for condition_id, condition in conditions.items():
        kind = condition["biological_status"]
        _require(kind in {"biotic", "abiotic"}, f"Invalid biological_status: {kind}")
        additives = _flask_additives(condition)
        metadata_rows[condition_id] = {
            "organism": common["organism_abbreviation"] if kind == "biotic" else "Ab",
            "kind": kind, "batch": str(condition["start_date"]), "medium": common["medium"],
            "incubation_temperature_C": common["incubation_temperature_C"],
            "vanadyl_interference": (condition.get("cerimetry_interpretation") == "vanadyl_interference"
                                     or any("vanadyl" in item["name"].lower()
                                            for item in condition["additives"])),
            **additives,
        }
        _require(len({item["replicate"] for item in condition["flasks"]}) == len(condition["flasks"]),
                 f"Duplicate replicate assignment: {condition_id}")
        for flask in condition["flasks"]:
            identity = flask["flask_ID"]
            _require(identity not in assignments, f"Flask assigned more than once: {identity}")
            culture = flask.get("inoculum_culture_id")
            _require((culture not in {None, "NA"}) == (kind == "biotic"),
                     f"Culture assignment disagrees with biological status: {identity}")
            assignments[identity] = (condition_id, flask)
    _require(set(frame.flask_ID) == set(assignments), "Flask identity coverage disagrees with descriptor")
    for identity, rows in frame.groupby("flask_ID", sort=False):
        condition_id, assignment = assignments[identity]
        _require(rows.condition_ID.eq(condition_id).all(), f"Condition mismatch for {identity}")
        _require(rows.replicate.eq(assignment["replicate"]).all(), f"Replicate mismatch for {identity}")
        _require(rows.inoculum_culture_id.map(lambda value: _equal(value, assignment["inoculum_culture_id"])).all(),
                 f"Culture mismatch for {identity}")
        count = conditions[condition_id]["sampling_events_per_flask"]
        _require(set(rows.sampling_event) == set(range(1, count + 1)), f"Missing/extra sampling events for {identity}")
        _require(rows.loc[rows.sampling_event.eq(1), "elapsed_days_mean"].eq(0).all(),
                 f"Initial flask event must be time zero: {identity}")
    _status_matches(frame, "cerimetric_FeII_equivalent_mM", "cerimetry_status", {"reported", "vanadyl_interference"})
    for fraction in ("planktonic", "nonplanktonic"):
        _status_matches(frame, f"measured_{fraction}_fluorescence_RFU", f"{fraction}_fluorescence_status", {"reported"})
    stored_od600 = "planktonic_OD600_calculation" in meta
    if stored_od600:
        _verify_od600(frame, raw, meta)
        od600_rule = meta["planktonic_OD600_calculation"]
    else:
        rule_path = Path(__file__).resolve().parents[1] / "metadata/flask_od600_calculation.yaml"
        with rule_path.open(encoding="utf-8") as handle:
            od600_rule = yaml.safe_load(handle)
        output, status = od600_rule["output_column"], od600_rule["status_column"]
        _require(output not in frame and status not in frame,
                 "Stored OD600 columns require a fluorescence calibration descriptor")
        frame[output], frame[status] = _od600_values(frame, od600_rule)
    joined = pd.DataFrame.from_dict(metadata_rows, orient="index")
    for column in joined:
        _require(column not in frame, f"Derived metadata overlaps an input column: {column}")
        frame[column] = frame.condition_ID.map(joined[column])
    flagged = frame.cerimetry_status.eq("vanadyl_interference")
    _require((~flagged | frame.vanadyl_interference).all(), "Vanadyl status lacks matching condition metadata")
    # The vanadyl condition flag applies to all measurements from that condition.
    frame["feii_eligible"] = frame.cerimetry_status.eq("reported") & ~frame.vanadyl_interference
    frame["eligible_feii_mM"] = frame.cerimetric_FeII_equivalent_mM.where(frame.feii_eligible)
    frame["derived_feiii_mM"] = frame.soluble_total_Fe_mM - frame.eligible_feii_mM
    frame["feiii_negative_difference"] = frame.derived_feiii_mM.lt(0)
    _validate_count(meta, "conditions", len(conditions))
    _validate_count(meta, "flasks", frame.flask_ID.nunique())
    _validate_count(meta, "flask_sampling_records", len(frame))
    _validate_count(meta, "inoculum_cultures", frame.inoculum_culture_id.nunique())
    fluorescence = meta.get("measurement_methods", {}).get("fluorescence", {})
    initial = meta.get("initial_chemistry_assignments", {})
    # Initial chemistry is shared when the descriptor specifies assignment columns
    # and does not identify the measurements as independent for each flask.
    assigned = bool(initial.get("columns")) and initial.get(
        "independently_measured_in_each_flask", False) is False
    if assigned:
        for column in initial["columns"]:
            _require(column in frame, f"Unknown initial chemistry assignment column: {column}")
            frame[f"{column}_is_initial_assignment"] = (
                frame.sampling_event.eq(initial["sampling_event"]) & frame[column].notna())
    frame.attrs.update(metadata=meta, source_file=f"{_FLASK_DIR}/{meta['files']['measurements']}",
                       fluorescence_background_corrected=fluorescence.get("background_corrected", False),
                       fluorescence_background_blank=fluorescence.get("background_blank"),
                       fluorescence_blank_subtracted_before_deposition=fluorescence.get(
                           "blank_subtracted_before_deposition", False),
                       od600_proxy_verified=stored_od600,
                       od600_proxy_origin="verified input values" if stored_od600 else "computed in analysis",
                       od600_calculation=od600_rule, experimental_unit="physical flask",
                       initial_chemistry_assignments=initial,
                       baseline_independence=("assigned starting-medium chemistry; not independently measured in each flask"
                                              if assigned else "not established"))
    return frame


def _tube_additives(treatment: dict) -> dict[str, float]:
    result = dict.fromkeys(_NOMINAL_COLUMNS, 0.0)
    if "ion" in treatment:
        element = treatment["ion"].split("(")[0]
        key = f"nominal_{element}_mM"
        _require(key in result, f"Unsupported tube treatment ion: {element}")
        result[key] = float(treatment["final_concentration_mM"])
    else:
        _require("Ga_final_mM" in treatment and "In_final_mM" in treatment,
                 "Treatment must explicitly define its additions")
        result["nominal_Ga_mM"] = float(treatment["Ga_final_mM"])
        result["nominal_In_mM"] = float(treatment["In_final_mM"])
    _require(all(np.isfinite(value) and value >= 0 for value in result.values()), "Invalid tube nominal concentration")
    return result


def load_tubes(root: str | Path) -> pd.DataFrame:
    """Load tube time courses and retain shared references without pseudoreplication.

    RFU is already background corrected and is copied unchanged. The loader
    does not derive OD600, Fe(II), Fe(III), or proton amounts from tube data.
    ``is_independent_measurement`` distinguishes observed rows from repeated
    reference assignments; it does NOT make repeated samples of one tube
    statistically independent experimental replicates.
    """
    meta = get_tube_metadata(root)
    frame, _ = _read_measurements(root, _TUBE_DIR, meta)
    conditions = _index_conditions(meta)
    _require(set(frame.condition_ID) == set(conditions), "Tube condition coverage disagrees with descriptor")
    assignments, metadata_rows, block_by_tube = {}, {}, {}
    for block_id, block in meta["timing_blocks"].items():
        prefix, first = block["first_tube_id"].rsplit("-", 1)
        last_prefix, last = block["last_tube_id"].rsplit("-", 1)
        _require(prefix == last_prefix, f"Inconsistent tube range for {block_id}")
        for number in range(int(first), int(last) + 1):
            identity = f"{prefix}-{number}"
            _require(identity not in block_by_tube, f"Overlapping timing blocks: {identity}")
            block_by_tube[identity] = block_id
    abiotic_code = meta["abiotic_code"]["code"]
    for condition_id, condition in conditions.items():
        organism = condition["organism_code"]
        _require(organism == abiotic_code or organism in meta["organisms"], f"Unknown organism {organism}")
        _require(condition["medium_ID"] in meta["media"], f"Unknown medium for {condition_id}")
        _require(condition["treatment_ID"] in meta["treatments"], f"Unknown treatment for {condition_id}")
        metadata_rows[condition_id] = {
            "organism": organism, "kind": "abiotic" if organism == abiotic_code else "biotic",
            "batch": meta["common_conditions"]["start_datetime"], "medium": condition["medium_ID"],
            **_tube_additives(meta["treatments"][condition["treatment_ID"]]),
        }
        identities, replicates = condition["tube_ids"], condition["replicate_order"]
        _require(len(identities) == len(replicates) and len(replicates) == len(set(replicates)),
                 f"Invalid tube replicate assignments for {condition_id}")
        for identity, replicate in zip(identities, replicates):
            _require(identity not in assignments, f"Tube assigned more than once: {identity}")
            _require(identity in block_by_tube, f"No timing block for {identity}")
            culture = None
            if organism != abiotic_code:
                cultures = meta["replication"]["culture_ids"][organism]
                _require(1 <= replicate <= len(cultures), f"Undefined culture replicate for {identity}")
                culture = cultures[replicate - 1]
            assignments[identity] = (condition_id, replicate, culture)
    _require(set(frame.tube_ID) == set(assignments) == set(block_by_tube), "Tube identity coverage disagrees with descriptor")
    for identity, rows in frame.groupby("tube_ID", sort=False):
        condition_id, replicate, culture = assignments[identity]
        condition = conditions[condition_id]
        for name in ("condition_ID", "organism_code", "medium_ID", "treatment_ID", "incubation_temperature_C"):
            expected = condition_id if name == "condition_ID" else condition[name]
            _require(rows[name].eq(expected).all(), f"{name} mismatch for {identity}")
        _require(rows.replicate.eq(replicate).all(), f"Replicate mismatch for {identity}")
        _require(rows.inoculum_culture_ID.map(lambda value: _equal(value, culture)).all(), f"Culture mismatch for {identity}")
        block_id = block_by_tube[identity]
        _require(rows.timing_block_ID.eq(block_id).all(), f"Timing block mismatch for {identity}")
        last_event = meta["timing_blocks"][block_id]["poststart_event_count"]
        _require(set(rows.sampling_event) == set(range(last_event + 1)), f"Missing/extra sampling events for {identity}")
        expected_roles = np.where(rows.sampling_event.eq(0), "shared_initial_reference",
                                  np.where(rows.sampling_event.eq(last_event), "endpoint", "measurement"))
        _require(rows.record_role.eq(expected_roles).all(), f"Incorrect record roles for {identity}")
    references = meta["shared_initial_references"]["values"]
    reference_index = {item["initial_reference_ID"]: item for item in references}
    _require(len(reference_index) == len(references), "Duplicate initial_reference_ID in descriptor")
    is_initial = frame.record_role.eq("shared_initial_reference")
    _require(frame.loc[~is_initial, "initial_reference_ID"].isna().all(), "Reference ID on a measured tube row")
    _require(set(frame.loc[is_initial, "initial_reference_ID"]) == set(reference_index), "Shared initial reference coverage mismatch")
    for row in frame.loc[is_initial].itertuples(index=False):
        reference = reference_index[row.initial_reference_ID]
        for field, expected in (("timing_block_ID", reference["timing_block_ID"]),
                                ("treatment_ID", reference["treatment_ID"]),
                                ("replicate", reference["reference_index"]),
                                ("pH", reference["pH"]), ("orp_mV", reference["orp_mV"])):
            _require(_equal(getattr(row, field), expected), f"Shared reference mismatch: {row.tube_ID}, {field}")
        _require(row.time_average_d == 0 and row.time_sd_d == 0 and row.datetime_basis == "assigned_start",
                 f"Invalid shared-reference timing: {row.tube_ID}")
    _require(frame.datetime_basis.isin({"recorded", "interpolated", "assigned_start"}).all(), "Invalid datetime_basis")
    _require(~frame.loc[~is_initial, "datetime_basis"].eq("assigned_start").any(), "Assigned start on measured row")
    _status_matches(frame, "planktonic_background_corrected_rfu", "fluorescence_status", {"reported"})
    endpoint = frame.record_role.eq("endpoint")
    _require(frame.fluorescence_status.eq("not_scheduled").eq(~endpoint).all(), "Fluorescence scheduling disagrees with record role")
    joined = pd.DataFrame.from_dict(metadata_rows, orient="index")
    for column in joined:
        _require(column not in frame, f"Derived metadata overlaps an input column: {column}")
        frame[column] = frame.condition_ID.map(joined[column])
    frame["is_shared_initial_reference"] = is_initial
    frame["is_independent_measurement"] = ~is_initial
    _validate_count(meta, "conditions", len(conditions))
    _validate_count(meta, "tubes", frame.tube_ID.nunique())
    _validate_count(meta, "total_csv_rows", len(frame))
    _validate_count(meta, "time_course_rows", int((~is_initial).sum()))
    _validate_count(meta, "shared_initial_reference_rows", int(is_initial.sum()))
    _validate_count(meta, "distinct_initial_reference_assignments", len(reference_index))
    _validate_count(meta, "independent_starter_cultures", frame.inoculum_culture_ID.nunique())
    frame.attrs.update(metadata=meta, source_file=f"{_TUBE_DIR}/{meta['files']['measurements']}",
                       fluorescence_background_corrected=True, initial_reference_count=len(reference_index),
                       experimental_unit="physical tube")
    return frame
