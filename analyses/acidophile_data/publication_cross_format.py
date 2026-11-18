"""Cross format."""
from __future__ import annotations
import json
import numpy as np
import pandas as pd
from .metadata import load_metadata


CHEMICAL_QUANTITIES = {
    "gallium_sulfate": ("Ga", 2), "gallium": ("Ga", 1),
    "aluminum_sulfate": ("Al", 2), "copper_ii_sulfate": ("Cu", 1),
    "vanadyl_sulfate": ("VO", 1), "sodium_chloride": ("NaCl", 1),
    "citric_acid": ("citric_acid", 1),
}


FLASK_QUANTITIES = {"gallium(III)": "Ga", "aluminium(III)": "Al",
                    "copper(II)": "Cu", "vanadyl(IV)": "VO",
                    "sodium chloride": "NaCl", "citric acid": "citric_acid",
                    "citrate": "citrate_unspecified_form"}


LIMITATION = (
    "Independent culture cohorts and batches; static 0.2-mL plates versus shaken 100-mL flasks; "
    "different actual observation windows; matching doses does not establish assay or vessel equivalence."
)


def _json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def harmonize_plate_quantities(values: dict, bases: dict) -> dict:
    """Map declared final-well mM onto explicit comparable chemical quantities."""
    result = {}
    for compound, concentration in values.items():
        if compound not in CHEMICAL_QUANTITIES:
            raise ValueError(f"Unreviewed cross-format chemical quantity: {compound}")
        quantity, stoichiometry = CHEMICAL_QUANTITIES[compound]
        basis = bases[compound]
        if basis == "compound_formula_units":
            factor = stoichiometry
        elif basis == quantity:
            factor = 1
        else:
            raise ValueError(f"Unreviewed concentration basis for {compound}: {basis}")
        result[quantity] = result.get(quantity, 0.0) + float(concentration) * factor
    return {key: value for key, value in result.items() if value != 0}


def _flask_registry(root):
    document = load_metadata(root / "data/04_flask_validation/flask_validation.yaml")
    common = document["common_experimental_conditions"]
    rows = []
    for item in document["conditions"]:
        quantities = {FLASK_QUANTITIES[addition["name"]]: float(addition["nominal_concentration_mM"])
                      for addition in item["additives"] if addition["nominal_concentration_mM"] != 0}
        rows.append(dict(flask_condition_ID=item["condition_ID"], flask_kind=item["biological_status"],
            flask_batch=item["start_date"], flask_quantities_mM_json=_json(quantities),
            flask_vessel_ids=_json([v["flask_ID"] for v in item["flasks"]]),
            flask_culture_ids=_json(sorted({v["inoculum_culture_id"] for v in item["flasks"]
                if v.get("inoculum_culture_id") not in (None, "NA")})),
            flask_temperature_C=common["incubation_temperature_C"],
            flask_initial_volume_mL=common["initial_liquid_volume_mL"],
            flask_mixing_rpm=common["shaking_speed_rpm"], flask_medium=common["medium"],
            flask_initial_target_OD600=(common.get("initial_target_OD600", np.nan)
                                        if item["biological_status"] == "biotic" else np.nan),
            flask_inoculum_OD600_status=("reported preparation target; distinct from fluorescence-derived OD600"
                if item["biological_status"] == "biotic" and "initial_target_OD600" in common else
                "not applicable: sterile medium added" if item["biological_status"] == "abiotic" else
                "not specified in deposited flask methods")))
    return pd.DataFrame(rows)


def build_correspondence(plate: pd.DataFrame, flasks: pd.DataFrame) -> pd.DataFrame:
    """Equal doses or explicit bracketing doses, plus complete unmatched coverage."""
    records, used = [], set()
    quantities = plate.comparable_quantities_mM_json.map(json.loads)
    for _, flask in flasks[flasks.flask_kind.eq("biotic")].iterrows():
        target = json.loads(flask.flask_quantities_mM_json)
        candidates = plate[quantities.map(lambda value: set(value) == set(target))]
        exact = candidates[candidates.comparable_quantities_mM_json.map(
            lambda value: all(np.isclose(json.loads(value)[k], v, rtol=0, atol=1e-10)
                              for k, v in target.items()))]
        selected = [(idx, "equal_reported_dose") for idx in exact.index]
        # Dose brackets are explicit context only. Multicomponent formulations
        # require exact vectors; no scalar distance hides an unmatched component.
        if not selected and len(target) == 1 and not candidates.empty:
            quantity, dose = next(iter(target.items()))
            candidate_doses = candidates.comparable_quantities_mM_json.map(lambda value: json.loads(value)[quantity])
            for direction, subset in (("lower", candidate_doses[candidate_doses < dose]),
                                      ("upper", candidate_doses[candidate_doses > dose])):
                if not subset.empty:
                    nearest = subset.max() if direction == "lower" else subset.min()
                    selected.extend((idx, f"nearest_{direction}_dose_context") for idx in subset[subset.eq(nearest)].index)
        abiotic = flasks[flasks.flask_kind.eq("abiotic") & flasks.flask_quantities_mM_json.eq(flask.flask_quantities_mM_json)]
        common = {**flask.to_dict(), "flask_abiotic_condition_ids": _json(abiotic.flask_condition_ID.tolist()),
                  "flask_abiotic_batches": _json(sorted(abiotic.flask_batch.unique())),
                  "shared_culture_or_batch_pairing": False, "format_equivalence_test": False,
                  "limitation": LIMITATION}
        if selected:
            for idx, status in selected:
                used.add(idx)
                records.append({**common, **plate.loc[idx].to_dict(), "match_status": status})
        else:
            records.append({**common, "match_status": "no_comparable_plate_formulation"})
    for idx, row in plate.iterrows():
        if idx not in used:
            records.append({**row.to_dict(), "match_status": "plate_dose_without_flask_correspondence",
                "shared_culture_or_batch_pairing": False, "format_equivalence_test": False, "limitation": LIMITATION})
    result = pd.DataFrame(records)
    result.insert(0, "correspondence_id", [f"CF{i:04d}" for i in range(1, len(result) + 1)])
    return result


def temporal_shape_error(plate_times, plate_values, flask_times, flask_values):
    """Time-weighted shape errors inside observed overlap; no extrapolation."""
    arrays = [np.asarray(value, dtype=float) for value in (plate_times, plate_values, flask_times, flask_values)]
    pt, py, ft, fy = arrays
    if len(pt) != len(py) or len(ft) != len(fy):
        raise ValueError("Observation times and values must have matching lengths")
    result = dict(status="insufficient_observations", overlap_start_d=np.nan, overlap_end_d=np.nan,
                  plate_source_acquisition_n=len(pt), flask_source_acquisition_n=len(ft),
                  plate_unavailable_source_acquisition_n=int((~np.isfinite(py)).sum()),
                  flask_unavailable_source_acquisition_n=int((~np.isfinite(fy)).sum()),
                  plate_observations_used=0, flask_observations_used=0,
                  unavailable_window_bracketing_acquisition_n=0,
                  time_weighted_mae=np.nan, time_weighted_rmse=np.nan,
                  plate_window_change=np.nan, flask_window_change=np.nan)
    if len(pt) < 2 or len(ft) < 2:
        return result
    if not np.isfinite(pt).all() or not np.isfinite(ft).all():
        result["status"] = "invalid_observation_times"
        return result
    if np.any(np.diff(pt) <= 0) or np.any(np.diff(ft) <= 0):
        raise ValueError("Temporal correspondence requires strictly increasing unique observation times")
    start, end = max(pt[0], ft[0]), min(pt[-1], ft[-1])
    result.update(overlap_start_d=start, overlap_end_d=end)
    if end <= start:
        result["status"] = "no_positive_common_window"
        return result
    # Include the adjacent observed values required for interpolation at an
    # overlap boundary, but do not count unused samples outside the window.
    for name, times, values in (("plate", pt, py), ("flask", ft, fy)):
        lower = max(0, int(np.searchsorted(times, start, side="right")) - 1)
        upper = min(len(times) - 1, int(np.searchsorted(times, end, side="left")))
        result[f"{name}_observations_used"] = upper - lower + 1
        result["unavailable_window_bracketing_acquisition_n"] += int(
            (~np.isfinite(values[lower:upper + 1])).sum())
    if result["unavailable_window_bracketing_acquisition_n"]:
        result["status"] = "unavailable_acquisition_in_common_window"
        return result
    grid = np.unique(np.r_[start, end, pt[(pt >= start) & (pt <= end)], ft[(ft >= start) & (ft <= end)]])
    pvalues, fvalues = np.interp(grid, pt, py), np.interp(grid, ft, fy)
    pchange, fchange = pvalues[-1] - pvalues[0], fvalues[-1] - fvalues[0]
    result.update(plate_window_change=pchange, flask_window_change=fchange)
    if pchange <= 0 or fchange <= 0:
        result["status"] = "nonpositive_observed_window_change"
        return result
    delta = (pvalues - pvalues[0]) / pchange - (fvalues - fvalues[0]) / fchange
    result.update(status="descriptive_shape_comparison", union_grid_n=len(grid),
                  time_weighted_mae=np.trapezoid(np.abs(delta), grid) / (end - start),
                  time_weighted_rmse=np.sqrt(np.trapezoid(delta ** 2, grid) / (end - start)))
    return result
