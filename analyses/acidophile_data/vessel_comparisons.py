"""Vessel comparisons."""
from __future__ import annotations
import json
import numpy as np
import pandas as pd


FLASK_METRICS = {
    "eligible_feii_mM": ("mM", "eligible cerimetric Fe(II) estimate"),
    "feii_depletion_fraction": ("fraction", "fractional Fe(II) depletion from assigned starting-medium reference; reference uncertainty not propagated"),
    "soluble_total_Fe_mM": ("mM", "measured soluble total Fe"),
    "derived_feiii_mM": ("mM", "same-row total Fe minus eligible Fe(II); not a mass balance"),
    "soluble_Ga_mM": ("mM", "measured soluble Ga"),
    "soluble_Cu_mM": ("mM", "measured soluble Cu"),
    "soluble_V_mM": ("mM", "measured soluble V"),
    "soluble_Al_mM": ("mM", "measured soluble Al"),
    "measured_pH": ("pH", "measured pH; not proton production"),
    "measured_ORP_mV": ("mV", "measured electrode potential versus Ag/AgCl with 3 M KCl"),
    "measured_planktonic_fluorescence_RFU": ("RFU", "TE Buffer blank-subtracted planktonic DNA-associated fluorescence"),
    "measured_nonplanktonic_fluorescence_RFU": ("RFU", "TE Buffer blank-subtracted nonplanktonic DNA-associated fluorescence"),
    "calculated_planktonic_OD600_equivalent": ("OD600-equivalent", "documented fluorescence-derived proxy; not viable biomass"),
}


def _ids(values: pd.Series) -> str:
    """Serialize identities without inventing identities for abiotic cultures."""
    return json.dumps(sorted(str(value) for value in values.dropna().unique()))


def _constant(group: pd.DataFrame, column: str):
    values = group[column].dropna().unique()
    if len(values) > 1:
        raise ValueError(f"Expected constant {column} within condition/event; found {values}")
    return values[0] if len(values) else np.nan


def flask_feii_depletion(flasks: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return per-observation depletion and per-flask first-90% brackets.

    Only a valid, positive, deposited time-zero Fe(II) value can establish the
    denominator. The descriptor assigns this value from starting-medium
    chemistry, not measured independently in each flask. Reference uncertainty
    is not propagated. Later measurements never stand in for a missing baseline.
    ``initial_reference_provenance`` identifies the source of ``initial_observed_feii_mM``.
    The bracket uses reported elapsed-time means; these are wider-series
    timing summaries, not individually recorded flask acquisition times.
    """
    observations = flasks.copy()
    # Source descriptors are returned separately by the loaders. Repeated
    # pandas slicing would otherwise deep-copy the full descriptor thousands
    # of times during numerical aggregation without adding provenance.
    observations.attrs = {}
    observations["initial_observed_feii_mM"] = np.nan
    observations["initial_reference_provenance"] = "no_initial_reference"
    observations["feii_depletion_fraction"] = np.nan
    observations["feii_depletion_status"] = "no_valid_initial_reference"
    intervals = []
    for flask_id, group in observations.groupby("flask_ID", sort=True):
        group = group.sort_values(["elapsed_days_mean", "sampling_event"])
        initial = group.loc[group["elapsed_days_mean"].eq(0)]
        if len(initial) > 1:
            raise ValueError(f"Multiple initial measurements for flask {flask_id}")
        baseline = initial["eligible_feii_mM"].iloc[0] if len(initial) else np.nan
        baseline_ok = pd.notna(baseline) and np.isfinite(baseline) and baseline > 0
        assignment_flag = "cerimetric_FeII_equivalent_mM_is_initial_assignment"
        assigned = assignment_flag in initial and initial[assignment_flag].any()
        reference_provenance = ("assigned_starting_medium; not independently measured per flask"
                                if assigned else "deposited initial; independence not established")
        observations.loc[group.index, "initial_reference_provenance"] = reference_provenance
        result = {
            "flask_ID": flask_id,
            "condition_ID": _constant(group, "condition_ID"),
            "inoculum_culture_id": _constant(group, "inoculum_culture_id"),
            "initial_observed_feii_mM": baseline if baseline_ok else np.nan,
            "initial_reference_provenance": reference_provenance,
            "depletion_threshold_fraction": 0.90,
            "interval_lower_d": np.nan,
            "interval_upper_d": np.nan,
            "lower_sampling_event": np.nan,
            "upper_sampling_event": np.nan,
            "last_valid_feii_time_d": np.nan,
            "last_observation_time_d": group["elapsed_days_mean"].max(),
            "censoring_status": "not_estimable_no_valid_initial",
            "timing_basis": "reported wider-series elapsed-time means; not exact flask times",
            "interpretation": "measured Fe(II) depletion; not growth lag or inferred biological fraction",
        }
        for column in ("batch", "kind", "organism", "nominal_citrate_mM", "nominal_Ga_mM"):
            if column in group:
                result[column] = _constant(group, column)
        if baseline_ok:
            observations.loc[group.index, "initial_observed_feii_mM"] = baseline
            values = pd.to_numeric(group["eligible_feii_mM"], errors="coerce")
            valid = values.notna() & np.isfinite(values)
            fractions = (baseline - values) / baseline
            observations.loc[group.index, "feii_depletion_fraction"] = fractions.where(valid)
            observations.loc[group.index, "feii_depletion_status"] = np.where(
                valid, "estimated_from_assigned_initial" if assigned else "estimated_from_deposited_initial",
                "ineligible_or_unavailable_FeII"
            )
            valid_group = group.loc[valid].copy()
            valid_group["fraction"] = fractions.loc[valid]
            result["last_valid_feii_time_d"] = valid_group["elapsed_days_mean"].max()
            reached = valid_group.loc[valid_group["fraction"].ge(0.90)]
            if not reached.empty:
                first = reached.iloc[0]
                earlier = valid_group.loc[valid_group["elapsed_days_mean"].lt(first["elapsed_days_mean"])]
                previous = earlier.iloc[-1] if not earlier.empty else None
                result.update(
                    interval_lower_d=previous["elapsed_days_mean"] if previous is not None else 0.0,
                    interval_upper_d=first["elapsed_days_mean"],
                    lower_sampling_event=previous["sampling_event"] if previous is not None else np.nan,
                    upper_sampling_event=first["sampling_event"],
                    censoring_status="interval_censored_first_observed_reaching_threshold",
                )
            elif len(valid_group) > 1:
                last = valid_group.iloc[-1]
                result.update(
                    interval_lower_d=last["elapsed_days_mean"],
                    lower_sampling_event=last["sampling_event"],
                    censoring_status="right_censored_not_reached_at_last_valid_observation",
                )
            else:
                result["censoring_status"] = "not_estimable_no_valid_followup"
        intervals.append(result)
    return observations, pd.DataFrame(intervals)


def _summarize(
    frame: pd.DataFrame,
    *,
    vessel_column: str,
    culture_column: str,
    time_column: str,
    time_sd_column: str,
    metrics: dict[str, tuple[str, str]],
) -> pd.DataFrame:
    """Long-form sample summaries, preserving the actual sampling schedule."""
    records = []
    for (condition, event), group in frame.groupby(["condition_ID", "sampling_event"], sort=True):
        if group[vessel_column].duplicated().any():
            raise ValueError(f"Duplicate vessel measurements in {condition}, event {event}")
        common = {
            "condition_ID": condition,
            "sampling_event": event,
            "elapsed_time_mean_d": _constant(group, time_column),
            "elapsed_time_sd_d": _constant(group, time_sd_column),
            "n_vessels_in_event": group[vessel_column].nunique(),
            "vessel_ids": _ids(group[vessel_column]),
            "inoculum_culture_ids": _ids(group[culture_column]),
        }
        for column in (
            "batch", "kind", "organism", "medium", "medium_ID", "treatment_ID",
            "incubation_temperature_C", "nominal_Ga_mM", "nominal_citrate_mM",
            "nominal_Cu_mM", "nominal_Al_mM", "nominal_V_mM", "nominal_NaCl_mM",
        ):
            if column in group:
                common[column] = _constant(group, column)
        for metric, (unit, interpretation) in metrics.items():
            if metric not in group:
                continue
            numeric = pd.to_numeric(group[metric], errors="coerce")
            usable = numeric.notna() & np.isfinite(numeric)
            values = numeric.loc[usable]
            records.append({
                **common,
                "metric": metric,
                "unit": unit,
                "mean": values.mean(),
                "sd": values.std(ddof=1),
                "n": len(values),
                "contributing_vessel_ids": _ids(group.loc[usable, vessel_column]),
                "contributing_culture_ids": _ids(group.loc[usable, culture_column]),
                "summary_status": "reported" if len(values) else "no_eligible_measurements",
                "interpretation": interpretation,
            })
    return pd.DataFrame(records)
