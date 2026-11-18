"""Plate comparisons."""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from .metadata import load_metadata


CITRATE_EXPERIMENT = "AF_AFM1_gallium_sulfate_OD0p001_citrate_ratios"


CROSS_FORMAT_EXPERIMENTS = frozenset(
    "AF_AFM1_" + compound + "_OD0p001" for compound in (
        "gallium_sulfate", "copper_ii_sulfate", "sodium_chloride",
        "aluminum_sulfate", "vanadyl_sulfate", "citric_acid",
    )
)


RATIO_DENOMINATOR_MIN_AU = 1e-6


CONDITION_KEYS = ["species", "experiment_id", "group_id", "column"]


SERIES_KEYS = CONDITION_KEYS + ["plate_id", "batch_id", "culture_id", "row", "signal"]


SAMPLE_KEYS = CONDITION_KEYS + ["plate_id", "batch_id", "source", "ordinal", "days"]


MEASURES = ["biotic_signal", "abiotic_signal", "biotic_change", "abiotic_change", "net_change"]


CONDITION_METADATA = [
    "condition_json", "component_units_json", "component_basis_json",
    "component_preparation_mM_json", "component_final_mM_json", "medium_id",
    "preparation_id", "readout", "initial_OD600_200_uL", "initial_OD600_final",
]


def selected_experiment_ids(root: str | Path) -> list[str]:
    """Select all pH experiments, all LF conditions and complete flask-bridge series."""
    base = Path(root) / "data/03_biological_microplates"
    ids = []
    for path in sorted(base.glob("*/*_experiments.yaml")):
        document = load_metadata(path)
        for experiment in document["experiments"]:
            identifier = experiment["experiment_id"]
            if (document["species_id"] == "LF" or experiment["readout_id"] == "ph"
                    or identifier == CITRATE_EXPERIMENT
                    or identifier in CROSS_FORMAT_EXPERIMENTS):
                ids.append(identifier)
    return ids


def _safe_channel(frame: pd.DataFrame, name: str) -> pd.Series:
    values = pd.to_numeric(frame[name], errors="coerce")
    # Per-channel masks are required: an invalid unrelated channel must not
    # exclude otherwise usable A380 or indicator measurements.
    flag = name + "_eligible"
    eligible = np.isfinite(values) & values.le(4.0)
    if flag in frame:
        eligible &= frame[flag].fillna(False).astype(bool)
    return values.where(eligible)


def _ratio(numerator: pd.Series, denominator: pd.Series, cutoff: float):
    status = pd.Series("eligible", index=numerator.index, dtype=object)
    status.loc[denominator.le(cutoff) & denominator.gt(0)] = "near_zero_denominator"
    status.loc[denominator.le(0)] = "nonpositive_denominator"
    status.loc[~np.isfinite(numerator) | ~np.isfinite(denominator)] = "missing_or_ineligible_channel"
    return (numerator / denominator.where(status.eq("eligible"))), status


def derive_plate_signals(observations: pd.DataFrame,
                         denominator_min: float = RATIO_DENOMINATOR_MIN_AU) -> pd.DataFrame:
    """Return one row per physical-well membership and applicable signal.

    The denominator minimum is a numerical-stability guard, not an empirical
    detection limit. Its value and exclusions travel with the derived table.
    """
    if not np.isfinite(denominator_min) or denominator_min < 0:
        raise ValueError("Ratio denominator minimum must be finite and nonnegative")
    frame = observations.copy()
    if frame.empty:
        return pd.DataFrame()
    unique = ["experiment_id", "group_id", "physical_observation_key"]
    if frame.duplicated(unique).any():
        raise ValueError("Duplicate physical observation within an experiment/group membership")
    pieces = []
    for readout, signal_names in [("iron", ["raw_A380", "A380_minus_A750"]),
                                  ("ph", ["raw_R545_485", "corrected_R545_485", "A545_minus_A750"])]:
        group = frame[frame.readout.eq(readout)].copy()
        if group.empty:
            continue
        values = {name: _safe_channel(group, name) for name in ("A380", "A750", "A545", "A485")}
        signals = {
            "raw_A380": values["A380"],
            "A380_minus_A750": values["A380"] - values["A750"],
            "A545_minus_A750": values["A545"] - values["A750"],
        }
        ratio_status = {}
        for name, numerator, denominator in [
            ("raw_R545_485", values["A545"], values["A485"]),
            ("corrected_R545_485", values["A545"] - values["A750"], values["A485"] - values["A750"]),
        ]:
            signals[name], ratio_status[name] = _ratio(numerator, denominator, denominator_min)
        keep = list(dict.fromkeys(SAMPLE_KEYS + ["culture_id", "row", "well", "kind", "physical_observation_key"]))
        if "physical_plate_id" in group:
            keep.append("physical_plate_id")
        for name in signal_names:
            result = group[keep].copy()
            result["signal"] = name
            result["value"] = signals[name]
            result["signal_status"] = ratio_status.get(name, pd.Series(
                np.where(np.isfinite(signals[name]), "eligible", "missing_or_ineligible_channel"), index=group.index))
            result["ratio_denominator_min_AU"] = denominator_min if "R545" in name else np.nan
            pieces.append(result)
    return pd.concat(pieces, ignore_index=True) if pieces else pd.DataFrame()


def compare_biotic_with_controls(signals: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Preserve each biotic culture; require all three matched abiotic wells."""
    if signals.empty:
        return pd.DataFrame(), pd.DataFrame()
    sample_keys = SAMPLE_KEYS + ["signal"]
    controls = signals[signals.kind.eq("abiotic")].groupby(sample_keys, dropna=False, as_index=False).agg(
        abiotic_available_mean=("value", "mean"), abiotic_sd=("value", "std"),
        abiotic_n=("value", "count"), abiotic_rows=("row", "nunique"),
        abiotic_invalid_n=("signal_status", lambda s: int(s.ne("eligible").sum())),
    )
    controls["abiotic_signal"] = controls.abiotic_available_mean.where(controls.abiotic_n.eq(3) & controls.abiotic_rows.eq(3))
    controls["controls_status"] = np.where(controls.abiotic_n.eq(3) & controls.abiotic_rows.eq(3), "complete", "incomplete_matched_abiotic")
    biotic = signals[signals.kind.eq("biotic")].rename(columns={"value": "biotic_signal"})
    result = biotic.merge(controls, on=sample_keys, how="left", validate="many_to_one")
    # The registry defines the initial acquisition by its measured elapsed time;
    # a deposited ordinal of 1 is not evidence that the initial read is absent.
    baseline = result[result.days.eq(0)].copy()
    if baseline.duplicated(SERIES_KEYS).any():
        raise ValueError("Ambiguous initial source for a physical well")
    baseline = baseline[SERIES_KEYS + ["biotic_signal", "abiotic_signal", "source"]].rename(columns={
        "biotic_signal": "biotic_baseline", "abiotic_signal": "abiotic_baseline", "source": "baseline_source"})
    result = result.merge(baseline, on=SERIES_KEYS, how="left", validate="many_to_one")
    result["baseline_status"] = np.select([
        result.baseline_source.isna(), result.biotic_baseline.isna(), result.abiotic_baseline.isna()],
        ["missing_initial_acquisition", "invalid_biotic_baseline", "invalid_abiotic_baseline"], default="eligible")
    result["biotic_change"] = result.biotic_signal - result.biotic_baseline
    result["abiotic_change"] = result.abiotic_signal - result.abiotic_baseline
    result["net_change"] = (result.biotic_change - result.abiotic_change).where(result.baseline_status.eq("eligible"))
    result["comparison_status"] = np.select([
        result.baseline_status.ne("eligible"), result.biotic_signal.isna(), result.abiotic_signal.isna()],
        ["invalid_or_missing_baseline", "invalid_biotic_signal", "incomplete_matched_abiotic"], default="eligible")
    return result.sort_values(SERIES_KEYS + ["days"]).reset_index(drop=True), controls


def _average_technical(frame: pd.DataFrame, keys: list[str], time_column: str) -> pd.DataFrame:
    """Average three technical biotic wells per independent physical microplate."""
    spec = {name + "_available_mean": (name, "mean") for name in MEASURES}
    spec.update({name + "_valid_well_n": (name, "count") for name in MEASURES})
    spec.update({"days_mean": (time_column, "mean"), "days_min": (time_column, "min"),
                 "days_max": (time_column, "max"), "technical_well_n": ("row", "nunique"),
                 "culture_ids": ("culture_id", lambda s: ";".join(sorted(s.dropna().unique()))),
                 "physical_plate_ids": ("plate_id", lambda s: ";".join(sorted(s.unique())))})
    if "physical_plate_id" in frame:
        if frame.groupby(keys, dropna=False).physical_plate_id.nunique().ne(1).any():
            raise ValueError("A microplate mean must have one canonical physical_plate_id")
        spec["physical_plate_id"] = ("physical_plate_id", "first")
    result = frame.groupby(keys, dropna=False, as_index=False).agg(**spec)
    for name in MEASURES:
        result[name] = result[name + "_available_mean"].where(
            result[name + "_valid_well_n"].eq(3) & result.technical_well_n.eq(3))
    result["net_valid_well_n"] = result.net_change_valid_well_n
    if "physical_plate_id" not in result:
        # Standalone synthetic inputs may use their local plate identity only.
        # Canonical deposited observations always carry the persistent identity.
        result["physical_plate_id"] = result.plate_id
    result["biological_replicate_id"] = result.physical_plate_id
    result["experimental_unit"] = "independent_microplate"
    # These column aliases refer to the same measured well cultures.
    result["culture_id"] = result.plate_id
    result["row"] = "within_microplate_mean"
    result["technical_plate_n"] = 1
    result["net_valid_plate_n"] = result.net_change.notna().astype(int)
    result["aggregation"] = "mean_of_three_technical_wells_within_independent_microplate; complete_well_set_required"
    return result
