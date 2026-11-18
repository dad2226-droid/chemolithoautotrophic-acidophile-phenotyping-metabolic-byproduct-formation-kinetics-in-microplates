"""Baseline sulfur trajectories, wavelength comparisons and response summaries."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .plate_comparisons import CONDITION_KEYS, CITRATE_EXPERIMENT


def sulfur_baselines(tables, registry_path=None):
    path = registry_path or Path(__file__).resolve().parents[1] / "metadata/sulfur_baseline_references.json"
    requested = pd.DataFrame(json.loads(Path(path).read_text())["references"])
    keys = ["experiment_id", "group_id", "column"]
    if requested.duplicated(keys).any():
        raise ValueError("Duplicate sulfur reference formulation")
    conditions = tables["context_plate_trajectories"]
    columns = keys + ["species", "medium_id", "preparation_id", "initial_OD600_200_uL",
                      "initial_OD600_final", "component_preparation_mM_json", "component_final_mM_json"]
    conditions = conditions[columns].drop_duplicates()
    registry = requested.merge(conditions, on=keys, how="left", validate="one_to_one", indicator=True)
    if not registry._merge.eq("both").all():
        raise ValueError("A registered sulfur reference is absent from canonical conditions")
    registry = registry.drop(columns="_merge")
    source = tables["sulfur_plate_trajectories"].copy()
    source["source_table"] = "sulfur_plate_trajectories"
    source["source_row"] = np.arange(1, len(source) + 1)
    curves = source.merge(registry[keys + ["reference_id", "source_panel", "legend_label"]],
                          on=keys, how="inner", validate="many_to_one")
    if set(curves.reference_id) != set(registry.reference_id):
        raise ValueError("A registered sulfur reference has no recorded trajectories")
    support = []
    for reference, g in curves[curves.signal.eq("raw_R545_485")].groupby("reference_id", sort=False):
        cultures = sorted({s for value in g.culture_ids.dropna() for s in str(value).split(";")})
        support.append(dict(reference_id=reference, batch_ids=";".join(sorted(g.batch_id.unique())),
            physical_plate_ids=";".join(sorted(g.plate_id.unique())), physical_plate_n=g.plate_id.nunique(),
            culture_ids=";".join(cultures), culture_n=len(cultures),
            first_time_d=g.days_mean.min(), last_time_d=g.days_mean.max(),
            recorded_plate_acquisition_n=len(g), eligible_net_acquisition_n=int(g.net_change.notna().sum())))
    registry = registry.merge(pd.DataFrame(support), on="reference_id", validate="one_to_one")
    for basis in ("preparation", "final"):
        registry[basis + "_components_mM"] = registry["component_" + basis + "_mM_json"].map(
            lambda value: "; ".join(f"{key}: {dose:g}" for key, dose in json.loads(value).items()))
    registry["reference_role"] = "design_selected_optical_baseline_not_universal_zero_addition_control"
    registry["uncertainty_scope"] = "sample_SD_between_independent_microplate_means; technical_wells_averaged_within_microplate"
    return registry, curves


def timing_comparisons(tables):
    summary = tables["consolidated_activity_comparisons"].copy()
    replicates = tables["consolidated_activity_replicates"].copy()
    rows = []
    for (comparison, fmt), g in replicates.groupby(["comparison_id", "format"], sort=False):
        # A reused flask cohort remains identified; comparisons are never pooled.
        g = g.drop_duplicates("replicate_id")
        rows.append(dict(comparison_id=comparison, format=fmt,
            final_read_maximum_n=int(g.treated_maximum_at_final_read.eq(True).sum()),
            candidate_replicate_n=len(g),
            replicate_ids=";".join(sorted(g.replicate_id.astype(str)))))
    flags = pd.DataFrame(rows)
    for fmt in ("plate", "flask"):
        q = flags[flags["format"].eq(fmt)].drop(columns="format")
        q = q.rename(columns={c: fmt + "_" + c for c in q if c != "comparison_id"})
        summary = summary.merge(q, on="comparison_id", how="left", validate="one_to_one")
    summary["final_read_maximum_present"] = (summary.plate_final_read_maximum_n.fillna(0)
                                              + summary.flask_final_read_maximum_n.fillna(0)).gt(0)
    summary["timing_available"] = np.isfinite(summary.time_to_max_difference_d)
    return summary


def citrate_wavelength_pairs(tables):
    keys = CONDITION_KEYS + ["batch_id", "plate_id", "ordinal"]
    values = ["days_mean", "biotic_change", "abiotic_change", "net_change", "sources", "culture_ids"]
    raw = tables["af_iron_plate_trajectories"]
    raw = raw[raw.experiment_id.eq(CITRATE_EXPERIMENT) & raw.signal.eq("raw_A380")]
    alt = tables["channel_trajectories"]
    alt = alt[alt.experiment_id.eq(CITRATE_EXPERIMENT) & alt.signal.eq("raw_A428")]
    paired = raw[keys + values].merge(alt[keys + values], on=keys, how="outer",
                                     suffixes=("_380", "_428"), validate="one_to_one")
    meta = tables["af_iron_condition_coverage"][CONDITION_KEYS + ["component_final_mM_json"]]
    paired = paired.merge(meta, on=CONDITION_KEYS, how="left", validate="many_to_one")
    paired["citrate_mM"] = paired.component_final_mM_json.map(lambda s: json.loads(s).get("citric_acid", 0))
    paired["Ga_mM"] = paired.component_final_mM_json.map(lambda s: json.loads(s).get("gallium", 0))
    paired["pair_available"] = np.isfinite(paired.net_change_380) & np.isfinite(paired.net_change_428)
    paired["pair_status"] = np.where(paired.pair_available, "same_plate_formulation_acquisition",
                                     "missing_or_ineligible_channel_or_matched_baseline")
    return paired
