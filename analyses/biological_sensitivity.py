"""Repeatability of treatment effects, monitoring-schedule sensitivity, and citrate endpoints.

The analyses use physical vessel, matched control, incubation batch, and acquisition
identities to group the measurements.
"""
from __future__ import annotations

import itertools
import json

import numpy as np
import pandas as pd
from scipy import stats

from acidophile_data.peak_response import observed_peak
from acidophile_data.publication_si_screening_validation import (
    CONTEXT, _classify_curve, _components, _context_metadata, _id, _joined,
    _percent, _quantile, _finite,
)
from acidophile_data.publication_si_biological_results import _mean, _mean_interval, _sd
from acidophile_data.publication_si_correspondence_validation import _median


def _effect_units(tables, metadata):
    """Keep treatment and reference design identities, then average physical units."""
    lookup = {tuple(r[c] for c in CONTEXT): r for r in metadata.to_dict("records")}
    rows = []
    for row in tables["effect_plate_units"].to_dict("records"):
        if row["no_addition"]:
            # These are the selected same-context normalization references.
            # Some citrate=0 references contain Ga; self-normalization is not precision.
            continue
        meta = lookup[tuple(row[c] for c in CONTEXT)]
        refs = []
        for reference in str(row["reference_ids"]).split(";"):
            parts = reference.split("/")
            refs.append(lookup[(*parts[:3], int(parts[3]))]["condition_id"])
        references = sorted(set(refs))
        condition = _id([meta["condition_id"], references])
        definition = json.loads(meta["condition_definition_json"])
        definition["components_mM"] = sorted(definition["components_mM"])
        # Stock pH can vary between recorded preparations and is part of the series.
        series = _id([definition, references])
        rows.append(dict(format="plate", species=row["species"], stratum=meta["stratum"],
            condition_id=condition, native_condition_id=meta["condition_id"], series_id=series,
            batch_id=row["batch_id"], unit_id=row["plate_id"], physical_location=str(row["column"]),
            response_percent=row["response_percent"], delay_d=row["delay_d"],
            source_response_status=row["response_status"], source_timing_status=row["timing_status"],
            reference_condition_ids=";".join(references), reference_ids=row["reference_ids"],
            final_read_maximum=row["final_read_maximum"],
            source_context="/".join(str(row[c]) for c in CONTEXT)))
    for row in tables["flask_effect_units"].to_dict("records"):
        if row["no_addition"]:
            continue
        definition = _components(row["formulation"])
        rows.append(dict(format="flask", species="AF", stratum="All flask additions",
            condition_id="flask_" + _id(definition), native_condition_id="flask_" + _id(definition),
            series_id="flask_" + _id(sorted(definition)), batch_id=row["batch_id"],
            unit_id=row["unit_id"], physical_location=row["unit_id"],
            response_percent=row["response_percent"], delay_d=row["delay_d"],
            source_response_status=row["response_status"], source_timing_status=row["timing_status"],
            reference_condition_ids="same_batch_no_additive", reference_ids=row["reference_ids"],
            final_read_maximum=row["final_read_maximum"], source_context=row["formulation"]))
    frame = pd.DataFrame(rows)
    keys = ["format", "species", "stratum", "condition_id", "series_id", "batch_id", "unit_id"]
    duplicate = frame[frame.duplicated(keys + ["physical_location"], keep=False)]
    for key, g in duplicate.groupby(keys + ["physical_location"], sort=False):
        if not np.allclose(g[["response_percent", "delay_d"]],
                           g[["response_percent", "delay_d"]].iloc[0], rtol=0, atol=0, equal_nan=True):
            raise ValueError(f"Conflicting repeated effect membership: {key}")
    frame = frame.drop_duplicates(keys + ["physical_location"])
    return frame.groupby(keys, as_index=False).agg(
        response_percent=("response_percent", "mean"), delay_d=("delay_d", "mean"),
        physical_location_n=("physical_location", "size"),
        reference_ids=("reference_ids", _joined), reference_condition_ids=("reference_condition_ids", _joined),
        source_context=("source_context", _joined), final_read_maximum=("final_read_maximum", "max"))


def _effect_results(units):
    batch_keys = ["format", "species", "stratum", "condition_id", "series_id", "batch_id"]
    batches = units.groupby(batch_keys, as_index=False).agg(
        response_percent=("response_percent", "mean"), delay_d=("delay_d", "mean"),
        physical_unit_n=("unit_id", "size"), unit_ids=("unit_id", _joined),
        reference_ids=("reference_ids", _joined), final_read_maximum_n=("final_read_maximum", "sum"))
    conditions, rank_pairs = [], []
    for metric, field, center, unit in [("relative AUC", "response_percent", 100., "percentage points"),
                                         ("peak-time delay", "delay_d", 0., "d")]:
        for key, g in batches.groupby(batch_keys[:-1], sort=True):
            values = _finite(g[field])
            effects = values-center
            sd = float(values.std(ddof=1)) if len(values) >= 2 else np.nan
            signs = np.sign(effects)
            conditions.append(dict(zip(batch_keys[:-1], key)) | dict(metric=metric, unit=unit,
                candidate_batch_n=len(g), eligible_batch_n=len(values), effect_mean=float(effects.mean()) if len(values) else np.nan,
                between_batch_sd=sd, signed_effect_min=float(effects.min()) if len(effects) else np.nan,
                signed_effect_max=float(effects.max()) if len(effects) else np.nan,
                direction_agreement=bool(np.all(signs == signs[0])) if len(signs) >= 2 else np.nan,
                complete_batch_support=len(values) == len(g), batch_ids=_joined(g.batch_id),
                unit_ids=_joined(g.unit_ids), reference_ids=_joined(g.reference_ids)))
        # Spearman correlations use identical canonical treatment conditions in
        # each pair, not unmatched concentration grids or technical well points.
        for key, g in batches.groupby(["format", "species", "stratum", "series_id"], sort=True):
            wide = g.pivot(index="condition_id", columns="batch_id", values=field)
            for a, b in itertools.combinations(wide.columns, 2):
                pair = wide[[a, b]].dropna()
                rho = float(stats.spearmanr(pair[a], pair[b]).statistic) if (
                    len(pair) >= 3 and pair[a].nunique() > 1 and pair[b].nunique() > 1) else np.nan
                rank_pairs.append(dict(zip(["format", "species", "stratum", "series_id"], key)) |
                    dict(metric=metric, batch_a=a, batch_b=b, common_condition_n=len(pair),
                         spearman_rho=rho, status="descriptive_rank_correlation" if np.isfinite(rho)
                         else "fewer_than_three_common_conditions_or_constant_values"))
    conditions = pd.DataFrame(conditions)
    rank_pairs = pd.DataFrame(rank_pairs, columns=[
        "format", "species", "stratum", "series_id", "metric", "batch_a", "batch_b",
        "common_condition_n", "spearman_rho", "status",
    ])
    rows = []
    for key, g in conditions.groupby(["format", "species", "stratum", "metric", "unit"], sort=True):
        eligible = g[g.eligible_batch_n.ge(2)]
        complete = eligible[eligible.complete_batch_support]
        ranks = rank_pairs
        for field, value in zip(["format", "species", "stratum", "metric"], [key[0], key[1], key[2], key[3]]):
            ranks = ranks[ranks[field].eq(value)]
        series_rho = ranks.groupby("series_id").spearman_rho.mean().dropna()
        directions = eligible.direction_agreement.astype(bool)
        rows.append(dict(zip(["format", "species", "stratum", "metric", "unit"], key)) | dict(
            candidate_condition_n=len(g), repeated_condition_n=int(g.candidate_batch_n.ge(2).sum()),
            estimable_sd_condition_n=len(eligible), unavailable_sd_condition_n=len(g)-len(eligible),
            three_or_more_batch_condition_n=int(g.eligible_batch_n.ge(3).sum()),
            median_between_batch_sd=_quantile(eligible.between_batch_sd, .5),
            between_batch_sd_q25=_quantile(eligible.between_batch_sd, .25),
            between_batch_sd_q75=_quantile(eligible.between_batch_sd, .75),
            between_batch_sd_min=_quantile(eligible.between_batch_sd, 0),
            between_batch_sd_max=_quantile(eligible.between_batch_sd, 1),
            direction_agreement_n=int(directions.sum()),
            direction_agreement_percent=_percent(directions.sum(), len(directions)),
            complete_repeated_condition_n=len(complete),
            complete_direction_agreement_percent=_percent(complete.direction_agreement.astype(bool).sum(), len(complete)),
            eligible_rank_series_n=len(series_rho), eligible_rank_pair_n=int(ranks.spearman_rho.notna().sum()),
            median_spearman_rho=_quantile(series_rho, .5),
            spearman_rho_q25=_quantile(series_rho, .25), spearman_rho_q75=_quantile(series_rho, .75),
            status="descriptive_shared_microplates_and_controls_no_population_test"))
    summary = pd.DataFrame(rows)
    summary["metric_label"] = summary.metric.map({"relative AUC": "relative AUC (pp)",
                                                   "peak-time delay": "peak-time delay (d)"})
    return summary, batches, conditions, rank_pairs


def _decision(curve):
    decision = _classify_curve(curve)["decision"]
    if pd.isna(decision):
        return "unevaluable"
    return "qualifying" if decision else "no_observed_qualifying_signal"


def _scenario(curve, scenario):
    g = curve.sort_values("days_mean", kind="stable").copy().reset_index(drop=True)
    if scenario == "remove final acquisition":
        return g.iloc[:-1].copy()
    if scenario == "retain first 75% of duration":
        limit = g.days_mean.iloc[0]+.75*(g.days_mean.iloc[-1]-g.days_mean.iloc[0])
        return g[g.days_mean.le(limit)].copy()
    if scenario == "every second acquisition":
        # Preserve source NaNs even if that event would otherwise be removed.
        mask = np.arange(len(g)) % 2 == 0
        mask[-1] = True
        mask |= g.biotic_change.isna().to_numpy()
        if "abiotic_change" in g:
            mask |= g.abiotic_change.isna().to_numpy()
        return g[mask].copy()
    if scenario == "withhold middle acquisition":
        if len(g) > 2:
            fields = [x for x in ["biotic_change", "abiotic_change"] if x in g]
            g.loc[len(g)//2, fields] = np.nan
        return g
    raise ValueError(scenario)


def _acquisition_sensitivity(tables):
    """Recompute physical-unit AUC, peak timing and AF plate screening under four schedules."""
    pc = tables["af_iron_plate_trajectories"]
    pc = pc[pc.signal.eq("raw_A380")].copy()
    extra = tables["channel_trajectories"]
    extra = extra[extra.signal.isin(["raw_A428", "raw_A600"])]
    frames = pd.concat([pc, extra], ignore_index=True)
    physical_keys = ["species", "batch_id", "plate_id", "column", "signal", "ordinal"]
    duplicates = frames[frames.duplicated(physical_keys, keep=False)]
    for _, g in duplicates.groupby(physical_keys, sort=False):
        for field in ["days_mean", "biotic_change", "abiotic_change"]:
            if not np.allclose(g[field], g[field].iloc[0], equal_nan=True, rtol=1e-12, atol=1e-12):
                raise ValueError("Conflicting repeated physical-plate membership in acquisition sensitivity")
    frames = frames.drop_duplicates(physical_keys)
    keys = ["experiment_id", "group_id", "column", "batch_id", "plate_id", "signal"]
    curves = [("plate", k, g) for k, g in frames.groupby(keys, sort=True)]
    fs = tables["validation_flask_points"]
    fs = fs[fs.kind.eq("biotic")]
    for identity, g in fs.groupby("flask_ID", sort=True):
        curves.append(("flask", (g.condition_ID.iloc[0], "", 0, g.batch.iloc[0], identity, "Fe(II) loss"), g))
    scenarios = ["remove final acquisition", "retain first 75% of duration", "every second acquisition", "withhold middle acquisition"]
    rows = []
    for form, key, g in curves:
        g = g.sort_values("days_mean", kind="stable")
        base = observed_peak(g)
        decision = _decision(g) if form == "plate" and key[-1] == "raw_A380" else "not_applicable"
        for scenario in scenarios:
            simulated = _scenario(g, scenario)
            p = observed_peak(simulated)
            primary, alternate = base["auc_to_max"], p["auc_to_max"]
            eligible = np.isfinite(primary) and np.isfinite(alternate)
            simulated_decision = _decision(simulated) if decision != "not_applicable" else "not_applicable"
            rows.append(dict(format=form, readout=key[-1], experiment_id=key[0], group_id=key[1], column=key[2],
                batch_id=key[3], physical_id=key[4], scenario=scenario, source_observation_n=len(g),
                retained_observation_n=len(simulated), primary_auc=primary, alternative_auc=alternate,
                primary_status=base["status"], alternative_status=p["status"],
                auc_change_percent=100*(alternate-primary)/abs(primary) if eligible and primary != 0 else np.nan,
                timing_change_d=p["time_to_max_d"]-base["time_to_max_d"] if eligible and base["positive_rise"] and p["positive_rise"] else np.nan,
                primary_final_read_maximum=base["maximum_at_final_read"],
                alternative_final_read_maximum=p["maximum_at_final_read"],
                primary_decision=decision, alternative_decision=simulated_decision))

    details = pd.DataFrame(rows)
    output = []
    for key, g in details.groupby(["format", "readout", "scenario"], sort=False):
        primary_valid = g.primary_auc.notna()
        paired = primary_valid & g.alternative_auc.notna()
        b = g.loc[paired].groupby("batch_id")[["auc_change_percent", "timing_change_d"]].median()
        distribution = []
        for batch, batch_rows in g.loc[paired].groupby("batch_id"):
            values = _finite(batch_rows.auc_change_percent)
            timing = _finite(batch_rows.timing_change_d)
            distribution.append(dict(batch=batch,
                auc_changed_percent=100*np.mean(np.abs(values) > 1e-9) if len(values) else np.nan,
                p90_absolute_auc_change_percent=float(np.quantile(np.abs(values), .9)) if len(values) else np.nan,
                timing_changed_percent=100*np.mean(np.abs(timing) > 1e-9) if len(timing) else np.nan))
        distribution = pd.DataFrame(distribution, columns=["batch", "auc_changed_percent", "p90_absolute_auc_change_percent", "timing_changed_percent"])
        decisions = g[g.primary_decision.ne("not_applicable")]
        output.append(dict(format=key[0], readout=key[1], scenario=key[2], trajectory_n=len(g),
            readout_label=("plate " + key[1].replace("raw_", "")) if key[0] == "plate" else "flask Fe(II) loss",
            source_batch_n=g.batch_id.nunique(), physical_vessel_n=g.physical_id.nunique(),
            primary_available_n=int(primary_valid.sum()), paired_auc_n=int(paired.sum()),
            new_unavailable_n=int((primary_valid & g.alternative_auc.isna()).sum()),
            median_auc_change_percent=_median(b.auc_change_percent), median_timing_change_d=_median(b.timing_change_d),
            mean_batch_auc_changed_percent=_mean(distribution.auc_changed_percent),
            median_batch_p90_absolute_auc_change_percent=_median(distribution.p90_absolute_auc_change_percent),
            mean_batch_timing_changed_percent=_mean(distribution.timing_changed_percent),
            primary_final_read_maximum_percent=100*g.primary_final_read_maximum.mean(),
            alternative_final_read_maximum_percent=100*g.alternative_final_read_maximum.mean(),
            source_qualifying_n=int(decisions.primary_decision.eq("qualifying").sum()) if len(decisions) else np.nan,
            qualifying_to_unevaluable_n=int((decisions.primary_decision.eq("qualifying") & decisions.alternative_decision.eq("unevaluable")).sum()) if len(decisions) else np.nan,
            qualifying_to_no_event_n=int((decisions.primary_decision.eq("qualifying") & decisions.alternative_decision.eq("no_observed_qualifying_signal")).sum()) if len(decisions) else np.nan,
            qualifying_lost_n=int((decisions.primary_decision.eq("qualifying") & decisions.alternative_decision.ne("qualifying")).sum()) if len(decisions) else np.nan,
            decision_changed_n=int(decisions.primary_decision.ne(decisions.alternative_decision).sum()) if len(decisions) else np.nan,
            alternative_unevaluable_n=int(decisions.alternative_decision.eq("unevaluable").sum()) if len(decisions) else np.nan,
            source_batch_ids=_joined(g.batch_id), status="deterministic_schedule_sensitivity_not_new_experiment"))
    return pd.DataFrame(output), details


def _citrate_endpoint_effects(tables):
    """Citric-dose minus zero-citric effects on strictly paired final support.

    Endpoint eligibility is fixed before looking at measured values. A missing
    final value, mismatched final schedule or missing matched reference cannot
    cause selection of an earlier acquisition.
    """
    flask_doses = (1., 6., 12., 18.)
    experiment = "AF_AFM1_gallium_sulfate_OD0p001_citrate_ratios"
    identity = ["species", "experiment_id", "group_id", "column"]
    p = tables["af_iron_plate_trajectories"]
    p = p[p.experiment_id.eq(experiment) & p.signal.eq("raw_A380")].copy()
    meta = tables["af_iron_condition_coverage"][identity+["component_final_mM_json"]]
    p = p.merge(meta, on=identity, validate="many_to_one")
    p["citric_acid_mM"] = p.component_final_mM_json.map(lambda s: json.loads(s).get("citric_acid", 0.))
    p["Ga_mM"] = p.component_final_mM_json.map(lambda s: json.loads(s).get("gallium", 0.))
    p = p[p.Ga_mM.isin((0., 1.))]
    # The two citrate-ratio plate layouts include additional measured doses beyond
    # the four shared with the two citrate flask experiments.
    plate_doses = sorted(p.loc[p.citric_acid_mM.gt(0), "citric_acid_mM"].unique())
    normalized = []
    for readout, value in [("plate net A380 (AU)", "net_change"),
                            ("plate abiotic A380 (AU)", "abiotic_change")]:
        g = p.rename(columns={"batch_id":"batch", "plate_id":"unit_id", "ordinal":"event",
                              "days_mean":"endpoint_d", value:"endpoint_value"}).copy()
        g["match_id"] = g.unit_id
        g["readout"] = readout
        normalized.append(g[["readout", "batch", "unit_id", "match_id", "event", "endpoint_d",
                             "citric_acid_mM", "Ga_mM", "endpoint_value"]])
    f = tables["validation_flask_points"].copy()
    citrate_batches = f.loc[f.nominal_citrate_mM.gt(0), "batch"].unique()
    f = f[f.kind.eq("biotic") & f.organism.eq("AF") & f.medium.eq("AFM1") &
          f.incubation_temperature_C.eq(30) & f.batch.isin(citrate_batches)]
    for column in ["nominal_Cu_mM", "nominal_Al_mM", "nominal_V_mM", "nominal_NaCl_mM", "nominal_In_mM"]:
        f = f[f[column].eq(0)]
    f = f[f.nominal_citrate_mM.isin((0., *flask_doses)) & f.nominal_Ga_mM.isin((0., 1.))]
    f = f.rename(columns={"flask_ID":"unit_id", "inoculum_culture_id":"match_id", "sampling_event":"event",
                          "days_mean":"endpoint_d", "nominal_citrate_mM":"citric_acid_mM",
                          "nominal_Ga_mM":"Ga_mM", "net_change":"endpoint_value"})
    f["readout"] = "flask apparent Fe(III) (mM)"
    normalized.insert(0, f[["readout", "batch", "unit_id", "match_id", "event", "endpoint_d",
                            "citric_acid_mM", "Ga_mM", "endpoint_value"]])
    q = pd.concat(normalized, ignore_index=True)
    keys = ["readout", "batch", "unit_id", "citric_acid_mM", "Ga_mM"]
    if q.duplicated(keys+["event"]).any():
        raise ValueError("Duplicated physical-unit citrate acquisition")
    endpoint = q[q.event.eq(q.groupby(keys).event.transform("max"))].copy()
    match_keys = ["readout", "batch", "match_id", "Ga_mM"]
    control = endpoint[endpoint.citric_acid_mM.eq(0)].drop(columns="citric_acid_mM")
    if control.duplicated(match_keys).any():
        raise ValueError("Ambiguous same-plate or same-starter zero-citric reference")
    treated = endpoint[endpoint.citric_acid_mM.gt(0)]
    paired = treated.merge(control, on=match_keys, how="left", suffixes=("", "_reference"),
                           validate="many_to_one")
    paired["pair_status"] = np.select([
        paired.match_id.isna() | paired.unit_id_reference.isna(),
        paired.event.ne(paired.event_reference),
        ~np.isclose(paired.endpoint_d, paired.endpoint_d_reference, rtol=0, atol=1e-9),
        ~np.isfinite(paired.endpoint_value) | ~np.isfinite(paired.endpoint_value_reference)],
        ["missing_matched_zero_citric_endpoint", "final_acquisition_mismatch", "final_time_mismatch",
         "nonfinite_selected_endpoint"], default="eligible")
    paired["paired_effect"] = (paired.endpoint_value-paired.endpoint_value_reference).where(paired.pair_status.eq("eligible"))
    rows = []
    for readout in ["flask apparent Fe(III) (mM)", "plate net A380 (AU)", "plate abiotic A380 (AU)"]:
        for ga in (0., 1.):
            for dose in flask_doses if readout.startswith("flask") else plate_doses:
                candidate = paired[paired.readout.eq(readout) & paired.Ga_mM.eq(ga) & paired.citric_acid_mM.eq(dose)]
                eligible = candidate[candidate.pair_status.eq("eligible")]
                batch = eligible.groupby("batch").paired_effect.mean()
                low, high = _mean_interval(batch, minimum_n=2)
                status = ("no_eligible_paired_endpoint" if len(batch)==0 else
                          "single_batch_descriptive" if len(batch)==1 else "marginal_paired_batch_t_interval")
                # Identity and complete selection reasons stay in the machine-readable result.
                rows.append(dict(readout=readout, citric_acid_mM=dose, Ga_mM=ga,
                    paired_batch_n=len(batch), paired_unit_n=len(eligible), candidate_unit_n=len(candidate),
                    mean_effect=_mean(batch), between_batch_effect_sd=_sd(batch),
                    within_batch_pair_sd=_sd(eligible.paired_effect) if len(batch)==1 else np.nan,
                    ci_low=low, ci_high=high, endpoint_min_d=_quantile(eligible.endpoint_d, 0),
                    endpoint_max_d=_quantile(eligible.endpoint_d, 1),
                    final_event_min=_quantile(eligible.event,0), final_event_max=_quantile(eligible.event,1),
                    candidate_batch_n=candidate.batch.nunique(), unavailable_unit_n=len(candidate)-len(eligible),
                    reference_citric_acid_mM=0., batch_ids=_joined(eligible["batch"]),
                    treated_unit_ids=_joined(eligible["unit_id"]), reference_unit_ids=_joined(eligible["unit_id_reference"]),
                    pair_match_ids=_joined(eligible["match_id"]),
                    unavailable_reasons="; ".join(sorted(candidate.loc[~candidate.pair_status.eq("eligible"),"pair_status"].unique())),
                    status=status))
    return pd.DataFrame(rows), paired


def run(tables: dict[str, pd.DataFrame]) -> dict[str, pd.DataFrame]:
    """Compute three manuscript analyses from the biological workflow's numerical inputs."""
    effect_units = _effect_units(tables, _context_metadata(tables))
    repeatability, batches, conditions, ranks = _effect_results(effect_units)
    sensitivity, acquisitions = _acquisition_sensitivity(tables)
    citrate, citrate_pairs = _citrate_endpoint_effects(tables)
    return {
        "treatment_effect_repeatability": repeatability,
        "treatment_effect_physical_units": effect_units,
        "treatment_effect_batches": batches,
        "treatment_effect_condition_results": conditions,
        "treatment_rank_correlations": ranks,
        "acquisition_sensitivity": sensitivity,
        "acquisition_sensitivity_trajectories": acquisitions,
        "citrate_endpoint_effects": citrate,
        "citrate_endpoint_pairs": citrate_pairs,
    }
