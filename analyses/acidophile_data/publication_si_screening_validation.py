"""Si screening validation."""
from __future__ import annotations
import hashlib
import json
import numpy as np
import pandas as pd
from .publication_screening_support import add_decision_support, add_agreement_support


CONTEXT = ["species", "experiment_id", "group_id", "column"]


PHYSICAL = ["species", "batch_id", "plate_id", "column"]


def _id(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                    allow_nan=False).encode()).hexdigest()[:20]


def _finite(values):
    x = np.asarray(values, float)
    return x[np.isfinite(x)]


def _quantile(values, p):
    x = _finite(values)
    return float(np.quantile(x, p)) if len(x) else np.nan


def _percent(numerator, denominator):
    return 100. * numerator / denominator if denominator else np.nan


def _joined(values):
    return ";".join(sorted({item for value in values if pd.notna(value)
                            for item in str(value).split(";") if item}))


def _components(value):
    return {k: float(v) for k, v in json.loads(value).items() if float(v) != 0.}


def _result(key, title, caption, frame, display_columns, sources, panels, claims,
            result_columns, support_tables=None):
    if frame.empty or any(c not in frame for c in display_columns + result_columns):
        raise ValueError(f"Invalid screening validation result: {key}")
    labels = {"candidate_plate_condition_n": "plate-conditions (n)",
              "unevaluable_decision_n": "unevaluable (n)", "event_n": "events (n)",
              "event_percent": "events (%)", "abiotic_threshold_exceeded_n": "Ab exceedance (n)",
              "first_event_net_min_AU": "minimum net (AU)",
              "first_event_net_median_AU": "median net (AU)",
              "event_without_net_1AU_n": "weak-net events (n)",
              "comparison_level": "comparison", "candidate_comparison_n": "candidate (n)",
              "evaluable_comparison_n": "evaluable (n)",
              "identical_event_fraction_percent": "identical rates (%)",
              "event_fraction_range_median_pp": "median range (pp)",
              "event_fraction_range_q75_pp": "Q3 range (pp)",
              "candidate_condition_n": "candidate (n)",
              "estimable_sd_condition_n": "SD estimable (n)",
              "three_or_more_batch_condition_n": "at least 3 batches (n)",
              "median_between_batch_sd": "median SD", "between_batch_sd_q25": "Q1 SD",
              "between_batch_sd_q75": "Q3 SD", "direction_agreement_percent": "same direction (%)",
              "eligible_rank_series_n": "rank series (n)", "median_spearman_rho": "median Spearman rho",
              "metric_label": "metric (SD units)"}
    return dict(key=key, title=title, caption=caption, frame=frame,
                display_columns=display_columns, source_tables=sources,
                figure_panels=panels, claims=claims, result_columns=result_columns,
                column_labels={c: labels[c] for c in display_columns if c in labels},
                support_tables=support_tables or {})


def _context_metadata(tables):
    """Identify AF experimental conditions by their chemistry and preparation."""
    units = tables["condition_repeatability_units"]
    units = units[units["format"].eq("plate") & units.metric.eq("AUC_to_own_max")]
    identities = {}
    for row in units.itertuples():
        for context in str(row.source_context_ids).split(";"):
            parts = context.split("/")
            key = (*parts[:3], int(parts[3]))
            value = (row.condition_id, row.condition_definition_json)
            if key in identities and identities[key] != value:
                raise ValueError(f"Conflicting canonical condition identity: {key}")
            identities[key] = value
    records = []
    for row in tables["af_iron_condition_coverage"].to_dict("records"):
        key = tuple(row[c] for c in CONTEXT)
        condition_id, definition = identities[key]
        components = _components(row["component_final_mM_json"])
        records.append({c: row[c] for c in CONTEXT} | dict(
            condition_id=condition_id, condition_definition_json=definition,
            stratum=row["additive_class"] if components else "No added compound",
            added_compound=bool(components), component_final_mM_json=row["component_final_mM_json"],
            medium_id=row["medium_id"]))
    lf = tables["validation_lf_plate_trajectories"]
    fields = CONTEXT + ["component_final_mM_json", "component_basis_json", "medium_id",
                        "preparation_id", "readout", "initial_OD600_final"]
    lf = lf[fields].drop_duplicates()
    if lf.duplicated(CONTEXT).any():
        raise ValueError("Conflicting LF condition metadata")
    for row in lf.to_dict("records"):
        components = _components(row["component_final_mM_json"])
        definition = {field: row[field] for field in fields if field not in CONTEXT}
        definition["components_mM"] = components
        del definition["component_final_mM_json"]
        # A zero concentration does not create a distinct no-addition chemistry.
        basis = json.loads(definition.pop("component_basis_json"))
        definition["concentration_bases"] = {k: basis.get(k, "unspecified") for k in components}
        definition["species"] = row["species"]
        records.append({c: row[c] for c in CONTEXT} | dict(
            condition_id="lf_" + _id(definition),
            condition_definition_json=json.dumps(definition, sort_keys=True),
            stratum="LF additions" if components else "No added compound",
            added_compound=bool(components), component_final_mM_json=row["component_final_mM_json"],
            medium_id=row["medium_id"]))
    metadata = pd.DataFrame(records)
    if metadata.duplicated(CONTEXT).any():
        raise ValueError("Duplicate condition metadata")
    return metadata


def _classify_curve(curve):
    """A non-event is established only on a complete recorded postbaseline set."""
    q = curve.sort_values(["days_mean", "ordinal"])
    initial = q[q.days_mean.eq(0)]
    initial_valid = bool(len(initial) == 1 and
                         np.isfinite(initial[["biotic_change", "abiotic_change"]].to_numpy(float)).all())
    # The recorded initial zero cannot by itself establish an absent event.
    q = q[q.days_mean.gt(0)]
    biotic = np.isfinite(q.biotic_change.to_numpy(float))
    abiotic = np.isfinite(q.abiotic_change.to_numpy(float))
    matched = biotic & abiotic
    event = matched & q.biotic_change.ge(1).to_numpy() & q.abiotic_change.lt(1).to_numpy()
    complete = bool(initial_valid and len(q) and matched.all())
    outcome = 1. if event.any() else 0. if complete else np.nan
    passing = q.loc[event]
    net = passing.biotic_change - passing.abiotic_change
    status = ("recorded_event" if event.any() else "complete_recorded_window_no_event" if complete
              else "no_postbaseline_acquisition" if not len(q)
              else "partial_recorded_window_no_event" if matched.any()
              else "no_matched_postbaseline_acquisition")
    return dict(decision=outcome, decision_status=status, event=bool(event.any()),
        complete_recorded_window=complete, initial_pair_valid=initial_valid, postbaseline_acquisition_n=len(q),
        eligible_matched_acquisition_n=int(matched.sum()),
        unpaired_biotic_threshold_acquisition_n=int((biotic & ~abiotic & q.biotic_change.ge(1)).sum()),
        abiotic_threshold_exceeded=bool((abiotic & q.abiotic_change.ge(1)).any()),
        first_event_net_AU=float(net.iloc[0]) if len(net) else np.nan,
        maximum_event_net_AU=float(net.max()) if len(net) else np.nan,
        event_without_net_1AU=bool(len(net) and net.max() < 1),
        first_event_time_d=float(passing.days_mean.iloc[0]) if len(passing) else np.nan,
        recorded_end_d=float(q.days_mean.max()) if len(q) else np.nan)


def _screen_units(tables, metadata):
    source = pd.concat([tables["af_iron_plate_trajectories"],
                        tables["validation_lf_plate_trajectories"]], ignore_index=True)
    source = source[source.signal.eq("raw_A380")]
    fields = CONTEXT + ["batch_id", "plate_id", "ordinal", "days_mean", "biotic_change",
                        "abiotic_change", "net_change", "culture_ids"]
    source = source[fields].merge(metadata, on=CONTEXT, validate="many_to_one")
    # Repeated metadata memberships must agree exactly before deduplication.
    duplicate = source[source.duplicated(PHYSICAL + ["ordinal"], keep=False)]
    for key, g in duplicate.groupby(PHYSICAL + ["ordinal"], sort=False):
        for field in ("days_mean", "biotic_change", "abiotic_change", "net_change"):
            values = g[field].to_numpy(float)
            if not np.allclose(values, values[0], rtol=0, atol=0, equal_nan=True):
                raise ValueError(f"Conflicting physical screening observation: {key}, {field}")
        if g.condition_id.nunique() != 1 or g.stratum.nunique() != 1:
            raise ValueError(f"Conflicting duplicate physical condition: {key}")
    source = source.drop_duplicates(PHYSICAL + ["ordinal"])
    rows = []
    for key, g in source.groupby(PHYSICAL, sort=True):
        first = g.iloc[0]
        rows.append(dict(zip(PHYSICAL, key)) | {
            field: first[field] for field in ("condition_id", "condition_definition_json", "stratum",
                                              "added_compound", "medium_id")}
            | dict(culture_ids=_joined(g.culture_ids), recorded_acquisition_n=len(g))
            | _classify_curve(g))
    output = pd.DataFrame(rows)
    if output.duplicated(["species", "condition_id", "batch_id", "plate_id"]).any():
        raise ValueError("Repeated equal-chemistry columns require physical-plate aggregation")
    # Event calls must agree between the trajectory and plate-level summaries.
    plate_events = tables["screen_plate_evidence"]
    joined = plate_events.merge(output[PHYSICAL + ["event"]], on=PHYSICAL, validate="many_to_one")
    if len(joined) != len(plate_events) or not joined.optical_activity_event.eq(joined.event).all():
        raise ValueError("Trajectory and plate-level event calls differ")
    return output


def _decision_performance(units):
    rows = []
    for (species, stratum), g in units.groupby(["species", "stratum"], sort=True):
        event = g[g.event]
        usable = int(g.decision.notna().sum())
        rows.append(dict(species=species, stratum=stratum, candidate_plate_condition_n=len(g),
            physical_plate_n=g.plate_id.nunique(), starter_batch_n=g.batch_id.nunique(),
            evaluable_decision_n=usable, unevaluable_decision_n=len(g)-usable,
            partial_no_event_n=int(g.decision_status.eq("partial_recorded_window_no_event").sum()),
            event_n=len(event), event_percent=_percent(len(event), usable),
            no_event_n=int(g.decision.eq(0).sum()),
            abiotic_threshold_exceeded_n=int(g.abiotic_threshold_exceeded.sum()),
            abiotic_threshold_exceeded_percent=_percent(g.abiotic_threshold_exceeded.sum(), len(g)),
            first_event_net_min_AU=_quantile(event.first_event_net_AU, 0),
            first_event_net_median_AU=_quantile(event.first_event_net_AU, .5),
            first_event_net_q25_AU=_quantile(event.first_event_net_AU, .25),
            first_event_net_q75_AU=_quantile(event.first_event_net_AU, .75),
            event_without_net_1AU_n=int(event.event_without_net_1AU.sum()),
            event_without_net_1AU_percent=_percent(event.event_without_net_1AU.sum(), len(event)),
            postbaseline_acquisition_n=int(g.postbaseline_acquisition_n.sum()),
            eligible_matched_acquisition_n=int(g.eligible_matched_acquisition_n.sum()),
            unpaired_biotic_threshold_acquisition_n=int(g.unpaired_biotic_threshold_acquisition_n.sum()),
            observed_window_median_d=_quantile(g.recorded_end_d, .5),
            observed_window_min_d=_quantile(g.recorded_end_d, 0),
            observed_window_max_d=_quantile(g.recorded_end_d, 1),
            batch_ids=_joined(g.batch_id), status="descriptive_recorded_screen_not_diagnostic_validation"))
    return pd.DataFrame(rows)


def _agreement_records(units):
    """Assess independent-microplate spread separately from batch event-rate spread."""
    batch_rows, within = [], []
    keys = ["species", "stratum", "condition_id", "batch_id"]
    for key, g in units.groupby(keys, sort=True):
        r = dict(zip(keys, key))
        values = _finite(g.decision)
        complete = len(values) == len(g)
        r.update(candidate_unit_n=len(g), evaluable_unit_n=len(values),
                 complete=complete, event_fraction=float(values.mean()) if complete else np.nan,
                 fraction_range=float(values.max()-values.min()) if complete else np.nan,
                 physical_plate_ids=_joined(g.plate_id), culture_ids=_joined(g.culture_ids))
        batch_rows.append(r)
        if len(g) >= 2:
            within.append(r | dict(level="plates within batch", candidate=True))
    batches = pd.DataFrame(batch_rows)
    between = []
    for key, g in batches.groupby(keys[:-1], sort=True):
        if len(g) < 2:
            continue
        complete = bool(g.complete.all())
        values = _finite(g.event_fraction)
        between.append(dict(zip(keys[:-1], key)) | dict(level="between batches", candidate=True,
            candidate_unit_n=len(g), evaluable_unit_n=len(values), complete=complete,
            event_fraction=float(values.mean()) if complete else np.nan,
            fraction_range=float(values.max()-values.min()) if complete else np.nan,
            mixed_batch_n=int((g.event_fraction.gt(0) & g.event_fraction.lt(1)).sum()),
            batch_ids=_joined(g.batch_id), physical_plate_ids=_joined(g.physical_plate_ids),
            culture_ids=_joined(g.culture_ids)))
    records = pd.DataFrame(within + between)
    rows = []
    # Include unavailable levels explicitly, including LF between-batch support.
    for (species, stratum), _ in units.groupby(["species", "stratum"], sort=True):
        for level in ("plates within batch", "between batches"):
            g = records[records.species.eq(species) & records.stratum.eq(stratum) & records.level.eq(level)]
            eligible = g[g.complete]
            ranges = _finite(eligible.fraction_range)
            rows.append(dict(species=species, stratum=stratum, comparison_level=level,
                candidate_comparison_n=len(g), evaluable_comparison_n=len(eligible),
                incomplete_comparison_n=len(g)-len(eligible),
                identical_event_fraction_n=int((ranges == 0).sum()),
                identical_event_fraction_percent=_percent((ranges == 0).sum(), len(ranges)),
                disagreement_percent=_percent((ranges > 0).sum(), len(ranges)),
                event_fraction_range_median_pp=100*_quantile(ranges, .5),
                event_fraction_range_q25_pp=100*_quantile(ranges, .25),
                event_fraction_range_q75_pp=100*_quantile(ranges, .75),
                event_fraction_range_max_pp=100*_quantile(ranges, 1),
                status="descriptive_consistency" if len(eligible) else "no_evaluable_repeated_groups"))
    return pd.DataFrame(rows), batches, records


def build(tables):
    metadata = _context_metadata(tables)
    screen_units = _screen_units(tables, metadata)
    decisions = _decision_performance(screen_units)
    agreement, screen_batches, agreement_records = _agreement_records(screen_units)
    screen_sources = ["af_iron_plate_trajectories", "validation_lf_plate_trajectories",
                      "screen_plate_evidence", "af_iron_condition_coverage",
                      "condition_repeatability_units"]
    results = [
        _result("screen_decision_performance", "Recorded screen decisions and biotic–abiotic separation",
            "The prespecified operational event is a physical-plate mean biotic "
            "$\\Delta A_{380}\\geq1$ AU with matched abiotic $\\Delta A_{380}<1$ AU at the same "
            "recorded postbaseline acquisition. Technical wells are averaged within physical plate "
            "using the complete-assigned-well rule. Repeated scientific memberships of the "
            "same plate and column are counted once. An event is evaluable when observed; absence "
            "is evaluable only with a finite initial pair and when every recorded postbaseline acquisition has both finite signals. "
            "Other no-event trajectories remain unevaluable, not negative. Event percentage uses "
            "evaluable decisions. Abiotic exceedance counts any recorded abiotic change at least "
            "1 AU. Net separation is biotic minus abiotic change at the same acquisition; the "
            "minimum and median refer to the first event per passing plate-condition. Weak-net "
            "events are passing plate-conditions with no passing acquisition reaching 1 AU net "
            "separation. No-added-compound rows retain all recorded nominal preparation/pH contexts. "
            "Other rows group the tested additions by chemical class, not by biological mechanism. "
            "Percentages describe this observed, duration-dependent screen; shared microplates and "
            "controls prevent independent binomial inference. These are not limits of detection, "
            "diagnostic sensitivity/specificity, or viable-cell classifications.",
            decisions, ["species", "stratum", "candidate_plate_condition_n", "unevaluable_decision_n",
                        "event_n", "event_percent", "abiotic_threshold_exceeded_n",
                        "first_event_net_min_AU", "first_event_net_median_AU", "event_without_net_1AU_n"],
            screen_sources, ["5a", "7g"],
            ["The operational screen criterion is distinct from a 1 AU net biotic–abiotic separation.",
             "Missing-control cases cannot establish an absent abiotic contribution or a negative result."],
            ["event_percent", "first_event_net_min_AU", "first_event_net_median_AU", "event_without_net_1AU_percent"],
            {"screen_operational_decisions": screen_units}),
        _result("screen_decision_repeatability", "Consistency of recorded screen classifications",
            "Classification is defined in the screen-decision table. Conditions retain chemical "
            "composition, concentration basis, medium, inoculum and recorded preparation/pH "
            "attributes. Within a starting batch, each comparison contains at least two distinct "
            "physical plates at the same condition; its event-fraction range is 0 or 100 percentage "
            "points. Between batches, each comparison contains the same condition in at least two "
            "batches; its range is the maximum minus minimum of the batch-specific fractions of "
            "passing independent microplates. Batches receive equal weight and differing microplate outcomes "
            "are retained as fractional rates. A comparison is evaluable only if every constituent "
            "plate decision is evaluable. Identical fractions do not establish equivalence or "
            "between-date reproducibility. Quantiles summarize condition-level ranges, "
            "without treating shared microplates, controls or formulations as independent replicates. "
            "LF has one recorded starting batch; its between-batch consistency is unavailable. "
            "Differences in observation windows remain part of the observed consistency.",
            agreement, ["species", "stratum", "comparison_level", "candidate_comparison_n",
                        "evaluable_comparison_n", "identical_event_fraction_percent",
                        "event_fraction_range_median_pp", "event_fraction_range_q75_pp"],
            screen_sources, ["5a", "7g"],
            ["Within-date agreement among independent microplates and between-date consistency are distinct validation quantities.",
             "The LF screen cannot estimate between-starting-batch repeatability."],
            ["identical_event_fraction_percent", "event_fraction_range_median_pp", "event_fraction_range_q75_pp"],
            {"screen_batch_decisions": screen_batches, "screen_agreement_comparisons": agreement_records})]
    results[0]['frame'], support = add_decision_support(tables, metadata, screen_units, decisions)
    results[0]['support_tables'].update(support)
    results[0]['caption'] += (
        ' Additional support reports event, complete no-event and unevaluable prevalence with explicit '
        'candidate and evaluable denominators. At each acquired pair the biotic gate margin is '
        'biotic change minus 1 AU (passing at or above zero), and the abiotic margin is 1 AU minus '
        'abiotic change (passing strictly above zero). First-event margins use the first qualifying '
        'acquisition. For complete no-event trajectories, nearest joint approach maximizes the '
        'smaller of the two margins at the same acquisition, breaking ties by earliest ordinal. '
        'Margins from different acquisitions are never combined; missing pairs remain unavailable.')
    results[1]['frame'], support = add_agreement_support(agreement, agreement_records)
    results[1]['support_tables'].update(support)
    results[1]['caption'] += (
        ' Additional counts partition candidate comparisons into all-positive, all-negative, '
        'identical mixed fractions, discordant and unevaluable outcomes. Identical mixed fractions '
        'apply only between batches. The first three categories sum to the identical '
        'fraction total; unevaluable comparisons are excluded from its percentage denominator.')
    return results
