"""Si correspondence validation."""
from __future__ import annotations
import json
import numpy as np
import pandas as pd
from .peak_response import observed_peak


def _finite(values):
    x = np.asarray(values, float)
    return x[np.isfinite(x)]


def _mean(values):
    x = _finite(values)
    return float(x.mean()) if len(x) else np.nan


def _median(values):
    x = _finite(values)
    return float(np.median(x)) if len(x) else np.nan


def _ids(values):
    return ";".join(sorted(set(str(x) for x in values if pd.notna(x))))


def _result(key, title, caption, frame, columns, sources, panels, claims, support=None):
    labels = {"plate_dose": "Plate dose (mM)", "flask_dose": "Flask dose (mM)",
        "dose_match": "Dose relation", "match_status": "Dose relation",
        "paired_context_n": "Paired contexts (n)", "unavailable_context_n": "Unavailable (n)",
        "raw_error": "A380 error", "corrected_error": "A*380 error",
        "paired_error_difference": "Corrected − raw error", "readout_label": "Readout",
        "paired_auc_n": "Paired AUC (n)", "new_unavailable_n": "New unavailable (n)",
        "mean_batch_auc_changed_percent": "Changed AUC (%)",
        "median_batch_p90_absolute_auc_change_percent": "90th percentile |ΔAUC| (%)",
        "median_timing_change_d": "Median Δtime (d)", "qualifying_lost_n": "Qualifying lost (n)",
        "A380_estimate": "A380 estimate", "A600_estimate": "A600 estimate",
        "paired_difference": "A600 − A380", "ci_low": "95% CI lower", "ci_high": "95% CI upper",
        "plate_dose_mM": "Plate dose (mM)", "flask_dose_mM": "Flask dose (mM)",
        "Ga_mM": "Ga (mM)", "plate_response_change_pp": "Plate Δresponse (pp)",
        "flask_response_change_pp": "Flask Δresponse (pp)",
        "format_difference_change_pp": "Δ(format difference) (pp)"}
    return dict(key=key, title=title, caption=caption, frame=frame,
        display_columns=columns, source_tables=sources, figure_panels=panels,
        column_labels={c: labels[c] for c in columns if c in labels},
        claims=claims, result_columns=[c for c in frame if any(s in c for s in
            ("error", "difference", "change", "sensitivity", "ci_", "unavailable", "loss", "qualifying", "median_", "estimate"))],
        support_tables=support or {})


def _dose(value):
    q = json.loads(value)
    return "; ".join(f"{k.replace('citric_acid', 'citric acid')} {v:g}" for k, v in sorted(q.items())) or "no additive"


def _strict_window(curve, start, end, value):
    """Observed knots and necessary boundary brackets, never across a missing read."""
    g = curve.sort_values("days_mean", kind="stable")
    x, y = g.days_mean.to_numpy(float), g[value].to_numpy(float)
    if not np.isfinite([start, end]).all() or len(x) < 2 or not np.isfinite(x).all() or np.any(np.diff(x) <= 0):
        return None
    if start < x[0] or end > x[-1] or end <= start or not np.isfinite(y[0]):
        return None
    lo = max(0, np.searchsorted(x, start, side="right")-1)
    hi = min(len(x)-1, np.searchsorted(x, end, side="left"))
    if not np.isfinite(y[lo:hi+1]).all():
        return None
    knots = np.unique(np.r_[start, x[(x >= start) & (x <= end)], end])
    return knots, np.interp(knots, x, y)


def _shape_error(p, f, start, end):
    from .publication_cross_format import temporal_shape_error
    a, b = _strict_window(p, start, end, "net_change"), _strict_window(f, start, end, "mean")
    if a is None or b is None:
        return dict(status="missing_recorded_common_window_support", time_weighted_mae=np.nan,
                    time_weighted_rmse=np.nan)
    result = temporal_shape_error(*a, *b)
    if result["status"] == "descriptive_shape_comparison":
        grid = np.unique(np.r_[a[0], b[0]])
        pv, fv = np.interp(grid, *a), np.interp(grid, *b)
        delta = (pv-pv[0])/(pv[-1]-pv[0])-(fv-fv[0])/(fv[-1]-fv[0])
        crosses = np.flatnonzero(delta[:-1]*delta[1:] < 0)
        roots = grid[crosses]-delta[crosses]*np.diff(grid)[crosses]/np.diff(delta)[crosses]
        absolute_grid = np.unique(np.r_[grid, roots])
        result["time_weighted_mae"] = float(np.trapezoid(np.abs(np.interp(absolute_grid, grid, delta)), absolute_grid)/(end-start))
        result["time_weighted_rmse"] = float(np.sqrt(np.sum(np.diff(grid)*(delta[:-1]**2+delta[:-1]*delta[1:]+delta[1:]**2)/3)/(end-start)))
    return result


def _common_readout(tables):
    correspondence = tables["cross_format_condition_correspondence"]
    trajectories = tables["cross_format_plate_timecourses"]
    flask_source = tables["cross_format_flask_timecourses"]
    keys = ["experiment_id", "group_id", "column", "batch_id", "signal", "ordinal"]
    means = trajectories.groupby(keys, as_index=False).agg(
        days_mean=("days_mean", "mean"), net_change=("net_change", "mean"),
        valid_n=("net_change", "count"), culture_n=("culture_id", "nunique"),
        physical_ids=("physical_plate_ids", _ids), culture_ids=("culture_id", _ids))
    curve_keys = keys[:-1]
    means["expected_n"] = means.groupby(curve_keys).culture_n.transform("max")
    means.loc[means.valid_n.ne(means.expected_n), "net_change"] = np.nan
    curves = dict(tuple(means.groupby(curve_keys)))
    fs = flask_source[flask_source.metric.eq("feii_depletion_fraction")].copy()
    fs["days_mean"] = fs.elapsed_time_mean_d
    fs.loc[fs.n.ne(fs.n_vessels_in_event), "mean"] = np.nan
    flasks = dict(tuple(fs.groupby("condition_ID")))
    rows = []
    selected = correspondence[correspondence.match_status.isin([
        "equal_reported_dose", "nearest_lower_dose_context", "nearest_upper_dose_context"])]
    for _, m in selected.iterrows():
        f = flasks.get(m.flask_condition_ID)
        if f is None:
            continue
        prefix = (m.experiment_id, m.group_id, m.column)
        for key, raw in curves.items():
            if key[:3] != prefix or key[-1] != "raw_A380":
                continue
            corrected = curves.get((*key[:-1], "A380_minus_A750"))
            if corrected is None:
                continue
            # The same originally assigned plate acquisitions define both readouts.
            r = raw.merge(corrected[["ordinal", "net_change"]], on="ordinal", suffixes=("", "_corrected"), validate="one_to_one")
            start = max(r.days_mean.min(), f.days_mean.min())
            end = min(r.days_mean.max(), f.days_mean.max())
            rr = _shape_error(r, f, start, end)
            cc = _shape_error(r.assign(net_change=r.net_change_corrected), f, start, end)
            valid = rr["status"] == cc["status"] == "descriptive_shape_comparison"
            rows.append(dict(correspondence_id=m.correspondence_id, flask_condition_ID=m.flask_condition_ID,
                plate_batch=key[3], flask_batch=m.flask_batch, plate_physical_ids=_ids(raw.physical_ids),
                plate_culture_ids=_ids(raw.culture_ids), flask_culture_ids=m.flask_culture_ids,
                plate_dose=_dose(m.comparable_quantities_mM_json), flask_dose=_dose(m.flask_quantities_mM_json),
                match_status=m.match_status, window_start_d=start, window_end_d=end,
                raw_error=rr["time_weighted_mae"], corrected_error=cc["time_weighted_mae"],
                paired_error_difference=cc["time_weighted_mae"]-rr["time_weighted_mae"] if valid else np.nan,
                raw_rmse=rr["time_weighted_rmse"], corrected_rmse=cc["time_weighted_rmse"],
                status="paired_common_window_estimated" if valid else rr["status"]+";"+cc["status"]))
    details = pd.DataFrame(rows).drop_duplicates([
        "flask_condition_ID", "plate_batch", "plate_physical_ids", "plate_dose", "flask_dose"])
    output = []
    for key, g in details.groupby(["plate_dose", "flask_dose", "match_status"], sort=True):
        good = g[g.paired_error_difference.notna()]
        # Equal weight for source-cohort combinations, not technical conditions.
        b = good.groupby(["flask_condition_ID", "plate_batch"])[["raw_error", "corrected_error", "paired_error_difference"]].mean()
        output.append(dict(plate_dose=key[0], flask_dose=key[1], dose_match=key[2].replace("_", " "),
            candidate_context_n=len(g), paired_context_n=len(good), cohort_pair_n=len(b),
            plate_batch_n=g.plate_batch.nunique(), flask_batch_n=g.flask_batch.nunique(),
            raw_error=_mean(b.raw_error), corrected_error=_mean(b.corrected_error),
            paired_error_difference=_mean(b.paired_error_difference),
            error_improved_context_percent=100*good.paired_error_difference.lt(0).mean() if len(good) else np.nan,
            unavailable_context_n=len(g)-len(good), window_min_d=g.window_end_d.min(), window_max_d=g.window_end_d.max(),
            correspondence_ids=_ids(g.correspondence_id), plate_batch_ids=_ids(g.plate_batch),
            flask_batch_ids=_ids(g.flask_batch), status="descriptive_reused_cohorts_no_equivalence_test"))
    return _result("common_window_readout_error", "Common-window optical–chemical trajectory discrepancies",
        "Raw and A750-corrected net A380 changes are compared with fractional flask Fe(II) depletion over "
        "the same observed overlap. Each trajectory is shifted to its overlap start and divided by its positive "
        "end-minus-start change. Error is the time-weighted mean absolute discrepancy of piecewise-linear "
        "normalized trajectories, integrating exactly across zero crossings; a difference below zero favors "
        "the corrected readout. The support export also retains exact piecewise-linear RMS discrepancies. This normalization "
        "assesses shape, not amplitude or concentration accuracy. Every assigned culture and recorded "
        "measurement needed in the overlap must be valid; missing initial/internal values are not bridged. "
        "Only contexts estimable for both readouts contribute paired means. Technical contexts are averaged "
        "within each flask-condition/plate-batch combination before equally averaging those combinations. "
        "Cohorts are reused across comparisons and no population CI, P value or format-equivalence claim is "
        "assigned. Doses are mM; neighboring doses remain context comparisons, not concentration matches. "
        "Shared cohort identities and unavailable results are retained in the context-level output.",
        pd.DataFrame(output), ["plate_dose", "flask_dose", "dose_match", "paired_context_n", "unavailable_context_n",
            "raw_error", "corrected_error", "paired_error_difference"],
        ["cross_format_condition_correspondence", "cross_format_plate_timecourses", "cross_format_flask_timecourses"],
        ["5a", "5b", "7a", "7b", "7c", "7d", "7e", "7f"],
        ["Optical–chemical correspondence depends on the readout, chemical context and observed overlap; amplitude agreement is not tested."],
        {"context_results": details})


def _window_response(tables):
    pc = tables["af_iron_plate_trajectories"]
    pc = pc[pc.signal.eq("raw_A380")]
    pkeys = ["species", "experiment_id", "group_id", "column", "batch_id", "plate_id"]
    plate = dict(tuple(pc.groupby(pkeys)))
    fs = tables["validation_flask_points"]
    flask = dict(tuple(fs[fs.kind.eq("biotic")].groupby("flask_ID")))
    primary = tables["consolidated_activity_replicates"]
    rows = []
    for comparison, group in primary.groupby("comparison_id", sort=True):
        selected, all_curves = [], []
        for _, r in group.iterrows():
            if r["format"] == "plate":
                key = (r.species, r.experiment_id, r.group_id, r.column, r.plate_batch, r.replicate_id)
                g = plate[key]
                refs = [plate.get((*key[:3], int(float(c)), *key[4:])) for c in str(r.reference_ids).split(";") if c and c != "nan"]
            else:
                g = flask[r.replicate_id]
                refs = [flask.get(x) for x in str(r.reference_ids).split(";") if x and x != "nan"]
            if not refs:
                refs = [None]  # A shared window requires all assigned control trajectories.
            # The established primary uses biotic A380 change and flask Fe(II) loss.
            check = observed_peak(g)
            if not np.isclose(check["auc_to_max"], r.treated_auc_to_max, equal_nan=True, rtol=1e-9, atol=1e-10):
                raise ValueError(f"Primary own-maximum AUC not reproduced: {comparison}/{r.replicate_id}")
            selected.append((r, g, refs))
            all_curves.extend([g, *refs])
        complete = all(c is not None for c in all_curves)
        start = max(c.days_mean.min() for c in all_curves) if complete else np.nan
        end = min(c.days_mean.max() for c in all_curves) if complete else np.nan
        for r, g, refs in selected:
            spans = [_strict_window(c, start, end, "biotic_change") if complete else None for c in [g, *refs]]
            areas = [float(np.trapezoid(y, x)) if span is not None else np.nan
                     for span in spans for x, y in ([span] if span is not None else [(None, None)])]
            fixed = 100*areas[0]/np.mean(areas[1:]) if np.isfinite(areas).all() and np.mean(areas[1:]) > 0 and pd.notna(r.activity_percent) else np.nan
            rows.append(dict(comparison_id=comparison, format=r["format"], additive=r.additive,
                plate_dose_mM=r.plate_dose_mM, flask_dose_mM=r.flask_dose_mM, Ga_mM=r.Ga_mM,
                match_status=r.match_status, plate_batch=r.plate_batch, flask_batch=r.flask_batch,
                flask_condition_ID=r.flask_condition_ID, physical_id=r.replicate_id,
                reference_ids=r.reference_ids, window_start_d=start, window_end_d=end,
                primary_response_percent=r.activity_percent, fixed_response_percent=fixed,
                response_change_pp=fixed-r.activity_percent, status="paired_estimated" if np.isfinite(fixed) else "missing_common_window_or_primary_eligibility"))
    details = pd.DataFrame(rows)
    output = []
    keys = ["additive", "plate_dose_mM", "flask_dose_mM", "Ga_mM", "match_status"]
    for key, g in details.groupby(keys, dropna=False, sort=True):
        reduced = []
        for (cid, form), vessels in g.groupby(["comparison_id", "format"]):
            candidates = vessels[vessels.primary_response_percent.notna()]
            complete = len(candidates) > 0 and candidates.fixed_response_percent.notna().all()
            reduced.append(dict(comparison_id=cid, format=form, **{field: _mean(candidates[field]) if complete else np.nan
                for field in ["primary_response_percent", "fixed_response_percent", "response_change_pp"]}))
        c = pd.DataFrame(reduced).set_index(["comparison_id", "format"])
        wide = c.reset_index().pivot(index="comparison_id", columns="format", values=["primary_response_percent", "fixed_response_percent", "response_change_pp"])
        paired = wide.dropna()
        change = paired["response_change_pp"]["plate"]-paired["response_change_pp"]["flask"]
        output.append(dict(zip(keys, key), candidate_context_n=g.comparison_id.nunique(), paired_context_n=len(paired),
            plate_response_change_pp=_mean(paired["response_change_pp"]["plate"]),
            flask_response_change_pp=_mean(paired["response_change_pp"]["flask"]),
            format_difference_change_pp=_mean(change),
            primary_format_difference_pp=_mean(paired["primary_response_percent"]["plate"]-paired["primary_response_percent"]["flask"]),
            fixed_format_difference_pp=_mean(paired["fixed_response_percent"]["plate"]-paired["fixed_response_percent"]["flask"]),
            window_min_d=g.window_end_d.min(), window_max_d=g.window_end_d.max(),
            plate_batch_n=g.plate_batch.nunique(), flask_batch_n=g.flask_batch.nunique(),
            unavailable_context_n=g.comparison_id.nunique()-len(paired), comparison_ids=_ids(g.comparison_id),
            status="descriptive_shared_controls_and_cohorts_no_population_test"))
    return _result("common_window_response_sensitivity", "Control-relative response sensitivity to a common observation window",
        "The primary response remains 100×treated AUC/mean separately integrated no-additive control AUC, "
        "with every trajectory ending at its own first observed maximum. The sensitivity calculation instead "
        "integrates all treated and reference trajectories to a single common end within each plate/flask "
        "comparison, retaining postmaximum observations within that window. This deliberately changes the "
        "integration estimand while holding the source vessels, same-plate/same-batch controls and primary "
        "eligibility fixed. Plate integration uses biotic A380 change and flask integration uses biotic Fe(II) "
        "loss; neither receives abiotic subtraction. The common end is the earliest final scheduled reading "
        "among both formats and all required controls. Linear boundary interpolation uses adjacent observed "
        "values only; missing initial/internal support withholds the result. Physical replicates are averaged "
        "within each correspondence, then correspondences are averaged within dose pair. Differences are "
        "percentage points, calculated as common-window minus primary; the final column is the change in "
        "plate-minus-flask response. Reused cohorts and reference sets are not independent replicates. "
        "Exact and neighboring doses remain distinguished, and no slope/intercept or equivalence test is "
        "assigned: these differently measured and control-normalized observational cohorts cannot identify "
        "method agreement or isolate a causal format effect. Full-precision primary, alternative, support and "
        "window values remain available per vessel outside the numbered SI.",
        pd.DataFrame(output), ["additive", "plate_dose_mM", "flask_dose_mM", "Ga_mM", "match_status",
            "paired_context_n", "plate_response_change_pp", "flask_response_change_pp", "format_difference_change_pp"],
        ["af_iron_plate_trajectories", "validation_flask_points", "consolidated_activity_replicates"],
        ["7a", "7b", "7c", "7d", "7e", "7f"],
        ["Plate–flask relative-response discrepancies depend on whether each trajectory has its own observed endpoint or shares a common monitoring window."],
        {"vessel_results": details})


def build(tables):
    """Return two single-purpose analysis-result tables and unnumbered support."""
    return [_common_readout(tables), _window_response(tables)]
