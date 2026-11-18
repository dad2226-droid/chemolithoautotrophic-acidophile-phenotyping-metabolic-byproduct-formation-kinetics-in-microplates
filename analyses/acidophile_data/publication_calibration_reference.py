"""Calibration reference."""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from .publication_calibration_stats import four_parameter_logistic


def fit_reference_4pl(x, y, *, fixed_bottom=False, decreasing=False, start=None):
    """Fit the same base-10, linear-predictor 4PL used by manuscript analyses.

    Numerical bounds are broad stability guards, not prior precision. Flagged
    fits remain available, with predictions and failures visible to the reader.
    """
    x, y = np.asarray(x, float), np.asarray(y, float)
    keep = np.isfinite(x) & np.isfinite(y)
    x, y = x[keep], y[keep]
    npar = 3 if fixed_bottom else 4
    record = dict(fit_success=False, fit_qc="insufficient_variable_levels", free_parameter_n=npar,
                  bottom_fixed_zero=bool(fixed_bottom), n_observations=len(x), n_levels=len(np.unique(x)))
    if len(np.unique(x)) <= npar or len(x) <= npar or np.ptp(x) <= 0 or np.ptp(y) <= 0:
        return record
    lo, scale, yr = float(x.min()), float(np.ptp(x)), float(np.ptp(y))
    z = (x - lo) / scale
    lb = np.array([float(y.min()-10*yr), 1e-12, -2., 1e-6 if decreasing else -100.])
    ub = np.array([float(y.max()+10*yr), 20*yr, 3., 100. if decreasing else -1e-6])
    sign = 1 if decreasing else -1
    if fixed_bottom:
        lb, ub = lb[1:], ub[1:]
    def unpack(p):
        return np.r_[0., p] if fixed_bottom else p
    def fun(p):
        return four_parameter_logistic(z, *unpack(p)) - y
    def jac(p):
        b, s, m, h = unpack(p)
        f = 1/(1+np.power(10., np.clip((z-m)*h, -300, 300)))
        j = np.column_stack([np.ones(len(z)), f, s*np.log(10)*h*f*(1-f),
                             -s*np.log(10)*(z-m)*f*(1-f)])
        return j[:, 1:] if fixed_bottom else j
    starts = []
    if start and start.get("fit_success"):
        p = np.array([start["bottom"], start["span"], (start["midpoint"]-lo)/scale,
                      start["slope"]*scale])
        starts.append(p[1:] if fixed_bottom else p)
    else:
        for slope in (2., 8.):
            p = np.array([float(y.min()), max(yr, y.max()) if fixed_bottom else yr, .5, sign*slope])
            starts.append(p[1:] if fixed_bottom else p)
    fits = []
    for p in starts:
        try:
            p = np.clip(p, lb+1e-10*(ub-lb), ub-1e-10*(ub-lb))
            fits.append(least_squares(fun, p, jac=jac, bounds=(lb, ub), max_nfev=800,
                                      ftol=1e-9, xtol=1e-9, gtol=1e-9))
        except (ValueError, FloatingPointError, np.linalg.LinAlgError):
            continue
    if not fits:
        return {**record, "fit_qc": "optimizer_failed"}
    best = min(fits, key=lambda a: (not a.success, np.sum(a.fun*a.fun)))
    b, s, m, h = unpack(best.x)
    norms = np.linalg.norm(best.jac, axis=0)
    cond = float(np.linalg.cond(best.jac/np.where(norms > 0, norms, 1)))
    distance = np.minimum((best.x-lb)/(ub-lb), (ub-best.x)/(ub-lb))
    flags = []
    if not best.success:
        flags.append("optimizer_not_converged")
    if np.any(distance < 1e-5):
        flags.append("parameter_at_numerical_bound")
    if not 0 <= m <= 1:
        flags.append("midpoint_outside_observed_range")
    if not np.isfinite(cond) or cond > 1e6:
        flags.append("weak_parameter_identification")
    sse = float(np.sum(best.fun**2))
    return {**record, "fit_success": bool(best.success), "fit_qc": ";".join(flags) or "no_numerical_flag",
            "bottom": float(b), "span": float(s), "midpoint": float(lo+scale*m), "slope": float(h/scale),
            "training_min": lo, "training_max": float(x.max()), "sse": sse,
            "r_squared": float(1-sse/np.sum((y-y.mean())**2)),
            "sy_x": float(np.sqrt(sse/(len(y)-npar))), "jacobian_condition_number": cond}


def inverse_reference_4pl(response, fit):
    """Invert inside model asymptotes, never clip to the training X range."""
    y = np.asarray(response, float)
    out = np.full(y.shape, np.nan)
    if not fit.get("fit_success"):
        return out
    f = (y-fit["bottom"])/fit["span"]
    valid = np.isfinite(f) & (f > 0) & (f < 1) & (abs(fit["slope"]) > 0)
    out[valid] = fit["midpoint"] + np.log10(1/f[valid]-1)/fit["slope"]
    return out


def _ci(values, requested):
    values = np.asarray(values, float)
    values = values[np.isfinite(values)]
    # A sparse set of successful bootstrap fits is not a reliable 95% interval.
    ok = requested >= 20 and len(values) >= max(20, int(np.ceil(.8*requested)))
    lo, hi = np.quantile(values, [.025, .975]) if ok else (np.nan, np.nan)
    return float(lo), float(hi), len(values)


def _means_from_groups(groups, rng=None):
    """Keep wavelength vectors paired; incomplete groups never become complete."""
    ans = []
    for matrix in groups:
        complete = np.isfinite(matrix).all(axis=0) & (len(matrix) == 4)
        selected = matrix if rng is None else matrix[rng.integers(0, len(matrix), len(matrix))]
        mean = np.mean(selected, axis=0)
        ans.append(np.where(complete, mean, np.nan))
    return np.asarray(ans)


def _attach_intervals(point, draws, keys, metrics, bootstrap_n):
    result = point.copy()
    result["bootstrap_requested_n"] = bootstrap_n
    index = pd.concat(draws, ignore_index=True).groupby(keys, dropna=False) if draws else None
    for metric in metrics:
        lower, upper, counts = [], [], []
        for row in result.to_dict("records"):
            key = tuple(row[k] for k in keys)
            try:
                values = index.get_group(key)[metric].to_numpy() if index is not None else []
            except KeyError:
                values = []
            lo, hi, n = _ci(values, bootstrap_n)
            lower.append(lo); upper.append(hi); counts.append(n)
        result[metric.removesuffix("_mM").removesuffix("_pH")+"_ci_low"] = lower
        result[metric.removesuffix("_mM").removesuffix("_pH")+"_ci_high"] = upper
        result[metric+"_bootstrap_valid_n"] = counts
    result["ci_method"] = "95%_percentile_technical_bootstrap; suppressed_if_B_below20_or_success_below80%"
    return result


def _ph_groups(wells):
    # Include every intended ratio and preserve missingness before resampling.
    data = wells[wells.signal.isin(("ratio_raw", "ratio_corrected")) & wells.dye_nominal_uM.gt(0)].copy()
    keys = ["dataset", "plate_id", "layout_group", "column"]
    records, matrices = [], []
    for key, g in data.groupby(keys, sort=True):
        first = g.iloc[0]
        records.append({**dict(zip(keys, key)), "acid_m": first.acid_m,
                        "reference_ph": first.reference_ph_after_dye, "source_file": first.source_file,
                        "reference_pH_measurement_stage": first.reference_pH_measurement_stage})
        wide = g.pivot(index="well", columns="signal", values="response").reindex(columns=["ratio_raw", "ratio_corrected"])
        matrices.append(wide.to_numpy(float))
    return pd.DataFrame(records), matrices


def _ph_predictions(meta, means, sd, warm=None):
    """Same two-plate training sets predict withheld acid and ferric matrices.

    All ferric points are evaluated by each acid holdout model. Fold replication
    is collapsed at the physical plate/level before accuracy summaries; it never
    creates additional experimental units.
    """
    rows, fits, cache = [], [], {}
    acid = meta.dataset.eq("acid_calibration").to_numpy()
    x = meta.reference_ph.to_numpy(float)
    plates = sorted(meta.loc[acid, "plate_id"].unique())
    for sidx, signal in enumerate(("ratio_raw", "ratio_corrected")):
        for plate in plates:
            train = acid & meta.plate_id.ne(plate).to_numpy()
            test = (~acid) | (acid & meta.plate_id.eq(plate).to_numpy())
            fit_id = f"{signal}|acid_holdout_{plate}"
            fit = fit_reference_4pl(x[train], means[train, sidx], decreasing=True,
                                    start=(warm or {}).get(fit_id))
            cache[fit_id] = fit
            fits.append({"fit_id": fit_id, "signal": signal, "held_out_acid_plate": plate, **fit})
            idx = np.flatnonzero(test)
            estimates = inverse_reference_4pl(means[idx, sidx], fit)
            for i, estimate in zip(idx, estimates):
                support = bool(fit.get("fit_success") and np.isfinite(x[i]) and fit["training_min"] <= x[i] <= fit["training_max"])
                if np.isfinite(estimate):
                    fraction = (means[i, sidx]-fit["bottom"])/fit["span"]
                    derivative = abs(fit["span"]*np.log(10)*fit["slope"]*fraction*(1-fraction))
                    technical_sd = sd[i, sidx]/derivative if derivative > 0 else np.nan
                else:
                    technical_sd = np.nan
                rows.append({**meta.iloc[i].to_dict(), "signal": signal, "fit_id": fit_id,
                    "held_out_acid_plate": plate,
                    "matrix": "acid_held_out_plate" if acid[i] else "ferric_transfer",
                    "response_mean": means[i, sidx], "estimated_pH": estimate, "error_pH": estimate-x[i],
                    "reference_in_training_range": support,
                    "estimate_in_training_range": bool(np.isfinite(estimate) and fit["training_min"] <= estimate <= fit["training_max"]),
                    "fit_qc": fit["fit_qc"], "local_delta_technical_sd_pH": technical_sd,
                    "prediction_status": "valid_inverse" if np.isfinite(estimate) else (
                        "fit_failed" if not fit.get("fit_success") else "incomplete_optical_mean" if not np.isfinite(means[i, sidx]) else "outside_model_asymptotes"),
                    "prediction_selection": "reference_support_only; estimated_pH_never_clipped_or_used_to_remove_error",
                    "comparator_design": "same_two_acid_plates_train_ferric_transfer_and_heldout_acid_predictions"})
    return pd.DataFrame(rows), pd.DataFrame(fits), cache


def _ph_collapse(predictions):
    keys = ["dataset", "plate_id", "layout_group", "column", "acid_m", "reference_ph", "source_file", "signal", "matrix"]
    rows = []
    for key, g in predictions.groupby(keys, dropna=False, sort=True):
        complete = g.error_pH.notna().all()
        rows.append({**dict(zip(keys, key)), "error_pH": g.error_pH.mean() if complete else np.nan,
                     "estimated_pH": g.estimated_pH.mean() if complete else np.nan,
                     "error_mse_pH2": float(np.mean(g.error_pH**2)) if complete else np.nan,
                     "local_delta_technical_sd_pH": np.sqrt(np.mean(g.local_delta_technical_sd_pH**2)) if complete else np.nan,
                     "reference_in_training_range": bool(g.reference_in_training_range.all()),
                     "fold_n": len(g), "valid_fold_n": int(g.error_pH.notna().sum()),
                     "fit_flagged": bool(g.fit_qc.ne("no_numerical_flag").any())})
    return pd.DataFrame(rows)


def _ph_support(collapsed, bin_width=.25):
    """Define comparable pH strata by REFERENCES, before looking at errors.

    Equal total weight is assigned to each common 0.25-pH stratum in both
    matrices, and equal weight to observations within a stratum. This avoids
    attributing different concentration-grid density to a matrix effect.
    """
    rows = []
    for acid_m in (.001, .01, .1, 1.):
        for signal in ("ratio_raw", "ratio_corrected"):
            acid = collapsed[collapsed.matrix.eq("acid_held_out_plate") & collapsed.signal.eq(signal)].copy()
            ferric = collapsed[collapsed.matrix.eq("ferric_transfer") & collapsed.signal.eq(signal) & collapsed.acid_m.eq(acid_m)].copy()
            valid_f = ferric[ferric.reference_in_training_range & np.isfinite(ferric.reference_ph)]
            valid_a = acid[acid.reference_in_training_range & np.isfinite(acid.reference_ph)]
            lo = max(valid_f.reference_ph.min(), valid_a.reference_ph.min()) if len(valid_f) and len(valid_a) else np.nan
            hi = min(valid_f.reference_ph.max(), valid_a.reference_ph.max()) if len(valid_f) and len(valid_a) else np.nan
            common = set()
            if np.isfinite(lo) and np.isfinite(hi) and lo <= hi:
                aa = valid_a[valid_a.reference_ph.between(lo, hi)]
                ff = valid_f[valid_f.reference_ph.between(lo, hi)]
                common = set(np.floor(aa.reference_ph/bin_width).astype(int)) & set(np.floor(ff.reference_ph/bin_width).astype(int))
            for matrix, d in (("acid_held_out_plate", acid), ("ferric_transfer", ferric)):
                for point in d.to_dict("records"):
                    ph = point["reference_ph"]
                    bin_id = int(np.floor(ph/bin_width)) if np.isfinite(ph) else -999
                    eligible = bool(point["reference_in_training_range"] and np.isfinite(lo) and lo <= ph <= hi and bin_id in common)
                    rows.append({**point, "comparison_acid_m": acid_m, "comparison_matrix": matrix,
                        "common_pH_min": lo, "common_pH_max": hi, "pH_bin": bin_id, "pH_bin_width": bin_width,
                        "common_pH_bin_n": len(common), "common_reference_support": eligible})
    return pd.DataFrame(rows)


def _ph_metrics(group):
    eligible = group[group.common_reference_support]
    complete = eligible[np.isfinite(eligible.error_pH)]
    row = dict(attempted_n=len(group), common_support_n=len(eligible), finite_prediction_n=len(complete),
               inverse_failure_n=len(eligible)-len(complete), common_pH_min=group.common_pH_min.iloc[0],
               common_pH_max=group.common_pH_max.iloc[0], common_pH_bin_n=int(group.common_pH_bin_n.iloc[0]))
    if complete.empty:
        return {**row, "bias_pH": np.nan, "rmse_pH": np.nan, "technical_sd_pH": np.nan,
                "status": "no_valid_common_reference_support_or_inverse_predictions"}
    binned = complete.groupby("pH_bin").agg(error=("error_pH", "mean"),
        mse=("error_mse_pH2", "mean"), variance=("local_delta_technical_sd_pH", lambda x: np.mean(x*x)))
    # A failed whole bin must not silently disappear from the equal-bin target.
    all_bins = set(eligible.pH_bin)
    full_bin_coverage = set(binned.index) == all_bins
    return {**row, "bias_pH": float(binned.error.mean()) if full_bin_coverage else np.nan,
            "rmse_pH": float(np.sqrt(binned.mse.mean())) if full_bin_coverage else np.nan,
            "technical_sd_pH": float(np.sqrt(binned.variance.mean())) if full_bin_coverage else np.nan,
            "status": "common_bin_weighted_finite_predictions; failures_reported" if full_bin_coverage else "failed_entire_common_pH_bin; metric_not_estimable"}


def _ph_summaries(predictions):
    collapsed = _ph_collapse(predictions)
    supported = _ph_support(collapsed)
    rows, contrasts = [], []
    keys = ["comparison_acid_m", "signal", "comparison_matrix"]
    for key, g in supported.groupby(keys, sort=True):
        rows.append(dict(acid_m=key[0], signal=key[1], matrix=key[2], **_ph_metrics(g),
                         weighting="equal_weight_per_common_0.25_pH_reference_bin; finite_points_equal_within_bin",
                         training_scope="two_acid_plates_per_fit; squared_errors_averaged_across_folds_before_physical_unit_summary",
                         uncertainty_scope="technical_whole_well_bootstrap_with_model_refitting; conditional_observed_plate_blocks_and_unquantified_pH_reference_error"))
    accuracy = pd.DataFrame(rows)
    for (acid_m, signal), g in accuracy.groupby(["acid_m", "signal"], sort=True):
        a = g[g.matrix.eq("acid_held_out_plate")].iloc[0]
        f = g[g.matrix.eq("ferric_transfer")].iloc[0]
        contrasts.append(dict(acid_m=acid_m, signal=signal, comparison="ferric_minus_heldout_acid",
            matrix_rmse_difference_pH=f.rmse_pH-a.rmse_pH, bias_difference_pH=f.bias_pH-a.bias_pH,
            technical_sd_ratio=f.technical_sd_pH/a.technical_sd_pH if a.technical_sd_pH > 0 else np.nan,
            acid_finite_n=a.finite_prediction_n, ferric_finite_n=f.finite_prediction_n,
            acid_failure_n=a.inverse_failure_n, ferric_failure_n=f.inverse_failure_n,
            common_pH_min=f.common_pH_min, common_pH_max=f.common_pH_max,
            interpretation="matrix_transfer_vs_heldout_acid_at_common_pH_support; technical_not_independent_preparation_comparison"))
    return accuracy, pd.DataFrame(contrasts), supported


def _ph_readout_contrasts(supported):
    """Compare corrections on the same physical pH-reference observations."""
    records = []
    ids = ["dataset", "plate_id", "layout_group", "column"]
    for (acid_m, matrix), group in supported.groupby(["comparison_acid_m", "comparison_matrix"]):
        raw = group[group.signal.eq("ratio_raw")].set_index(ids)
        corr = group[group.signal.eq("ratio_corrected")].set_index(ids)
        cols = ["error_pH", "error_mse_pH2", "local_delta_technical_sd_pH", "common_reference_support", "pH_bin"]
        pairs = raw[cols].join(corr[cols], lsuffix="_raw", rsuffix="_corrected")
        intended = pairs.common_reference_support_raw & pairs.common_reference_support_corrected
        finite = intended & pairs.error_pH_raw.notna() & pairs.error_pH_corrected.notna()
        p = pairs[finite]
        complete_bins = set(p.pH_bin_raw) == set(pairs.loc[intended, "pH_bin_raw"])
        result = dict(acid_m=acid_m, matrix=matrix, comparison="corrected_minus_raw",
                      attempted_pair_n=int(intended.sum()), paired_n=len(p),
                      raw_failure_n=int((intended & pairs.error_pH_raw.isna()).sum()),
                      corrected_failure_n=int((intended & pairs.error_pH_corrected.isna()).sum()),
                      rmse_difference_pH=np.nan, bias_difference_pH=np.nan, technical_sd_ratio=np.nan,
                      comparison_scope="identical_physical_plate_standard_pairs; common_reference_bins; failures_reported")
        if len(p) and complete_bins:
            b = p.groupby("pH_bin_raw").mean(numeric_only=True)
            result["rmse_difference_pH"] = float(np.sqrt(b.error_mse_pH2_corrected.mean())-np.sqrt(b.error_mse_pH2_raw.mean()))
            result["bias_difference_pH"] = float((b.error_pH_corrected-b.error_pH_raw).mean())
            # Pool squared technical SD within bins, then weight bins equally.
            raw_var = p.assign(v=p.local_delta_technical_sd_pH_raw**2).groupby("pH_bin_raw").v.mean().mean()
            corr_var = p.assign(v=p.local_delta_technical_sd_pH_corrected**2).groupby("pH_bin_raw").v.mean().mean()
            result["technical_sd_ratio"] = float(np.sqrt(corr_var/raw_var)) if raw_var > 0 else np.nan
        records.append(result)
    return pd.DataFrame(records)


def analyze_ph_matrix(wells, bootstrap_n=199, seed=1731):
    meta, groups = _ph_groups(wells)
    means = _means_from_groups(groups)
    sd = np.asarray([np.std(g, axis=0, ddof=1) if len(g) > 1 else [np.nan, np.nan] for g in groups])
    pred, fits, warm = _ph_predictions(meta, means, sd)
    accuracy, contrasts, support = _ph_summaries(pred)
    readouts = _ph_readout_contrasts(support)
    rng = np.random.default_rng(seed)
    draws_a, draws_c, draws_r = [], [], []
    for draw in range(bootstrap_n):
        selected = [g[rng.integers(0, len(g), len(g))] for g in groups]
        bmeans = _means_from_groups(selected)
        bmeans[~np.isfinite(means)] = np.nan
        bsd = np.asarray([np.std(g, axis=0, ddof=1) for g in selected])
        bp, _, _ = _ph_predictions(meta, bmeans, bsd, warm=warm)
        ba, bc, bs = _ph_summaries(bp)
        draws_a.append(ba); draws_c.append(bc)
        draws_r.append(_ph_readout_contrasts(bs))
        if (draw + 1) % 100 == 0:
            print(f"pH matrix reference bootstrap: {draw + 1}/{bootstrap_n}", flush=True)
    accuracy = _attach_intervals(accuracy, draws_a, ["acid_m", "signal", "matrix"],
                                 ["bias_pH", "rmse_pH", "technical_sd_pH"], bootstrap_n)
    contrasts = _attach_intervals(contrasts, draws_c, ["acid_m", "signal", "comparison"],
                                  ["matrix_rmse_difference_pH", "bias_difference_pH", "technical_sd_ratio"], bootstrap_n)
    readouts = _attach_intervals(readouts, draws_r, ["acid_m", "matrix", "comparison"],
                                 ["rmse_difference_pH", "bias_difference_pH", "technical_sd_ratio"], bootstrap_n)
    return {"ph_matrix_accuracy": accuracy, "ph_matrix_comparisons": contrasts,
            "ph_matrix_predictions": pred, "ph_matrix_common_support": support, "ph_matrix_fits": fits,
            "ph_matrix_readout_comparisons": readouts}
