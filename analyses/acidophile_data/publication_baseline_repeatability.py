"""Baseline repeatability."""
from __future__ import annotations
import json
import numpy as np
import pandas as pd
from scipy.stats import f, t


METRICS = ("AUC_to_own_max", "time_to_observed_max")


METHOD_URL = "https://www.itl.nist.gov/div898/handbook/eda/section3/eda359.htm"


LOGNORMAL_URL = "https://www.itl.nist.gov/div898/handbook/eda/section3/eda3669.htm"


PAIRED_METHOD_URL = "https://search.r-project.org/CRAN/refmans/PairedData/html/Var.test.html"


FLASK_PAIRS = ((0., 1.), (0., 100.), (1., 100.))


def baseline_plate_mask(frame):
    """Standard no-additive AFM1 controls, excluding deliberately altered pH.

    The condition table represents the standard preparation pH by a
    null nominal-pH field; recorded 1.0, 1.2 and 1.4 preparations are separate
    conditions. Accept an explicit 1.8 representation of the same target too.
    No additive, stock-pH, medium or species is inferred from a Ga value alone.
    """
    def eligible(row):
        if row["format"] != "plate" or not bool(row["no_addition"]):
            return False
        d = json.loads(row["condition_definition_json"])
        required = dict(species="AF", medium_id="AFM1", preparation_id="iron_no_dye",
                        signal="biotic_delta_A380", readout="iron")
        if any(d.get(key) != value for key, value in required.items()):
            return False
        if d.get("components_mM") or json.loads(row["components_mM_json"]):
            return False
        for key, target in (("initial_OD600_final", .001),
                            ("final_well_volume_uL", 200.),
                            ("incubation_temperature_C", 30.)):
            if d.get(key) is None or not np.isclose(float(d[key]), target):
                return False
        ph = d.get("nominal_initial_pH_json")
        return ph is None or np.isclose(float(ph), 1.8)
    return frame.apply(eligible, axis=1).astype(bool)


def baseline_batches(batches, metric):
    """Return one equally weighted observation per starting batch."""
    q = batches[batches.metric.eq(metric)].copy()
    plate = q[baseline_plate_mask(q)].copy()
    flask = q[q["format"].eq("flask") & q.no_addition].copy()
    for name, group in (("plate", plate), ("flask", flask)):
        if group.condition_id.nunique() > 1:
            raise ValueError(f"Multiple baseline {name} designs require explicit stratification")
        if group.batch_id.duplicated().any():
            raise ValueError(f"Duplicate baseline {name} batch: pseudo-replication is not allowed")
    return plate, flask


def _joined(frame, field):
    return ";".join(sorted({part for value in frame[field].dropna()
                            for part in str(value).split(";") if part}))


def flask_ga_batches(batches, metric, ga):
    """Select only the specified Ga-only AFM1 flask design, preserving batches."""
    q = batches[batches.metric.eq(metric) & batches["format"].eq("flask")].copy()
    def eligible(row):
        design = json.loads(row.condition_definition_json)
        if any(design.get(key) != value for key, value in
               dict(species="AF", medium_id="AFM1", signal="biotic_apparent_FeIII").items()):
            return False
        components = json.loads(row.components_mM_json)
        expected = {} if ga == 0 else {"Ga": ga}
        return (components == expected and design.get("components_mM") == expected
                and np.isclose(float(design.get("incubation_temperature_C", np.nan)), 30.))
    q = q[q.apply(eligible, axis=1)].copy()
    if q.condition_id.nunique() > 1:
        raise ValueError(f"Multiple {ga:g} mM Ga flask designs require stratification")
    if q.batch_id.duplicated().any():
        raise ValueError(f"Duplicate {ga:g} mM Ga flask batch: pseudo-replication is not allowed")
    return q


def _holm(result, family_size):
    """Retain unavailable planned comparisons in the prespecified family size."""
    result = result.copy()
    result["p_holm"] = np.nan
    order = result.p_value.dropna().sort_values().index
    adjusted = np.maximum.accumulate([(family_size-rank) * result.loc[index, "p_value"]
                                     for rank, index in enumerate(order)])
    result.loc[order, "p_holm"] = np.minimum(adjusted, 1.)
    result["reject_at_005_holm"] = result.p_holm.lt(.05)
    result["multiplicity_family_size"] = family_size
    return result


def _paired_log_variance(a, b):
    """Pitman-Morgan test and CI obtained by inversion on paired log outcomes.

    With X=log(a), Y=log(b), equal variances imply corr(X+Y,X-Y)=0.
    The exact t reference requires independent bivariate-normal batch pairs.
    For a hypothesized variance ratio r, apply the same test to X/sqrt(r), Y;
    solving its quadratic acceptance inequality yields the marginal ratio CI.
    """
    a, b = np.asarray(a, float), np.asarray(b, float)
    row = dict(status="estimated", p_value=np.nan, log_variance_ratio=np.nan,
               ci_low=np.nan, ci_high=np.nan, statistic=np.nan,
               df=len(a)-2 if len(a)>=3 else np.nan)
    if len(a) != len(b):
        raise ValueError("Paired variance test requires equal paired input lengths")
    if len(a) < 3:
        return {**row, "status": "fewer_than_three_matched_starting_batches"}
    if np.any(~np.isfinite(a)) or np.any(~np.isfinite(b)):
        return {**row, "status": "nonfinite_native_batch_metric"}
    if np.any(a <= 0) or np.any(b <= 0):
        return {**row, "status": "nonpositive_native_batch_metric_log_test_unavailable"}
    x, y = np.log(a), np.log(b)
    va, vb = np.var(x, ddof=1), np.var(y, ddof=1)
    cov = np.cov(x, y, ddof=1)[0, 1]
    det = va*vb - cov*cov
    row.update(numerator_log_variance=va, denominator_log_variance=vb,
               paired_log_covariance=cov)
    if va <= 0 or vb <= 0 or det <= np.finfo(float).eps*max(va*vb, np.finfo(float).tiny):
        return {**row, "status": "degenerate_paired_log_covariance"}
    ratio = va/vb
    statistic = np.sqrt(len(a)-2)*(va-vb)/(2*np.sqrt(det))
    spread = t.ppf(.975, len(a)-2)*np.sqrt(det/(len(a)-2))/vb
    upper = (np.sqrt(ratio+spread*spread)+spread)**2
    lower = ratio*ratio/upper
    row.update(log_variance_ratio=ratio, statistic=statistic,
               p_value=2*t.sf(abs(statistic), len(a)-2), ci_low=lower, ci_high=upper)
    return row


def flask_pairwise_repeatability_tests(batches):
    """Six exploratory metric-by-Ga contrasts, retaining shared-batch pairing.

    The Ga groups share starter cultures and acquisition batches, so an
    independent-samples F test would be inappropriate. Complete matched batch
    pairs support Pitman-Morgan tests of natural-log variance; incomplete
    batches remain identified and are not silently treated as independent.
    """
    rows = []
    for metric in METRICS:
        for ga_a, ga_b in FLASK_PAIRS:
            groups = [flask_ga_batches(batches, metric, ga) for ga in (ga_a, ga_b)]
            usable = [g[g.eligible & g.raw_value.notna()].set_index("batch_id", drop=False)
                      for g in groups]
            common = sorted(set(usable[0].index) & set(usable[1].index))
            paired = [g.loc[common] for g in usable]
            row = dict(metric=metric, comparison=f"AFM1 flask {ga_a:g} / AFM1 flask {ga_b:g}",
                       numerator_Ga_mM=ga_a, denominator_Ga_mM=ga_b,
                       method="two-sided paired Pitman-Morgan test of natural-log native batch variance",
                       method_url=PAIRED_METHOD_URL, alpha=.05, confidence_level=.95,
                       test_target="between-batch multiplicative dispersion within matched starting batches",
                       multiplicity_family="six exploratory flask pair by metric contrasts",
                       matched_batch_n=len(common), matched_batch_ids=";".join(common),
                       assumptions="independent starting-batch pairs; bivariate-normal log native batch outcomes",
                       limitation="three or fewer matched batches; unmatched batches excluded from paired inference only; failure to reject is not equivalence")
            for name, candidate, use, match in zip(("numerator", "denominator"), groups, usable, paired):
                x = match.raw_value.to_numpy(float)
                row.update({name+"_candidate_batch_n": len(candidate),
                            name+"_usable_batch_n": len(use),
                            name+"_matched_batch_n": len(match),
                            name+"_unmatched_batch_ids": ";".join(sorted(set(use.index)-set(common))),
                            name+"_condition_id": _joined(candidate, "condition_id"),
                            name+"_native_unit": _joined(candidate, "native_unit"),
                            name+"_unit_ids": _joined(match, "unit_ids"),
                            name+"_culture_ids": _joined(match, "culture_ids"),
                            name+"_source_record_ids": _joined(match, "source_record_ids"),
                            name+"_physical_unit_n": int(match.physical_unit_n.sum()),
                            name+"_native_mean": float(np.mean(x)) if len(x) else np.nan,
                            name+"_native_sd": float(np.std(x, ddof=1)) if len(x)>1 else np.nan,
                            name+"_cv_percent": 100*float(np.std(x, ddof=1)/np.mean(x))
                                if len(x)>1 and np.mean(x)>0 else np.nan})
            row.update(_paired_log_variance(*(g.raw_value.to_numpy(float) for g in paired)))
            rows.append(row)
    return _holm(pd.DataFrame(rows), len(METRICS)*len(FLASK_PAIRS))


def baseline_repeatability_tests(batches):
    """Two planned, two-sided log-variance tests with Holm-adjusted P values.

    Native nonpositive measurements cause an unavailable test, rather than being
    silently removed. The 95% variance-ratio intervals are marginal intervals;
    only the P values are adjusted across the two planned metrics.
    """
    rows = []
    for metric in METRICS:
        groups = baseline_batches(batches, metric)
        row = dict(metric=metric, comparison="AFM1 plate 0 / AFM1 flask 0",
                   method="two-sided F test of natural-log native batch variance",
                   test_target="between-batch multiplicative dispersion",
                   alpha=.05, confidence_level=.95,
                   assumptions="independent starting batches; positive lognormal native batch outcomes",
                   limitation="observational formats differ in vessel, readout and acquisition schedule; small batch samples; failure to reject is not equivalence",
                   method_url=METHOD_URL, lognormal_url=LOGNORMAL_URL,
                   status="estimated", p_value=np.nan,
                   log_variance_ratio=np.nan, ci_low=np.nan, ci_high=np.nan)
        values = []
        for name, group in zip(("plate", "flask"), groups):
            valid = group[group.eligible & group.raw_value.notna()]
            x = valid.raw_value.to_numpy(float)
            values.append(x)
            row.update({name + "_batch_n": len(x),
                        name + "_candidate_batch_n": len(group),
                        name + "_physical_unit_n": int(valid.physical_unit_n.sum()),
                        name + "_eligible_physical_unit_n": int(valid.eligible_unit_n.sum()),
                        name + "_batch_ids": _joined(valid, "batch_id"),
                        name + "_unit_ids": _joined(valid, "unit_ids"),
                        name + "_culture_ids": _joined(valid, "culture_ids"),
                        name + "_source_record_ids": _joined(valid, "source_record_ids"),
                        name + "_condition_id": _joined(group, "condition_id"),
                        name + "_native_unit": _joined(group, "native_unit"),
                        name + "_native_mean": float(np.mean(x)) if len(x) else np.nan,
                        name + "_native_sd": float(np.std(x, ddof=1)) if len(x) > 1 else np.nan,
                        name + "_cv_percent": 100 * float(np.std(x, ddof=1) / np.mean(x))
                             if len(x) > 1 and np.mean(x) > 0 else np.nan})
        a, b = values
        if min(len(a), len(b)) < 2:
            row["status"] = "fewer_than_two_usable_batches"
        elif np.any(a <= 0) or np.any(b <= 0):
            row["status"] = "nonpositive_native_batch_metric_log_test_unavailable"
        else:
            va, vb = np.var(np.log(a), ddof=1), np.var(np.log(b), ddof=1)
            row.update(plate_log_variance=va, flask_log_variance=vb,
                       numerator_df=len(a)-1, denominator_df=len(b)-1)
            if va <= 0 or vb <= 0:
                row["status"] = "degenerate_log_variance"
            else:
                ratio = float(va / vb)
                df = (len(a)-1, len(b)-1)
                row.update(log_variance_ratio=ratio,
                           p_value=min(1., 2 * min(f.cdf(ratio, *df), f.sf(ratio, *df))),
                           ci_low=ratio / f.ppf(.975, *df),
                           ci_high=ratio / f.ppf(.025, *df))
        rows.append(row)
    result = pd.DataFrame(rows)
    result["p_holm"] = np.nan
    # Keep both planned comparisons in the family even if one is unavailable.
    order = result.p_value.dropna().sort_values().index
    adjusted = np.maximum.accumulate([(len(METRICS)-rank) * result.loc[index, "p_value"]
                                     for rank, index in enumerate(order)])
    result.loc[order, "p_holm"] = np.minimum(adjusted, 1.)
    result["reject_at_005_holm"] = result.p_holm.lt(.05)
    return result
