"""Ratio support."""
from __future__ import annotations
import numpy as np
import pandas as pd


DENOMINATOR_GUARD_AU = 1e-6


JOINT_CATEGORIES = ("both_valid", "raw_only", "corrected_only", "neither_valid")


INVALIDITY_PRIORITY = (
    "instrument_over_range", "instrument_parse_error", "required_channel_absent",
    "missing_or_nonfinite_required_channel", "required_channel_ineligible",
    "nonpositive_denominator", "near_zero_denominator",
)


def _true(values):
    return values.eq(True) | values.astype(str).str.lower().eq("true")


def _status_reason(status):
    tokens = set(str(status).split(";"))
    mapping = {
        "instrument_over_range": {"over_range", "censored_high"},
        "instrument_parse_error": {"parse_error"},
        "required_channel_absent": {"channel_absent"},
        "missing_or_nonfinite_required_channel": {"missing", "nonfinite"},
        "nonpositive_denominator": {"nonpositive_denominator"},
        "near_zero_denominator": {"near_zero_denominator"},
    }
    for reason in INVALIDITY_PRIORITY:
        if tokens & mapping.get(reason, set()):
            return reason
    return "eligible" if status == "eligible" else "required_channel_ineligible"


def classify_channel_ratios(channels):
    """Classify acquired observations without rejecting negative numerators.

    Required A545/A485/A750 channel values and optional per-channel eligibility,
    overrange, censored_high, and parse_error fields follow context_sulfur_channels.
    An explicit ``acquired=False`` row is excluded before computing any yield.
    Required absent channels differ from acquired but nonfinite channel values.
    """
    q = channels.copy()
    if "acquired" in q:
        q = q.loc[_true(q.acquired)].copy()
    absent = {}
    for wave in (545, 485, 750):
        field = f"A{wave}"
        absent[wave] = field not in q
        q[field] = pd.to_numeric(q[field], errors="coerce") if field in q else np.nan
    for prefix, waves in (("raw", (545, 485)), ("corrected", (545, 485, 750))):
        flags = {key: pd.Series(False, index=q.index) for key in INVALIDITY_PRIORITY[:5]}
        channel_ok = pd.Series(True, index=q.index)
        for wave in waves:
            field = f"A{wave}"
            finite = np.isfinite(q[field])
            eligible = _true(q[field + "_eligible"]) if field + "_eligible" in q else finite
            over = finite & q[field].gt(4.0)
            parse = pd.Series(False, index=q.index)
            for suffix in ("_overrange", "_censored_high"):
                if field + suffix in q:
                    over |= _true(q[field + suffix])
            if field + "_parse_error" in q:
                parse = _true(q[field + "_parse_error"])
            flags["instrument_over_range"] |= over
            flags["instrument_parse_error"] |= parse
            flags["required_channel_absent"] |= absent[wave]
            flags["missing_or_nonfinite_required_channel"] |= ~finite
            flags["required_channel_ineligible"] |= ~eligible
            channel_ok &= eligible & finite & ~over & ~parse
        denominator = q.A485 if prefix == "raw" else q.A485 - q.A750
        q[prefix + "_denominator"] = denominator
        q[prefix + "_numerator"] = q.A545 if prefix == "raw" else q.A545 - q.A750
        q[prefix + "_ok"] = channel_ok & denominator.gt(DENOMINATOR_GUARD_AU)
        flags["nonpositive_denominator"] = denominator.le(0)
        flags["near_zero_denominator"] = denominator.gt(0) & denominator.le(DENOMINATOR_GUARD_AU)
        q[prefix + "_primary_reason"] = np.select(
            [flags[reason] for reason in INVALIDITY_PRIORITY], INVALIDITY_PRIORITY,
            default="eligible",
        )
        if not q[prefix + "_ok"].eq(q[prefix + "_primary_reason"].eq("eligible")).all():
            raise ValueError(f"Unclassified {prefix} ratio eligibility")
    return _complete_pair_fields(q)


def _complete_pair_fields(q):
    q = q.copy()
    raw, corrected = q.raw_ok.to_numpy(bool), q.corrected_ok.to_numpy(bool)
    q["joint_category"] = np.select(
        [raw & corrected, raw & ~corrected, ~raw & corrected],
        JOINT_CATEGORIES[:3], default="neither_valid",
    )
    q["denominator_retained_fraction"] = (
        q.corrected_denominator.abs() / q.raw_denominator.abs().where(q.raw_denominator.ne(0))
    )
    q["denominator_guard_AU"] = DENOMINATOR_GUARD_AU
    q["eligibility_scope"] = "same_acquired_well_ratio_only"
    q["repeated_observation_unit"] = "physical plate/well/acquisition; not independent trials"
    return q


def pair_canonical_ratios(observations):
    """Pair optical ratios using the recorded channel eligibility flags."""
    q = observations.loc[observations.assay.eq("pH")].copy()
    if "acquired" in q:
        q = q.loc[_true(q.acquired)].copy()
    ids = ["acquisition_id", "physical_plate_id", "well"]
    raw = q.loc[q.signal.eq("ratio_raw")]
    corrected = q.loc[q.signal.eq("ratio_corrected")]
    if raw.duplicated(ids).any() or corrected.duplicated(ids).any():
        raise ValueError("Duplicate physical optical ratio observation")
    fields = ids + ["eligible", "denominator", "eligibility_status"]
    pair = raw.merge(corrected[fields], on=ids, suffixes=("_raw", "_corrected"),
                     validate="one_to_one")
    if len(pair) != len(raw) or len(pair) != len(corrected):
        raise ValueError("Unmatched acquired optical ratio")
    for prefix in ("raw", "corrected"):
        pair[prefix + "_ok"] = _true(pair["eligible_" + prefix])
        pair[prefix + "_denominator"] = pair["denominator_" + prefix]
        pair[prefix + "_primary_reason"] = pair["eligibility_status_" + prefix].map(_status_reason)
        if not pair[prefix + "_ok"].eq(pair[prefix + "_primary_reason"].eq("eligible")).all():
            raise ValueError(f"Inconsistent canonical {prefix} ratio status")
    pair["raw_numerator"] = pair.A545
    pair["corrected_numerator"] = pair.A545 - pair.A750
    pair["unit"] = pair.physical_plate_id
    pair["physical_observation_key"] = pair.acquisition_id.astype(str) + "|" + pair.physical_plate_id.astype(str) + "|" + pair.well.astype(str)
    return _complete_pair_fields(pair)


def summarize_ratio_yield(frame):
    """Summarize ratio availability and joint counts with equal weight per plate."""
    units = []
    for plate, p in frame.groupby("unit", sort=True):
        raw, corr = p.raw_ok.to_numpy(bool), p.corrected_ok.to_numpy(bool)
        cancel = pd.to_numeric(p.denominator_retained_fraction, errors="coerce")
        row = dict(unit=plate, attempt_n=len(p), raw_percent=100 * raw.mean(),
                   corrected_percent=100 * corr.mean(), loss_percent=100 * np.mean(raw & ~corr),
                   gain_percent=100 * np.mean(~raw & corr),
                   difference_pp=100 * np.mean(corr.astype(int) - raw.astype(int)),
                   raw_valid_n=int(raw.sum()), corrected_valid_n=int(corr.sum()),
                   cancellation_median=cancel[np.isfinite(cancel)].median())
        for category in JOINT_CATEGORIES:
            mask = p.joint_category.eq(category)
            row[category + "_n"] = int(mask.sum())
            row[category + "_percent"] = 100 * mask.mean()
        units.append(row)
    u = pd.DataFrame(units)
    if u.empty:
        raise ValueError("No acquired ratio observations")
    row = dict(plate_n=len(u), attempt_n=int(u.attempt_n.sum()))
    for key in ["raw_percent", "corrected_percent", "loss_percent", "gain_percent", "difference_pp", "cancellation_median"] + [x + "_percent" for x in JOINT_CATEGORIES]:
        row[key + "_mean"] = u[key].mean()
        row[key + "_min"] = u[key].min()
        row[key + "_max"] = u[key].max()
    for key in ["raw_valid_n", "corrected_valid_n"] + [x + "_n" for x in JOINT_CATEGORIES]:
        row[key] = int(u[key].sum())
    return row, u


def _coverage_details(observations, prefix, eligible=None):
    """Observed reference coverage, with levels retained as preparation identities."""
    group = observations if eligible is None else observations.loc[eligible]
    levels = []
    for level, wells in group.groupby("level_id", sort=True, dropna=False):
        reference = pd.to_numeric(wells.predictor, errors="coerce").to_numpy(float)
        finite = reference[np.isfinite(reference)]
        if len(np.unique(finite)) > 1:
            raise ValueError("Inconsistent reference within one media preparation level")
        mean_response = pd.to_numeric(wells.response, errors="coerce").mean()
        levels.append(dict(level_id=level, well_n=len(wells),
                           predictor=float(finite[0]) if len(finite) else np.nan,
                           response_mean=mean_response))
    level_frame = pd.DataFrame(levels, columns=["level_id", "well_n", "predictor", "response_mean"])
    represented = level_frame[np.isfinite(level_frame.predictor)]
    if eligible is not None:
        represented = represented[np.isfinite(represented.response_mean)]
    replicated = represented[represented.well_n.ge(2)]
    row = {prefix + "_well_n": len(group), prefix + "_level_n": len(represented),
           prefix + "_singleton_level_n": int(represented.well_n.eq(1).sum()),
           prefix + "_replicated_level_n": len(replicated)}
    for name, frame in ((prefix, represented), (prefix + "_replicated", replicated)):
        row[name + "_reference_min"] = frame.predictor.min()
        row[name + "_reference_max"] = frame.predictor.max()
        row[name + "_level_ids"] = ";".join(frame.level_id.astype(str))
    return row


def _media_acid_support(observations):
    acid = pd.to_numeric(observations.acid_M, errors="coerce")
    iron = observations.assay.eq("iron").all()
    if iron and not (np.isfinite(acid).all() and np.allclose(acid, .5, atol=1e-12, rtol=0)):
        raise ValueError("Ferric media source does not support the documented 0.5 M added-acid condition")
    return dict(
        nominal_added_h2so4_M_min=acid.min() if iron else np.nan,
        nominal_added_h2so4_M_max=acid.max() if iron else np.nan,
        acid_condition_scope=("0.5 M nominal added H2SO4 after mixing; not untreated medium" if iron else
                              "not_applicable_to_variable_pH_series; use measured pH reference; source acid_M is retained without interpreting it as reference-level acid concentration"),
        reference_quantity=("added ferric-stock Fe (mM); excludes medium iron background" if iron else "measured medium-specific pH before dye"),
    )


def _matching_media_observations(observations, row, signal):
    q = observations.loc[
        observations.acquisition_id.eq(row.acquisition_id)
        & observations.physical_plate_id.eq(row.physical_plate_id)
        & observations.condition.eq(row.condition) & observations.signal.eq(signal)
    ].copy()
    if q.empty or q.well.duplicated().any():
        raise ValueError("Missing or duplicate medium wells for the specified observation")
    return q


def build_media_analysis_support(tables):
    """Describe the standards and wells contributing to each medium analysis.

    The replicated-standard span requires at least two standards with replicate
    wells. A singleton standard can contribute to the common-well span. At
    least one replicated common standard is required for the pooled sample SD.
    """
    observations = tables["pub_media_observations"]
    performance = tables["pub_media_performance"]
    performance = performance.loc[
        (performance.assay.eq("pH") & performance.condition.eq("with_dye"))
        | (performance.assay.eq("iron") & performance.signal.isin(["A380", "A380_corrected"]))
    ]
    performance_rows = []
    for source_index, row in performance.iterrows():
        q = _matching_media_observations(observations, row, row.signal)
        scheduled = _coverage_details(q, "scheduled_reference")
        valid = _coverage_details(q, "individually_eligible", _true(q.eligible))
        if scheduled["scheduled_reference_well_n"] != row.scheduled_well_n or valid["individually_eligible_well_n"] != row.eligible_well_n:
            raise ValueError("Media performance source counts disagree with acquired observations")
        if valid["individually_eligible_replicated_level_n"] != row.eligible_level_n:
            raise ValueError("Media performance replicated-level support differs from source fit population")
        performance_rows.append(row.to_dict() | scheduled | valid | _media_acid_support(q) | dict(
            support_source_table="pub_media_performance", source_row_index=source_index,
            span_eligibility_status=("eligible_replicated_reference_range" if row.eligible_level_n >= 2 else
                                     "unavailable_fewer_than_two_replicated_levels"),
            numbered_table_span=(row.observed_response_range if row.eligible_level_n >= 2 else np.nan),
            numbered_table_span_to_sd=(row.range_to_technical_sd if row.eligible_level_n >= 2 else np.nan),
            numbered_table_rule="The span requires at least two standard levels with at least two eligible wells per level.",
        ))
    comparison = tables["pub_media_readout_comparison"]
    comparison = comparison.loc[
        (comparison.assay.eq("pH") & comparison.dye_uM.gt(0)
         & comparison.reference_signal.eq("ratio_raw") & comparison.candidate_signal.eq("ratio_corrected"))
        | (comparison.assay.eq("iron") & comparison.dye_uM.eq(0)
           & comparison.reference_signal.eq("A380") & comparison.candidate_signal.eq("A380_corrected"))
    ]
    paired_rows = []
    for source_index, row in comparison.iterrows():
        raw = _matching_media_observations(observations, row, row.reference_signal)
        corrected = _matching_media_observations(observations, row, row.candidate_signal)
        joined = raw.merge(corrected[["well", "eligible", "level_id", "predictor"]], on="well",
                           suffixes=("", "_candidate"), validate="one_to_one")
        if len(joined) != len(raw) or len(joined) != len(corrected):
            raise ValueError("Unmatched media correction source wells")
        if not joined.level_id.eq(joined.level_id_candidate).all() or not np.allclose(joined.predictor, joined.predictor_candidate, equal_nan=True):
            raise ValueError("Paired media wells use different reference levels")
        common = _true(joined.eligible) & _true(joined.eligible_candidate)
        scheduled = _coverage_details(joined, "scheduled_reference")
        ref_details = _coverage_details(raw, "reference_individually_eligible", _true(raw.eligible))
        cand_details = _coverage_details(corrected, "candidate_individually_eligible", _true(corrected.eligible))
        common_details = _coverage_details(joined, "common_eligible", common)
        expected = [row.scheduled_well_n, row.reference_eligible_n, row.candidate_eligible_n, row.common_eligible_n]
        actual = [len(joined), int(_true(raw.eligible).sum()), int(_true(corrected.eligible).sum()), int(common.sum())]
        if expected != actual:
            raise ValueError("Media correction source counts disagree with paired observations")
        ref_q, cand_q = row.reference_range_to_sd_common_wells, row.candidate_range_to_sd_common_wells
        ratio = cand_q / ref_q if np.isfinite(ref_q) and ref_q != 0 else np.nan
        estimable = bool(ref_q > 0 and cand_q > 0 and np.isfinite(ratio))
        if estimable:
            status = "estimable_positive_finite_Q_pair"
        elif common_details["common_eligible_well_n"] == 0:
            status = "unavailable_no_common_wells"
        elif common_details["common_eligible_level_n"] < 2:
            status = "unavailable_fewer_than_two_common_levels"
        elif common_details["common_eligible_replicated_level_n"] == 0:
            status = "unavailable_no_replicated_common_level_for_SD"
        else:
            status = "unavailable_nonpositive_or_nonfinite_Q"
        paired_rows.append(row.to_dict() | scheduled | ref_details | cand_details | common_details | _media_acid_support(joined) | dict(
            support_source_table="pub_media_readout_comparison", source_row_index=source_index,
            contrast_eligible=estimable, contrast_eligibility_status=status,
            contrast_log2_Q_ratio=float(np.log2(ratio)) if estimable else np.nan,
            common_variance_level_n=common_details["common_eligible_replicated_level_n"],
            numbered_table_rule="media_correction_effects uses identical common wells; singleton levels may enter span but do not supply within-level sample variance",
        ))
    paired = pd.DataFrame(paired_rows)
    if len(paired):
        paired["contrast_group_eligible_plate_n"] = paired.groupby(["medium_id", "reference_signal"]).contrast_eligible.transform("sum").astype(int)
        paired["contrast_group_inference_status"] = np.select(
            [paired.contrast_group_eligible_plate_n.eq(0), paired.contrast_group_eligible_plate_n.eq(1)],
            ["unavailable_no_eligible_plate", "unavailable_single_eligible_plate"],
            default="multiple_eligible_plates; group interval/test reported separately in media_correction_effects",
        )
    return {"performance_plate_support": pd.DataFrame(performance_rows), "paired_plate_support": paired}
