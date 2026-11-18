"""Sulfur endpoint support."""
from __future__ import annotations
import numpy as np
import pandas as pd


SIGNALS = {"raw": "raw_R545_485", "corrected": "corrected_R545_485"}


PLATE_KEYS = ["species", "experiment_id", "group_id", "column", "batch_id", "plate_id"]


WELL_KEYS = PLATE_KEYS + ["ordinal", "well"]


STAGES = ("channel", "ratio", "baseline", "control", "net", "noise", "snr")


PRIMARY_REASONS = (
    "missing_derived_endpoint", "missing_or_ineligible_channel",
    "nonpositive_denominator", "near_zero_denominator", "invalid_ratio_signal",
    "missing_initial_acquisition", "invalid_biotic_baseline", "invalid_abiotic_baseline",
    "invalid_or_missing_baseline", "incomplete_matched_abiotic",
    "missing_or_nonfinite_net_contrast", "nonfinite_technical_noise",
    "zero_technical_noise", "negative_technical_noise", "nonfinite_technical_snr",
)


def _unique(frame, keys, label):
    if frame.duplicated(keys).any():
        raise ValueError(f"Duplicate {label}")


def _true(values):
    return values.eq(True) | values.astype(str).str.lower().eq("true")


def _joint(raw, corrected):
    raw, corrected = np.asarray(raw, dtype=bool), np.asarray(corrected, dtype=bool)
    return np.select([raw & corrected, raw & ~corrected, ~raw & corrected],
                     ["both", "raw_only", "corrected_only"], default="neither")


def _counts(values, prefix):
    return {f"{prefix}_{label}_n": int(values.eq(label).sum())
            for label in ("both", "raw_only", "corrected_only", "neither")}


def endpoint_well_support(channels, evidence):
    """One row per acquired endpoint biotic well and ratio readout.

    Eligibility stages are cumulative. The first failed stage provides one
    mutually exclusive primary reason; source statuses remain as diagnostics.
    The source evidence supplies the minimum eligible ratio denominator.
    Noise must be finite and strictly positive; technical scores must also be
    finite (a finite numerator/noise can still overflow).
    """
    acquired = channels[channels.kind.eq("biotic")].copy()
    _unique(acquired, WELL_KEYS, "acquired sulfur well/acquisition")
    acquired = acquired[acquired.ordinal.eq(
        acquired.groupby(PLATE_KEYS, dropna=False).ordinal.transform("max"))].copy()
    keep = WELL_KEYS + ["culture_id", "physical_observation_key", "medium_id",
                        "A545", "A485", "A750", "A545_eligible", "A485_eligible", "A750_eligible"]
    acquired = acquired[keep]
    source = evidence[evidence.signal.isin(SIGNALS.values())].copy()
    _unique(source, WELL_KEYS + ["signal"], "sulfur derived well/readout")
    diagnostics = ["signal_status", "ratio_denominator_min_AU", "baseline_status",
                   "controls_status", "biotic_signal", "biotic_baseline", "abiotic_baseline",
                   "abiotic_signal", "abiotic_n", "abiotic_rows", "net_change", "change_noise",
                   "technical_snr_signed", "technical_snr_magnitude",
                   "technical_snr_signed_common", "technical_snr_magnitude_common"]
    pieces = []
    for readout, signal in SIGNALS.items():
        q = acquired.merge(source[source.signal.eq(signal)][WELL_KEYS + diagnostics],
                           on=WELL_KEYS, how="left", validate="one_to_one", indicator=True)
        q["signal"] = signal
        q["readout"] = readout
        present = q.pop("_merge").eq("both")
        names = ["A545", "A485"] + (["A750"] if readout == "corrected" else [])
        channel_ok = pd.Series(True, index=q.index)
        for name in names:
            channel_ok &= np.isfinite(q[name]) & _true(q[name + "_eligible"])
        denominator = q.A485 - q.A750 if readout == "corrected" else q.A485
        q["endpoint_ratio_denominator_AU"] = denominator
        q["channel_eligible"] = present & channel_ok
        q["ratio_eligible"] = (q.channel_eligible & q.signal_status.eq("eligible") &
                               np.isfinite(q.biotic_signal) & np.isfinite(denominator) &
                               denominator.gt(q.ratio_denominator_min_AU))
        baseline_ok = (q.baseline_status.eq("eligible") & np.isfinite(q.biotic_baseline) &
                       np.isfinite(q.abiotic_baseline))
        q["baseline_eligible"] = q.ratio_eligible & baseline_ok
        control_ok = (q.controls_status.eq("complete") & q.abiotic_n.eq(3) &
                      q.abiotic_rows.eq(3) & np.isfinite(q.abiotic_signal))
        q["control_eligible"] = q.baseline_eligible & control_ok
        q["net_eligible"] = q.control_eligible & np.isfinite(q.net_change)
        q["noise_eligible"] = q.net_eligible & np.isfinite(q.change_noise) & q.change_noise.gt(0)
        q["snr_eligible"] = (q.noise_eligible & np.isfinite(q.technical_snr_signed) &
                            np.isfinite(q.technical_snr_magnitude))
        q["common_snr_eligible"] = (np.isfinite(q.technical_snr_signed_common) &
                                            np.isfinite(q.technical_snr_magnitude_common))
        reason = pd.Series("eligible", index=q.index)
        reason.loc[~q.snr_eligible] = "nonfinite_technical_snr"
        reason.loc[~q.noise_eligible & q.change_noise.lt(0)] = "negative_technical_noise"
        reason.loc[~q.noise_eligible & q.change_noise.eq(0)] = "zero_technical_noise"
        reason.loc[~q.noise_eligible & ~np.isfinite(q.change_noise)] = "nonfinite_technical_noise"
        reason.loc[~q.net_eligible] = "missing_or_nonfinite_net_contrast"
        reason.loc[~q.control_eligible] = "incomplete_matched_abiotic"
        reason.loc[~q.baseline_eligible] = q.baseline_status.where(
            q.baseline_status.isin(PRIMARY_REASONS), "invalid_or_missing_baseline")
        reason.loc[~q.ratio_eligible] = q.signal_status.where(
            q.signal_status.isin(PRIMARY_REASONS), "invalid_ratio_signal")
        reason.loc[~q.channel_eligible] = "missing_or_ineligible_channel"
        reason.loc[~present] = "missing_derived_endpoint"
        q["primary_exclusion_reason"] = reason
        if not q.snr_eligible.eq(reason.eq("eligible")).all():
            raise ValueError("Sulfur primary exclusion reasons do not partition eligibility")
        # Common S/N values require eligibility at every stage.
        if (q.common_snr_eligible & ~q.snr_eligible).any():
            raise ValueError("Common sulfur S/N includes a failed eligibility stage")
        pieces.append(q)
    return pd.concat(pieces, ignore_index=True)


def paired_endpoint_wells(wells):
    """Join actual raw/corrected endpoint memberships, retaining every candidate."""
    fields = [f"{stage}_eligible" for stage in STAGES] + [
        "net_change", "technical_snr_signed_common", "technical_snr_magnitude_common",
        "common_snr_eligible", "primary_exclusion_reason"]
    raw = wells[wells.readout.eq("raw")][WELL_KEYS + ["culture_id", "physical_observation_key", "medium_id"] + fields]
    corrected = wells[wells.readout.eq("corrected")][WELL_KEYS + fields]
    q = raw.merge(corrected, on=WELL_KEYS, how="outer", suffixes=("_raw", "_corrected"),
                  validate="one_to_one")
    for stage in STAGES:
        q[f"{stage}_joint_support"] = _joint(q[f"{stage}_eligible_raw"], q[f"{stage}_eligible_corrected"])
    q["paired_endpoint_well_eligible"] = (q.common_snr_eligible_raw &
                                         q.common_snr_eligible_corrected)
    if not q.paired_endpoint_well_eligible.eq(q.snr_joint_support.eq("both")).all():
        raise ValueError("Common sulfur S/N does not match joint finite endpoint support")
    return q


def correction_endpoint_support(wells, plate_evidence):
    """Return paired endpoint support by plate and by species and medium.

    Readout-marginal plate support requires at least one eligible well; common
    plate support requires at least one identical paired endpoint well. Those
    counts can differ when marginal well sets are disjoint. Net signs classify
    positive/negative agreement, opposition, and any-zero separately.
    """
    pair = paired_endpoint_wells(wells)
    reference = plate_evidence[plate_evidence.signal.isin(SIGNALS.values())]
    _unique(reference, PLATE_KEYS + ["signal"], "sulfur plate endpoint by readout")
    rows = []
    for identity, g in pair.groupby(PLATE_KEYS, sort=True, dropna=False):
        row = dict(zip(PLATE_KEYS, identity))
        row.update(medium_id=g.medium_id.iloc[0], final_ordinal=g.ordinal.iloc[0],
                   candidate_endpoint_well_n=len(g))
        for stage in STAGES:
            row.update(_counts(g[f"{stage}_joint_support"], f"endpoint_well_{stage}"))
            row[f"{stage}_joint_support"] = _joint(
                bool(g[f"{stage}_eligible_raw"].any()),
                bool(g[f"{stage}_eligible_corrected"].any())).item()
        common = g[g.paired_endpoint_well_eligible]
        row["common_eligible_endpoint_well_n"] = len(common)
        row["common_eligible_endpoint_plate"] = bool(len(common))
        for readout in SIGNALS:
            row[f"mean_net_response_{readout}_common"] = common[f"net_change_{readout}"].mean()
            row[f"mean_signed_snr_{readout}_common"] = common[f"technical_snr_signed_common_{readout}"].mean()
            row[f"mean_magnitude_snr_{readout}_common"] = common[f"technical_snr_magnitude_common_{readout}"].mean()
            row[f"marginal_eligible_endpoint_well_n_{readout}"] = int(g[f"snr_eligible_{readout}"].sum())
            for reason in PRIMARY_REASONS:
                row[f"endpoint_well_{readout}_{reason}_n"] = int(g[f"primary_exclusion_reason_{readout}"].eq(reason).sum())
        raw, corrected = (row[f"mean_net_response_{r}_common"] for r in SIGNALS)
        row["mean_net_response_sign_category"] = (
            "unavailable" if not len(common) else "zero" if raw == 0 or corrected == 0
            else "same" if np.sign(raw) == np.sign(corrected) else "opposite")
        row["mean_net_response_zero_category"] = (
            "both_zero" if len(common) and raw == 0 and corrected == 0 else
            "raw_zero_only" if len(common) and raw == 0 else
            "corrected_zero_only" if len(common) and corrected == 0 else "no_zero")
        rows.append(row)
    plates = pd.DataFrame(rows)
    # Compare plate-level results with the common wells used for the paired scores.
    for readout, signal in SIGNALS.items():
        check = plates.merge(reference[reference.signal.eq(signal)], on=PLATE_KEYS, how="outer",
                             suffixes=("_support", "_reference"), validate="one_to_one", indicator=True)
        if not check._merge.eq("both").all():
            raise ValueError("Acquired sulfur endpoint plate population differs from the plate evidence")
        if not check.final_ordinal_support.eq(check.final_ordinal_reference).all():
            raise ValueError("Acquired sulfur endpoint ordinal differs from the plate evidence")
        if not check.common_eligible_endpoint_well_n.eq(check.final_snr_common_n).all():
            raise ValueError("Common sulfur endpoint well count differs from the plate evidence")
        for kind in ("signed", "magnitude"):
            if not np.allclose(check[f"mean_{kind}_snr_{readout}_common"],
                               check[f"final_snr_{kind}_common"],
                               equal_nan=True, rtol=1e-10, atol=1e-12):
                raise ValueError(f"Common sulfur {kind} S/N differs from the plate evidence")
    summaries = []
    for (species, medium), g in plates.groupby(["species", "medium_id"], sort=True):
        common = g[g.common_eligible_endpoint_plate]
        row = dict(species=species, medium=medium,
                   candidate_endpoint_plate_formulation_n=len(g),
                   candidate_endpoint_well_n=int(g.candidate_endpoint_well_n.sum()),
                   common_eligible_endpoint_plate_formulation_n=len(common),
                   common_eligible_endpoint_well_n=int(g.common_eligible_endpoint_well_n.sum()))
        for readout in SIGNALS:
            row[f"marginal_eligible_endpoint_plate_formulation_n_{readout}"] = int(
                g[f"marginal_eligible_endpoint_well_n_{readout}"].gt(0).sum())
        for stage in STAGES:
            row.update(_counts(g[f"{stage}_joint_support"], f"endpoint_plate_{stage}"))
        for field in g:
            if field.startswith("endpoint_well_") or field.startswith("marginal_eligible_endpoint_well_n_"):
                row[field] = int(g[field].sum())
        for category in ("same", "opposite", "zero"):
            n = int(common.mean_net_response_sign_category.eq(category).sum())
            row[f"mean_net_response_{category}_sign_n"] = n
            row[f"mean_net_response_{category}_sign_percent"] = 100 * n / len(common) if len(common) else np.nan
        for category in ("both_zero", "raw_zero_only", "corrected_zero_only"):
            row[f"mean_net_response_{category}_n"] = int(common.mean_net_response_zero_category.eq(category).sum())
        row["mean_net_response_sign_denominator_n"] = len(common)
        row["mean_net_response_nonzero_sign_denominator_n"] = int(common.mean_net_response_sign_category.isin(["same", "opposite"]).sum())
        agreements = np.sign(common.mean_signed_snr_raw_common).eq(np.sign(common.mean_signed_snr_corrected_common))
        row["mean_signed_snr_direction_agreement_percent"] = 100 * agreements.mean() if len(common) else np.nan
        summaries.append(row)
    return plates, pd.DataFrame(summaries)


def baseline_endpoint_support(trajectories, registry):
    """Compute union, marginal and common endpoint-plate counts."""
    q = trajectories[trajectories.signal.isin(SIGNALS.values())].copy()
    q = q.merge(registry[["reference_id", "initial_OD600_200_uL"]],
                on="reference_id", validate="many_to_one")
    q["cohort"] = q.apply(lambda r: (
        f"AF thiosulfate OD {r.initial_OD600_200_uL:g}" if str(r.reference_id).startswith("af_ts_")
        else "AF tetrathionate 5 mM" if r.reference_id == "af_tt_5"
        else "AF tetrathionate 15 mM" if r.reference_id == "af_tt_15"
        else "AT tetrathionate 20 mM" if r.reference_id == "at_tt" else "AT thiosulfate 20 mM"), axis=1)
    keys = ["reference_id", "batch_id", "plate_id"]
    q = q[q.ordinal.eq(q.groupby(keys + ["signal"], dropna=False).ordinal.transform("max"))]
    _unique(q, keys + ["signal"], "baseline sulfur plate endpoint")
    q["endpoint_eligible"] = np.isfinite(q.biotic_change) & np.isfinite(q.abiotic_change)
    fields = ["ordinal", "endpoint_eligible"]
    raw = q[q.signal.eq(SIGNALS["raw"])][keys + ["cohort"] + fields]
    corr = q[q.signal.eq(SIGNALS["corrected"])][keys + ["cohort"] + fields]
    paired = raw.merge(corr, on=keys + ["cohort"], how="outer", suffixes=("_raw", "_corrected"), validate="one_to_one")
    for readout in SIGNALS:
        paired[f"endpoint_eligible_{readout}"] = paired[f"endpoint_eligible_{readout}"].fillna(False).astype(bool)
    paired["endpoint_joint_support"] = _joint(paired.endpoint_eligible_raw, paired.endpoint_eligible_corrected)
    paired["common_eligible_endpoint_plate"] = (paired.endpoint_joint_support.eq("both") &
                                                paired.ordinal_raw.eq(paired.ordinal_corrected))
    rows = []
    for cohort, g in paired.groupby("cohort", sort=True):
        for readout, label in (("raw", "R"), ("corrected", "R*")):
            row = dict(cohort=cohort, readout=label, candidate_endpoint_plate_n=len(g),
                       marginal_eligible_endpoint_plate_n=int(g[f"endpoint_eligible_{readout}"].sum()),
                       common_eligible_endpoint_plate_n=int(g.common_eligible_endpoint_plate.sum()),
                       endpoint_ordinal_mismatch_n=int((g.endpoint_joint_support.eq("both") &
                                                        ~g.ordinal_raw.eq(g.ordinal_corrected)).sum()))
            row.update(_counts(g.endpoint_joint_support, "endpoint_plate"))
            rows.append(row)
    return paired, pd.DataFrame(rows)


def build_tables(tables):
    """Build well- and plate-level eligibility tables for baseline and paired-correction endpoints."""
    wells = endpoint_well_support(tables["context_sulfur_channels"],
                                  tables["validation_sulfur_trajectory_evidence"])
    paired = paired_endpoint_wells(wells)
    plates, correction = correction_endpoint_support(wells, tables["sulfur_plate_evidence"])
    baseline_plates, baseline = baseline_endpoint_support(
        tables["sulfur_baseline_trajectories"], tables["sulfur_baseline_registry"])
    return {"sulfur_endpoint_well_support": wells,
            "sulfur_endpoint_paired_wells": paired,
            "sulfur_endpoint_plate_support": plates,
            "sulfur_correction_endpoint_support": correction,
            "sulfur_baseline_endpoint_plates": baseline_plates,
            "sulfur_baseline_endpoint_support": baseline}
