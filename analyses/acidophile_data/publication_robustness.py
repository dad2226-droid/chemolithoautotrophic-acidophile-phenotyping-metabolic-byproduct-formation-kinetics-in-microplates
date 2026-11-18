"""Robustness."""
from __future__ import annotations
import hashlib
import json
import numpy as np
import pandas as pd
from openpyxl import load_workbook
from .optical import load_media_references
from .plate_reader import read_plate


ABSORBANCE_LIMIT_AU = 4.0


DENOMINATOR_GUARD_AU = 1e-6


SCOPE = ("Conditional optical/product-associated response; reference uncertainty is not "
         "quantified. Sequential acquisition, thermal history and handling are confounded; "
         "no equivalence or isolated causal effect is established.")


CURVE_KEYS = ["acquisition_id", "condition", "signal"]


META_KEYS = ["acquisition_id", "physical_plate_id", "reference_set_id", "source_path",
             "source_sheet", "experiment_id", "assay", "family", "condition", "signal",
             "run_number", "plate_type", "temperature_C", "acid_M", "dye_uM",
             "reading_label", "step_id", "step_order", "read_index", "recorded_time",
             "medium_id", "condition_role", "group_id", "independent_preparations"]


def _number(value):
    return float(value) if value is not None else np.nan


def _channels(plate, well):
    values, flags = {}, {}
    for wavelength in (380, 428, 485, 545, 750):
        try:
            value = plate.value(well, wavelength)
            flag = ("over_range" if plate.is_censored_high(well, wavelength) or value > ABSORBANCE_LIMIT_AU else
                    "parse_error" if plate.is_parse_error(well, wavelength) else
                    "ok" if np.isfinite(value) else "missing")
        except KeyError:
            value, flag = np.nan, "channel_absent"
        values[wavelength], flags[wavelength] = (value if flag == "ok" else np.nan), flag
    return values, flags


def _signal_records(values, flags, assay):
    specifications = ({"ratio_raw": (545, 485, False), "ratio_corrected": (545, 485, True)}
                      if assay == "pH" else
                      {"A380": (380, None, False), "A380_corrected": (380, None, True),
                       "A428": (428, None, False), "A428_corrected": (428, None, True)})
    result = []
    for name, (numerator, denominator, corrected) in specifications.items():
        required = [numerator] + ([denominator] if denominator else []) + ([750] if corrected else [])
        errors = sorted({flags[w] for w in required if flags[w] != "ok"})
        a = values[numerator] - (values[750] if corrected else 0)
        b = (values[denominator] - (values[750] if corrected else 0)) if denominator else 1.
        if denominator and not errors and b <= DENOMINATOR_GUARD_AU:
            errors = ["nonpositive_denominator" if b <= 0 else "near_zero_denominator"]
        value = a / b if not errors else np.nan
        result.append({"signal": name, "response": value, "eligible": not errors,
                       "eligibility_status": ";".join(errors) if errors else "eligible",
                       "denominator": b if denominator else np.nan,
                       "denominator_guard_AU": DENOMINATOR_GUARD_AU if denominator else np.nan})
    return result


def _lid_cells(path):
    """Read exact worksheet geometry, retaining OVER and source-cell provenance."""
    workbook = load_workbook(path, data_only=True, read_only=True)
    result = {}
    try:
        for sheet in workbook:
            rows = list(sheet.values)
            if not any(any(x is not None for x in row) for row in rows):
                continue
            cells = {}
            for start, row in enumerate(rows):
                if not row or row[0] != "Measurement Wavelength":
                    continue
                wavelength = int(row[4])
                headers = [i for i in range(start + 1, min(start + 26, len(rows)))
                           if rows[i] and rows[i][0] == "<>"]
                if len(headers) != 1:
                    raise ValueError(f"Ambiguous lid matrix {sheet.title}: {wavelength}")
                header = headers[0]
                if list(rows[header][1:13]) != list(range(1, 13)):
                    raise ValueError(f"Invalid lid matrix columns: {sheet.title}")
                for i, letter in enumerate("ABCDEFGH", 1):
                    data = rows[header + i]
                    if data[0] != letter:
                        raise ValueError(f"Invalid lid matrix row: {sheet.title}")
                    for col in range(1, 13):
                        raw = data[col]
                        value = float(raw) if isinstance(raw, (int, float)) and not isinstance(raw, bool) else np.nan
                        if str(raw).upper() == "OVER" or value > ABSORBANCE_LIMIT_AU:
                            flag, value = "over_range", np.nan
                        elif np.isfinite(value):
                            flag = "ok"
                        else:
                            flag, value = "missing_or_nonnumeric", np.nan
                        cells[(letter + str(col), wavelength)] = (value, flag, f"{chr(65 + col)}{header + i + 1}")
            if len(cells) != 576:
                raise ValueError(f"Incomplete lid worksheet: {sheet.title}")
            result[sheet.title] = cells
    finally:
        workbook.close()
    return result


def build_robustness_observations(root, registry):
    records = []
    selected = registry[(registry.experiment_id.str.contains("configuration")) |
                        registry.experiment_id.isin(["final_indicator_acid", "final_indicator_ferric", "paired_lid_comparison"])].copy()
    workbook = {}
    for row in selected.to_dict("records"):
        acquisition = json.loads(row["acquisition_metadata_json"])
        if row["reading_label"] == "spectrum" or "spectrum" in row["readout_id"]:
            continue
        experiment = json.loads(row["experiment_metadata_json"])
        preparation = experiment["preparation"]
        physical_plate = json.loads(row.get("physical_plate_metadata_json", "{}"))
        reference = json.loads(row["reference_metadata_json"])[row["reference_set_id"]]
        assay = "pH" if row["family"] == "ratiometric_pH_indicator" or row.get("acquisition_family") == "acid" else "iron"
        replication = experiment.get("replication", {})
        independent_preparations = replication.get("independent_standard_preparations_by_family", {}).get(
            row.get("acquisition_family"), replication.get("independent_standard_preparations", 1)
        )
        if row["source_sheet"]:
            if row["source_path"] not in workbook:
                workbook[row["source_path"]] = _lid_cells(root / row["source_path"])
            cells = workbook[row["source_path"]][row["source_sheet"]]
            times = acquisition["channel_start_times"]
            timestamp = max(x["measurement_start"] for x in times)
            plate = None
        else:
            plate = read_plate(root / row["source_path"])
            timestamp = acquisition.get("date_last_saved", "")
        base = {key: row.get(key, "") for key in ("acquisition_id", "physical_plate_id", "reference_set_id", "source_path", "source_sheet", "experiment_id", "family", "run_number", "reading_label")}
        base.update(assay=assay, plate_type=physical_plate.get("plate_type", preparation.get("plate_type", "UV_transparent")),
                    temperature_C=preparation.get("condition_temperature_C", 23),
                    acid_M=preparation.get("declared_h2so4_M", _number(preparation.get("declared_h2so4_mM")) / 1000),
                    step_id=acquisition.get("incubation_step_id", ""), step_order=acquisition.get("incubation_step_order", np.nan),
                    read_index=acquisition.get("read_index", np.nan), recorded_time=timestamp,
                    medium_id="", group_id="", condition_role="calibration", independent_preparations=independent_preparations,
                    preparation_series_id=physical_plate.get("preparation_series_id", ""),
                    source_sha256=hashlib.sha256((root / row["source_path"]).read_bytes()).hexdigest())
        for group_index, group in enumerate(reference["row_groups"]):
            condition = group["condition"]
            dye = preparation.get("calculated_dyed_well_thymol_blue_uM", 9.523809523809524) if condition == "thymol_blue" else 0
            # Ratios of undyed acid matrix controls are retained in raw channels,
            # not fitted as an indicator calibration.
            if assay == "pH" and condition != "thymol_blue":
                continue
            report = group["reported_concentration"]
            for position, col in enumerate(group["plate_columns"]):
                if assay == "pH":
                    predictors = group.get("reference_pH_by_column")
                    predictor = _number(predictors[position]) if predictors else np.nan
                    basis = group.get("reference_pH_status", "not_reported")
                else:
                    predictors = report.get("values_by_column")
                    predictor = _number(predictors[position]) if predictors else np.nan
                    predictor *= report.get("final_well_concentration_multiplier", 1)
                    predictor *= 1000 if report.get("unit") == "M" else 1
                    basis = "total_Fe_mM_final_well; " + report.get("status", "")
                for letter in group["rows"]:
                    well = letter + str(col)
                    if plate is None:
                        values = {w: cells[(well, w)][0] for w in (380, 428, 485, 545, 750)}
                        flags = {w: cells[(well, w)][1] for w in values}
                    else:
                        values, flags = _channels(plate, well)
                    for signal in _signal_records(values, flags, assay):
                        # Configuration/storage A428 comparisons duplicate Figure 2;
                        # retain these channels but the present curve models use A380.
                        if assay == "iron" and signal["signal"].startswith("A428"):
                            continue
                        records.append({**base, **signal, "condition": condition, "dye_uM": dye,
                                        "level_id": f"{group_index}:{col}", "row_group": "".join(group["rows"]),
                                        "well": well, "row": letter, "column": col,
                                        "predictor": predictor, "predictor_unit": "pH" if assay == "pH" else "mM_total_Fe",
                                        "predictor_basis": basis, "final_well_volume_uL": group["final_well_volume_uL"],
                                        **{f"A{w}": value for w, value in values.items()}})
    result = pd.DataFrame(records)
    if result.duplicated(["acquisition_id", "well", "signal"]).any():
        raise ValueError("Duplicate canonical robustness well assignments")
    return result


def build_media_observations(root, registry):
    references = load_media_references(root).set_index("reference_ID")
    records = []
    for row in registry[registry.family.eq("media_validation")].to_dict("records"):
        preparation = json.loads(row["preparation_metadata_json"])
        experiment = json.loads(row["experiment_metadata_json"])
        groups = json.loads(row["group_metadata_json"])
        refs = json.loads(row["reference_metadata_json"])
        plate = read_plate(root / row["source_path"])
        for group_id, group in groups.items():
            assay = "pH" if group["assay"] == "pH" else "iron"
            prep = preparation["well_preparations"][group["well_preparation_id"]]
            volume = prep["final_well_volume_uL"]
            for i, letter in enumerate(group["plate_rows"]):
                fraction = [1, .5, .25, .125, .0625, .03125, .015625, 0][i]
                reference_id = (refs[group["reference_set_id"]]["reference_ids_by_row"][i] if assay == "pH" else group["stock_reference_id"])
                predictor = float(references.loc[reference_id, "value"])
                if assay == "iron":
                    predictor *= fraction * preparation["calibration_solution_volume_uL"] / volume * 1000
                for column in group["technical_replicate_columns"]:
                    well = letter + str(column)
                    values, flags = _channels(plate, well)
                    for signal in _signal_records(values, flags, assay):
                        records.append({**signal, "acquisition_id": row["acquisition_id"] + ":" + group_id,
                                        "physical_plate_id": row["physical_plate_id"], "reference_set_id": group.get("reference_set_id", reference_id),
                                        "reference_id": reference_id, "source_path": row["source_path"], "source_sheet": "", "experiment_id": row["experiment_id"],
                                        "assay": assay, "family": "media_validation", "condition": group["well_preparation_id"],
                                        "run_number": row["plate_replicate_number"], "plate_type": "non_UV_transparent",
                                        "temperature_C": 23, "acid_M": .5 if assay == "iron" else fraction * 100 / volume,
                                        "dye_uM": prep.get("calculated_thymol_blue_concentration_uM", 0),
                                        "reading_label": "single_endpoint", "step_id": "", "step_order": np.nan, "read_index": 0,
                                        "recorded_time": row["recorded_datetime"], "medium_id": group["medium_id"],
                                        "group_id": group_id, "condition_role": group["condition_role"], "independent_preparations": 1,
                                        "level_id": letter, "row_group": "ABCDEFGH", "well": well, "row": letter, "column": column,
                                        "predictor": predictor, "predictor_unit": "pH" if assay == "pH" else "mM_added_Fe",
                                        "predictor_basis": "measured_medium_specific_pH_before_dye" if assay == "pH" else "measured_total_Fe_in_ferric_stock; added_Fe_excludes_medium_background",
                                        "final_well_volume_uL": volume, "relative_stock_fraction": fraction,
                                        "medium_strength_relative_to_standard": 200 / volume,
                                        **{f"A{w}": value for w, value in values.items()}})
    return pd.DataFrame(records)


def summarize_levels(observations):
    records = []
    for _, group in observations.groupby(CURVE_KEYS + ["level_id"], dropna=False, sort=True):
        first = group.iloc[0]
        good = group[group.eligible]
        records.append({**{key: first[key] for key in META_KEYS}, "level_id": first.level_id,
                        "row_group": first.row_group, "column": first.column,
                        "predictor": first.predictor, "predictor_unit": first.predictor_unit,
                        "predictor_basis": first.predictor_basis,
                        "response_mean": good.response.mean(), "technical_sd": good.response.std(ddof=1),
                        "technical_n": len(good), "scheduled_well_n": len(group), "complete_level": len(good) == len(group),
                        "excluded_well_n": len(group) - len(good),
                        "exclusion_reasons": ";".join(sorted(set(group.loc[~group.eligible, "eligibility_status"])))})
    return pd.DataFrame(records)


def build_matches(observations):
    """Prespecified same-plate pairs, plus same-preparation plate-material pairs."""
    reads = observations.drop_duplicates("acquisition_id").copy()
    pairs = []
    for _, candidate in reads.iterrows():
        common = reads[reads.assay.eq(candidate.assay) & reads.physical_plate_id.eq(candidate.physical_plate_id)]
        kind = None
        if candidate.family == "lid_comparison" and candidate.reading_label == "lid_present":
            reference = common[common.reading_label.eq("lid_absent")]
            kind = "lid_present_minus_absent"
        elif candidate.experiment_id.startswith("final_") and candidate.read_index > 0:
            reference = common[common.read_index.eq(0) & common.family.eq(candidate.family)]
            kind = "storage_minus_initial"
        elif "configuration" in candidate.experiment_id and candidate.step_id != "initial_23_C":
            ref_step = "overnight_23_C" if candidate.step_id == "membrane_23_C" else "initial_23_C"
            reference = common[common.step_id.eq(ref_step)]
            kind = "membrane_minus_overnight" if candidate.step_id == "membrane_23_C" else "configuration_minus_initial"
        elif "configuration" in candidate.experiment_id and candidate.plate_type == "non_UV_transparent":
            reference = reads[reads.assay.eq(candidate.assay) & reads.reference_set_id.eq(candidate.reference_set_id)
                              & reads.plate_type.eq("UV_transparent") & reads.step_id.eq(candidate.step_id)
                              & reads.dye_uM.eq(candidate.dye_uM)]
            kind = "nonUV_minus_UV"
        if kind:
            pairs.append({"comparison_id": f"robust_pair_{len(pairs) + 1:04d}", "comparison": kind,
                          "candidate_acquisition_id": candidate.acquisition_id,
                          "reference_acquisition_id": reference.acquisition_id.iloc[0] if len(reference) == 1 else "",
                          "pair_status": "matched" if len(reference) == 1 else "missing_reference" if not len(reference) else "ambiguous_reference",
                          "candidate_physical_plate_id": candidate.physical_plate_id,
                          "reference_physical_plate_id": reference.physical_plate_id.iloc[0] if len(reference) == 1 else "",
                          "assay": candidate.assay, "step_id": candidate.step_id,
                          "reference_set_id": candidate.reference_set_id,
                          "same_physical_plate": bool(len(reference) == 1 and reference.physical_plate_id.iloc[0] == candidate.physical_plate_id),
                          "scope": SCOPE})
    # Material comparison is independent of the within-plate sequential pairing.
    for _, candidate in reads[reads.experiment_id.str.contains("configuration") & reads.plate_type.eq("non_UV_transparent") & ~reads.step_id.eq("initial_23_C")].iterrows():
        reference = reads[reads.assay.eq(candidate.assay) & reads.reference_set_id.eq(candidate.reference_set_id)
                          & reads.plate_type.eq("UV_transparent") & reads.step_id.eq(candidate.step_id) & reads.dye_uM.eq(candidate.dye_uM)]
        pairs.append({"comparison_id": f"robust_pair_{len(pairs) + 1:04d}", "comparison": "nonUV_minus_UV",
                      "candidate_acquisition_id": candidate.acquisition_id,
                      "reference_acquisition_id": reference.acquisition_id.iloc[0] if len(reference) == 1 else "",
                      "pair_status": "matched" if len(reference) == 1 else "missing_or_ambiguous_reference",
                      "candidate_physical_plate_id": candidate.physical_plate_id,
                      "reference_physical_plate_id": reference.physical_plate_id.iloc[0] if len(reference) == 1 else "",
                      "assay": candidate.assay, "step_id": candidate.step_id, "reference_set_id": candidate.reference_set_id,
                      "same_physical_plate": False, "scope": SCOPE})
    return pd.DataFrame(pairs)


def paired_observations(observations, matches):
    output = []
    grouped = {name: group for name, group in observations.groupby("acquisition_id")}
    for pair in matches[matches.pair_status.eq("matched")].to_dict("records"):
        left, right = grouped[pair["reference_acquisition_id"]], grouped[pair["candidate_acquisition_id"]]
        joined = left.merge(right, on=["well", "condition", "signal"], suffixes=("_reference", "_candidate"), validate="one_to_one")
        for _, row in joined.iterrows():
            same_predictor = bool(np.isclose(row.predictor_reference, row.predictor_candidate, equal_nan=True))
            valid = bool(row.eligible_reference and row.eligible_candidate and same_predictor)
            delay = (pd.Timestamp(row.recorded_time_candidate) - pd.Timestamp(row.recorded_time_reference)).total_seconds() / 3600
            output.append({**pair, "condition": row.condition, "signal": row.signal, "well": row.well,
                           "level_id": row.level_id_candidate, "predictor": row.predictor_candidate,
                           "predictor_unit": row.predictor_unit_candidate, "response_reference": row.response_reference,
                           "response_candidate": row.response_candidate,
                           "response_difference": row.response_candidate - row.response_reference if valid else np.nan,
                           "pair_eligible": valid, "same_predictor": same_predictor,
                           "elapsed_hours": delay, "recorded_time_reference": row.recorded_time_reference,
                           "recorded_time_candidate": row.recorded_time_candidate,
                           "plate_type": row.plate_type_candidate, "temperature_C": row.temperature_C_candidate,
                           "acid_M": row.acid_M_candidate, "dye_uM": row.dye_uM_candidate,
                           "read_index": row.read_index_candidate, "step_order": row.step_order_candidate,
                           "independent_preparations": row.independent_preparations_candidate,
                           "reference_status": row.eligibility_status_reference, "candidate_status": row.eligibility_status_candidate})
    return pd.DataFrame(output)


def summarize_pairs(pairs):
    rows = []
    for _, group in pairs.groupby(["comparison_id", "condition", "signal"], dropna=False, sort=True):
        first = group.iloc[0]
        valid = group[group.pair_eligible]
        # Complete paired levels only; technical wells are averaged first.
        level = valid.groupby("level_id").agg(reference=("response_reference", "mean"), difference=("response_difference", "mean"))
        width = level.reference.max() - level.reference.min()
        rows.append({**{key: first[key] for key in ("comparison_id", "comparison", "candidate_acquisition_id", "reference_acquisition_id", "candidate_physical_plate_id", "reference_physical_plate_id", "assay", "step_id", "reference_set_id", "condition", "signal", "elapsed_hours", "plate_type", "temperature_C", "acid_M", "dye_uM", "read_index", "step_order", "independent_preparations")},
                     "paired_well_n": len(valid), "scheduled_pair_n": len(group), "paired_level_n": len(level),
                     "mean_response_difference": level.difference.mean(), "level_difference_sd": level.difference.std(ddof=1),
                     "reference_observed_range": width,
                     "mean_difference_fraction_of_reference_range": level.difference.mean() / width if width > 0 else np.nan,
                     "maximum_absolute_difference": valid.response_difference.abs().max(),
                     "rms_difference": np.sqrt(np.mean(valid.response_difference ** 2)),
                     "scope": SCOPE})
    return pd.DataFrame(rows)


def signal_comparisons(observations):
    rows = []
    for identity, group in observations.groupby(["acquisition_id", "condition"], dropna=False):
        first = group.iloc[0]
        comparisons = [("ratio_raw", "ratio_corrected")] if first.assay == "pH" else [("A380", "A380_corrected"), ("A428", "A428_corrected"), ("A380", "A428")]
        for reference, candidate in comparisons:
            a = group[group.signal.eq(reference)]
            b = group[group.signal.eq(candidate)]
            if a.empty or b.empty:
                continue
            joined = a.merge(b, on="well", suffixes=("_reference", "_candidate"), validate="one_to_one")
            valid = joined[joined.eligible_reference & joined.eligible_candidate]
            record = {**{key: first[key] for key in META_KEYS if key != "signal"},
                      "reference_signal": reference, "candidate_signal": candidate,
                      "scheduled_well_n": len(joined), "reference_eligible_n": int(joined.eligible_reference.sum()),
                      "candidate_eligible_n": int(joined.eligible_candidate.sum()), "common_eligible_n": len(valid),
                      "eligibility_lost_to_candidate": int((joined.eligible_reference & ~joined.eligible_candidate).sum()),
                      "denominator_guard_is_detection_limit": False}
            for role in ("reference", "candidate"):
                level = valid.groupby("level_id_" + role)["response_" + role].agg(["mean", "std", "count"])
                df = (level["count"] - 1).sum()
                sd = np.sqrt(((level["count"] - 1) * level["std"] ** 2).sum() / df) if df > 0 else np.nan
                width = level["mean"].max() - level["mean"].min()
                record.update({role + "_range_common_wells": width, role + "_pooled_sd_common_wells": sd,
                               role + "_range_to_sd_common_wells": width / sd if sd > 0 else np.nan})
            rows.append(record)
    return pd.DataFrame(rows)


def media_interval_diagnostics(levels):
    rows = []
    for _, curve in levels[levels.condition_role.ne("no_dye_matrix_control")].groupby(CURVE_KEYS):
        curve = curve.sort_values("predictor")
        for i in range(len(curve) - 1):
            a, b = curve.iloc[i], curve.iloc[i + 1]
            noise = np.sqrt(a.technical_sd ** 2 + b.technical_sd ** 2)
            change = b.response_mean - a.response_mean
            spacing = b.predictor - a.predictor
            rows.append({**{key: a[key] for key in META_KEYS}, "predictor_unit": a.predictor_unit,
                         "predictor_low": a.predictor, "predictor_high": b.predictor,
                         "observed_response_change": change, "observed_secant_sensitivity": change / spacing if spacing > 0 else np.nan,
                         "change_to_combined_well_sd": abs(change) / noise if noise > 0 else np.nan,
                         "complete_interval": bool(a.complete_level and b.complete_level),
                         "quantification_limit_established": False})
    return pd.DataFrame(rows)


def media_background(observations):
    rows = []
    wide = observations[observations.assay.eq("pH")].drop_duplicates(["acquisition_id", "well"])
    for (plate, medium, level), group in wide.groupby(["physical_plate_id", "medium_id", "level_id"]):
        for wavelength in (485, 545, 750):
            a = group[group.condition.eq("no_dye")][f"A{wavelength}"]
            b = group[group.condition.eq("with_dye")][f"A{wavelength}"]
            rows.append({"physical_plate_id": plate, "medium_id": medium, "level_id": level,
                         "wavelength_nm": wavelength, "no_dye_mean_AU": a.mean(), "with_dye_mean_AU": b.mean(),
                         "with_dye_minus_no_dye_AU": b.mean() - a.mean(), "no_dye_n": int(a.count()), "with_dye_n": int(b.count()),
                         "preparation_comparison": "200uL_undyed_vs210uL_dyed; medium_dilution_also_differs"})
    return pd.DataFrame(rows)
