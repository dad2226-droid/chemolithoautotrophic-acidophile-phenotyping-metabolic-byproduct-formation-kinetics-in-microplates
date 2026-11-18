"""Export the numerical results supporting the manuscript's Technical Validation."""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from html import escape
import json
from pathlib import Path
import re
from typing import Mapping

import pandas as pd


@dataclass(frozen=True)
class Table:
    title: str
    panels: str
    claim: str
    sources: str
    fields: str = ""


FIT = ("plate_n plate_mean_n level_n fit_rmse r_squared identified converged identification_flags "
       "ci_df bottom bottom_ci_low bottom_ci_high span span_ci_low span_ci_high midpoint midpoint_ci_low "
       "midpoint_ci_high slope slope_ci_low slope_ci_high")
TEST = "plate_n estimate sd ci_low ci_high p_value p_holm multiplicity_family_n status"
DOSE = ("additive plate_dose_mM flask_dose_mM Ga_mM match_status comparison_context_n defined_context_n "
        "plate_batch_n flask_batch_n plate_vessel_n flask_vessel_n median_log2_fold q25_log2_fold "
        "q75_log2_fold median_fold minimum_log2_fold maximum_log2_fold final_read_maximum_context_n status")

TABLES = (
    Table("pH indicator calibration", "1b–d",
          "Four-parameter logistic fits describe the measured and scatter-corrected ratio response where replicate groups are complete. Bottom and span are in ratio units, midpoint is pH and slope is per pH unit.",
          "ph_response_fit", "dye_uM readout "+FIT+" unavailable_plate_mean_n"),
    Table("Scatter correction and pH-indicator precision", "1c–d",
          "At 10 µM thymol blue, correction increases response span and within-standard SD. The span-to-SD contrast has Holm-adjusted p = 0.0505. The estimates and their intervals use log₂(corrected/measured). The geometric ratios give the corresponding folds.",
          "initial_correction_precision", "metric raw_mean corrected_mean "+TEST+" geometric_ratio ratio_ci_low ratio_ci_high"),
    Table("pH prediction in acid standards and ferric mixtures", "1b–d",
          "Scatter correction reduces withheld-acid prediction RMSE. Both readouts have higher error in ferric mixtures over the common pH support. The paired SD comparison is a ratio. Paired bias and RMSE contrasts are corrected minus measured.",
          "ph_transfer_accuracy ph_transfer_correction"),
    Table("Acid-specific ferric calibrations", "1f",
          "The fitted A380 response approaches a plateau. The acid-specific model gives baseline and span in AU, half-response concentration in mM total Fe and a dimensionless shape parameter.",
          "ferric_response_fit", "acid_M "+FIT),
    Table("Cerimetry against measured total iron", "1g",
          "The ferrous-standard slopes are below one in 0.001 and 1 M sulfuric acid. The ferric-standard series includes a constant cerimetry response. Tests compare the slope with zero and one and the intercept with zero.",
          "chemical_reference_regression", "series_ID n x_min_M x_max_M slope slope_ci_low slope_ci_high intercept_M intercept_ci_low_M intercept_ci_high_M sy_x_M residual_df slope_vs_zero_p slope_zero_p_holm slope_vs_one_p slope_one_p_holm intercept_vs_zero_p intercept_zero_p_holm"),
    Table("Ferric–ferrous separation at 380 and 428 nm", "1e",
          "The ferric–ferrous contrast in AU is larger at 380 nm than at 428 nm across the same supported concentrations. The wavelength-difference row retains within-plate wavelength pairing.",
          "ferric_wavelength_selectivity", "readout ferric_plate_n ferrous_plate_n common_level_n reference_min_M reference_max_M estimate ci_low ci_high welch_df"),
    Table("Dye and scatter correction in the ferric standards", "1h",
          "With-versus-without-dye comparisons give Holm-adjusted p = 1.00 for both A380 and scatter-corrected A380. The four planned contrasts use matched standard positions within physical plates and report effects in AU.",
          "iron_dye_correction", "contrast context "+TEST),
    Table("Lid effects on indicator response", "2a, 3a",
          "With the lid, measured and corrected pH ratios are lower. The corresponding ferric indicator contrasts have Holm-adjusted p = 0.140 and 0.193. Effects are lid-present minus lid-absent, expressed as a percentage of the paired reference span.",
          "lid_effect", "assay dye readout paired_well_n paired_level_n_min paired_level_n_max "+TEST),
    Table("Indicator drift during benchtop storage", "2b–c, 3b–e",
          "Mean final ferric-indicator changes are small relative to initial spans, while some individual concentrations show larger drift. Deviations also occur among standards measured throughout storage. The interval type distinguishes confidence limits from observed ranges.",
          "storage_drift storage_trajectory_drift storage_trajectory_drift__common_support_drift_summary"),
    Table("Edge–interior differences in storage drift", "2d, 3f",
          "None of the 12 endpoint and time-averaged edge–interior contrasts reaches significance after Holm adjustment. Effects are percentages of the initial observed response span.",
          "edge_drift", "assay dye readout metric mean_endpoint_d maximum_absolute_observed_percent "+TEST),
    Table("Medium-specific indicator response and precision", "4a–d",
          "Response span, within-standard SD, their ratio and adjacent-standard reversals vary by medium. AFM4 has too few corrected standards with replicates to estimate a span in the primary comparison.",
          "media_response_performance", "medium readout plate_n eligible_well_n scheduled_well_n span_plate_n span_mean span_ci_low span_ci_high sd_plate_n technical_sd_mean sd_ci_low sd_ci_high span_sd_mean median_rank_rho plates_with_adjacent_reversal"),
    Table("Paired scatter-correction effects in growth media", "4d",
          "Across seven estimable comparisons, Holm-adjusted p ≥ 0.891 for changes in span divided by within-standard SD. Each comparison uses the same eligible wells for both readouts. The effect is log₂(corrected/measured).",
          "media_correction_effects", "medium readout scheduled_plate_n unavailable_plate_n common_well_n raw_Q_mean corrected_Q_mean "+TEST),
    Table("Between-batch variance in plates and flasks", "5c–d",
          "Both plate–flask comparisons give Holm-adjusted p = 0.617. All four estimable paired flask comparisons give p = 1.00 after adjustment. The independent-format and paired-flask tests have separate families of two and six planned comparisons, respectively.",
          "bio_baseline_repeatability bio_flask_repeatability"),
    Table("Batch repeatability of integrated response and timing", "5a–d",
          "The batch distributions quantify baseline consistency and show how AUC and time to first observed maximum vary among independent batches. Native AUC units are specific to each readout.",
          "bio_repeatability_distribution", "metric cohort native_unit batch_n vessel_n native_mean native_sd cv_percent log_outcome_sd centered_deviation_iqr_pp final_read_maximum_vessels status"),
    Table("Gallium response across turbidity and fluorescence readouts", "6a–d",
          "The same-batch AUC and timing responses differ between plate A600, corrected A600 and fluorescence-derived apparent planktonic OD600. The AUC effect is percentage points relative to the untreated response, with zero denoting no effect. Timing is in days.",
          "od_ga_response od_ga_timing", "format readout Ga_mM unit candidate_batch_n batch_n physical_unit_n effect between_batch_sd within_batch_unit_sd ci_low ci_high final_read_maximum_n status"),
    Table("Net 380 versus 428 nm response at citrate endpoints", "7d",
          "None of 18 citrate–gallium formulations shows a significant net wavelength difference after Holm adjustment. The paired contrast is net A380 minus net A428 in AU.",
          "citrate_wavelength_contrast", "citric_acid_mM Ga_mM paired_plate_n batch_n net_380_mean_AU net_428_mean_AU paired_380_minus_428_AU between_batch_difference_sd_AU ci_low ci_high statistic df p_value p_holm days_min days_max status"),
    Table("Plate–flask correspondence in AUC and timing", "7e–f",
          "Correspondence varies by treatment for both within-format control-normalized AUC and time to first observed maximum. The folds compare the plate quantity with its flask counterpart. The match status distinguishes exact and neighboring doses.",
          "product_format_response product_format_timing", DOSE),
    Table("L. ferriphilum plate response and time-matched A. ferrooxidans benchmark", "7g",
          "L. ferriphilum A380 remains close to the matched abiotic response, while the independent A. ferrooxidans benchmark shows a larger change at the matched duration.",
          "lf_plate_benchmark", "LF_plate_context LF_Ga_mM LF_In_mM LF_batch_n LF_plate_n LF_inoculum_culture_n LF_biotic_mean_AU LF_abiotic_mean_AU LF_net_change_AU LF_net_between_microplate_sd_AU LF_endpoint_min_d LF_endpoint_max_d AF_reference_day_d AF_reference_mean_AU AF_reference_batch_sd_AU AF_batch_n AF_plate_n LF_biotic_minus_AF_reference_AU status"),
    Table("Fresh-starter culture-tube endpoint contrasts", "7h–i, 8g–h",
          "Fresh L. ferriphilum starters give positive ORP contrasts under some conditions and higher apparent OD600. A. ferrooxidans and A. thiooxidans cultures supplied with 5 mM tetrathionate have lower final pH and higher apparent OD600 than their abiotic controls.",
          "tube_context_contrasts sulfur_tube_contrasts", "organism Ga_mM In_mM metric unit endpoint_d independent_batch_n biotic_starter_n biotic_n abiotic_n biotic_mean biotic_sd abiotic_mean abiotic_sd biotic_minus_abiotic status"),
    Table("Soluble-sulfur endpoint ratio contrasts", "8a–f",
          "Measured biotic-minus-abiotic ratio changes are small or variable across the tested thiosulfate and tetrathionate cultures. Counts distinguish all candidate plates from plates with eligible scheduled endpoints.",
          "sulfur_baseline_contrast", "cohort readout candidate_plate_n paired_endpoint_plate_n batch_n inoculum_culture_n net_endpoint_change net_between_microplate_sd minimum_net_change maximum_net_change endpoint_min_d endpoint_max_d marginal_eligible_endpoint_plate_n common_eligible_endpoint_plate_n status"),
    Table("Scatter correction and sulfur endpoint separation", "8a–f",
          "Scatter correction does not consistently increase endpoint separation scaled by abiotic SD. Signed and absolute score changes use the paired plate–formulation observations.",
          "sulfur_paired_correction", "species medium candidate_plate_formulation_pairs paired_n unavailable_pair_n batch_n inoculum_culture_n median_signed_snr_change signed_snr_q25 signed_snr_q75 median_magnitude_change magnitude_q25 magnitude_q75 magnitude_increased_n magnitude_increased_percent direction_agreement_percent status"),
    Table("Apparent iron-oxidation screen outcomes", "5a",
          "Most no-additive A. ferrooxidans controls satisfy the simultaneous screening criteria. Some positive screens have net biotic changes below 1 AU. The evaluable denominator includes positive events and complete negative observations.",
          "screen_decision_performance", "species stratum candidate_plate_condition_n physical_plate_n starter_batch_n evaluable_decision_n unevaluable_decision_n event_n event_percent no_event_n partial_no_event_n first_event_net_min_AU first_event_net_q25_AU first_event_net_median_AU first_event_net_q75_AU event_without_net_1AU_n event_without_net_1AU_percent observed_window_median_d status"),
    Table("Within- and between-batch agreement of iron screens", "5a",
          "Event fractions are generally consistent across batches. Counts distinguish agreement, discordance and insufficient information at each comparison level.",
          "screen_decision_repeatability", "species stratum comparison_level candidate_comparison_n evaluable_comparison_n incomplete_comparison_n identical_event_fraction_n identical_event_fraction_percent event_fraction_range_median_pp event_fraction_range_q25_pp event_fraction_range_q75_pp event_fraction_range_max_pp all_positive_comparison_n all_negative_comparison_n identical_mixed_comparison_n discordant_comparison_n disagreement_percent status"),
    Table("Inverse ferric concentration predictions", "1f",
          "Near saturation, inverse concentration predictions show larger errors or cannot be evaluated. Rows report the attempted and successful predictions for each withheld reference window.",
          "ferric_inverse_validation", "acid_M validation nominal_reference_window_mM reference_min_mM reference_max_mM plate_n attempted_n predicted_n failed_n outside_reference_support_n response_or_inverse_failure_n converged_unidentified_fit_n bias_mM rmse_mM median_absolute_error_mM median_absolute_relative_error_percent relative_error_n plate_rmse_min_mM plate_rmse_max_mM"),
    Table("Withheld prediction by shared and acid-specific ferric models", "1f",
          "Acid-specific fits do not consistently reduce prediction error at withheld microplates or concentrations. RMSE and bias describe absorbance in AU.",
          "ferric_predictive_model_comparison", "acid_M validation plate_n attempted_n paired_prediction_n shared_rmse_AU shared_bias_AU acid_specific_rmse_AU acid_specific_bias_AU rmse_change_AU squared_error_reduction_percent"),
    Table("Separation between adjacent indicator standards", "1b–d, 1f, 4a–c",
          "Local separation varies across the tested standards, with reversed ordering in some conditions. Each interval reports its response change, slope per reference unit and response change divided by local measurement noise.",
          "local_response_resolution", "family cohort readout assay reference_unit interval_rank scheduled_plate_n reference_low reference_high response_change_n response_change response_change_ci_low response_change_ci_high secant_sensitivity_n secant_sensitivity secant_sensitivity_ci_low secant_sensitivity_ci_high change_over_noise_n change_over_noise change_over_noise_ci_low change_over_noise_ci_high opposite_direction_plate_n interval_status"),
    Table("Response span, precision and separation across handling conditions", "2e–h, 3g–l",
          "A similar span across temperature, plate and membrane conditions can coexist with a different pooled within-standard SD or signed adjacent-standard separation. Folds summarize span and SD. Local separation is a signed score difference.",
          "handling_discrimination", "assay readout dye material contrast step_id acid_M paired_plate_n span_plate_n span_fold_estimate span_fold_ci_low span_fold_ci_high sd_plate_n sd_fold_estimate sd_fold_ci_low sd_fold_ci_high local_plate_n local_estimate local_ci_low local_ci_high reference_reversal_plate_n candidate_reversal_plate_n span_status sd_status local_status"),
    Table("Eligibility of measured and scatter-corrected ratios", "2b–c, 4a–b, 8a–f",
          "Corrected ratios fail eligibility more often in some optical configurations and A. ferrooxidans sulfur treatments. Counts include each acquired well and read. Percentages are averaged across physical plates with equal weight.",
          "optical_ratio_availability sulfur_read_availability", "experiment_id plate_n attempt_n raw_valid_n corrected_valid_n both_valid_n raw_only_n corrected_only_n neither_valid_n raw_percent_mean corrected_percent_mean difference_pp_mean difference_pp_min difference_pp_max"),
    Table("Agreement with flask chemistry over common observation windows", "7b, 7e–f",
          "Scatter correction does not consistently reduce the error between the plate and chemical trajectories scaled to their changes over the shared interval. Differences give corrected minus measured error.",
          "common_window_readout_error", "plate_dose flask_dose dose_match candidate_context_n paired_context_n unavailable_context_n cohort_pair_n plate_batch_n flask_batch_n raw_error corrected_error paired_error_difference error_improved_context_percent window_min_d window_max_d status"),
    Table("AUC comparisons over a common observation window", "7e–f",
          "Using the common observation interval changes plate and flask control-relative AUC estimates. Changes are percentage points relative to the corresponding time-to-maximum calculation.",
          "common_window_response_sensitivity", "additive plate_dose_mM flask_dose_mM Ga_mM match_status candidate_context_n paired_context_n unavailable_context_n plate_batch_n flask_batch_n plate_response_change_pp flask_response_change_pp format_difference_change_pp primary_format_difference_pp fixed_format_difference_pp window_min_d window_max_d status"),
    Table("Matched citrate contrasts in plate and flask endpoints", "7a–c",
          "Citrate-dependent changes occur in the abiotic plate signal as well as the plate net response and apparent ferric generation in flasks. Each contrast uses zero citrate at the same gallium concentration. Plate endpoints are in AU and flask endpoints in mM.",
          "citrate_endpoint_effects", "readout citric_acid_mM reference_citric_acid_mM Ga_mM candidate_batch_n paired_batch_n candidate_unit_n paired_unit_n unavailable_unit_n mean_effect between_batch_effect_sd within_batch_pair_sd ci_low ci_high endpoint_min_d endpoint_max_d status"),
    Table("Direction of sulfur indicator changes after scatter correction", "8a–f",
          "The sign of the net sulfur indicator change sometimes reverses after scatter correction. Counts identify the common observations and the nonzero comparisons used to assess sign changes.",
          "sulfur_correction_direction", "experiment_id cohort plate_n postinitial_condition_acquisition_n common_valid_n nonzero_common_n opposite_sign_n opposite_sign_percent opposite_sign_percent_min opposite_sign_percent_max raw_comparison_yield_percent corrected_comparison_yield_percent"),
)

# The shared schemas give one measurement or statistical contrast per row.
MEASURES = {
    "od_ga_response": "AUC relative to same-batch control",
    "od_ga_timing": "time to maximum minus same-batch control",
    "product_format_response": "control-relative AUC",
    "product_format_timing": "time to first observed maximum",
    "tube_context_contrasts": "L. ferriphilum fresh-starter tubes",
    "sulfur_tube_contrasts": "5 mM tetrathionate tubes",
}


def _select(frame, columns, rename=None):
    result = frame.loc[:, columns.split()].copy()
    if rename:
        result.rename(columns=rename, inplace=True)
    return result


def _transfer(tables):
    rows = []
    groups = (
        ("ph_transfer_accuracy", "readout", "attempted_n", "finite_prediction_n",
         (("bias", "pH", "bias_pH", "bias"), ("RMSE", "pH", "rmse_pH", "rmse"),
          ("technical SD", "pH", "technical_sd_pH", "technical_sd"))),
        ("ph_transfer_correction", "comparison", "attempted_pair_n", "paired_n",
         (("paired bias difference", "pH", "bias_difference_pH", "bias_difference"),
          ("paired RMSE difference", "pH", "rmse_difference_pH", "rmse_difference"),
          ("paired SD ratio", "fold", "technical_sd_ratio", "technical_sd_ratio"))),
    )
    for key, readout, attempted, available, metrics in groups:
        data = tables[key]
        for metric, unit, estimate, interval in metrics:
            values = _select(data,
                f"acid_m matrix_display {readout} {attempted} {available} support_reference_pH_min "
                f"support_reference_pH_max support_shared_bin_n {estimate} {interval}_ci_low "
                f"{interval}_ci_high bootstrap_requested_n {estimate}_bootstrap_valid_n",
                {"acid_m": "acid_M", "matrix_display": "matrix", readout: "readout_or_comparison",
                 attempted: "attempted_n", available: "eligible_n", estimate: "estimate",
                 interval+"_ci_low": "ci_low", interval+"_ci_high": "ci_high",
                 estimate+"_bootstrap_valid_n": "bootstrap_success_n"})
            values.insert(3, "metric", metric)
            values.insert(4, "unit", unit)
            rows.append(values)
    return pd.concat(rows, ignore_index=True)


def _storage(tables):
    """Report final drift, drift across acquisitions and fixed-standard sensitivity."""
    final = _select(tables["storage_drift"],
        "assay dye readout metric mean_endpoint_d plate_n estimate sd ci_low ci_high p_value p_holm multiplicity_family_n status",
        {"mean_endpoint_d": "endpoint_day", "ci_low": "interval_low", "ci_high": "interval_high"})
    final["support_basis"] = "available final values"
    final["interval_type"] = "95% Student-t CI"
    trajectory = tables["storage_trajectory_drift"]
    pieces = [final]
    for field in ("mean_observed_drift", "max_absolute_plate_drift", "max_absolute_level_drift", "rms_observed_level_drift"):
        summary = _select(trajectory,
            f"assay dye readout last_recorded_day_mean plate_n {field}_mean {field}_min {field}_max",
            {"last_recorded_day_mean": "endpoint_day", field+"_mean": "estimate",
             field+"_min": "interval_low", field+"_max": "interval_high"})
        summary["support_basis"] = "available level and acquisition values"
        summary["metric"] = field+" (% initial span)"
        summary["interval_type"] = "observed plate range"
        pieces.append(summary)
    sensitivity = _select(tables["storage_trajectory_drift__common_support_drift_summary"],
        "assay condition signal metric support_basis endpoint_day_mean eligible_plate_n "
        "estimate interval_type interval_low interval_high fixed_common_level_n status",
        {"signal":"readout", "endpoint_day_mean":"endpoint_day", "eligible_plate_n":"plate_n"})
    sensitivity["dye"] = sensitivity.pop("condition").map({"no_thymol_blue":"0 µM", "thymol_blue":"10 µM"})
    sensitivity["readout"] = sensitivity["readout"].replace({
        "ratio_raw":"R", "ratio_corrected":"R*", "A380_corrected":"A*380"})
    sensitivity["metric"] += " (% initial span)"
    pieces.append(sensitivity)
    return pd.concat(pieces, ignore_index=True)


def _variance(tables):
    common = "metric comparison method log_outcome_variance_ratio ci_low ci_high p_value p_holm status"
    plate = _select(tables["bio_baseline_repeatability"], common
                    +" plate_batch_n flask_batch_n numerator_df denominator_df",
                    {"plate_batch_n":"numerator_batch_n", "flask_batch_n":"denominator_batch_n"})
    plate["family_planned_n"] = 2
    flask = _select(tables["bio_flask_repeatability"], common
                    +" numerator_matched_batch_n denominator_matched_batch_n matched_batch_n statistic df multiplicity_family_size",
                    {"numerator_matched_batch_n":"numerator_batch_n", "denominator_matched_batch_n":"denominator_batch_n",
                     "multiplicity_family_size":"family_planned_n", "df":"paired_test_df"})
    return pd.concat([plate, flask], ignore_index=True)


def _build(table, inputs):
    keys = table.sources.split()
    if keys[0] == "ph_transfer_accuracy":
        return _transfer(inputs)
    if keys[0] == "storage_drift":
        return _storage(inputs)
    if keys[0] == "bio_baseline_repeatability":
        return _variance(inputs)
    pieces = []
    for key in keys:
        result = _select(inputs[key], table.fields)
        if key in MEASURES:
            result.insert(0, "comparison_quantity" if "tube" not in key else "cohort", MEASURES[key])
        if key == "optical_ratio_availability":
            source = inputs[key]
            result.insert(0, "cohort", "optical standards")
            result.insert(2, "plate_type", source["plate_type"])
            result.insert(3, "medium_id", source["medium_id"])
        if key == "sulfur_read_availability":
            result.insert(0, "cohort", "sulfur microplate cultures")
            result.insert(2, "culture_type", inputs[key]["kind"])
        pieces.append(result)
    return pd.concat(pieces, ignore_index=True)


def _input_frames(results):
    """Load only the necessary exports, preserving each source decimal as written."""
    from hashlib import sha256

    keys = {k for table in TABLES for k in table.sources.split()}
    frames = {}
    for section in ("optical", "biological"):
        directory = results / section
        with (directory / "table_index.csv").open(encoding="utf-8", newline="") as stream:
            for row in csv.DictReader(stream):
                key = row["table"]
                if key not in keys:
                    continue
                if row["file"] != key+".csv" or key in frames:
                    raise ValueError(f"Invalid indexed entry for {key}")
                path = directory / row["file"]
                if sha256(path.read_bytes()).hexdigest() != row["sha256"]:
                    raise ValueError(f"Checksum mismatch for {path}")
                frames[key] = pd.read_csv(path, dtype=str, keep_default_na=False)
    if missing := keys - frames.keys():
        raise ValueError(f"Analysis outputs missing {', '.join(sorted(missing))}")
    return frames


def _markdown_cell(value):
    # Preserve scientific identifiers and punctuation as literal text in GitHub.
    value = escape(value, quote=False)
    for char in ("\\", "|", "*", "[", "]", "`", "~"):
        value = value.replace(char, f"&#{ord(char)};")
    return value.replace("_", "&#95;<wbr>").replace("\r\n", "\n").replace("\n", "<br>")


def _write_captioned_table(path, number, table, definitions, frame):
    columns = list(frame.columns)
    rows = list(frame.fillna("").astype(str).itertuples(index=False, name=None))
    # Numeric columns align at the right edge. Blank estimates retain their source cells.
    decimal = re.compile(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?")
    separators = []
    for index in range(len(columns)):
        values = [r[index] for r in rows if r[index]]
        separators.append("---:" if values and all(decimal.fullmatch(v) for v in values) else "---")
    lines = [
        f"# Supplementary Table S{number}. {table.title}", "",
        table.claim, "",
        f"**Figure panels** {table.panels} · [All Supplementary tables](README.md)", "",
        "**Column definitions**", "",
        *(f"- `{column}` — " + escape(definitions[column], quote=False).replace("*", "&#42;")
          for column in columns),
        "",
    ]
    if any("" in row for row in rows):
        lines.extend(("Blank numeric cells indicate quantities that were unavailable or inapplicable for that row.", ""))
    lines.extend([
        "| " + " | ".join(_markdown_cell(c) for c in columns) + " |",
        "| " + " | ".join(separators) + " |",
    ])
    lines.extend("| " + " | ".join(_markdown_cell(v) for v in row) + " |" for row in rows)
    path.write_text("\n".join(lines)+"\n", encoding="utf-8")


def write_tables(tables: Path | str | Mapping[str, pd.DataFrame], output_dir: Path):
    """Write the full-precision CSVs, captioned Supplementary tables and manuscript claim guide."""
    inputs = _input_frames(Path(tables)) if isinstance(tables, (str, Path)) else tables
    description_file = Path(__file__).with_name("supplementary_column_descriptions.json")
    definitions = json.loads(description_file.read_text(encoding="utf-8"))
    results = [_build(table, inputs) for table in TABLES]
    if set(definitions) != {str(i) for i in range(1, len(TABLES) + 1)}:
        raise ValueError("Supplementary column descriptions must cover each table")
    for number, result in enumerate(results, 1):
        columns = definitions[str(number)]
        if result.empty or not result.columns.is_unique:
            raise ValueError(f"Invalid Supplementary table S{number}")
        if (set(result.columns) != set(columns)
                or any(not isinstance(value, str) or not value.strip() for value in columns.values())):
            raise ValueError(f"Column descriptions must match every header in Supplementary Table S{number}")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    guide = [
        "# Supplementary tables", "",
        "Each table caption defines every column in the order it appears. Select a table title "
        "to read the caption and the pipe-delimited results in one Markdown file. "
        "The figure column identifies the associated panels, including nearby panels when the "
        "statistic is reported in the surrounding text. Table rows retain the precision of the analysis exports.",
        "",
        "A blank numeric cell indicates that a quantity was unavailable or inapplicable. Counts "
        "identify the contributing physical plates, vessels, or independent batches. CI columns "
        "are 95% confidence limits unless the table identifies an observed range. The columns "
        "`p_value` and `p_holm` give the unadjusted and Holm-adjusted p-values, respectively. "
        "Percentages and percentage-point changes are distinguished by the units or metric. "
        "In cultures, A380 and the thymol blue pH ratio are indicators of the metabolic byproducts Fe³⁺ and H⁺, respectively. "
        "R denotes the measured pH indicator ratio and R* denotes its scatter-corrected counterpart. "
        "An asterisk on an absorbance readout denotes subtraction of A750.",
        "",
        "| Supplementary table | Figure panels | Technical Validation claim and interpretation |",
        "|---|---|---|",
    ]
    for number, (table, result) in enumerate(zip(TABLES, results), 1):
        path = output_dir / f"Supplementary_Table_S{number}.md"
        _write_captioned_table(path, number, table, definitions[str(number)], result)
        guide.append(f"| [S{number}. {table.title}]({path.name}) | {table.panels} | {table.claim} |")
    guide.extend(("", "## Analysis outputs", "",
                  "The source names below refer to the CSV files written by "
                  "[the analysis workflow](../../analyses/README.md).",
                  "", "| Supplementary table | Analysis CSVs |", "|---|---|"))
    for number, table in enumerate(TABLES, 1):
        sources = ", ".join(f"`{key}.csv`" for key in table.sources.split())
        guide.append(f"| S{number} | {sources} |")
    (output_dir / "README.md").write_text("\n".join(guide)+"\n", encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", required=True, type=Path, help="Completed optical and biological analysis run")
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent / "supplementary_tables")
    args = parser.parse_args()
    write_tables(args.results, args.output)
