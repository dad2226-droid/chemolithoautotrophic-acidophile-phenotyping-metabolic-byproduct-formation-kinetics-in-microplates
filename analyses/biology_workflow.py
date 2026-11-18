"""Calculate manuscript Figures 5–8 and their numerical validation tables."""

import json
from pathlib import Path

import numpy as np

import pandas as pd

from acidophile_data.biological import load_biological_observations
from acidophile_data.metadata import load_metadata
from acidophile_data.peak_response import plate_means
from acidophile_data.plate_comparisons import (
    CONDITION_KEYS, CONDITION_METADATA, SERIES_KEYS, SAMPLE_KEYS, CITRATE_EXPERIMENT, CROSS_FORMAT_EXPERIMENTS,
    _average_technical, selected_experiment_ids,
    compare_biotic_with_controls, derive_plate_signals,
)
from acidophile_data.publication_additional_channels import extract
from acidophile_data.iron_screen_and_sulfur_trajectories import screen, sulfur
from acidophile_data.publication_biology import _condition_fields, control_plate_means
from acidophile_data.publication_condition_repeatability import build_tables as repeatability_tables
from acidophile_data.publication_consolidation import activity_tables
from acidophile_data.publication_cross_format import _flask_registry, _json, build_correspondence, harmonize_plate_quantities
from acidophile_data.sulfur_and_treatment_comparisons import (
    citrate_wavelength_pairs, sulfur_baselines, timing_comparisons,
)
from acidophile_data.publication_ratio_support import classify_channel_ratios, summarize_ratio_yield
from acidophile_data.optical_and_treatment_contrasts import (
    effect_repeatability, fold_values, od_dose_tables,
)
from acidophile_data.publication_response_validation import (
    flask_tables, inhibition_tables, sulfur_signal_to_noise, tube_points,
)
from acidophile_data.indicator_and_tube_summaries import tube_od
from acidophile_data.publication_si_biological_results import build as biological_results
from acidophile_data.publication_si_correspondence_validation import build as correspondence_results
from acidophile_data.publication_si_robustness_validation import _sulfur_direction
from acidophile_data.publication_si_screening_validation import build as screening_results
from acidophile_data.vessel_comparisons import FLASK_METRICS, _summarize, flask_feii_depletion
from acidophile_data.vessels import load_flasks, load_tubes
from biological_sensitivity import run as sensitivity_results


def _iron_plates(root):
    registry = next((root / "data/03_biological_microplates").glob("*/AF_experiments.yaml"))
    ids = [e["experiment_id"] for e in load_metadata(registry)["experiments"] if e["readout_id"] == "iron"]
    design = CONDITION_METADATA + [
        "final_well_volume_uL", "initial_pH_status", "initial_pH_measured",
        "initial_pH_reported_plus_minus", "initial_pH_scope", "initial_pH_uncertainty_type",
        "nominal_initial_pH_json", "incubation_temperature_C", "components_json", "additive_stock_pH",
    ]
    coverage, curves = [], []
    # Each experiment is processed in full so that shared physical controls, wells,
    # acquisition ordinals and same-column matches remain intact.
    for start in range(0, len(ids), 8):
        raw = load_biological_observations(root, species="AF", experiment_ids=ids[start:start + 8])
        fields = raw[CONDITION_KEYS + [c for c in design if c in raw]].drop_duplicates()
        descriptions = pd.DataFrame([_condition_fields(row) for _, row in fields.iterrows()])
        fields = pd.concat([fields.reset_index(drop=True), descriptions], axis=1)
        support = raw[raw.kind.eq("biotic")].groupby(CONDITION_KEYS, as_index=False).agg(
            physical_culture_n=("culture_id", "nunique"), batch_n=("batch_id", "nunique"),
            physical_plate_n=("plate_id", "nunique"), physical_observation_n=("physical_observation_key", "nunique"),
            actual_first_days=("days", "min"), actual_last_days=("days", "max"),
        )
        coverage.append(fields.merge(support, on=CONDITION_KEYS, validate="one_to_one"))
        required = set(CONDITION_KEYS + SERIES_KEYS + SAMPLE_KEYS + [
            "physical_observation_key", "well", "kind", "source", "days", "ordinal", "timestamp",
            "physical_plate_id", "readout", "A380", "A380_eligible", "A750", "A750_eligible",
            "A485", "A485_eligible", "A545", "A545_eligible",
        ])
        raw = raw[[column for column in raw if column in required]]
        trajectories, _ = compare_biotic_with_controls(derive_plate_signals(raw))
        curves.append(plate_means(trajectories))
    conditions = pd.concat(coverage, ignore_index=True)
    conditions["independent_plate_n"] = conditions.physical_plate_n
    conditions["experimental_unit"] = "independent_microplate"
    trajectories = pd.concat(curves, ignore_index=True)
    return {
        "af_iron_condition_coverage": conditions,
        "af_iron_plate_trajectories": trajectories,
        "af_iron_control_plate_trajectories": control_plate_means(trajectories, conditions),
    }


def _correspondence(root, plate_tables, source_flasks):
    chosen = set(CROSS_FORMAT_EXPERIMENTS) | {CITRATE_EXPERIMENT}
    conditions = plate_tables["plate_conditions"]
    conditions = conditions[conditions.experiment_id.isin(chosen)].copy().reset_index(drop=True)
    if chosen - set(conditions.experiment_id):
        raise ValueError("Complete general and citrate/Ga plate series are required for publication")
    conditions["comparable_quantities_mM_json"] = conditions.apply(lambda row: _json(harmonize_plate_quantities(
        json.loads(row.component_final_mM_json), json.loads(row.component_basis_json))), axis=1)
    definitions = load_metadata(root / "data/03_biological_microplates/biological_definitions.yaml")
    conditions["plate_temperature_C"] = conditions.species.map(
        lambda species: definitions["species"][species]["incubation_temperature_C"])
    conditions["plate_final_volume_mL"] = conditions.preparation_id.map(
        lambda preparation: definitions["preparations"][preparation]["final_well_volume_uL"] / 1000)
    conditions["plate_mixing"] = definitions["incubation"]["mixing"]
    conditions["plate_abiotic_controls"] = "same-plate same-column B/C/D; per-acquisition eligibility is required"
    conditions["plate_concentration_determination"] = np.where(conditions.experiment_id.eq(CITRATE_EXPERIMENT),
        "nominal citrate and explicitly elemental Ga", "reported preparation concentrations; general stock-assay conventions")
    trajectory = plate_tables["plate_culture_trajectories"]
    trajectory = trajectory[trajectory.experiment_id.isin(chosen)].copy()
    windows = trajectory.groupby(CONDITION_KEYS, as_index=False).agg(
        plate_first_acquisition_d=("days_min", "min"), plate_last_acquisition_d=("days_max", "max"),
        plate_n=("plate_id", "nunique"), plate_culture_n=("plate_id", "nunique"), plate_batch_n=("batch_id", "nunique"),
        plate_culture_ids=("plate_id", lambda s: _json(sorted(s.dropna().unique()))),
        plate_batch_ids=("batch_id", lambda s: _json(sorted(s.unique()))),
        plate_ids=("physical_plate_ids", lambda s: _json(sorted({plate for value in s for plate in value.split(";")}))))
    conditions = conditions.merge(windows, on=CONDITION_KEYS, validate="one_to_one")
    registry = _flask_registry(root)
    flask_windows = source_flasks.groupby("condition_ID", as_index=False).agg(
        flask_first_acquisition_d=("elapsed_days_mean", "min"), flask_last_acquisition_d=("elapsed_days_mean", "max"))
    registry = registry.merge(flask_windows, left_on="flask_condition_ID", right_on="condition_ID", validate="one_to_one").drop(columns="condition_ID")
    correspondence = build_correspondence(conditions, registry)
    flask_course = _summarize(
        source_flasks, vessel_column="flask_ID", culture_column="inoculum_culture_id",
        time_column="elapsed_days_mean", time_sd_column="elapsed_days_sd",
        metrics={"feii_depletion_fraction": FLASK_METRICS["feii_depletion_fraction"]},
    )
    return {
        "cross_format_plate_formulations": conditions,
        "cross_format_condition_correspondence": correspondence,
        "cross_format_plate_timecourses": trajectory.merge(
            conditions[CONDITION_KEYS + ["comparable_quantities_mM_json"]], on=CONDITION_KEYS, validate="many_to_one"),
        "cross_format_flask_timecourses": flask_course,
    }


def _plate_and_vessel_inputs(root):
    ids = selected_experiment_ids(root)
    conditions, sulfur_channels, well_trajectories = [], [], []
    cols = CONDITION_KEYS + ["plate_id", "batch_id", "culture_id", "kind", "row", "well", "source",
                             "ordinal", "days", "physical_observation_key"]
    cols += [f"A{wave}{suffix}" for wave in (545, 485, 750)
             for suffix in ("", "_eligible", "_overrange", "_parse_error", "_censored_high")]
    for start in range(0, len(ids), 8):
        raw = load_biological_observations(root, experiment_ids=ids[start:start + 8])
        conditions.append(raw[CONDITION_KEYS + CONDITION_METADATA].drop_duplicates())
        sulfur_channels.append(raw.loc[raw.readout.eq("ph"), cols].copy())
        trajectories, _ = compare_biotic_with_controls(derive_plate_signals(raw))
        well_trajectories.append(trajectories)
    conditions = pd.concat(conditions, ignore_index=True)
    trajectories = pd.concat(well_trajectories, ignore_index=True)
    chosen = set(CROSS_FORMAT_EXPERIMENTS) | {CITRATE_EXPERIMENT}
    comparable = trajectories[trajectories.experiment_id.isin(chosen)]
    means = _average_technical(comparable, CONDITION_KEYS + ["batch_id", "plate_id", "signal", "ordinal"], "days")
    flask_points, _ = flask_feii_depletion(load_flasks(root))
    tables = _correspondence(root, {"plate_conditions": conditions, "plate_culture_trajectories": means}, flask_points)
    sulfur_conditions = conditions[conditions.readout.eq("ph")]
    sulfur_rows = trajectories[trajectories.experiment_id.isin(sulfur_conditions.experiment_id) &
                               trajectories.signal.isin(["raw_R545_485", "corrected_R545_485"])]
    tables["validation_sulfur_trajectory_evidence"] = sulfur_signal_to_noise(sulfur_rows)
    tables["context_sulfur_eligibility"] = sulfur_conditions
    tables["context_plate_trajectories"] = conditions
    tables["context_sulfur_channels"] = pd.concat(sulfur_channels, ignore_index=True).merge(
        conditions, on=CONDITION_KEYS, validate="many_to_one")
    lf = plate_means(trajectories[trajectories.species.eq("LF")])
    lf = lf.merge(conditions, on=CONDITION_KEYS, validate="many_to_one")
    lf["culture_id"] = lf.plate_id
    lf["experimental_unit"] = "independent_microplate"
    tables["validation_lf_plate_trajectories"] = lf
    tables["validation_tube_points"] = tube_points(load_tubes(root))
    tables.update(flask_tables(tables, {"flask_feii_timecourse": flask_points}))
    return tables


def _treatment_responses(root, tables):
    tables.update(inhibition_tables(tables))
    tables.update(repeatability_tables(tables))
    tables.update(screen(tables))
    tables.update(sulfur(tables))
    registry, curves = sulfur_baselines(tables, root / "analyses/metadata/sulfur_baseline_references.json")
    tables["sulfur_baseline_registry"] = registry
    tables["sulfur_baseline_trajectories"] = curves
    tables["tube_od600"] = tube_od(tables)
    tables.update(activity_tables(tables))
    timing = timing_comparisons(tables)
    tables["timing_comparisons"] = timing
    folds = timing.assign(readout="raw_A380")
    for label, numerator, denominator in [("response", "plate_mean", "flask_mean"),
                                          ("time", "plate_time_to_max_mean", "flask_time_to_max_mean")]:
        folds[label + "_fold"], folds[label + "_log2_fold"], folds[label + "_fold_status"] = (
            fold_values(folds[numerator], folds[denominator]))
    tables["format_folds"] = folds
    tables.update(extract(root))
    tables.update(od_dose_tables(tables))
    tables["citrate_wavelength_pairs"] = citrate_wavelength_pairs(tables)
    tables.update(effect_repeatability(tables))


def _sulfur_availability(tables):
    classified = classify_channel_ratios(tables["context_sulfur_channels"])
    classified["unit"] = classified.batch_id.astype(str) + "|" + classified.plate_id.astype(str)
    summary, plates = [], []
    for (experiment, kind), group in classified.groupby(["experiment_id", "kind"]):
        row, physical = summarize_ratio_yield(group)
        keys = {"experiment_id": experiment, "kind": kind}
        summary.append(keys | row)
        plates.append(physical.assign(**keys))
    return {"sulfur_read_availability": pd.DataFrame(summary),
            "sulfur_read_availability_plates": pd.concat(plates, ignore_index=True)}


# Physical-unit curves and contrasts required to redraw the manuscript panels.
FIGURE_TABLES = {
    "af_iron_condition_coverage": ("5", "Iron plate formulations, dose bases, controls and physical coverage."),
    "af_iron_plate_trajectories": ("5–7", "AF iron changes by acquisition, formulation and physical plate, including technical-well support."),
    "af_iron_control_plate_trajectories": ("5a, 7g", "No-addition AF iron control changes by starting batch and physical plate."),
    "condition_repeatability_batches": ("5c–d", "Observed response and timing means by standardized condition and independent starting batch."),
    "validation_plate_peak_metrics": ("5–7", "Observed optical maxima, corresponding times and integrated areas for individual plate conditions."),
    "validation_flask_points": ("5b, 6b–d, 7a, 7e–f", "Individual flask acquisitions, including eligible iron measurements, measurable depletion and calculated apparent planktonic OD."),
    "od_plate_trajectories": ("6a, 6c–d", "Observed gallium plate time courses for A600, A600 minus A750 and A750."),
    "od_dose_units": ("6c–d", "Matched gallium response and timing for each microplate or flask."),
    "format_folds": ("7e–f", "Matched cross-format response and time ratios for the iron product readout."),
    "consolidated_activity_replicates": ("7e–f", "AUC ratios, control identities and timing for matched physical units in both formats."),
    "cross_format_condition_correspondence": ("7e–f", "Factual plate-to-flask dose matching and recorded experimental differences."),
    "cross_format_plate_timecourses": ("7b–c", "Per-plate biotic and matched abiotic signal, change and net change at every recorded acquisition, including citrate series."),
    "citrate_wavelength_pairs": ("7d", "Within-unit paired citrate detector readings with readout eligibility and acquisition identity."),
    "validation_lf_plate_trajectories": ("7g", "Observed LF iron time courses by formulation and independent microplate."),
    "sulfur_baseline_registry": ("8a–f", "Explicit soluble-sulfur reference condition for each displayed cohort."),
    "sulfur_baseline_trajectories": ("8a–f", "Measured and scatter-corrected ratio time courses for the registered soluble-sulfur formulations and their matched abiotic wells."),
    "tube_od600": ("7h–i, 8g–h", "Tube pH, ORP, change from shared initial reference and apparent planktonic OD for the fresh starter cohort."),
}

SUPPORT = {
    "treatment_effect_repeatability": ("5", "Across-batch repeatability of matched treatment responses and timing."),
    "treatment_effect_physical_units": ("5", "Matched effect estimates by physical plate or flask."),
    "treatment_effect_batches": ("5", "Batch means and rank support for matched effects."),
    "treatment_effect_condition_results": ("5", "Condition-level treatment contrast repeatability."),
    "treatment_rank_correlations": ("5", "Batch-pair rank concordance on shared treatments."),
    "acquisition_sensitivity": ("5–7", "Sensitivity of observed response and timing to removing the final acquisition, retaining the first 75% of the window, retaining alternate acquisitions or withholding the middle acquisition."),
    "acquisition_sensitivity_trajectories": ("5–7", "Full and thinned observation grids within each physical-unit trajectory."),
    "citrate_endpoint_effects": ("7", "Matched citrate endpoint effects by readout and format."),
    "citrate_endpoint_pairs": ("7", "Citrate condition and reference endpoint pairs used for each contrast."),
    "sulfur_read_availability": ("8", "Equal-plate availability of paired raw and corrected ratios by sulfur experiment and well role."),
    "sulfur_read_availability_plates": ("8", "Per-plate paired acquisition counts, valid fractions and corrected-denominator retention."),
}


def _results(tables):
    outputs, guide = {}, []
    specs = [*biological_results(tables), *screening_results(tables),
             *correspondence_results(tables), _sulfur_direction(tables)]
    for spec in specs:
        key = spec["key"]
        outputs[key] = spec["frame"]
        guide.append({"table": key, "title": spec["title"], "caption": spec["caption"],
                      "panels": ", ".join(spec["figure_panels"]), "source_tables": ", ".join(spec["source_tables"])})
        for name, frame in spec.get("support_tables", {}).items():
            name = key + "__" + name
            outputs[name] = frame
            guide.append({"table": name, "title": name.replace("_", " "), "caption": spec["caption"],
                          "panels": ", ".join(spec["figure_panels"]), "source_tables": ", ".join(spec["source_tables"])})
    extra = {key: tables[key] for key in FIGURE_TABLES}
    extra.update(sensitivity_results(tables))
    extra.update(_sulfur_availability(tables))
    for key, frame in extra.items():
        panels, caption = (FIGURE_TABLES | SUPPORT)[key]
        outputs[key] = frame
        guide.append({"table": key, "title": key.replace("_", " "), "caption": caption, "panels": panels, "source_tables": ""})
    outputs["result_guide"] = pd.DataFrame(guide)
    return outputs


def run(root: Path) -> dict[str, pd.DataFrame]:
    root = Path(root)
    tables = _iron_plates(root)
    tables.update(_plate_and_vessel_inputs(root))
    _treatment_responses(root, tables)
    return _results(tables)
