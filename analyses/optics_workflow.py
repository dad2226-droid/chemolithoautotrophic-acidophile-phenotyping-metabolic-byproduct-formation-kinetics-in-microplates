"""Calculate manuscript optical calibrations and handling, storage, and medium results."""
from pathlib import Path

import numpy as np
import pandas as pd

from acidophile_data.optical import load_optical_registry
from acidophile_data.publication_calibration import (
    IRON_SIGNALS, _chemical_references, _core, _level_means, _read_sources, signal_value,
)
from acidophile_data.publication_calibration_reference import analyze_ph_matrix
from acidophile_data.publication_calibration_reference_adapter import canonical_ph_adapter
from acidophile_data.optical_and_treatment_contrasts import optical_folds
from acidophile_data.publication_ratio_support import pair_canonical_ratios, summarize_ratio_yield
from acidophile_data.indicator_and_tube_summaries import (
    calibration, iron_spectra, reference_regression, spectra, storage,
)
from acidophile_data.publication_robustness import (
    CURVE_KEYS, META_KEYS, build_matches, build_media_observations, build_robustness_observations,
    media_background, media_interval_diagnostics, paired_observations, signal_comparisons,
    summarize_levels, summarize_pairs,
)
from acidophile_data.publication_si_calibration_validation import build as validate_calibrations
from acidophile_data.publication_si_optical_results import (
    MEDIUM, _holm, _summary, build as optical_results,
)
from acidophile_data.publication_si_robustness_validation import (
    _handling as handling_discrimination, _storage as storage_sensitivity,
)


def _precision_pairs(levels):
    """Compare the same complete standards before and after scatter correction."""
    core = _core(levels)
    core = core[core.assay.eq("pH") & core.nominal_dye_uM.eq(10)]
    rows = []
    for (plate, condition), group in core.groupby(["physical_plate_id", "condition"]):
        pairs = group.pivot(index="condition_id", columns="signal", values=["mean", "technical_sd"])
        pairs = pairs.dropna()
        if pairs.empty:
            continue
        span = pairs["mean"].max() - pairs["mean"].min()
        # Every included standard has four eligible technical wells.
        sd = np.sqrt(pairs["technical_sd"].pow(2).mean())
        metrics = {"response_span": span, "pooled_technical_sd": sd,
                   "range_over_technical_sd": span / sd.where(sd.gt(0))}
        for metric, values in metrics.items():
            rows.append(dict(physical_plate_id=plate, condition=condition, assay="pH",
                             paired_level_n=len(pairs), metric=metric, raw_value=values.ratio_raw,
                             corrected_value=values.ratio_corrected))
    return pd.DataFrame(rows)


def _media_precision(levels):
    rows = []
    relevant = levels[levels.condition_role.ne("no_dye_matrix_control")]
    for _, group in relevant.groupby(CURVE_KEYS, sort=True, dropna=False):
        first = group.iloc[0]
        curve = group[group.technical_n.ge(2) & np.isfinite(group.response_mean)
                      & np.isfinite(group.predictor)].sort_values("predictor")
        weights = curve.technical_n - 1
        sd = np.sqrt((weights * curve.technical_sd.pow(2)).sum() / weights.sum()) if len(curve) else np.nan
        span = curve.response_mean.max() - curve.response_mean.min() if len(curve) >= 2 else np.nan
        # Adjacent mean differences need one eligible well at each standard. A
        # missing mean leaves both neighboring differences unavailable.
        changes = group.sort_values("predictor").response_mean.diff()
        reversals = changes.gt(0) if first.assay == "pH" else changes.lt(0)
        rank = (curve.predictor.rank().corr(curve.response_mean.rank())
                if len(curve) >= 3 and curve.response_mean.nunique() > 1 else np.nan)
        row = {key: first[key] for key in META_KEYS}
        row.update(eligible_level_n=len(curve), scheduled_well_n=int(group.scheduled_well_n.sum()),
                   eligible_well_n=int(group.technical_n.sum()), observed_response_range=span,
                   pooled_technical_sd=sd, range_to_technical_sd=span / sd if sd > 0 else np.nan,
                   spearman_rho=rank, adjacent_direction_reversals=int(reversals.sum()))
        rows.append(row)
    return pd.DataFrame(rows)


def _calibration_tables(root, wells, spectral_levels):
    levels = _level_means(wells)
    chemical = _chemical_references(root)
    iron = levels[levels.series_ID.ne("") & levels.signal.isin(IRON_SIGNALS)].merge(
        chemical, on=["series_ID", "plate_column"], validate="many_to_one")
    for analyte in ("total_Fe", "FeII", "FeIII"):
        iron[analyte + "_mM"] = iron[analyte + "_M"] * 1000 * 200 / iron.final_well_volume_uL
    tables = {"calibration_spectral_summary": spectral_levels, "calibration_level_means": levels,
              "calibration_correction_effects": _precision_pairs(levels),
              "iron_chemical_reference_summary": chemical, "iron_optical_reference_levels": iron}
    tables.update(spectra(tables))
    tables.update(calibration(tables))
    tables.update(reference_regression(tables))
    # A zero residual variance leaves regression uncertainty and tests undefined.
    regressions = tables["feii_regression"]
    inferential = [c for c in regressions if "_ci_" in c or c.endswith("_p")]
    regressions.loc[regressions.sy_x_M.eq(0), inferential] = np.nan
    tables["iron_spectral_replicates"] = iron_spectra(tables)
    return tables


def _handling_tables(root):
    registry = load_optical_registry(root)
    handling = build_robustness_observations(root, registry)
    # Figures 2–3 use the two ferric channels at 380 nm and the two pH ratios.
    handling = handling[handling.signal.isin(["ratio_raw", "ratio_corrected", "A380", "A380_corrected"])]
    media = build_media_observations(root, registry)
    media = media[media.signal.isin(["ratio_raw", "ratio_corrected", "A380", "A380_corrected"])]
    levels = summarize_levels(media)
    pairs = paired_observations(handling, build_matches(handling))
    tables = {"pub_robustness_observations": handling, "pub_robustness_levels": summarize_levels(handling),
              "pub_robustness_paired_wells": pairs, "pub_robustness_paired_summary": summarize_pairs(pairs),
              "pub_media_observations": media, "pub_media_levels": levels, "pub_media_performance": _media_precision(levels),
              "pub_media_readout_comparison": signal_comparisons(media), "pub_media_background": media_background(media),
              "pub_media_interval_diagnostics": media_interval_diagnostics(levels)}
    tables.update(storage(tables))
    tables.update(optical_folds(tables))
    return tables


def _handling_offsets(tables):
    """Separate Holm families for pH and ferric contrasts, equal weight per plate."""
    pairs = tables["pub_robustness_paired_summary"]
    pairs = pairs[~pairs.comparison.isin(["lid_present_minus_absent", "storage_minus_initial"])]
    pairs = pairs[~pairs.comparison.eq("nonUV_minus_UV") | pairs.step_id.eq("initial_23_C")]
    families = []
    for assay, key, panels in (("pH", "ph_handling_effect", "2e–h"), ("iron", "iron_handling_effect", "3g–l")):
        fields = ["comparison", "plate_type", "step_id", "signal"]
        if assay == "iron":
            fields += ["acid_M", "condition"]
        rows = []
        for identity, group in pairs[pairs.assay.eq(assay)].groupby(fields, dropna=False):
            if group.candidate_physical_plate_id.duplicated().any():
                raise ValueError("Duplicate physical plate in a handling contrast")
            row = dict(zip(fields, identity))
            row.update(_summary(100 * group.mean_difference_fraction_of_reference_range),
                       paired_well_n=int(group.paired_well_n.sum()),
                       plate_ids=";".join(group.candidate_physical_plate_id.astype(str)))
            rows.append(row)
        families.append(dict(key=key, title=f"{assay} indicator offsets during handling",
            caption="Matched-standard differences divided by the comparator span, averaged within each physical plate "
                    "and expressed as percentages. The initial acquisition is the temperature comparator, matched "
                    "preparations compare plate materials, and the preceding overnight acquisition is the membrane "
                    "comparator. Two-sided t tests against zero use Holm adjustment across estimable contrasts for this indicator.",
            frame=_holm(pd.DataFrame(rows)), figure_panels=[panels], source_tables=["pub_robustness_paired_summary"]))
    return families


def _dye_background(tables):
    rows = []
    for (medium, wavelength), group in tables["pub_media_background"].groupby(["medium_id", "wavelength_nm"]):
        plate_means = group.groupby("physical_plate_id").with_dye_minus_no_dye_AU.mean()
        rows.append(dict(medium=MEDIUM.get(medium, medium), wavelength_nm=wavelength, **_summary(plate_means)))
    return dict(key="media_background_contrast", title="Effect of dye addition on medium background",
                caption="Dyed minus undyed absorbance in AU at 485, 545 and 750 nm. Matched standard differences "
                        "are averaged per plate. Holm adjustment spans all estimable wavelength-by-medium tests. "
                        "Dye addition also changes the final volume from 200 to 210 µL.",
                frame=_holm(pd.DataFrame(rows)), figure_panels=["4a–b"], source_tables=["pub_media_background"])


def _availability(tables):
    summaries, plates = [], []
    for source in ("pub_robustness_observations", "pub_media_observations"):
        paired = pair_canonical_ratios(tables[source])
        if source == "pub_media_observations":
            paired = paired[paired.condition.eq("with_dye")]
            fields = ["experiment_id", "medium_id"]
        else:
            fields = ["experiment_id", "plate_type"]
        for identity, group in paired.groupby(fields, sort=True, dropna=False):
            stats, units = summarize_ratio_yield(group)
            cohort = dict(zip(fields, identity))
            summaries.append(cohort | stats)
            plates.append(units.assign(**cohort))
    return dict(key="optical_ratio_availability", title="Availability of paired optical pH ratios",
                caption="Measured and corrected ratios are paired by acquired well. Percentages are averaged "
                        "equally across physical plates within each experiment and material or medium.",
                frame=pd.DataFrame(summaries), support_tables={"physical_plates": pd.concat(plates, ignore_index=True)},
                source_tables=["pub_robustness_observations", "pub_media_observations"], figure_panels=["2", "4"])


def run(root: Path) -> dict[str, pd.DataFrame]:
    """Read the repository data/ tree and return full-precision optical analysis results."""
    wells, spectral_levels, _ = _read_sources(root)
    tables = _calibration_tables(root, wells, spectral_levels)
    ph = canonical_ph_adapter(wells, signal_value)
    tables.update(analyze_ph_matrix(ph, bootstrap_n=199, seed=1731))
    tables.update(_handling_tables(root))
    return _results(tables)


def _results(tables):
    specifications = [*optical_results(tables), *validate_calibrations(tables), handling_discrimination(tables),
                      storage_sensitivity(tables), *_handling_offsets(tables), _dye_background(tables), _availability(tables)]
    # These numerical summaries supply the individual figure panels, including fitted lines and uncertainty.
    panel_tables = {
        "spectral_replicates": ("1a", "Native wavelength absorbance by physical plate and pH standard"),
        "spectral_means": ("1a", "Native wavelength spectra by measured pH"),
        "calibration_level_means": ("1b–d, 1f", "Initial standard means and technical-well SDs"),
        "response_curves": ("1b–c, 1f", "Fitted responses and pointwise plate-clustered intervals"),
        "iron_spectral_replicates": ("1e", "Ferric minus ferrous spectra at matched concentrations"),
        "iron_chemical_reference_summary": ("1g", "Total and ferrous iron reference concentrations"),
        "feii_regression_curves": ("1g", "Cerimetry regression predictions and confidence intervals"),
        "iron_optical_reference_levels": ("1e, 1h", "Absorbance by wavelength and measured chemical reference"),
        "pub_robustness_levels": ("2a, 2e–f, 3a, 3g–h", "Lid and handling series by standard and plate"),
        "storage_folds": ("2b–c, 3b–e", "Storage changes at each acquisition relative to baseline"),
        "edge_replicates": ("2d, 3f", "Concentration-matched changes in edge versus interior wells"),
        "ph_handling_folds": ("2g–h", "Ratio changes against the common pH indicator reference"),
        "iron_handling_folds": ("3i–l", "Ferric-indicator changes against the interpolated common reference"),
        "pub_media_levels": ("4a–c", "Medium-specific standard means and technical SDs"),
    }
    exports = {key: tables[key] for key in panel_tables}
    guide = [dict(table=key, title=title, caption=title + ".", panels=panels, source_tables="")
             for key, (panels, title) in panel_tables.items()]
    for spec in specifications:
        key = spec["key"]
        exports[key] = spec["frame"]
        guide.append(dict(table=key, title=spec["title"], caption=spec["caption"],
                          panels=", ".join(spec["figure_panels"]), source_tables=", ".join(spec["source_tables"])))
        for name, frame in spec.get("support_tables", {}).items():
            filename = key + "__" + name
            exports[filename] = frame
            guide.append(dict(table=filename, title=name.replace("_", " "), caption=f"Supporting values for {spec['title']}.",
                              panels=", ".join(spec["figure_panels"]), source_tables=key))
    exports["result_guide"] = pd.DataFrame(guide)
    return exports
