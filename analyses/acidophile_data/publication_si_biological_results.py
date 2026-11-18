"""Si biological results."""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy import stats
from .optical_and_treatment_contrasts import wavelength_endpoint_mask


def _finite(values):
    a = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(float)
    return a[np.isfinite(a)]


def _mean(values):
    a = _finite(values)
    return float(np.mean(a)) if len(a) else np.nan


def _sd(values):
    a = _finite(values)
    return float(np.std(a, ddof=1)) if len(a) > 1 else np.nan


def _quantile(values, p):
    a = _finite(values)
    return float(np.quantile(a, p)) if len(a) else np.nan


def _mean_interval(values, minimum_n=3):
    a = _finite(values)
    if len(a) < minimum_n:
        return np.nan, np.nan
    half = stats.t.ppf(.975, len(a)-1) * np.std(a, ddof=1) / np.sqrt(len(a))
    return float(np.mean(a)-half), float(np.mean(a)+half)


def _holms(values):
    """Retain unavailable planned slots in the multiplicity family."""
    a = np.asarray(values, float)
    output = np.full(len(a), np.nan)
    order = np.flatnonzero(np.isfinite(a))
    order = order[np.argsort(a[order], kind="stable")]
    ceiling = 0.
    for rank, index in enumerate(order):
        ceiling = max(ceiling, min(1., (len(a)-rank)*a[index]))
        output[index] = ceiling
    return output


def _ids(frame, field):
    if field not in frame:
        return ""
    return ";".join(sorted({part for value in frame[field].dropna()
                            for part in str(value).split(";") if part}))


def _count_ids(frame, field):
    value = _ids(frame, field)
    return len(value.split(";")) if value else 0


def _result(key, title, caption, frame, display_columns, sources, panels, claims):
    if not len(frame):
        raise ValueError(f"Empty analysis table: {key}")
    if any(column not in frame for column in display_columns):
        raise ValueError(f"Missing display column in {key}")
    results = [column for column in frame if any(token in column for token in
        ("effect", "difference", "minus", "variance_ratio", "cv_percent", "outcome_sd", "deviation_iqr",
         "log2_fold", "net_change", "net_endpoint", "snr_change", "magnitude_change", "agreement_percent",
         "increased_percent", "p_holm", "ci_low", "ci_high"))]
    return dict(key=key, title=title, caption=caption, frame=frame, results_columns=results,
                display_columns=display_columns, source_tables=sources,
                figure_panels=panels, claims=claims)


def _repeatability(tables):
    from .publication_baseline_repeatability import (
        baseline_batches, flask_ga_batches, METRICS, baseline_repeatability_tests,
        flask_pairwise_repeatability_tests,
    )
    source = tables["condition_repeatability_batches"]
    baseline = baseline_repeatability_tests(source)
    paired = flask_pairwise_repeatability_tests(source)
    baseline = baseline.rename(columns={"log_variance_ratio": "log_outcome_variance_ratio"})
    paired = paired.rename(columns={"log_variance_ratio": "log_outcome_variance_ratio"})
    output = [_result("bio_baseline_repeatability", "Plate–flask baseline dispersion tests",
        "Two planned, two-sided F tests compare between-starting-batch variance of natural-log native "
        "AUC and time to observed maximum for standard no-additive AFM1 plates versus flasks. "
        "Each physical vessel is reduced to its batch mean before testing; the 345 plates represent "
        "eight batches and the 15 flasks represent five. AUC is the signed trapezoidal integral to "
        "each trajectory's own first observed maximum, without later readings. Unit conversion does "
        "not affect log variance. The effect is Var[ln(Yplate)]/Var[ln(Yflask)], not its logarithm. "
        "The F distribution assumes independent lognormal batch outcomes; intervals are marginal 95% "
        "intervals and Holm P covers both planned tests. Format, readout, cohort and schedule are "
        "confounded. Nonsignificance does not establish equal repeatability or equivalence.",
        baseline, ["metric", "plate_batch_n", "flask_batch_n", "plate_cv_percent", "flask_cv_percent",
                   "log_outcome_variance_ratio", "ci_low", "ci_high", "numerator_df", "denominator_df",
                   "p_value", "p_holm"],
        ["condition_repeatability_batches"], ["5c", "5d"],
        ["Baseline multiplicative dispersion is not detectably different between these observational plate and flask cohorts."]),
        _result("bio_flask_repeatability", "Paired flask dispersion tests across Ga conditions",
        "Exploratory two-sided Pitman–Morgan tests compare the variances of natural-log native batch "
        "AUC or time to observed maximum at 0, 1 and 100 mM Ga. Shared starting batches are paired. "
        "The 0/1 and 0/100 mM comparisons each use three paired batches; 1/100 mM has one shared "
        "batch and its variance test is unavailable. The effect is the numerator/denominator log-outcome "
        "variance ratio. Marginal 95% intervals invert the paired test. Holm correction includes all "
        "six planned endpoint-by-condition contrasts, including the two unavailable slots. The "
        "bivariate-normal log-outcome assumption is not assessable at three pairs. Results are "
        "conditional, exploratory and do not establish equivalence; unmatched controls remain in Figure 5.",
        paired, ["metric", "numerator_Ga_mM", "denominator_Ga_mM", "matched_batch_n",
                 "log_outcome_variance_ratio", "ci_low", "ci_high", "statistic", "df", "p_value", "p_holm", "status"],
        ["condition_repeatability_batches"], ["5c", "5d"],
        ["The paired flask comparisons do not establish a Ga-associated change in repeatability; one comparison lacks sufficient paired batches."])]
    rows = []
    for metric in METRICS:
        plate, flask = baseline_batches(tables["condition_repeatability_batches"], metric)
        for label, group in [("plate 0", plate), ("flask 0", flask),
                             ("flask 1", flask_ga_batches(tables["condition_repeatability_batches"], metric, 1.)),
                             ("flask 100", flask_ga_batches(tables["condition_repeatability_batches"], metric, 100.))]:
            x = _finite(group.raw_value)
            mean, sd = _mean(x), _sd(x)
            rows.append(dict(metric=metric, cohort=label, batch_n=len(x),
                vessel_n=_count_ids(group, "unit_ids"), native_mean=mean, native_sd=sd,
                native_unit=group.native_unit.iloc[0], cv_percent=100*sd/abs(mean) if mean else np.nan,
                log_outcome_sd=_sd(np.log(x[x > 0])),
                centered_deviation_iqr_pp=100*(_quantile(x,.75)-_quantile(x,.25))/abs(mean) if mean else np.nan,
                final_read_maximum_vessels=int(group.final_read_maximum_n.sum()),
                batch_ids=_ids(group, "batch_id"), status="descriptive_between_starting_batch_dispersion"))
    output.append(_result("bio_repeatability_distribution", "Native batch repeatability of the displayed cohorts",
        "Precision estimates use the exact native batch AUC and timing values underlying Figure 5c–d. "
        "Physical vessels are equally averaged within starting batch, then batches are equally weighted. "
        "CV is 100×sample SD/absolute mean; log-outcome SD measures multiplicative dispersion. "
        "The interquartile range of 100×(Ybatch/mean(Y)−1) describes the plotted centered distribution. "
        "Its mean is zero by construction and is not a testable between-group effect. AUC units differ "
        "between formats, so native means are not compared. Final-read maxima do not establish a plateau. "
        "Only standard no-additive AFM1 plate controls and Ga-only AFM1 flask designs are included; "
        "deliberately adjusted-pH and additive-treated plate conditions are excluded from these plotted cohorts.",
        pd.DataFrame(rows), ["metric", "cohort", "batch_n", "vessel_n", "cv_percent", "log_outcome_sd",
                            "centered_deviation_iqr_pp", "final_read_maximum_vessels"],
        ["condition_repeatability_batches"], ["5a", "5b", "5c", "5d"],
        ["The displayed repeatability distributions quantify between-batch dispersion rather than differences in their constructed zero means."]))
    return output


def _od(tables):
    output = []
    source = tables["od_dose_batch_means"]
    summaries = tables["od_dose_summary"]
    for key, metric, title, panels, offset, unit in [
        ("od_ga_response", "response_percent", "Ga effects on control-relative optical integrated response", ["6c"], 100., "percentage points"),
        ("od_ga_timing", "delay_d", "Ga-associated optical timing delays", ["6d"], 0., "d")]:
        rows = []
        for (form, signal, ga), group in source[~source.is_reference].groupby(["format", "signal", "Ga_mM"]):
            values = _finite(group[metric])-offset
            if not len(values):
                continue
            low, high = _mean_interval(values)
            summary = summaries[summaries["format"].eq(form)&summaries.signal.eq(signal)&summaries.Ga_mM.eq(ga)].iloc[0]
            prefix = "response" if metric=="response_percent" else "delay"
            rows.append(dict(format=form, readout=signal, Ga_mM=ga, effect=_mean(values), unit=unit,
                batch_n=len(values), candidate_batch_n=len(group), physical_unit_n=int(summary[prefix+"_unit_n"]),
                between_batch_sd=_sd(values), within_batch_unit_sd=float(summary[prefix+"_unit_sd"]) if len(values)==1 else np.nan,
                ci_low=low, ci_high=high,
                final_read_maximum_n=int(group.final_read_maximum_n.sum()),
                status="descriptive_fewer_than_three_batches" if len(values)<3 else "conditional_batch_mean_t_interval"))
        caption = ("Effects use physical-vessel trajectories reduced to equal-weight starting-batch means, each "
            "with its own same-batch no-Ga reference. "+
            ("The response effect is 100×AUC_Ga/mean(AUC_0)−100 percentage points, where each signed trapezoidal "
             "area stops at that trajectory's own first observed maximum. Negative effects indicate a smaller "
             "integrated response, not a calibrated inhibition or growth-rate estimate. " if metric=="response_percent" else
             "The timing effect is tmax,Ga−mean(tmax,0), in days. It is not a fold of a delay that can be zero or negative. "
             "Non-rising trajectories have unavailable timing; final-read maxima are observed endpoints, not established peaks. ")+
            "All Ga doses with an available estimate are retained, including doses outside the compact plotted subset. "
            "Raw and A750-corrected plate readouts use paired wells. Flask apparent OD is calculated per eligible "
            "TE-buffer-blank-subtracted RFU observation as (RFU/74930)^(1/1.386)−0.001 before aggregation; "
            "missing initial RFU is not imputed. This AF-derived relation is not a universal biomass calibration. "
            "Marginal Student-t 95% intervals are provided only for at least three independent starting batches; "
            "they assume approximately normal batch effects and do not remove shared-control covariance across doses. "
            "Within-batch physical-unit SD is shown only when a single batch contributes; it is not biological uncertainty. "
            "No cross-format neighboring-dose quotient or significance test is calculated.")
        frame = pd.DataFrame(rows)
        display = ["format", "readout", "Ga_mM", "batch_n", "physical_unit_n", "effect", "within_batch_unit_sd",
                   "final_read_maximum_n"]
        if frame.ci_low.notna().any():
            display += ["between_batch_sd", "ci_low", "ci_high"]
        output.append(_result(key, title, caption, frame, display,
            ["od_dose_batch_means", "od_dose_summary"], panels,
            ["Optical response and timing estimates remain readout-specific and incompletely supported across Ga doses."]))
    return output


def _wavelength(tables):
    source = tables["citrate_wavelength_pairs"]
    q = source[wavelength_endpoint_mask(source)].copy()
    if q.duplicated(["experiment_id", "group_id", "column", "batch_id", "plate_id", "ordinal"]).any():
        raise ValueError("Duplicated physical-plate wavelength endpoint")
    q["paired_difference"] = q.net_change_380-q.net_change_428
    rows = []
    for (citric, ga), group in q.groupby(["citrate_mM", "Ga_mM"]):
        finite = group[np.isfinite(group.net_change_380)&np.isfinite(group.net_change_428)]
        batch = finite.groupby("batch_id")[["net_change_380", "net_change_428", "paired_difference"]].mean()
        delta = batch.paired_difference
        low, high = _mean_interval(delta)
        test = stats.ttest_1samp(delta, 0.) if len(delta)>=3 else None
        rows.append(dict(citric_acid_mM=citric, Ga_mM=ga, paired_plate_n=len(finite), batch_n=len(batch),
            net_380_mean_AU=_mean(batch.net_change_380), net_428_mean_AU=_mean(batch.net_change_428),
            paired_380_minus_428_AU=_mean(delta), between_batch_difference_sd_AU=_sd(delta),
            ci_low=low, ci_high=high, statistic=float(test.statistic) if test else np.nan,
            df=len(delta)-1 if test else np.nan, p_value=float(test.pvalue) if test else np.nan,
            days_min=float(finite.days_mean_380.min()), days_max=float(finite.days_mean_380.max()),
            batch_ids=_ids(finite, "batch_id"), plate_ids=_ids(finite, "plate_id"),
            status="exploratory_paired_starting_batch_t_test" if test else "insufficient_starting_batches"))
    frame = pd.DataFrame(rows)
    frame["p_holm"] = _holms(frame.p_value)
    return _result("citrate_wavelength_contrast", "Paired wavelength contrasts at the final acquisition",
        "For every citric-acid/Ga formulation in Figure 7d, ΔAλ,net=[Aλ,AF(t)−Aλ,AF(0)]−"
        "[Aλ,Ab(t)−Aλ,Ab(0)] is evaluated at each physical plate's actual final acquisition using only plates with both "
        "380 and 428 nm values. The contrast is ΔA380,net−ΔA428,net within a physical plate. Paired differences "
        "are averaged within starting batch, then batches receive equal weight. Zero citric acid combines "
        "the compatible preparations before batch aggregation. The 60 paired plate endpoints are not "
        "60 independent batches; each formulation has three or four batch results. Exploratory two-sided "
        "one-sample t tests of batch differences assume independent normally distributed batch contrasts. "
        "Intervals are marginal 95%; Holm adjustment spans every formulation in this table. Different "
        "formulations share batches. A larger optical response at 380 nm does not establish better "
        "concentration resolution, biomass accuracy or universal wavelength superiority.", frame,
        ["citric_acid_mM", "Ga_mM", "paired_plate_n", "batch_n", "paired_380_minus_428_AU",
         "ci_low", "ci_high", "df", "p_value", "p_holm"],
        ["citrate_wavelength_pairs"], ["7d"],
        ["The paired wavelength difference is quantified on common scheduled plate endpoints without treating plates or doses as independent batches."])


def _format_effects(tables):
    q = tables["format_folds"]
    q = q[q.readout.eq("raw_A380")].copy()
    groups = ["additive", "plate_dose_mM", "flask_dose_mM", "Ga_mM", "match_status"]
    output = []
    for key, field, title, panel in [
        ("product_format_response", "response_log2_fold", "Plate–flask integrated-response correspondence", "7e"),
        ("product_format_timing", "time_log2_fold", "Plate–flask time-to-maximum correspondence", "7f")]:
        rows = []
        for labels, group in q.groupby(groups, dropna=False):
            valid = group[np.isfinite(group[field])]
            values = _finite(valid[field])
            row = dict(zip(groups, labels))
            row.update(comparison_context_n=len(group), defined_context_n=len(values),
                plate_batch_n=valid.plate_batch.nunique(), flask_batch_n=valid.flask_batch.nunique(),
                plate_vessel_n=_count_ids(valid, "plate_replicate_ids"), flask_vessel_n=_count_ids(valid, "flask_replicate_ids"),
                median_log2_fold=_quantile(values,.5), q25_log2_fold=_quantile(values,.25),
                q75_log2_fold=_quantile(values,.75), minimum_log2_fold=_quantile(values,0), maximum_log2_fold=_quantile(values,1),
                median_fold=2**_quantile(values,.5),
                final_read_maximum_context_n=int(valid.final_read_maximum_present.eq(True).sum()),
                status="descriptive_shared_cohorts_and_references" if len(values) else "unavailable_nonpositive_or_missing_metric")
            rows.append(row)
        caption = ("Comparison-level effects from the metabolic byproduct indicators are summarized within each actual "
            "plate/flask dose pairing; no new regression, independence assumption or neighboring-dose matching is imposed. "+
            ("The effect is log2[(AUC_treated/AUC_control)_plate/(AUC_treated/AUC_control)_flask]. Each signed trapezoidal "
             "AUC ends at its trajectory's own first observed maximum. " if field=="response_log2_fold" else
             "The effect is log2(tmax,plate/tmax,flask), where times are measured from the recorded initial reading "
             "to each trajectory's first observed maximum; final-read maxima do not establish an eventual plateau. ")+
            "Zero means equal, +1 twofold and −1 half. Nonpositive or missing quantities have no finite log2 fold. "
            "The median, quartiles and range are descriptive distributions across attempted cohort pairings; they are "
            "not confidence intervals. A given flask cohort/reference is reused for multiple plate batches and doses; "
            "the context count is not an independent sample size. Exact and neighboring reported doses remain explicit, "
            "and medium, inoculum, readout, vessel and schedule differences remain confounded. No population CI or P "
            "value is assigned to these cross-format effects. Ga_mM is the additional Ga factor for citric acid; "
            "for additive=Ga, the two dose columns contain the Ga exposure.")
        output.append(_result(key, title, caption, pd.DataFrame(rows),
            ["additive", "plate_dose_mM", "flask_dose_mM", "Ga_mM", "match_status", "defined_context_n",
             "median_log2_fold", "q25_log2_fold", "q75_log2_fold", "final_read_maximum_context_n"],
            ["format_folds"], [panel],
            ["Cross-format effect estimates preserve exact versus neighboring dose contexts and do not establish assay interchangeability."]))
    return output


def _lf_benchmark(tables):
    from .measurement_summaries import lf_baseline_and_nearby_treatments
    from .publication_lf_benchmarks import at_exact_endpoint
    q = lf_baseline_and_nearby_treatments(tables["validation_lf_plate_trajectories"])
    target = float(q.days_mean.max())
    _, batches, af = at_exact_endpoint(tables, target)
    af = af.iloc[0]
    rows = []
    for label, group in q.groupby("display_condition", sort=False):
        endpoint = group[group.ordinal.eq(group.groupby(["batch_id", "plate_id", "column"]).ordinal.transform("max"))]
        endpoint = endpoint.drop_duplicates(["batch_id", "plate_id", "column", "ordinal"])
        finite = endpoint[np.isfinite(endpoint.biotic_change)&np.isfinite(endpoint.abiotic_change)]
        net = finite.biotic_change-finite.abiotic_change
        rows.append(dict(LF_plate_context=label, LF_plate_n=len(finite), LF_batch_n=finite.batch_id.nunique(),
            LF_inoculum_culture_n=_count_ids(finite,"culture_ids"), LF_biotic_mean_AU=_mean(finite.biotic_change),
            LF_abiotic_mean_AU=_mean(finite.abiotic_change), LF_net_change_AU=_mean(net),
            LF_net_between_microplate_sd_AU=_sd(net), AF_reference_mean_AU=float(af["mean"]),
            AF_reference_batch_sd_AU=float(af.sd), AF_batch_n=int(af.batch_n), AF_plate_n=int(af.physical_plate_n),
            LF_biotic_minus_AF_reference_AU=_mean(finite.biotic_change)-float(af["mean"]),
            AF_reference_day_d=target, LF_endpoint_min_d=float(finite.days_mean.min()), LF_endpoint_max_d=float(finite.days_mean.max()),
            LF_Ga_mM=142.8 if label=="Ga" else 0., LF_In_mM=6.22 if label=="In" else 0.,
            status="descriptive_species_cohort_and_medium_confounding"))
    result = _result("lf_plate_benchmark", "LF plate separation relative to matched abiotic and AF controls",
        "Each LF physical plate is evaluated at its final scheduled acquisition. A missing endpoint is unevaluable. The LF net effect is its biotic initial-subtracted A380 change minus "
        "the corresponding abiotic change. Technical wells are averaged within each microplate, and SD "
        "describes variation among independent microplates in the same acquisition-date group. The AF reference "
        "is interpolated to the exact last LF acquisition, using only directly bracketing finite AF readings "
        "without extrapolation or missing-value bridges. Eligible AF plates are averaged within starting "
        "batch, then eight batches are equally weighted. The contextual LF-biotic-minus-AF gap compares "
        "different organisms, cohorts and media; no species-effect CI or test is assigned. LF Ga and In "
        "plate concentrations are 142.8 and 6.22 mM, respectively, distinct from the 100 and 10 mM tube "
        "treatments. Small biotic–abiotic separation does not demonstrate starter nonviability.",
        pd.DataFrame(rows), ["LF_plate_context", "LF_plate_n", "LF_batch_n", "LF_net_change_AU",
                            "LF_net_between_microplate_sd_AU", "AF_reference_mean_AU", "AF_reference_batch_sd_AU",
                            "LF_biotic_minus_AF_reference_AU"],
        ["validation_lf_plate_trajectories", "af_iron_control_plate_trajectories"], ["7g"],
        ["LF baseline and selected additive plate responses have a small matched abiotic separation relative to the AF contextual benchmark."])
    result["support_tables"] = {"batch_estimates": batches}
    return result


def _tube_effects(tables, sulfur=False):
    q = tables["tube_od600"]
    temperature = 30 if sulfur else 37
    q = q[q.record_role.eq("endpoint") & q.incubation_temperature_C.eq(temperature)].copy()
    if sulfur:
        q = q[q.medium.eq("SM4_tetrathionate") & q.nominal_Ga_mM.eq(0)]
        organisms = ["AF", "AT"]
        fields = [("pH", "pH units"), ("apparent_OD600", "AU")]
    else:
        q = q[q.medium.eq("HH") & q.organism.isin(["LF", "Ab"])]
        organisms = ["LF"]
        fields = [("orp_mV", "mV"), ("apparent_OD600", "AU")]
    rows = []
    for (ga, indium), group in q.groupby(["nominal_Ga_mM", "nominal_In_mM"]):
        ab = group[group.organism.eq("Ab")]
        for organism in organisms:
            biotic = group[group.organism.eq(organism)]
            matching = ["batch", "medium_ID", "treatment_ID", "incubation_temperature_C", "sampling_event"]
            if len(group[matching].drop_duplicates()) != 1:
                raise ValueError("Tube endpoint contrasts require one explicitly matched cohort and acquisition")
            if ab.tube_ID.duplicated().any() or biotic.tube_ID.duplicated().any():
                raise ValueError("Repeated tube endpoints cannot be counted as independent observations")
            for metric, unit in fields:
                a, b = _finite(biotic[metric]), _finite(ab[metric])
                rows.append(dict(organism=organism, Ga_mM=ga, In_mM=indium, metric=metric, unit=unit,
                    biotic_n=len(a), abiotic_n=len(b), biotic_mean=_mean(a), abiotic_mean=_mean(b),
                    biotic_minus_abiotic=_mean(a)-_mean(b), biotic_sd=_sd(a), abiotic_sd=_sd(b),
                    biotic_starter_n=biotic.inoculum_culture_ID.nunique(),
                    endpoint_d=_mean(biotic.time_average_d), independent_batch_n=biotic.batch.nunique(),
                    status="descriptive_one_starting_cohort_shared_initial_reference"))
    return _result("sulfur_tube_contrasts" if sulfur else "tube_context_contrasts",
        "Tetrathionate tube endpoint contrasts" if sulfur else "Fresh LF tube endpoint contrasts",
        f"Endpoint biotic-minus-abiotic contrasts use the recorded final sampling event at {temperature} °C, with three independent "
        "physical tubes per available condition; a tube is never multiplied by its repeated observations. "
        "Separate group SDs describe observed dispersion, not uncertainty in a species or treatment effect. "
        "Initial medium references are shared and are not additional tube measurements. No population test "
        "or CI is assigned to the single starting cohort. Apparent OD600 is calculated per eligible RFU "
        "observation as (RFU/74930)^(1/1.386), after the deposited TE-buffer subtraction, without the "
        "flask-specific −0.001 offset; negative/flagged RFU remain unavailable. This AF-derived relation "
        "is not an organism-specific biomass calibration. "+
        ("Only baseline tetrathionate tubes are analyzed here, matching Figure 8g–h. AF and AT contrasts "
         "share an abiotic group in medium containing 5 mM tetrathionate."
         if sulfur else
         "ORP is relative to Ag/AgCl, 3 M KCl. Baseline, 100 mM Ga and 10 mM In are matched separately to "
         "their abiotic tubes. The tube cultures use fresh LF starters and describe a separate cohort from the microplate cultures."),
        pd.DataFrame(rows), ["organism", "Ga_mM", "In_mM", "metric", "biotic_n", "abiotic_n",
                            "biotic_minus_abiotic", "biotic_sd", "abiotic_sd", "endpoint_d"],
        ["tube_od600"], ["8g", "8h"] if sulfur else ["7h", "7i"],
        ["Tube biotic–abiotic endpoint contrasts provide independent experimental context without establishing plate–tube equivalence or viability."])


def _sulfur(tables):
    from .publication_sulfur_endpoint_support import build_tables as build_endpoint_support
    q = tables["sulfur_baseline_trajectories"].copy()
    registry = tables["sulfur_baseline_registry"][["reference_id", "initial_OD600_200_uL"]]
    q = q.merge(registry, on="reference_id", validate="many_to_one")
    q["cohort"] = q.apply(lambda r: (f"AF thiosulfate OD {r.initial_OD600_200_uL:g}" if str(r.reference_id).startswith("af_ts_")
        else "AF tetrathionate 5 mM" if r.reference_id=="af_tt_5"
        else "AF tetrathionate 15 mM" if r.reference_id=="af_tt_15"
        else "AT tetrathionate 20 mM" if r.reference_id=="at_tt" else "AT thiosulfate 20 mM"), axis=1)
    identifiers = ["reference_id", "signal", "batch_id", "plate_id"]
    q = q[q.ordinal.eq(q.groupby(identifiers).ordinal.transform("max"))]
    if q.duplicated(identifiers).any():
        raise ValueError("Duplicate baseline plate endpoint")
    rows = []
    for (cohort, signal), group in q.groupby(["cohort", "signal"]):
        g = group[np.isfinite(group.biotic_change)&np.isfinite(group.abiotic_change)].copy()
        difference = g.biotic_change-g.abiotic_change
        rows.append(dict(cohort=cohort, readout="R*" if signal.startswith("corrected") else "R",
            candidate_plate_n=len(group), paired_endpoint_plate_n=len(g), batch_n=g.batch_id.nunique(),
            inoculum_culture_n=_count_ids(g,"culture_ids"),
            net_endpoint_change=_mean(difference), net_between_microplate_sd=_sd(difference),
            minimum_net_change=_quantile(difference,0), maximum_net_change=_quantile(difference,1),
            endpoint_min_d=float(g.days_mean.min()), endpoint_max_d=float(g.days_mean.max()),
            status="descriptive_independent_microplates_single_acquisition_date"))
    output = [_result("sulfur_baseline_contrast", "Soluble-sulfur baseline biotic–abiotic separation",
        "At each baseline physical plate's final scheduled acquisition, net ratio change is "
        "(Rbiotic,t−Rbiotic,0)−(Rabiotic,t−Rabiotic,0), and analogously for R*. R=A545/A485 and "
        "R*=(A545−A750)/(A485−A750). Both members must be finite at the scheduled endpoint. Compatible AF thiosulfate no-addition series are combined "
        "within nominal initial OD. Each microplate used its own set of three inoculum cultures; their "
        "well responses are averaged within microplate. SD and range describe variation among independent "
        "microplates from one acquisition-date group. No population CI or test is assigned in this descriptive "
        "comparison. Concentrations retain the deposited nominal basis; dye dilution "
        "is retained in the source registry. The 15 mM AF tetrathionate comparator is retained as an "
        "analysis despite its omission from the compact figure. A net ratio change is an optical "
        "contrast, not calibrated pH or direct evidence of growth or viability.",
        pd.DataFrame(rows), ["cohort", "readout", "paired_endpoint_plate_n", "batch_n", "net_endpoint_change",
                            "net_between_microplate_sd", "minimum_net_change", "maximum_net_change"],
        ["sulfur_baseline_trajectories", "sulfur_baseline_registry"], ["8a", "8b", "8c", "8d", "8e", "8f"],
        ["Raw and corrected baseline ratio separation differs across the retained soluble-sulfur references without demonstrating universal growth or incompatibility."])]
    source = tables["sulfur_plate_evidence"]
    index = ["species", "experiment_id", "group_id", "column", "batch_id", "plate_id", "medium_id"]
    fields = ["final_snr_signed_common", "final_snr_magnitude_common"]
    raw = source[source.signal.eq("raw_R545_485")][index+fields+['culture_ids']]
    corrected = source[source.signal.eq("corrected_R545_485")][index+fields]
    paired = raw.merge(corrected, on=index, how="outer", suffixes=("_raw", "_corrected"), validate="one_to_one")
    rows = []
    for (species, medium), group in paired.groupby(["species", "medium_id"]):
        finite = group[np.isfinite(group.final_snr_signed_common_raw)&np.isfinite(group.final_snr_signed_common_corrected)&
                       np.isfinite(group.final_snr_magnitude_common_raw)&np.isfinite(group.final_snr_magnitude_common_corrected)]
        signed = finite.final_snr_signed_common_corrected-finite.final_snr_signed_common_raw
        magnitude = finite.final_snr_magnitude_common_corrected-finite.final_snr_magnitude_common_raw
        agrees = np.sign(finite.final_snr_signed_common_corrected).eq(np.sign(finite.final_snr_signed_common_raw))
        rows.append(dict(species=species, medium=medium, candidate_plate_formulation_pairs=len(group), paired_n=len(finite),
            unavailable_pair_n=len(group)-len(finite), batch_n=finite.batch_id.nunique(), inoculum_culture_n=_count_ids(finite,'culture_ids'),
            median_signed_snr_change=_quantile(signed,.5), signed_snr_q25=_quantile(signed,.25), signed_snr_q75=_quantile(signed,.75),
            median_magnitude_change=_quantile(magnitude,.5), magnitude_q25=_quantile(magnitude,.25), magnitude_q75=_quantile(magnitude,.75),
            magnitude_increased_n=int(magnitude.gt(0).sum()), magnitude_increased_percent=100*float(magnitude.gt(0).mean()) if len(finite) else np.nan,
            direction_agreement_percent=100*float(agrees.mean()) if len(finite) else np.nan,
            status="descriptive_paired_indicator_effect_shared_microplates_across_conditions"))
    output.append(_result("sulfur_paired_correction", "Paired correction effects on sulfur technical signal-to-noise",
        "Raw and corrected ratio readouts are paired within the same physical plate/formulation at the final "
        "scheduled acquisition using common finite wells. For each contributing starter/well observation, "
        "the signed technical score is net ratio change divided by sqrt(SD_abiotic,t^2+SD_abiotic,0^2), "
        "an independent-read noise approximation. The signed and absolute scores are averaged separately "
        "within physical plate; the absolute value is taken before averaging, so positive and negative "
        "responses cannot cancel in the magnitude score. Effects are corrected minus raw scores. "
        "Medians and quartiles summarize "
        "paired effects; an increase count uses a strict positive magnitude difference. These are "
        "technical sensitivity summaries across conditions that can share a microplate, not independent "
        "biological tests or a calibrated limit of detection. No binomial CI or P is assigned to pair "
        "fractions. A larger absolute score may accompany a negative response and does not establish "
        "growth or universal correction benefit. Missing common support remains unavailable.",
        pd.DataFrame(rows), ["species", "medium", "paired_n", "unavailable_pair_n", "median_signed_snr_change",
                            "median_magnitude_change", "magnitude_q25", "magnitude_q75", "magnitude_increased_percent",
                            "direction_agreement_percent"],
        ["sulfur_plate_evidence"], ["8a", "8b", "8c", "8d", "8e", "8f"],
        ["A750 correction changes signed and magnitude technical signal-to-noise heterogeneously rather than providing universal improvement."]))
    support = build_endpoint_support(tables)
    baseline = support['sulfur_baseline_endpoint_support']
    correction = support['sulfur_correction_endpoint_support']
    for result, additions, keys in [(output[0], baseline, ['cohort','readout']),
                                    (output[1], correction, ['species','medium'])]:
        overlap = set(result['frame']).intersection(additions) - set(keys)
        if overlap:
            raise ValueError(f'Sulfur endpoint support would overwrite result fields: {overlap}')
        original = result['frame']
        result['frame'] = original.merge(additions, on=keys, how='left', validate='one_to_one', sort=False)
        if len(result['frame']) != len(original):
            raise ValueError('Sulfur endpoint support changes the number of summary rows')
    if not output[0]['frame'].paired_endpoint_plate_n.eq(
            output[0]['frame'].marginal_eligible_endpoint_plate_n).all():
        raise ValueError('Sulfur baseline marginal support differs from the baseline endpoint counts')
    if not output[1]['frame'].paired_n.eq(
            output[1]['frame'].common_eligible_endpoint_plate_formulation_n).all():
        raise ValueError('Sulfur paired endpoint support differs from the paired S/N counts')
    output[0]['support_tables'] = {
        name.removeprefix('sulfur_'): frame for name, frame in support.items()
        if name.startswith('sulfur_baseline_')}
    output[1]['support_tables'] = {
        name.removeprefix('sulfur_'): frame for name, frame in support.items()
        if not name.startswith('sulfur_baseline_')}
    output[1]['source_tables'] += ['context_sulfur_channels','validation_sulfur_trajectory_evidence']
    output[0]['caption'] += (
        ' Additional support identifies union candidate, marginally eligible raw/corrected and '
        'commonly eligible endpoint plates, retaining any final-ordinal mismatch. The '
        'raw and corrected summaries retain their respective marginal populations; their '
        'difference is not a paired correction estimate.')
    output[1]['caption'] += (
        ' Endpoint support traces acquired biotic wells through required channels, the '
        'ratio-denominator guard, initial references, complete matched abiotic controls, finite '
        'net contrast and finite positive technical noise. Joint availability is reported at '
        'well and plate/formulation levels with mutually exclusive primary exclusions per readout. '
        'Marginally available readouts need not share eligible wells. Mean net-response signs are '
        'computed only on the exact common endpoint wells admitted to this paired score analysis; '
        'same, opposite and any-zero outcomes have explicit denominators. The direction '
        'agreement statistic uses the sign of mean signed S/N. These '
        'endpoint comparisons make no temporal-persistence claim.')
    return output


def build(tables):
    """Return the analysis-only results for Figures 5–8."""
    output = _repeatability(tables)
    output.extend(_od(tables))
    output.append(_wavelength(tables))
    output.extend(_format_effects(tables))
    output.append(_lf_benchmark(tables))
    output.append(_tube_effects(tables))
    output.extend(_sulfur(tables))
    output.append(_tube_effects(tables, sulfur=True))
    if len({entry["key"] for entry in output}) != len(output):
        raise ValueError("Duplicate biological analysis key")
    return output
