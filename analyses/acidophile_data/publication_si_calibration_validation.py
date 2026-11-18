"""Si calibration validation."""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, t as student_t
from .publication_calibration_fourpl import fit_ferric_response


SIGNALS = {'ratio_raw': 'R', 'ratio_corrected': 'R*', 'a380': 'A380',
           'A380': 'A380', 'A380_corrected': 'A*380', 'A428': 'A428',
           'A428_corrected': 'A*428'}


MEDIA = {'Starkey_brown': 'St-BS', 'Starkey_colloidal': 'St-YS'}


def _ci(mean, low, high):
    if not np.isfinite(mean):
        return 'NA'
    return (f'{mean:.3g} [{low:.3g}, {high:.3g}]' if np.isfinite(low + high)
            else f'{mean:.3g} (NA)')


def _mean_ci(values, positive=False):
    a = np.asarray(values, float)
    a = a[np.isfinite(a)]
    if positive:
        a = a[a > 0]
        a = np.log(a)
    n = len(a)
    mean = float(a.mean()) if n else np.nan
    lo = hi = np.nan
    if n > 1:
        delta = student_t.ppf(.975, n-1) * np.std(a, ddof=1) / np.sqrt(n)
        lo, hi = mean-delta, mean+delta
    if positive:
        mean, lo, hi = np.exp([mean, lo, hi])
    return n, float(mean), float(lo), float(hi)


def _table(key, title, caption, frame, columns, sources, panels, claims, results):
    labels = {'acid_M': 'H2SO4 (M)', 'validation': 'Withholding',
              'nominal_reference_window_mM': 'Nominal window (mM)',
              'attempted_n': 'Attempted', 'predicted_n': 'Predicted',
              'outside_reference_support_n': 'Reference outside range',
              'response_or_inverse_failure_n': 'Inverse unavailable',
              'bias_mM': 'Bias (mM)', 'rmse_mM': 'RMSE (mM)',
              'median_absolute_relative_error_percent': 'Median absolute error (%)',
              'paired_prediction_n': 'Paired predictions', 'shared_rmse_AU': 'Shared RMSE (AU)',
              'acid_specific_rmse_AU': 'Acid-specific RMSE (AU)',
              'rmse_change_AU': 'RMSE change (AU)',
              'shared_max_absolute_level_bias_AU': 'Shared max. level bias (AU)',
              'acid_specific_max_absolute_level_bias_AU': 'Acid-specific max. level bias (AU)',
              'family': 'Family', 'cohort': 'Cohort', 'readout': 'Readout',
              'reference_unit': 'Reference unit', 'reference_interval': 'Reference interval',
              'response_change_n': 'Plates', 'response_change_ci': 'Change [95% CI]',
              'secant_sensitivity_ci': 'Sensitivity [95% CI]',
              'change_over_noise_ci': 'Change/noise [95% CI]',
              'opposite_direction_plate_n': 'Reversals'}
    return dict(key=key, title=title, caption=caption, frame=frame.reset_index(drop=True),
                display_columns=columns, source_tables=sources, figure_panels=panels,
                claims=claims, result_columns=results,
                column_labels={c: labels[c] for c in columns if c in labels})


def _ferric(tables):
    q = tables['calibration_level_means']
    q = q[q.assay.eq('iron') & q.condition.eq('no_thymol_blue') & q.signal.eq('a380')
          & q.experiment_id.str.startswith('iron_configuration_')
          & q.experiment_id.str.endswith('_23_C_UV')].copy()
    q = q.sort_values(['physical_plate_id', 'axis_value']).reset_index(drop=True)
    if q.duplicated(['physical_plate_id', 'axis_value']).any():
        raise ValueError('Ferric validation requires one mean per physical plate/reference.')
    q['dilution_rank'] = q.groupby('physical_plate_id').cumcount()
    q['observation_id'] = np.arange(len(q))
    if set(q.groupby('physical_plate_id').size()) != {12}:
        raise ValueError('The ferric calibration requires twelve common dilution ranks.')
    # Validate that rank follows the same source well-column map across plates.
    if q.groupby('dilution_rank').plate_column.nunique().max() != 1:
        raise ValueError('Common dilution-rank/well-column correspondence changed.')
    return q


def _valid_fit(fit):
    return bool(fit.converged and np.isfinite(fit.estimates).all())


def _inverse(value, fit, low, high):
    if not _valid_fit(fit):
        return np.nan, 'fit_failed'
    if not np.isfinite(value):
        return np.nan, 'response_unavailable'
    bottom, span, midpoint, shape = fit.estimates
    supported_y = fit.predict([low, high])
    if value < supported_y[0] or value > supported_y[1]:
        return np.nan, 'response_outside_fitted_training_range'
    scaled = (value-bottom) / span
    if not 0 <= scaled < 1:
        return np.nan, 'inverse_domain_failure'
    prediction = midpoint * (-np.log1p(-scaled) / np.log(2)) ** (1/shape)
    if not np.isfinite(prediction) or not low <= prediction <= high:
        return np.nan, 'prediction_outside_training_range'
    return float(prediction), 'predicted'


def _validation_records(q):
    rows = []
    for scheme in ('held-out plate', 'held-out interior level'):
        units = sorted(q.physical_plate_id.unique()) if scheme == 'held-out plate' else range(1, 11)
        for unit in units:
            test_mask = q.physical_plate_id.eq(unit) if scheme == 'held-out plate' else q.dilution_rank.eq(unit)
            training = q[~test_mask]
            test = q[test_mask]
            shared = fit_ferric_response(training.axis_value, training['mean'])
            specific = {acid: fit_ferric_response(g.axis_value, g['mean'])
                        for acid, g in training.groupby('acid_M') if acid in test.acid_M.unique()}
            for row in test.itertuples():
                acid_train = training[training.acid_M.eq(row.acid_M)]
                low, high = float(acid_train.axis_value.min()), float(acid_train.axis_value.max())
                supported_reference = low <= row.axis_value <= high
                # The same acid-specific concentration support is used for both models.
                for model, fit in [('shared', shared), ('acid-specific', specific[row.acid_M])]:
                    predicted_A = float(fit.predict([row.axis_value])[0]) if supported_reference and _valid_fit(fit) else np.nan
                    if supported_reference:
                        predicted_C, status = _inverse(row.mean, fit, low, high)
                    else:
                        predicted_C, status = np.nan, 'reference_outside_training_range'
                    rows.append(dict(scheme=scheme, fold=str(unit), acid_M=row.acid_M,
                                     physical_plate_id=row.physical_plate_id,
                                     observation_id=row.observation_id, dilution_rank=row.dilution_rank,
                                     reference_mM=row.axis_value, observed_A380=row.mean,
                                     model=model, training_plate_n=training.physical_plate_id.nunique()
                                     if model == 'shared' else acid_train.physical_plate_id.nunique(),
                                     fit_converged=_valid_fit(fit), fit_identified=fit.identified,
                                     training_low_mM=low, training_high_mM=high,
                                     reference_in_training_range=supported_reference,
                                     predicted_A380=predicted_A,
                                     response_error_AU=predicted_A-row.mean,
                                     predicted_mM=predicted_C,
                                     concentration_error_mM=predicted_C-row.axis_value,
                                     prediction_status=status))
    return pd.DataFrame(rows)


def _inverse_table(records):
    rows = []
    strata = []
    for (scheme, acid), full in records[records.model.eq('acid-specific')].groupby(['scheme', 'acid_M']):
        strata.append((scheme, acid, 'all', full))
        for label, low, high in [('0–1', 0, 2), ('>1–10', 3, 5),
                                 ('>10–100', 6, 8), ('>100–1000', 9, 11)]:
            strata.append((scheme, acid, label, full[full.dilution_rank.between(low, high)]))
    for scheme, acid, window, g in strata:
        finite = np.isfinite(g.concentration_error_mM)
        e = g.loc[finite, 'concentration_error_mM']
        plate_rmse = g[finite].groupby('physical_plate_id').concentration_error_mM.apply(lambda a: np.sqrt(np.mean(a*a)))
        relative = g.loc[finite & g.reference_mM.gt(0)].copy()
        relative_error = 100*relative.concentration_error_mM/relative.reference_mM
        rows.append(dict(acid_M=acid, validation=scheme, nominal_reference_window_mM=window,
                         reference_min_mM=float(g.reference_mM.min()), reference_max_mM=float(g.reference_mM.max()),
                         plate_n=g.physical_plate_id.nunique(),
                         attempted_n=len(g), predicted_n=int(finite.sum()),
                         failed_n=int((~finite).sum()),
                         outside_reference_support_n=int((~g.reference_in_training_range).sum()),
                         response_or_inverse_failure_n=int((~finite & g.reference_in_training_range).sum()),
                         converged_unidentified_fit_n=int((g.fit_converged & ~g.fit_identified).sum()),
                         bias_mM=float(e.mean()), rmse_mM=float(np.sqrt(np.mean(e*e))),
                         median_absolute_error_mM=float(np.median(np.abs(e))) if len(e) else np.nan,
                         median_absolute_relative_error_percent=float(np.median(np.abs(relative_error))) if len(relative_error) else np.nan,
                         relative_error_n=len(relative_error),
                         plate_rmse_min_mM=float(plate_rmse.min()), plate_rmse_max_mM=float(plate_rmse.max()),
                         uncertainty_scope='conditional_CV_descriptive_no_independent_fold_CI'))
    return _table('ferric_inverse_validation', 'Held-out ferric concentration recovery',
        'Figure 1f is refitted as A380=b+s[1−exp{−ln(2)(C/K)^h}], using four-well physical-plate means '
        'and recorded total-Fe references C (mM). Each complete physical plate is withheld in turn; '
        'the second scheme withholds each of the ten interior dilution levels from all physical plates '
        'simultaneously. Level matching follows source well columns, not slightly different individually measured reference '
        'concentrations. The two endpoint levels are never extrapolated in the level-withholding scheme. '
        'Inverse predictions outside the fitted response over the training concentration range are '
        'unavailable; held-out references outside that range are separately counted. Bias is predicted '
        'minus reference concentration, and error summaries describe successful predictions only. '
        'Median absolute relative error excludes zero reference. Failed values are never clipped to '
        'a boundary or assigned zero error. Alongside the overall result, prespecified nominal reference '
        'windows group corresponding source dilution ranks (0, 0.3, 1; 3, 6, 10; 30, 60, 100; '
        '300, 600, 1000 mM); actual recorded extrema are retained in the CSV. Windows describe local '
        'prediction performance and are not declared quantification limits. Overlapping training folds are not independent replicates; '
        'these descriptive cross-validation errors have no fold-based confidence interval or P value. '
        'The table reports the physical-plate RMSE range. Concentration recovery is conditional on '
        'the measured reference values and excludes uncertainty from independent standard preparation.',
        pd.DataFrame(rows), ['acid_M', 'validation', 'nominal_reference_window_mM', 'attempted_n', 'predicted_n', 'outside_reference_support_n',
                             'response_or_inverse_failure_n', 'bias_mM', 'rmse_mM', 'median_absolute_relative_error_percent'],
        ['calibration_level_means'], ['1f'],
        ['The tested ferric range is not automatically an accurately invertible range.',
         'Conditional concentration recovery includes explicit inverse failures.'],
        ['bias_mM', 'rmse_mM', 'failed_n', 'median_absolute_relative_error_percent'])


def _predictive_table(records):
    rows = []
    for (scheme, acid), group in records.groupby(['scheme', 'acid_M']):
        wide = group.pivot(index=['observation_id', 'dilution_rank', 'reference_mM'], columns='model', values='response_error_AU')
        paired = wide.dropna(subset=['shared', 'acid-specific']).reset_index()
        r = dict(acid_M=acid, validation=scheme, plate_n=group.physical_plate_id.nunique(),
                 attempted_n=group.observation_id.nunique(), paired_prediction_n=len(paired))
        for model, prefix in [('shared', 'shared'), ('acid-specific', 'acid_specific')]:
            e = paired[model].to_numpy()
            by_level = paired.groupby('dilution_rank')[model].mean()
            r[prefix+'_rmse_AU'] = float(np.sqrt(np.mean(e*e))) if len(e) else np.nan
            r[prefix+'_bias_AU'] = float(np.mean(e)) if len(e) else np.nan
            r[prefix+'_max_absolute_level_bias_AU'] = float(by_level.abs().max())
            r[prefix+'_residual_reference_rho'] = (float(spearmanr(paired.reference_mM, e).statistic)
                                                  if len(e)>2 and np.ptp(e)>0 else np.nan)
        r['rmse_change_AU'] = r['acid_specific_rmse_AU']-r['shared_rmse_AU']
        r['squared_error_reduction_percent'] = 100*(1-r['acid_specific_rmse_AU']**2/r['shared_rmse_AU']**2)
        rows.append(r)
    return _table('ferric_predictive_model_comparison', 'Predictive benefit of acid-specific calibration',
        'The shared and four acid-specific models use the Figure 1f stretched-exponential '
        'equation with respectively four and sixteen total parameters. Models are refitted after withholding '
        'one physical plate or one common interior dilution level, as in the concentration-recovery analysis. '
        'The shared model can use other acid backgrounds; the corresponding acid-specific model uses only '
        'the training plates in that acid background. Both are evaluated on exactly the same supported '
        'held-out four-well means. Forward response error is predicted minus measured A380. Negative RMSE '
        'change means lower held-out error for acid-specific calibration; squared-error reduction may be '
        'negative. Maximum absolute level bias averages residuals over each nominal dilution level before '
        'taking the largest absolute value. It diagnoses localized misfit that overall RMSE can conceal. '
        'Machine-readable exports also retain residual–reference Spearman correlations as descriptive '
        'checks, without treating shared-fold residuals as independent tests. No training-SSE improvement '
        'is presented as predictive evidence. These errors are conditional on recorded references and '
        'the observed calibration plates, with no independent-fold confidence interval or P value.',
        pd.DataFrame(rows), ['acid_M', 'validation', 'paired_prediction_n', 'shared_rmse_AU',
                             'acid_specific_rmse_AU', 'rmse_change_AU', 'shared_max_absolute_level_bias_AU',
                             'acid_specific_max_absolute_level_bias_AU'],
        ['calibration_level_means'], ['1f'],
        ['Acid-specific training fit improvement must be distinguished from held-out predictive benefit.',
         'Held-out residual bias varies over the calibration range.'],
        ['rmse_change_AU', 'squared_error_reduction_percent', 'acid_specific_max_absolute_level_bias_AU'])


def _calibration_intervals(tables):
    rows = []
    groups = []
    ph = tables['ph_replicates'].copy()
    for (dose, signal, plate), g in ph.groupby(['nominal_dye_uM', 'signal', 'physical_plate_id']):
        # Multiple source preparation levels can have the same plotted mean pH.
        # Merge equal recorded x before defining an interval, retaining pooled within-level SD.
        collapsed = []
        for x, h in g.groupby('mean_measured_pH'):
            finite = h[np.isfinite(h['mean'])]
            collapsed.append(dict(x=x, y=float(finite['mean'].mean()),
                                  sd=float(np.sqrt(np.mean(finite.technical_sd**2))),
                                  source_levels=len(h)))
        groups.append((f'{dose:g} µM dye', SIGNALS[signal], 'pH', 'pH', plate,
                       pd.DataFrame(collapsed).sort_values('x')))
    fe = _ferric(tables)
    for (acid, plate), g in fe.groupby(['acid_M', 'physical_plate_id']):
        groups.append((f'{acid:g} M H2SO4', 'A380', 'iron', 'mM Fe', plate,
                       g.rename(columns={'axis_value': 'x', 'mean': 'y', 'technical_sd': 'sd'}).sort_values('x')))
    for cohort, readout, assay, unit, plate, g in groups:
        values = g[['x', 'y', 'sd']].to_numpy()
        for rank, (low, high) in enumerate(zip(values[:-1], values[1:])):
            change = high[1]-low[1]
            noise = np.sqrt(high[2]**2+low[2]**2)
            rows.append(dict(family='initial calibration', cohort=cohort, readout=readout,
                             assay=assay, reference_unit=unit, physical_plate_id=plate,
                             interval_rank=rank, reference_low=low[0], reference_high=high[0],
                             response_change=change, secant_sensitivity=change/(high[0]-low[0]),
                             change_over_noise=abs(change)/noise if noise>0 else np.nan))
    return rows


def _local_table(tables):
    rows = _calibration_intervals(tables)
    q = tables['pub_media_interval_diagnostics'].copy()
    # The upstream table already computes adjacent means with the same physical plate.
    for _, cohort in q.groupby(['medium_id', 'signal']):
        grids = [g.sort_values('predictor_low')[['predictor_low', 'predictor_high']].to_numpy()
                 for _, g in cohort.groupby('physical_plate_id')]
        if any(grid.shape != grids[0].shape or not np.allclose(grid, grids[0], rtol=1e-12, atol=1e-12)
               for grid in grids[1:]):
            raise ValueError('Media local-resolution intervals require a common ordered reference grid.')
    for (medium, signal, plate), g in q.groupby(['medium_id', 'signal', 'physical_plate_id']):
        g = g.sort_values('predictor_low')
        for rank, r in enumerate(g.itertuples()):
            rows.append(dict(family='medium validation', cohort=MEDIA.get(medium, medium),
                             readout=SIGNALS[signal], assay=r.assay,
                             reference_unit='pH' if r.assay=='pH' else 'mM added Fe',
                             physical_plate_id=plate, interval_rank=rank,
                             reference_low=r.predictor_low, reference_high=r.predictor_high,
                             response_change=r.observed_response_change,
                             secant_sensitivity=r.observed_secant_sensitivity,
                             change_over_noise=r.change_to_combined_well_sd))
    q = pd.DataFrame(rows)
    if q.duplicated(['family', 'cohort', 'readout', 'physical_plate_id', 'interval_rank']).any():
        raise ValueError('Local interval analysis accidentally duplicates a physical plate.')
    summaries = []
    for keys, g in q.groupby(['family', 'cohort', 'readout', 'assay', 'reference_unit', 'interval_rank']):
        family, cohort, readout, assay, unit, rank = keys
        row = dict(family=family, cohort=cohort, readout=readout, assay=assay, reference_unit=unit,
                   interval_rank=rank, scheduled_plate_n=g.physical_plate_id.nunique(),
                   reference_low=float(g.reference_low.mean()), reference_high=float(g.reference_high.mean()),
                   reference_low_min=float(g.reference_low.min()), reference_low_max=float(g.reference_low.max()),
                   reference_high_min=float(g.reference_high.min()), reference_high_max=float(g.reference_high.max()))
        row['reference_interval'] = f'{row["reference_low"]:.3g}–{row["reference_high"]:.3g}'
        for name, positive in [('response_change', False), ('secant_sensitivity', False), ('change_over_noise', True)]:
            n, estimate, low, high = _mean_ci(g[name], positive)
            row[name+'_n'] = n
            row[name] = estimate
            row[name+'_ci_low'], row[name+'_ci_high'] = low, high
            row[name+'_ci'] = _ci(estimate, low, high)
        direction = -1 if assay=='pH' else 1
        finite = np.isfinite(g.response_change)
        row['opposite_direction_plate_n'] = int((direction*g.loc[finite, 'response_change']<0).sum())
        row['interval_status'] = 'estimated' if row['response_change_n']>1 else ('single_plate' if row['response_change_n'] else 'unavailable')
        row['formal_detection_or_quantification_limit'] = False
        summaries.append(row)
    result = _table('local_response_resolution', 'Local response discrimination across calibration intervals',
        'Each row is an analysis of an adjacent reference interval, not an inventory of standards. '
        'Within each physical plate, response change is the higher-reference mean minus the lower-reference '
        'mean; secant sensitivity divides that change by the recorded reference interval. The diagnostic '
        'change/noise is the absolute change divided by sqrt(SDlow²+SDhigh²), using within-level well SDs '
        'rather than standard errors. This is not a formal detection or quantification limit. Signed change '
        'and sensitivity are averaged over physical plates with two-sided 95% Student-t intervals; '
        'change/noise uses a geometric mean and log-scale interval over positive finite ratios. Sample '
        'sizes and full-precision numeric limits for each metric are retained in the CSV. Intervals are '
        'pointwise descriptive intervals, without a multiple-testing claim that a complete range is validated. '
        'Opposite-direction counts identify increasing pH-ratio intervals or decreasing iron-response '
        'intervals; an interval with no eligible complete readings remains unavailable. The initial pH '
        'analysis uses the mean recorded pH plotted in Figure 1b–c; equal-x preparation levels are merged '
        'within plate before defining intervals, pooling their within-level SDs. Ferric initial intervals '
        'match the common source dilution rank and report mean recorded concentration bounds; original '
        'bounds for each plate remain in the source tables and their extrema in the CSV. These results '
        'are conditional on the recorded reference values and do not propagate the unquantified uncertainty or covariance of measured reference concentrations. '
        'Initial iron and medium results apply to Figure 1f and Figure 4a–c respectively.',
        pd.DataFrame(summaries), ['family', 'cohort', 'readout', 'reference_unit', 'reference_interval',
                                  'response_change_n', 'response_change_ci', 'secant_sensitivity_ci',
                                  'change_over_noise_ci', 'opposite_direction_plate_n'],
        ['ph_replicates', 'calibration_level_means', 'pub_media_interval_diagnostics'],
        ['1b', '1c', '1f', '4a', '4b', '4c'],
        ['Local discrimination and reversals vary across the tested range.',
         'Overall response span does not establish resolution at every reference interval.'],
        ['response_change', 'secant_sensitivity', 'change_over_noise', 'opposite_direction_plate_n'])
    from .publication_ph_support import summarize_local_resolution
    result['support_tables'] = {'physical_plate_intervals': q,
                               'cohort_summary': summarize_local_resolution(result['frame'])}
    result['caption'] += (' A separate cohort summary retains every cohort/readout, including unavailable '
                          'ones, and reports tested bounds, evaluable and reversed interval counts and the '
                          'minimum finite geometric separation with its reference interval. This selected '
                          'minimum is descriptive; it is not a global interval or pass criterion. Medium '
                          'iron intervals apply to 0.5 M nominal added sulfuric acid after mixing.')
    return result


def build(tables):
    records = _validation_records(_ferric(tables))
    inverse = _inverse_table(records)
    inverse['support_tables'] = {'held_out_predictions': records}
    return [inverse, _predictive_table(records), _local_table(tables)]
