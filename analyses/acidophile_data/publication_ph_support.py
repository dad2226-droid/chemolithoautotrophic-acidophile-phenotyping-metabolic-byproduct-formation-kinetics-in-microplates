"""Ph support."""
from __future__ import annotations
import json
import numpy as np
import pandas as pd


_UNIT = ['dataset', 'plate_id', 'layout_group', 'column']


_STRATUM = ['comparison_acid_m', 'comparison_matrix']


_LOCAL = ['family', 'cohort', 'readout', 'assay', 'reference_unit']


def _require(frame, columns, name):
    missing = set(columns) - set(frame.columns)
    if missing:
        raise ValueError(f'{name} lacks required columns: {sorted(missing)}')


def _json_values(values):
    return json.dumps(sorted(set(values.dropna().astype(str))), ensure_ascii=False)


def _finite_bounds(values):
    values = pd.to_numeric(values, errors='raise')
    finite = values[np.isfinite(values)]
    return (float(finite.min()), float(finite.max())) if len(finite) else (np.nan, np.nan)


def _coverage(records):
    """Count physical reference observations, never resampling folds or wells."""
    eligible = records[records.common_reference_support]
    predicted = records[records.prediction_available]
    low, high = _finite_bounds(eligible.reference_ph)
    predicted_low, predicted_high = _finite_bounds(predicted.reference_ph)
    bins = sorted(eligible.pH_bin.unique().tolist())
    predicted_bins = sorted(predicted.pH_bin.unique().tolist())
    return dict(
        support_attempted_physical_reference_n=len(records),
        support_eligible_physical_reference_n=len(eligible),
        support_prediction_physical_reference_n=len(predicted),
        support_failed_physical_reference_n=len(eligible)-len(predicted),
        support_reference_pH_min=low,
        support_reference_pH_max=high,
        support_unique_reference_pH_n=int(eligible.reference_ph.nunique()),
        support_prediction_reference_pH_min=predicted_low,
        support_prediction_reference_pH_max=predicted_high,
        support_prediction_unique_reference_pH_n=int(predicted.reference_ph.nunique()),
        support_shared_bin_n=len(bins),
        support_prediction_bin_n=len(predicted_bins),
        support_failed_entire_bin_n=len(set(bins)-set(predicted_bins)),
        support_bin_ids=json.dumps([int(x) for x in bins]),
        support_bin_width_pH=float(records.pH_bin_width.iloc[0]) if len(records) else np.nan,
        support_complete_bin_coverage=bool(bins and bins == predicted_bins),
        support_attempted_dataset_ids=_json_values(records.dataset),
        support_eligible_dataset_ids=_json_values(eligible.dataset),
        support_attempted_plate_ids=_json_values(records.plate_id),
        support_eligible_plate_ids=_json_values(eligible.plate_id),
        support_prediction_plate_ids=_json_values(predicted.plate_id),
        support_attempted_source_files=_json_values(records.source_file),
        support_eligible_source_files=_json_values(eligible.source_file),
        support_eligible_plate_n=int(eligible.plate_id.nunique()),
        support_prediction_plate_n=int(predicted.plate_id.nunique()),
        support_status=('no_common_reference_support' if not len(eligible) else
                        'no_available_prediction' if not len(predicted) else
                        'failed_entire_common_bin' if bins != predicted_bins else
                        'all_common_bins_represented'),
        support_count_unit='physical_plate_preparation_level; repeated_fits_are_not_replicates',
    )


def _append(frame, records):
    additional = pd.DataFrame(records, index=frame.index)
    overlap = set(frame.columns) & set(additional.columns)
    if overlap:
        raise ValueError(f'Coverage would overwrite result columns: {sorted(overlap)}')
    return pd.concat([frame.copy(deep=True), additional], axis=1)


def _check_count(row, column, actual):
    if column in row and int(row[column]) != int(actual):
        raise ValueError(f'pH transfer {column} disagrees with original support: '
                         f'{row[column]} != {actual}')


def _reference_records(group, population, acid, matrix, signal):
    keep = _UNIT + ['source_file', 'reference_ph', 'pH_bin', 'pH_bin_width',
                    'common_reference_support', 'prediction_available']
    q = group[keep].copy()
    q.insert(0, 'support_population', population)
    q.insert(1, 'comparison_acid_m', acid)
    q.insert(2, 'comparison_matrix', matrix)
    q.insert(3, 'signal', signal)
    return q


def augment_ph_transfer(accuracy, correction, common_support):
    """Summarize physical reference coverage for pH accuracy and correction.

    Accuracy uses the common-reference support for each readout. The correction
    comparison pairs measured and corrected ratios from the same physical unit
    and requires both readouts to meet the reference-support criteria. Bounds
    describe the pH values of the retained observations. Prediction errors
    must be finite to contribute to the accuracy statistics.
    """
    required = (_UNIT + _STRATUM + ['source_file', 'signal', 'reference_ph',
                'error_pH', 'pH_bin', 'pH_bin_width', 'common_reference_support'])
    _require(common_support, required, 'ph_matrix_common_support')
    _require(accuracy, ['acid_m', 'matrix', 'signal'], 'pH accuracy result')
    _require(correction, ['acid_m', 'matrix'], 'pH correction result')
    supported = common_support.copy(deep=True)
    if supported.duplicated(_STRATUM + ['signal'] + _UNIT).any():
        raise ValueError('pH support duplicates a physical plate/preparation unit.')
    if supported.common_reference_support.isna().any():
        raise ValueError('pH common-reference support flags must be explicit.')
    if len(supported) and not supported.pH_bin_width.eq(.25).all():
        raise ValueError('pH transfer coverage requires the existing 0.25-pH bins.')
    if (supported.error_pH.notna() & ~np.isfinite(supported.error_pH)).any():
        raise ValueError('pH errors must be finite or explicitly unavailable.')
    supported['common_reference_support'] = supported.common_reference_support.eq(True)
    training = supported[supported.comparison_matrix.eq('acid_held_out_plate')]
    provenance = {
        'support_training_dataset_ids': _json_values(training.dataset),
        'support_training_plate_ids': _json_values(training.plate_id),
        'support_training_source_files': _json_values(training.source_file),
        'support_training_scope': 'acid_cohort_pool; each_fit_withholds_one_acid_plate',
    }
    arows, crows, arefs, crefs = [], [], [], []
    for row in accuracy.to_dict('records'):
        group = supported[supported.comparison_acid_m.eq(row['acid_m']) &
                          supported.comparison_matrix.eq(row['matrix']) &
                          supported.signal.eq(row['signal'])].copy()
        group['prediction_available'] = group.common_reference_support & np.isfinite(group.error_pH)
        coverage = _coverage(group)
        for count_field, coverage_field in [('attempted_n', 'support_attempted_physical_reference_n'),
                         ('common_support_n', 'support_eligible_physical_reference_n'),
                         ('finite_prediction_n', 'support_prediction_physical_reference_n'),
                         ('inverse_failure_n', 'support_failed_physical_reference_n')]:
            _check_count(row, count_field, coverage[coverage_field])
        arows.append({**coverage, **provenance})
        arefs.append(_reference_records(group, 'single_readout', row['acid_m'], row['matrix'], row['signal']))
    for row in correction.to_dict('records'):
        group = supported[supported.comparison_acid_m.eq(row['acid_m']) &
                          supported.comparison_matrix.eq(row['matrix'])]
        raw = group[group.signal.eq('ratio_raw')].set_index(_UNIT)
        corr = group[group.signal.eq('ratio_corrected')].set_index(_UNIT)
        cols = ['source_file', 'reference_ph', 'pH_bin', 'pH_bin_width',
                'common_reference_support', 'error_pH']
        pairs = raw[cols].join(corr[cols], lsuffix='_raw', rsuffix='_corrected')
        intended = pairs.common_reference_support_raw.eq(True) & pairs.common_reference_support_corrected.eq(True)
        if (pairs.loc[intended, 'reference_ph_raw'].ne(pairs.loc[intended, 'reference_ph_corrected']).any()
                or pairs.loc[intended, 'pH_bin_raw'].ne(pairs.loc[intended, 'pH_bin_corrected']).any()):
            raise ValueError('Paired pH readouts do not share the same physical reference.')
        records = pairs.reset_index()[_UNIT].copy()
        for name in ['source_file', 'reference_ph', 'pH_bin', 'pH_bin_width']:
            records[name] = pairs[name+'_raw'].to_numpy()
        records['common_reference_support'] = intended.to_numpy()
        raw_valid = pairs.error_pH_raw.notna()
        corr_valid = pairs.error_pH_corrected.notna()
        records['prediction_available'] = (intended & raw_valid & corr_valid).to_numpy()
        coverage = _coverage(records)
        _check_count(row, 'attempted_pair_n', coverage['support_eligible_physical_reference_n'])
        _check_count(row, 'paired_n', coverage['support_prediction_physical_reference_n'])
        _check_count(row, 'raw_failure_n', int((intended & ~raw_valid).sum()))
        _check_count(row, 'corrected_failure_n', int((intended & ~corr_valid).sum()))
        coverage.update(
            support_raw_only_reference_n=int((pairs.common_reference_support_raw.eq(True) &
                                              ~pairs.common_reference_support_corrected.eq(True)).sum()),
            support_corrected_only_reference_n=int((~pairs.common_reference_support_raw.eq(True) &
                                                    pairs.common_reference_support_corrected.eq(True)).sum()),
            support_raw_prediction_failure_n=int((intended & ~raw_valid).sum()),
            support_corrected_prediction_failure_n=int((intended & ~corr_valid).sum()),
        )
        crows.append({**coverage, **provenance})
        refs = _reference_records(records, 'paired_raw_corrected', row['acid_m'], row['matrix'], 'corrected_minus_raw')
        refs['raw_common_reference_support'] = pairs.common_reference_support_raw.eq(True).to_numpy()
        refs['corrected_common_reference_support'] = pairs.common_reference_support_corrected.eq(True).to_numpy()
        refs['raw_prediction_available'] = raw_valid.to_numpy()
        refs['corrected_prediction_available'] = corr_valid.to_numpy()
        crefs.append(refs)
    a = _append(accuracy, arows)
    c = _append(correction, crows)
    exports = {
        'accuracy_support': a[['acid_m', 'matrix', 'signal'] + list(arows[0])].copy() if arows else a.copy(),
        'correction_support': c[['acid_m', 'matrix'] + list(crows[0])].copy() if crows else c.copy(),
        'accuracy_references': pd.concat(arefs, ignore_index=True) if arefs else pd.DataFrame(),
        'correction_references': pd.concat(crefs, ignore_index=True) if crefs else pd.DataFrame(),
    }
    return a, c, exports


def summarize_local_resolution(intervals):
    """Describe every original cohort/readout without selecting by significance.

    One row is returned per existing family/cohort/readout/assay/unit, including
    all-unavailable readouts. Change/noise and reversal eligibility are separate
    denominators. The minimum is over positive finite interval geometric means;
    ties use the original interval rank then numeric reference bounds. No CI for
    this selected minimum is promoted to a guarantee for the complete range.
    Actual extrema over recorded plate-specific bounds are retained separately
    from the selected interval's mean reference bounds.
    """
    required = _LOCAL + ['interval_rank', 'reference_low', 'reference_high',
        'reference_low_min', 'reference_low_max', 'reference_high_min', 'reference_high_max',
        'response_change_n', 'response_change', 'change_over_noise_n', 'change_over_noise',
        'opposite_direction_plate_n']
    _require(intervals, required, 'local response interval results')
    if intervals.duplicated(_LOCAL + ['interval_rank']).any():
        raise ValueError('Local resolution contains duplicate cohort intervals.')
    if ((intervals.opposite_direction_plate_n < 0) |
            (intervals.opposite_direction_plate_n > intervals.response_change_n)).any():
        raise ValueError('Local reversal count exceeds evaluated physical plates.')
    rows = []
    for identity, group in intervals.groupby(_LOCAL, dropna=False, sort=True):
        change_ok = group.response_change_n.gt(0) & np.isfinite(group.response_change)
        separation_ok = (group.change_over_noise_n.gt(0) & np.isfinite(group.change_over_noise)
                         & group.change_over_noise.gt(0))
        usable = group[separation_ok].sort_values(
            ['change_over_noise', 'interval_rank', 'reference_low', 'reference_high'], kind='stable')
        minimum = usable.iloc[0] if len(usable) else None
        result = dict(zip(_LOCAL, identity))
        result.update(
            reference_min=_finite_bounds(group.reference_low_min)[0],
            reference_max=_finite_bounds(group.reference_high_max)[1],
            total_interval_n=len(group),
            change_evaluable_interval_n=int(change_ok.sum()),
            change_unavailable_interval_n=int((~change_ok).sum()),
            separation_evaluable_interval_n=int(separation_ok.sum()),
            separation_unavailable_interval_n=int((~separation_ok).sum()),
            intervals_with_any_reversed_plate_n=int((change_ok & group.opposite_direction_plate_n.gt(0)).sum()),
            reversal_evaluable_interval_n=int(change_ok.sum()),
            evaluated_plate_interval_n=int(group.loc[change_ok, 'response_change_n'].sum()),
            reversed_plate_interval_n=int(group.loc[change_ok, 'opposite_direction_plate_n'].sum()),
            minimum_interval_geometric_change_over_noise=float(minimum.change_over_noise) if minimum is not None else np.nan,
            minimum_interval_rank=int(minimum.interval_rank) if minimum is not None else np.nan,
            minimum_interval_reference_low=float(minimum.reference_low) if minimum is not None else np.nan,
            minimum_interval_reference_high=float(minimum.reference_high) if minimum is not None else np.nan,
            minimum_interval_reference_low_min=float(minimum.reference_low_min) if minimum is not None else np.nan,
            minimum_interval_reference_low_max=float(minimum.reference_low_max) if minimum is not None else np.nan,
            minimum_interval_reference_high_min=float(minimum.reference_high_min) if minimum is not None else np.nan,
            minimum_interval_reference_high_max=float(minimum.reference_high_max) if minimum is not None else np.nan,
            minimum_interval_plate_n=int(minimum.change_over_noise_n) if minimum is not None else 0,
            summary_status='no_estimable_separation' if minimum is None else 'descriptive_interval_summary',
            reversal_denominator='intervals_with_at_least_one_finite_physical_plate_response_change',
            minimum_selection='minimum_positive_finite_interval_geometric_mean; ties_by_interval_rank',
            uncertainty_scope='selected_minimum_is_descriptive; no_overall_range_CI_or_detection_limit',
            formal_detection_or_quantification_limit=False,
        )
        rows.append(result)
    return pd.DataFrame(rows)
