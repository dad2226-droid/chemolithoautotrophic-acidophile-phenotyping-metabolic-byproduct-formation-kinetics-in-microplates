"""Storage support."""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy.stats import t as student_t


PLATE_KEYS = ['assay', 'physical_plate_id', 'condition']


LEVEL_KEYS = ['row_group', 'column']


SIGNAL_PAIRS = {'pH': ('ratio_raw', 'ratio_corrected'),
                'iron': ('A380', 'A380_corrected')}


SUPPORT_BASES = ('available_case', 'fixed_common_levels')


REASONS = ('missing_raw_acquisition', 'missing_corrected_acquisition',
           'nonfinite_raw_signal_or_initial_reference',
           'nonfinite_corrected_signal_or_initial_reference',
           'invalid_raw_original_initial_span',
           'invalid_corrected_original_initial_span',
           'nonfinite_raw_normalized_drift', 'nonfinite_corrected_normalized_drift')


def _finite_mean(values):
    a = np.asarray(values, float)
    a = a[np.isfinite(a)]
    return float(a.mean()) if len(a) else np.nan


def _finite_max(values):
    a = np.asarray(values, float)
    a = a[np.isfinite(a)]
    return float(a.max()) if len(a) else np.nan


def _level_eligibility(group, scheduled_indices, raw_signal, corrected_signal):
    """One row per physical reference level, with explicit dropout precedence."""
    result = []
    for identity, level in group.groupby(LEVEL_KEYS, dropna=False, sort=True):
        row = dict(zip(LEVEL_KEYS, identity))
        flags = {}
        for role, signal in [('raw', raw_signal), ('corrected', corrected_signal)]:
            observations = level[level.signal.eq(signal)]
            row[role + '_acquisition_n'] = len(observations)
            missing = set(scheduled_indices) - set(observations.read_index)
            row[role + '_missing_acquisition_n'] = len(missing)
            row[role + '_missing_acquisition_indices'] = ';'.join(map(str, sorted(missing)))
            signal_good = (np.isfinite(observations['mean']) &
                           np.isfinite(observations.initial_mean))
            span_good = (np.isfinite(observations.initial_response_range) &
                         observations.initial_response_range.gt(0))
            drift_good = np.isfinite(observations.normalized_difference)
            row[role + '_unavailable_signal_acquisition_n'] = int((~signal_good).sum())
            row[role + '_unavailable_drift_acquisition_n'] = int((~drift_good).sum())
            flags['missing_' + role + '_acquisition'] = bool(missing)
            flags['nonfinite_' + role + '_signal_or_initial_reference'] = bool((~signal_good).any())
            flags['invalid_' + role + '_original_initial_span'] = bool((~span_good).any())
            flags['nonfinite_' + role + '_normalized_drift'] = bool((~drift_good).any())
        row.update({'failure_' + reason: flags[reason] for reason in REASONS})
        row['fixed_common_level'] = not any(flags.values())
        row['primary_dropout_reason'] = next((r for r in REASONS if flags[r]), 'retained')
        result.append(row)
    return pd.DataFrame(result)


def _plate_metrics(acquisitions, selected, identity, basis, signal, span):
    post = acquisitions[acquisitions.read_index.gt(0)].sort_values('read_index')
    latest = acquisitions.sort_values('read_index').iloc[-1]
    drift = selected.loc[selected.read_index.gt(0), 'drift_percent_initial_span'].to_numpy(float)
    finite = drift[np.isfinite(drift)]
    signs = np.sign(post.mean_drift_percent_initial_span.to_numpy(float))
    adjacent = (np.diff(post.read_index.to_numpy()) == 1)
    eligible_pairs = adjacent & np.isfinite(signs[:-1]) & np.isfinite(signs[1:])
    row = dict(identity, signal=signal, support_basis=basis,
               original_initial_span=span,
               endpoint_read_index=int(latest.read_index), endpoint_day=latest.days,
               endpoint_drift_percent_initial_span=(latest.mean_drift_percent_initial_span
                                                    if latest.read_index > 0 else np.nan),
               endpoint_eligible_level_n=(int(latest.eligible_level_n) if latest.read_index > 0 else 0),
               endpoint_attempted_level_n=int(latest.attempted_level_n),
               post_initial_acquisition_n=len(post),
               eligible_post_initial_acquisition_n=int(np.isfinite(post.mean_drift_percent_initial_span).sum()),
               attempted_post_initial_level_acquisition_n=int(post.attempted_level_n.sum()),
               eligible_post_initial_level_acquisition_n=int(post.eligible_level_n.sum()),
               unavailable_post_initial_level_acquisition_n=int((post.attempted_level_n-post.eligible_level_n).sum()),
               mean_observed_drift=_finite_mean(post.mean_drift_percent_initial_span),
               max_absolute_plate_drift=_finite_max(np.abs(post.mean_drift_percent_initial_span)),
               max_absolute_level_drift=_finite_max(np.abs(finite)),
               rms_observed_level_drift=float(np.sqrt(np.mean(finite ** 2))) if len(finite) else np.nan,
               consecutive_pair_n=int(eligible_pairs.sum()),
               same_direction_pair_n=int(((signs[:-1]*signs[1:] > 0) & eligible_pairs).sum()))
    row['endpoint_status'] = ('estimated' if np.isfinite(row['endpoint_drift_percent_initial_span'])
                              else 'no_eligible_level_at_last_acquired_index')
    return row


def _summary(metrics):
    rows = []
    group_keys = ['assay', 'condition', 'signal', 'support_basis']
    fields = ['endpoint_drift_percent_initial_span',
              'max_absolute_plate_drift', 'max_absolute_level_drift']
    for identity, group in metrics.groupby(group_keys, dropna=False, sort=True):
        for field in fields:
            values = group[field].to_numpy(float)
            values = values[np.isfinite(values)]
            n = len(values)
            estimate = float(values.mean()) if n else np.nan
            lower = upper = np.nan
            if field == 'endpoint_drift_percent_initial_span':
                interval_type = '95_percent_Student_t_CI_across_physical_plates'
                if n >= 2:
                    margin = student_t.ppf(.975, n-1)*values.std(ddof=1)/np.sqrt(n)
                    lower, upper = estimate-margin, estimate+margin
                status = 'estimated' if n >= 2 else ('single_plate_no_interval' if n else 'unavailable')
            else:
                interval_type = 'observed_physical_plate_range_not_confidence_interval'
                if n:
                    lower, upper = float(values.min()), float(values.max())
                status = 'descriptive' if n else 'unavailable'
            row = dict(zip(group_keys, identity))
            row.update(metric=field, unit='percent_original_initial_span',
                       attempted_plate_n=len(group), eligible_plate_n=n, unavailable_plate_n=len(group)-n,
                       endpoint_day_mean=_finite_mean(group.endpoint_day),
                       endpoint_day_min=group.endpoint_day.min(), endpoint_day_max=group.endpoint_day.max(),
                       endpoint_read_index_min=int(group.endpoint_read_index.min()),
                       endpoint_read_index_max=int(group.endpoint_read_index.max()),
                       plates_with_acquisition_gaps_n=int(group.absent_index_n.gt(0).sum()),
                       absent_index_n=int(group.absent_index_n.sum()),
                       observed_schedule_scope='acquired_indices_only_no_interpolation_or_missing_index_bridging',
                       estimate=estimate, interval_type=interval_type, interval_low=lower,
                       interval_high=upper, status=status)
            for count in ['attempted_reference_level_n', 'fixed_common_level_n',
                          'endpoint_attempted_level_n', 'endpoint_eligible_level_n',
                          'post_initial_acquisition_n', 'eligible_post_initial_acquisition_n',
                          'attempted_post_initial_level_acquisition_n',
                          'eligible_post_initial_level_acquisition_n',
                          'unavailable_post_initial_level_acquisition_n']:
                row[count] = int(group[count].sum())
            rows.append(row)
    return pd.DataFrame(rows)


def storage_common_support(tables):
    """Return five reproducible support exports for existing storage analyses.

    ``level_eligibility`` gives mutually exclusive primary dropout reasons plus
    explicitly nonexclusive failure flags. ``physical_plate_eligibility`` counts
    reference levels once, irrespective of readout or repeated acquisitions.
    ``acquisition_means`` and ``physical_plate_metrics`` show both original
    available-case and fixed-common-level drift. ``drift_summary`` uses only
    physical plates for endpoint CIs; selected maxima have observed ranges only.

    Missingness can be diagnosed only to the level represented in this aggregated
    source: absent acquisition rows, nonfinite signals/references, invalid initial
    spans, or nonfinite normalized drift. It does not infer instrument-flag or
    ratio-denominator causes absent from these rows.
    """
    q = tables['storage_folds'].copy()
    required = PLATE_KEYS + LEVEL_KEYS + ['signal', 'read_index', 'days', 'mean',
                                         'initial_mean', 'initial_response_range', 'normalized_difference']
    missing = sorted(set(required)-set(q.columns))
    if missing:
        raise ValueError('Storage support is missing columns: ' + ', '.join(missing))
    if q.empty:
        raise ValueError('Storage support requires acquired storage observations')
    if q[PLATE_KEYS + LEVEL_KEYS].isna().any().any():
        raise ValueError('Storage physical-plate and level identities must be present')
    if q.duplicated(PLATE_KEYS + LEVEL_KEYS + ['signal', 'read_index']).any():
        raise ValueError('Duplicate storage level acquisition')
    if (~np.isfinite(q.read_index) | q.read_index.lt(0) | q.read_index.mod(1).ne(0)).any():
        raise ValueError('Storage acquisition indices must be nonnegative integers')
    q['read_index'] = q.read_index.astype(int)
    q['drift_percent_initial_span'] = 100*q.normalized_difference
    level_frames, eligibility_rows, acquisition_rows, metric_rows = [], [], [], []
    for key, group in q.groupby(PLATE_KEYS, dropna=False, sort=True):
        identity = dict(zip(PLATE_KEYS, key))
        if identity['assay'] not in SIGNAL_PAIRS:
            raise ValueError('Unsupported storage assay: ' + str(identity['assay']))
        signals = SIGNAL_PAIRS[identity['assay']]
        if not set(group.signal).issubset(signals):
            raise ValueError('Unexpected storage readout for assay: ' + str(identity['assay']))
        acquired_indices = sorted(group.read_index.unique())
        scheduled_indices = sorted(set(acquired_indices) | {0})
        absent_indices = sorted(set(range(max(acquired_indices)+1))-set(acquired_indices))
        eligibility = _level_eligibility(group, scheduled_indices, *signals)
        for name, value in identity.items():
            eligibility[name] = value
        level_frames.append(eligibility)
        plate_eligibility = dict(identity,
            attempted_reference_level_n=len(eligibility), fixed_common_level_n=int(eligibility.fixed_common_level.sum()),
            excluded_reference_level_n=int((~eligibility.fixed_common_level).sum()),
            acquired_index_n=len(acquired_indices), acquired_indices=';'.join(map(str, acquired_indices)),
            absent_index_n=len(absent_indices), absent_indices=';'.join(map(str, absent_indices)),
            required_index_n=len(scheduled_indices),
            fixed_support_status='estimated' if eligibility.fixed_common_level.any() else 'no_common_levels')
        for reason in REASONS:
            plate_eligibility['primary_' + reason + '_level_n'] = int(eligibility.primary_dropout_reason.eq(reason).sum())
            plate_eligibility['nonexclusive_' + reason + '_level_n'] = int(eligibility['failure_' + reason].sum())
        eligibility_rows.append(plate_eligibility)
        annotated = group.merge(eligibility[LEVEL_KEYS + ['fixed_common_level']],
                                on=LEVEL_KEYS, validate='many_to_one')
        for signal in signals:
            source = annotated[annotated.signal.eq(signal)]
            spans = source.initial_response_range.drop_duplicates()
            if len(spans) > 1:
                raise ValueError('Original initial span changed within one plate/readout')
            span = float(spans.iloc[0]) if len(spans) else np.nan
            for basis in SUPPORT_BASES:
                selected = source if basis == 'available_case' else source[source.fixed_common_level]
                this_acquisition_rows = []
                for index in acquired_indices:
                    current = selected[selected.read_index.eq(index)]
                    values = current.drift_percent_initial_span.to_numpy(float)
                    days = group.loc[group.read_index.eq(index), 'days']
                    row = dict(identity, signal=signal, support_basis=basis, read_index=index,
                        days=_finite_mean(days), attempted_level_n=len(eligibility),
                        acquired_level_observation_n=int(source.read_index.eq(index).sum()),
                        selected_level_observation_n=len(current), eligible_level_n=int(np.isfinite(values).sum()),
                        fixed_common_level_n=plate_eligibility['fixed_common_level_n'],
                        mean_drift_percent_initial_span=_finite_mean(values), original_initial_span=span)
                    acquisition_rows.append(row)
                    this_acquisition_rows.append(row)
                metrics = _plate_metrics(pd.DataFrame(this_acquisition_rows), selected,
                                         identity, basis, signal, span)
                metrics.update(attempted_reference_level_n=len(eligibility),
                               fixed_common_level_n=plate_eligibility['fixed_common_level_n'],
                               absent_index_n=len(absent_indices),
                               absent_indices=';'.join(map(str, absent_indices)))
                metric_rows.append(metrics)
    metrics = pd.DataFrame(metric_rows)
    return {'level_eligibility': pd.concat(level_frames, ignore_index=True),
            'physical_plate_eligibility': pd.DataFrame(eligibility_rows),
            'acquisition_means': pd.DataFrame(acquisition_rows),
            'physical_plate_metrics': metrics,
            'drift_summary': _summary(metrics)}
