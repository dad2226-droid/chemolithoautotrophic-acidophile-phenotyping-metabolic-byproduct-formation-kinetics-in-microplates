"""Lf benchmarks."""
from __future__ import annotations
import json
import numpy as np
import pandas as pd


def at_exact_endpoint(tables, target_day_d):
    """Interpolate adjacent recorded acquisitions; never extrapolate or bridge NA."""
    source = tables['af_iron_control_plate_trajectories'].copy()
    source['source_row'] = np.arange(1, len(source) + 1)
    source = source[source.signal.eq('raw_A380') & source.control_cohort_definition.map(
        lambda value: json.loads(value).get('nominal_initial_pH_json') == 'null')]
    units = []
    for (batch, plate), group in source.groupby(['batch_id', 'plate_id'], sort=True):
        group = group.sort_values('days_mean', kind='stable')
        before = group[group.days_mean.le(target_day_d)]
        after = group[group.days_mean.ge(target_day_d)]
        row = dict(batch_id=batch, plate_id=plate, target_day_d=target_day_d,
                   source_table='af_iron_control_plate_trajectories',
                   benchmark_value=np.nan, eligible=False, status='outside_recorded_support')
        if len(before) and len(after):
            left, right = before.iloc[-1], after.iloc[0]
            row.update(left_source_row=int(left.source_row), right_source_row=int(right.source_row),
                       left_time_d=left.days_mean, right_time_d=right.days_mean,
                       left_value=left.biotic_change, right_value=right.biotic_change,
                       status='missing_bracketing_measurement')
            if np.isfinite([left.biotic_change, right.biotic_change]).all():
                weight = ((target_day_d-left.days_mean)/(right.days_mean-left.days_mean)
                          if right.days_mean > left.days_mean else 0.)
                row.update(interpolation_weight=weight,
                           benchmark_value=(1-weight)*left.biotic_change+weight*right.biotic_change,
                           eligible=True, status='linear_interpolation_of_adjacent_acquisitions'
                           if weight else 'exact_recorded_acquisition')
        units.append(row)
    units = pd.DataFrame(units)
    batches = units[units.eligible].groupby('batch_id', as_index=False).agg(
        mean=('benchmark_value', 'mean'), within_batch_plate_sd=('benchmark_value', 'std'),
        physical_plate_n=('benchmark_value', 'count'))
    batches['target_day_d'] = target_day_d
    summary = pd.DataFrame([dict(target_day_d=target_day_d, mean=batches['mean'].mean(),
        sd=batches['mean'].std(), batch_n=len(batches), physical_plate_n=int(units.eligible.sum()),
        candidate_plate_n=len(units), candidate_batch_n=units.batch_id.nunique(),
        selection='adjacent_acquisition_linear_interpolation_no_extrapolation_no_missing_bridge',
        aggregation='equal_plate_weight_within_batch_equal_batch_mean_and_between_batch_sample_SD')])
    return units, batches, summary
