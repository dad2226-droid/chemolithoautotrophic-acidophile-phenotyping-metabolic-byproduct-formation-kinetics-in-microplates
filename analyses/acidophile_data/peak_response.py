"""Peak response."""
from __future__ import annotations
import numpy as np
import pandas as pd
from .plate_comparisons import CONDITION_KEYS


def plate_means(frame: pd.DataFrame) -> pd.DataFrame:
    """One mean per physical plate, condition, signal and recorded acquisition.

Shared scientific memberships retain their condition identity but cannot add
duplicate wells within a condition. Missing/invalid wells are retained in the
counts; the primary response requires the complete assigned set at that read.
"""
    keys = CONDITION_KEYS + ['batch_id', 'plate_id', 'signal', 'ordinal']
    f = frame.copy()
    identity = 'physical_observation_key' if 'physical_observation_key' in f else 'well'
    f = f.drop_duplicates(keys + [identity])
    measures = [v for v in ('biotic_signal', 'abiotic_signal', 'biotic_change',
                            'abiotic_change', 'net_change') if v in f]
    spec = {v: (v, 'mean') for v in measures}
    spec.update({v + '_sd': (v, 'std') for v in measures})
    spec.update({v + '_n': (v, 'count') for v in measures})
    spec.update(days_mean=('days', 'mean'), days_min=('days', 'min'), days_max=('days', 'max'),
                well_n=(identity, 'nunique'),
                culture_ids=('culture_id', lambda s: ';'.join(sorted(set(s.dropna())))),
                sources=('source', lambda s: ';'.join(sorted(set(s.dropna())))))
    out = f.groupby(keys, dropna=False, as_index=False).agg(**spec)
    series = CONDITION_KEYS + ['batch_id', 'plate_id', 'signal']
    out['expected_well_n'] = out.groupby(series).well_n.transform('max')
    for v in measures:
        out[v + '_available_mean'] = out[v]
        out[v] = out[v].where(out[v + '_n'].eq(out.expected_well_n))
    out['aggregation'] = 'mean_of_recorded_wells_within_physical_plate_condition_acquisition'
    return out


def observed_peak(curve: pd.DataFrame, *, value='biotic_change', time='days_mean') -> dict:
    """Integrate from the recorded initial read to this curve's first maximum.

Only adjacent recorded, finite points before and including that maximum enter
the trapezoid sum. A missing initial value, duplicate time, or internal gap
withholds the integral. Missing data after the maximum do not affect it.
Negative increments/areas are never clipped. A non-rising curve has a zero
area at its initial maximum; its time of zero is not interpreted as a fast rise.
"""
    out = dict(auc_to_max=np.nan, initial_time_d=np.nan, maximum_time_d=np.nan,
               time_to_max_d=np.nan, maximum_response=np.nan, observed_end_d=np.nan,
               observation_n=len(curve), valid_observation_n=0, integrated_observation_n=0,
               candidate_premax_observation_n=0,
               excluded_post_max_n=0, maximum_tie_n=0, maximum_at_final_read=False,
               positive_rise=False, maximum_status='unavailable', status='no_observations')
    if curve.empty:
        return out
    g = curve.sort_values(time, kind='stable')
    x, y = g[time].to_numpy(float), g[value].to_numpy(float)
    finite = np.isfinite(x) & np.isfinite(y)
    out.update(initial_time_d=x[0], observed_end_d=x[-1], valid_observation_n=int(finite.sum()))
    if not np.isfinite(x).all() or len(np.unique(x)) != len(x):
        out['status'] = 'invalid_or_duplicate_recorded_time'
        return out
    if not finite.any():
        out['status'] = 'no_valid_response'
        return out
    maximum = float(np.max(y[finite]))
    indices = np.flatnonzero(finite & (y == maximum))
    peak = int(indices[0])
    out.update(maximum_time_d=float(x[peak]), time_to_max_d=float(x[peak] - x[0]),
               maximum_response=maximum, candidate_premax_observation_n=peak + 1,
               excluded_post_max_n=len(g) - peak - 1, maximum_tie_n=len(indices),
               maximum_at_final_read=peak == len(g) - 1,
               positive_rise=bool(finite[0] and maximum > y[0]),
               maximum_status='final_read_maximum_plateau_unconfirmed' if peak == len(g) - 1
               else 'first_of_tied_observed_maxima' if len(indices) > 1 else 'interior_observed_maximum')
    if not finite[0]:
        out['status'] = 'missing_initial_response'
        return out
    if not finite[:peak + 1].all():
        out['status'] = 'missing_recorded_value_before_maximum'
        return out
    if finite.sum() < 2:
        out['status'] = 'insufficient_acquisitions'
        return out
    out['auc_to_max'] = float(np.trapezoid(y[:peak + 1], x[:peak + 1])) if peak else 0.
    out['integrated_observation_n'] = peak + 1
    out['status'] = 'estimated' if out['positive_rise'] else 'no_positive_rise'
    if not out['positive_rise']:
        out['maximum_status'] = 'initial_maximum_no_positive_rise'
    return out


def compare_peaks(treated: dict, references: list[dict | None]) -> dict:
    """Normalize to the mean of separately integrated control replicates.

The denominator is never the area of an averaged control trajectory. Every
recorded reference is required, avoiding silent replacement by surviving ones.
"""
    out = {f'treated_{key}': value for key, value in treated.items()}
    out.update(control_integral=np.nan, control_integral_sd=np.nan, control_n=len(references),
               control_time_to_max_d=np.nan, control_time_to_max_sd_d=np.nan,
               treated_integral=treated['auc_to_max'], integral_inhibition_percent=np.nan,
               activity_percent=np.nan, time_to_max_delay_d=np.nan,
               status='missing_no_additive_control')
    if not references:
        return out
    if any(r is None for r in references):
        out['status'] = 'missing_assigned_control'
        return out
    times = np.array([r['time_to_max_d'] for r in references], float)
    if treated['positive_rise'] and all(r['positive_rise'] for r in references):
        out['time_to_max_delay_d'] = treated['time_to_max_d'] - times.mean()
    if any(not np.isfinite(r['auc_to_max']) for r in references):
        out['status'] = 'unavailable_reference_integral'
        return out
    areas = np.array([r['auc_to_max'] for r in references], float)
    out.update(control_integral=float(areas.mean()),
               control_integral_sd=float(areas.std(ddof=1)) if len(areas) > 1 else np.nan,
               control_time_to_max_d=float(times.mean()),
               control_time_to_max_sd_d=float(times.std(ddof=1)) if len(times) > 1 else np.nan,
               control_final_read_maximum_n=sum(r['maximum_at_final_read'] for r in references),
               control_maximum_times_d=';'.join(str(r['maximum_time_d']) for r in references))
    if areas.mean() <= 0:
        out['status'] = 'nonpositive_reference_integral'
        return out
    if not np.isfinite(treated['auc_to_max']):
        out['status'] = treated['status']
        return out
    out['activity_percent'] = 100 * treated['auc_to_max'] / areas.mean()
    out['integral_inhibition_percent'] = 100 - out['activity_percent']
    out['status'] = 'estimated'
    return out
