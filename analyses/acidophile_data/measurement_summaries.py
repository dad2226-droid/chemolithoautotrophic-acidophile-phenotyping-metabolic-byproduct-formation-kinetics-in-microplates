"""Measurement summaries."""
from __future__ import annotations
import json
import numpy as np
import pandas as pd


def storage_reference_ranges(storage):
    """Initial observed level span within each physical plate/readout, never pooled."""
    keys = ['assay', 'physical_plate_id', 'condition', 'signal']
    initial = storage[storage.read_index.eq(0)]
    return initial.groupby(keys, dropna=False, as_index=False).agg(
        initial_min=('mean', 'min'), initial_max=('mean', 'max'),
        initial_level_n=('mean', 'count')).assign(
        initial_response_range=lambda d: d.initial_max-d.initial_min)


def storage_drift(storage):
    """Subtract the same physical plate/level's first acquisition; retain missingness."""
    keys = ['assay', 'physical_plate_id', 'condition', 'signal', 'row_group', 'column']
    initial = storage[storage.read_index.eq(0)][keys+['mean']].rename(columns={'mean':'initial_mean'})
    result = storage.merge(initial, on=keys, how='left', validate='many_to_one')
    result = result.merge(storage_reference_ranges(storage), on=keys[:4],
                          how='left', validate='many_to_one')
    result['difference'] = result['mean']-result.initial_mean
    result['normalized_difference'] = result.difference.div(
        result.initial_response_range.where(result.initial_response_range.gt(0)))
    result['normalization_reference'] = 'same_physical_plate_initial_observed_response_range'
    return result


def common_ph_handling(levels):
    """Compare every readout to UV/raw R on the same preparation level and run.

    Temperature columns use 23 C; both membrane columns use the overnight
    minus-membrane reference. The full reference R span supplies common units.
    Equal run numbers identify matched solution preparations, not extra cultures.
    Missing level references stay missing and no interpolation is performed.
    """
    f=levels[levels.assay.eq('pH') & levels.step_order.between(1,6) &
             levels.condition.eq('thymol_blue')].copy()
    f['reference_step_order']=np.where(f.step_order.le(4),1.,5.)
    reference=f[f.plate_type.eq('UV_transparent') & f.signal.eq('ratio_raw') &
                f.step_order.isin([1,5])].copy()
    keys=['run_number','condition','row_group','column','reference_step_order']
    reference['reference_step_order']=reference.step_order
    spans=reference.groupby(['run_number','reference_step_order']).response_mean.agg(
        lambda s:s.max()-s.min()).rename('reference_span')
    reference=reference.join(spans,on=['run_number','reference_step_order'])
    reference=reference[keys+['response_mean','reference_span','physical_plate_id','acquisition_id']].rename(
        columns={'response_mean':'reference_response','physical_plate_id':'reference_physical_plate_id',
                 'acquisition_id':'reference_acquisition_id'})
    result=f.merge(reference,on=keys,how='left',validate='many_to_one')
    result['difference']=result.response_mean-result.reference_response
    result['normalized_difference']=result.difference.div(result.reference_span.where(result.reference_span.gt(0)))
    result['normalization_reference']=np.where(result.reference_step_order.eq(1),
        'UV_raw_R_initial_23C_same_preparation_level_and_run',
        'UV_raw_R_overnight_minus_membrane_same_preparation_level_and_run')
    return result


def no_addition(value):return all(float(v)==0 for v in json.loads(value).values())


def lf_baseline_and_nearby_treatments(plates):
    """Select nearest elemental Ga/In dose, retaining all recorded no-addition plates."""
    raw=plates[plates.signal.eq('raw_A380')].copy()
    baseline=raw[raw.component_final_mM_json.map(no_addition)].drop_duplicates(
        ['batch_id','plate_id','column','ordinal']).copy()
    baseline['display_condition']='HH'
    selections=[baseline]
    for experiment,column,label,ion_mM in [
        ('LF_HH_indium_sulfate_OD0p001',5,'In',6.22),
        ('LF_HH_gallium_sulfate_OD0p001',7,'Ga',142.8)]:
        g=raw[raw.experiment_id.eq(experiment)&raw.column.eq(column)].copy()
        g['display_condition']=label;g['display_ion_mM']=ion_mM
        selections.append(g)
    return pd.concat(selections,ignore_index=True)
