"""Calculate optical fold changes, treatment repeatability and dose responses."""
from __future__ import annotations
import json
import numpy as np
import pandas as pd
from .measurement_summaries import storage_drift, common_ph_handling
from .peak_response import compare_peaks


def fold_values(numerator, denominator):
    """Conventional positive-signal folds; never offset or take absolute values."""
    a, b = np.asarray(numerator, float), np.asarray(denominator, float)
    good = np.isfinite(a) & np.isfinite(b) & (a > 0) & (b > 0)
    fold = np.full(a.shape, np.nan)
    np.divide(a, b, out=fold, where=good)
    status = np.where(~np.isfinite(a) | ~np.isfinite(b), 'missing_signal_or_reference',
                      np.where(good, 'positive_signal_and_reference', 'nonpositive_signal_or_reference'))
    return fold, np.log2(fold), status


def optical_folds(tables):
    storage = storage_drift(tables['storage_replicates'])
    storage['fold'], storage['log2_fold'], storage['fold_status'] = fold_values(storage['mean'], storage.initial_mean)
    storage['normalization_reference'] = 'same_physical_plate_readout_level_initial_signal'
    ph = common_ph_handling(tables['pub_robustness_levels'])
    ph['fold'], ph['log2_fold'], ph['fold_status'] = fold_values(ph.response_mean, ph.reference_response)
    levels = tables['pub_robustness_levels']
    iron = levels[levels.assay.eq('iron') & levels.step_order.between(1, 6) &
                  levels.signal.isin(['A380', 'A380_corrected'])].copy()
    iron['eligible_response'] = iron.response_mean.where(iron.complete_level)
    keys = ['plate_type', 'condition', 'acid_M', 'signal', 'step_order', 'column']
    iron = iron.groupby(keys, as_index=False).agg(Fe_mM=('predictor', 'mean'),
        mean_response=('eligible_response', 'mean'), sd_response=('eligible_response', 'std'),
        plate_n=('eligible_response', 'count'))
    reference = iron[iron.plate_type.eq('UV_transparent') & iron.condition.eq('no_thymol_blue') &
                     iron.acid_M.eq(1) & iron.signal.eq('A380') & iron.step_order.eq(1)]
    reference = reference.dropna(subset=['mean_response']).sort_values('Fe_mM')
    if len(reference) < 2:
        raise ValueError('Common UV/no-dye/1 M acid iron reference is unavailable')
    iron['reference_response'] = np.interp(iron.Fe_mM, reference.Fe_mM, reference.mean_response,
                                          left=np.nan, right=np.nan)
    iron['fold'], iron['log2_fold'], iron['fold_status'] = fold_values(iron.mean_response, iron.reference_response)
    iron['normalization_reference'] = 'initial_23C_UV_raw_A380_no_dye_1M_H2SO4_all_steps_including_membrane'
    iron['reference_matching'] = 'plate_level_means; adjacent_Fe_interpolation_within_support_only'
    return dict(storage_folds=storage, ph_handling_folds=ph, iron_handling_folds=iron)


def wavelength_endpoint_mask(pairs):
    """Select each physical plate's final acquisition before checking eligibility.

    Chronological indices belong to individual plates. Different acquisition
    counts do not justify dropping a plate's actual final measurement, and an
    invalid final wavelength pair must never be replaced by an earlier pair.
    """
    keys = ['species', 'experiment_id', 'citrate_mM', 'Ga_mM', 'batch_id', 'plate_id']
    final = pairs.groupby(keys, dropna=False).ordinal.transform('max')
    return pairs.ordinal.eq(final)


def timing_contrast(treated, references):
    """Observed-maximum delay, independently eligible from integrated response."""
    result = dict(delay_d=np.nan, control_time_d=np.nan,
        timing_status='missing_no_additive_control', control_timing_n=len(references))
    if not references:
        return result
    if any(r is None for r in references):
        result['timing_status'] = 'missing_assigned_control'
        return result
    times = np.asarray([r['time_to_max_d'] for r in references], float)
    if not np.isfinite(times).all() or not all(r['positive_rise'] for r in references):
        result['timing_status'] = 'unavailable_or_nonrising_reference_time'
        return result
    result['control_time_d'] = float(times.mean())
    if not treated['positive_rise'] or not np.isfinite(treated['time_to_max_d']):
        result['timing_status'] = 'unavailable_or_nonrising_treated_time'
        return result
    result.update(delay_d=float(treated['time_to_max_d'] - times.mean()),
                  timing_status='estimated_independently_of_integral')
    return result


def _joined(values):
    return ';'.join(sorted({str(x) for x in values if pd.notna(x)}))


def _effect_batches(units, identity):
    return units.groupby(identity + ['batch_id'], as_index=False, dropna=False).agg(
        response_percent=('response_percent', 'mean'), delay_d=('delay_d', 'mean'),
        absolute_time_d=('absolute_time_d', 'mean'), area=('area', 'mean'),
        response_n=('response_percent', 'count'), delay_n=('delay_d', 'count'),
        unit_n=('unit_id', 'nunique'), unit_ids=('unit_id', _joined),
        reference_ids=('reference_ids', _joined),
        final_read_maximum_n=('final_read_maximum', 'sum'),
        control_final_read_maximum_n=('control_final_read_maximum_n', 'sum'))


def _effect_summary(batches, identity):
    return batches.groupby(identity, as_index=False, dropna=False).agg(
        response_mean_percent=('response_percent', 'mean'), response_sd_pp=('response_percent', 'std'),
        response_batch_n=('response_percent', 'count'), delay_mean_d=('delay_d', 'mean'),
        delay_sd_d=('delay_d', 'std'), delay_batch_n=('delay_d', 'count'),
        absolute_time_mean_d=('absolute_time_d', 'mean'), absolute_time_sd_d=('absolute_time_d', 'std'),
        area_mean=('area', 'mean'), area_sd=('area', 'std'),
        batch_n=('batch_id', 'nunique'), batch_ids=('batch_id', _joined), unit_n=('unit_n', 'sum'),
        unit_ids=('unit_ids', _joined), final_read_maximum_n=('final_read_maximum_n', 'sum'),
        control_final_read_maximum_n=('control_final_read_maximum_n', 'sum'))


def effect_repeatability(tables):
    """Repeatability of treatment effects, with absolute controls retained separately."""
    from .plate_comparisons import CONDITION_KEYS
    metrics = tables['validation_plate_peak_metrics']
    keycols = CONDITION_KEYS + ['batch_id', 'plate_id']
    lookup = {tuple(r[c] for c in keycols): r for r in metrics.to_dict('records')}
    rows = []
    for source_row, r in enumerate(tables['validation_inhibition_plates'].to_dict('records'), 1):
        key = tuple(r[c] for c in keycols)
        treated = lookup[key]
        columns = [int(float(c)) for c in str(r['control_columns']).split(';') if c and c != 'nan']
        reference_keys = [(*key[:3], column, *key[4:]) for column in columns]
        references = [lookup.get(k) for k in reference_keys]
        timing = timing_contrast(treated, references)
        control = r['principal_final_mM'] == 0
        rows.append({**{c:r[c] for c in keycols}, 'unit_id':r['plate_id'],
            'additive_label':r['additive_label'], 'principal_final_mM':r['principal_final_mM'],
            'no_addition':control, 'response_percent':np.nan if control else r['activity_percent'],
            'delay_d':np.nan if control else timing['delay_d'],
            'absolute_time_d':treated['time_to_max_d'] if treated['positive_rise'] else np.nan,
            'area':treated['auc_to_max'], 'final_read_maximum':treated['maximum_at_final_read'],
            'control_final_read_maximum_n':sum(x['maximum_at_final_read'] for x in references if x is not None),
            'response_status':'self_normalized_reference_not_precision' if control else r['status'],
            'timing_status':'self_subtracted_reference_not_precision' if control else timing['timing_status'],
            'control_time_d':timing['control_time_d'],
            'reference_ids':';'.join('/'.join(str(x) for x in k) for k in reference_keys),
            'source_table':'validation_inhibition_plates', 'source_row':source_row})
    plate_units = pd.DataFrame(rows)
    identity = CONDITION_KEYS + ['additive_label', 'principal_final_mM', 'no_addition']
    plate_batches = _effect_batches(plate_units, identity)
    plate_summary = _effect_summary(plate_batches[~plate_batches.no_addition], identity)
    distributions = []
    for metric, field, n, unit in [('integrated_response','response_sd_pp','response_batch_n','percentage_points'),
                                  ('maximum_time','delay_sd_d','delay_batch_n','days')]:
        g = plate_summary.copy()
        g['metric'],g['between_batch_sd'],g['eligible_batch_n'],g['unit'] = metric,g[field],g[n],unit
        g['eligible'] = np.isfinite(g[field]) & g[n].ge(2)
        g['eligibility_reason'] = np.where(g.eligible,'finite_SD_at_least_two_usable_batch_results',
                                         'fewer_than_two_usable_batch_results')
        g['estimand'] = 'control_relative_integrated_response' if metric=='integrated_response' else 'control_relative_time_to_maximum_delay'
        distributions.append(g)
    fmetrics = tables['validation_flask_peak_metrics'].merge(
        tables['validation_flask_points'][['flask_ID','kind','no_addition']].drop_duplicates('flask_ID'),
        on='flask_ID',validate='one_to_one')
    fmetrics = fmetrics[fmetrics.kind.eq('biotic')]
    rows=[]
    for batch,g in fmetrics.groupby('batch',sort=True):
        refs=g[g.no_addition].to_dict('records')
        for r in g.to_dict('records'):
            contrast=compare_peaks(r,refs);timing=timing_contrast(r,refs)
            rows.append(dict(unit_id=r['flask_ID'],batch_id=batch,formulation=r['quantities_json'],
                no_addition=r['no_addition'],response_percent=np.nan if r['no_addition'] else contrast['activity_percent'],
                delay_d=np.nan if r['no_addition'] else timing['delay_d'],
                absolute_time_d=r['time_to_max_d'] if r['positive_rise'] else np.nan,
                area=r['auc_to_max'],final_read_maximum=r['maximum_at_final_read'],
                control_final_read_maximum_n=sum(x['maximum_at_final_read'] for x in refs),
                reference_ids=';'.join(x['flask_ID'] for x in refs),control_time_d=timing['control_time_d'],
                response_status='self_normalized_reference_not_precision' if r['no_addition'] else contrast['status'],
                timing_status='self_subtracted_reference_not_precision' if r['no_addition'] else timing['timing_status']))
    flask_units=pd.DataFrame(rows)
    flask_batches=_effect_batches(flask_units,['formulation','no_addition'])
    flask_summary=_effect_summary(flask_batches,['formulation','no_addition'])
    flask_summary=flask_summary[flask_summary.batch_n.ge(2)].copy()
    flask_summary['label']=flask_summary.formulation.map(lambda s:'0 mM Ga' if not json.loads(s) else '; '.join(f'{v:g} mM {k}' for k,v in json.loads(s).items()))
    flask_summary['response_status']=np.where(flask_summary.no_addition,'self_normalized_reference_not_precision','observed_batch_SD')
    flask_summary['timing_status']=np.where(flask_summary.no_addition,'self_subtracted_reference_not_precision','observed_batch_SD')
    combined=[]
    for metric,field,count,mean,unit in [
        ('integrated_response','response_sd_pp','response_batch_n','response_mean_percent','percentage_points'),
        ('maximum_time','delay_sd_d','delay_batch_n','delay_mean_d','days')]:
        for fmt,summary in [('plate',plate_summary),('flask',flask_summary)]:
            for r in summary.to_dict('records'):
                label=r['additive_label']+'; '+f"{r['principal_final_mM']:g} mM compound" if fmt=='plate' else r['label']
                reference=bool(r['no_addition'])
                combined.append(dict(format=fmt,condition=label,metric=metric,mean=r[mean],
                    between_batch_sd=r[field],eligible_batch_n=r[count],unit=unit,
                    candidate_batch_n=r['batch_n'],batch_ids=r['batch_ids'],unit_ids=r['unit_ids'],
                    status='self_reference_not_precision' if reference else
                    'estimated' if np.isfinite(r[field]) and r[count]>=2 else 'fewer_than_two_usable_batch_results',
                    source_formulation=r.get('formulation',r.get('group_id','')),
                    source_column=r.get('column',np.nan),source_experiment=r.get('experiment_id',''),
                    final_read_maximum_n=r['final_read_maximum_n'],
                    control_final_read_maximum_n=r['control_final_read_maximum_n']))
    return dict(effect_repeatability=pd.DataFrame(combined), effect_plate_units=plate_units, effect_plate_batch_means=plate_batches,
        effect_repeatability_distribution=pd.concat(distributions,ignore_index=True),
        flask_effect_units=flask_units,flask_effect_batch_means=flask_batches,
        flask_effect_repeatability=flask_summary)


def od_plate_trajectories(tables, *, ga_only=True):
    """Recorded same-well A600-A750; no substitution from other channels."""
    from .plate_comparisons import CONDITION_KEYS, compare_biotic_with_controls
    from .peak_response import plate_means
    source=tables['channel_observations']
    # Only the Ga-only formulation series supplies this dose assessment.
    formulations=tables['cross_format_plate_formulations'].copy()
    quantities=formulations.comparable_quantities_mM_json.map(json.loads)
    if ga_only:
        selected=quantities.map(lambda q: set(q).issubset({'Ga'})) & formulations.experiment_id.eq('AF_AFM1_gallium_sulfate_OD0p001')
        formulations=formulations[selected].copy()
    formulations['Ga_mM']=formulations.comparable_quantities_mM_json.map(lambda s:json.loads(s).get('Ga',0.))
    f=source.merge(formulations[CONDITION_KEYS+['Ga_mM']],on=CONDITION_KEYS,validate='many_to_one')
    identities=CONDITION_KEYS+['batch_id','plate_id','culture_id','row','well','source','ordinal','days','kind','physical_observation_key']
    a=pd.to_numeric(f.A600,errors='coerce').where(f.A600_eligible.eq(True))
    b=pd.to_numeric(f.A750,errors='coerce').where(f.A750_eligible.eq(True))
    frames=[]
    for signal,value in [('raw_A600',a),('A600_minus_A750',a-b),('raw_A750',b)]:
        q=f[identities].copy();q['signal']=signal;q['value']=value
        q['signal_status']=np.where(np.isfinite(value),'eligible','missing_or_ineligible_channel')
        wells,_=compare_biotic_with_controls(q)
        trajectories=plate_means(wells)
        frames.append(trajectories)
    out=pd.concat(frames,ignore_index=True).merge(formulations[CONDITION_KEYS+['Ga_mM']],on=CONDITION_KEYS,validate='many_to_one')
    out['source_table']='channel_observations'
    out['scatter_correction']='same_well_A600_minus_A750_before_initial_subtraction'
    return out


def od_dose_tables(tables):
    """Unmatched-dose-free A600/corrected-A600/flask-OD treatment effects."""
    from .plate_comparisons import CONDITION_KEYS
    from .peak_response import observed_peak
    curves=od_plate_trajectories(tables)
    keycols=CONDITION_KEYS+['batch_id','plate_id','signal']
    grouped=dict(tuple(curves[curves.signal.ne('raw_A750')].groupby(keycols)))
    peaks={key:observed_peak(g) for key,g in grouped.items()}
    rows=[]
    for key,g in grouped.items():
        first=g.iloc[0]; control=first.Ga_mM==0
        reference_keys=[k for k,v in grouped.items() if k[:3]==key[:3] and k[4:]==key[4:] and v.Ga_mM.iloc[0]==0]
        refs=[peaks[k] for k in reference_keys]
        contrast=compare_peaks(peaks[key],refs);timing=timing_contrast(peaks[key],refs)
        response=contrast['activity_percent'];delay=timing['delay_d']
        rows.append({**dict(zip(keycols,key)), 'unit_id':key[-2], 'format':'plate',
            'Ga_mM':first.Ga_mM,'is_reference':control,'response_percent':response,'delay_d':delay,
            'response_status':contrast['status'],'timing_status':timing['timing_status'],
            'absolute_time_d':peaks[key]['time_to_max_d'] if peaks[key]['positive_rise'] else np.nan,
            'area':peaks[key]['auc_to_max'],'control_area':contrast['control_integral'],
            'control_time_d':timing['control_time_d'],
            'final_read_maximum':peaks[key]['maximum_at_final_read'],
            'control_final_read_maximum_n':sum(x['maximum_at_final_read'] for x in refs),
            'reference_ids':';'.join('/'.join(str(x) for x in k) for k in reference_keys),
            'source_ids':_joined(g.sources), 'source_table':'od_plate_trajectories',
            'reference_interpretation':'normalization_reference_not_precision' if control else 'same_plate_batch_and_readout_no_Ga'})
    f=tables['validation_flask_points'].copy()
    q=f.quantities_json.map(lambda s:json.loads(s) if isinstance(s,str) else {})
    f=f[f.kind.eq('biotic') & q.map(lambda q:set(q).issubset({'Ga'}))].copy()
    value='calculated_planktonic_OD600_equivalent'
    # Keep the actual first measurement, including NA; never use first nonmissing.
    first=f[f.days_mean.eq(0)].drop_duplicates('flask_ID').set_index('flask_ID')[value]
    f['biotic_change']=f[value]-f.flask_ID.map(first)
    grouped=dict(tuple(f.groupby('flask_ID')))
    peaks={key:observed_peak(g) for key,g in grouped.items()}
    for flask,g in grouped.items():
        r=g.iloc[0]; control=bool(r.no_addition)
        reference_ids=[k for k,v in grouped.items() if v.batch.iloc[0]==r.batch and v.no_addition.iloc[0]]
        refs=[peaks[k] for k in reference_ids]
        contrast=compare_peaks(peaks[flask],refs);timing=timing_contrast(peaks[flask],refs)
        rows.append(dict(unit_id=flask,format='flask',signal='flask_OD600_apparent',batch_id=r.batch,
            condition_ID=r.condition_ID,Ga_mM=r.nominal_Ga_mM,is_reference=control,
            response_percent=contrast['activity_percent'],delay_d=timing['delay_d'],
            response_status=contrast['status'],timing_status=timing['timing_status'],
            absolute_time_d=peaks[flask]['time_to_max_d'] if peaks[flask]['positive_rise'] else np.nan,
            area=peaks[flask]['auc_to_max'],control_area=contrast['control_integral'],
            control_time_d=timing['control_time_d'],final_read_maximum=peaks[flask]['maximum_at_final_read'],
            control_final_read_maximum_n=sum(x['maximum_at_final_read'] for x in refs),
            reference_ids=';'.join(reference_ids),source_ids=flask,source_table='validation_flask_points',
            reference_interpretation='normalization_reference_not_precision' if control else 'same_batch_no_additive_flasks'))
    units=pd.DataFrame(rows)
    batches=_effect_batches(units,['format','signal','Ga_mM','is_reference'])
    summary=_effect_summary(batches,['format','signal','Ga_mM','is_reference'])
    point_summaries=[]
    for key,g in units.groupby(['format','signal','Ga_mM','is_reference'],dropna=False):
        row=dict(zip(['format','signal','Ga_mM','is_reference'],key))
        for label,value in [('response','response_percent'),('delay','delay_d')]:
            valid=g[g[value].notna()]
            row[label+'_unit_ids']=_joined(valid.unit_id)
            row[label+'_unit_n']=len(valid)
            row[label+'_unit_sd']=valid[value].std() if valid.batch_id.nunique()==1 and not key[-1] else np.nan
        point_summaries.append(row)
    summary=summary.merge(pd.DataFrame(point_summaries),on=['format','signal','Ga_mM','is_reference'],validate='one_to_one')
    summary['uncertainty_scope']='between_batch_SD_if_multiple_usable_batches; within_batch_physical_unit_SD_if_one; self_reference_SD_withheld'
    # Never plot self-normalized reference dispersion as independent precision.
    summary.loc[summary.is_reference,['response_sd_pp','delay_sd_d']]=np.nan
    for label,sd,count in [('response','response_sd_pp','response_batch_n'),('delay','delay_sd_d','delay_batch_n')]:
        summary[label+'_plotted_sd']=summary[sd].where(summary[count].gt(1),summary[label+'_unit_sd']).where(~summary.is_reference)
        summary[label+'_sd_scope']=np.where(summary.is_reference,'self_reference_SD_withheld',
            np.where(summary[count].gt(1),'between_starting_batch_SD','within_starting_batch_physical_unit_SD'))
    return dict(od_plate_trajectories=curves,od_dose_units=units,
                od_dose_batch_means=batches,od_dose_summary=summary)
