"""Si robustness validation."""
from __future__ import annotations
import numpy as np
import pandas as pd
from .publication_si_optical_results import SIGNAL, _summary, _ci
from .publication_storage_support import storage_common_support


def _spec(key,title,caption,frame,display,sources,panels,claims,results,support=None):
    assert len(frame) and all(c in frame for c in display),key
    labels={'assay':'Assay','acid_M':'H2SO4 (M)','dye':'Dye','readout':'Readout','material':'Plate',
            'contrast':'Candidate / reference','eligible_plate_counts':'Eligible n (span/SD/local)',
            'span_display':'Span fold [95% CI]','sd_display':'SD fold [95% CI]',
            'local_display':'Δ local separation [95% CI]','plate_n':'Plates',
            'last_recorded_day_mean':'Final day','mean_observed_drift_display':'Mean drift [range]',
            'max_absolute_plate_drift_display':'Maximum |plate drift| [range]',
            'max_absolute_level_drift_display':'Maximum |level drift| [range]',
            'rms_observed_level_drift_display':'RMS drift [range]',
            'same_direction_pair_percent':'Persistent sign (%)','cohort':'Cohort','kind':'Wells',
            'attempt_n':'Acquired well reads','raw_percent_mean':'R valid (%)',
            'corrected_percent_mean':'R* valid (%)','loss_percent_mean':'R-only loss (%)',
            'gain_percent_mean':'R*-only gain (%)','difference_pp_mean':'Yield Δ (pp)',
            'cancellation_median_mean':'Denominator retained fraction',
            'raw_comparison_yield_percent':'R net valid (%)',
            'corrected_comparison_yield_percent':'R* net valid (%)',
            'comparison_loss_percent':'Net yield loss (%)','opposite_sign_percent':'Opposite sign (%)',
}
    return dict(key=key,title=title,caption=caption,frame=frame,
                display_columns=display,source_tables=sources,figure_panels=panels,
                claims=claims,result_columns=results,column_labels={c:labels[c] for c in display if c in labels},
                support_tables=support or {})


def _true(values):
    return values.eq(True) | values.astype(str).str.lower().eq('true')


def _curve_metrics(pair,field):
    """Common complete levels; local contrasts never bridge an omitted level."""
    rows=[]
    for level,g in pair.groupby('level_id',sort=False):
        usable=_true(g.pair_eligible)&_true(g.same_predictor)&np.isfinite(g.predictor)
        complete=bool(usable.all() and len(g)>=2)
        values=g.loc[usable,field].to_numpy(float)
        rows.append(dict(level=level,predictor=g.predictor.iloc[0],complete=complete,
                         mean=values.mean() if complete else np.nan,
                         sd=values.std(ddof=1) if complete else np.nan,n=len(values)))
    levels=pd.DataFrame(rows).sort_values('predictor',kind='stable')
    eligible=levels[levels.complete & np.isfinite(levels.predictor)]
    span=float(eligible['mean'].max()-eligible['mean'].min()) if len(eligible)>=2 else np.nan
    denom=float((eligible.n-1).sum())
    sd=float(np.sqrt(((eligible.n-1)*eligible.sd**2).sum()/denom)) if denom>0 else np.nan
    direction=-1 if pair.assay.iloc[0]=='pH' else 1
    local=[]
    for j in range(1,len(levels)):
        a,b=levels.iloc[j-1],levels.iloc[j]
        if not a.complete or not b.complete or not np.isfinite(a.predictor+b.predictor) or b.predictor<=a.predictor:
            continue
        pooled=float(np.sqrt((a.sd**2+b.sd**2)/2))
        if pooled>0 and np.isfinite(pooled):local.append(direction*(b['mean']-a['mean'])/pooled)
    return dict(span=span,sd=sd,local=float(np.median(local)) if local else np.nan,
                level_n=len(eligible),adjacent_n=len(local),wrong_direction_n=int(np.sum(np.asarray(local)<0)))


def _prefixed_summary(values,prefix,fold=False):
    s=_summary(values,test=False)
    row={prefix+'_'+key:value for key,value in s.items()}
    if fold:
        for field in ['estimate','ci_low','ci_high']:
            row[prefix+'_fold_'+field]=float(2**s[field]) if np.isfinite(s[field]) else np.nan
        row[prefix+'_display']=_ci(row[prefix+'_fold_estimate'],row[prefix+'_fold_ci_low'],row[prefix+'_fold_ci_high'])
    else:row[prefix+'_display']=s['estimate_ci']
    return row


def _handling(tables):
    q=tables['pub_robustness_paired_wells'].copy()
    q=q[~q.comparison.eq('storage_minus_initial')]
    q=q[~q.comparison.eq('nonUV_minus_UV')|q.step_id.eq('initial_23_C')]
    # Keep the lid comparison used in Figure 3a (undyed ferric signal).
    q=q[~(q.comparison.eq('lid_present_minus_absent')&q.assay.eq('iron')&q.condition.ne('no_thymol_blue'))]
    unit_rows=[]
    for identity,g in q.groupby(['comparison_id','condition','signal'],dropna=False,sort=True):
        if g.well.duplicated().any():raise ValueError(f'Duplicate matched well: {identity}')
        ref=_curve_metrics(g,'response_reference');cand=_curve_metrics(g,'response_candidate')
        row=g.iloc[0][['assay','comparison','step_id','plate_type','acid_M','condition','signal',
                      'candidate_physical_plate_id','reference_physical_plate_id']].to_dict()
        for metric in ['span','sd']:
            r,c=ref[metric],cand[metric]
            row[metric+'_effect']=float(np.log2(c/r)) if np.isfinite(c+r) and c>0 and r>0 else np.nan
            row[metric+'_reference']=r;row[metric+'_candidate']=c
        row['local_effect']=cand['local']-ref['local']
        row.update(local_reference=ref['local'],local_candidate=cand['local'],
                   common_complete_level_n=ref['level_n'],common_adjacent_n=ref['adjacent_n'],
                   reference_wrong_direction_n=ref['wrong_direction_n'],candidate_wrong_direction_n=cand['wrong_direction_n'])
        unit_rows.append(row)
    units=pd.DataFrame(unit_rows);rows=[]
    groups=['assay','comparison','step_id','plate_type','acid_M','condition','signal']
    for identity,g in units.groupby(groups,dropna=False,sort=True):
        if g.candidate_physical_plate_id.duplicated().any():raise ValueError(f'Duplicate inferential plate: {identity}')
        row=dict(zip(groups,identity));row['paired_plate_n']=len(g)
        for metric in ['span','sd','local']:
            row.update(_prefixed_summary(g[metric+'_effect'],metric,fold=metric!='local'))
            row[metric+'_reference_mean']=g[metric+'_reference'].mean()
            row[metric+'_candidate_mean']=g[metric+'_candidate'].mean()
        row['eligible_plate_counts']='/'.join(str(row[metric+'_plate_n']) for metric in ['span','sd','local'])
        for metric in ['common_complete_level_n','common_adjacent_n']:
            row[metric+'_min']=int(g[metric].min());row[metric+'_max']=int(g[metric].max())
        row['reference_reversal_plate_n']=int(g.reference_wrong_direction_n.gt(0).sum())
        row['candidate_reversal_plate_n']=int(g.candidate_wrong_direction_n.gt(0).sum())
        row['plate_ids']=';'.join(g.candidate_physical_plate_id)
        row['readout']=SIGNAL[row['signal']]
        row['dye']='10 µM' if row['condition']=='thymol_blue' else '0 µM'
        row['material']={'UV_transparent':'UV','non_UV_transparent':'non-UV'}.get(row['plate_type'],'—')
        row['contrast']={'incubation_4_C':'4°C / initial','incubation_30_C':'30°C / initial',
                         'incubation_75_C':'75°C / initial','overnight_23_C':'overnight / initial',
                         'membrane_23_C':'membrane / overnight','initial_23_C':'non-UV / UV'}.get(row['step_id'],'lid+ / lid−')
        rows.append(row)
    frame=pd.DataFrame(rows)
    return _spec('handling_discrimination','Handling effects on span, scatter and local separation',
        'Each acquisition/reference pair is reduced to one estimate per physical-plate pair. Span is max−min '
        'of common complete preparation-level means; pooled technical SD uses within-level sample variances. '
        'The same eligible physical wells and complete levels enter both sides. Reported span and SD folds '
        'are geometric means of candidate/reference, with 95% Student-t intervals on log2 folds over plate '
        'pairs. Local separation is the median signed difference between neighboring level means divided by '
        'sqrt[(SD1²+SD2²)/2], positive in the expected direction (decreasing pH signal, increasing ferric signal). '
        'Omitted levels are never bridged. Its effect is candidate minus reference, with an arithmetic-mean '
        'Student-t interval. This empirical standardized separation is not a detection limit or proof of '
        'concentration recovery. Missing/zero reference metrics produce NA rather than a pseudocount. '
        'Temperature conditions use their own material’s initial 23°C acquisition; membrane uses the preceding '
        'overnight acquisition; material comparison is only at initial acquisition. Lid pairs are actual '
        'lid-present/lid-absent readings. Typically three paired plates summarize differences among '
        'plate estimates, including preparation variation where plates use independently prepared '
        'standard series. Repeated measurements within a plate do not add preparations. '
        'Time, thermal history and handling order are confounded; intervals are '
        'descriptive pointwise intervals, with no confirmatory significance or equivalence claim. Eligible '
        'plate counts are given in span/SD/local order; scheduled pair counts are retained in the CSV.',
        frame,['assay','acid_M','dye','readout','material','contrast','eligible_plate_counts','span_display','sd_display','local_display'],
        ['pub_robustness_paired_wells'],['2a','2e','2f','2g','2h','3a','3g','3h','3i','3j','3k','3l'],
        ['Handling can change optical span, technical scatter and adjacent-reference separation independently.'],
        ['span_fold_estimate','sd_fold_estimate','local_estimate'],{'physical_plate_metrics':units})


def _storage(tables):
    q=tables['storage_folds'].copy();unit_rows=[]
    keys=['assay','condition','signal','physical_plate_id']
    for identity,g in q.groupby(keys,dropna=False,sort=True):
        if g.duplicated(['row_group','column','read_index']).any():raise ValueError('Duplicate storage level acquisition')
        g=g[g.read_index.gt(0)].copy();g['drift']=100*g.normalized_difference
        acquisitions=g.groupby('read_index',as_index=False).agg(value=('drift','mean'),days=('days','mean'))
        finite=acquisitions[np.isfinite(acquisitions.value)]
        row=dict(zip(keys,identity));row.update(
            acquisition_n=len(acquisitions),eligible_acquisition_n=len(finite),
            last_recorded_day=float(acquisitions.days.max()),
            unavailable_level_acquisition_n=int((~np.isfinite(g.drift)).sum()),
            available_level_acquisition_n=int(np.isfinite(g.drift).sum()),
            mean_observed_drift=finite.value.mean(),
            max_absolute_plate_drift=finite.value.abs().max(),
            max_absolute_level_drift=g.drift.abs().max(),
            rms_observed_level_drift=float(np.sqrt(np.nanmean(g.drift**2))) if np.isfinite(g.drift).any() else np.nan)
        common=acquisitions.sort_values('read_index'); signs=np.sign(common.value.to_numpy())
        eligible=np.isfinite(signs[:-1])&np.isfinite(signs[1:])
        consecutive=(np.diff(common.read_index.to_numpy())==1)&eligible
        row['consecutive_pair_n']=int(consecutive.sum())
        row['same_direction_pair_n']=int(((signs[:-1]*signs[1:]>0)&consecutive).sum())
        unit_rows.append(row)
    units=pd.DataFrame(unit_rows);rows=[]
    for identity,g in units.groupby(['assay','condition','signal'],sort=True):
        row=dict(zip(['assay','condition','signal'],identity))
        row['readout']=SIGNAL[row['signal']];row['dye']='10 µM' if row['condition']=='thymol_blue' else '0 µM'
        row['plate_n']=len(g);row['last_recorded_day_mean']=g.last_recorded_day.mean()
        for field in ['mean_observed_drift','max_absolute_plate_drift','max_absolute_level_drift','rms_observed_level_drift']:
            row[field+'_mean']=g[field].mean();row[field+'_min']=g[field].min();row[field+'_max']=g[field].max()
            row[field+'_display']=f'{g[field].mean():.3g} [{g[field].min():.3g}, {g[field].max():.3g}]'
        for field in ['acquisition_n','eligible_acquisition_n','available_level_acquisition_n','unavailable_level_acquisition_n','consecutive_pair_n','same_direction_pair_n']:
            row[field]=int(g[field].sum())
        row['same_direction_pair_percent']=100*row['same_direction_pair_n']/row['consecutive_pair_n'] if row['consecutive_pair_n'] else np.nan
        rows.append(row)
    result = _spec('storage_trajectory_drift','Drift throughout the observed storage sequence',
        'At every recorded post-initial acquisition, drift is 100×(Xt−X0)/span(X0), using the same physical '
        'plate, readout and preparation level as its initial reference. Level drifts are averaged within '
        'acquisition before obtaining each plate’s mean observed drift and maximum absolute plate-average '
        'drift. Maximum absolute level drift and RMS of all observed level-by-acquisition drifts additionally '
        'retain effects that cancel in a level average. Brackets are the observed minimum and maximum over '
        'physical plates, not confidence intervals. The observed-acquisition mean weights acquisitions equally '
        'and makes no interpolation or stability claim between them. Consecutive-pair persistence counts '
        'adjacent recorded indices with the same nonzero sign of plate-average drift; missing indices or '
        'undefined values are not bridged. Maxima are selected descriptive statistics without P tests. '
        'All drift measures are percent of initial span. Missing corrected values remain unavailable; the '
        'analysis does not establish an acceptable storage duration without a prespecified tolerance.',
        pd.DataFrame(rows),['assay','dye','readout','plate_n','last_recorded_day_mean','mean_observed_drift_display',
                           'max_absolute_plate_drift_display','max_absolute_level_drift_display','rms_observed_level_drift_display',
                           'same_direction_pair_percent'],
        ['storage_folds'],['2b','2c','3b','3c','3d','3e'],
        ['Endpoint drift alone can understate larger excursions during the recorded storage interval.'],
        ['mean_observed_drift_mean','max_absolute_plate_drift_mean','rms_observed_level_drift_mean'],
        {'physical_plate_metrics':units})
    support = storage_common_support(tables)
    result['support_tables'].update({'common_support_' + name: frame
                                     for name, frame in support.items()})
    result['caption'] += (' Additional support exports compare available-case drift with a fixed set of '
                          'reference levels valid for both readouts at every acquired index, including '
                          'the initial read. Selection retains each original plate/readout initial span. '
                          'They report level exclusions, acquisition means, actual final indices/days and '
                          'acquisition gaps. Fixed-support signed endpoint means use 95% Student-t '
                          'intervals across physical plates; maximum absolute plate and level drifts '
                          'remain descriptive means with observed plate ranges. Empty support is '
                          'unavailable, and missing indices are never interpolated or bridged.')
    return result


def _experiment_labels(channels):
    # Printed labels summarize an analytical cohort; machine outputs retain IDs.
    labels={}
    for experiment,g in channels.groupby('experiment_id',sort=True):
        text=experiment
        for a,b in [('Acidithiobacillus_',''),('potassium_tetrathionate','tetrathionate'),
                    ('sodium_thiosulfate','thiosulfate'),('aluminum_sulfate','Al'),('cobalt_ii_sulfate','Co'),
                    ('copper_ii_sulfate','Cu'),('gallium_sulfate','Ga'),('potassium_sulfate','K'),
                    ('zinc_sulfate','Zn'),('sulfuric_acid','acid'),('threefold_dye','3× dye'),('and_','+'),
                    ('OD0p001','OD 10⁻³'),('OD0p01','OD 10⁻²'),('OD0p1','OD 10⁻¹')]:
            text=text.replace(a,b)
        labels[experiment]=text.replace('_',' ')
    return labels


def _sulfur_direction(tables):
    q=tables['validation_sulfur_trajectory_evidence'];ids=['physical_observation_key']
    raw=q[q.signal.eq('raw_R545_485')].copy();corr=q[q.signal.eq('corrected_R545_485')].copy()
    if raw.duplicated(ids).any() or corr.duplicated(ids).any():raise ValueError('Duplicate sulfur ratio trajectory observation')
    pair=raw.merge(corr[ids+['net_change','comparison_status']],on=ids,validate='one_to_one',suffixes=('_raw','_corrected'))
    if len(pair)!=len(raw) or len(pair)!=len(corr):raise ValueError('Missing paired sulfur acquisition')
    pair['raw_ok']=pair.comparison_status_raw.eq('eligible')&np.isfinite(pair.net_change_raw)
    pair['corrected_ok']=pair.comparison_status_corrected.eq('eligible')&np.isfinite(pair.net_change_corrected)
    labels=_experiment_labels(tables['context_sulfur_channels']);plate_rows=[]
    # Mean net change per condition/physical plate/acquisition; require all
    # recorded biotic replicates and complete controls, and pair on common wells.
    condition_rows=[]
    group=['experiment_id','batch_id','plate_id','group_id','column','ordinal','days']
    for identity,g in pair.groupby(group,sort=True,dropna=False):
        if len(g)!=3 or g.well.nunique()!=3:
            raise ValueError(f'Incomplete scheduled sulfur biotic triplicate: {identity}')
        valid=g.raw_ok&g.corrected_ok
        row=dict(zip(group,identity));row.update(raw_ok=bool(g.raw_ok.all()),corrected_ok=bool(g.corrected_ok.all()),
            common_ok=bool(valid.all()),raw_net=float(g.net_change_raw.mean()) if valid.all() else np.nan,
            corrected_net=float(g.net_change_corrected.mean()) if valid.all() else np.nan,
            scheduled_biotic_well_n=len(g))
        condition_rows.append(row)
    conditions=pd.DataFrame(condition_rows)
    for identity,g in conditions.groupby(['experiment_id','batch_id','plate_id'],sort=True):
        after=g[g.days.gt(0)];common=after[after.common_ok]
        nonzero=common.raw_net.ne(0)&common.corrected_net.ne(0)
        eligible_sign=common[nonzero]
        opposite=eligible_sign.raw_net*eligible_sign.corrected_net<0
        row=dict(zip(['experiment_id','batch_id','plate_id'],identity))
        row.update(postinitial_condition_acquisition_n=len(after),common_valid_n=len(common),
                   opposite_sign_n=int(opposite.sum()),nonzero_common_n=len(eligible_sign),
                   opposite_sign_percent=float(100*opposite.mean()) if len(opposite) else np.nan,
                   raw_comparison_yield_percent=100*after.raw_ok.mean(),corrected_comparison_yield_percent=100*after.corrected_ok.mean(),
                   comparison_loss_percent=100*(after.raw_ok&~after.corrected_ok).mean())
        plate_rows.append(row)
    plates=pd.DataFrame(plate_rows);rows=[]
    for experiment,g in plates.groupby('experiment_id',sort=True):
        row=dict(experiment_id=experiment,cohort=labels[experiment],plate_n=len(g))
        for field in ['postinitial_condition_acquisition_n','common_valid_n','opposite_sign_n','nonzero_common_n']:
            row[field]=int(g[field].sum())
        for field in ['opposite_sign_percent','raw_comparison_yield_percent','corrected_comparison_yield_percent',
                      'comparison_loss_percent']:
            row[field]=g[field].mean();row[field+'_min']=g[field].min();row[field+'_max']=g[field].max()
        rows.append(row)
    return _spec('sulfur_correction_direction','Correction effects on sulfur-signal direction',
        'Net change is the change from the initial biotic ratio minus the corresponding change in the '
        'condition-matched abiotic ratio. A condition/plate/acquisition is valid only when all recorded '
        'biotic replicates have valid signal, baseline and complete matched abiotic controls. Raw/corrected '
        'direction comparisons use the identical valid wells; condition-level means are computed within '
        'physical plate before signs are compared. The initial zero-change acquisition is excluded. '
        'Opposite direction is a negative product of the raw and corrected nonzero mean net changes. '
        'Percentages are calculated within physical plate, then averaged equally over '
        'plates. Reused controls, concentrations and repeated acquisitions are not independent trials; '
        'no binomial inference is made. A sign near zero may change through measurement noise, and this '
        'threshold-free direction diagnostic is not a biological activity or acidification criterion. '
        'Missing controls or invalid denominators reduce evaluable yield; they are never negative results.',
        pd.DataFrame(rows),['cohort','plate_n','raw_comparison_yield_percent','corrected_comparison_yield_percent',
                           'comparison_loss_percent','opposite_sign_percent'],
        ['validation_sulfur_trajectory_evidence','context_sulfur_channels'],['8a','8b','8c','8d','8e','8f'],
        ['A750 correction can change the direction of matched sulfur-associated net signals.'],
        ['comparison_loss_percent','opposite_sign_percent'],
        {'physical_plate_metrics':plates,'paired_condition_acquisitions':conditions})
