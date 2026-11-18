"""Si optical results."""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy.stats import t as student_t
from .publication_calibration_fourpl import fit_fourpl, fit_ferric_response
from .publication_ph_support import augment_ph_transfer
from .publication_ratio_support import build_media_analysis_support


SIGNAL = {'ratio_raw':'R', 'ratio_corrected':'R*', 'raw_ratio':'R',
          'A750_corrected_ratio':'R*','a380':'A380','a380_corrected':'A*380',
          'A380':'A380','A380_corrected':'A*380','a428':'A428','a428_corrected':'A*428'}


MEDIUM = {'Starkey_brown':'St-BS','Starkey_colloidal':'St-YS'}


TECHNICAL = ('Intervals are two-sided 95% Student-t intervals over physical-plate '
             'estimates (usually n=3), conditional on the recorded preparations and references; '
             'they include between-preparation variation where each plate uses an independent '
             'standard series and otherwise describe technical plate variation. They do not '
             'quantify reference-measurement uncertainty or biological variation. ')


def _ci(value, low, high, digits=3):
    if not np.isfinite(value): return 'NA'
    if not np.isfinite(low+high): return f'{value:.{digits}g} (NA)'
    return f'{value:.{digits}g} [{low:.{digits}g}, {high:.{digits}g}]'


def _summary(values, test=True):
    a=np.asarray(values,float); a=a[np.isfinite(a)]; n=len(a)
    m=float(a.mean()) if n else np.nan
    sd=float(a.std(ddof=1)) if n>1 else np.nan
    lo=hi=p=np.nan
    if n>1 and np.isfinite(sd):
        margin=student_t.ppf(.975,n-1)*sd/np.sqrt(n)
        lo,hi=m-margin,m+margin
        if test and sd>0: p=float(2*student_t.sf(abs(m)/(sd/np.sqrt(n)),n-1))
    return {'plate_n':n,'estimate':m,'sd':sd,'ci_low':lo,'ci_high':hi,
            'p_value':p,'estimate_ci':_ci(m,lo,hi),
            'status':'estimated' if n>1 else ('single_plate_no_interval' if n else 'unavailable')}


def _holm(frame, col='p_value'):
    frame=frame.copy(); values=pd.to_numeric(frame[col],errors='coerce').to_numpy()
    eligible=np.flatnonzero(np.isfinite(values)); adjusted=np.full(len(values),np.nan)
    if len(eligible):
        order=eligible[np.argsort(values[eligible])]
        adjusted[order]=np.minimum(1,np.maximum.accumulate(values[order]*(len(order)-np.arange(len(order)))))
    frame['p_holm']=adjusted
    frame['multiplicity_family_n']=len(eligible)
    return frame


def _table(key,title,caption,frame,columns,sources,panels,claims):
    frame=frame.reset_index(drop=True)
    return dict(key=key,title=title,caption=caption,frame=frame,
                display_columns=columns,source_tables=sources,figure_panels=panels,claims=claims,
                result_columns=[c for c in frame if c in {
                    'estimate','p_holm','fit_rmse','r_squared','bottom','span','midpoint','slope',
                    'bias_pH','rmse_pH','technical_sd_pH','rmse_difference_pH','bias_difference_pH',
                    'technical_sd_ratio','sse_reduction_percent','acid_specific_rmse','intercept_M',
                    'residual_sd_mM','span_mean','technical_sd_mean','span_sd_mean',
                    'median_rank_rho','plates_with_adjacent_reversal'}])


def _cluster_parameters(group,x,y,iron=False):
    q=group[np.isfinite(group[x])&np.isfinite(group[y])].copy()
    if iron: q=q[q[x].ge(0)]
    fit=fit_ferric_response(q[x],q[y]) if iron else fit_fourpl(q[x],q[y])
    ids=q.physical_plate_id.to_numpy(); units=np.unique(ids); g=len(units)
    lows=highs=np.full(4,np.nan)
    if g>1 and fit.converged and np.isfinite(fit.influence).all():
        scores=np.array([fit.influence[ids==u].sum(axis=0) for u in units])
        covariance=scores.T@scores*g/(g-1)
        delta=student_t.ppf(.975,g-1)*np.sqrt(np.maximum(0,np.diag(covariance)))
        lows,highs=fit.estimates-delta,fit.estimates+delta
    row=dict(plate_n=g,plate_mean_n=len(q),level_n=q[x].nunique(),
             fit_rmse=float(np.sqrt(fit.sse/len(q))) if len(q) else np.nan,
             r_squared=fit.r_squared,identified=fit.identified,converged=fit.converged,
             identification_flags='; '.join(fit.flags),
             ci_method='physical_plate_clustered_HC3_delta_t',ci_df=g-1)
    for i,name in enumerate(('bottom','span','midpoint','slope')):
        row[name]=fit.estimates[i]; row[name+'_ci_low']=lows[i];row[name+'_ci_high']=highs[i]
        row[name+'_ci']=_ci(fit.estimates[i],lows[i],highs[i])
    return row


def _models(t):
    result=[]; rows=[]
    q=t['ph_replicates']
    for (dose,signal),g in q.groupby(['nominal_dye_uM','signal'],sort=True):
        row=dict(dye_uM=dose,readout=SIGNAL[signal],**_cluster_parameters(g,'mean_measured_pH','mean'))
        row['unavailable_plate_mean_n']=int((~np.isfinite(g['mean'])).sum())
        rows.append(row)
    result.append(_table('ph_response_fit','Fitted thymol-blue response',
        'Figure 1b–c uses Y=b+s/[1+10^((pH−m)h)], where b is the lower asymptote, s the span, '
        'm the midpoint and h the slope. Fits use four-well physical-plate means against the mean '
        'recorded pH for each preparation level. Parameter intervals use the same plate-clustered '
        'HC3 linearized influence, G/(G−1) correction and t(G−1) convention as the plotted response bands. '
        'RMSE is in ratio units. The 1 µM corrected curve remains unestimable; denominator failures are '
        'not zero response. No inverse-pH prediction accuracy follows from these fit residuals. '
        'Parameter identification is the existing Jacobian/coverage gate, not a new test of accuracy. '
        'Each dye-dose cohort uses three independently prepared standard series on three UV-transparent '
        'physical plates. The 10 µM cohort is the initial configuration acquisition; the 1 µM cohort uses '
        'separate initial-calibration plates. Technical wells within a preparation level share its recorded '
        'reference pH. These independent preparations do not provide independent analytical validation '
        'of the reference values. The non-UV acid-final plates are a separate cohort.',
        pd.DataFrame(rows),['dye_uM','readout','plate_n','plate_mean_n','bottom_ci','span_ci','midpoint_ci','slope_ci','fit_rmse','identified'],
        ['ph_replicates'],['1b','1c'],['Dose and correction change fitted ratio response.','The 1 µM corrected curve is not supported.']))
    q=t['calibration_level_means']; q=q[q.assay.eq('iron')&q.condition.eq('no_thymol_blue')&q.signal.eq('a380')]
    # Use the same acid-configuration cohort for fitting and precision estimates.
    q=q[q.experiment_id.str.startswith('iron_configuration_') & q.experiment_id.str.endswith('_23_C_UV')]
    rows=[]
    for acid,g in q.groupby('acid_M'):
        row=dict(acid_M=acid,**_cluster_parameters(g,'axis_value','mean',True));rows.append(row)
    result.append(_table('ferric_response_fit','Acid-specific ferric absorbance response',
        'The displayed model is A380=b+s[1−exp{−ln(2)(C/K)^h}], with C the recorded total-Fe reference '
        'in mM, including zero; b and s are AU, K is mM and h is dimensionless. Four free parameters '
        'are estimated separately for each sulfuric-acid background without dye. Intervals use the '
        'physical-plate-clustered HC3/t convention described for the pH response. RMSE describes fitted '
        'optical response, not concentration recovery. Reference concentrations are treated as fixed; '
        'reference measurement uncertainty is not propagated.',pd.DataFrame(rows),
        ['acid_M','plate_n','plate_mean_n','bottom_ci','span_ci','midpoint_ci','slope_ci','fit_rmse','identified'],
        ['calibration_level_means'],['1f'],['Acid-specific empirical curves describe a nonlinear ferric response.']))
    return result


def _initial_precision(t):
    q=t['calibration_correction_effects'];q=q[q.assay.eq('pH')].copy()
    rows=[]
    for metric,g in q.groupby('metric'):
        ratio=g.corrected_value/g.raw_value
        logs=np.log2(ratio.where((ratio>0)&np.isfinite(ratio)))
        row=dict(metric={'response_span':'response span','pooled_technical_sd':'technical SD','range_over_technical_sd':'span/technical SD'}[metric],
                 raw_mean=g.raw_value.mean(),corrected_mean=g.corrected_value.mean(),**_summary(logs))
        row['geometric_ratio']=2**row['estimate'];row['ratio_ci_low']=2**row['ci_low'];row['ratio_ci_high']=2**row['ci_high']
        row['corrected_over_raw_ci']=_ci(row['geometric_ratio'],row['ratio_ci_low'],row['ratio_ci_high']);rows.append(row)
    return _table('initial_correction_precision','Effect of A750 correction on initial ratio precision',
        'Paired results for the 10 µM thymol-blue cohort in Figure 1d. Each physical plate contributes '
        'one corrected/raw ratio of response span, pooled within-preparation-level four-well sample SD, '
        'or their quotient. Span is max−min of level means; pooled SD=sqrt[Σ(nj−1)sj²/Σ(nj−1)]. '
        'Intervals and two-sided zero tests are on log2(corrected/raw); reported folds back-transform '
        'the interval. Holm adjustment covers the three metrics. '+TECHNICAL+
        'An increase in span accompanied by increased technical SD is not automatically improved precision.',
        _holm(pd.DataFrame(rows)),['metric','plate_n','raw_mean','corrected_mean','corrected_over_raw_ci','p_holm'],
        ['calibration_correction_effects'],['1d'],['Correction changes optical span and well-to-well scatter differently.'])


def _transfer(t):
    q=t['ph_matrix_accuracy'].copy();q['readout']=q.signal.map(SIGNAL)
    q['matrix_display']=q.matrix.map({'acid_held_out_plate':'held-out acid','ferric_transfer':'ferric transfer'})
    q['analysis_status']=np.where(q.finite_prediction_n.gt(0),'estimated','unavailable')
    q['bias_ci']=[_ci(*r) for r in q[['bias_pH','bias_ci_low','bias_ci_high']].to_numpy()]
    q['rmse_ci']=[_ci(*r) for r in q[['rmse_pH','rmse_ci_low','rmse_ci_high']].to_numpy()]
    q['technical_sd_ci']=[_ci(*r) for r in q[['technical_sd_pH','technical_sd_ci_low','technical_sd_ci_high']].to_numpy()]
    a=_table('ph_transfer_accuracy','Conditional pH prediction and ferric-matrix transfer',
        'The acid calibration is trained on two physical plates and evaluated on the held-out acid plate '
        'or transferred to ferric mixtures. Only shared 0.25-pH bins are used; bins receive equal weight. '
        'Bias is predicted−reference pH, RMSE is root-mean-square error, and technical SD is within-reference '
        'scatter of inferred pH. The 199-resample whole-well bootstrap refits the model and conditions on '
        'the observed plate blocks and unquantified reference-pH error. Training uses independently '
        'prepared acid standard series measured in UV-transparent plates (acid_calibration). '
        'The transfer evaluation uses ferric mixtures (ferric_matrix). '
        'The recorded reference values are shared within preparation levels and are not independently '
        'validated by this prediction analysis. All attempted strata remain, including no-common-support '
        'rows; an inverse failure is unavailable, never zero error.',q,
        ['acid_m','matrix_display','readout','finite_prediction_n','inverse_failure_n','bias_ci','rmse_ci','technical_sd_ci','analysis_status'],
        ['ph_matrix_accuracy','ph_matrix_common_support'],['1b','1c','3h'],['An acid-only pH response need not transfer accurately into ferric matrices.'])
    q=t['ph_matrix_readout_comparisons'].copy()
    q['matrix_display']=q.matrix.map({'acid_held_out_plate':'held-out acid','ferric_transfer':'ferric transfer'})
    for out,v,l,h in [('rmse_change_ci','rmse_difference_pH','rmse_difference_ci_low','rmse_difference_ci_high'),
                      ('bias_change_ci','bias_difference_pH','bias_difference_ci_low','bias_difference_ci_high'),
                      ('sd_ratio_ci','technical_sd_ratio','technical_sd_ratio_ci_low','technical_sd_ratio_ci_high')]:
        q[out]=[_ci(*r) for r in q[[v,l,h]].to_numpy()]
    b=_table('ph_transfer_correction','Paired correction effects on conditional pH prediction',
        'Corrected-minus-raw RMSE and bias use identical physical-plate/standard pairs and shared reference '
        'bins from the transfer analysis. The technical-SD result is corrected/raw. Intervals use the same '
        '199-resample refitting bootstrap and preserve raw or corrected inverse failures. Positive RMSE '
        'change means worse prediction. Technical uncertainty conditions on recorded plates and references; '
        'the comparison does not establish performance for a newly prepared standard or culture.',q,
        ['acid_m','matrix_display','paired_n','raw_failure_n','corrected_failure_n','rmse_change_ci','bias_change_ci','sd_ratio_ci'],
        ['ph_matrix_readout_comparisons','ph_matrix_common_support'],['1b','1c','3h'],['Background subtraction can affect accuracy and technical precision in different directions.'])
    a['frame'], b['frame'], support = augment_ph_transfer(
        a['frame'], b['frame'], t['ph_matrix_common_support'])
    a['support_tables'] = {name: support[name] for name in ('accuracy_support','accuracy_references')}
    b['support_tables'] = {name: support[name] for name in ('correction_support','correction_references')}
    a['caption'] += (' Full exports report the actual retained discrete pH bounds, shared 0.25-pH-bin '
                     'counts, physical reference counts and source/plate identities separately for each '
                     'readout; these bounds do not imply continuous-range validation.')
    b['caption'] += (' Paired coverage is derived separately from single-readout accuracy coverage and '
                     'reported with actual pH bounds, shared bins, physical reference counts and '
                     'source/plate identities in the full and support exports.')
    return [a,b]


def _chemical_reference(t):
    q=t['feii_regression'].copy()
    q['slope_ci']=[_ci(*r) for r in q[['slope','slope_ci_low','slope_ci_high']].to_numpy()]
    q['intercept_ci_mM']=[_ci(*(1000*r)) for r in q[['intercept_M','intercept_ci_low_M','intercept_ci_high_M']].to_numpy()]
    q['residual_sd_mM']=1000*q.sy_x_M
    ps=q[['slope_vs_zero_p','slope_vs_one_p','intercept_vs_zero_p']].to_numpy().ravel()
    adj=_holm(pd.DataFrame({'p_value':ps})).p_holm.to_numpy().reshape((-1,3))
    for i,c in enumerate(['slope_zero_p_holm','slope_one_p_holm','intercept_zero_p_holm']):q[c]=adj[:,i]
    b=_table('chemical_reference_regression','Ferrous titration versus the recorded total-Fe reference',
        'Ordinary least squares regresses the mean FeII titration on the recorded total-Fe reference '
        'within each chemical series. The twelve reference levels, each summarized from three titrations, '
        'define the fitted relation; individual concentration measurements do not make these levels independent standard preparations. '
        'Slope tests are against 0 and 1 and intercept tests against 0; Holm correction covers all finite '
        'tests across the three series. Intervals condition on fixed references and the OLS residual model; '
        'uncertainty and covariance of the individually measured reference concentrations are not propagated. The constant ferric-series reading has '
        'zero fitted residual spread, so its inferential intervals and tests are withheld.',q,
        ['series_ID','n','slope_ci','intercept_ci_mM','residual_sd_mM','slope_zero_p_holm','slope_one_p_holm','intercept_zero_p_holm'],
        ['feii_regression'],['1g'],['Ferrous readings track the recorded total-Fe reference differently across chemical series.'])
    return b


def _iron_contrasts(t):
    q=t['iron_optical_reference_levels'];q=q[q.series_ID.eq('ferric') & q.signal.isin(['a380','a380_corrected'])].copy()
    rows=[]
    for condition,g in q.groupby('condition'):
        p=g.pivot(index=['physical_plate_id','plate_column'],columns='signal',values='mean')
        v=(p.a380_corrected-p.a380).groupby(level=0).mean()
        rows.append(dict(contrast='A*380 − A380',context='with dye' if condition=='thymol_blue' else 'without dye',**_summary(v)))
    for signal,g in q.groupby('signal'):
        p=g.pivot(index=['physical_plate_id','plate_column'],columns='condition',values='mean')
        v=(p.thymol_blue-p.no_thymol_blue).groupby(level=0).mean()
        rows.append(dict(contrast='dyed − undyed',context=SIGNAL[signal],**_summary(v)))
    dye_table=_table('iron_dye_correction','Paired ferric-signal changes with dye and A750 subtraction',
        'At each scheduled ferric reference level, paired differences are averaged equally within '
        'physical plate before computing the across-plate effect (AU). Dyed-minus-undyed contrasts '
        'compare recorded preparation positions; the dyed wells also have a different final volume '
        'and row block, so this is a combined formulation/volume/position contrast, not an isolated dye '
        'effect at identical concentration. '+TECHNICAL+'Two-sided tests against zero receive Holm '
        'adjustment across four contrasts. Negative absorbance differences are retained.',_holm(pd.DataFrame(rows)),
        ['contrast','context','plate_n','estimate_ci','p_holm'],['iron_optical_reference_levels'],['1h'],
        ['A750 subtraction and dyed-well preparation produce different paired optical offsets.'])
    # Same reference grid and interpolation as the Fig. 1e endpoint contrast.
    q=t['iron_optical_reference_levels'];q=q[q.condition.eq('no_thymol_blue')&q.signal.isin(['a380','a428'])]
    ferric=q[q.series_ID.eq('ferric')];ferrous=q[q.series_ID.eq('ferrous_1M')]
    target=np.sort(ferric.total_Fe_M.unique())
    vals={}
    for chemical,g in [('ferric',ferric),('ferrous',ferrous)]:
        for signal,h in g.groupby('signal'):
            r=[]
            for plate,v in h.groupby('physical_plate_id'):
                v=v.sort_values('total_Fe_M'); z=np.interp(target,v.total_Fe_M,v['mean'],left=np.nan,right=np.nan)
                r.append((plate,z))
            vals[(chemical,signal)]=r
    all_vectors=[v for a in vals.values() for _,v in a]
    keep=np.isfinite(np.asarray(all_vectors)).all(axis=0)
    rows=[]
    for signal in ['a380','a428','a380_minus_a428']:
        a=np.array([v[keep].mean() for _,v in vals[('ferric','a380' if signal=='a380_minus_a428' else signal)]])
        b=np.array([v[keep].mean() for _,v in vals[('ferrous','a380' if signal=='a380_minus_a428' else signal)]])
        if signal=='a380_minus_a428':
            a-=np.array([v[keep].mean() for _,v in vals[('ferric','a428')]])
            b-=np.array([v[keep].mean() for _,v in vals[('ferrous','a428')]])
        va=a.var(ddof=1)/len(a);vb=b.var(ddof=1)/len(b);se=np.sqrt(va+vb)
        df=(va+vb)**2/(va**2/(len(a)-1)+vb**2/(len(b)-1)) if va+vb>0 else np.nan
        effect=a.mean()-b.mean(); margin=student_t.ppf(.975,df)*se
        rows.append(dict(readout=SIGNAL.get(signal,'A380 − A428'),ferric_plate_n=len(a),ferrous_plate_n=len(b),
                         common_level_n=int(keep.sum()),estimate=effect,ci_low=effect-margin,ci_high=effect+margin,
                         estimate_ci=_ci(effect,effect-margin,effect+margin),reference_min_M=target[keep].min(),
                         reference_max_M=target[keep].max(),welch_df=df))
    b=_table('ferric_wavelength_selectivity','Ferric-minus-ferrous endpoint absorbance contrast',
        'For the undyed 1 M sulfuric-acid references in Figure 1e, each chemical family is reduced '
        'to one mean per physical plate across exactly the same supported ferric reference levels. '
        'Ferrous responses are interpolated only between recorded adjacent concentration levels, '
        'never in wavelength. Contrasts are ferric minus ferrous at 380 and 428 nm; the final row is '
        'the difference between those contrasts. Wavelength pairing is preserved within each plate; '
        'ferric and ferrous preparations are independent groups, so 95% Welch intervals combine '
        'their between-plate variances rather than pairing equal run numbers. Intervals condition '
        'on the chosen reference grid and shared prepared stocks. Effects are AU and summarize '
        'the entire common grid; they are not a detection limit or a local concentration sensitivity.',
        pd.DataFrame(rows),['readout','ferric_plate_n','ferrous_plate_n','common_level_n','estimate_ci'],
        ['iron_optical_reference_levels'],['1e'],['The native endpoint channels separate ferric and ferrous optical response.'])
    return [b,dye_table]


def _paired_effects(t):
    q=t['pub_robustness_paired_summary'].copy()
    q['normalized_effect']=100*q.mean_difference_fraction_of_reference_range
    q=q[~q.comparison.eq('nonUV_minus_UV') | q.step_id.eq('initial_23_C')].copy()
    specs=[('lid_effect','Lid-associated paired response changes',
            q[q.comparison.eq('lid_present_minus_absent') &
              (q.assay.eq('pH')|q.condition.eq('no_thymol_blue'))],
            ['assay','condition','signal'],['2a','3a'])]
    out=[]
    for key,title,data,groupcols,panels in specs:
        rows=[]
        for identity,g in data.groupby(groupcols,dropna=False,sort=True):
            meta=dict(zip(groupcols,identity if isinstance(identity,tuple) else (identity,)))
            # One summary per physical plate/preparation pair. Never let repeated wells become n.
            if g.candidate_physical_plate_id.duplicated().any():
                raise ValueError(f'duplicate physical plate in {key}: {meta}')
            row={**meta,**_summary(g.normalized_effect),
                 'paired_well_n':int(g.paired_well_n.sum()),'paired_level_n_min':int(g.paired_level_n.min()),
                 'paired_level_n_max':int(g.paired_level_n.max()),
                 'plate_ids':';'.join(g.candidate_physical_plate_id.astype(str))}
            row['readout']=SIGNAL[meta['signal']]
            if 'condition' in meta:row['dye']='10 µM' if meta['condition']=='thymol_blue' else '0 µM'
            rows.append(row)
        frame=_holm(pd.DataFrame(rows))
        display=[c for c in ['assay','acid_M','dye','readout','plate_material','contrast_label','plate_n','estimate_ci','p_holm'] if c in frame]
        caption=('Effects are 100×mean(lid-present−lid-absent)/span(lid-absent), in percent '
                 'of each paired reference’s optical range. Eligible well responses are paired at '
                 'the same preparation level; level contrasts receive equal weight within each '
                 'physical plate and do not increase independent n. These signed additive effects '
                 'are distinct from the figures’ shared-reference log2 heatmaps. '+TECHNICAL+
                 'Two-sided tests against zero use Holm correction over the four finite lid '
                 'contrasts. Small n limits inference.')
        out.append(_table(key,title,caption,frame,display,['pub_robustness_paired_summary'],panels,
                          ['Paired optical response changes with lid addition.']))
    return out


def _storage_edge(t):
    out=[];q=t['storage_folds'].copy();rows=[]
    keys=['assay','physical_plate_id','condition','signal']
    last=q.groupby(keys).read_index.transform('max');q=q[q.read_index.eq(last)]
    for identity,g in q.groupby(['assay','condition','signal']):
        assay,condition,signal=identity
        for metric,field,mult in [('signed change (% initial span)','normalized_difference',100),('log2(final/initial)','log2_fold',1)]:
            p=g.groupby('physical_plate_id',as_index=False).agg(value=(field,'mean'),days=('days','mean'))
            s=_summary(mult*p.value)
            rows.append(dict(assay=assay,dye='10 µM' if condition=='thymol_blue' else '0 µM',readout=SIGNAL[signal],
                             metric=metric,mean_endpoint_d=p.days.mean(),eligible_level_plate_n=int(np.isfinite(g[field]).sum()),
                             unavailable_level_plate_n=int((~np.isfinite(g[field])).sum()),**s))
    out.append(_table('storage_drift','Final recorded storage drift',
        'The final scheduled storage acquisition is compared with each physical plate/readout/level’s '
        'initial observation. Signed change is 100×(Xt−X0)/span(X0), averaged equally over eligible '
        'preparation levels within plate. Log2(final/initial) uses only positive finite operands and '
        'is averaged within plate before the across-plate summary; its exponent is a geometric fold. '
        'Both estimands are shown because a valid negative corrected signal cannot enter a logarithm. '
        'No last-available substitution, absolute-value transform or pseudocount is used. '+TECHNICAL+
        'Holm adjustment covers all finite tests in this table. This is an endpoint drift analysis, '
        'not proof of stability throughout the interval or equivalence.',_holm(pd.DataFrame(rows)),
        ['assay','dye','readout','metric','mean_endpoint_d','plate_n','estimate_ci','unavailable_level_plate_n','p_holm'],
        ['storage_folds'],['2b','2c','3b','3c','3d','3e'],['Storage changes optical signals and can invalidate positive-reference log folds.']))
    edge=t['edge_replicates'].copy();initial=t['storage_replicates'];initial=initial[initial.read_index.eq(0)]
    spans=initial.groupby(keys,as_index=False)['mean'].agg(lambda a:a.max()-a.min()).rename(columns={'mean':'initial_span'})
    edge=edge.merge(spans,on=keys,validate='many_to_one')
    edge['percent_drift']=100*edge.edge_minus_interior/edge.initial_span.where(edge.initial_span.gt(0))
    units=edge.groupby(keys+['read_index'],as_index=False).agg(value=('percent_drift','mean'),days=('days','mean'),eligible_level_n=('percent_drift','count'))
    rows=[]
    for identity,g in units.groupby(['assay','condition','signal']):
        assay,condition,signal=identity
        for metric in ['final acquisition','time-average over observed window']:
            values=[];endpoints=[];maximum=[]
            for plate,p in g.groupby('physical_plate_id'):
                p=p.sort_values('read_index');endpoints.append(p.days.iloc[-1]);maximum.append(p.value.abs().max())
                if metric=='final acquisition':v=p.value.iloc[-1]
                else:
                    complete=np.isfinite(p.value).all() and np.isfinite(p.days).all() and len(p)>1 and p.days.iloc[-1]>p.days.iloc[0]
                    v=np.trapezoid(p.value,p.days)/(p.days.iloc[-1]-p.days.iloc[0]) if complete else np.nan
                values.append(v)
            rows.append(dict(assay=assay,dye='10 µM' if condition=='thymol_blue' else '0 µM',readout=SIGNAL[signal],metric=metric,
                             mean_endpoint_d=np.mean(endpoints),maximum_absolute_observed_percent=np.nanmax(maximum),**_summary(values)))
    out.append(_table('edge_drift','Concentration-matched spatial drift relative to initial span',
        'The exact displayed estimand is 100×[(Xedge,t−Xedge,0)−(Xinterior,t−Xinterior,0)]/span(X0). '
        'Each level uses edge and interior wells at the same preparation concentration; levels without '
        'both positions are excluded from that contrast. Complete level contrasts are averaged within '
        'physical plate. The initial span uses all available initial level means in the same plate, '
        'dye and readout. Results are percent of initial span for both pH and iron. Final acquisition '
        'and actual-time trapezoidal mean drift summarize each plate before inference. The maximum '
        'absolute observed plate drift is descriptive and has no selection-based P test. '+TECHNICAL+
        'Holm correction covers both summaries and all readouts. Spatial drift does not identify '
        'evaporation as its cause or validate a perimeter-water intervention.',_holm(pd.DataFrame(rows)),
        ['assay','dye','readout','metric','plate_n','estimate_ci','maximum_absolute_observed_percent','p_holm'],
        ['edge_replicates','storage_replicates'],['2d','3f'],['Matched edge and interior signal changes differ relative to initial calibration span.']))
    return out


def _media(t):
    q=t['pub_media_performance'];q=q[(q.assay.eq('pH')&q.condition.eq('with_dye')) | (q.assay.eq('iron')&q.signal.isin(['A380','A380_corrected']))].copy()
    # A single eligible reference level has no estimable response span. The
    # upstream max-minus-min value of zero must not imply a valid flat response.
    q.loc[q.eligible_level_n.lt(2),['observed_response_range','range_to_technical_sd']]=np.nan
    rows=[]
    for (medium,signal),g in q.groupby(['medium_id','signal']):
        row=dict(medium=MEDIUM.get(medium,medium),readout=SIGNAL[signal],plate_n=g.physical_plate_id.nunique(),
                 span_mean=g.observed_response_range.mean(),technical_sd_mean=g.pooled_technical_sd.mean(),
                 span_sd_mean=g.range_to_technical_sd.mean(),median_rank_rho=g.spearman_rho.dropna().median() if g.spearman_rho.notna().any() else np.nan,
                 plates_with_adjacent_reversal=int(g.adjacent_direction_reversals.gt(0).sum()),
                 eligible_well_n=int(g.eligible_well_n.sum()),scheduled_well_n=int(g.scheduled_well_n.sum()))
        # A symmetric t interval can put an impossible negative limit on an
        # arithmetic mean of positive plate spans/SDs. Delta-method inference on
        # log(mean), followed by exponentiation, retains the arithmetic estimand.
        # This is an approximate small-sample interval, not a lognormal/geometric
        # mean fit; neither measurements nor limits are clipped.
        for prefix,values in [('span',g.observed_response_range),('sd',g.pooled_technical_sd)]:
            values=np.asarray(values,float);values=values[np.isfinite(values)]
            n=len(values);mean=float(values.mean()) if n else np.nan
            low=high=np.nan
            if n>1 and mean>0:
                se=float(values.std(ddof=1)/np.sqrt(n))
                half_width=float(student_t.ppf(.975,n-1)*se/mean)
                low=float(mean*np.exp(-half_width));high=float(mean*np.exp(half_width))
            row[prefix+'_plate_n']=n
            row[prefix+'_ci_low']=low;row[prefix+'_ci_high']=high
            row[prefix+'_ci']=_ci(mean,low,high)
        row['positive_mean_ci_method']='approximate_log_delta_t_on_arithmetic_mean'
        rows.append(row)
    a=_table('media_response_performance','Medium-specific optical response and technical scatter',
        'Results summarize each physical plate’s response span (max−min of replicated eligible level means; '
        'at least two eligible wells per level), '
        'pooled within-level technical SD and span/SD. Span and SD are in AU for A380 and in ratio '
        'units for R; the quotient is dimensionless. Plates receive equal weight. Median Spearman '
        'rho and the number of plates with a wrong-direction adjacent step describe monotonicity '
        'against the recorded reference (decreasing pH response, increasing iron response); no '
        'population-level monotonicity test is inferred. Span and SD retain their arithmetic plate means. '
        'Their approximate 95% intervals are mean×exp[±t(n−1)×SE(mean)/mean], using the delta method '
        'on log(mean). This positive-scale interval does not clip limits or change the estimand to a '
        'geometric mean. With typically three plates, it is exploratory and its coverage is not '
        'guaranteed; a single eligible plate has no interval. Numeric interval limits and the eligible '
        'plate count for each metric are retained in the machine-readable table. Reference/preparation '
        'uncertainty remains unquantified. Ferric-addition wells contain 0.5 M nominal added sulfuric '
        'acid after mixing; their response does not describe untreated medium. '
        'A response span and span/SD require at least two replicated eligible reference levels; the one-level '
        'corrected AFM4 acquisitions therefore have unavailable spans rather than zero-span responses. '
        'Eligibility requires finite unflagged channels and a corrected ratio denominator >10⁻⁶ AU. '
        'The count is a diagnostic accompanying analytical results, not evidence of analytical recovery.',pd.DataFrame(rows),
        ['medium','readout','plate_n','span_ci','sd_ci','span_sd_mean','median_rank_rho','plates_with_adjacent_reversal','eligible_well_n'],
        ['pub_media_performance','pub_media_observations'],['4a','4b','4c'],['Medium composition changes signal span, scatter and monotonicity.'])
    q=t['pub_media_readout_comparison'];q=q[(q.assay.eq('pH')&q.dye_uM.gt(0)&q.reference_signal.eq('ratio_raw')&q.candidate_signal.eq('ratio_corrected')) | (q.assay.eq('iron')&q.dye_uM.eq(0)&q.reference_signal.eq('A380')&q.candidate_signal.eq('A380_corrected'))].copy()
    ratio=q.candidate_range_to_sd_common_wells/q.reference_range_to_sd_common_wells
    q['log2_ratio']=np.log2(ratio.where(q.candidate_range_to_sd_common_wells.gt(0)&q.reference_range_to_sd_common_wells.gt(0)&np.isfinite(ratio)))
    rows=[]
    for (medium,signal),g in q.groupby(['medium_id','reference_signal']):
        s=_summary(g.log2_ratio); rows.append(dict(medium=MEDIUM.get(medium,medium),readout=SIGNAL[signal],
            scheduled_plate_n=g.physical_plate_id.nunique(),unavailable_plate_n=int((~np.isfinite(g.log2_ratio)).sum()),
            common_well_n=int(g.common_eligible_n.sum()),raw_Q_mean=g.reference_range_to_sd_common_wells.mean(),
            corrected_Q_mean=g.candidate_range_to_sd_common_wells.mean(),**s))
    b=_table('media_correction_effects','Paired correction effects on medium-specific span-to-SD',
        'The effect is Δlog2(Q)=log2(Qcorrected/Qraw), exactly matching Figure 4d; Q=span/pooled '
        'within-level sample SD. Each raw/corrected pair uses exactly the same eligible physical wells. '
        'One contrast per physical plate enters the mean, interval and two-sided zero test; Holm correction '
        'covers the seven estimable comparisons. '+TECHNICAL+'AFM4 has only one defined physical-plate '
        'contrast and therefore no CI or P value; the other two remain unavailable. The defined AFM4 '
        'contrast uses one eligible well at pH 1.39 and three at pH 3.79 on media plate 08. Its range '
        'therefore includes a singleton level, while its pooled SD uses the triplicate level alone. '
        'This sparse common-well contrast is distinct from the replicated-level performance analysis, '
        'which cannot estimate an AFM4 corrected span. Positive effects '
        'mean a larger optical span-to-scatter ratio and do not establish better pH or Fe recovery. '
        'Ferric-addition contrasts apply to medium with 0.5 M nominal added sulfuric acid after mixing.',
        _holm(pd.DataFrame(rows)),['medium','readout','plate_n','unavailable_plate_n','raw_Q_mean','corrected_Q_mean','estimate_ci','p_holm'],
        ['pub_media_readout_comparison','pub_media_observations'],['4d'],['A750 correction changes span-to-SD differently across media.'])
    support = build_media_analysis_support(t)
    a['support_tables'] = {'physical_plate_support': support['performance_plate_support']}
    b['support_tables'] = {'physical_plate_support': support['paired_plate_support']}
    a['caption'] += (' Physical-plate support exports give scheduled, individually eligible and '
                     'replicated reference-level counts and their actual bounds; a singleton does '
                     'not qualify as a replicated level for this analysis.')
    b['caption'] += (' Physical-plate support exports retain common-level bounds, singleton and '
                     'variance-contributing level counts, and unavailable contrasts separately from '
                     'the replicated-level performance population.')
    return [a,b]


def build(tables):
    """Return complete Figure 1–4 analyses-only SI specifications."""
    models=_models(tables)
    out=[models[0],_initial_precision(tables),*_transfer(tables),models[1],
         _chemical_reference(tables),*_iron_contrasts(tables),*_paired_effects(tables),
         *_storage_edge(tables),*_media(tables)]
    assert len({x['key'] for x in out})==len(out)
    for item in out:
        assert len(item['frame']) and all(c in item['frame'] for c in item['display_columns']),item['key']
    return out
