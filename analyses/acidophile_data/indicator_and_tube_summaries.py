"""Summarize indicator calibration, storage stability and tube measurements."""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy.stats import t as student_t, linregress
from .publication_calibration import _core
from .publication_calibration_fourpl import fit_fourpl, _jacobian, fit_ferric_response, ferric_response_jacobian


def fitted_response(frame, x, y, identity, *, fixed_bottom=None):
    q = frame[np.isfinite(frame[x]) & np.isfinite(frame[y])].copy()
    iron_model = identity.get('scope', '').startswith('iron_')
    if iron_model:
        q = q[q[x].ge(0)].copy()
    fit = fit_ferric_response(q[x], q[y]) if iron_model else fit_fourpl(q[x], q[y], fixed_bottom=fixed_bottom)
    grid = np.linspace(q[x].min(), q[x].max(), 181) if len(q) else np.array([np.nan])
    if iron_model and q[x].gt(0).any():
        grid = np.r_[0., np.geomspace(q.loc[q[x].gt(0), x].min(), q[x].max(), 181)]
    mean = fit.predict(grid)
    low = high = np.full(len(grid), np.nan)
    ids = q.physical_plate_id.to_numpy() if 'physical_plate_id' in q else np.arange(len(q))
    groups = np.unique(ids)
    if len(groups) > 1 and np.isfinite(fit.influence).all() and fit.converged:
        influences = np.array([fit.influence[ids == unit].sum(axis=0) for unit in groups])
        covariance = influences.T @ influences * len(groups)/(len(groups)-1)
        j = ferric_response_jacobian(grid, fit.estimates) if iron_model else _jacobian(grid, fit.estimates)
        se = np.sqrt(np.maximum(0, np.einsum('ij,jk,ik->i', j,covariance,j)))
        margin = student_t.ppf(.975,len(groups)-1)*se
        low, high = mean-margin, mean+margin
    curves = pd.DataFrame(dict(x=grid, fitted=mean, ci_low=low, ci_high=high))
    for k,v in identity.items(): curves[k] = v
    model = {**identity, **fit.as_dict(), 'physical_plate_n':len(groups),
             'predictor_unit':'mM total Fe' if iron_model else 'measured pH',
             'model_family':'blank_plus_stretched_exponential' if iron_model else 'base10_pH_4PL',
             'slope_parameter_unit':'dimensionless' if iron_model else '1/pH',
             'reference_scope':'individually_measured_total_Fe_uncertainty_not_propagated' if iron_model else 'mean_of_three_recorded_pH_references',
             # The plotting helper draws the fitted-response coordinates even
             # for weak models, but never invents coordinates for missing fits.
             'displayed':bool(np.any(np.isfinite(grid) & np.isfinite(mean))),
             'response_interval_method':'physical_plate_clustered_HC3_delta_t',
             'interval_df':len(groups)-1, 'fixed_zero_lower_asymptote':False if iron_model else fixed_bottom == 0}
    return curves, model


def spectra(tables):
    f = tables['calibration_spectral_summary']
    f = f[f.experiment_id.eq('pH_UV_200uM_stock_initial_spectra')].copy()
    refs = f[['physical_plate_id','row_block','plate_column','reference_pH']].drop_duplicates()
    ph = refs.groupby(['row_block','plate_column'],as_index=False).reference_pH.mean().rename(columns={'reference_pH':'mean_measured_pH'})
    f = f.merge(ph,on=['row_block','plate_column'],validate='many_to_one')
    s = f.groupby(['row_block','plate_column','mean_measured_pH','wavelength_nm'], as_index=False).agg(
        mean=('mean_absorbance','mean'), sd=('mean_absorbance','std'), plate_n=('mean_absorbance','count'))
    return {'spectral_replicates':f, 'spectral_means':s}


def calibration(tables):
    levels = _core(tables['calibration_level_means'])
    ph = levels[levels.assay.eq('pH') & levels.condition.eq('thymol_blue')].copy()
    ph['mean_measured_pH'] = ph.groupby(['experiment_id','signal','row_block','plate_column']).reference_pH.transform('mean')
    models, curves = [], []
    for (experiment,signal),g in ph.groupby(['experiment_id','signal']):
        c,m = fitted_response(g,'mean_measured_pH','mean',dict(scope='pH',experiment_id=experiment,
            signal=signal,nominal_dye_uM=g.nominal_dye_uM.iloc[0]))
        curves.append(c);models.append(m)
    iron = levels[levels.assay.eq('iron') & levels.signal.eq('a380') & levels.condition.eq('no_thymol_blue')].copy()
    # The manuscript reports the acid-specific calibration curves directly.
    # Their held-out predictive benefit is assessed separately; training-error
    # hypothesis tests do not select which calibration is displayed.
    for acid, group in iron.groupby('acid_M'):
        c,m = fitted_response(group, 'axis_value', 'mean', dict(scope='iron_acid_specific',
            signal='a380', experiment_id='iron_all_acid_matrices', acid_M=acid))
        curves.append(c); models.append(m)
    return {'ph_replicates':ph, 'response_curves':pd.concat(curves,ignore_index=True),
            'response_models':pd.DataFrame(models)}


def iron_spectra(tables):
    f = tables['calibration_spectral_summary']
    ferric = f[f.experiment_id.eq('final_indicator_ferric') & f.condition.eq('no_thymol_blue')]
    ferrous = f[f.experiment_id.eq('final_indicator_ferrous_1M') & f.condition.eq('no_thymol_blue')]
    rows = []
    refs = dict(tuple(ferrous.groupby('wavelength_nm')))
    for _,r in ferric.iterrows():
        g = refs.get(r.wavelength_nm)
        if g is None:continue
        target = r.reference_concentration_M
        reference_values=[]
        for _,preparation in g.groupby('physical_plate_id'):
            preparation=preparation.sort_values('reference_concentration_M')
            x,y=preparation.reference_concentration_M.to_numpy(),preparation.mean_absorbance.to_numpy()
            reference_values.append(np.interp(target,x,y) if len(x) and np.isfinite(x).all()
                                    and x[0]<=target<=x[-1] else np.nan)
        reference=np.mean(reference_values) if np.isfinite(reference_values).all() else np.nan
        reference_sd=np.std(reference_values,ddof=1) if len(reference_values)>1 else np.nan
        rows.append(dict(physical_plate_id=r.physical_plate_id,run_number=r.run_number,
            plate_column=r.plate_column,total_Fe_M=target,wavelength_nm=r.wavelength_nm,
            ferric_mean=r.mean_absorbance,ferrous_reference=reference,
            ferrous_reference_sd=reference_sd,ferrous_plate_n=len(reference_values),
            difference=r.mean_absorbance-reference,source_path=r.source_path,
            ferrous_source_path=';'.join(sorted(g.source_path.unique()))))
    return pd.DataFrame(rows)


def reference_regression(tables):
    refs = tables['iron_chemical_reference_summary']; models,curves=[],[]
    for series,g in refs.groupby('series_ID'):
        q = g.dropna(subset=['total_Fe_M','FeII_M'])
        x,y = q.total_Fe_M.to_numpy(),q.FeII_M.to_numpy()
        fit = linregress(x,y); df = len(x)-2; crit=student_t.ppf(.975,df)
        residual = y-fit.intercept-fit.slope*x; syx=np.sqrt(np.sum(residual**2)/df)
        models.append(dict(series_ID=series,n=len(x),slope=fit.slope,intercept_M=fit.intercept,
            slope_se=fit.stderr,intercept_se_M=fit.intercept_stderr,
            slope_ci_low=fit.slope-crit*fit.stderr,slope_ci_high=fit.slope+crit*fit.stderr,
            intercept_ci_low_M=fit.intercept-crit*fit.intercept_stderr,
            intercept_ci_high_M=fit.intercept+crit*fit.intercept_stderr,
            slope_vs_zero_p=fit.pvalue,
            slope_vs_one_p=2*student_t.sf(abs((fit.slope-1)/fit.stderr),df),
            intercept_vs_zero_p=2*student_t.sf(abs(fit.intercept/fit.intercept_stderr),df),
            residual_df=df,sy_x_M=syx,x_min_M=x.min(),x_max_M=x.max()))
        grid=np.linspace(x.min(),x.max(),101); mean=fit.intercept+fit.slope*grid
        margin=crit*syx*np.sqrt(1/len(x)+(grid-x.mean())**2/np.sum((x-x.mean())**2))
        curves.append(pd.DataFrame(dict(series_ID=series,x=grid,fitted=mean,ci_low=mean-margin,ci_high=mean+margin)))
    return {'feii_regression':pd.DataFrame(models),
            'feii_regression_curves':pd.concat(curves,ignore_index=True)}


def storage(tables):
    f=tables['pub_robustness_observations']
    f=f[f.experiment_id.isin(['final_indicator_ferric','final_indicator_acid'])].copy()
    f['recorded_datetime']=pd.to_datetime(f.recorded_time)
    f['days']=(f.recorded_datetime-f.groupby('physical_plate_id').recorded_datetime.transform('min')).dt.total_seconds()/86400
    keys=['assay','physical_plate_id','condition','signal','row_group','column','read_index']
    level=f.groupby(keys,as_index=False).agg(days=('days','mean'),predictor=('predictor','mean'),
        mean=('response','mean'),sd=('response','std'),well_n=('response','count'))
    base=f[f.read_index.eq(0)][['physical_plate_id','signal','well','response']].rename(columns={'response':'initial'})
    f=f.merge(base,on=['physical_plate_id','signal','well'],how='left',validate='many_to_one')
    f['drift']=f.response-f.initial
    f['edge']=f.row.isin(['A','H']) | f.column.isin([1,12])
    rows=[]
    for key,g in f.groupby(keys):
        e=g.loc[g.edge,'drift'];i=g.loc[~g.edge,'drift']
        rows.append(dict(zip(keys,key)) | dict(days=g.days.mean(),predictor=g.predictor.mean(),
            edge_n=e.count(),interior_n=i.count(),edge_minus_interior=e.mean()-i.mean(),
            status='estimated' if e.count() and i.count() else 'no_same_level_interior_comparator'))
    return {'storage_replicates':level,'edge_replicates':pd.DataFrame(rows)}


def tube_od(tables):
    f=tables['validation_tube_points'].copy()
    rfu=pd.to_numeric(f.planktonic_background_corrected_rfu,errors='coerce')
    valid=rfu.ge(0) & np.isfinite(rfu)
    # Only source-observed fluorescence is eligible. A shared initial record
    # is retained, but cannot become an independent endpoint in a boxplot.
    if 'fluorescence_status' in f:
        valid &= f.fluorescence_status.fillna('').str.lower().isin(['reported','eligible','observed','measured'])
    f['apparent_OD600']=(rfu.where(valid)/74930)**(1/1.386)
    f['od_status']=np.where(valid,'approximate_AF_calibration',np.where(rfu.lt(0),'negative_RFU_no_real_calibrated_OD','missing_or_ineligible_RFU'))
    return f
