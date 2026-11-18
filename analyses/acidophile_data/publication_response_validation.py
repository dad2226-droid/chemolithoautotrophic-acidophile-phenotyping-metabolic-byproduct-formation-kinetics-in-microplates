"""Response validation."""
from __future__ import annotations
import json
import numpy as np
import pandas as pd
from .plate_comparisons import CONDITION_KEYS
from .peak_response import observed_peak, compare_peaks


def crossing(curve, start, end):
    g = curve[curve.days_mean.between(start, end)].sort_values('days_mean')
    valid = g.biotic_change.notna() & g.abiotic_change.notna()
    hit = valid & g.biotic_change.ge(1) & g.abiotic_change.lt(1)
    if not hit.any():
        return np.nan, np.nan, False
    index = int(np.flatnonzero(hit)[0])
    return (float(g.days_mean.iloc[index-1]) if index and valid.iloc[index-1] else np.nan,
            float(g.days_mean.iloc[index]), True)


def inhibition_tables(tables):
    conditions = tables['af_iron_condition_coverage']
    curve = tables['af_iron_plate_trajectories']
    curve = curve[curve.signal.eq('raw_A380')]
    metadata = conditions.set_index(CONDITION_KEYS)
    keys = CONDITION_KEYS + ['batch_id', 'plate_id']
    grouped = dict(tuple(curve.groupby(keys)))
    controls = {}
    for key, row in metadata.iterrows():
        if row.principal_final_mM != 0:
            continue
        values = json.loads(row.component_final_mM_json)
        values.pop(row.principal_component, None)
        controls.setdefault((*key[:3], json.dumps(values, sort_keys=True)), []).append(key[3])
    peaks = {key: observed_peak(g) for key, g in grouped.items()}
    records, peak_rows = [], []
    for key, t in grouped.items():
        row = metadata.loc[key[:4]]
        base = {**dict(zip(keys, key)),
                **row[['additive_label','additive_class','principal_final_mM']].to_dict(),
                'culture_ids': ';'.join(sorted(set(';'.join(t.culture_ids).split(';')))),
                'signal': 'biotic_delta_A380', 'integral_unit': 'AU day',
                'abiotic_subtraction': False}
        peak_rows.append({**base, **peaks[key]})
        values = json.loads(row.component_final_mM_json)
        values.pop(row.principal_component, None)
        columns = controls.get((*key[:3], json.dumps(values, sort_keys=True)), [])
        reference_keys = [(*key[:3], col, *key[4:]) for col in columns]
        result = compare_peaks(peaks[key], [peaks.get(k) for k in reference_keys])
        active = bool(reference_keys) and all(
            k in grouped and crossing(grouped[k], peaks[k]['initial_time_d'],
                                      peaks[k]['maximum_time_d'])[2]
            for k in reference_keys)
        records.append({**base, **result, 'control_active': active,
                        'control_columns': ';'.join(str(k[3]) for k in reference_keys),
                        'reference_plate_ids': ';'.join(sorted({k[-1] for k in reference_keys}))})
    plates = pd.DataFrame(records)
    plates['rising_time_to_max_d'] = plates.treated_time_to_max_d.where(plates.treated_positive_rise)
    batches = plates.groupby(CONDITION_KEYS + ['batch_id','additive_label','additive_class',
        'principal_final_mM'], as_index=False).agg(
        mean=('integral_inhibition_percent','mean'), sd=('integral_inhibition_percent','std'),
        n=('integral_inhibition_percent','count'),
        time_to_max_mean_d=('rising_time_to_max_d','mean'),
        time_to_max_sd_d=('rising_time_to_max_d','std'),
        time_to_max_n=('rising_time_to_max_d','count'),
        time_to_max_delay_mean_d=('time_to_max_delay_d','mean'),
        time_to_max_delay_sd_d=('time_to_max_delay_d','std'),
        timing_n=('time_to_max_delay_d','count'),
        control_active_n=('control_active','sum'),
        final_read_maximum_n=('treated_maximum_at_final_read','sum'),
        plate_ids=('plate_id', lambda x: ';'.join(sorted(set(x)))))
    return dict(validation_inhibition_plates=plates, validation_inhibition_batches=batches,
                validation_plate_peak_metrics=pd.DataFrame(peak_rows))


def flask_tables(tables, vessel_tables):
    f = vessel_tables['flask_feii_timecourse'].copy()
    f['net_change'] = f.initial_observed_feii_mM - f.eligible_feii_mM
    f['biotic_change'] = f['net_change']
    f['days_mean'] = f.elapsed_days_mean
    quantities = tables['cross_format_condition_correspondence'].dropna(
        subset=['flask_condition_ID']).drop_duplicates('flask_condition_ID').set_index(
        'flask_condition_ID').flask_quantities_mM_json.to_dict()
    f['quantities_json'] = f.condition_ID.map(quantities)
    additions = [c for c in f if c.startswith('nominal_') and c.endswith('_mM')]
    f['no_addition'] = f[additions].fillna(0).eq(0).all(axis=1)
    grouped = dict(tuple(f[f.kind.eq('biotic')].groupby('flask_ID')))
    peaks = {identity: observed_peak(g) for identity, g in grouped.items()}
    records, peak_rows = [], []
    for flask, t in grouped.items():
        first = t.iloc[0]
        base = dict(condition_ID=first.condition_ID, flask_ID=flask, batch=first.batch,
                    signal='biotic_apparent_FeIII', integral_unit='mM day', abiotic_subtraction=False,
                    quantities_json=first.quantities_json, nominal_citrate_mM=first.nominal_citrate_mM,
                    nominal_Ga_mM=first.nominal_Ga_mM)
        peak_rows.append({**base, **peaks[flask]})
        for scope in ['general', 'citrate']:
            if scope == 'citrate' and first.condition_ID not in [f'FV{i:02d}' for i in range(21,34)]:
                continue
            reference = f[f.kind.eq('biotic') & f.batch.eq(first.batch)]
            reference = (reference[reference.no_addition] if scope == 'general' else
                reference[reference.nominal_citrate_mM.eq(0) & reference.nominal_Ga_mM.eq(first.nominal_Ga_mM)])
            ids = sorted(reference.flask_ID.unique())
            result = compare_peaks(peaks[flask], [peaks[k] for k in ids])
            records.append({**base, **result, 'scope': scope, 'window': 'each_trajectory_own_observed_maximum',
                            'reference_ids': ';'.join(ids),
                            'reference_condition_ids': ';'.join(sorted(reference.condition_ID.unique()))})
    return dict(validation_flask_points=f, validation_flask_inhibition=pd.DataFrame(records),
                validation_flask_peak_metrics=pd.DataFrame(peak_rows))


def tube_points(tubes):
    f=tubes.copy();f.attrs={}
    f=f.sort_values(['tube_ID','sampling_event'])
    initial=f.groupby('tube_ID').head(1).set_index('tube_ID')
    f['delta_pH']=f.pH-f.tube_ID.map(initial.pH)
    f['delta_ORP']=f.orp_mV-f.tube_ID.map(initial.orp_mV)
    keys=['medium_ID','treatment_ID','incubation_temperature_C','sampling_event']
    abiotic=f[f.kind.eq('abiotic')].groupby(keys,as_index=False).agg(abiotic_delta_ORP=('delta_ORP','mean'),abiotic_n=('delta_ORP','count'))
    f=f.merge(abiotic,on=keys,how='left',validate='many_to_one')
    f['delta_delta_ORP']=f.delta_ORP-f.abiotic_delta_ORP
    f['observed']=~f.is_shared_initial_reference
    return f


def sulfur_signal_to_noise(trajectories):
    """Keep direction, magnitude, eligibility and shared readout support distinct.

    Absolute values are taken for each starter/plate/acquisition observation,
    never after averaging signed responses that could cancel one another.
    """
    p = trajectories.copy()
    keys = CONDITION_KEYS + ['plate_id', 'batch_id', 'culture_id', 'row', 'signal']
    initial = (p.loc[p.days.eq(0), keys + ['abiotic_sd']]
               .rename(columns={'abiotic_sd': 'initial_abiotic_sd'}))
    p = p.merge(initial, on=keys, how='left', validate='many_to_one')
    p['change_noise'] = np.hypot(p.abiotic_sd, p.initial_abiotic_sd)
    noise_ok = np.isfinite(p.change_noise) & p.change_noise.gt(0)
    response_ok = np.isfinite(p.net_change)
    p['technical_snr_signed'] = p.net_change.where(response_ok) / p.change_noise.where(noise_ok)
    p['technical_snr_magnitude'] = p.technical_snr_signed.abs()
    # The unqualified technical signal-to-noise field contains the signed value.
    p['technical_snr'] = p.technical_snr_signed
    p['noise_status'] = np.select(
        [~noise_ok, ~response_ok],
        ['missing_zero_or_nonfinite_abiotic_noise', 'missing_or_nonfinite_net_change'],
        default='independent_read_noise_approximation')
    observation_keys = [key for key in keys if key != 'signal'] + ['ordinal']
    paired = p.pivot(index=observation_keys, columns='signal', values='technical_snr_signed')
    paired = paired.reindex(columns=['raw_R545_485', 'corrected_R545_485'])
    common = np.isfinite(paired).all(axis=1).rename('technical_snr_pair_eligible').reset_index()
    p = p.merge(common, on=observation_keys, how='left', validate='many_to_one')
    for definition in ['signed', 'magnitude']:
        p[f'technical_snr_{definition}_common'] = p[f'technical_snr_{definition}'].where(
            p.technical_snr_pair_eligible)
    return p


def consecutive_snr_summary(frame, value):
    """Selected adjacent deposited-acquisition minimum; exclude invalid pairs."""
    g = frame.sort_values('ordinal')
    values = g[value].to_numpy(float)
    adjacent = np.diff(g.ordinal.to_numpy()) == 1
    valid = adjacent & np.isfinite(values[:-1]) & np.isfinite(values[1:])
    pairs = np.minimum(values[:-1], values[1:])[valid]
    return (float(pairs.max()) if len(pairs) else np.nan, int(len(pairs)))
