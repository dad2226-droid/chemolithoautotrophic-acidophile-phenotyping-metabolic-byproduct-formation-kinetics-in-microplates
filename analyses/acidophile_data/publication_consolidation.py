"""Consolidation."""
from __future__ import annotations
import json
import numpy as np
import pandas as pd
from .plate_comparisons import CONDITION_KEYS
from .peak_response import observed_peak, compare_peaks


def no_addition(value):
    return all(float(v) == 0 for v in json.loads(value).values())


def activity_tables(tables):
    pc = tables['af_iron_plate_trajectories']
    pc = pc[pc.signal.eq('raw_A380')]
    plate = dict(tuple(pc.groupby(CONDITION_KEYS + ['batch_id','plate_id'])))
    ppeaks = {k: observed_peak(g) for k, g in plate.items()}
    formulations = tables['cross_format_plate_formulations']
    zero = formulations[formulations.comparable_quantities_mM_json.map(no_addition)]
    raw_flask = tables['validation_flask_points']
    flask = dict(tuple(raw_flask[raw_flask.kind.eq('biotic')].groupby('flask_ID')))
    fpeaks = {k: observed_peak(g, value='net_change') for k, g in flask.items()}
    maps = tables['cross_format_condition_correspondence']
    maps = maps[maps.flask_kind.eq('biotic') & maps.match_status.isin(
        ['equal_reported_dose','nearest_lower_dose_context','nearest_upper_dose_context'])]
    replicate_rows, comparison_rows = [], []
    for _, m in maps.iterrows():
        key = tuple(m[k] for k in CONDITION_KEYS)
        group = [k for k in plate if k[:4] == key]
        if not group:
            continue
        fgroups = [g for g in flask.values() if g.condition_ID.iloc[0] == m.flask_condition_ID]
        fcontrols = [g for g in flask.values() if g.batch.iloc[0] == m.flask_batch and g.no_addition.iloc[0]]
        pq, fq = json.loads(m.comparable_quantities_mM_json), json.loads(m.flask_quantities_mM_json)
        additive = 'citric_acid' if 'citric_acid' in pq else next(iter(pq), 'control')
        if not pq:
            continue
        z = zero[(zero.experiment_id == m.experiment_id) & (zero.group_id == m.group_id)]
        for batch in sorted({k[4] for k in group}):
            ts = [k for k in group if k[4] == batch]
            cid = f'{m.correspondence_id}/{batch}'
            base = dict(comparison_id=cid, **dict(zip(CONDITION_KEYS,key)), plate_batch=batch,
                        flask_batch=m.flask_batch, flask_condition_ID=m.flask_condition_ID,
                        additive=additive, plate_dose_mM=pq.get(additive,0),
                        flask_dose_mM=fq.get(additive,np.nan),
                        Ga_mM=pq.get('Ga',0) if additive == 'citric_acid' else 0,
                        match_status=m.match_status, integration_endpoint='each_trajectory_own_observed_maximum',
                        abiotic_subtraction=False)
            rows = []
            for pk in ts:
                refs = [(*key[:3], col, batch, pk[-1]) for col in z.column.unique()]
                result = compare_peaks(ppeaks[pk], [ppeaks.get(k) for k in refs])
                rows.append({**base, **result, 'format':'plate', 'replicate_id':pk[-1],
                    'inoculum_culture_id':';'.join(sorted(set(';'.join(plate[pk].culture_ids).split(';')))),
                    'reference_ids':';'.join(str(k[3]) for k in refs),
                    'abiotic_status':'biotic_only_matched_abiotic_retained_for_validation'})
            for t in fgroups:
                identity = t.flask_ID.iloc[0]
                result = compare_peaks(fpeaks[identity], [fpeaks[g.flask_ID.iloc[0]] for g in fcontrols])
                rows.append({**base, **result, 'format':'flask', 'replicate_id':identity,
                    'inoculum_culture_id':t.inoculum_culture_id.iloc[0],
                    'reference_ids':';'.join(g.flask_ID.iloc[0] for g in fcontrols),
                    'abiotic_status':'biotic_only_no_abiotic_subtraction'})
            replicate_rows.extend(rows)
            for fmt in ('plate', 'flask'):
                q = pd.DataFrame([r for r in rows if r['format'] == fmt])
                for value, prefix in [('activity_percent',fmt), ('treated_time_to_max_d',fmt+'_time_to_max'),
                                      ('time_to_max_delay_d',fmt+'_delay')]:
                    values = q[value].dropna() if len(q) else pd.Series(dtype=float)
                    if value == 'treated_time_to_max_d' and len(q):
                        values = q.loc[q.treated_positive_rise, value].dropna()
                    base.update({prefix+'_mean':values.mean(), prefix+'_sd':values.std(), prefix+'_n':len(values)})
            base['activity_difference_pp'] = base['plate_mean'] - base['flask_mean']
            base['time_to_max_difference_d'] = base['plate_time_to_max_mean'] - base['flask_time_to_max_mean']
            base['delay_difference_d'] = base['plate_delay_mean'] - base['flask_delay_mean']
            base['status'] = 'estimated' if np.isfinite(base['activity_difference_pp']) else 'activity_unavailable'
            base['unavailable_reasons'] = ';'.join(sorted({r['status'] for r in rows if r['status'] != 'estimated'}))
            comparison_rows.append(base)
    return {'consolidated_activity_replicates':pd.DataFrame(replicate_rows),
            'consolidated_activity_comparisons':pd.DataFrame(comparison_rows)}
