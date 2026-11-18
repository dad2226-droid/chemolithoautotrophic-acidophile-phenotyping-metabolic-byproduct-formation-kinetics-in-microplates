"""Additional channels."""
from __future__ import annotations
import pandas as pd
from .biological import load_biological_observations
from .plate_comparisons import CONDITION_KEYS, CROSS_FORMAT_EXPERIMENTS, CITRATE_EXPERIMENT, derive_plate_signals, compare_biotic_with_controls
from .peak_response import plate_means


def extract(repo_root):
    ids = sorted(set(CROSS_FORMAT_EXPERIMENTS) | {CITRATE_EXPERIMENT})
    f = load_biological_observations(repo_root, species='AF', experiment_ids=ids)
    # The comparison library needs four base channels; substitute the recorded
    # alternative channel only in a working copy, never in the source observation table.
    pieces = []
    for wave, names in ((428, ('raw_A428','A428_minus_A750')), (600, ('raw_A600',))):
        q = f.copy()
        q['A380'], q['A380_eligible'] = q[f'A{wave}'], q[f'A{wave}_eligible']
        sig = derive_plate_signals(q)
        mapping = {'raw_A380':names[0]}
        if wave == 428: mapping['A380_minus_A750'] = names[1]
        sig = sig[sig.signal.isin(mapping)].copy(); sig['signal'] = sig.signal.map(mapping)
        trajectories, _ = compare_biotic_with_controls(sig)
        pieces.append(plate_means(trajectories))
    cols = list(dict.fromkeys(CONDITION_KEYS+['batch_id','plate_id','culture_id','row','well',
        'source','ordinal','days','kind','physical_observation_key'] +
        [k for k in f if k.startswith(('A380','A428','A600','A750'))]))
    return {'channel_observations':f[cols],
            'channel_trajectories':pd.concat(pieces, ignore_index=True)}
