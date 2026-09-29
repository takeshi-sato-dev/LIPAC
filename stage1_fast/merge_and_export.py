"""Merge the chunk files of fast_stage1.py and write Stage 2 input tables.

Outputs in <out>/:
  stage1_fast_all.npz            every per-frame array, frames, time, residue ids, lipid types
  causal_data_legacy_subset.csv  published definition (legacy filters, residue subset): reproduces
                                 dpg3_causal_data_all_proteins.csv of the published Stage 1
                                 when the same frames, leaflet and subset are used
  causal_data_legacy_full.csv    published definition over the full protein (GitHub version)
  causal_data_exact_full.csv     exact contacts over the full protein, correct molecule counts
  causal_data_exact_subset.csv   exact contacts over the residue subset
Column layout follows the published file (frame, protein, target_lipid_bound, <L>_contacts,
<L>_unique_molecules for every non-target lipid), with extra columns time_ps,
<TARGET>_contacts and <TARGET>_molecules (continuous measures of target-lipid contact).
"""
import os, sys, glob, json, argparse
import numpy as np, pandas as pd

ap = argparse.ArgumentParser(); ap.add_argument('--out', required=True); a = ap.parse_args()
meta = json.load(open(os.path.join(a.out, 'meta.json')))
files = sorted(glob.glob(os.path.join(a.out, 'chunks', 'chunk_*.npz')))
if not files:
    sys.exit('no chunks')
D = {}
for f in files:
    z = np.load(f)
    for k in z.files:
        D.setdefault(k, []).append(z[k])
D = {k: np.concatenate(v) for k, v in D.items()}
o = np.argsort(D['frames']); D = {k: v[o] for k, v in D.items()}
types = meta['lipids']; tgt = meta['target']; ti = types.index(tgt) if tgt in types else None
resids = np.array(meta['resids'])
subset = None
if meta['subset']:
    lo, hi = map(int, meta['subset'].split(':')); subset = (resids >= lo) & (resids <= hi)
np.savez_compressed(os.path.join(a.out, 'stage1_fast_all.npz'), resids=resids, lipids=np.array(types),
                    segids=np.array(meta['segids']), **D)
nf, nc = D['res_exact'].shape[:2]
expected = meta['n_frames_analysed']
print(f'{nf} frames merged (expected {expected}); copies {nc}; residues {len(resids)}')


def table(res, mol, mask):
    r = res[:, :, mask, :].sum(2)                        # frames x copies x types
    rows = []
    for c in range(nc):
        d = {'frame': D['frames'], 'time_ps': D['time_ps'], 'protein': f'Protein_{c+1}'}
        bound = (r[:, c, ti] > 0) if ti is not None else np.zeros(nf, bool)
        d['target_lipid_bound'] = bound
        for k, t in enumerate(types):
            if k == ti:
                continue
            d[f'{t}_contacts'] = r[:, c, k].astype(float)
            d[f'{t}_unique_molecules'] = mol[:, c, k].astype(int)
        if ti is not None:
            d[f'{tgt}_contacts'] = r[:, c, ti].astype(float)
            d[f'{tgt}_molecules'] = mol[:, c, ti].astype(int)
        rows.append(pd.DataFrame(d))
    return pd.concat(rows).sort_values(['protein', 'frame'])


full = np.ones(len(resids), bool)
out = {'legacy_full': table(D['res_legacy'], D['mol_legacy'] if subset is None else D['mol_full'], full),
       'exact_full': table(D['res_exact'], D['mol_full'], full)}
if subset is not None:
    out['legacy_subset'] = table(D['res_legacy'], D['mol_legacy'], subset)
    out['exact_subset'] = table(D['res_exact'], D['mol_subset'], subset)
for k, df in out.items():
    p = os.path.join(a.out, f'causal_data_{k}.csv'); df.to_csv(p, index=False)
    print(f'{p}: bound fraction {df.target_lipid_bound.mean():.3f}')
