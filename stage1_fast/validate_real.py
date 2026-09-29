"""Check the fast Stage 1 against an existing published-Stage-1 table on the same frames.

Example (EGFR system 1; run fast_stage1.py first with --frames on the published frames):
  python fast_stage1.py ... --frames $(seq 20000 20 20400) --out val_EGFR1 --workers 4 --chunk 5
  python merge_and_export.py --out val_EGFR1
  python validate_real.py --new val_EGFR1/causal_data_legacy_subset.csv \
      --old .../with_target_lipid/dpg3_causal_data_all_proteins.csv
The legacy table must agree exactly (binding state, every <L>_contacts, every
<L>_unique_molecules). The script prints every disagreement it finds.
"""
import argparse, sys
import pandas as pd

ap = argparse.ArgumentParser(); ap.add_argument('--new', required=True); ap.add_argument('--old', required=True)
a = ap.parse_args()
new = pd.read_csv(a.new); old = pd.read_csv(a.old)
m = new.merge(old, on=['frame', 'protein'], suffixes=('_new', '_old'))
if m.empty:
    sys.exit('no common (frame, protein) rows')
cols = [c for c in old.columns if c.endswith('_contacts') or c.endswith('_unique_molecules') or c == 'target_lipid_bound']
bad = 0
for c in cols:
    if c + '_new' not in m:
        print('column missing in new table:', c); continue
    x = m[c + '_new'].astype(float); y = m[c + '_old'].astype(float)
    d = m[x != y]
    if len(d):
        bad += len(d)
        print(f'{c}: {len(d)} of {len(m)} rows differ, e.g.')
        print(d[['frame', 'protein', c + '_new', c + '_old']].head(5).to_string(index=False))
print(f'{len(m)} rows compared, {bad} disagreements')
