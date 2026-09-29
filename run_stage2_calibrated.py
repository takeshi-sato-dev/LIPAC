#!/usr/bin/env python3
"""Calibrated Stage 2 for a Stage 1 table (dpg3_causal_data_all_proteins.csv format).

    python run_stage2_calibrated.py --input with_target_lipid/dpg3_causal_data_all_proteins.csv \
        --lipids CHOL DPSM DIPC --n-null 200 --out stage2_calibrated.csv

For every lipid type, prints and writes the block-bootstrap effect of binding for each
copy, the time-shift null p-values, and the classification of
stage2_contact_analysis/analysis/calibration.py. Use this classification in place of
the uncalibrated WAIC rule (delta WAIC > 2) of the Bayesian models.
"""
import os, sys, argparse
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'stage2_contact_analysis'))
from analysis.calibration import calibrated_classification

ap = argparse.ArgumentParser()
ap.add_argument('--input', required=True)
ap.add_argument('--lipids', nargs='+', default=None)
ap.add_argument('--n-null', type=int, default=200)
ap.add_argument('--n-blocks', type=int, default=None, help='default: chosen from the autocorrelation time')
ap.add_argument('--out', default='stage2_calibrated.csv')
a = ap.parse_args()

d = pd.read_csv(a.input).sort_values(['protein', 'frame'])
copies = sorted(d.protein.unique()); P = d.protein.map({c: i for i, c in enumerate(copies)}).values
T = d.target_lipid_bound.astype(int).values
lipids = a.lipids or [c[:-9] for c in d.columns if c.endswith('_contacts') and not c.startswith('DPG3')]
rows = []
for L in lipids:
    Y = d[f'{L}_contacts'].values.astype(float)
    r = calibrated_classification(Y, T, P, len(copies), n_null=a.n_null, n_blocks=a.n_blocks)
    bb, ns = r['block_bootstrap'], r['time_shift_null']
    print(f"\n{L}: {r['classification']}  (mean effect {bb['mean_effect']:+.2f}, 95% block CI "
          f"[{bb['mean_ci'][0]:+.2f}, {bb['mean_ci'][1]:+.2f}], p_effect {ns['p_effect']:.3f}, "
          f"dAIC {ns['observed']['d_aic']:.1f} vs null 95th pct {ns['d_aic_null_95']:.1f}, p_mixture {ns['p_mixture']:.3f})")
    for c, name in enumerate(copies):
        print(f"   {name}: {bb['effect'][c]:+.2f} [{bb['ci_low'][c]:+.2f}, {bb['ci_high'][c]:+.2f}]  "
              f"(block {bb['block_frames'][c]} frames, autocorrelation time {bb['tau_frames'][c]:.0f} frames)")
        if bb['block_frames'][c] < 2 * bb['tau_frames'][c]:
            print(f"   warning: block shorter than 2 autocorrelation times; the trajectory holds too few independent blocks for a reliable interval")
    rows.append(dict(lipid=L, classification=r['classification'], mean_effect=bb['mean_effect'],
                     ci_low=bb['mean_ci'][0], ci_high=bb['mean_ci'][1], p_effect=ns['p_effect'],
                     d_aic=ns['observed']['d_aic'], d_aic_null_95=ns['d_aic_null_95'], p_mixture=ns['p_mixture']))
pd.DataFrame(rows).to_csv(a.out, index=False)
print(f"\nwritten {a.out}")
