#!/usr/bin/env python3
"""
Run Stage 2 Hierarchical Bayesian Causal Analysis

This script performs hierarchical Bayesian analysis on Stage 1 causal data
to assess population-level effects of target lipid binding on lipid contacts.

Usage:
    python run_stage2_hierarchical.py

Author: Takeshi Sato, PhD (with Claude Code)
Kyoto Pharmaceutical University
2024
"""

import os
import sys

# Add stage2_contact_analysis to path
stage2_dir = os.path.join(os.path.dirname(__file__), 'stage2_contact_analysis')
sys.path.insert(0, stage2_dir)

from analysis.causal_analysis import load_causal_data
from analysis.hierarchical_causal_analysis import perform_hierarchical_causal_analysis


def main():
    """Main function to run hierarchical analysis"""

    print("=" * 80)
    print("STAGE 2: HIERARCHICAL BAYESIAN CAUSAL ANALYSIS")
    print("=" * 80)

    # Configuration
    BASE_DIR = os.path.dirname(__file__)
    STAGE1_OUTPUT = os.path.join(BASE_DIR, 'stage1_contact_analysis/output')
    STAGE2_OUTPUT = os.path.join(BASE_DIR, 'stage2_output')
    TARGET_LIPID_NAME = 'gm3'  # Change this to match your target lipid

    # Check if stage1 output exists
    if not os.path.exists(STAGE1_OUTPUT):
        print(f"Error: Stage 1 output directory not found: {STAGE1_OUTPUT}")
        print("Please run Stage 1 analysis first")
        return

    # Create output directory
    os.makedirs(STAGE2_OUTPUT, exist_ok=True)

    print(f"\nConfiguration:")
    print(f"  Stage 1 output: {STAGE1_OUTPUT}")
    print(f"  Stage 2 output: {STAGE2_OUTPUT}")
    print(f"  Target lipid: {TARGET_LIPID_NAME.upper()}")

    # Load causal data from Stage 1
    print(f"\n{'=' * 80}")
    print("LOADING CAUSAL DATA FROM STAGE 1")
    print(f"{'=' * 80}")

    causal_data = load_causal_data(STAGE1_OUTPUT, target_lipid_name=TARGET_LIPID_NAME)

    if causal_data is None or len(causal_data) == 0:
        print("Error: No causal data found. Cannot proceed with hierarchical analysis.")
        return

    print(f"\nLoaded data for {len(causal_data)} proteins:")
    for protein_name, df in causal_data.items():
        n_bound = df['target_lipid_bound'].sum()
        n_total = len(df)
        print(f"  {protein_name}: {n_total} frames, {n_bound} with {TARGET_LIPID_NAME.upper()} bound ({n_bound/n_total*100:.1f}%)")

    # Perform hierarchical Bayesian analysis
    print(f"\n{'=' * 80}")
    print("RUNNING HIERARCHICAL BAYESIAN ANALYSIS")
    print(f"{'=' * 80}")

    hierarchical_results = perform_hierarchical_causal_analysis(
        causal_data,
        output_dir=STAGE2_OUTPUT,
        target_lipid_name=TARGET_LIPID_NAME
    )

    if hierarchical_results is None or len(hierarchical_results) == 0:
        print("\nWarning: No hierarchical results generated")
        return

    # Print summary
    print(f"\n{'=' * 80}")
    print("HIERARCHICAL ANALYSIS COMPLETED")
    print(f"{'=' * 80}")

    print(f"\nAnalyzed {len(hierarchical_results)} lipid types:")
    for lipid_type, results in hierarchical_results.items():
        mu_beta = results['mu_beta_mean']
        ci = results['mu_beta_ci_94']
        rhat = results['mu_beta_rhat']
        ess = results['mu_beta_ess']

        # Convergence status
        converged = "✓" if rhat < 1.05 and ess > 400 else "✗"

        print(f"\n  {lipid_type}:")
        print(f"    Population effect (μβ): {mu_beta:.3f} [{ci[0]:.3f}, {ci[1]:.3f}]")
        print(f"    Convergence {converged}: r̂ = {rhat:.3f}, ESS = {ess:.0f}")

        if rhat > 1.1:
            print(f"    WARNING: Poor convergence (r̂ > 1.1)")
        elif rhat > 1.05:
            print(f"    CAUTION: Marginal convergence (r̂ > 1.05)")

        if ess < 400:
            print(f"    WARNING: Low effective sample size (ESS < 400)")

    print(f"\n{'=' * 80}")
    print(f"Results saved to: {STAGE2_OUTPUT}")
    print(f"  - Traces: {os.path.join(STAGE2_OUTPUT, 'hierarchical_analysis/traces')}")
    print(f"  - Plots: {os.path.join(STAGE2_OUTPUT, 'hierarchical_analysis')}")
    print(f"  - Summary: {os.path.join(STAGE2_OUTPUT, 'hierarchical_analysis/hierarchical_summary.csv')}")
    print(f"{'=' * 80}")

    print("\n✓ Stage 2 hierarchical analysis completed successfully!")


if __name__ == '__main__':
    main()
