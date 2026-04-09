"""
Example: Mixture Causal Model Analysis

This script demonstrates how to apply the Bayesian mixture causal model
to Stage 1 output data for detecting cooperative lipid dynamics.

Author: Takeshi Sato
Kyoto Pharmaceutical University
"""

import sys
import os

# Add parent directory to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from stage2_contact_analysis.analysis.mixture_causal_analysis import (
    load_stage1_data,
    analyze_lipid_type
)


def main():
    """Run mixture causal analysis on example data"""

    # Path to Stage 1 output CSV (replace with your actual path)
    input_csv = "path/to/stage1_output.csv"

    # Output directory for results
    output_dir = "mixture_results"

    # Load Stage 1 data
    # Expected columns: frame, protein_copy, binding_state, lipid_type, contact_count
    print("Loading data...")
    df = load_stage1_data(input_csv)

    print(f"Loaded {len(df)} observations")
    print(f"Lipid types: {', '.join(df['lipid_type'].unique())}")

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Analyze each lipid type
    results = []
    for lipid_name in df['lipid_type'].unique():
        df_lipid = df[df['lipid_type'] == lipid_name].copy()

        # Skip if insufficient data
        if len(df_lipid) < 100:
            print(f"\nSkipping {lipid_name}: insufficient data (N={len(df_lipid)})")
            continue

        # Run mixture analysis
        try:
            result = analyze_lipid_type(df_lipid, lipid_name, output_dir)
            results.append(result)

            # Print summary
            print(f"\n{lipid_name}: {result['classification'].upper()}")
            print(f"  ΔWAIC = {result['delta_waic']:.1f}")
            print(f"  β_coop = {result['beta_coop_mean']:.2f} "
                  f"[{result['beta_coop_hdi'][0]:.2f}, {result['beta_coop_hdi'][1]:.2f}]")
            print(f"  π = {result['pi_mean']:.2f} "
                  f"[{result['pi_hdi'][0]:.2f}, {result['pi_hdi'][1]:.2f}]")

        except Exception as e:
            print(f"\nError analyzing {lipid_name}: {e}")

    # Save summary
    if results:
        import pandas as pd
        results_df = pd.DataFrame(results)
        results_df.to_csv(os.path.join(output_dir, 'summary.csv'), index=False)
        print(f"\nResults saved to {output_dir}/summary.csv")


if __name__ == '__main__':
    main()
