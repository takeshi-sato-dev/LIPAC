"""
Save main results with MAIN_ prefix for easy identification
"""

import os
import pandas as pd


def save_main_causal_effects_summary(output_dir, protein_effect=None, lipid_effect=None,
                                    hierarchical_results=None):
    """
    Save main causal effects summary to a MAIN_ prefixed file

    Parameters
    ----------
    output_dir : str
        Output directory path
    protein_effect : dict
        Protein model effect results with 'beta', 'hdi', 'prob_negative'
    lipid_effect : dict
        Lipid model effect results with 'beta', 'hdi', 'prob_positive'
    hierarchical_results : dict
        Hierarchical model results with protein-specific effects
    """

    summary_path = os.path.join(output_dir, 'MAIN_causal_effects_summary.txt')

    with open(summary_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("                 MAIN CAUSAL EFFECTS SUMMARY\n")
        f.write("           GM3/DPG3 Effects on Molecular Interactions\n")
        f.write("=" * 80 + "\n\n")

        if protein_effect:
            f.write("PROTEIN-PROTEIN CONTACT EFFECTS:\n")
            f.write("-" * 40 + "\n")
            f.write(f"Causal Effect (β): {protein_effect.get('beta', 'N/A'):.4f}\n")

            if 'hdi' in protein_effect:
                f.write(f"95% HDI: [{protein_effect['hdi'][0]:.4f}, {protein_effect['hdi'][1]:.4f}]\n")

            if 'prob_negative' in protein_effect:
                f.write(f"P(β < 0) = {protein_effect['prob_negative']:.4f}\n")

            f.write("\nInterpretation: ")
            if protein_effect.get('beta', 0) < 0:
                f.write("Target lipid INHIBITS protein-protein contacts\n")
            else:
                f.write("Target lipid PROMOTES protein-protein contacts\n")
            f.write("\n")

        if lipid_effect:
            f.write("LIPID-PROTEIN CONTACT EFFECTS:\n")
            f.write("-" * 40 + "\n")
            f.write(f"Causal Effect (β): {lipid_effect.get('beta', 'N/A'):.4f}\n")

            if 'hdi' in lipid_effect:
                f.write(f"95% HDI: [{lipid_effect['hdi'][0]:.4f}, {lipid_effect['hdi'][1]:.4f}]\n")

            if 'prob_positive' in lipid_effect:
                f.write(f"P(β > 0) = {lipid_effect['prob_positive']:.4f}\n")

            f.write("\nInterpretation: ")
            if lipid_effect.get('beta', 0) > 0:
                f.write("Target lipid INCREASES lipid-protein contacts\n")
            else:
                f.write("Target lipid DECREASES lipid-protein contacts\n")
            f.write("\n")

        if hierarchical_results:
            f.write("PROTEIN-SPECIFIC EFFECTS:\n")
            f.write("-" * 40 + "\n")

            # Sort by effect size
            if 'protein_effects' in hierarchical_results:
                sorted_effects = sorted(hierarchical_results['protein_effects'].items(),
                                      key=lambda x: x[1].get('beta', 0))

                for protein, effect_data in sorted_effects:
                    beta = effect_data.get('beta', 'N/A')
                    prob = effect_data.get('prob_negative', 'N/A')

                    if isinstance(beta, (int, float)):
                        f.write(f"{protein}: β = {beta:.4f}, P(β < 0) = {prob:.4f}\n")
                    else:
                        f.write(f"{protein}: β = {beta}, P(β < 0) = {prob}\n")

            if 'global_effect' in hierarchical_results:
                f.write(f"\nGlobal Effect (μ_β): {hierarchical_results['global_effect']:.4f}\n")

        f.write("\n" + "=" * 80 + "\n")
        f.write("Note: These causal effects are from between-trajectory comparison\n")
        f.write("(with target lipid vs without target lipid systems)\n")
        f.write("=" * 80 + "\n")

    print(f"Main causal effects summary saved to: {summary_path}")

    # Also save as CSV for easier processing
    if protein_effect or lipid_effect:
        csv_data = []

        if protein_effect:
            csv_data.append({
                'Model': 'Protein Contact',
                'Beta': protein_effect.get('beta', None),
                'HDI_Lower': protein_effect.get('hdi', [None, None])[0],
                'HDI_Upper': protein_effect.get('hdi', [None, None])[1],
                'Probability': protein_effect.get('prob_negative', None),
                'Direction': 'Negative (Inhibitory)'
            })

        if lipid_effect:
            csv_data.append({
                'Model': 'Lipid Contact',
                'Beta': lipid_effect.get('beta', None),
                'HDI_Lower': lipid_effect.get('hdi', [None, None])[0],
                'HDI_Upper': lipid_effect.get('hdi', [None, None])[1],
                'Probability': lipid_effect.get('prob_positive', None),
                'Direction': 'Positive (Promoting)'
            })

        if csv_data:
            csv_path = os.path.join(output_dir, 'MAIN_causal_effects.csv')
            pd.DataFrame(csv_data).to_csv(csv_path, index=False)
            print(f"Main causal effects CSV saved to: {csv_path}")