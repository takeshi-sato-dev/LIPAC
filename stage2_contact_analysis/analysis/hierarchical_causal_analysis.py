"""
Hierarchical Bayesian Causal Analysis for Stage 2 Lipid Contacts

This module implements hierarchical Bayesian models to assess population-level
causal effects across multiple protein copies.

Author: Takeshi Sato, PhD (with Claude Code)
Kyoto Pharmaceutical University
2024
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pymc as pm
import arviz as az
import warnings
warnings.filterwarnings('ignore')


def perform_hierarchical_causal_analysis(causal_data, output_dir, target_lipid_name='TARGET_LIPID'):
    """Perform hierarchical Bayesian analysis across all proteins for each lipid type

    Parameters
    ----------
    causal_data : dict
        Dictionary with protein names as keys and DataFrames as values
    output_dir : str
        Output directory for hierarchical model results
    target_lipid_name : str
        Name of the target lipid

    Returns
    -------
    dict
        Hierarchical analysis results for each lipid type
    """
    print("\n===== Hierarchical Bayesian Causal Analysis =====")
    print(f"Analyzing population-level effects of {target_lipid_name} binding")

    # Create output directories
    hierarchical_dir = os.path.join(output_dir, 'hierarchical_analysis')
    os.makedirs(hierarchical_dir, exist_ok=True)

    traces_dir = os.path.join(hierarchical_dir, 'traces')
    os.makedirs(traces_dir, exist_ok=True)

    # Get all lipid types from the first protein
    first_protein_df = list(causal_data.values())[0]
    target_lipid_col = f'{target_lipid_name}_contacts'

    residue_contact_cols = [col for col in first_protein_df.columns
                           if col.endswith('_contacts') and col != target_lipid_col]
    unique_molecule_cols = [col for col in first_protein_df.columns
                           if col.endswith('_unique_molecules')]

    lipid_cols = residue_contact_cols + unique_molecule_cols

    if not lipid_cols:
        print("Error: No lipid contact columns found!")
        return None

    print(f"Found {len(lipid_cols)} lipid types to analyze:")
    for col in lipid_cols:
        print(f"  - {col}")

    hierarchical_results = {}

    # Analyze each lipid type with hierarchical model
    for lipid_col in lipid_cols:
        if lipid_col.endswith('_contacts'):
            lipid_name = lipid_col.replace('_contacts', '')
            analysis_type = 'contacts'
        elif lipid_col.endswith('_unique_molecules'):
            lipid_name = lipid_col.replace('_unique_molecules', '')
            analysis_type = 'unique'
        else:
            continue

        print(f"\n--- Hierarchical Analysis: {lipid_name} {analysis_type} ---")

        # Prepare data from all proteins
        all_target_bound = []
        all_lipid_contacts = []
        protein_indices = []
        protein_names = []

        for protein_idx, (protein_name, df) in enumerate(causal_data.items()):
            target_lipid_bound = df['target_lipid_bound'].astype(int).values
            lipid_contacts = df[lipid_col].values

            # Remove NaN values
            valid_mask = ~np.isnan(lipid_contacts)
            target_lipid_bound = target_lipid_bound[valid_mask]
            lipid_contacts = lipid_contacts[valid_mask]

            if len(lipid_contacts) < 20:
                print(f"  Skipping {protein_name}: insufficient data (n={len(lipid_contacts)})")
                continue

            all_target_bound.extend(target_lipid_bound)
            all_lipid_contacts.extend(lipid_contacts)
            protein_indices.extend([protein_idx] * len(lipid_contacts))
            protein_names.append(protein_name)

            print(f"  {protein_name}: {len(lipid_contacts)} frames, "
                  f"{np.sum(target_lipid_bound)} bound ({np.mean(target_lipid_bound)*100:.1f}%)")

        if len(protein_names) < 2:
            print(f"  Insufficient proteins for hierarchical analysis (need ≥2, have {len(protein_names)})")
            continue

        # Convert to numpy arrays
        all_target_bound = np.array(all_target_bound)
        all_lipid_contacts = np.array(all_lipid_contacts)
        protein_indices = np.array(protein_indices)
        n_proteins = len(protein_names)

        # Check for zero variance
        std_lipid = np.std(all_lipid_contacts)
        if std_lipid < 1e-6:
            print(f"  Skipping {lipid_name}: zero variance")
            continue

        print(f"  Building hierarchical model for {n_proteins} proteins...")
        print(f"  Total data points: {len(all_lipid_contacts)}")

        # Build hierarchical Bayesian model
        try:
            with pm.Model() as hierarchical_model:
                # Hyperpriors for population-level parameters
                mu_alpha = pm.Normal('mu_alpha', mu=np.mean(all_lipid_contacts), sigma=std_lipid*2)
                sigma_alpha = pm.HalfNormal('sigma_alpha', sigma=std_lipid)

                mu_beta = pm.Normal('mu_beta', mu=0, sigma=std_lipid*2)
                sigma_beta = pm.HalfNormal('sigma_beta', sigma=std_lipid)

                # Protein-specific parameters (hierarchical)
                alpha = pm.Normal('alpha', mu=mu_alpha, sigma=sigma_alpha, shape=n_proteins)
                beta = pm.Normal('beta', mu=mu_beta, sigma=sigma_beta, shape=n_proteins)

                # Residual error
                sigma = pm.HalfNormal('sigma', sigma=std_lipid)

                # Linear model for each data point
                mu = alpha[protein_indices] + beta[protein_indices] * all_target_bound

                # Likelihood
                y = pm.Normal('y', mu=mu, sigma=sigma, observed=all_lipid_contacts)

                # Sample
                print(f"  Sampling posterior distribution...")
                trace = pm.sample(2000, tune=1000, chains=4, random_seed=42,
                                progressbar=True, return_inferencedata=True)

            print(f"  ✓ Sampling completed")

            # Extract results
            mu_beta_samples = trace.posterior['mu_beta'].values.flatten()
            sigma_beta_samples = trace.posterior['sigma_beta'].values.flatten()

            # Individual protein effects
            beta_samples = trace.posterior['beta'].values  # shape: (chains, draws, n_proteins)

            # Calculate statistics for population-level effect
            mu_beta_mean = np.mean(mu_beta_samples)
            mu_beta_std = np.std(mu_beta_samples)
            mu_beta_ci = np.percentile(mu_beta_samples, [3, 97])  # 94% HDI
            mu_beta_prob_negative = (mu_beta_samples < 0).mean()
            mu_beta_prob_positive = (mu_beta_samples > 0).mean()

            # Convergence diagnostics
            summary = az.summary(trace, var_names=['mu_beta', 'sigma_beta', 'alpha', 'beta'])
            mu_beta_rhat = summary.loc['mu_beta', 'r_hat']
            mu_beta_ess = summary.loc['mu_beta', 'ess_bulk']

            print(f"\n  Population-level effect (μβ):")
            print(f"    Mean: {mu_beta_mean:.3f} ± {mu_beta_std:.3f}")
            print(f"    94% HDI: [{mu_beta_ci[0]:.3f}, {mu_beta_ci[1]:.3f}]")
            print(f"    P(μβ < 0): {mu_beta_prob_negative:.3f}")
            print(f"    r̂: {mu_beta_rhat:.3f}")
            print(f"    ESS: {mu_beta_ess:.0f}")

            # Individual protein effects
            individual_effects = {}
            for i, protein_name in enumerate(protein_names):
                beta_i_samples = beta_samples[:, :, i].flatten()
                beta_i_mean = np.mean(beta_i_samples)
                beta_i_ci = np.percentile(beta_i_samples, [3, 97])
                beta_i_rhat = summary.loc[f'beta[{i}]', 'r_hat']
                beta_i_ess = summary.loc[f'beta[{i}]', 'ess_bulk']

                individual_effects[protein_name] = {
                    'mean': beta_i_mean,
                    'ci_94': beta_i_ci,
                    'rhat': beta_i_rhat,
                    'ess': beta_i_ess
                }

                print(f"    {protein_name}: β = {beta_i_mean:.3f} [{beta_i_ci[0]:.3f}, {beta_i_ci[1]:.3f}], "
                      f"r̂ = {beta_i_rhat:.3f}, ESS = {beta_i_ess:.0f}")

            # Store results
            hierarchical_results[f'{lipid_name}_{analysis_type}'] = {
                'mu_beta_mean': mu_beta_mean,
                'mu_beta_std': mu_beta_std,
                'mu_beta_ci_94': mu_beta_ci,
                'mu_beta_prob_negative': mu_beta_prob_negative,
                'mu_beta_prob_positive': mu_beta_prob_positive,
                'mu_beta_rhat': mu_beta_rhat,
                'mu_beta_ess': mu_beta_ess,
                'sigma_beta_mean': np.mean(sigma_beta_samples),
                'sigma_beta_ci_94': np.percentile(sigma_beta_samples, [3, 97]),
                'individual_effects': individual_effects,
                'trace': trace,
                'protein_names': protein_names,
                'n_proteins': n_proteins,
                'convergence_summary': summary
            }

            # Save trace
            trace_path = os.path.join(traces_dir, f'{lipid_name}_{analysis_type}_hierarchical_trace.nc')
            trace.to_netcdf(trace_path)
            print(f"  ✓ Trace saved: {trace_path}")

            # Generate diagnostic plots
            _generate_hierarchical_plots(trace, lipid_name, analysis_type, protein_names,
                                        hierarchical_dir, target_lipid_name)

        except Exception as e:
            print(f"  Error in hierarchical analysis for {lipid_name}: {e}")
            import traceback
            traceback.print_exc()
            continue

    # Save summary
    _save_hierarchical_summary(hierarchical_results, hierarchical_dir, target_lipid_name)

    return hierarchical_results


def _generate_hierarchical_plots(trace, lipid_name, analysis_type, protein_names,
                                 output_dir, target_lipid_name):
    """Generate diagnostic plots for hierarchical model"""

    full_name = f'{lipid_name}_{analysis_type}'

    # 1. Trace plot for mu_beta
    try:
        fig, axes = plt.subplots(2, 1, figsize=(12, 8))

        # Posterior distribution
        mu_beta_samples = trace.posterior['mu_beta'].values.flatten()
        axes[0].hist(mu_beta_samples, bins=50, density=True, alpha=0.7,
                    color='steelblue', edgecolor='black')
        axes[0].axvline(mu_beta_samples.mean(), color='red', linestyle='--',
                       linewidth=2, label=f'Mean = {mu_beta_samples.mean():.3f}')
        axes[0].set_xlabel('μβ (population-level effect)')
        axes[0].set_ylabel('Density')
        axes[0].set_title(f'{target_lipid_name} → {full_name}: Posterior Distribution')
        axes[0].legend()

        # Add convergence diagnostics
        summary = az.summary(trace, var_names=['mu_beta'])
        rhat = summary.loc['mu_beta', 'r_hat']
        ess = summary.loc['mu_beta', 'ess_bulk']
        axes[0].text(0.05, 0.95, f'r̂ = {rhat:.2f}\nESS = {ess:.0f}',
                    transform=axes[0].transAxes, verticalalignment='top',
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

        # Trace plot
        for chain in range(trace.posterior['mu_beta'].shape[0]):
            chain_samples = trace.posterior['mu_beta'].values[chain, :]
            axes[1].plot(chain_samples, alpha=0.7, linewidth=1)
        axes[1].set_xlabel('Iteration')
        axes[1].set_ylabel('μβ')
        axes[1].set_title('MCMC Trace')
        axes[1].axhline(mu_beta_samples.mean(), color='red', linestyle='--',
                       linewidth=1, alpha=0.5)

        plt.tight_layout()
        plot_path_png = os.path.join(output_dir, f'{full_name}_mu_beta_trace.png')
        plot_path_svg = os.path.join(output_dir, f'{full_name}_mu_beta_trace.svg')
        plt.savefig(plot_path_png, dpi=300, bbox_inches='tight')
        plt.savefig(plot_path_svg, format='svg', bbox_inches='tight')
        plt.close()
        print(f"  ✓ mu_beta trace plot saved: {plot_path_png}")
        print(f"  ✓ mu_beta trace plot saved: {plot_path_svg}")
    except Exception as e:
        print(f"  Warning: Could not generate mu_beta trace plot: {e}")

    # 2. Forest plot for individual effects
    try:
        fig, ax = plt.subplots(figsize=(10, 6))
        az.plot_forest(trace, var_names=['beta'], combined=False, ax=ax)
        ax.set_title(f'{target_lipid_name} → {full_name}: Individual Protein Effects')
        ax.set_yticklabels(protein_names)
        plt.tight_layout()

        forest_path_png = os.path.join(output_dir, f'{full_name}_forest.png')
        forest_path_svg = os.path.join(output_dir, f'{full_name}_forest.svg')
        plt.savefig(forest_path_png, dpi=300, bbox_inches='tight')
        plt.savefig(forest_path_svg, format='svg', bbox_inches='tight')
        plt.close()
        print(f"  ✓ Forest plot saved: {forest_path_png}")
        print(f"  ✓ Forest plot saved: {forest_path_svg}")
    except Exception as e:
        print(f"  Warning: Could not generate forest plot: {e}")

    # 3. Full trace plot (all parameters)
    try:
        az.plot_trace(trace, var_names=['mu_beta', 'sigma_beta', 'alpha', 'beta'])
        plt.suptitle(f'{target_lipid_name} → {full_name}: Full Hierarchical Model', y=1.02)
        plt.tight_layout()

        full_trace_path_png = os.path.join(output_dir, f'{full_name}_full_trace.png')
        full_trace_path_svg = os.path.join(output_dir, f'{full_name}_full_trace.svg')
        plt.savefig(full_trace_path_png, dpi=300, bbox_inches='tight')
        plt.savefig(full_trace_path_svg, format='svg', bbox_inches='tight')
        plt.close()
        print(f"  ✓ Full trace plot saved: {full_trace_path_png}")
        print(f"  ✓ Full trace plot saved: {full_trace_path_svg}")
    except Exception as e:
        print(f"  Warning: Could not generate full trace plot: {e}")


def _save_hierarchical_summary(results, output_dir, target_lipid_name):
    """Save hierarchical analysis summary"""

    if not results:
        print("No hierarchical results to save")
        return

    # Save as CSV
    summary_data = []
    for lipid_type, res in results.items():
        row = {
            'lipid_type': lipid_type,
            'mu_beta_mean': res['mu_beta_mean'],
            'mu_beta_std': res['mu_beta_std'],
            'mu_beta_ci_lower': res['mu_beta_ci_94'][0],
            'mu_beta_ci_upper': res['mu_beta_ci_94'][1],
            'mu_beta_prob_negative': res['mu_beta_prob_negative'],
            'mu_beta_prob_positive': res['mu_beta_prob_positive'],
            'mu_beta_rhat': res['mu_beta_rhat'],
            'mu_beta_ess': res['mu_beta_ess'],
            'sigma_beta_mean': res['sigma_beta_mean'],
            'n_proteins': res['n_proteins']
        }
        summary_data.append(row)

        # Add individual effects
        for protein_name, eff in res['individual_effects'].items():
            row_ind = {
                'lipid_type': lipid_type,
                'protein': protein_name,
                'beta_mean': eff['mean'],
                'beta_ci_lower': eff['ci_94'][0],
                'beta_ci_upper': eff['ci_94'][1],
                'beta_rhat': eff['rhat'],
                'beta_ess': eff['ess']
            }
            summary_data.append(row_ind)

    summary_df = pd.DataFrame(summary_data)
    summary_path = os.path.join(output_dir, 'hierarchical_summary.csv')
    summary_df.to_csv(summary_path, index=False)
    print(f"\n✓ Hierarchical summary saved: {summary_path}")

    # Save as text
    text_path = os.path.join(output_dir, 'hierarchical_summary.txt')
    with open(text_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write(f"Hierarchical Bayesian Analysis Summary - {target_lipid_name}\n")
        f.write("=" * 80 + "\n\n")

        for lipid_type, res in results.items():
            f.write(f"\n{'=' * 80}\n")
            f.write(f"{lipid_type.upper()}\n")
            f.write(f"{'=' * 80}\n\n")

            f.write("Population-level effect (μβ):\n")
            f.write(f"  Mean: {res['mu_beta_mean']:.4f} ± {res['mu_beta_std']:.4f}\n")
            f.write(f"  94% HDI: [{res['mu_beta_ci_94'][0]:.4f}, {res['mu_beta_ci_94'][1]:.4f}]\n")
            f.write(f"  P(μβ < 0): {res['mu_beta_prob_negative']:.4f}\n")
            f.write(f"  P(μβ > 0): {res['mu_beta_prob_positive']:.4f}\n")
            f.write(f"  r̂: {res['mu_beta_rhat']:.4f}\n")
            f.write(f"  ESS: {res['mu_beta_ess']:.0f}\n\n")

            f.write("Individual protein effects:\n")
            for protein_name, eff in res['individual_effects'].items():
                f.write(f"  {protein_name}:\n")
                f.write(f"    β = {eff['mean']:.4f} [{eff['ci_94'][0]:.4f}, {eff['ci_94'][1]:.4f}]\n")
                f.write(f"    r̂ = {eff['rhat']:.4f}, ESS = {eff['ess']:.0f}\n")

    print(f"✓ Hierarchical text summary saved: {text_path}")
