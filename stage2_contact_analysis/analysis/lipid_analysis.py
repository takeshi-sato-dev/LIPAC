"""
Lipid-specific analysis functions for BayesianLipidAnalysis

Author: Takeshi Sato, PhD
Kyoto Pharmaceutical University
2024
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import pymc as pm
import arviz as az

# Use absolute import instead of relative import
from visualization.plots import plot_convergence_diagnostics, plot_lipid_comparison


def analyze_lipid_specific_effects(combined_df, matched_df, output_dir):
    """Analyze the effect of target lipid on interactions with specific lipid types
    
    Parameters
    ----------
    combined_df : pd.DataFrame
        Combined dataframe with target lipid and without target lipid data
    matched_df : pd.DataFrame
        Matched dataframe with calculated differences
    output_dir : str
        Output directory
    
    Returns
    -------
    bool
        Success status
    """
    print("\n===== Lipid-Specific Target Lipid Effect Analysis =====")
    
    # Create output directory
    lipid_dir = os.path.join(output_dir, "lipid_specific_analysis")
    os.makedirs(lipid_dir, exist_ok=True)
    print(f"Created directory: {lipid_dir}")
    
    # Check for empty dataframes
    if combined_df is None or matched_df is None:
        print("Error: Input dataframes are empty or None.")
        return False
    
    # Split combined data into with_lipid and without_lipid
    with_lipid_df = combined_df[combined_df['has_target_lipid'] == True].copy()
    without_lipid_df = combined_df[combined_df['has_target_lipid'] == False].copy()
    
    print(f"With target lipid data: {len(with_lipid_df)} rows")
    print(f"Without target lipid data: {len(without_lipid_df)} rows")
    
    # Identify lipid contact columns - hard-coded list to ensure they're found
    lipid_columns = ['CHOL_contact', 'DOPS_contact', 'DIPC_contact', 'DPSM_contact']
    available_lipids = [col for col in lipid_columns if col in with_lipid_df.columns and col in without_lipid_df.columns]
    
    if not available_lipids:
        print("Error: No lipid contact columns found in both datasets.")
        # Debug information
        print(f"Columns in with_lipid_df: {with_lipid_df.columns.tolist()}")
        print(f"Columns in without_lipid_df: {without_lipid_df.columns.tolist()}")
        return False
    
    print(f"Available lipid contact columns: {available_lipids}")
    
    # Create enhanced matched dataframe with lipid-specific differences
    enhanced_matched = matched_df.copy()
    
    # Add lipid-specific differences
    print("Creating lipid-specific difference columns")
    for protein in matched_df['protein'].unique():
        # Get residues for this protein
        protein_residues = matched_df[matched_df['protein'] == protein]['residue'].unique()
        
        # Extract with_lipid and without_lipid data for this protein
        with_lipid_protein = with_lipid_df[with_lipid_df['protein'] == protein]
        without_lipid_protein = without_lipid_df[without_lipid_df['protein'] == protein]
        
        # For each residue, calculate lipid-specific differences
        for residue in protein_residues:
            # Get rows for this residue
            with_lipid_row = with_lipid_protein[with_lipid_protein['residue'] == residue]
            without_lipid_row = without_lipid_protein[without_lipid_protein['residue'] == residue]
            
            if len(with_lipid_row) == 0 or len(without_lipid_row) == 0:
                continue
            
            # For each lipid type, calculate the difference
            for lipid_col in available_lipids:
                with_lipid_val = with_lipid_row[lipid_col].values[0]
                without_lipid_val = without_lipid_row[lipid_col].values[0]
                diff = with_lipid_val - without_lipid_val
                
                # Add to enhanced matched dataframe
                mask = (enhanced_matched['protein'] == protein) & (enhanced_matched['residue'] == residue)
                enhanced_matched.loc[mask, f"{lipid_col}_with_lipid"] = with_lipid_val
                enhanced_matched.loc[mask, f"{lipid_col}_without_lipid"] = without_lipid_val
                enhanced_matched.loc[mask, f"{lipid_col}_diff"] = diff
    
    # Verify target lipid contact column exists
    if 'target_lipid_contact' not in enhanced_matched.columns:
        # Try to find any target lipid column
        if 'CHOL_contact' in enhanced_matched.columns:
            enhanced_matched['target_lipid_contact'] = enhanced_matched['CHOL_contact']
        elif 'DPG3_contact' in enhanced_matched.columns:
            enhanced_matched['target_lipid_contact'] = enhanced_matched['DPG3_contact']
        else:
            enhanced_matched['target_lipid_contact'] = 0.0
    
    # Ensure target_lipid_contact column exists
    if 'target_lipid_contact' not in enhanced_matched.columns:
        print("Error: Target lipid contact column not found in matched data.")
        return False
    
    # Analysis results storage
    lipid_results = []
    
    # Analyze each lipid type using Bayesian regression
    for lipid_col in available_lipids:
        lipid_type = lipid_col.replace('_contact', '')
        print(f"\nAnalyzing {lipid_type} interactions")
        
        # Check if difference column exists
        diff_col = f"{lipid_col}_diff"
        if diff_col not in enhanced_matched.columns:
            print(f"Warning: {diff_col} column not found, skipping {lipid_type}")
            continue
        
        # Get X and Y data
        X_target_lipid = enhanced_matched['target_lipid_contact'].values
        Y_lipid_diff = enhanced_matched[diff_col].values
        
        # Print data stats
        print(f"  Data points: {len(X_target_lipid)}")
        print(f"  Target lipid contact range: {np.min(X_target_lipid):.3f} to {np.max(X_target_lipid):.3f}")
        print(f"  {lipid_type} diff range: {np.min(Y_lipid_diff):.3f} to {np.max(Y_lipid_diff):.3f}")
        
        # Correlation analysis
        corr, p_val = stats.pearsonr(X_target_lipid, Y_lipid_diff)
        print(f"  Correlation: r={corr:.3f}, p={p_val:.3e}")
        
        # Bayesian regression model
        with pm.Model() as lipid_model:
            # Prior distributions
            alpha = pm.Normal('alpha', mu=0, sigma=1)  # Intercept
            beta = pm.Normal('beta', mu=0, sigma=1)    # Target lipid contact coefficient
            sigma = pm.HalfNormal('sigma', sigma=1)    # Error standard deviation
            
            # Linear model
            mu = alpha + beta * X_target_lipid
            
            # Likelihood
            Y_obs = pm.Normal('Y_obs', mu=mu, sigma=sigma, observed=Y_lipid_diff)
            
            # Sampling
            try:
                import sys
                sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
                import config
                
                lipid_trace = pm.sample(
                    config.MCMC_SAMPLES,
                    tune=config.TUNE_SAMPLES,
                    chains=config.CHAINS,
                    random_seed=config.RANDOM_SEED,
                    target_accept=config.TARGET_ACCEPT
                )
                
                # Save trace netcdf
                trace_path = os.path.join(lipid_dir, f"{lipid_type}_trace.nc")
                az.to_netcdf(lipid_trace, trace_path)
                
                # Get summary statistics
                summary = az.summary(lipid_trace)
                summary.to_csv(os.path.join(lipid_dir, f"{lipid_type}_summary.csv"))
                
                # Add convergence diagnostics
                plot_convergence_diagnostics(lipid_trace, f'{lipid_type}_model', lipid_dir)
                
                # Extract posterior samples for beta
                beta_samples = lipid_trace.posterior['beta'].values.flatten()
                beta_mean = beta_samples.mean()
                beta_hdi = az.hdi(beta_samples, hdi_prob=0.95)
                
                # Determine effect direction and probability
                if beta_mean < 0:
                    prob_effect = (beta_samples < 0).mean()
                    effect_direction = "negative"
                else:
                    prob_effect = (beta_samples > 0).mean()
                    effect_direction = "positive"
                
                print(f"  Effect size (β): {beta_mean:.3f}")
                print(f"  95% HDI: [{beta_hdi[0]:.3f}, {beta_hdi[1]:.3f}]")
                print(f"  P({effect_direction[0]}): {prob_effect:.3f}")
                
                # Store results
                lipid_results.append({
                    'lipid_type': lipid_type,
                    'mean_effect': beta_mean,
                    'hdi_lower': beta_hdi[0],
                    'hdi_upper': beta_hdi[1],
                    'prob_effect': prob_effect,
                    'effect_direction': effect_direction
                })
                
                # Create posterior distribution plot
                plt.figure(figsize=(10, 6))
                sns.histplot(beta_samples, kde=True, stat="density")
                plt.axvline(x=0, color='r', linestyle='--')
                plt.axvline(x=beta_hdi[0], color='k', linestyle='--', alpha=0.5)
                plt.axvline(x=beta_hdi[1], color='k', linestyle='--', alpha=0.5)
                
                # Statistical information
                plt.text(0.05, 0.95, f'95% HDI: [{beta_hdi[0]:.3f}, {beta_hdi[1]:.3f}]', 
                        transform=plt.gca().transAxes, bbox=dict(facecolor='white', alpha=0.7))
                plt.text(0.05, 0.85, f'P({effect_direction[0]}): {prob_effect:.3f}', 
                        transform=plt.gca().transAxes, bbox=dict(facecolor='white', alpha=0.7))
                
                plt.title(f'Posterior Distribution of Target Lipid Effect on {lipid_type} Contacts')
                plt.xlabel('Effect Size (β)')
                plt.ylabel('Density')
                plt.grid(True, alpha=0.3)
                plt.savefig(os.path.join(lipid_dir, f'{lipid_type}_effect.png'), dpi=300)
                plt.savefig(os.path.join(lipid_dir, f'{lipid_type}_effect.ps'))
                plt.close()
                
                # Create regression scatter plot
                plt.figure(figsize=(10, 6))
                plt.scatter(X_target_lipid, Y_lipid_diff, alpha=0.7)
                
                # Add regression line
                x_range = np.linspace(0, np.max(X_target_lipid), 100)
                alpha_mean = lipid_trace.posterior['alpha'].values.flatten().mean()
                plt.plot(x_range, alpha_mean + beta_mean * x_range, 'r-', linewidth=2, 
                        label=f'β = {beta_mean:.3f}')
                
                # Reference line
                plt.axhline(y=0, color='k', linestyle='--', alpha=0.5)
                
                # Add residue labels
                for i, row in enhanced_matched.iterrows():
                    if row['target_lipid_contact'] > np.percentile(X_target_lipid, 75) or abs(row[diff_col]) > np.percentile(abs(Y_lipid_diff), 75):
                        plt.annotate(
                            f"{row['residue']}", 
                            (row['target_lipid_contact'], row[diff_col]),
                            xytext=(5, 5), textcoords='offset points',
                            fontsize=8
                        )
                
                plt.title(f'Target Lipid Contact vs Change in {lipid_type} Contact')
                plt.xlabel('Target Lipid Contact Frequency')
                plt.ylabel(f'Change in {lipid_type} Contact')
                plt.grid(True, alpha=0.3)
                plt.legend()
                plt.savefig(os.path.join(lipid_dir, f'{lipid_type}_regression.png'), dpi=300)
                plt.savefig(os.path.join(lipid_dir, f'{lipid_type}_regression.ps'))
                plt.close()
                
            except Exception as e:
                print(f"  Error in Bayesian analysis: {str(e)}")
                import traceback
                traceback.print_exc()
    
    # Save consolidated results if any successful analyses
    if lipid_results:
        # Convert to DataFrame
        results_df = pd.DataFrame(lipid_results)
        
        # Sort by effect size magnitude
        results_df['abs_effect'] = results_df['mean_effect'].abs()
        results_df = results_df.sort_values('abs_effect', ascending=False)
        results_df = results_df.drop('abs_effect', axis=1)
        
        # Save results
        results_df.to_csv(os.path.join(lipid_dir, 'lipid_specific_effects.csv'), index=False)
        
        # Create comparison plot
        plot_lipid_comparison(results_df, lipid_dir)
        
        # Generate summary report
        _generate_lipid_summary_report(results_df, lipid_dir)
        
        print(f"\nLipid-specific analysis complete. Results saved to {lipid_dir}")
        return True
    else:
        print("\nNo successful lipid-specific analyses completed.")
        return False


def _generate_lipid_summary_report(results_df, output_dir):
    """Generate summary report for lipid-specific analysis
    
    Parameters
    ----------
    results_df : pd.DataFrame
        Results dataframe
    output_dir : str
        Output directory
    """
    with open(os.path.join(output_dir, 'lipid_specific_analysis_summary.txt'), 'w') as f:
        f.write("===== Target Lipid Effects on Specific Lipid Types =====\n\n")
        
        for _, row in results_df.iterrows():
            f.write(f"Lipid Type: {row['lipid_type']}\n")
            f.write(f"  Effect Size (β): {row['mean_effect']:.3f}\n")
            f.write(f"  95% HDI: [{row['hdi_lower']:.3f}, {row['hdi_upper']:.3f}]\n")
            
            if row['effect_direction'] == 'positive':
                f.write(f"  Probability of positive effect: {row['prob_effect']:.3f}\n")
            else:
                f.write(f"  Probability of negative effect: {row['prob_effect']:.3f}\n")
            
            # Add interpretation
            f.write("  Interpretation: ")
            if row['prob_effect'] > 0.95:
                if row['effect_direction'] == 'positive':
                    f.write("Strong evidence that target lipid increases contact with this lipid type.\n")
                else:
                    f.write("Strong evidence that target lipid decreases contact with this lipid type.\n")
            elif row['prob_effect'] > 0.75:
                if row['effect_direction'] == 'positive':
                    f.write("Moderate evidence that target lipid increases contact with this lipid type.\n")
                else:
                    f.write("Moderate evidence that target lipid decreases contact with this lipid type.\n")
            else:
                f.write("Weak evidence for target lipid effect on this lipid type.\n")
            f.write("\n")
        
        # Overall conclusions
        f.write("\n===== Overall Conclusions =====\n")
        
        # Count effect directions
        strong_positive = sum(1 for _, r in results_df.iterrows() 
                            if r['effect_direction'] == 'positive' and r['prob_effect'] > 0.95)
        strong_negative = sum(1 for _, r in results_df.iterrows() 
                            if r['effect_direction'] == 'negative' and r['prob_effect'] > 0.95)
        
        f.write(f"Strong positive effects: {strong_positive} lipid types\n")
        f.write(f"Strong negative effects: {strong_negative} lipid types\n")
        
        if strong_positive > 0 and strong_negative > 0:
            f.write("\nTarget lipid shows differential effects on different lipid types, suggesting ")
            f.write("specific lipid-lipid interactions that modulate protein-membrane contacts.\n")
        elif strong_positive > strong_negative:
            f.write("\nTarget lipid primarily increases contacts with other lipid types, suggesting ")
            f.write("cooperative lipid-lipid interactions.\n")
        elif strong_negative > strong_positive:
            f.write("\nTarget lipid primarily decreases contacts with other lipid types, suggesting ")
            f.write("competitive or repulsive interactions.\n")
        else:
            f.write("\nTarget lipid shows weak or mixed effects on other lipid types.\n")


def analyze_per_peptide_effects(combined_df, matched_df, output_dir):
    """Analyze lipid-specific effects for each peptide separately
    
    Parameters
    ----------
    combined_df : pd.DataFrame
        Combined dataframe with target lipid and without target lipid data
    matched_df : pd.DataFrame
        Matched dataframe with calculated differences
    output_dir : str
        Output directory for analysis results
    
    Returns
    -------
    bool
        Success status
    """
    print("\n===== Per-Peptide Lipid-Specific Analysis =====")
    
    # Create per-peptide output directory
    per_peptide_dir = os.path.join(output_dir, "per_peptide_analysis")
    os.makedirs(per_peptide_dir, exist_ok=True)
    
    # Get unique proteins from the data
    if 'protein' not in combined_df.columns:
        print("Error: 'protein' column not found in data. Cannot perform per-peptide analysis.")
        return False
    
    unique_proteins = combined_df['protein'].unique()
    print(f"Found {len(unique_proteins)} unique proteins: {list(unique_proteins)}")
    
    successful_analyses = 0
    
    for protein in unique_proteins:
        print(f"\nAnalyzing protein: {protein}")
        
        # Filter data for this protein
        protein_combined_df = combined_df[combined_df['protein'] == protein].copy()
        protein_matched_df = matched_df[matched_df['protein'] == protein].copy()
        
        print(f"  Data points - Combined: {len(protein_combined_df)}, Matched: {len(protein_matched_df)}")
        
        if len(protein_matched_df) < 10:
            print(f"  Warning: Too few data points ({len(protein_matched_df)}) for protein {protein}. Skipping.")
            continue
        
        # Create protein-specific output directory
        protein_output_dir = os.path.join(per_peptide_dir, f"protein_{protein}")
        os.makedirs(protein_output_dir, exist_ok=True)
        
        try:
            # Run the same lipid-specific analysis for this protein
            protein_success = analyze_lipid_specific_effects(protein_combined_df, protein_matched_df, protein_output_dir)
            
            if protein_success:
                successful_analyses += 1
                print(f"  ✓ Successfully analyzed protein {protein}")
                
                # Rename the output files to include protein name
                _rename_output_files(protein_output_dir, protein)
                
            else:
                print(f"  ✗ Failed to analyze protein {protein}")
                
        except Exception as e:
            print(f"  ✗ Error analyzing protein {protein}: {str(e)}")
            continue
    
    print(f"\n✓ Per-peptide analysis completed: {successful_analyses}/{len(unique_proteins)} proteins analyzed successfully")
    print(f"Results saved to: {per_peptide_dir}")
    
    return successful_analyses > 0


def _rename_output_files(protein_output_dir, protein_name):
    """Rename output files to include protein name"""
    import glob
    import shutil
    
    # Rename lipid_comparison files
    for ext in ['png', 'svg', 'ps']:
        old_path = os.path.join(protein_output_dir, f'lipid_comparison.{ext}')
        new_path = os.path.join(protein_output_dir, f'lipid_comparison_{protein_name}.{ext}')
        if os.path.exists(old_path):
            shutil.move(old_path, new_path)
    
    # Rename other files if needed
    files_to_rename = [
        'lipid_specific_analysis_summary.txt',
        'lipid_specific_effects.csv'
    ]
    
    for filename in files_to_rename:
        old_path = os.path.join(protein_output_dir, filename)
        if os.path.exists(old_path):
            name_part, ext = os.path.splitext(filename)
            new_filename = f"{name_part}_{protein_name}{ext}"
            new_path = os.path.join(protein_output_dir, new_filename)
            shutil.move(old_path, new_path)