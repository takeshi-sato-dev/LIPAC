"""
Visualization functions for BayesianLipidAnalysis
All plotting and visualization utilities

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
import arviz as az


def plot_convergence_diagnostics(trace, model_name, output_dir):
    """Generate and save convergence diagnostic plots
    
    Parameters
    ----------
    trace : arviz.InferenceData
        MCMC trace from PyMC model
    model_name : str
        Name of the model for file naming
    output_dir : str
        Directory to save plots
    
    Returns
    -------
    pd.DataFrame or None
        Summary statistics dataframe
    """
    try:
        # 1. Trace plot
        fig = az.plot_trace(trace)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'{model_name}_trace_plot.png'), dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(output_dir, f'{model_name}_trace_plot.svg'), format='svg', bbox_inches='tight')
        plt.close()
        
        # 2. Autocorrelation plot
        fig = az.plot_autocorr(trace)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'{model_name}_autocorr.png'), dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(output_dir, f'{model_name}_autocorr.svg'), format='svg', bbox_inches='tight')
        plt.close()
        
        # 3. Summary statistics
        summary = az.summary(trace)
        print(f"\n{model_name} - Convergence Diagnostics:")
        print(f"R-hat values: {summary['r_hat'].to_dict()}")
        print(f"ESS bulk: {summary['ess_bulk'].to_dict()}")
        print(f"ESS tail: {summary['ess_tail'].to_dict()}")
        
        # Check convergence criteria
        r_hat_ok = (summary['r_hat'] < 1.01).all()
        ess_bulk_ok = (summary['ess_bulk'] > 400).all()
        ess_tail_ok = (summary['ess_tail'] > 400).all()
        
        if r_hat_ok and ess_bulk_ok and ess_tail_ok:
            print(f"✔ {model_name}: All convergence criteria met")
        else:
            print(f"⚠ {model_name}: Some convergence criteria not met")
            if not r_hat_ok:
                print(f"  - R-hat > 1.01 for some parameters")
            if not ess_bulk_ok:
                print(f"  - ESS bulk < 400 for some parameters")
            if not ess_tail_ok:
                print(f"  - ESS tail < 400 for some parameters")
        
        return summary
        
    except Exception as e:
        print(f"Error generating convergence diagnostics for {model_name}: {e}")
        import traceback
        traceback.print_exc()
        return None


def plot_residue_profile(merged_data, protein, partner, output_dir):
    """Create residue profile plot for protein pair
    
    Parameters
    ----------
    merged_data : pd.DataFrame
        Merged data with residue information
    protein : str
        Protein name
    partner : str
        Partner protein name
    output_dir : str
        Directory to save plots
    """
    # Sort by residue ID
    sorted_data = merged_data.sort_values('residue').reset_index(drop=True)
    
    # Create figure
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10), sharex=True)
    
    # X-axis: residue numbers
    x = sorted_data['residue'].values
    x_idx = np.arange(len(x))
    
    # First subplot: Contact frequencies without GM3
    ax1.bar(x_idx, sorted_data['protein_contact_no_gm3'], width=0.4, alpha=0.6, 
           label='Protein Contact (No GM3)', color='blue')
    ax1.axhline(y=np.percentile(sorted_data['protein_contact_no_gm3'], 75), color='blue', 
                linestyle='--', alpha=0.7, label='Interface Threshold')
    
    ax1.set_ylabel('Contact Frequency', fontsize=12)
    ax1.set_title(f'{protein}-{partner}: Residue Profile', fontsize=14)
    ax1.legend(loc='upper right')
    ax1.grid(True, alpha=0.3)
    
    # Second subplot: GM3 contact and protein contact difference
    ax2.bar(x_idx, sorted_data['protein_contact_diff'], width=0.4, alpha=0.6,
           color=np.where(sorted_data['protein_contact_diff'] < 0, 'red', 'green'),
           label='Protein Contact Change')
    ax2.axhline(y=0, color='black', linestyle='-', alpha=0.5)
    
    # Add GM3 contact line on shared axis
    ax3 = ax2.twinx()
    ax3.plot(x_idx, sorted_data['DPG3_contact'], 'r-', linewidth=2, label='GM3 Contact')
    ax3.axhline(y=np.percentile(sorted_data['DPG3_contact'], 75), color='red', 
                linestyle='--', alpha=0.7, label='GM3 Binding Threshold')
    
    # Highlight competitive binding sites
    for i, row in sorted_data.iterrows():
        if row.get('is_competition', False):
            ax2.axvspan(i-0.5, i+0.5, color='purple', alpha=0.2)
    
    # Formatting
    ax2.set_xlabel('Residue', fontsize=12)
    ax2.set_ylabel('Contact Change', fontsize=12)
    ax3.set_ylabel('GM3 Contact Frequency', fontsize=12, color='red')
    
    # Set x-tick labels to residue IDs
    plt.xticks(x_idx, x, rotation=90)
    
    # Add a combined legend for subplot 2
    lines_ax2, labels_ax2 = ax2.get_legend_handles_labels()
    lines_ax3, labels_ax3 = ax3.get_legend_handles_labels()
    ax2.legend(lines_ax2 + lines_ax3, labels_ax2 + labels_ax3, loc='upper right')
    
    # Add text for competitive sites
    if 'is_competition' in sorted_data.columns:
        comp_count = sorted_data['is_competition'].sum()
        if comp_count > 0:
            ax2.text(0.02, 0.95, f'Competitive sites: {comp_count}', transform=ax2.transAxes,
                    va='top', ha='left', bbox=dict(facecolor='white', alpha=0.7))
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{protein}_{partner}_residue_profile.png'), dpi=300)
    plt.savefig(os.path.join(output_dir, f'{protein}_{partner}_residue_profile.ps'))
    plt.close()


def plot_protein_contact_distribution(combined_df, figures_dir):
    """Plot protein-protein contact distribution"""
    plt.figure(figsize=(10, 6))
    sns.histplot(data=combined_df, x='protein_contact', hue='has_gm3', kde=True, bins=30, alpha=0.5)
    plt.title('Distribution of Protein-Protein Contacts')
    plt.xlabel('Protein-Protein Contact Frequency')
    plt.savefig(os.path.join(figures_dir, 'protein_contact_dist.png'), dpi=300)
    plt.close()


def plot_lipid_contact_distribution(combined_df, figures_dir):
    """Plot lipid-protein contact distribution"""
    plt.figure(figsize=(10, 6))
    sns.histplot(data=combined_df, x='lipid_contact', hue='has_gm3', kde=True, bins=30, alpha=0.5)
    plt.title('Distribution of Lipid-Protein Contacts')
    plt.xlabel('Lipid-Protein Contact Frequency')
    plt.savefig(os.path.join(figures_dir, 'lipid_contact_dist.png'), dpi=300)
    plt.close()


def plot_ratio_distribution(combined_df, figures_dir):
    """Plot lipid-to-protein contact ratio distribution"""
    plt.figure(figsize=(10, 6))
    sns.histplot(data=combined_df, x='ratio', hue='has_gm3', kde=True, bins=30, alpha=0.5)
    plt.title('Distribution of Lipid-to-Protein Contact Ratio')
    plt.xlabel('Lipid-to-Protein Contact Ratio')
    plt.savefig(os.path.join(figures_dir, 'ratio_dist.png'), dpi=300)
    plt.close()


def plot_gm3_correlation(matched_df, figures_dir):
    """Plot GM3 contact correlations"""
    # GM3 contact vs protein contact difference
    plt.figure(figsize=(10, 6))
    sns.regplot(data=matched_df, x='gm3_contact', y='protein_contact_diff', scatter_kws={'alpha': 0.5})
    plt.title('GM3 Contact vs Change in Protein-Protein Contact')
    plt.xlabel('GM3 Contact Frequency')
    plt.ylabel('Change in Protein-Protein Contact (With GM3 - Without GM3)')
    plt.axhline(y=0, color='r', linestyle='--')
    
    corr = matched_df[['gm3_contact', 'protein_contact_diff']].corr().iloc[0, 1]
    plt.text(0.05, 0.95, f'Correlation: {corr:.3f}', transform=plt.gca().transAxes, 
            bbox=dict(facecolor='white', alpha=0.7))
    
    plt.savefig(os.path.join(figures_dir, 'gm3_vs_protein_contact_diff.png'), dpi=300)
    plt.close()
    
    # GM3 contact vs lipid contact difference
    plt.figure(figsize=(10, 6))
    sns.regplot(data=matched_df, x='gm3_contact', y='lipid_contact_diff', scatter_kws={'alpha': 0.5})
    plt.title('GM3 Contact vs Change in Lipid-Protein Contact')
    plt.xlabel('GM3 Contact Frequency')
    plt.ylabel('Change in Lipid-Protein Contact (With GM3 - Without GM3)')
    plt.axhline(y=0, color='r', linestyle='--')
    
    corr = matched_df[['gm3_contact', 'lipid_contact_diff']].corr().iloc[0, 1]
    plt.text(0.05, 0.95, f'Correlation: {corr:.3f}', transform=plt.gca().transAxes,
            bbox=dict(facecolor='white', alpha=0.7))
    
    plt.savefig(os.path.join(figures_dir, 'gm3_vs_lipid_contact_diff.png'), dpi=300)
    plt.close()
    
    # GM3 contact vs ratio difference
    plt.figure(figsize=(10, 6))
    sns.regplot(data=matched_df, x='gm3_contact', y='ratio_diff', scatter_kws={'alpha': 0.5})
    plt.title('GM3 Contact vs Change in Contact Ratio')
    plt.xlabel('GM3 Contact Frequency')
    plt.ylabel('Change in Ratio (With GM3 - Without GM3)')
    plt.axhline(y=0, color='r', linestyle='--')
    
    corr = matched_df[['gm3_contact', 'ratio_diff']].corr().iloc[0, 1]
    plt.text(0.05, 0.95, f'Correlation: {corr:.3f}', transform=plt.gca().transAxes,
            bbox=dict(facecolor='white', alpha=0.7))
    
    plt.savefig(os.path.join(figures_dir, 'gm3_vs_ratio_diff.png'), dpi=300)
    plt.close()


def plot_gm3_effect_violin(gm3_high, gm3_low, ttest_result, figures_dir):
    """Plot violin plot comparing GM3 effects"""
    plt.figure(figsize=(12, 6))
    plot_data = pd.concat([
        pd.DataFrame({'GM3 Contact': 'High', 'Protein Contact Diff': gm3_high['protein_contact_diff']}),
        pd.DataFrame({'GM3 Contact': 'Low/None', 'Protein Contact Diff': gm3_low['protein_contact_diff']})
    ])
    sns.violinplot(data=plot_data, x='GM3 Contact', y='Protein Contact Diff')
    plt.title(f'Effect of GM3 on Protein-Protein Contacts\np = {ttest_result.pvalue:.2e}')
    plt.axhline(y=0, color='r', linestyle='--')
    plt.savefig(os.path.join(figures_dir, 'gm3_effect_violin.png'), dpi=300)
    plt.close()


def plot_bayesian_regression(X, Y, trace, param_name, output_dir, model_name, ylabel, color='r'):
    """Plot Bayesian regression results
    
    Parameters
    ----------
    X : np.ndarray
        Predictor variable
    Y : np.ndarray
        Response variable
    trace : arviz.InferenceData
        MCMC trace
    param_name : str
        Name of parameter to plot
    output_dir : str
        Output directory
    model_name : str
        Model name for file naming
    ylabel : str
        Y-axis label
    color : str
        Color for regression line
    """
    try:
        # Extract parameters
        alpha_samples = trace.posterior['alpha'].values.flatten()
        beta_samples = trace.posterior[param_name].values.flatten()
        
        plt.figure(figsize=(10, 6))
        plt.scatter(X, Y, alpha=0.5)
        
        # Average regression line
        alpha_mean = alpha_samples.mean()
        beta_mean = beta_samples.mean()
        
        x_range = np.linspace(0, np.max(X) * 1.1, 100)
        plt.plot(x_range, alpha_mean + beta_mean * x_range, f'{color}-', linewidth=2, 
                label=f'y = {alpha_mean:.3f} + {beta_mean:.3f}x')
        
        plt.axhline(y=0, color='k', linestyle='--', alpha=0.3)
        plt.xlabel('GM3 Contact Frequency')
        plt.ylabel(ylabel)
        plt.title(f'Bayesian Regression: {model_name}')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.savefig(os.path.join(output_dir, f'{model_name}_regression.png'), dpi=300)
        plt.savefig(os.path.join(output_dir, f'{model_name}_regression.svg'), dpi=300)
        plt.close()
        
        # Effect size posterior distribution
        plt.figure(figsize=(8, 6))
        plt.hist(beta_samples, bins=30, alpha=0.7, density=True)
        plt.axvline(x=0, color='r', linestyle='--')
        plt.xlabel(f'Effect Size ({param_name})')
        plt.ylabel('Density')
        plt.title(f'Posterior Distribution: {model_name}')
        
        # 95% HDI
        beta_hdi = np.percentile(beta_samples, [2.5, 97.5])
        plt.axvline(x=beta_hdi[0], color='k', linestyle='--', alpha=0.5)
        plt.axvline(x=beta_hdi[1], color='k', linestyle='--', alpha=0.5)
        
        # Statistical information
        prob_effect = (beta_samples < 0).mean() if beta_mean < 0 else (beta_samples > 0).mean()
        direction = '<' if beta_mean < 0 else '>'
        
        plt.text(0.05, 0.95, f'95% HDI: [{beta_hdi[0]:.3f}, {beta_hdi[1]:.3f}]', 
                transform=plt.gca().transAxes, bbox=dict(facecolor='white', alpha=0.7))
        plt.text(0.05, 0.85, f'P(β {direction} 0) = {prob_effect:.3f}', 
                transform=plt.gca().transAxes, bbox=dict(facecolor='white', alpha=0.7))
        
        plt.grid(True, alpha=0.3)
        plt.savefig(os.path.join(output_dir, f'{model_name}_effect.png'), dpi=300)
        plt.savefig(os.path.join(output_dir, f'{model_name}_effect.svg'), dpi=300)
        plt.close()
        
    except Exception as e:
        print(f"Error in Bayesian regression plot: {str(e)}")


def plot_lipid_comparison(results_df, output_dir):
    """Create comparison bar plot for lipid-specific effects
    
    Parameters
    ----------
    results_df : pd.DataFrame
        DataFrame with lipid analysis results
    output_dir : str
        Output directory
    """
    plt.figure(figsize=(12, 8))
    
    # Create error bars for each lipid type
    y_pos = np.arange(len(results_df))
    
    # Color based on effect direction and certainty
    colors = []
    for _, row in results_df.iterrows():
        if row['effect_direction'] == 'positive':
            colors.append('darkgreen' if row['prob_effect'] > 0.95 else 'lightgreen')
        else:
            colors.append('darkred' if row['prob_effect'] > 0.95 else 'lightcoral')
    
    # Plot bars
    plt.barh(y_pos, results_df['mean_effect'], 
            xerr=np.vstack([
                results_df['mean_effect'] - results_df['hdi_lower'], 
                results_df['hdi_upper'] - results_df['mean_effect']
            ]), 
            align='center', color=colors, alpha=0.7)
    
    # Reference line
    plt.axvline(x=0, color='k', linestyle='-', alpha=0.3)
    
    # Labels
    plt.yticks(y_pos, results_df['lipid_type'])
    plt.xlabel('GM3 Effect Size (β)')
    plt.title('GM3 Effects on Different Lipid Types')
    
    # Add probability annotations
    for i, (_, row) in enumerate(results_df.iterrows()):
        direction = "+" if row['effect_direction'] == 'positive' else "-"
        plt.text(0, i, f"  P({direction}): {row['prob_effect']:.3f}", 
                va='center', ha='left' if row['mean_effect'] < 0 else 'right',
                fontsize=9, backgroundcolor='white', alpha=0.7)
    
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'lipid_comparison.png'), dpi=300)
    plt.savefig(os.path.join(output_dir, 'lipid_comparison.svg'), format='svg', dpi=300)
    plt.savefig(os.path.join(output_dir, 'lipid_comparison.ps'))
    plt.close()


def plot_interface_enrichment(pair_stats_df, output_dir):
    """Plot enrichment of GM3 binding sites at interface residues
    
    Parameters
    ----------
    pair_stats_df : pd.DataFrame
        Statistics for each protein pair
    output_dir : str
        Output directory
    """
    plt.figure(figsize=(10, 6))
    bar_colors = ['green' if p < 0.05 else 'gray' for p in pair_stats_df['p_value']]
    
    bars = plt.bar(
        pair_stats_df['pair'], 
        pair_stats_df['enrichment'],
        color=bar_colors, alpha=0.7
    )
    
    # Reference line for no enrichment
    plt.axhline(y=1.0, linestyle='--', color='red', alpha=0.7)
    
    # Mark statistical significance
    for i, p in enumerate(pair_stats_df['p_value']):
        sig = ''
        if p < 0.001:
            sig = '***'
        elif p < 0.01:
            sig = '**'
        elif p < 0.05:
            sig = '*'
        
        if sig:
            plt.text(
                i, 
                pair_stats_df['enrichment'].iloc[i] + 0.1, 
                sig, 
                ha='center', 
                fontsize=12
            )
    
    plt.ylabel('Enrichment of GM3 Binding Sites at Interface Residues', fontsize=12)
    plt.xlabel('Protein Pair', fontsize=12)
    plt.title('Analysis of GM3 Binding Site and Interface Residue Overlap', fontsize=14)
    plt.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'gm3_interface_enrichment.png'), dpi=300)
    plt.savefig(os.path.join(output_dir, 'gm3_interface_enrichment.ps'))
    plt.close()


def plot_residue_distribution_pie(pair_stats_df, output_dir):
    """Create pie charts showing residue distribution
    
    Parameters
    ----------
    pair_stats_df : pd.DataFrame
        Statistics for each protein pair
    output_dir : str
        Output directory
    """
    plt.figure(figsize=(15, 5))
    for i, (_, row) in enumerate(pair_stats_df.iterrows()):
        plt.subplot(1, len(pair_stats_df), i+1)
        
        # Competition sites, non-competition interface, non-competition GM3, other
        competition = row['competition_sites']
        interface_only = row['interface_residues'] - competition
        gm3_only = row['gm3_binding_residues'] - competition
        other = row['total_residues'] - interface_only - gm3_only - competition
        
        plt.pie(
            [competition, interface_only, gm3_only, other],
            labels=[f'Competition ({competition:.0f})', 
                   f'Interface Only ({interface_only:.0f})', 
                   f'GM3 Binding Only ({gm3_only:.0f})', 
                   f'Other ({other:.0f})'],
            autopct='%1.1f%%',
            colors=['red', 'blue', 'green', 'lightgray'],
            startangle=90
        )
        plt.title(row['pair'], fontsize=12)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'residue_distribution.png'), dpi=300)
    plt.savefig(os.path.join(output_dir, 'residue_distribution.ps'))
    plt.close()


def plot_hierarchical_forest(hierarchical_trace, proteins, output_dir):
    """Create forest plot for hierarchical model results
    
    Parameters
    ----------
    hierarchical_trace : arviz.InferenceData
        Hierarchical model trace
    proteins : list
        List of protein names
    output_dir : str
        Output directory
    """
    plt.figure(figsize=(10, 8))
    
    # Calculate effect for each protein
    means = []
    lower_cis = []
    upper_cis = []
    protein_names = []
    
    for i, protein in enumerate(proteins):
        beta_samples = hierarchical_trace.posterior['beta'][:, :, i].values.flatten()
        mean = beta_samples.mean()
        ci = np.percentile(beta_samples, [2.5, 97.5])
        
        means.append(mean)
        lower_cis.append(ci[0])
        upper_cis.append(ci[1])
        protein_names.append(f'Protein {i+1}: {protein}')
    
    # Add global effect
    mu_beta_samples = hierarchical_trace.posterior['mu_beta'].values.flatten()
    mu_beta_mean = mu_beta_samples.mean()
    mu_beta_ci = np.percentile(mu_beta_samples, [2.5, 97.5])
    
    means.append(mu_beta_mean)
    lower_cis.append(mu_beta_ci[0])
    upper_cis.append(mu_beta_ci[1])
    protein_names.append('Global')
    
    # Display bottom to top
    y_pos = np.arange(len(protein_names))
    
    # Create plot
    plt.errorbar(means, y_pos, 
                xerr=[np.array(means)-np.array(lower_cis), 
                      np.array(upper_cis)-np.array(means)], 
                fmt='o', capsize=5, elinewidth=2, markeredgewidth=2)
    
    # Zero line
    plt.axvline(x=0, color='r', linestyle='--', alpha=0.5)
    
    # Label each protein
    plt.yticks(y_pos, protein_names)
    
    plt.xlabel('Effect Size (β)')
    plt.title('GM3 Effect on Protein-Protein Contacts by Protein')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'hierarchical_model_forest.png'), dpi=300)
    plt.close()