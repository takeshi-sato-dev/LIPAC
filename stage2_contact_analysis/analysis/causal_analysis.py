"""
Causal Bayesian analysis for target lipid-induced lipid composition changes

This module implements causal inference to determine how target lipid binding
affects the composition of other lipids around proteins.

Author: Takeshi Sato, PhD
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


def load_causal_data(data_dir, target_lipid_name='target_lipid'):
    """Load target lipid causal data from stage1 output
    
    Parameters
    ----------
    data_dir : str
        Directory containing target lipid causal data files
    target_lipid_name : str
        Name of the target lipid (e.g., 'dpg3', 'chol', etc.)
    
    Returns
    -------
    dict
        Dictionary with protein names as keys and DataFrames as values
    """
    print(f"\n===== Loading {target_lipid_name.upper()} Causal Data =====")
    
    causal_data = {}
    
    # Load individual protein data
    import glob
    pattern = os.path.join(data_dir, f'{target_lipid_name.lower()}_causal_data_Protein_*.csv')
    files = glob.glob(pattern)
    
    if not files:
        # Try loading the combined file
        combined_file = os.path.join(data_dir, f'{target_lipid_name.lower()}_causal_data_all_proteins.csv')
        if os.path.exists(combined_file):
            print(f"Loading combined data from: {combined_file}")
            all_df = pd.read_csv(combined_file)
            
            # Split by protein
            for protein in all_df['protein'].unique():
                causal_data[protein] = all_df[all_df['protein'] == protein].copy()
                print(f"  Loaded {protein}: {len(causal_data[protein])} frames")
        else:
            print(f"Error: No {target_lipid_name.upper()} causal data files found!")
            return None
    else:
        for file in files:
            protein_name = os.path.basename(file).replace(f'{target_lipid_name.lower()}_causal_data_', '').replace('.csv', '')
            df = pd.read_csv(file)
            causal_data[protein_name] = df
            print(f"  Loaded {protein_name}: {len(df)} frames, "
                  f"{target_lipid_name.upper()} bound in {df['target_lipid_bound'].sum()} frames ({df['target_lipid_bound'].mean()*100:.1f}%)")
    
    return causal_data


def perform_causal_bayesian_analysis(causal_data, output_dir, target_lipid_name='TARGET_LIPID'):
    """Perform causal Bayesian analysis for target lipid effects on lipid composition
    
    Parameters
    ----------
    causal_data : dict
        Dictionary with protein DataFrames
    output_dir : str
        Output directory for results
    target_lipid_name : str
        Name of the target lipid for analysis
    
    Returns
    -------
    dict
        Analysis results
    """
    print("\n===== Causal Bayesian Analysis =====")
    print(f"Analyzing causal effects of {target_lipid_name} binding on lipid composition")
    
    # Create output directory
    causal_dir = os.path.join(output_dir, 'causal_analysis')
    os.makedirs(causal_dir, exist_ok=True)
    
    results = {}
    
    # Analyze each protein separately
    for protein_name, df in causal_data.items():
        print(f"\n--- Analyzing {protein_name} ---")
        
        # Get lipid columns (exclude the target lipid itself)
        target_lipid_col = f'{target_lipid_name}_contacts'
        residue_contact_cols = [col for col in df.columns 
                               if col.endswith('_contacts') and col != target_lipid_col]
        unique_molecule_cols = [col for col in df.columns 
                               if col.endswith('_unique_molecules')]
        
        # Combine both types of columns
        lipid_cols = residue_contact_cols + unique_molecule_cols
        
        if not lipid_cols:
            print(f"  Warning: No lipid contact columns found for {protein_name}")
            continue
        
        protein_results = {}
        
        # Analyze each lipid type
        for lipid_col in lipid_cols:
            if lipid_col.endswith('_contacts'):
                lipid_name = lipid_col.replace('_contacts', '') + '_contacts'
                analysis_type = 'residue_contacts'
            elif lipid_col.endswith('_unique_molecules'):
                lipid_name = lipid_col.replace('_unique_molecules', '') + '_unique'
                analysis_type = 'unique_molecules'
            else:
                continue
            print(f"  Analyzing {lipid_name}...")
            
            # Prepare data
            target_lipid_bound = df['target_lipid_bound'].astype(int).values
            lipid_contacts = df[lipid_col].values
            
            # Remove NaN values
            valid_mask = ~np.isnan(lipid_contacts)
            target_lipid_bound = target_lipid_bound[valid_mask]
            lipid_contacts = lipid_contacts[valid_mask]
            
            if len(lipid_contacts) < 20:
                print(f"    Insufficient data for {lipid_name} (n={len(lipid_contacts)})")
                continue
            
            # Check for zero variance (all values are the same)
            std_lipid = np.std(lipid_contacts)
            if std_lipid < 1e-6:  # Essentially zero variance
                print(f"    Skipping {lipid_name}: zero variance (all values ≈ {np.mean(lipid_contacts):.3f})")
                continue
            
            # Bayesian causal model
            with pm.Model() as causal_model:
                # Priors with minimum sigma values
                alpha = pm.Normal('alpha', mu=np.mean(lipid_contacts), sigma=max(std_lipid, 1.0))
                beta_target_lipid = pm.Normal('beta_target_lipid', mu=0, sigma=max(std_lipid, 1.0))
                sigma = pm.HalfNormal('sigma', sigma=max(std_lipid, 1.0))
                
                # Linear model: contacts = alpha + beta_target_lipid * target_lipid_bound
                mu = alpha + beta_target_lipid * target_lipid_bound
                
                # Likelihood
                y = pm.Normal('y', mu=mu, sigma=sigma, observed=lipid_contacts)
                
                # Sample
                trace = pm.sample(2000, tune=1000, chains=4, random_seed=42, 
                                progressbar=False, return_inferencedata=True)
            
            # Extract results
            beta_samples = trace.posterior['beta_target_lipid'].values.flatten()
            
            # Calculate statistics
            effect_mean = np.mean(beta_samples)
            effect_std = np.std(beta_samples)
            ci_95 = np.percentile(beta_samples, [2.5, 97.5])
            prob_positive = (beta_samples > 0).mean()
            prob_negative = (beta_samples < 0).mean()
            
            # Calculate average lipid contacts with and without target lipid
            avg_with_target_lipid = lipid_contacts[target_lipid_bound == 1].mean() if np.any(target_lipid_bound == 1) else np.nan
            avg_without_target_lipid = lipid_contacts[target_lipid_bound == 0].mean() if np.any(target_lipid_bound == 0) else np.nan
            
            protein_results[lipid_name] = {
                'effect_mean': effect_mean,
                'effect_std': effect_std,
                'ci_95': ci_95,
                'prob_positive': prob_positive,
                'prob_negative': prob_negative,
                'avg_with_target_lipid': avg_with_target_lipid,
                'avg_without_target_lipid': avg_without_target_lipid,
                'n_with_target_lipid': np.sum(target_lipid_bound == 1),
                'n_without_target_lipid': np.sum(target_lipid_bound == 0),
                'trace': trace,
                'beta_samples': beta_samples
            }
            
            print(f"    Effect: {effect_mean:.2f} ± {effect_std:.2f}")
            print(f"    95% CI: [{ci_95[0]:.2f}, {ci_95[1]:.2f}]")
            print(f"    P(effect > 0): {prob_positive:.3f}")
        
        results[protein_name] = protein_results
        
        # Generate protein-specific plots
        _generate_causal_plots(protein_results, protein_name, causal_dir, target_lipid_name)
    
    # Generate comparative analysis
    _generate_comparative_causal_plots(results, causal_dir, target_lipid_name)
    
    # Generate horizontal bar chart
    _generate_causal_bar_chart(results, causal_dir, target_lipid_name)
    
    # Save summary
    _save_causal_summary(results, causal_dir, target_lipid_name)
    
    return results


def _generate_causal_plots(protein_results, protein_name, output_dir, target_lipid_name='TARGET_LIPID'):
    """Generate causal effect plots for a single protein"""
    
    if not protein_results:
        return
    
    # Create figure with subplots
    n_lipids = len(protein_results)
    fig, axes = plt.subplots(1, n_lipids, figsize=(4*n_lipids, 5))
    
    if n_lipids == 1:
        axes = [axes]
    
    for idx, (lipid_name, results) in enumerate(protein_results.items()):
        ax = axes[idx]
        
        # Plot posterior distribution
        beta_samples = results['beta_samples']
        ax.hist(beta_samples, bins=50, alpha=0.7, color='blue', edgecolor='black')
        ax.axvline(0, color='red', linestyle='--', alpha=0.5)
        ax.axvline(results['effect_mean'], color='green', linestyle='-', linewidth=2)
        
        ax.set_xlabel(f'{target_lipid_name} Effect on {lipid_name}')
        ax.set_ylabel('Density')
        ax.set_title(f'{lipid_name}\nEffect: {results["effect_mean"]:.2f} ± {results["effect_std"]:.2f}')
        
        # Add probability text
        prob_text = f'P(>0): {results["prob_positive"]:.3f}\nP(<0): {results["prob_negative"]:.3f}'
        ax.text(0.05, 0.95, prob_text, transform=ax.transAxes, 
               verticalalignment='top', fontsize=10,
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.suptitle(f'{target_lipid_name} Causal Effects - {protein_name}', fontsize=14)
    plt.tight_layout()
    
    # Save figure
    plot_path = os.path.join(output_dir, f'causal_effects_{protein_name}.png')
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.savefig(plot_path.replace('.png', '.svg'), format='svg', bbox_inches='tight')
    plt.close()


def _generate_comparative_causal_plots(all_results, output_dir, target_lipid_name='TARGET_LIPID'):
    """Generate separate comparative plots for contacts and unique molecules"""
    
    # Prepare data for comparison, separating contacts and unique molecules
    contact_data = []
    unique_data = []
    
    for protein_name, protein_results in all_results.items():
        for lipid_name, results in protein_results.items():
            data_point = {
                'protein': protein_name,
                'lipid': lipid_name,
                'effect': results['effect_mean'],
                'ci_lower': results['ci_95'][0],
                'ci_upper': results['ci_95'][1],
                'prob_positive': results['prob_positive']
            }
            
            if lipid_name.endswith('_contacts'):
                # Remove _contacts suffix for cleaner display
                data_point['lipid'] = lipid_name.replace('_contacts', '')
                contact_data.append(data_point)
            elif lipid_name.endswith('_unique'):
                # Remove _unique suffix for cleaner display
                data_point['lipid'] = lipid_name.replace('_unique', '')
                unique_data.append(data_point)
    
    # Generate contacts plot
    if contact_data:
        _generate_single_comparison_plot(contact_data, output_dir, target_lipid_name, 
                                       'contacts', 'Residue Contacts')
    
    # Generate unique molecules plot
    if unique_data:
        _generate_single_comparison_plot(unique_data, output_dir, target_lipid_name, 
                                       'unique_molecules', 'Unique Molecules')


def _generate_single_comparison_plot(comparison_data, output_dir, target_lipid_name, 
                                   plot_type, plot_label):
    """Generate a single comparison plot for either contacts or unique molecules"""
    
    df = pd.DataFrame(comparison_data)
    
    # Create comparison plot
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Group by lipid type
    lipids = df['lipid'].unique()
    proteins = sorted(df['protein'].unique())  # Sort proteins to ensure 1,2,3,4 order
    
    x = np.arange(len(lipids))
    width = 0.8 / len(proteins)
    
    colors = plt.cm.Set3(np.linspace(0, 1, len(proteins)))
    
    for i, protein in enumerate(proteins):
        protein_df = df[df['protein'] == protein]
        effects = []
        errors = []
        
        for lipid in lipids:
            lipid_df = protein_df[protein_df['lipid'] == lipid]
            if not lipid_df.empty:
                effects.append(lipid_df['effect'].values[0])
                errors.append([
                    lipid_df['effect'].values[0] - lipid_df['ci_lower'].values[0],
                    lipid_df['ci_upper'].values[0] - lipid_df['effect'].values[0]
                ])
            else:
                effects.append(0)
                errors.append([0, 0])
        
        errors = np.array(errors).T
        positions = x + i * width - (len(proteins) - 1) * width / 2
        
        ax.bar(positions, effects, width, label=protein, color=colors[i], alpha=0.8)
        ax.errorbar(positions, effects, yerr=errors, fmt='none', color='black', capsize=3)
    
    ax.set_xlabel('Lipid Type')
    ax.set_ylabel(f'{target_lipid_name} Causal Effect ({plot_label})')
    ax.set_title(f'{target_lipid_name} Causal Effects on {plot_label} Across Proteins')
    ax.set_xticks(x)
    ax.set_xticklabels(lipids)
    ax.axhline(y=0, color='red', linestyle='--', alpha=0.5)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save figure with specific naming
    plot_path = os.path.join(output_dir, f'causal_effects_comparison_{plot_type}.png')
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.savefig(plot_path.replace('.png', '.svg'), format='svg', bbox_inches='tight')
    plt.close()
    print(f"Comparative {plot_label} plot saved: {plot_path}")


def _save_causal_summary(results, output_dir, target_lipid_name='TARGET_LIPID'):
    """Save causal analysis summary to text file"""
    
    summary_path = os.path.join(output_dir, 'causal_analysis_summary.txt')
    
    with open(summary_path, 'w') as f:
        f.write(f"===== {target_lipid_name} Causal Bayesian Analysis Summary =====\n")
        f.write("=" * 60 + "\n\n")
        
        for protein_name, protein_results in results.items():
            f.write(f"\n--- {protein_name} ---\n")
            
            for lipid_name, res in protein_results.items():
                f.write(f"\n{lipid_name}:\n")
                f.write(f"  Causal Effect: {res['effect_mean']:.3f} ± {res['effect_std']:.3f}\n")
                f.write(f"  95% CI: [{res['ci_95'][0]:.3f}, {res['ci_95'][1]:.3f}]\n")
                f.write(f"  P(effect > 0): {res['prob_positive']:.3f}\n")
                f.write(f"  P(effect < 0): {res['prob_negative']:.3f}\n")
                f.write(f"  Average contacts with {target_lipid_name}: {res['avg_with_target_lipid']:.2f} (n={res['n_with_target_lipid']})\n")
                f.write(f"  Average contacts without {target_lipid_name}: {res['avg_without_target_lipid']:.2f} (n={res['n_without_target_lipid']})\n")
                
                # Interpretation
                if res['prob_positive'] > 0.95:
                    f.write(f"  → Strong evidence that {target_lipid_name} increases {lipid_name} contacts\n")
                elif res['prob_negative'] > 0.95:
                    f.write(f"  → Strong evidence that {target_lipid_name} decreases {lipid_name} contacts\n")
                elif res['prob_positive'] > 0.75:
                    f.write(f"  → Moderate evidence that {target_lipid_name} increases {lipid_name} contacts\n")
                elif res['prob_negative'] > 0.75:
                    f.write(f"  → Moderate evidence that {target_lipid_name} decreases {lipid_name} contacts\n")
                else:
                    f.write(f"  → Weak or no evidence for {target_lipid_name} effect on {lipid_name}\n")
        
        f.write("\n" + "=" * 60 + "\n")
        f.write("CAUSAL INTERPRETATION:\n")
        f.write(f"This analysis shows the direct causal effect of {target_lipid_name} binding on lipid composition.\n")
        f.write(f"Positive effects indicate {target_lipid_name} binding increases contacts with that lipid type.\n")
        f.write(f"Negative effects indicate {target_lipid_name} binding decreases contacts with that lipid type.\n")
    
    print(f"\nSummary saved to: {summary_path}")


def _generate_causal_bar_chart(results, output_dir, target_lipid_name='TARGET_LIPID'):
    """Generate horizontal bar chart of causal effects"""
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    
    # Collect data for bar chart
    proteins = []
    lipids = []
    effects = []
    ci_lower = []
    ci_upper = []
    colors = []
    
    for protein_name, protein_results in results.items():
        for lipid_name, lipid_results in protein_results.items():
            proteins.append(protein_name)
            lipids.append(lipid_name)
            effects.append(lipid_results['effect_mean'])
            ci_lower.append(lipid_results['ci_95'][0])
            ci_upper.append(lipid_results['ci_95'][1])
            
            # Color based on effect direction and significance
            if lipid_results['prob_positive'] > 0.95:  # Strong positive effect
                colors.append('green')
            elif lipid_results['prob_negative'] > 0.95:  # Strong negative effect
                colors.append('red')
            elif lipid_results['prob_positive'] > 0.8:  # Moderate positive
                colors.append('lightgreen')
            elif lipid_results['prob_negative'] > 0.8:  # Moderate negative
                colors.append('lightcoral')
            else:  # Uncertain
                colors.append('gray')
    
    if not effects:
        return
    
    # Create labels
    labels = [f"{p}_{l}" for p, l in zip(proteins, lipids)]
    
    # Create horizontal bar chart
    fig, ax = plt.subplots(figsize=(10, max(6, len(labels) * 0.4)))
    
    y_pos = range(len(labels))
    bars = ax.barh(y_pos, effects, color=colors, alpha=0.7, edgecolor='black')
    
    # Add error bars (confidence intervals)
    xerr_lower = [e - c_l for e, c_l in zip(effects, ci_lower)]
    xerr_upper = [c_u - e for e, c_u in zip(ci_upper, effects)]
    ax.errorbar(effects, y_pos, xerr=[xerr_lower, xerr_upper], 
                fmt='none', color='black', capsize=3, capthick=1)
    
    # Add vertical line at x=0
    ax.axvline(0, color='black', linestyle='--', alpha=0.5)
    
    # Formatting
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels)
    ax.set_xlabel(f'{target_lipid_name} Causal Effect (Contact Count Change)')
    ax.set_title(f'{target_lipid_name} Causal Effects on Lipid Contacts\n(with 95% Confidence Intervals)')
    
    # Add legend
    legend_elements = [
        Patch(facecolor='green', alpha=0.7, label='Strong Positive (P>0.95)'),
        Patch(facecolor='red', alpha=0.7, label='Strong Negative (P<0.05)'),
        Patch(facecolor='lightgreen', alpha=0.7, label='Moderate Positive (P>0.8)'),
        Patch(facecolor='lightcoral', alpha=0.7, label='Moderate Negative (P<0.2)'),
        Patch(facecolor='gray', alpha=0.7, label='Uncertain')
    ]
    ax.legend(handles=legend_elements, loc='lower right')
    
    # Add text annotations for values
    for i, (effect, bar) in enumerate(zip(effects, bars)):
        ax.text(effect + 0.1 if effect >= 0 else effect - 0.1, i, f'{effect:.1f}', 
                va='center', ha='left' if effect >= 0 else 'right', fontsize=9)
    
    plt.tight_layout()
    plot_path = os.path.join(output_dir, f'{target_lipid_name.lower()}_causal_effects_barplot.png')
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.savefig(plot_path.replace('.png', '.svg'), format='svg', bbox_inches='tight')
    plt.close()
    print(f"Causal effects bar chart saved: {plot_path} and {plot_path.replace('.png', '.svg')}")