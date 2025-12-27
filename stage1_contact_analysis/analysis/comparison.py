#!/usr/bin/env python3
"""Comparison analysis functions for with/without target lipid systems."""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# Import the new competition plotting functions
from ..visualization.competition_plots import (
    plot_energy_comparison,
    plot_target_lipid_protein_interface_competition,
    plot_posterior_distribution_target_lipid_on_protein_contacts,
    plot_bayesian_regression_target_lipid_effect,
    plot_protein_model_effect_simple
)

def compare_with_without_target_lipid(with_lipid_df, without_lipid_df, output_dir, target_lipid_col='DPG3_contact'):
    """Compare systems with and without target lipid
    
    Parameters
    ----------
    with_lipid_df : pd.DataFrame
        Contact data from system with target lipid
    without_lipid_df : pd.DataFrame
        Contact data from system without target lipid
    output_dir : str
        Output directory for comparison plots
    target_lipid_col : str
        Column name for target lipid contact data
    """
    print("\n===== Comparing Systems With and Without Target Lipid =====")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Check for missing data
    if with_lipid_df is None or without_lipid_df is None:
        print("Missing data for comparison. Aborting.")
        return
    
    # Add lipid presence flags
    with_lipid_df = with_lipid_df.copy()
    without_lipid_df = without_lipid_df.copy()
    with_lipid_df['has_target_lipid'] = True
    without_lipid_df['has_target_lipid'] = False
    
    # Merge keys
    merge_keys = ['protein', 'residue']
    
    # Prepare subsets for merging
    with_lipid_subset = with_lipid_df[merge_keys + ['protein_contact', 'lipid_contact', 'ratio', 'has_target_lipid', target_lipid_col]]
    without_lipid_subset = without_lipid_df[merge_keys + ['protein_contact', 'lipid_contact', 'ratio', 'has_target_lipid']]
    
    # Rename columns to distinguish with/without target lipid
    with_lipid_subset = with_lipid_subset.rename(columns={
        'protein_contact': 'protein_contact_with_lipid',
        'lipid_contact': 'lipid_contact_with_lipid',
        'ratio': 'ratio_with_lipid'
    })
    
    without_lipid_subset = without_lipid_subset.rename(columns={
        'protein_contact': 'protein_contact_without_lipid',
        'lipid_contact': 'lipid_contact_without_lipid',
        'ratio': 'ratio_without_lipid'
    })
    
    # Outer merge to keep residues from both datasets
    merged_df = pd.merge(with_lipid_subset, without_lipid_subset, on=merge_keys, how='outer')
    
    # Fill NaN values with 0
    merged_df = merged_df.fillna(0)
    
    # Calculate differences
    merged_df['ratio_diff'] = merged_df['ratio_with_lipid'] - merged_df['ratio_without_lipid']
    merged_df['protein_contact_diff'] = merged_df['protein_contact_with_lipid'] - merged_df['protein_contact_without_lipid']
    merged_df['lipid_contact_diff'] = merged_df['lipid_contact_with_lipid'] - merged_df['lipid_contact_without_lipid']
    
    # Target lipid contact classification
    if target_lipid_col in merged_df.columns:
        merged_df['has_target_contact'] = merged_df[target_lipid_col] > 0.04
    
    # ----- 1. Target lipid contact vs ratio difference scatter plot -----
    if target_lipid_col in merged_df.columns:
        plt.figure(figsize=(10, 8))
        
        scatter = plt.scatter(
            merged_df[target_lipid_col],
            merged_df['ratio_diff'],
            c=merged_df['protein_contact_with_lipid'],
            cmap='viridis',
            alpha=0.7,
            s=50
        )
        
        # Colorbar
        cbar = plt.colorbar(scatter)
        cbar.set_label('Protein-Protein Contact Frequency')
        
        # Linear regression line
        try:
            # Only consider residues with target lipid contact
            target_contact_df = merged_df[merged_df[target_lipid_col] > 0]
            
            slope, intercept, r_value, p_value, std_err = stats.linregress(
                target_contact_df[target_lipid_col],
                target_contact_df['ratio_diff']
            )
            
            # Plot regression line and R² value
            x = np.array([0, target_contact_df[target_lipid_col].max()])
            y = slope * x + intercept
            plt.plot(x, y, 'r--', linewidth=2)
            
            plt.text(
                0.05, 0.95,
                f'R² = {r_value**2:.3f}\np = {p_value:.3e}',
                transform=plt.gca().transAxes,
                bbox=dict(facecolor='white', alpha=0.7)
            )
        except Exception as e:
            print(f"Warning: Could not calculate regression line: {str(e)}")
        
        # Labels and title
        plt.xlabel('Target Lipid Contact Frequency')
        plt.ylabel('Ratio Difference (With Target Lipid - Without Target Lipid)')
        plt.title('Effect of Target Lipid Contact on Contact Ratio')
        plt.grid(True, alpha=0.3)
        
        # Add zero line
        plt.axhline(y=0, color='k', linestyle='-', alpha=0.3)
        
        # Save
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "target_lipid_contact_ratio_diff.png"), dpi=300)
        plt.savefig(os.path.join(output_dir, "target_lipid_contact_ratio_diff.ps"))
        plt.close()
    
    # ----- 2. Protein contact change vs lipid contact change scatter plot -----
    plt.figure(figsize=(10, 8))
    
    # Color and size based on target lipid contact
    colors = []
    sizes = []
    
    if target_lipid_col in merged_df.columns:
        for target_contact in merged_df[target_lipid_col]:
            if target_contact > 1.0:  # Strong target lipid contact
                colors.append('red')
                sizes.append(100)
            elif target_contact > 0.04:  # Weak target lipid contact
                colors.append('orange')
                sizes.append(70)
            else:  # No target lipid contact
                colors.append('blue')
                sizes.append(30)
    else:
        colors = ['blue'] * len(merged_df)
        sizes = [50] * len(merged_df)
    
    # Scatter plot
    plt.scatter(
        merged_df['protein_contact_diff'],
        merged_df['lipid_contact_diff'],
        c=colors,
        alpha=0.7,
        s=sizes
    )
    
    # Linear regression line
    try:
        slope, intercept, r_value, p_value, std_err = stats.linregress(
            merged_df['protein_contact_diff'],
            merged_df['lipid_contact_diff']
        )
        
        # Plot regression line and R² value
        x = np.array([merged_df['protein_contact_diff'].min(), merged_df['protein_contact_diff'].max()])
        y = slope * x + intercept
        plt.plot(x, y, 'k--', linewidth=2)
        
        plt.text(
            0.05, 0.95,
            f'R² = {r_value**2:.3f}\np = {p_value:.3e}',
            transform=plt.gca().transAxes,
            bbox=dict(facecolor='white', alpha=0.7)
        )
    except Exception as e:
        print(f"Warning: Could not calculate regression line: {str(e)}")
    
    # Labels and title
    plt.xlabel('Change in Protein-Protein Contact (With Target Lipid - Without Target Lipid)')
    plt.ylabel('Change in Lipid-Protein Contact (With Target Lipid - Without Target Lipid)')
    plt.title('Impact of Target Lipid on Contact Distribution')
    plt.grid(True, alpha=0.3)
    
    # Add zero lines
    plt.axhline(y=0, color='k', linestyle='-', alpha=0.3)
    plt.axvline(x=0, color='k', linestyle='-', alpha=0.3)
    
    # Legend (if target lipid contact data available)
    if target_lipid_col in merged_df.columns:
        from matplotlib.lines import Line2D
        legend_elements = [
            Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=12, label='Strong Target Lipid Contact'),
            Line2D([0], [0], marker='o', color='w', markerfacecolor='orange', markersize=10, label='Weak Target Lipid Contact'),
            Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=8, label='No Target Lipid Contact')
        ]
        plt.legend(handles=legend_elements, loc='best')
    
    # Save
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "contact_change_scatter.png"), dpi=300)
    plt.savefig(os.path.join(output_dir, "contact_change_scatter.ps"))
    plt.close()
    
    # ----- 3. Per-protein ratio change heatmaps -----
    for protein in merged_df['protein'].unique():
        protein_df = merged_df[merged_df['protein'] == protein]
        
        # Skip if too few residues
        if len(protein_df) < 5:
            continue
        
        # Prepare heatmap data
        heatmap_df = protein_df.sort_values('residue')
        
        # Create pivot for plotting
        heatmap_data = pd.DataFrame({
            'residue': heatmap_df['residue'],
            'Ratio Difference': heatmap_df['ratio_diff'],
            'Protein Contact Difference': heatmap_df['protein_contact_diff'],
            'Lipid Contact Difference': heatmap_df['lipid_contact_diff']
        })
        
        if target_lipid_col in heatmap_df.columns:
            heatmap_data['Target Lipid Contact'] = heatmap_df[target_lipid_col]
        
        # Set residue ID as index
        heatmap_data = heatmap_data.set_index('residue')
        
        # Plot
        plt.figure(figsize=(12, 8))
        
        # Colormap
        cmap = sns.diverging_palette(220, 10, as_cmap=True)
        
        # Heatmap
        sns.heatmap(
            heatmap_data, 
            cmap=cmap,
            center=0,
            linewidths=0.5,
            annot=False,
            cbar_kws={'label': 'Difference'}
        )
        
        # Title and labels
        plt.title(f'Target Lipid Impact on Contact Distribution for {protein}')
        plt.ylabel('Residue ID')
        plt.tight_layout()
        
        # Save
        plt.savefig(os.path.join(output_dir, f"{protein}_difference_heatmap.png"), dpi=300)
        plt.savefig(os.path.join(output_dir, f"{protein}_difference_heatmap.ps"))
        plt.close()
    
    # ----- 4. Target lipid competitive inhibition mechanism analysis -----
    if target_lipid_col in merged_df.columns:
        # Compare residues with strong vs weak/no target lipid contact
        target_high = merged_df[merged_df[target_lipid_col] > 1.0]
        target_low = merged_df[merged_df[target_lipid_col] <= 0.04]
        
        # Statistical analysis
        print("\nStatistical Analysis of Target Lipid Impact:")
        
        try:
            # Protein contact change
            ttest_protein = stats.ttest_ind(
                target_high['protein_contact_diff'],
                target_low['protein_contact_diff'],
                equal_var=False
            )
            print(f"Protein Contact Change (target lipid high vs low): t={ttest_protein.statistic:.4f}, p={ttest_protein.pvalue:.4e}")
            
            # Lipid contact change
            ttest_lipid = stats.ttest_ind(
                target_high['lipid_contact_diff'],
                target_low['lipid_contact_diff'],
                equal_var=False
            )
            print(f"Lipid Contact Change (target lipid high vs low): t={ttest_lipid.statistic:.4f}, p={ttest_lipid.pvalue:.4e}")
            
            # Ratio change
            ttest_ratio = stats.ttest_ind(
                target_high['ratio_diff'],
                target_low['ratio_diff'],
                equal_var=False
            )
            print(f"Ratio Difference (target lipid high vs low): t={ttest_ratio.statistic:.4f}, p={ttest_ratio.pvalue:.4e}")
            
            # Violin plots to visualize target lipid impact
            plt.figure(figsize=(12, 10))
            
            # Prepare data
            plot_data = pd.DataFrame({
                'Target Lipid Contact Level': ['High'] * len(target_high) + ['Low/None'] * len(target_low),
                'Protein Contact Change': list(target_high['protein_contact_diff']) + list(target_low['protein_contact_diff']),
                'Lipid Contact Change': list(target_high['lipid_contact_diff']) + list(target_low['lipid_contact_diff']),
                'Ratio Difference': list(target_high['ratio_diff']) + list(target_low['ratio_diff'])
            })
            
            # Subplots
            fig, axes = plt.subplots(1, 3, figsize=(18, 6))
            
            # Protein contact change
            sns.violinplot(x='Target Lipid Contact Level', y='Protein Contact Change', data=plot_data, ax=axes[0])
            axes[0].set_title(f'Protein Contact Change\np={ttest_protein.pvalue:.2e}')
            axes[0].axhline(y=0, color='r', linestyle='--')
            
            # Lipid contact change
            sns.violinplot(x='Target Lipid Contact Level', y='Lipid Contact Change', data=plot_data, ax=axes[1])
            axes[1].set_title(f'Lipid Contact Change\np={ttest_lipid.pvalue:.2e}')
            axes[1].axhline(y=0, color='r', linestyle='--')
            
            # Ratio change
            sns.violinplot(x='Target Lipid Contact Level', y='Ratio Difference', data=plot_data, ax=axes[2])
            axes[2].set_title(f'Ratio Difference\np={ttest_ratio.pvalue:.2e}')
            axes[2].axhline(y=0, color='r', linestyle='--')
            
            # Overall title
            plt.suptitle('Impact of Target Lipid on Protein-Protein and Lipid-Protein Interactions', fontsize=16)
            plt.tight_layout()
            
            # Save
            plt.savefig(os.path.join(output_dir, "target_lipid_impact_violin.png"), dpi=300)
            plt.savefig(os.path.join(output_dir, "target_lipid_impact_violin.ps"))
            plt.close()
            
        except Exception as e:
            print(f"Warning: Could not perform statistical analysis: {str(e)}")
    
    # ----- 5. Generate essential competition plots -----
    print("\n===== Generating Essential Competition Plots =====")
    
    try:
        # 1. Energy comparison plot
        plot_energy_comparison(with_lipid_df, without_lipid_df, output_dir)
        
        # 2. Target lipid vs protein interface competition for each protein
        for protein in with_lipid_df['protein'].unique():
            if protein in without_lipid_df['protein'].unique():
                plot_target_lipid_protein_interface_competition(
                    with_lipid_df, without_lipid_df, output_dir, protein
                )
        
        # 3. Posterior distribution of target lipid effect on protein contacts
        plot_posterior_distribution_target_lipid_on_protein_contacts(
            with_lipid_df, without_lipid_df, output_dir
        )
        
        # 4. Bayesian regression: effect of target lipid on protein-protein contacts
        plot_bayesian_regression_target_lipid_effect(
            with_lipid_df, without_lipid_df, output_dir
        )
        
        # 5. protein_model_effect_simple plot (exact reproduction from bayesian_models.py)
        plot_protein_model_effect_simple(
            with_lipid_df, without_lipid_df, output_dir
        )
        
        print("All essential competition plots completed!")
        
    except Exception as e:
        print(f"Error generating essential competition plots: {str(e)}")
        import traceback
        traceback.print_exc()
    
    print(f"All comparison plots generated and saved to {output_dir}")