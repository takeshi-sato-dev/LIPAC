"""
Exploratory data analysis functions for BayesianLipidAnalysis
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats


def exploratory_analysis(combined_df, matched_df, output_dir):
    """Exploratory data analysis"""
    print("\n===== Exploratory Data Analysis =====")
    
    # Create output directory
    figures_dir = os.path.join(output_dir, "figures")
    os.makedirs(figures_dir, exist_ok=True)
    
    # 1. Compare target lipid-containing system with target lipid-free system
    if combined_df is not None:
        # Compare protein-protein contact distribution between target lipid-present and target lipid-absent systems
        plt.figure(figsize=(10, 6))
        sns.histplot(data=combined_df, x='protein_contact', hue='has_target_lipid', kde=True, bins=30, alpha=0.5)
        plt.title('Distribution of Protein-Protein Contacts')
        plt.xlabel('Protein-Protein Contact Frequency')
        plt.savefig(os.path.join(figures_dir, 'protein_contact_dist.png'), dpi=300)
        plt.close()
        
        # Compare lipid-protein contact distribution between target lipid-present and target lipid-absent systems
        plt.figure(figsize=(10, 6))
        sns.histplot(data=combined_df, x='lipid_contact', hue='has_target_lipid', kde=True, bins=30, alpha=0.5)
        plt.title('Distribution of Lipid-Protein Contacts')
        plt.xlabel('Lipid-Protein Contact Frequency')
        plt.savefig(os.path.join(figures_dir, 'lipid_contact_dist.png'), dpi=300)
        plt.close()
        
        # Compare ratio distribution between target lipid-present and target lipid-absent systems
        plt.figure(figsize=(10, 6))
        sns.histplot(data=combined_df, x='ratio', hue='has_target_lipid', kde=True, bins=30, alpha=0.5)
        plt.title('Distribution of Lipid-to-Protein Contact Ratio')
        plt.xlabel('Lipid-to-Protein Contact Ratio')
        plt.savefig(os.path.join(figures_dir, 'ratio_dist.png'), dpi=300)
        plt.close()
        
        # Statistical testing
        with_lipid_subset = combined_df[combined_df['has_target_lipid']]
        without_lipid_subset = combined_df[~combined_df['has_target_lipid']]
        
        # Protein-protein contact test
        ttest_protein = stats.ttest_ind(
            with_lipid_subset['protein_contact'],
            without_lipid_subset['protein_contact'],
            equal_var=False
        )
        print(f"Protein-protein contact difference (with target lipid vs without): t={ttest_protein.statistic:.4f}, p={ttest_protein.pvalue:.4e}")
        
        # Lipid-protein contact test
        ttest_lipid = stats.ttest_ind(
            with_lipid_subset['lipid_contact'],
            without_lipid_subset['lipid_contact'],
            equal_var=False
        )
        print(f"Lipid-protein contact difference (with target lipid vs without): t={ttest_lipid.statistic:.4f}, p={ttest_lipid.pvalue:.4e}")
        
        # Ratio test
        ttest_ratio = stats.ttest_ind(
            with_lipid_subset['ratio'],
            without_lipid_subset['ratio'],
            equal_var=False
        )
        print(f"Ratio difference (with target lipid vs without): t={ttest_ratio.statistic:.4f}, p={ttest_ratio.pvalue:.4e}")
    
    # 2. Relationship between target lipid contact and contact differences
    if matched_df is not None and 'target_lipid_contact' in matched_df.columns:
        # Scatter plot: Target lipid contact vs protein contact difference
        plt.figure(figsize=(10, 6))
        sns.regplot(data=matched_df, x='target_lipid_contact', y='protein_contact_diff', scatter_kws={'alpha': 0.5})
        plt.title('Target Lipid Contact vs Change in Protein-Protein Contact')
        plt.xlabel('Target Lipid Contact Frequency')
        plt.ylabel('Change in Protein-Protein Contact (With Target Lipid - Without Target Lipid)')
        plt.axhline(y=0, color='r', linestyle='--')
        
        # Calculate and display correlation coefficient
        corr = matched_df[['target_lipid_contact', 'protein_contact_diff']].corr().iloc[0, 1]
        plt.text(0.05, 0.95, f'Correlation: {corr:.3f}', transform=plt.gca().transAxes, 
                bbox=dict(facecolor='white', alpha=0.7))
        
        plt.savefig(os.path.join(figures_dir, 'target_lipid_vs_protein_contact_diff.png'), dpi=300)
        plt.close()
        
        # Scatter plot: Target lipid contact vs lipid contact difference
        plt.figure(figsize=(10, 6))
        sns.regplot(data=matched_df, x='target_lipid_contact', y='lipid_contact_diff', scatter_kws={'alpha': 0.5})
        plt.title('Target Lipid Contact vs Change in Lipid-Protein Contact')
        plt.xlabel('Target Lipid Contact Frequency')
        plt.ylabel('Change in Lipid-Protein Contact (With Target Lipid - Without Target Lipid)')
        plt.axhline(y=0, color='r', linestyle='--')
        
        # Calculate and display correlation coefficient
        corr = matched_df[['target_lipid_contact', 'lipid_contact_diff']].corr().iloc[0, 1]
        plt.text(0.05, 0.95, f'Correlation: {corr:.3f}', transform=plt.gca().transAxes,
                bbox=dict(facecolor='white', alpha=0.7))
        
        plt.savefig(os.path.join(figures_dir, 'target_lipid_vs_lipid_contact_diff.png'), dpi=300)
        plt.close()
        
        # Scatter plot: Target lipid contact vs ratio difference
        plt.figure(figsize=(10, 6))
        sns.regplot(data=matched_df, x='target_lipid_contact', y='ratio_diff', scatter_kws={'alpha': 0.5})
        plt.title('Target Lipid Contact vs Change in Contact Ratio')
        plt.xlabel('Target Lipid Contact Frequency')
        plt.ylabel('Change in Ratio (With Target Lipid - Without Target Lipid)')
        plt.axhline(y=0, color='r', linestyle='--')
        
        # Calculate and display correlation coefficient
        corr = matched_df[['target_lipid_contact', 'ratio_diff']].corr().iloc[0, 1]
        plt.text(0.05, 0.95, f'Correlation: {corr:.3f}', transform=plt.gca().transAxes,
                bbox=dict(facecolor='white', alpha=0.7))
        
        plt.savefig(os.path.join(figures_dir, 'target_lipid_vs_ratio_diff.png'), dpi=300)
        plt.close()
        
        # Compare residues with strong target lipid contact to those with weak/no contact
        target_lipid_high = matched_df[matched_df['target_lipid_contact'] > 0.5]
        target_lipid_low = matched_df[matched_df['target_lipid_contact'] <= 0.1]
        
        print("\nComparison by target lipid contact strength:")
        print(f"Strong target lipid contact: {len(target_lipid_high)} residues")
        print(f"Weak/no target lipid contact: {len(target_lipid_low)} residues")
        
        # Compare protein contact differences
        ttest_protein_diff = stats.ttest_ind(
            target_lipid_high['protein_contact_diff'],
            target_lipid_low['protein_contact_diff'],
            equal_var=False
        )
        print(f"Protein contact difference (strong vs weak/no target lipid): t={ttest_protein_diff.statistic:.4f}, p={ttest_protein_diff.pvalue:.4e}")
        print(f"  Strong target lipid mean: {target_lipid_high['protein_contact_diff'].mean():.4f}")
        print(f"  Weak/no target lipid mean: {target_lipid_low['protein_contact_diff'].mean():.4f}")
        
        # Violin plot visualization of comparison
        plt.figure(figsize=(12, 6))
        plot_data = pd.concat([
            pd.DataFrame({'Target Lipid Contact': 'High', 'Protein Contact Diff': target_lipid_high['protein_contact_diff']}),
            pd.DataFrame({'Target Lipid Contact': 'Low/None', 'Protein Contact Diff': target_lipid_low['protein_contact_diff']})
        ])
        sns.violinplot(data=plot_data, x='Target Lipid Contact', y='Protein Contact Diff')
        plt.title(f'Effect of Target Lipid on Protein-Protein Contacts\np = {ttest_protein_diff.pvalue:.2e}')
        plt.axhline(y=0, color='r', linestyle='--')
        plt.savefig(os.path.join(figures_dir, 'target_lipid_effect_violin.png'), dpi=300)
        plt.close()
        
        # Save results to text file
        with open(os.path.join(output_dir, 'exploratory_analysis_results.txt'), 'w') as f:
            f.write("===== Target Lipid-Protein Interaction Exploratory Analysis =====\n\n")
            f.write(f"With target lipid: {len(with_lipid_subset)} residues\n")
            f.write(f"Without target lipid: {len(without_lipid_subset)} residues\n\n")
            
            f.write("--- Statistical Test Results ---\n")
            f.write(f"Protein-protein contact difference (with target lipid vs without): t={ttest_protein.statistic:.4f}, p={ttest_protein.pvalue:.4e}\n")
            f.write(f"Lipid-protein contact difference (with target lipid vs without): t={ttest_lipid.statistic:.4f}, p={ttest_lipid.pvalue:.4e}\n")
            f.write(f"Ratio difference (with target lipid vs without): t={ttest_ratio.statistic:.4f}, p={ttest_ratio.pvalue:.4e}\n\n")
            
            f.write("--- Comparison by Target Lipid Contact Strength ---\n")
            f.write(f"Strong target lipid contact: {len(target_lipid_high)} residues\n")
            f.write(f"Weak/no target lipid contact: {len(target_lipid_low)} residues\n")
            f.write(f"Protein contact difference (strong vs weak/no target lipid): t={ttest_protein_diff.statistic:.4f}, p={ttest_protein_diff.pvalue:.4e}\n")
            f.write(f"  Strong target lipid mean: {target_lipid_high['protein_contact_diff'].mean():.4f}\n")
            f.write(f"  Weak/no target lipid mean: {target_lipid_low['protein_contact_diff'].mean():.4f}\n")
        
        # 173行目（間違い）
        print(f"Exploratory analysis results saved to {os.path.join(output_dir, 'exploratory_analysis_results.txt')}")