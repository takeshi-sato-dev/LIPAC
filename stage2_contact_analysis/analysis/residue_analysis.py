"""
Residue-level analysis functions for BayesianLipidAnalysis

Author: Takeshi Sato, PhD
Kyoto Pharmaceutical University
2024
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats


def analyze_residue_contacts(with_gm3_df, no_gm3_df, output_dir):
    """Analyze GM3 contact and protein contact at residue level (data-driven approach)
    
    Parameters
    ----------
    with_gm3_df : pd.DataFrame
        Data with GM3
    no_gm3_df : pd.DataFrame
        Data without GM3
    output_dir : str
        Output directory
    
    Returns
    -------
    bool
        Success status
    """
    print("\n===== Residue-Level GM3/Protein Contact Analysis =====")
    
    # Create output directory
    residue_dir = os.path.join(output_dir, "residue_analysis")
    os.makedirs(residue_dir, exist_ok=True)
    
    # Check for empty dataframes
    if len(with_gm3_df) == 0 or len(no_gm3_df) == 0:
        print("Error: Dataframes are empty")
        return False
    
    # Display protein information in CSV
    proteins = with_gm3_df['protein'].unique()
    print(f"Proteins in CSV: {proteins}")
    
    # Extract actual protein pairs from data
    unique_pairs = []
    for protein in proteins:
        # Find all partners for this protein
        protein_data = with_gm3_df[with_gm3_df['protein'] == protein]
        if 'partner_protein' in protein_data.columns:
            partners = protein_data['partner_protein'].unique()
            for partner in partners:
                if partner != 'none' and partner != protein:
                    # Sort pair alphabetically to avoid duplicates
                    pair = tuple(sorted([protein, partner]))
                    if pair not in unique_pairs:
                        unique_pairs.append(pair)
    
    if not unique_pairs and 'partner_protein' not in with_gm3_df.columns:
        print("No partner protein information found. Performing single protein analysis.")
        # Analyze each protein individually
        for protein in proteins:
            analyze_single_protein(with_gm3_df, no_gm3_df, protein, residue_dir)
    else:
        # Analyze each pair
        print(f"Detected protein pairs:")
        for protein, partner in unique_pairs:
            print(f"  - {protein} and {partner}")
            analyze_protein_pair(with_gm3_df, no_gm3_df, protein, partner, residue_dir)
    
    # Quantitative analysis of GM3 binding site and interface overlap
    if len(unique_pairs) > 0:
        analyze_binding_interface_overlap(with_gm3_df, no_gm3_df, unique_pairs, residue_dir)
    
    print(f"Residue-level analysis results saved to {residue_dir}")
    return True


def analyze_single_protein(with_gm3_df, no_gm3_df, protein, output_dir):
    """Single protein residue analysis
    
    Parameters
    ----------
    with_gm3_df : pd.DataFrame
        Data with GM3
    no_gm3_df : pd.DataFrame
        Data without GM3
    protein : str
        Protein name
    output_dir : str
        Output directory
    """
    print(f"\nAnalysis: {protein} (single protein)")
    
    # Extract protein data
    with_gm3_data = with_gm3_df[with_gm3_df['protein'] == protein]
    no_gm3_data = no_gm3_df[no_gm3_df['protein'] == protein]
    
    if len(with_gm3_data) == 0 or len(no_gm3_data) == 0:
        print(f"  Warning: Insufficient data for {protein}")
        return
    
    # Correlation analysis of GM3 contact and protein contact
    merged_data = pd.merge(
        with_gm3_data[['residue', 'protein_contact', 'DPG3_contact']],
        no_gm3_data[['residue', 'protein_contact']],
        on='residue', suffixes=('_with_gm3', '_no_gm3')
    )
    
    # Calculate differences
    merged_data['protein_contact_diff'] = (
        merged_data['protein_contact_with_gm3'] - merged_data['protein_contact_no_gm3']
    )
    
    # Correlation between GM3 contact and protein contact change
    corr, p_value = stats.pearsonr(
        merged_data['DPG3_contact'], 
        merged_data['protein_contact_diff']
    )
    
    print(f"  Correlation between GM3 contact and protein contact change: r={corr:.3f}, p={p_value:.3e}")
    
    # Identify statistically significant GM3 binding residues
    gm3_threshold = np.percentile(merged_data['DPG3_contact'], 75)
    high_gm3 = merged_data[merged_data['DPG3_contact'] >= gm3_threshold]
    low_gm3 = merged_data[merged_data['DPG3_contact'] < gm3_threshold]
    
    # T-test comparing protein contact change between high GM3 and low GM3 binding residues
    if len(high_gm3) > 0 and len(low_gm3) > 0:
        t_stat, p_val = stats.ttest_ind(
            high_gm3['protein_contact_diff'],
            low_gm3['protein_contact_diff'],
            equal_var=False
        )
        print(f"  Protein contact difference (high GM3 vs low GM3): t={t_stat:.3f}, p={p_val:.3e}")
        print(f"  High GM3 mean change: {high_gm3['protein_contact_diff'].mean():.3f}")
        print(f"  Low GM3 mean change: {low_gm3['protein_contact_diff'].mean():.3f}")
    
    # Visualization: Correlation between GM3 contact and protein contact change
    plt.figure(figsize=(10, 8))
    plt.scatter(
        merged_data['DPG3_contact'],
        merged_data['protein_contact_diff'],
        alpha=0.7, s=50
    )
    
    # Residue labels
    for i, row in merged_data.iterrows():
        if row['DPG3_contact'] >= gm3_threshold or abs(row['protein_contact_diff']) > 0.1:
            plt.annotate(
                f"{row['residue']}", 
                (row['DPG3_contact'], row['protein_contact_diff']),
                xytext=(5, 5), textcoords='offset points',
                fontsize=9
            )
    
    # Regression line
    if len(merged_data) > 2:  # Enough data points for regression
        slope, intercept, r_value, p_value, std_err = stats.linregress(
            merged_data['DPG3_contact'], 
            merged_data['protein_contact_diff']
        )
        x_range = np.linspace(0, merged_data['DPG3_contact'].max() * 1.1, 100)
        plt.plot(
            x_range, 
            intercept + slope * x_range, 
            'r--', 
            label=f'y = {slope:.3f}x + {intercept:.3f} (r={r_value:.3f}, p={p_value:.3e})'
        )
    
    # Graph decoration
    plt.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    plt.xlabel('GM3 Contact Frequency', fontsize=12)
    plt.ylabel('Protein Contact Change (With GM3 - Without GM3)', fontsize=12)
    plt.title(f'{protein}: Relationship Between GM3 Binding and Protein Contact Change', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{protein}_gm3_effect.png'), dpi=300)
    plt.savefig(os.path.join(output_dir, f'{protein}_gm3_effect.ps'))
    plt.close()
    
    # Save results
    merged_data.to_csv(os.path.join(output_dir, f'{protein}_residue_analysis.csv'))


def analyze_protein_pair(with_gm3_df, no_gm3_df, protein, partner, output_dir):
    """Protein pair residue analysis
    
    Parameters
    ----------
    with_gm3_df : pd.DataFrame
        Data with GM3
    no_gm3_df : pd.DataFrame
        Data without GM3
    protein : str
        Protein name
    partner : str
        Partner protein name
    output_dir : str
        Output directory
    """
    from visualization.plots import plot_residue_profile
    
    print(f"\nAnalysis: {protein}-{partner} pair")
    
    # Extract pair data
    with_gm3_data = with_gm3_df[(with_gm3_df['protein'] == protein) & 
                               (with_gm3_df['partner_protein'] == partner)]
    
    no_gm3_data = no_gm3_df[(no_gm3_df['protein'] == protein) & 
                           (no_gm3_df['partner_protein'] == partner)]
    
    # Check reverse pair
    if len(with_gm3_data) == 0:
        with_gm3_data = with_gm3_df[(with_gm3_df['protein'] == partner) & 
                                  (with_gm3_df['partner_protein'] == protein)]
    
    if len(no_gm3_data) == 0:
        no_gm3_data = no_gm3_df[(no_gm3_df['protein'] == partner) & 
                               (no_gm3_df['partner_protein'] == protein)]
    
    if len(with_gm3_data) == 0 or len(no_gm3_data) == 0:
        print(f"  Warning: Insufficient data for {protein}-{partner} pair")
        return
    
    # Merge residues and calculate differences
    merged_data = pd.merge(
        with_gm3_data[['residue', 'protein_contact', 'DPG3_contact', 'lipid_contact']],
        no_gm3_data[['residue', 'protein_contact', 'lipid_contact']],
        on='residue', suffixes=('_with_gm3', '_no_gm3')
    )
    
    # Calculate differences
    merged_data['protein_contact_diff'] = (
        merged_data['protein_contact_with_gm3'] - merged_data['protein_contact_no_gm3']
    )
    merged_data['lipid_contact_diff'] = (
        merged_data['lipid_contact_with_gm3'] - merged_data['lipid_contact_no_gm3']
    )
    
    # Sort residues by ID
    merged_data = merged_data.sort_values('residue')
    
    # Identify interface residues and GM3 binding residues
    protein_threshold = np.percentile(merged_data['protein_contact_no_gm3'], 75)
    gm3_threshold = np.percentile(merged_data['DPG3_contact'], 75)
    
    merged_data['is_interface'] = merged_data['protein_contact_no_gm3'] >= protein_threshold
    merged_data['is_gm3_binding'] = merged_data['DPG3_contact'] >= gm3_threshold
    merged_data['is_competition'] = (
        merged_data['is_interface'] & 
        merged_data['is_gm3_binding'] & 
        (merged_data['protein_contact_diff'] < 0)
    )
    
    # Calculate number of competitive binding sites
    n_competition = merged_data['is_competition'].sum()
    n_interface = merged_data['is_interface'].sum()
    n_gm3_binding = merged_data['is_gm3_binding'].sum()
    
    print(f"  Interface residues: {n_interface}")
    print(f"  GM3 binding residues: {n_gm3_binding}")
    print(f"  Competitive binding sites: {n_competition}")
    
    # Correlation between GM3 contact and protein contact change
    corr, p_value = stats.pearsonr(
        merged_data['DPG3_contact'], 
        merged_data['protein_contact_diff']
    )
    print(f"  Correlation between GM3 contact and protein contact change: r={corr:.3f}, p={p_value:.3e}")
    
    # Relationship between interface residues and GM3 binding
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Scatter plot: GM3 contact vs protein contact change
    scatter = ax.scatter(
        merged_data['DPG3_contact'],
        merged_data['protein_contact_diff'],
        c=merged_data['is_interface'].map({True: 'blue', False: 'gray'}),
        alpha=0.7, s=70
    )
    
    # Highlight competitive binding sites
    competition_sites = merged_data[merged_data['is_competition']]
    if len(competition_sites) > 0:
        ax.scatter(
            competition_sites['DPG3_contact'],
            competition_sites['protein_contact_diff'],
            color='red', s=100, alpha=0.8, edgecolor='black',
            marker='o', label='Competitive Binding Sites'
        )
    
    # Residue labels
    for i, row in merged_data.iterrows():
        if row['is_competition'] or (row['is_interface'] and row['DPG3_contact'] > 0.2):
            ax.annotate(
                f"{row['residue']}", 
                (row['DPG3_contact'], row['protein_contact_diff']),
                xytext=(5, 5), textcoords='offset points',
                fontsize=10, fontweight='bold'
            )
    
    # Reference lines
    ax.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    ax.axvline(x=gm3_threshold, color='r', linestyle='--', alpha=0.5, 
               label=f'GM3 Threshold ({gm3_threshold:.2f})')
    
    # Regression line
    if len(merged_data) > 2:
        slope, intercept, r_value, p_value, std_err = stats.linregress(
            merged_data['DPG3_contact'], 
            merged_data['protein_contact_diff']
        )
        x_range = np.linspace(0, merged_data['DPG3_contact'].max() * 1.1, 100)
        ax.plot(
            x_range, 
            intercept + slope * x_range, 
            'r-', linewidth=2, 
            label=f'Regression: y = {slope:.3f}x + {intercept:.3f}'
        )
    
    # Graph decoration
    ax.set_xlabel('GM3 Contact Frequency', fontsize=12)
    ax.set_ylabel('Protein Contact Change (With GM3 - Without GM3)', fontsize=12)
    ax.set_title(f'{protein}-{partner}: Relationship Between GM3 Binding and Protein Contact Change', fontsize=14)
    
    # Custom legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', 
               markersize=10, label='Interface Residues'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', 
               markersize=10, label='Non-Interface Residues'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='red', 
               markersize=10, label='Competitive Binding Sites')
    ]
    ax.legend(handles=legend_elements, loc='best')
    
    # Statistical information
    stat_text = (
        f"Correlation: {corr:.3f} (p={p_value:.3e})\n"
        f"Interface residues: {n_interface} ({n_interface/len(merged_data)*100:.1f}%)\n"
        f"GM3 binding residues: {n_gm3_binding} ({n_gm3_binding/len(merged_data)*100:.1f}%)\n"
        f"Competitive binding sites: {n_competition} ({n_competition/len(merged_data)*100:.1f}%)"
    )
    ax.text(0.02, 0.98, stat_text, transform=ax.transAxes, 
            va='top', ha='left', bbox=dict(facecolor='white', alpha=0.7))
    
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{protein}_{partner}_competition_analysis.png'), dpi=300)
    plt.savefig(os.path.join(output_dir, f'{protein}_{partner}_competition_analysis.ps'))
    plt.close()
    
    # Create residue profile
    plot_residue_profile(merged_data, protein, partner, output_dir)
    
    # Save results
    merged_data.to_csv(os.path.join(output_dir, f'{protein}_{partner}_residue_analysis.csv'))
    
    # Save competitive binding site details
    if n_competition > 0:
        competition_sites.to_csv(
            os.path.join(output_dir, f'{protein}_{partner}_competition_sites.csv')
        )


def analyze_binding_interface_overlap(with_gm3_df, no_gm3_df, protein_pairs, output_dir):
    """Analyze overlap between GM3 binding sites and interface regions
    
    Parameters
    ----------
    with_gm3_df : pd.DataFrame
        Data with GM3
    no_gm3_df : pd.DataFrame
        Data without GM3
    protein_pairs : list
        List of protein pairs
    output_dir : str
        Output directory
    """
    from visualization.plots import plot_interface_enrichment, plot_residue_distribution_pie
    
    print("\n=== Analysis of GM3 Binding Site and Interface Region Overlap ===")
    
    # Compile data for all protein pairs
    all_data = []
    
    for protein, partner in protein_pairs:
        # Extract pair data
        with_gm3_data = with_gm3_df[(with_gm3_df['protein'] == protein) & 
                                  (with_gm3_df['partner_protein'] == partner)]
        
        no_gm3_data = no_gm3_df[(no_gm3_df['protein'] == protein) & 
                               (no_gm3_df['partner_protein'] == partner)]
        
        # Check reverse pair
        if len(with_gm3_data) == 0:
            with_gm3_data = with_gm3_df[(with_gm3_df['protein'] == partner) & 
                                      (with_gm3_df['partner_protein'] == protein)]
        
        if len(no_gm3_data) == 0:
            no_gm3_data = no_gm3_df[(no_gm3_df['protein'] == partner) & 
                                   (no_gm3_df['partner_protein'] == protein)]
        
        if len(with_gm3_data) == 0 or len(no_gm3_data) == 0:
            continue
        
        # Merge residues and calculate differences
        merged_data = pd.merge(
            with_gm3_data[['residue', 'protein', 'protein_contact', 'DPG3_contact']],
            no_gm3_data[['residue', 'protein_contact']],
            on=['residue', 'protein'], suffixes=('_with_gm3', '_no_gm3')
        )
        
        # Calculate differences
        merged_data['protein_contact_diff'] = (
            merged_data['protein_contact_with_gm3'] - merged_data['protein_contact_no_gm3']
        )
        
        # Identify interface residues and GM3 binding residues
        protein_threshold = np.percentile(merged_data['protein_contact_no_gm3'], 75)
        gm3_threshold = np.percentile(merged_data['DPG3_contact'], 75)
        
        merged_data['is_interface'] = merged_data['protein_contact_no_gm3'] >= protein_threshold
        merged_data['is_gm3_binding'] = merged_data['DPG3_contact'] >= gm3_threshold
        merged_data['is_competition'] = (
            merged_data['is_interface'] & 
            merged_data['is_gm3_binding'] & 
            (merged_data['protein_contact_diff'] < 0)
        )
        
        # Add pair information
        merged_data['pair'] = f"{protein}-{partner}"
        
        # Accumulate results
        all_data.append(merged_data)
    
    # Combine all data
    if not all_data:
        print("  No analyzable pair data available")
        return
    
    all_pairs_data = pd.concat(all_data, ignore_index=True)
    
    # Calculate statistics
    pair_stats = []
    for pair, group in all_pairs_data.groupby('pair'):
        n_total = len(group)
        n_interface = group['is_interface'].sum()
        n_gm3_binding = group['is_gm3_binding'].sum()
        n_competition = group['is_competition'].sum()
        
        # Expected overlap (probabilistic)
        expected_overlap = (n_interface / n_total) * (n_gm3_binding / n_total) * n_total
        
        # Ratio of observed to expected (enrichment)
        enrichment = n_competition / expected_overlap if expected_overlap > 0 else float('nan')
        
        # Fisher's exact test
        table = np.array([
            [n_competition, n_interface - n_competition],
            [n_gm3_binding - n_competition, n_total - n_interface - n_gm3_binding + n_competition]
        ])
        
        oddsratio, p_value = stats.fisher_exact(table)
        
        # Average effect size
        competition_effect = group.loc[group['is_competition'], 'protein_contact_diff'].mean() if n_competition > 0 else 0
        
        pair_stats.append({
            'pair': pair,
            'total_residues': n_total,
            'interface_residues': n_interface,
            'gm3_binding_residues': n_gm3_binding,
            'competition_sites': n_competition,
            'expected_overlap': expected_overlap,
            'enrichment': enrichment,
            'odds_ratio': oddsratio,
            'p_value': p_value,
            'avg_competition_effect': competition_effect
        })
    
    pair_stats_df = pd.DataFrame(pair_stats)
    
    # Display results
    print("  Competitive binding site analysis by protein pair:")
    for _, row in pair_stats_df.iterrows():
        print(f"  {row['pair']}:")
        print(f"    Competitive binding sites: {int(row['competition_sites'])} (out of {int(row['total_residues'])} residues)")
        print(f"    Interface residues: {int(row['interface_residues'])}")
        print(f"    GM3 binding residues: {int(row['gm3_binding_residues'])}")
        print(f"    Enrichment: {row['enrichment']:.2f}x")
        print(f"    Statistical significance: p={row['p_value']:.3e}")
        print(f"    Average effect size: {row['avg_competition_effect']:.3f}")
    
    # Create plots
    plot_interface_enrichment(pair_stats_df, output_dir)
    plot_residue_distribution_pie(pair_stats_df, output_dir)
    
    # Save results
    pair_stats_df.to_csv(os.path.join(output_dir, 'gm3_interface_overlap_statistics.csv'))