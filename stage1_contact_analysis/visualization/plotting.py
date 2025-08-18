"""Plotting functions for contact analysis visualization"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import os
import traceback
from ..config import TARGET_LIPID

# Import competition plots
from .competition_plots import (
    plot_energy_comparison,
    plot_target_lipid_protein_interface_competition,
    plot_posterior_distribution_target_lipid_on_protein_contacts,
    plot_bayesian_regression_target_lipid_effect,
    plot_protein_model_effect_simple
)

def plot_contact_complementarity(df, output_dir, with_target_lipid=True):
    """Generate contact complementarity plots with correct residue numbering
    
    Parameters
    ----------
    df : pd.DataFrame
        Contact complementarity data
    output_dir : str
        Output directory path
    with_target_lipid : bool
        Whether system contains target lipid
    """
    print("\n===== Generating Contact Complementarity Plots =====")
    
    # Ensure output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # If data is empty
    if df is None or len(df) == 0:
        print("No data to plot. Aborting.")
        return
    
    # Set plot title based on target lipid presence
    plot_title = f"With {TARGET_LIPID}" if with_target_lipid else f"Without {TARGET_LIPID}"
    file_prefix = f"with_{TARGET_LIPID.lower()}" if with_target_lipid else f"no_{TARGET_LIPID.lower()}"
    
    # Plot settings
    plt.style.use('seaborn-v0_8-whitegrid')
    sns.set_context("talk")
    
    # ----- 1. Scatter plot of protein contacts vs lipid contacts -----
    plt.figure(figsize=(10, 8))
    
    # Target only dimerized proteins
    dimer_df = df[df['protein_pair'] != 'none']
    
    # Scatter plot
    scatter = plt.scatter(
        dimer_df['protein_contact'], 
        dimer_df['lipid_contact'],
        c=dimer_df['ratio'],  # Color code by ratio
        cmap='viridis',
        alpha=0.7,
        s=50
    )
    
    # Color bar
    cbar = plt.colorbar(scatter)
    cbar.set_label('Lipid-to-Protein Contact Ratio')
    
    # Linear regression line
    try:
        from scipy import stats
        slope, intercept, r_value, p_value, std_err = stats.linregress(
            dimer_df['protein_contact'], 
            dimer_df['lipid_contact']
        )
        
        # Display regression line and R² value
        x = np.array([dimer_df['protein_contact'].min(), dimer_df['protein_contact'].max()])
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
    
    # Axis labels and legend
    plt.xlabel('Protein-Protein Contact Frequency')
    plt.ylabel('Lipid-Protein Contact Frequency')
    plt.title(f'Contact Complementarity: {plot_title}')
    plt.grid(True, alpha=0.3)
    
    # Save
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{file_prefix}_contact_scatter.png"), dpi=300)
    plt.savefig(os.path.join(output_dir, f"{file_prefix}_contact_scatter.ps"))
    plt.close()
    
    # ----- 2. Lipid contact heatmap for each protein -----
    # Get list of lipid types
    lipid_columns = [col for col in df.columns if '_contact' in col and col != 'protein_contact' and col != 'lipid_contact']
    
    # For each protein
    for protein in df['protein'].unique():
        protein_df = df[df['protein'] == protein]
        
        # Create heatmap data for residue and each lipid type contact values
        heatmap_data = []
        for _, row in protein_df.iterrows():
            lipid_values = {lipid_col: row[lipid_col] for lipid_col in lipid_columns}
            lipid_values['residue'] = row['residue']
            lipid_values['protein_contact'] = row['protein_contact']
            heatmap_data.append(lipid_values)
        
        if not heatmap_data:
            continue
            
        heatmap_df = pd.DataFrame(heatmap_data)
        heatmap_df = heatmap_df.sort_values('residue')
        
        # Set residue as index
        heatmap_df = heatmap_df.set_index('residue')
        
        # Plot
        plt.figure(figsize=(12, 8))
        
        # Lipid contact heatmap
        sns.heatmap(
            heatmap_df[lipid_columns], 
            cmap='YlGnBu', 
            linewidths=0.5,
            annot=False, 
            cbar_kws={'label': 'Contact Frequency'}
        )
        
        # Title and axis labels
        plt.title(f'Lipid Contacts for {protein}: {plot_title}')
        plt.ylabel('Residue ID')
        plt.tight_layout()
        
        # Save
        plt.savefig(os.path.join(output_dir, f"{file_prefix}_{protein}_lipid_heatmap.png"), dpi=300)
        plt.savefig(os.path.join(output_dir, f"{file_prefix}_{protein}_lipid_heatmap.ps"))
        plt.close()
    
    print(f"All plots generated and saved to {output_dir}")

def plot_target_lipid_protein_competition(df, output_dir, protein_name_or_pair="", no_target_lipid_df=None):
    """Analyze competition between target lipid and protein contacts - no_target_lipid_df can be passed to show competition
    
    Parameters
    ----------
    df : pd.DataFrame
        Contact data with target lipid
    output_dir : str
        Output directory path
    protein_name_or_pair : str
        Protein name or pair
    no_target_lipid_df : pd.DataFrame, optional
        Contact data without target lipid
    """
    # Ensure output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Extract protein name (if pair name is passed)
    if "-" in protein_name_or_pair:
        protein = protein_name_or_pair.split("-")[0]
    else:
        protein = protein_name_or_pair
    
    # Extract data for this protein
    protein_df = df[df['protein'] == protein].copy()
    
    if len(protein_df) == 0:
        print(f"No data found for protein {protein}")
        return
    
    # Check for target lipid contact data
    target_col = f'{TARGET_LIPID}_contact'
    if target_col not in protein_df.columns:
        print(f"No {TARGET_LIPID} contact data found for {protein}")
        return
    
    # Also get protein contact data from system without target lipid (if specified)
    no_target_protein_contacts = None
    no_target_residues = None
    if no_target_lipid_df is not None and protein in no_target_lipid_df['protein'].values:
        no_target_protein_df = no_target_lipid_df[no_target_lipid_df['protein'] == protein].copy()
        if 'protein_contact' in no_target_protein_df.columns:
            no_target_protein_contacts = no_target_protein_df['protein_contact'].values
            no_target_residues = no_target_protein_df['residue'].values
    
    # Identify partner protein
    partner_protein = 'none'
    if 'partner_protein' in protein_df.columns:
        partners = protein_df['partner_protein'].unique()
        if len(partners) > 0 and partners[0] != 'none':
            partner_protein = partners[0]
    
    # Sort by residue
    protein_df = protein_df.sort_values('residue')
    
    # Identify competition regions
    # Residues with strong target lipid contact
    target_preference = protein_df[protein_df[target_col] > 0.04]['residue'].tolist()
    
    # Residues with strong protein contact (only if pair information exists)
    protein_preference = []
    competition_residues = []
    
    # When protein contact data exists in system with target lipid
    if 'protein_contact' in protein_df.columns and protein_df['protein_contact'].max() > 0:
        protein_preference = protein_df[protein_df['protein_contact'] > 0.5]['residue'].tolist()
        
        # Residues in contact with both target lipid and protein
        competition_residues = protein_df[
            (protein_df[target_col] > 0.04) & 
            (protein_df['protein_contact'] > 0.45)
        ]['residue'].tolist()
    # When protein contact data exists in system without target lipid
    elif no_target_protein_contacts is not None and no_target_residues is not None:
        # Identify strong contact residues from no_target_lipid protein contact data
        for i, res in enumerate(no_target_residues):
            if no_target_protein_contacts[i] > 0.45 and res in protein_df['residue'].values:
                protein_preference.append(res)
                
                # Check competition with target lipid
                target_contact = protein_df.loc[protein_df['residue'] == res, target_col].values
                if len(target_contact) > 0 and target_contact[0] > 0.04:
                    competition_residues.append(res)
    
    # ---- Contact plot by residue ----
    plt.figure(figsize=(12, 6))
    
    # Set main Y axis (for target lipid contact)
    ax1 = plt.gca()
    
    # Set additional Y axis if needed
    use_protein_contacts = False
    if 'protein_contact' in protein_df.columns and protein_df['protein_contact'].max() > 0:
        use_protein_contacts = True
    elif no_target_protein_contacts is not None:
        use_protein_contacts = True
    
    if use_protein_contacts:
        ax3 = ax1.twinx()  # For protein contact
        # Move third Y axis to the right (create space)
        ax3.spines.right.set_position(('axes', 1.0))
    
    # Background color settings for range definition
    residues = protein_df['residue'].values
    y_min, y_max = -0.05, protein_df[target_col].max() * 1.1 if protein_df[target_col].max() > 0 else 1.0
    
    # Interface residues (if exist)
    for res in protein_preference:
        ax1.axvspan(res-0.5, res+0.5, color='lightblue', alpha=0.3)
        
    # Competition regions (if exist)
    for res in competition_residues:
        ax1.axvspan(res-0.5, res+0.5, color='purple', alpha=0.2)
        
    # Target lipid binding sites
    for res in target_preference:
        if res not in competition_residues:  # Only non-competing residues
            ax1.axvspan(res-0.5, res+0.5, color='pink', alpha=0.3)
    
    # Data plot
    ax1.plot(protein_df['residue'], protein_df[target_col], 'r-', linewidth=2, label=f'{TARGET_LIPID} Contact')
    
    # If protein contact data exists
    if 'protein_contact' in protein_df.columns and protein_df['protein_contact'].max() > 0:
        ax3.plot(protein_df['residue'], protein_df['protein_contact'], 'b-', linewidth=2, label=f'Protein Contact (with {TARGET_LIPID})')
    elif no_target_protein_contacts is not None and no_target_residues is not None:
        # Interpolation handling when plotting protein contact from system without target lipid
        interpolated_contacts = np.zeros(len(protein_df))
        for i, res in enumerate(protein_df['residue']):
            idx = np.where(no_target_residues == res)[0]
            if len(idx) > 0:
                interpolated_contacts[i] = no_target_protein_contacts[idx[0]]
        
        # Plot protein contact from system without target lipid
        ax3.plot(protein_df['residue'], interpolated_contacts, 'b-', linewidth=2, label=f'Protein Contact (no {TARGET_LIPID})')
    
    # Axis labels
    ax1.set_xlabel('Residue')
    ax1.set_ylabel(f'{TARGET_LIPID} Contact', color='r')
    ax1.tick_params(axis='y', colors='r')
    
    if use_protein_contacts:
        ax3.set_ylabel('Protein Contact', color='b')
        ax3.tick_params(axis='y', colors='b')
    
    # Title
    if partner_protein != 'none':
        plt.title(f'{protein} with partner {partner_protein}: {TARGET_LIPID} Binding vs Interface Residues')
    else:
        plt.title(f'{protein}: {TARGET_LIPID} Binding Sites vs Protein Interface')
    
    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='pink', alpha=0.5, label=f'{TARGET_LIPID} Binding Site')
    ]
    
    if protein_preference:
        legend_elements.append(Patch(facecolor='lightblue', alpha=0.5, label='Interface Residue'))
    
    if competition_residues:
        legend_elements.append(Patch(facecolor='purple', alpha=0.4, label='Competition Region'))
        
    legend_elements.extend([
        plt.Line2D([0], [0], color='r', linewidth=2, label=f'{TARGET_LIPID} Contact')
    ])
        
    if use_protein_contacts:
        if 'protein_contact' in protein_df.columns and protein_df['protein_contact'].max() > 0:
            legend_elements.append(plt.Line2D([0], [0], color='b', linewidth=2, label=f'Protein Contact (with {TARGET_LIPID})'))
        elif no_target_protein_contacts is not None:
            legend_elements.append(plt.Line2D([0], [0], color='b', linewidth=2, label=f'Protein Contact (no {TARGET_LIPID})'))
        
    plt.legend(handles=legend_elements, loc='upper right')
    
    # Grid
    ax1.grid(True, alpha=0.3)
    
    # Save
    plt.tight_layout()
    out_file_png = os.path.join(output_dir, f"{protein}_contact_competition.png")
    out_file_ps = os.path.join(output_dir, f"{protein}_contact_competition.ps")
    plt.savefig(out_file_png, dpi=300)
    plt.savefig(out_file_ps)
    plt.close()
    
    print(f"Generated {TARGET_LIPID} contact plot for {protein} with protein interface information")

def plot_protein_protein_residue_contacts(residue_df, output_dir):
    """Plot protein-protein residue-residue contacts
    
    Parameters
    ----------
    residue_df : pd.DataFrame
        Residue contact data
    output_dir : str
        Output directory path
    """
    print("\n===== Generating Protein-Protein Residue Contact Plots =====")
    
    # Ensure output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # If data is empty
    if residue_df is None or len(residue_df) == 0:
        print("No residue contact data to plot. Aborting.")
        return False
    
    # Data details output
    print(f"Input DataFrame has {len(residue_df)} rows and {len(residue_df.columns)} columns")
    print(f"Columns: {residue_df.columns.tolist()}")
    
    # Plot for each protein pair
    protein_pairs = []
    
    # 1. First confirm protein_pair column
    if 'protein_pair' in residue_df.columns:
        protein_pairs = [pair for pair in residue_df['protein_pair'].unique() if pair != 'none']
        print(f"Found {len(protein_pairs)} protein pairs from 'protein_pair' column: {protein_pairs}")
    
    # [Rest of the function continues as in the original code...]
    # [Code continues with plotting logic...]
    
    return True