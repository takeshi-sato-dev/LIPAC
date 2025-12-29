"""
Data loading and preprocessing functions for BayesianLipidAnalysis
"""

import pandas as pd
import numpy as np


def load_data(with_lipid_path, without_lipid_path):
    """Load and preprocess lipid present/absent data"""
    print("\n===== Data Loading and Preprocessing =====")
    
    # Load data
    try:
        with_lipid = pd.read_csv(with_lipid_path)
        without_lipid = pd.read_csv(without_lipid_path)
        
        print(f"With target lipid: {len(with_lipid)} rows loaded")
        print(f"Without target lipid: {len(without_lipid)} rows loaded")
        
        # Preprocessing: Add flags
        with_lipid['has_target_lipid'] = True
        without_lipid['has_target_lipid'] = False
        
        # Find target lipid contact column (could be any lipid type)
        target_lipid_col = None
        lipid_columns = [col for col in with_lipid.columns if col.endswith('_contact') and col not in ['protein_contact', 'lipid_contact']]
        
        # Find the lipid column with the highest total contact (most relevant target lipid)
        max_contact_sum = 0
        for col in lipid_columns:
            contact_sum = with_lipid[col].sum()
            print(f"  {col}: total contact = {contact_sum:.2f}")
            if contact_sum > max_contact_sum:
                max_contact_sum = contact_sum
                target_lipid_col = col
        
        # Force use DPG3 as target lipid if available and has contacts
        if 'DPG3_contact' in lipid_columns and with_lipid['DPG3_contact'].sum() > 0:
            target_lipid_col = 'DPG3_contact'
            max_contact_sum = with_lipid['DPG3_contact'].sum()
            print(f"  Forcing DPG3_contact as target lipid: total contact = {max_contact_sum:.2f}")
        
        if target_lipid_col is None or max_contact_sum == 0:
            # Fallback: create a placeholder column
            target_lipid_col = 'target_lipid_contact'
            with_lipid[target_lipid_col] = 0.0
            without_lipid[target_lipid_col] = 0.0
            print(f"No specific lipid contact column found with non-zero values, using placeholder: {target_lipid_col}")
        else:
            print(f"Using {target_lipid_col} as target lipid column (total contact: {max_contact_sum:.2f})")
        
        # Ensure target lipid column exists in both datasets
        if target_lipid_col not in with_lipid.columns:
            with_lipid[target_lipid_col] = 0.0
        if target_lipid_col not in without_lipid.columns:
            without_lipid[target_lipid_col] = 0.0
        
        # Combine data
        combined = pd.concat([with_lipid, without_lipid], ignore_index=True)
        
        # Create dataset for joining by key columns (protein, residue)
        with_lipid_key = with_lipid[['protein', 'residue', 'protein_contact', 'lipid_contact', 'ratio', target_lipid_col]]
        with_lipid_key = with_lipid_key.rename(columns={
            'protein_contact': 'protein_contact_with_lipid',
            'lipid_contact': 'lipid_contact_with_lipid',
            'ratio': 'ratio_with_lipid',
            target_lipid_col: 'target_lipid_contact'
        })
        
        without_lipid_key = without_lipid[['protein', 'residue', 'protein_contact', 'lipid_contact', 'ratio']]
        without_lipid_key = without_lipid_key.rename(columns={
            'protein_contact': 'protein_contact_without_lipid',
            'lipid_contact': 'lipid_contact_without_lipid',
            'ratio': 'ratio_without_lipid'
        })
        
        # Dataset comparable for each protein residue
        matched = pd.merge(with_lipid_key, without_lipid_key, on=['protein', 'residue'], how='inner')
        
        # Calculate differences
        matched['protein_contact_diff'] = matched['protein_contact_with_lipid'] - matched['protein_contact_without_lipid']
        matched['lipid_contact_diff'] = matched['lipid_contact_with_lipid'] - matched['lipid_contact_without_lipid']
        matched['ratio_diff'] = matched['ratio_with_lipid'] - matched['ratio_without_lipid']
        
        print(f"Matched data: {len(matched)} rows")
        print("Matched data summary:")
        print(matched.describe())
        
        return combined, matched
        
    except Exception as e:
        print(f"Data loading error: {str(e)}")
        return None, None


def prepare_residue_analysis_data(matched_df):
    """Prepare dataframe for residue-level analysis"""
    # Create copy for analysis
    analysis_df = matched_df.copy()
    
    # Required column mapping settings
    column_mappings = {
        'target_lipid_contact': ['target_lipid_contact', 'CHOL_contact', 'GM3_contact', 'DOPS_contact', 'DIPC_contact'],
        'protein_contact': ['protein_contact_with_lipid', 'protein_contact']
    }
    
    # Set target lipid contact data
    if 'target_lipid_contact' not in analysis_df.columns:
        for alt_col in column_mappings['target_lipid_contact']:
            if alt_col in analysis_df.columns:
                analysis_df['target_lipid_contact'] = analysis_df[alt_col]
                break
        else:
            analysis_df['target_lipid_contact'] = 0.0
    
    # Set protein contact data
    if 'protein_contact' not in analysis_df.columns:
        for alt_col in column_mappings['protein_contact']:
            if alt_col in analysis_df.columns:
                analysis_df['protein_contact'] = analysis_df[alt_col]
                break
    
    # Extract current protein pair relationships from data
    current_pairs = {}
    if 'protein_pair' in analysis_df.columns and 'protein' in analysis_df.columns and 'partner_protein' in analysis_df.columns:
        for _, row in analysis_df.iterrows():
            protein = row['protein']
            partner = row['partner_protein']
            if protein != 'none' and partner != 'none':
                current_pairs[protein] = partner
    
    # Set partner_protein column
    if 'partner_protein' not in analysis_df.columns or len(current_pairs) == 0:
        # Define exact pair relationships (based on data)
        protein_partners = {
            'Protein_1': 'Protein_4',
            'Protein_4': 'Protein_1',
            'Protein_2': 'Protein_3',
            'Protein_3': 'Protein_2'
        }
        analysis_df['partner_protein'] = analysis_df['protein'].map(
            lambda x: protein_partners.get(x, 'none')
        )
    
    # Set protein_pair column
    if 'protein_pair' not in analysis_df.columns:
        analysis_df['protein_pair'] = analysis_df.apply(
            lambda row: f"{min(str(row['protein']), str(row['partner_protein']))}-{max(str(row['protein']), str(row['partner_protein']))}" 
            if row['partner_protein'] != 'none' else 'none', 
            axis=1
        )
    
    # Check and set required columns
    required_columns = ['protein', 'residue', 'protein_contact', 'target_lipid_contact', 'protein_pair', 'partner_protein']
    for col in required_columns:
        if col not in analysis_df.columns:
            if col in ['protein_contact', 'target_lipid_contact']:
                analysis_df[col] = 0.0
            else:
                analysis_df[col] = 'none'
    
    return analysis_df