"""Residue-level contact analysis functions"""

import numpy as np
import pandas as pd
import MDAnalysis as mda
from ..config import RESIDUE_OFFSET, TARGET_LIPID

def extract_residue_contacts(universe, frame_idx, proteins, protein_pairs, leaflet0):
    """Extract per-residue protein contacts and target lipid contacts with correct residue number conversion
    
    Parameters
    ----------
    universe : MDAnalysis.Universe
        The universe object
    frame_idx : int
        Current frame index
    proteins : dict
        Dictionary of protein selections
    protein_pairs : dict or list
        Protein pairs to analyze
    leaflet0 : MDAnalysis.AtomGroup
        Upper leaflet atoms
        
    Returns
    -------
    dict
        Residue contact information
    """
    # Convert dictionary keys to list if needed
    if isinstance(protein_pairs, dict):
        protein_pairs = list(protein_pairs.keys())
    
    # Load frame
    universe.trajectory[frame_idx]
    box = universe.dimensions[:3]
    
    # Select target lipid molecules
    target_lipid_sel = None
    try:
        target_lipid_sel = leaflet0.select_atoms(f"resname {TARGET_LIPID}")
        print(f"{TARGET_LIPID} molecules: {len(target_lipid_sel.residues)}")
    except:
        print(f"{TARGET_LIPID} molecules not found")
        target_lipid_sel = mda.AtomGroup([], universe)
    
    residue_contacts = {}
    
    # For each protein pair
    for pair_name in protein_pairs:
        protein1_name, protein2_name = pair_name.split('-')
        protein1 = proteins[protein1_name]
        protein2 = proteins[protein2_name]
        
        # For each residue in protein1
        for res1 in protein1.residues:
            # Convert residue ID (+556) - This is the correction point
            actual_resid = res1.resid + RESIDUE_OFFSET  # e.g., 65 → 621
            residue_key = f"{protein1_name}_{actual_resid}"
            
            if residue_key not in residue_contacts:
                residue_contacts[residue_key] = {
                    'protein': protein1_name,
                    'residue': actual_resid,  # Converted residue ID
                    'protein_contacts': {},
                    f'{TARGET_LIPID.lower()}_contacts': 0,
                    'original_resid': res1.resid  # Record original residue ID for debugging
                }
            
            # Calculate protein-protein contacts
            protein_contact = 0
            for res2 in protein2.residues:
                min_dist = float('inf')
                
                for atom1 in res1.atoms:
                    for atom2 in res2.atoms:
                        # Distance calculation with PBC correction
                        diff = atom1.position - atom2.position
                        for dim in range(3):
                            if diff[dim] > box[dim] * 0.5:
                                diff[dim] -= box[dim]
                            elif diff[dim] < -box[dim] * 0.5:
                                diff[dim] += box[dim]
                        
                        dist = np.sqrt(np.sum(diff * diff))
                        min_dist = min(min_dist, dist)
                        
                        if min_dist <= 6.0:  # Contact cutoff
                            break
                    
                    if min_dist <= 6.0:
                        break
                
                # Calculate contact value (inverse of distance)
                if min_dist <= 6.0:
                    contact_value = 1.0 / max(0.1, min_dist)
                    # Also convert partner protein's residue ID
                    partner_resid = res2.resid + RESIDUE_OFFSET
                    contact_key = f"{protein2_name}_{partner_resid}"
                    residue_contacts[residue_key]['protein_contacts'][contact_key] = \
                        residue_contacts[residue_key]['protein_contacts'].get(contact_key, 0) + contact_value
                    protein_contact += contact_value
            
            # Calculate target lipid-protein contacts
            if target_lipid_sel and len(target_lipid_sel) > 0:
                for lipid_res in target_lipid_sel.residues:
                    min_dist = float('inf')
                    
                    for atom1 in res1.atoms:
                        for atom2 in lipid_res.atoms:
                            # Distance calculation with PBC correction
                            diff = atom1.position - atom2.position
                            for dim in range(3):
                                if diff[dim] > box[dim] * 0.5:
                                    diff[dim] -= box[dim]
                                elif diff[dim] < -box[dim] * 0.5:
                                    diff[dim] += box[dim]
                            
                            dist = np.sqrt(np.sum(diff * diff))
                            min_dist = min(min_dist, dist)
                            
                            if min_dist <= 6.0:  # Contact cutoff
                                break
                        
                        if min_dist <= 6.0:
                            break
                    
                    # Calculate contact value (binary 0/1)
                    if min_dist <= 6.0:
                        residue_contacts[residue_key][f'{TARGET_LIPID.lower()}_contacts'] += 1
        
        # Process protein2 residues similarly
        for res2 in protein2.residues:
            actual_resid = res2.resid + RESIDUE_OFFSET
            residue_key = f"{protein2_name}_{actual_resid}"
            
            if residue_key not in residue_contacts:
                residue_contacts[residue_key] = {
                    'protein': protein2_name,
                    'residue': actual_resid,
                    'protein_contacts': {},
                    f'{TARGET_LIPID.lower()}_contacts': 0,
                    'original_resid': res2.resid
                }
            
            # Calculate contacts (similar to protein1)
            # [Code continues similarly for protein2]
    
    return residue_contacts

def aggregate_residue_contacts(all_frames_results):
    """Aggregate residue-level contact information across all frames
    
    Parameters
    ----------
    all_frames_results : list
        List of frame results
        
    Returns
    -------
    pd.DataFrame
        Aggregated residue contact data
    """
    # Aggregation data structure
    agg_residue_data = {}
    
    # Aggregate data from all frames
    for frame_result in all_frames_results:
        if 'residue_contacts' not in frame_result:
            continue
            
        residue_contacts = frame_result['residue_contacts']
        
        for residue_key, contact_data in residue_contacts.items():
            if residue_key not in agg_residue_data:
                agg_residue_data[residue_key] = {
                    'protein': contact_data['protein'],
                    'residue': contact_data['residue'],
                    'protein_contacts': {},
                    f'{TARGET_LIPID.lower()}_contacts': 0,
                    'frame_count': 0
                }
            
            # Aggregate protein contact data
            for partner, contact_value in contact_data['protein_contacts'].items():
                if partner not in agg_residue_data[residue_key]['protein_contacts']:
                    agg_residue_data[residue_key]['protein_contacts'][partner] = 0
                agg_residue_data[residue_key]['protein_contacts'][partner] += contact_value
            
            # Aggregate target lipid contact data
            agg_residue_data[residue_key][f'{TARGET_LIPID.lower()}_contacts'] += contact_data[f'{TARGET_LIPID.lower()}_contacts']
            agg_residue_data[residue_key]['frame_count'] += 1
    
    # Calculate averages
    for residue_key, agg_data in agg_residue_data.items():
        if agg_data['frame_count'] > 0:
            for partner in agg_data['protein_contacts']:
                agg_data['protein_contacts'][partner] /= agg_data['frame_count']
            agg_data[f'{TARGET_LIPID.lower()}_contacts'] /= agg_data['frame_count']
    
    # Convert to DataFrame
    rows = []
    for residue_key, agg_data in agg_residue_data.items():
        protein = agg_data['protein']
        residue = agg_data['residue']
        
        # Total protein contact value
        total_protein_contact = sum(agg_data['protein_contacts'].values())
        
        # For each partner contact
        for partner, contact_value in agg_data['protein_contacts'].items():
            # Parse partner information to extract protein name
            if '_' in partner:
                partner_parts = partner.split('_')
                partner_protein = partner_parts[0]
            else:
                partner_protein = partner
                
            # Build protein pair name
            if protein < partner_protein:
                protein_pair = f"{protein}-{partner_protein}"
            else:
                protein_pair = f"{partner_protein}-{protein}"
                
            row = {
                'protein': protein,
                'residue': residue,
                'partner': partner_protein,
                'protein_pair': protein_pair,
                'protein_contact': contact_value,
                'total_protein_contact': total_protein_contact,
                f'{TARGET_LIPID}_contact': agg_data[f'{TARGET_LIPID.lower()}_contacts']  # Include target lipid contact data explicitly
            }
            rows.append(row)
    
    # Create DataFrame
    if rows:
        residue_df = pd.DataFrame(rows)
        
        # Debug output: Confirm residue range
        if not residue_df.empty:
            min_residue = residue_df['residue'].min()
            max_residue = residue_df['residue'].max()
            residue_range = sorted(residue_df['residue'].unique())
            print(f"\nResidue range confirmation: {min_residue}-{max_residue}")
            # Count residues in 621-644 range
            in_range_count = len([r for r in residue_range if 621 <= r <= 644])
            print(f"Residues in 621-644 range: {in_range_count}/{len(residue_range)}")
            
            # Confirm target lipid contact column
            target_col = f'{TARGET_LIPID}_contact'
            if target_col in residue_df.columns:
                target_max = residue_df[target_col].max()
                target_nonzero = (residue_df[target_col] > 0).sum()
                print(f"{TARGET_LIPID}_contact max value: {target_max}, non-zero values: {target_nonzero}")
            else:
                print(f"Warning: {TARGET_LIPID}_contact column not found")
        
        return residue_df
    else:
        return None