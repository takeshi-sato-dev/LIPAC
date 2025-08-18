"""Contact complementarity analysis functions"""

import numpy as np
import pandas as pd

try:
    from config import TARGET_LIPID
except ImportError:
    TARGET_LIPID = 'DPG3'

def analyze_contact_complementarity(protein_contact_results, lipid_contact_results):
    """Analyze relationship between protein-protein and protein-lipid contacts
    
    Parameters
    ----------
    protein_contact_results : dict
        Protein-protein contact results
    lipid_contact_results : dict
        Lipid-protein contact results
        
    Returns
    -------
    pd.DataFrame
        Contact complementarity analysis results
    """
    print("\nAnalyzing contact complementarity...")
    
    # Initialize results
    results = []
    
    # Maximum ratio value cap and minimum denominator
    MAX_RATIO = 10.0
    MIN_DENOMINATOR = 0.1
    
    # Get all protein names
    all_proteins = set()
    for lipid_type in lipid_contact_results:
        all_proteins.update(lipid_contact_results[lipid_type].keys())
    
    # Track processed proteins
    processed_proteins = set()
    
    # Process only dimerized protein pairs
    dimerized_pairs = {}
    for pair in protein_contact_results.keys():
        # Filter out invalid types
        if isinstance(pair, str) and "-" in pair:
            dimerized_pairs[pair] = True
    
    print(f"Found {len(dimerized_pairs)} dimerized protein pairs with contacts")
    
    # Process when protein-protein contacts exist
    if protein_contact_results and dimerized_pairs:
        # Process each protein pair
        for pair_name in dimerized_pairs:
            try:
                protein1_name, protein2_name = pair_name.split('-')
                processed_proteins.add(protein1_name)
                processed_proteins.add(protein2_name)
                
                # Get protein contact data for this pair
                if 'protein1' not in protein_contact_results[pair_name] or 'protein2' not in protein_contact_results[pair_name]:
                    print(f"Warning: Missing contact data for {pair_name}. Skipping.")
                    continue
                
                pp_contacts_1 = protein_contact_results[pair_name]['protein1']
                pp_contacts_2 = protein_contact_results[pair_name]['protein2']
                residue_ids1 = protein_contact_results[pair_name]['residue_ids1']
                residue_ids2 = protein_contact_results[pair_name]['residue_ids2']
                
                # Get minimum distance information if available
                min_distances1 = None
                min_distances2 = None
                if 'min_distances1' in protein_contact_results[pair_name]:
                    min_distances1 = protein_contact_results[pair_name]['min_distances1']
                if 'min_distances2' in protein_contact_results[pair_name]:
                    min_distances2 = protein_contact_results[pair_name]['min_distances2']
                
                # Process each residue in protein1
                for res_idx, pp_contact in enumerate(pp_contacts_1):
                    res_id = residue_ids1[res_idx]
                    
                    # Get minimum distance information
                    min_distance = 999.0
                    if min_distances1 is not None:
                        min_distance = min_distances1[res_idx]
                    
                    # Get lipid-protein contact data
                    lp_contacts = {}
                    for lipid_type in lipid_contact_results:
                        if protein1_name in lipid_contact_results[lipid_type]:
                            try:
                                # Search with converted residue ID
                                lp_idx = np.where(np.array(lipid_contact_results[lipid_type][protein1_name]['residue_ids']) == res_id)[0]
                                if len(lp_idx) > 0:
                                    lp_contacts[lipid_type] = lipid_contact_results[lipid_type][protein1_name]['contacts'][lp_idx[0]]
                                else:
                                    lp_contacts[lipid_type] = 0
                            except (ValueError, KeyError, IndexError) as e:
                                print(f"Warning: Error accessing lipid contact data for {protein1_name}, residue {res_id}, lipid {lipid_type}: {e}")
                                lp_contacts[lipid_type] = 0
                        else:
                            lp_contacts[lipid_type] = 0
                    
                    # Calculate total lipid contact
                    total_lp_contact = sum(lp_contacts.values())
                    
                    # Calculate ratio (with safety handling)
                    if pp_contact < 0.01:
                        if total_lp_contact < 0.01:
                            ratio = 0.0
                        else:
                            ratio = total_lp_contact / MIN_DENOMINATOR
                            ratio = min(ratio, MAX_RATIO)
                    else:
                        ratio = total_lp_contact / pp_contact
                        ratio = min(ratio, MAX_RATIO)
                    
                    # Save result with minimum distance information
                    result = {
                        'protein_pair': pair_name,
                        'protein': protein1_name,
                        'residue': res_id,
                        'protein_contact': pp_contact,
                        'min_distance': min_distance,
                        'lipid_contact': total_lp_contact,
                        'ratio': ratio,
                        'partner_protein': protein2_name
                    }
                    
                    # Add contact value for each lipid type
                    for lipid_type, contact in lp_contacts.items():
                        result[f'{lipid_type}_contact'] = contact
                    
                    # Ensure target lipid contact data is included
                    if f'{TARGET_LIPID}_contact' not in result:
                        result[f'{TARGET_LIPID}_contact'] = 0.0
                    
                    results.append(result)
                
                # Process each residue in protein2 similarly
                for res_idx, pp_contact in enumerate(pp_contacts_2):
                    res_id = residue_ids2[res_idx]
                    
                    # Get minimum distance information
                    min_distance = 999.0
                    if min_distances2 is not None:
                        min_distance = min_distances2[res_idx]
                    
                    # Get lipid-protein contact data
                    lp_contacts = {}
                    for lipid_type in lipid_contact_results:
                        if protein2_name in lipid_contact_results[lipid_type]:
                            try:
                                # Search with converted residue ID
                                lp_idx = np.where(np.array(lipid_contact_results[lipid_type][protein2_name]['residue_ids']) == res_id)[0]
                                if len(lp_idx) > 0:
                                    lp_contacts[lipid_type] = lipid_contact_results[lipid_type][protein2_name]['contacts'][lp_idx[0]]
                                else:
                                    lp_contacts[lipid_type] = 0
                            except (ValueError, KeyError, IndexError) as e:
                                print(f"Warning: Error accessing lipid contact data for {protein2_name}, residue {res_id}, lipid {lipid_type}: {e}")
                                lp_contacts[lipid_type] = 0
                        else:
                            lp_contacts[lipid_type] = 0
                    
                    # Calculate total lipid contact
                    total_lp_contact = sum(lp_contacts.values())
                    
                    # Calculate ratio (with safety handling)
                    if pp_contact < 0.01:
                        if total_lp_contact < 0.01:
                            ratio = 0.0
                        else:
                            ratio = total_lp_contact / MIN_DENOMINATOR
                            ratio = min(ratio, MAX_RATIO)
                    else:
                        ratio = total_lp_contact / pp_contact
                        ratio = min(ratio, MAX_RATIO)
                    
                    # Save result with minimum distance information
                    result = {
                        'protein_pair': pair_name,
                        'protein': protein2_name,
                        'residue': res_id,
                        'protein_contact': pp_contact,
                        'min_distance': min_distance,
                        'lipid_contact': total_lp_contact,
                        'ratio': ratio,
                        'partner_protein': protein1_name
                    }
                    
                    # Add contact value for each lipid type
                    for lipid_type, contact in lp_contacts.items():
                        result[f'{lipid_type}_contact'] = contact
                    
                    # Ensure target lipid contact data is included
                    if f'{TARGET_LIPID}_contact' not in result:
                        result[f'{TARGET_LIPID}_contact'] = 0.0
                    
                    results.append(result)
            except Exception as e:
                print(f"Error processing pair {pair_name}: {str(e)}")
    
    # Process all proteins including those without protein-protein contacts
    for protein_name in all_proteins:
        processed_lipid_residues = set()  # Track processed residue IDs
        
        for lipid_type in lipid_contact_results:
            if protein_name in lipid_contact_results[lipid_type]:
                lipid_data = lipid_contact_results[lipid_type][protein_name]
                
                if 'contacts' not in lipid_data or 'residue_ids' not in lipid_data:
                    continue
                
                contacts = lipid_data['contacts']
                residue_ids = lipid_data['residue_ids']
                
                for res_idx, lipid_contact in enumerate(contacts):
                    res_id = residue_ids[res_idx]
                    
                    # Skip already processed residues
                    res_key = f"{protein_name}_{res_id}"
                    if res_key in processed_lipid_residues:
                        continue
                    processed_lipid_residues.add(res_key)
                    
                    # Get all lipid type contacts for this residue
                    all_lipid_contacts = {}
                    total_contact = 0
                    
                    for inner_lipid_type in lipid_contact_results:
                        if protein_name in lipid_contact_results[inner_lipid_type]:
                            inner_lipid_data = lipid_contact_results[inner_lipid_type][protein_name]
                            
                            if 'contacts' not in inner_lipid_data or 'residue_ids' not in inner_lipid_data:
                                all_lipid_contacts[inner_lipid_type] = 0
                                continue
                            
                            # Search for this residue ID
                            try:
                                inner_res_idx = np.where(np.array(inner_lipid_data['residue_ids']) == res_id)[0]
                                if len(inner_res_idx) > 0:
                                    contact_val = inner_lipid_data['contacts'][inner_res_idx[0]]
                                    all_lipid_contacts[inner_lipid_type] = contact_val
                                    total_contact += contact_val
                                else:
                                    all_lipid_contacts[inner_lipid_type] = 0
                            except Exception as e:
                                print(f"Error finding residue {res_id} in {inner_lipid_type} data: {e}")
                                all_lipid_contacts[inner_lipid_type] = 0
                        else:
                            all_lipid_contacts[inner_lipid_type] = 0
                    
                    # Get protein contact data (0 if not available)
                    protein_contact = 0.0
                    min_distance = 999.0
                    partner_protein = 'none'
                    protein_pair = 'none'
                    
                    # Check if this protein is in any dimer pair
                    for pair_name in dimerized_pairs:
                        # Validate pair name
                        if not isinstance(pair_name, str) or '-' not in pair_name:
                            continue
                            
                        try:
                            pair_parts = pair_name.split('-')
                            if len(pair_parts) != 2:
                                continue
                                
                            protein1_name, protein2_name = pair_parts
                            
                            if protein_name == protein1_name:
                                if 'protein1' in protein_contact_results[pair_name] and 'residue_ids1' in protein_contact_results[pair_name]:
                                    residue_indices = np.where(np.array(protein_contact_results[pair_name]['residue_ids1']) == res_id)[0]
                                    if len(residue_indices) > 0:
                                        idx = residue_indices[0]
                                        protein_contact = protein_contact_results[pair_name]['protein1'][idx]
                                        if 'min_distances1' in protein_contact_results[pair_name]:
                                            min_distance = protein_contact_results[pair_name]['min_distances1'][idx]
                                        partner_protein = protein2_name
                                        protein_pair = pair_name
                            elif protein_name == protein2_name:
                                if 'protein2' in protein_contact_results[pair_name] and 'residue_ids2' in protein_contact_results[pair_name]:
                                    residue_indices = np.where(np.array(protein_contact_results[pair_name]['residue_ids2']) == res_id)[0]
                                    if len(residue_indices) > 0:
                                        idx = residue_indices[0]
                                        protein_contact = protein_contact_results[pair_name]['protein2'][idx]
                                        if 'min_distances2' in protein_contact_results[pair_name]:
                                            min_distance = protein_contact_results[pair_name]['min_distances2'][idx]
                                        partner_protein = protein1_name
                                        protein_pair = pair_name
                        except Exception as e:
                            print(f"Error checking protein contact in pair {pair_name}: {str(e)}")
                    
                    # Calculate ratio (with safety handling)
                    if protein_contact < 0.01:
                        if total_contact < 0.01:
                            ratio = 0.0
                        else:
                            ratio = total_contact / MIN_DENOMINATOR
                            ratio = min(ratio, MAX_RATIO)
                    else:
                        ratio = total_contact / protein_contact
                        ratio = min(ratio, MAX_RATIO)
                    
                    # Create lipid contact data (including cases with protein contacts)
                    result = {
                        'protein_pair': protein_pair,
                        'protein': protein_name,
                        'residue': res_id,
                        'protein_contact': protein_contact,
                        'min_distance': min_distance,
                        'lipid_contact': total_contact,
                        'ratio': ratio,
                        'partner_protein': partner_protein
                    }
                    
                    # Add contact value for each lipid type
                    for l_type, contact_value in all_lipid_contacts.items():
                        result[f'{l_type}_contact'] = contact_value
                    
                    # Ensure target lipid contact data is included
                    if f'{TARGET_LIPID}_contact' not in result:
                        result[f'{TARGET_LIPID}_contact'] = 0.0
                    
                    results.append(result)
    
    # Convert to DataFrame
    if results:
        df = pd.DataFrame(results)
        
        # Debug output: Target lipid contact data confirmation
        if f'{TARGET_LIPID}_contact' in df.columns:
            target_max = df[f'{TARGET_LIPID}_contact'].max()
            target_nonzero = (df[f'{TARGET_LIPID}_contact'] > 0).sum()
            print(f"{TARGET_LIPID}_contact max value: {target_max}, non-zero values: {target_nonzero}")
        else:
            print(f"Warning: {TARGET_LIPID}_contact column not found")
            
        # Check and add required columns
        required_columns = ['protein_pair', 'protein', 'residue', 'protein_contact', 
                           'min_distance', 'lipid_contact', 'ratio', 'partner_protein']
        
        for col in required_columns:
            if col not in df.columns:
                print(f"Warning: {col} column not found. Setting default value.")
                if col == 'protein_pair':
                    # Generate protein_pair if partner_protein exists
                    if 'partner_protein' in df.columns:
                        df[col] = df.apply(
                            lambda row: f"{min(str(row['protein']), str(row['partner_protein']))}-{max(str(row['protein']), str(row['partner_protein']))}" 
                            if row['partner_protein'] != 'none' else 'none', 
                            axis=1
                        )
                    else:
                        df[col] = 'none'
                elif col == 'partner_protein':
                    df[col] = 'none'
                elif col in ['protein_contact', 'min_distance', 'lipid_contact', 'ratio']:
                    df[col] = 0.0
                else:
                    df[col] = 'unknown'
        
        # Residue range confirmation
        print(f"Residue range: {df['residue'].min()}-{df['residue'].max()}")
        
        return df
    else:
        print("No complementarity data collected")
        return None