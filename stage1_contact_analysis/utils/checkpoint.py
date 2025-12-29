"""Checkpoint management utilities for resuming long calculations"""

import os
import pickle
import traceback

def save_checkpoint(checkpoint_file, processed_frames, all_protein_contacts, 
                   all_lipid_contacts, all_residue_contacts, selected_leaflet0=None,
                   all_target_lipid_binding_states=None, all_lipid_compositions=None,
                   all_unique_lipid_compositions=None):
    """Save checkpoint data to file
    
    Parameters
    ----------
    checkpoint_file : str
        Path to checkpoint file
    processed_frames : list
        List of processed frame indices
    all_protein_contacts : dict
        All protein contact data
    all_lipid_contacts : dict
        All lipid contact data
    all_residue_contacts : list
        All residue contact data
    selected_leaflet0 : MDAnalysis.AtomGroup, optional
        Selected leaflet information
    all_target_lipid_binding_states : dict, optional
        Target lipid binding states for causal analysis
    all_lipid_compositions : dict, optional
        Lipid compositions for causal analysis
    all_unique_lipid_compositions : dict, optional
        Unique molecule compositions for causal analysis
        
    Returns
    -------
    bool
        True if successful
    """
    try:
        checkpoint_data = {
            'processed_frames': processed_frames,
            'all_protein_contacts': all_protein_contacts,
            'all_lipid_contacts': all_lipid_contacts,
            'all_residue_contacts': all_residue_contacts,
            'selected_leaflet0': selected_leaflet0,  # Also save leaflet information
            'all_target_lipid_binding_states': all_target_lipid_binding_states or {},
            'all_lipid_compositions': all_lipid_compositions or {},
            'all_unique_lipid_compositions': all_unique_lipid_compositions or {}
        }
        
        with open(checkpoint_file, 'wb') as f:
            pickle.dump(checkpoint_data, f)
        print(f"Checkpoint saved successfully to {checkpoint_file}")
        return True
    except Exception as e:
        print(f"Warning: Failed to save checkpoint: {str(e)}")
        traceback.print_exc()
        return False

def load_checkpoint(checkpoint_file):
    """Load checkpoint data from file
    
    Parameters
    ----------
    checkpoint_file : str
        Path to checkpoint file
        
    Returns
    -------
    dict or None
        Checkpoint data if successful, None otherwise
    """
    if not os.path.exists(checkpoint_file):
        return None
        
    print(f"Found checkpoint file: {checkpoint_file}")
    print("Attempting to resume from checkpoint...")
    
    try:
        with open(checkpoint_file, 'rb') as f:
            checkpoint_data = pickle.load(f)
            
        # Get information from checkpoint
        processed_frames = checkpoint_data.get('processed_frames', [])
        all_protein_contacts = checkpoint_data.get('all_protein_contacts', {})
        all_lipid_contacts = checkpoint_data.get('all_lipid_contacts', {})
        all_residue_contacts = checkpoint_data.get('all_residue_contacts', [])
        
        # Also restore leaflet information (if exists)
        selected_leaflet0 = None
        if 'selected_leaflet0' in checkpoint_data:
            selected_leaflet0 = checkpoint_data['selected_leaflet0']
            if selected_leaflet0 is not None:
                print("Successfully restored leaflet information from checkpoint")
                print(f"Restored leaflet contains {len(selected_leaflet0.residues) if selected_leaflet0 else 0} residues")
                
                # Verify leaflet information - check target lipid molecules
                from ..config import TARGET_LIPID
                try:
                    target_count = len(selected_leaflet0.select_atoms(f"resname {TARGET_LIPID}").residues)
                    print(f"Checkpoint leaflet contains {target_count} {TARGET_LIPID} molecules")
                except Exception as e:
                    print(f"Error checking {TARGET_LIPID} in checkpoint leaflet: {str(e)}")
            else:
                print("Leaflet information in checkpoint is None")
        else:
            print("No leaflet information key in checkpoint")
        
        print(f"Successfully loaded checkpoint with {len(processed_frames)} processed frames")
        
        # Also restore causal analysis data
        all_target_lipid_binding_states = checkpoint_data.get('all_target_lipid_binding_states', {})
        all_lipid_compositions = checkpoint_data.get('all_lipid_compositions', {})
        
        print(f"Checkpoint causal data: {len(all_target_lipid_binding_states)} binding states, {len(all_lipid_compositions)} compositions")
        
        return {
            'processed_frames': processed_frames,
            'all_protein_contacts': all_protein_contacts,
            'all_lipid_contacts': all_lipid_contacts,
            'all_residue_contacts': all_residue_contacts,
            'selected_leaflet0': selected_leaflet0,
            'all_target_lipid_binding_states': all_target_lipid_binding_states,
            'all_lipid_compositions': all_lipid_compositions,
            'resume_from_frame': max(processed_frames) + 1 if processed_frames else 0
        }
        
    except Exception as e:
        print(f"Error loading checkpoint: {str(e)}")
        traceback.print_exc()
        return None

def create_intermediate_csv(all_protein_contacts, all_lipid_contacts, output_dir):
    """Create intermediate CSV file from current data
    
    Parameters
    ----------
    all_protein_contacts : dict
        All protein contact data
    all_lipid_contacts : dict
        All lipid contact data
    output_dir : str
        Output directory path
        
    Returns
    -------
    bool
        True if successful
    """
    from ..analysis.complementarity import analyze_contact_complementarity
    from ..config import TARGET_LIPID
    
    print("Generating intermediate complementarity CSV file...")
    try:
        # Perform complementarity analysis from currently collected data
        temp_complementarity_df = analyze_contact_complementarity(all_protein_contacts, all_lipid_contacts)
        
        # Save results as CSV
        if temp_complementarity_df is not None and len(temp_complementarity_df) > 0:
            temp_csv_path = os.path.join(output_dir, "contact_complementarity_temp.csv")
            temp_complementarity_df.to_csv(temp_csv_path, index=False)
            print(f"SUCCESS: Saved intermediate complementarity data to {temp_csv_path} ({len(temp_complementarity_df)} rows)")
            
            # Validate target lipid data
            target_col = f'{TARGET_LIPID}_contact'
            if target_col in temp_complementarity_df.columns:
                max_target = temp_complementarity_df[target_col].max()
                count_target = (temp_complementarity_df[target_col] > 0).sum()
                print(f"{TARGET_LIPID} contact data: max value = {max_target}, positive values = {count_target}")
            else:
                print(f"WARNING: No {TARGET_LIPID}_contact column in complementarity data!")
            
            return True
        else:
            print("WARNING: No complementarity data could be generated yet")
            return False
    except Exception as e:
        print(f"ERROR: Failed to generate intermediate complementarity file: {str(e)}")
        traceback.print_exc()
        return False