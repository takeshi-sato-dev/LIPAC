#!/usr/bin/env python3
"""Main entry point for contact analysis."""

import sys
import os
import time
import traceback
import argparse
from datetime import datetime

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import with correct function names
from stage1_contact_analysis.core.trajectory_loader import (
    load_universe, 
    identify_lipid_leaflets, 
    select_proteins,
    select_lipids
)
from stage1_contact_analysis.core.frame_processor import (
    process_frame,
    process_frames_serially
)
from stage1_contact_analysis.core.contact_calculator import (
    calculate_protein_com_distances,
    calculate_protein_protein_contacts,
    calculate_lipid_protein_contacts,
    check_tm_helix_interactions
)
from stage1_contact_analysis.analysis.complementarity import analyze_contact_complementarity
from stage1_contact_analysis.analysis.residue_contacts import (
    aggregate_residue_contacts,
    extract_residue_contacts
)
from stage1_contact_analysis.utils.checkpoint import (
    save_checkpoint, 
    load_checkpoint,
    create_intermediate_csv
)
from stage1_contact_analysis.utils.parallel import (
    process_batch,
    test_multiprocessing
)
from stage1_contact_analysis.visualization.plotting import (
    plot_contact_complementarity,
    plot_target_lipid_protein_competition,  
    plot_protein_protein_residue_contacts
)
from stage1_contact_analysis.analysis.comparison import compare_with_without_target_lipid
from stage1_contact_analysis.config import *


def save_target_lipid_causal_data(target_lipid_binding_states, lipid_compositions, output_dir, target_lipid_name='TARGET_LIPID', unique_lipid_compositions=None):
    """Save target lipid binding states and lipid composition data for causal analysis
    
    Parameters
    ----------
    target_lipid_binding_states : dict
        Frame-wise target lipid binding states for each protein
    lipid_compositions : dict
        Lipid compositions with and without target lipid
    output_dir : str
        Output directory path
    target_lipid_name : str
        Name of the target lipid for output files
    """
    print(f"\n----- Debug: save_target_lipid_causal_data called -----")
    print(f"Output directory: {output_dir}")
    print(f"Target lipid name: {target_lipid_name}")
    print(f"Binding states: {len(target_lipid_binding_states)} frames")
    print(f"Compositions: {lipid_compositions.keys() if lipid_compositions else 'None'}")
    import pandas as pd
    
    print(f"\n----- Saving {target_lipid_name} Causal Analysis Data -----")
    
    # Handle empty data gracefully
    if not target_lipid_binding_states:
        print(f"WARNING: No target lipid binding states provided - creating empty causal data files")
        target_lipid_binding_states = {}
    
    if not lipid_compositions:
        print(f"WARNING: No lipid compositions provided - creating empty causal data files")
        lipid_compositions = {'with_target_lipid': {}, 'without_target_lipid': {}}
    
    # Create DataFrames for each protein
    # Get protein names from binding states or lipid compositions
    protein_names = set()
    if target_lipid_binding_states:
        protein_names.update(set().union(*[set(frame_data.keys()) for frame_data in target_lipid_binding_states.values()]))
    
    # Also get protein names from lipid compositions if available
    for comp_type in ['with_target_lipid', 'without_target_lipid']:
        if comp_type in lipid_compositions:
            protein_names.update(lipid_compositions[comp_type].keys())
    
    # If still no protein names, create default protein names
    if not protein_names:
        print(f"WARNING: No protein names found - creating placeholder files")
        protein_names = {'Protein_1', 'Protein_2', 'Protein_3', 'Protein_4'}  # Default protein names
    
    for protein_name in protein_names:
        protein_data = []
        
        # Collect all frames for this protein
        if target_lipid_binding_states:
            for frame_idx, frame_states in target_lipid_binding_states.items():
                if protein_name in frame_states:
                    row = {
                        'frame': frame_idx,
                        'protein': protein_name,
                        'target_lipid_bound': frame_states[protein_name]
                    }
                    
                    # Add lipid composition data (residue contacts)
                    if frame_states[protein_name]:  # Target lipid bound
                        for comp_data in lipid_compositions.get('with_target_lipid', {}).get(protein_name, []):
                            if comp_data['frame'] == frame_idx:
                                for lipid_type, count in comp_data['composition'].items():
                                    row[f'{lipid_type}_contacts'] = count
                                break
                        # Also add unique molecule data if available
                        if unique_lipid_compositions:
                            for comp_data in unique_lipid_compositions.get('with_target_lipid', {}).get(protein_name, []):
                                if comp_data['frame'] == frame_idx:
                                    for lipid_type, count in comp_data['composition'].items():
                                        row[f'{lipid_type}_unique_molecules'] = count
                                    break
                    else:  # Target lipid not bound
                        for comp_data in lipid_compositions.get('without_target_lipid', {}).get(protein_name, []):
                            if comp_data['frame'] == frame_idx:
                                for lipid_type, count in comp_data['composition'].items():
                                    row[f'{lipid_type}_contacts'] = count
                                break
                        # Also add unique molecule data if available
                        if unique_lipid_compositions:
                            for comp_data in unique_lipid_compositions.get('without_target_lipid', {}).get(protein_name, []):
                                if comp_data['frame'] == frame_idx:
                                    for lipid_type, count in comp_data['composition'].items():
                                        row[f'{lipid_type}_unique_molecules'] = count
                                    break
                    
                    protein_data.append(row)
        
        # Always save CSV file, even if empty
        if protein_data:
            df = pd.DataFrame(protein_data)
            df = df.sort_values('frame')
            bound_count = df['target_lipid_bound'].sum()
        else:
            # Create empty dataframe with required columns
            df = pd.DataFrame(columns=['frame', 'protein', 'target_lipid_bound'])
            bound_count = 0
        
        # Save to CSV
        csv_path = os.path.join(output_dir, f'{target_lipid_name.lower()}_causal_data_{protein_name}.csv')
        df.to_csv(csv_path, index=False)
        print(f"  Saved {protein_name}: {len(df)} frames, {target_lipid_name} bound in {bound_count} frames")
    
    # Also save combined data for all proteins
    all_data = []
    if target_lipid_binding_states:
        for frame_idx, frame_states in target_lipid_binding_states.items():
            for protein_name, target_lipid_bound in frame_states.items():
                row = {
                    'frame': frame_idx,
                    'protein': protein_name,
                    'target_lipid_bound': target_lipid_bound
                }
                
                if target_lipid_bound:
                    for comp_data in lipid_compositions.get('with_target_lipid', {}).get(protein_name, []):
                        if comp_data['frame'] == frame_idx:
                            for lipid_type, count in comp_data['composition'].items():
                                row[f'{lipid_type}_contacts'] = count
                            break
                    # Also add unique molecule data if available
                    if unique_lipid_compositions:
                        for comp_data in unique_lipid_compositions.get('with_target_lipid', {}).get(protein_name, []):
                            if comp_data['frame'] == frame_idx:
                                for lipid_type, count in comp_data['composition'].items():
                                    row[f'{lipid_type}_unique_molecules'] = count
                                break
                else:
                    for comp_data in lipid_compositions.get('without_target_lipid', {}).get(protein_name, []):
                        if comp_data['frame'] == frame_idx:
                            for lipid_type, count in comp_data['composition'].items():
                                row[f'{lipid_type}_contacts'] = count
                            break
                    # Also add unique molecule data if available
                    if unique_lipid_compositions:
                        for comp_data in unique_lipid_compositions.get('without_target_lipid', {}).get(protein_name, []):
                            if comp_data['frame'] == frame_idx:
                                for lipid_type, count in comp_data['composition'].items():
                                    row[f'{lipid_type}_unique_molecules'] = count
                                break
                
                all_data.append(row)
    
    # Always save combined data file, even if empty
    if all_data:
        all_df = pd.DataFrame(all_data)
        all_df = all_df.sort_values(['protein', 'frame'])
    else:
        # Create empty dataframe with required columns
        all_df = pd.DataFrame(columns=['frame', 'protein', 'target_lipid_bound'])
    
    all_csv_path = os.path.join(output_dir, f'{target_lipid_name.lower()}_causal_data_all_proteins.csv')
    all_df.to_csv(all_csv_path, index=False)
    print(f"  Saved combined data: {len(all_df)} total entries")
    
    print(f"✓ {target_lipid_name} causal analysis data saved successfully")

# Variable name mapping for compatibility
WITH_TARGET_LIPID_PSF = DEFAULT_WITH_LIPID_PSF
WITH_TARGET_LIPID_XTC = DEFAULT_WITH_LIPID_XTC
WITHOUT_TARGET_LIPID_PSF = DEFAULT_WITHOUT_LIPID_PSF
WITHOUT_TARGET_LIPID_XTC = DEFAULT_WITHOUT_LIPID_XTC
WITH_TARGET_LIPID_OUTPUT = DEFAULT_WITH_LIPID_OUTPUT
WITHOUT_TARGET_LIPID_OUTPUT = DEFAULT_WITHOUT_LIPID_OUTPUT

def run_analysis(top_file, traj_file, output_dir, start_frame=START_FRAME, stop_frame=STOP_FRAME, 
               step_frame=STEP_FRAME, batch_size=BATCH_SIZE, debug_mode=False, force_parallel=True,
               mp_context=None, contact_cutoff=CONTACT_CUTOFF, min_cores=MIN_CORES):
    """Analyze system and generate contact_complementarity.csv and residue_contacts.csv
    
    Parameters
    ----------
    top_file : str
        Path to topology file
    traj_file : str
        Path to trajectory file
    output_dir : str
        Output directory path
    start_frame : int
        Starting frame
    stop_frame : int
        Ending frame
    step_frame : int
        Frame interval
    batch_size : int
        Batch size for processing
    debug_mode : bool
        Debug mode flag
    force_parallel : bool
        Force parallel processing
    mp_context : str
        Multiprocessing context
    contact_cutoff : float
        Contact cutoff distance
    min_cores : int
        Minimum number of cores
        
    Returns
    -------
    tuple
        (complementarity_df, residue_df)
    """
    print("\n===== Analyzing System =====")
    
    # Record start time
    start_time = time.time()
    
    # Ensure and create output directory - use absolute path
    abs_output_dir = os.path.abspath(output_dir)
    os.makedirs(abs_output_dir, exist_ok=True)
    print(f"Topology: {top_file}")
    print(f"Trajectory: {traj_file}")
    print(f"Output Directory: {abs_output_dir} (Absolute Path)")
    
    # Checkpoint file path
    checkpoint_file = os.path.join(abs_output_dir, "analysis_checkpoint.pkl")
    
    # Leaflet information file paths - use config setting
    leaflet_info_file = os.path.join(TEMP_FILES_DIR, "leaflet_info.pickle")
    leaflet_info_txt = os.path.join(TEMP_FILES_DIR, "leaflet_info.txt")
    
    # Selected leaflet
    selected_leaflet0 = None
    
    # Check for existing checkpoint
    checkpoint_data = load_checkpoint(checkpoint_file)
    if checkpoint_data:
        processed_frames = checkpoint_data.get('processed_frames', [])
        all_protein_contacts = checkpoint_data.get('all_protein_contacts', {})
        all_lipid_contacts = checkpoint_data.get('all_lipid_contacts', {})
        all_residue_contacts = checkpoint_data.get('all_residue_contacts', [])
        selected_leaflet0 = checkpoint_data.get('selected_leaflet0')
        # NEW: Initialize target lipid tracking variables from checkpoint
        all_target_lipid_binding_states = checkpoint_data.get('all_target_lipid_binding_states', {})
        all_lipid_compositions = checkpoint_data.get('all_lipid_compositions', {
            'with_target_lipid': {},
            'without_target_lipid': {}
        })
        all_unique_lipid_compositions = checkpoint_data.get('all_unique_lipid_compositions', {
            'with_target_lipid': {},
            'without_target_lipid': {}
        })
        print(f"Restored from checkpoint: {len(all_target_lipid_binding_states)} binding states, {len(all_lipid_compositions.get('with_target_lipid', {}))} with_target compositions")
        
        if processed_frames:
            resume_from_frame = max(processed_frames) + step_frame
            print(f"Resuming from frame {resume_from_frame}")
        else:
            resume_from_frame = start_frame
    else:
        print("No checkpoint found. Starting from scratch.")
        processed_frames = []
        all_protein_contacts = {}
        all_lipid_contacts = {}
        all_residue_contacts = []
        all_target_lipid_binding_states = {}  # NEW: Track target lipid binding states
        all_lipid_compositions = {    # NEW: Track lipid compositions
            'with_target_lipid': {},
            'without_target_lipid': {}
        }
        all_unique_lipid_compositions = {    # NEW: Track unique molecule compositions
            'with_target_lipid': {},
            'without_target_lipid': {}
        }
        resume_from_frame = start_frame
    
    # Create temporary files directory (for leaflet info, etc.)
    temp_files_dir = os.path.abspath(TEMP_FILES_DIR)
    os.makedirs(temp_files_dir, exist_ok=True)
    
    # Process leaflet information files
    if selected_leaflet0 is None:
        if os.path.exists(leaflet_info_file):
            print(f"No leaflet information in checkpoint - removing existing leaflet_info.pickle file")
            try:
                os.remove(leaflet_info_file)
            except Exception as e:
                print(f"Warning: Could not remove leaflet file: {str(e)}")
        
        if os.path.exists(leaflet_info_txt):
            print(f"No leaflet information in checkpoint - removing existing leaflet_info.txt file")
            try:
                os.remove(leaflet_info_txt)
            except Exception as e:
                print(f"Warning: Could not remove leaflet text file: {str(e)}")
    else:
        print("Using leaflet information from checkpoint - not regenerating")
    
    # Load trajectory and check
    try:
        universe = load_universe(top_file, traj_file)
        if universe is None:
            return None, None
            
        print(f"Successfully loaded trajectory with {len(universe.trajectory)} frames.")
        print(f"System contains {len(universe.atoms)} atoms and {len(universe.residues)} residues.")
        
        # Check lipid composition
        try:
            if hasattr(config, 'LIPID_TYPES'):
                lipid_resnames = " ".join(LIPID_TYPES)
            else:
                lipid_resnames = "CHOL DPSM DIPC DPG3 DOPS"
            
            lipid_atoms = universe.select_atoms(f"resname {lipid_resnames}")
            
            print("\nSystem lipid composition:")
            for lipid_type in lipid_resnames.split():
                try:
                    count = len(universe.select_atoms(f"resname {lipid_type}").residues)
                    print(f"  {lipid_type}: {count} molecules")
                except Exception as e:
                    print(f"  Error counting {lipid_type}: {str(e)}")
        except Exception as e:
            print(f"Warning: Could not check lipid composition: {str(e)}")
    except Exception as e:
        print(f"Failed to load universe: {str(e)}")
        traceback.print_exc()
        return None, None

    # Set frame range
    n_frames = len(universe.trajectory)
    stop_frame_val = min(stop_frame, n_frames) if stop_frame > 0 else n_frames
    frames = list(range(resume_from_frame, stop_frame_val, step_frame))
    
    print(f"Processing {len(frames)} frames from {resume_from_frame} to {stop_frame_val-1}, step {step_frame}")
    print(f"Using batch size: {batch_size}")
    print(f"Parallel processing: {'Enabled' if force_parallel else 'Disabled'}")
    print(f"Minimum cores: {min_cores}")
    
    # Batch processing
    frame_batches = [frames[i:i+batch_size] for i in range(0, len(frames), batch_size)]
    print(f"Created {len(frame_batches)} batches")
    
    # Process each batch
    results = []
    
    # Checkpoint interval (number of frames)
    checkpoint_interval = 200
    next_checkpoint_frame = frames[0] + checkpoint_interval if frames else 0
    
    # Process leaflet information
    if selected_leaflet0 is None and len(frames) > 0:
        try:
            universe.trajectory[frames[0]]
            print("\nNo valid leaflet from checkpoint. Detecting leaflets for the first frame...")
            leaflet0, _ = identify_lipid_leaflets(universe)
            selected_leaflet0 = leaflet0
            print("Leaflet detection completed for the first frame")
            if selected_leaflet0:
                print(f"Detected {len(selected_leaflet0.residues)} residues in upper leaflet")
        except Exception as e:
            print(f"Error detecting leaflets for the first frame: {str(e)}")
            traceback.print_exc()
    elif selected_leaflet0 is not None:
        print(f"\nUsing existing leaflet with {len(selected_leaflet0.residues)} residues")
    
    # Intermediate file creation flag
    created_intermediate_file = False
    
    for batch_idx, batch_frames in enumerate(frame_batches):
        print(f"\n----- Processing Batch {batch_idx+1}/{len(frame_batches)} -----")
        print(f"Frames: {batch_frames[0]} to {batch_frames[-1]}")
        
        # Display leaflet information status
        if selected_leaflet0 is not None:
            print(f"Using leaflet with {len(selected_leaflet0.residues)} residues for batch processing")
        else:
            print("Warning: No leaflet information available for batch processing!")
        
        # Batch processing (parallel or serial)
        batch_results = process_batch(
            batch_frames, 
            top_file, 
            traj_file,
            mp_context=mp_context,
            batch_size=batch_size,
            force_parallel=force_parallel,
            contact_cutoff=contact_cutoff,
            protein_cutoff=PROTEIN_CONTACT_CUTOFF,
            min_cores=min_cores,
            selected_leaflet0=selected_leaflet0
        )
        
        # Accumulate results
        if batch_results:
            # Process results from each frame and add to existing data
            for frame_result in batch_results:
                if not isinstance(frame_result, dict) or 'frame' not in frame_result:
                    print(f"Warning: Invalid frame result skipped: {type(frame_result)}")
                    continue
                    
                frame_idx = frame_result['frame']
                processed_frames.append(frame_idx)
                
                # Process protein-protein contacts
                if 'protein_contacts' in frame_result:
                    all_protein_contacts[frame_idx] = frame_result['protein_contacts']
                
                # Process lipid-protein contacts
                if 'lipid_contacts' in frame_result:
                    all_lipid_contacts[frame_idx] = frame_result['lipid_contacts']
                    
                # Process residue-level contact information
                if 'residue_contacts' in frame_result:
                    all_residue_contacts.append(frame_result['residue_contacts'])
                
                # NEW: Process target lipid binding states and lipid compositions
                if 'target_lipid_binding_state' in frame_result:
                    all_target_lipid_binding_states[frame_idx] = frame_result['target_lipid_binding_state']
                
                if 'lipid_composition_with_target_lipid' in frame_result:
                    for protein_name, composition in frame_result['lipid_composition_with_target_lipid'].items():
                        if protein_name not in all_lipid_compositions['with_target_lipid']:
                            all_lipid_compositions['with_target_lipid'][protein_name] = []
                        all_lipid_compositions['with_target_lipid'][protein_name].append({
                            'frame': frame_idx,
                            'composition': composition
                        })
                
                if 'lipid_composition_without_target_lipid' in frame_result:
                    for protein_name, composition in frame_result['lipid_composition_without_target_lipid'].items():
                        if protein_name not in all_lipid_compositions['without_target_lipid']:
                            all_lipid_compositions['without_target_lipid'][protein_name] = []
                        all_lipid_compositions['without_target_lipid'][protein_name].append({
                            'frame': frame_idx,
                            'composition': composition
                        })
                
                # NEW: Process unique molecule lipid compositions
                if 'unique_lipid_composition_with_target_lipid' in frame_result:
                    for protein_name, composition in frame_result['unique_lipid_composition_with_target_lipid'].items():
                        if protein_name not in all_unique_lipid_compositions['with_target_lipid']:
                            all_unique_lipid_compositions['with_target_lipid'][protein_name] = []
                        all_unique_lipid_compositions['with_target_lipid'][protein_name].append({
                            'frame': frame_idx,
                            'composition': composition
                        })
                
                if 'unique_lipid_composition_without_target_lipid' in frame_result:
                    for protein_name, composition in frame_result['unique_lipid_composition_without_target_lipid'].items():
                        if protein_name not in all_unique_lipid_compositions['without_target_lipid']:
                            all_unique_lipid_compositions['without_target_lipid'][protein_name] = []
                        all_unique_lipid_compositions['without_target_lipid'][protein_name].append({
                            'frame': frame_idx,
                            'composition': composition
                        })
                
                # Update leaflet information (only if not yet available)
                if selected_leaflet0 is None and 'leaflet0' in frame_result and frame_result['leaflet0'] is not None:
                    selected_leaflet0 = frame_result['leaflet0']
                    print("Updated selected leaflet from batch results")
                
                # Check checkpoint condition
                if frame_idx >= next_checkpoint_frame:
                    print(f"\n----- Creating Checkpoint at Frame {frame_idx} -----")
                    
                    # Create intermediate CSV file
                    created_intermediate_file = create_intermediate_csv(
                        all_protein_contacts, all_lipid_contacts, abs_output_dir
                    ) or created_intermediate_file
                    
                    # Save checkpoint with causal analysis data
                    save_checkpoint(
                        checkpoint_file, processed_frames, all_protein_contacts,
                        all_lipid_contacts, all_residue_contacts, selected_leaflet0,
                        all_target_lipid_binding_states, all_lipid_compositions,
                        all_unique_lipid_compositions
                    )
                    
                    # Set next checkpoint
                    next_checkpoint_frame = frame_idx + checkpoint_interval
            
            print(f"Batch {batch_idx+1} processed successfully, collected {len(batch_results)} frame results")
            
            # Create intermediate file after first batch if not yet created
            if batch_idx == 0 and not created_intermediate_file and processed_frames:
                print("\n----- Creating First Batch Checkpoint -----")
                created_intermediate_file = create_intermediate_csv(
                    all_protein_contacts, all_lipid_contacts, abs_output_dir
                ) or created_intermediate_file
                
                # Also save checkpoint after first batch with causal analysis data
                save_checkpoint(
                    checkpoint_file, processed_frames, all_protein_contacts,
                    all_lipid_contacts, all_residue_contacts, selected_leaflet0,
                    all_target_lipid_binding_states, all_lipid_compositions,
                    all_unique_lipid_compositions
                )
        else:
            print(f"Warning: Batch {batch_idx+1} returned no results")
    
    print(f"\nProcessed {len(processed_frames)} frames in total")
    
    # If results are empty
    if not processed_frames:
        print("No results collected. Analysis failed.")
        return None, None
    
    # Aggregate data across frames for complementarity analysis
    print("\n----- Aggregating Data Across Frames -----")
    
    # Aggregate protein-protein contacts
    aggregated_protein_contacts = {}
    for frame_idx, frame_contacts in all_protein_contacts.items():
        for pair_name, contact_data in frame_contacts.items():
            if pair_name not in aggregated_protein_contacts:
                aggregated_protein_contacts[pair_name] = {
                    'protein1': [],
                    'protein2': [],
                    'residue_ids1': None,
                    'residue_ids2': None,
                    'min_distances1': [],
                    'min_distances2': []
                }
            
            # Accumulate contact arrays
            if 'protein1' in contact_data:
                aggregated_protein_contacts[pair_name]['protein1'].append(contact_data['protein1'])
            if 'protein2' in contact_data:
                aggregated_protein_contacts[pair_name]['protein2'].append(contact_data['protein2'])
            
            # Store residue IDs (should be same across frames)
            if aggregated_protein_contacts[pair_name]['residue_ids1'] is None and 'residue_ids1' in contact_data:
                aggregated_protein_contacts[pair_name]['residue_ids1'] = contact_data['residue_ids1']
            if aggregated_protein_contacts[pair_name]['residue_ids2'] is None and 'residue_ids2' in contact_data:
                aggregated_protein_contacts[pair_name]['residue_ids2'] = contact_data['residue_ids2']
            
            # Accumulate min distances
            if 'min_distances1' in contact_data:
                aggregated_protein_contacts[pair_name]['min_distances1'].append(contact_data['min_distances1'])
            if 'min_distances2' in contact_data:
                aggregated_protein_contacts[pair_name]['min_distances2'].append(contact_data['min_distances2'])
    
    # Average the accumulated data
    import numpy as np
    for pair_name in aggregated_protein_contacts:
        if aggregated_protein_contacts[pair_name]['protein1']:
            aggregated_protein_contacts[pair_name]['protein1'] = np.mean(
                aggregated_protein_contacts[pair_name]['protein1'], axis=0
            )
        if aggregated_protein_contacts[pair_name]['protein2']:
            aggregated_protein_contacts[pair_name]['protein2'] = np.mean(
                aggregated_protein_contacts[pair_name]['protein2'], axis=0
            )
        if aggregated_protein_contacts[pair_name]['min_distances1']:
            aggregated_protein_contacts[pair_name]['min_distances1'] = np.mean(
                aggregated_protein_contacts[pair_name]['min_distances1'], axis=0
            )
        if aggregated_protein_contacts[pair_name]['min_distances2']:
            aggregated_protein_contacts[pair_name]['min_distances2'] = np.mean(
                aggregated_protein_contacts[pair_name]['min_distances2'], axis=0
            )
    
    # Aggregate lipid-protein contacts and restructure
    aggregated_lipid_contacts = {}
    for frame_idx, frame_lipid_contacts in all_lipid_contacts.items():
        for protein_name, protein_lipid_data in frame_lipid_contacts.items():
            for lipid_type, contact_data in protein_lipid_data.items():
                # Initialize structure: lipid_type -> protein_name
                if lipid_type not in aggregated_lipid_contacts:
                    aggregated_lipid_contacts[lipid_type] = {}
                if protein_name not in aggregated_lipid_contacts[lipid_type]:
                    aggregated_lipid_contacts[lipid_type][protein_name] = {
                        'contacts': [],
                        'residue_ids': None
                    }
                
                # Accumulate contacts
                if 'contacts' in contact_data:
                    aggregated_lipid_contacts[lipid_type][protein_name]['contacts'].append(contact_data['contacts'])
                
                # Store residue IDs (should be same across frames)
                if aggregated_lipid_contacts[lipid_type][protein_name]['residue_ids'] is None and 'residue_ids' in contact_data:
                    aggregated_lipid_contacts[lipid_type][protein_name]['residue_ids'] = contact_data['residue_ids']
    
    # Average the lipid contacts
    for lipid_type in aggregated_lipid_contacts:
        for protein_name in aggregated_lipid_contacts[lipid_type]:
            if aggregated_lipid_contacts[lipid_type][protein_name]['contacts']:
                aggregated_lipid_contacts[lipid_type][protein_name]['contacts'] = np.mean(
                    aggregated_lipid_contacts[lipid_type][protein_name]['contacts'], axis=0
                )
    
    # Save final results
    complementarity_df = analyze_contact_complementarity(aggregated_protein_contacts, aggregated_lipid_contacts)
    
    if complementarity_df is not None:
        complementarity_csv_path = os.path.join(abs_output_dir, "contact_complementarity.csv")
        complementarity_df.to_csv(complementarity_csv_path, index=False)
        print(f"Saved contact complementarity data to {complementarity_csv_path}")
    
    # NEW: Save target lipid binding states and lipid compositions for causal analysis
    print(f"\n----- Debug: Checking causal data before saving -----")
    print(f"Target lipid binding states collected: {len(all_target_lipid_binding_states)} frames")
    print(f"Lipid compositions collected: with_target_lipid={len(all_lipid_compositions.get('with_target_lipid', {}))}, without_target_lipid={len(all_lipid_compositions.get('without_target_lipid', {}))}")
    
    # Always save causal data regardless of binding states
    print(f"Target lipid binding states collected: {len(all_target_lipid_binding_states)} frames")
    if all_target_lipid_binding_states:
        print(f"Sample binding states: {list(all_target_lipid_binding_states.items())[:3]}")
    
    # Save causal data - this should work even if no target lipid binding occurred
    save_target_lipid_causal_data(all_target_lipid_binding_states, all_lipid_compositions, abs_output_dir, TARGET_LIPID, all_unique_lipid_compositions)
    
    # Generate plots
    if complementarity_df is not None and not debug_mode:
        try:
            # Check if target lipid column exists
            target_lipid_col = 'DPG3_contact'  # Default target lipid
            has_target = target_lipid_col in complementarity_df.columns and complementarity_df[target_lipid_col].max() > 0
            plot_contact_complementarity(complementarity_df, abs_output_dir, with_target_lipid=has_target)
        except Exception as e:
            print(f"Error generating plots: {str(e)}")
            traceback.print_exc()
    
    # Aggregate residue-level contact data
    residue_df = None
    if all_residue_contacts:
        print("\n----- Aggregating Residue-Level Contact Data -----")
        residue_df = aggregate_residue_contacts(results)
        
        # Save results as CSV
        if residue_df is not None:
            residue_csv_path = os.path.join(abs_output_dir, "residue_contacts.csv")
            residue_df.to_csv(residue_csv_path, index=False)
            print(f"Saved residue contact data to {residue_csv_path}")
    
    # Display processing time
    end_time = time.time()
    execution_time = end_time - start_time
    hours = int(execution_time // 3600)
    minutes = int((execution_time % 3600) // 60)
    seconds = execution_time % 60
    
    print(f"\n----- Analysis Complete -----")
    print(f"Results saved to {abs_output_dir}")
    print(f"Total execution time: {hours}h {minutes}m {seconds:.1f}s")
    
    # Remove checkpoint file (on successful completion)
    try:
        if os.path.exists(checkpoint_file):
            os.remove(checkpoint_file)
            print(f"Checkpoint file removed after successful completion")
    except Exception as e:
        print(f"Warning: Could not remove checkpoint file: {str(e)}")
    
    return complementarity_df, residue_df

def main():
    """Main execution function - parse command line arguments and run analysis"""
    
    # Execution start time
    start_time = time.time()
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    print(f"===== Contact Analysis Stage 1 =====")
    print(f"Started at: {timestamp}")
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Protein-protein and lipid-protein contact analysis')
    
    # File paths
    parser.add_argument('--with-target-lipid-psf', default=WITH_TARGET_LIPID_PSF, help='PSF file for system with target lipid')
    parser.add_argument('--with-target-lipid-xtc', default=WITH_TARGET_LIPID_XTC, help='Trajectory file for system with target lipid')
    parser.add_argument('--without-target-lipid-psf', default=WITHOUT_TARGET_LIPID_PSF, help='PSF file for system without target lipid')
    parser.add_argument('--without-target-lipid-xtc', default=WITHOUT_TARGET_LIPID_XTC, help='Trajectory file for system without target lipid')
    parser.add_argument('--with-target-lipid-output', default=WITH_TARGET_LIPID_OUTPUT, help='Output directory for system with target lipid')
    parser.add_argument('--without-target-lipid-output', default=WITHOUT_TARGET_LIPID_OUTPUT, help='Output directory for system without target lipid')
    parser.add_argument('--comparison-dir', default=DEFAULT_COMPARISON_OUTPUT, help='Output directory for comparison results')
    
    # Analysis parameters
    parser.add_argument('--start', type=int, default=START_FRAME, help='Starting frame')
    parser.add_argument('--stop', type=int, default=STOP_FRAME, help='Ending frame')
    parser.add_argument('--step', type=int, default=STEP_FRAME, help='Frame interval')
    parser.add_argument('--batch-size', type=int, default=BATCH_SIZE, help='Batch size')
    parser.add_argument('--cores', type=int, default=MIN_CORES, help='Minimum number of cores')
    parser.add_argument('--dimer-cutoff', type=float, default=DIMER_CUTOFF, help='Dimer detection cutoff distance')
    parser.add_argument('--contact-cutoff', type=float, default=CONTACT_CUTOFF, help='Contact detection cutoff distance')
    parser.add_argument('--leaflet-cutoff', type=float, default=LEAFLET_CUTOFF, help='Leaflet detection cutoff distance')
    
    # Optional flags
    parser.add_argument('--debug', action='store_true', help='Debug mode')
    parser.add_argument('--skip-with-target-lipid', action='store_true', help='Skip analysis of system with target lipid')
    parser.add_argument('--skip-without-target-lipid', action='store_true', help='Skip analysis of system without target lipid')
    parser.add_argument('--skip-plots', action='store_true', help='Skip plot generation')
    parser.add_argument('--no-parallel', dest='force_parallel', action='store_false', help='Disable parallel processing')
    parser.add_argument('--frames', type=int, default=0, help='Number of frames to process (0=all)')
    parser.add_argument('--mp-context', choices=['fork', 'spawn'], help='Multiprocessing context')
    parser.add_argument('--clear-leaflet', action='store_true', help='Remove existing leaflet information files')
    parser.add_argument('--force-leaflet-detection', action='store_true', help='Force leaflet detection (ignore checkpoint)')
    
    # Set default values
    parser.set_defaults(force_parallel=True)
    
    args = parser.parse_args()
    
    # Clear leaflet files if specified
    if args.clear_leaflet:
        leaflet_info_file = os.path.join(TEMP_FILES_DIR, "leaflet_info.pickle")
        leaflet_info_txt = os.path.join(TEMP_FILES_DIR, "leaflet_info.txt")
        
        if os.path.exists(leaflet_info_file):
            print(f"Removing leaflet info file: {leaflet_info_file}")
            os.remove(leaflet_info_file)
        
        if os.path.exists(leaflet_info_txt):
            print(f"Removing leaflet text report: {leaflet_info_txt}")
            os.remove(leaflet_info_txt)
    
    # Clear checkpoint files if forcing leaflet detection
    if args.force_leaflet_detection:
        print("Force leaflet detection mode - clearing checkpoint files")
        abs_with_target_lipid_output = os.path.abspath(args.with_target_lipid_output)
        abs_without_target_lipid_output = os.path.abspath(args.without_target_lipid_output)
        
        with_target_lipid_checkpoint = os.path.join(abs_with_target_lipid_output, "analysis_checkpoint.pkl")
        without_target_lipid_checkpoint = os.path.join(abs_without_target_lipid_output, "analysis_checkpoint.pkl")
        
        if os.path.exists(with_target_lipid_checkpoint):
            print(f"Removing with_target_lipid checkpoint: {with_target_lipid_checkpoint}")
            os.remove(with_target_lipid_checkpoint)
        
        if os.path.exists(without_target_lipid_checkpoint):
            print(f"Removing without_target_lipid checkpoint: {without_target_lipid_checkpoint}")
            os.remove(without_target_lipid_checkpoint)
    
    # Prepare output directories
    abs_with_target_lipid_output = os.path.abspath(args.with_target_lipid_output)
    abs_without_target_lipid_output = os.path.abspath(args.without_target_lipid_output)
    abs_comparison_dir = os.path.abspath(args.comparison_dir)
    
    os.makedirs(abs_with_target_lipid_output, exist_ok=True)
    os.makedirs(abs_without_target_lipid_output, exist_ok=True)
    os.makedirs(abs_comparison_dir, exist_ok=True)
    
    # Display execution settings
    print(f"\n[Execution Settings]")
    print(f"- Minimum cores: {args.cores}")
    print(f"- Parallel processing: {'Enabled' if args.force_parallel else 'Disabled'}")
    print(f"- Multiprocessing context: {args.mp_context if args.mp_context else 'Auto-select'}")
    print(f"- Batch size: {args.batch_size}")
    print(f"- Leaflet cutoff: {args.leaflet_cutoff} Å")
    
    # Set frame range
    if args.frames > 0:
        if args.stop > 0:
            args.stop = args.start + args.frames * args.step
            print(f"- Frame count specified, setting end frame to {args.stop}")
    
    print(f"- Frame range: {args.start} to {args.stop}, interval {args.step}")
    print(f"- Contact cutoff distance: {args.contact_cutoff} Å")
    print(f"- Dimer detection cutoff distance: {args.dimer_cutoff} Å")
    
    # Display directory information
    print(f"\nOutput directories:")
    print(f"  With {TARGET_LIPID}: {abs_with_target_lipid_output}")
    print(f"  Without {TARGET_LIPID}: {abs_without_target_lipid_output}")
    print(f"  Comparison results: {abs_comparison_dir}")
    
    # Initialize dataframes
    with_target_lipid_df = None
    without_target_lipid_df = None
    with_target_lipid_residue_df = None
    without_target_lipid_residue_df = None
    
    # Analyze system with target lipid
    if not args.skip_with_target_lipid:
        print(f"\n===== Analyzing System With {TARGET_LIPID} =====")
        with_target_lipid_df, with_target_lipid_residue_df = run_analysis(
            args.with_target_lipid_psf,
            args.with_target_lipid_xtc,
            abs_with_target_lipid_output,
            start_frame=args.start,
            stop_frame=args.stop,
            step_frame=args.step,
            batch_size=args.batch_size,
            debug_mode=args.debug,
            force_parallel=args.force_parallel,
            mp_context=args.mp_context,
            min_cores=args.cores,
            contact_cutoff=args.contact_cutoff
        )
        
        # Generate plots for with target lipid system
        if not args.skip_plots and with_target_lipid_df is not None:
            try:
                # Regular contact complementarity plots
                has_target_lipid = f'{TARGET_LIPID}_contact' in with_target_lipid_df.columns and with_target_lipid_df[f'{TARGET_LIPID}_contact'].max() > 0
                plot_contact_complementarity(with_target_lipid_df, abs_with_target_lipid_output, with_target_lipid=has_target_lipid)
                
                # Special plots for target lipid system
                for protein in with_target_lipid_df['protein'].unique():
                    # Target lipid competition plots
                    plot_target_lipid_protein_competition(with_target_lipid_df, abs_with_target_lipid_output, protein)
                
                # Protein-protein residue contact plots
                if with_target_lipid_residue_df is not None:
                    try:
                        print(f"\nSpecial plot: Protein-protein residue contacts (with {TARGET_LIPID})")
                        plot_protein_protein_residue_contacts(with_target_lipid_residue_df, abs_with_target_lipid_output)
                    except Exception as e:
                        print(f"Error in protein-protein residue contact plot: {str(e)}")
                        traceback.print_exc()
            except Exception as e:
                print(f"Error generating plots: {str(e)}")
                traceback.print_exc()
    
    # Analyze system without target lipid
    if not args.skip_without_target_lipid:
        print(f"\n===== Analyzing System Without {TARGET_LIPID} =====")
        without_target_lipid_df, without_target_lipid_residue_df = run_analysis(
            args.without_target_lipid_psf,
            args.without_target_lipid_xtc,
            abs_without_target_lipid_output,
            start_frame=args.start,
            stop_frame=args.stop,
            step_frame=args.step,
            batch_size=args.batch_size,
            debug_mode=args.debug,
            force_parallel=args.force_parallel,
            mp_context=args.mp_context,
            min_cores=args.cores,
            contact_cutoff=args.contact_cutoff
        )
        
        # Generate plots for without target lipid system
        if not args.skip_plots and without_target_lipid_df is not None:
            try:
                # Regular contact complementarity plots
                plot_contact_complementarity(without_target_lipid_df, abs_without_target_lipid_output, with_target_lipid=False)
                
                # Protein-protein residue contact plots
                if without_target_lipid_residue_df is not None:
                    try:
                        print(f"\nSpecial plot: Protein-protein residue contacts (without {TARGET_LIPID})")
                        plot_protein_protein_residue_contacts(without_target_lipid_residue_df, abs_without_target_lipid_output)
                    except Exception as e:
                        print(f"Error in protein-protein residue contact plot: {str(e)}")
                        traceback.print_exc()
            except Exception as e:
                print(f"Error generating plots: {str(e)}")
                traceback.print_exc()
    
    # Comparison analysis between with/without target lipid
    if not args.skip_plots and with_target_lipid_df is not None and without_target_lipid_df is not None:
        print("\n===== Comparison Analysis Between With/Without Target Lipid =====")
        
        try:
            # Check for target lipid data
            target_lipid_col = f'{TARGET_LIPID}_contact'
            if target_lipid_col in with_target_lipid_df.columns:
                target_max = with_target_lipid_df[target_lipid_col].max()
                target_nonzero = (with_target_lipid_df[target_lipid_col] > 0).sum()
                print(f"Target lipid contact data for comparison: max = {target_max}, positive values = {target_nonzero}")
                
                # Run comprehensive comparison analysis
                compare_with_without_target_lipid(
                    with_target_lipid_df, 
                    without_target_lipid_df, 
                    abs_comparison_dir, 
                    target_lipid_col=target_lipid_col
                )
                
                # Also generate individual protein comparison plots
                for protein in with_target_lipid_df['protein'].unique():
                    if protein in without_target_lipid_df['protein'].unique():
                        # Target lipid competition plots with comparison
                        plot_target_lipid_protein_competition(with_target_lipid_df, abs_comparison_dir, protein, without_target_lipid_df)
            else:
                print("WARNING: No target lipid contact data in with_lipid_df. Skipping comparison analysis.")
        except Exception as e:
            print(f"Error in comparison analysis: {str(e)}")
            traceback.print_exc()
    
    # Calculate processing time
    end_time = time.time()
    elapsed_time = end_time - start_time
    hours, remainder = divmod(elapsed_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    
    print(f"\n===== Analysis Complete =====")
    print(f"Total processing time: {int(hours)}h {int(minutes)}m {seconds:.1f}s")
    print(f"Results saved to:")
    print(f"  With {TARGET_LIPID}: {abs_with_target_lipid_output}")
    print(f"  Without {TARGET_LIPID}: {abs_without_target_lipid_output}")
    print(f"  Comparison results: {abs_comparison_dir}")
    
    return 0

if __name__ == "__main__":
    main()