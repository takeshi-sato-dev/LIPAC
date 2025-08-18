"""Parallel processing utilities"""

import os
import pickle
import traceback
import numpy as np
import multiprocessing as mp
from multiprocessing import Pool, cpu_count
import MDAnalysis as mda

from ..core.trajectory_loader import select_lipids, select_proteins, identify_lipid_leaflets
from ..core.contact_calculator import (
    calculate_protein_protein_contacts,
    calculate_lipid_protein_contacts,
    calculate_unique_lipid_protein_contacts,
    calculate_protein_com_distances
)
from ..analysis.residue_contacts import extract_residue_contacts
from ..config import CONTACT_CUTOFF, PROTEIN_CONTACT_CUTOFF, TEMP_FILES_DIR, TARGET_LIPID

def test_multiprocessing():
    """Test multiprocessing functionality
    
    Returns
    -------
    tuple
        (success, context_method)
    """
    print("\n===== Testing Multiprocessing Functionality =====")
    try:
        # Check available CPU cores
        n_cores = cpu_count()
        print(f"System has {n_cores} CPU cores available")
        
        # Display platform information
        import platform
        print(f"Operating System: {platform.system()} {platform.release()}")
        print(f"Python Version: {platform.python_version()}")
        print(f"Machine: {platform.machine()}")
        
        # Detect Apple Silicon
        is_apple_silicon = False
        if platform.system() == 'Darwin' and platform.machine() == 'arm64':
            is_apple_silicon = True
            print("Detected Apple Silicon (M1/M2/M3)")
        
        # Select optimal context
        if platform.system() == 'Darwin':
            ctx_method = 'fork'  # Use 'fork' on macOS
            print("Using 'fork' context for macOS")
        else:
            ctx_method = 'spawn'  # Use 'spawn' on other OS
            print("Using 'spawn' context for non-macOS platform")
        
        # Simple test function
        def square(x):
            import os
            import time
            # Display process ID to confirm execution in separate process
            print(f"Process {os.getpid()} computing square of {x}")
            # Add short sleep to make CPU usage more apparent
            time.sleep(0.5)
            return x * x
        
        # Explicitly start processes and check for setup issues
        print("Starting test pool with 4 processes (or max available)")
        n_test_processes = min(4, n_cores)
        
        # Explicitly set context
        ctx = mp.get_context(ctx_method)
        with ctx.Pool(processes=n_test_processes) as pool:
            test_data = [1, 2, 3, 4]
            print(f"Mapping {test_data} to worker processes...")
            
            # Start execution in worker processes
            print("Starting parallel execution...")
            
            # Measure start time
            import time
            start_time = time.time()
            
            # Execute parallel processing
            results = pool.map(square, test_data)
            
            # Measure end time
            end_time = time.time()
            duration = end_time - start_time
            
            print(f"Parallel execution completed in {duration:.2f} seconds")
            print(f"Results: {results}")
            
            # Check if multiple CPUs are actually being used
            expected_duration = 0.5 * len(test_data) / n_test_processes
            parallel_speedup = (0.5 * len(test_data)) / duration
            
            print(f"Expected sequential duration: ~{0.5 * len(test_data):.2f} seconds")
            print(f"Expected parallel duration: ~{expected_duration:.2f} seconds")
            print(f"Actual duration: {duration:.2f} seconds")
            print(f"Approximate parallel speedup: {parallel_speedup:.2f}x")
            
            if parallel_speedup > 1.5:
                print("PASSED: Multiple CPUs are being utilized!")
            else:
                print("WARNING: Parallel execution doesn't show expected speedup. May be using only one CPU.")
            
        print("Multiprocessing test completed successfully!")
        return True, ctx_method
    except Exception as e:
        print(f"Error during multiprocessing test: {str(e)}")
        traceback.print_exc()
        return False, 'fork'  # Return 'fork' as default

def _frame_processor_worker(args):
    """Worker function for multiprocessing - processes one frame
    
    Parameters
    ----------
    args : tuple
        Arguments for frame processing
        
    Returns
    -------
    dict
        Processing results
    """
    try:
        # Get process ID
        worker_pid = os.getpid()
        
        # Unpack arguments (including leaflet information)
        if len(args) >= 7:  # If leaflet information is included
            frame_idx, top_file, traj_file, contact_cutoff, protein_cutoff, is_first, selected_leaflet0 = args
        else:  # For backward compatibility
            frame_idx, top_file, traj_file, contact_cutoff, protein_cutoff, is_first = args
            selected_leaflet0 = None
        
        print(f"Worker {worker_pid}: Processing frame {frame_idx}")
                
        # Load trajectory
        universe = mda.Universe(top_file, traj_file)
        
        # Load frame
        universe.trajectory[frame_idx]
        box = universe.dimensions[:3]
        print(f"Worker {worker_pid}: Loaded frame {frame_idx}, box dimensions: {box}")
        
        # Get leaflet information
        leaflet0 = selected_leaflet0  # Use passed leaflet information
        
        # If no leaflet information, load or generate
        if leaflet0 is None:
            # Path to leaflet information file
            leaflet_info_file = os.path.join(TEMP_FILES_DIR, "leaflet_info.pickle")
            
            # Try to load from file
            if os.path.exists(leaflet_info_file):
                print(f"Worker {worker_pid}: Loading leaflet information from {leaflet_info_file}")
                try:
                    with open(leaflet_info_file, 'rb') as f:
                        leaflet_info = pickle.load(f)
                    
                    # Get upper leaflet lipid resids
                    upper_leaflet_resids = leaflet_info['upper_leaflet_resids']
                    
                    # Reconstruct leaflet
                    leaflet0 = mda.AtomGroup([], universe)
                    
                    # Select only lipid atoms first
                    lipid_atoms = universe.select_atoms("resname CHOL DPSM DIPC DPG3 DOPS")
                    
                    # Reconstruct leaflet with batch processing
                    batch_size = 100
                    for i in range(0, len(upper_leaflet_resids), batch_size):
                        batch = upper_leaflet_resids[i:i+batch_size]
                        batch_sel_str = " or ".join([f"resid {r}" for r in batch])
                        try:
                            # Select from lipid atoms
                            batch_atoms = lipid_atoms.select_atoms(batch_sel_str)
                            leaflet0 = leaflet0.union(batch_atoms)
                        except Exception as e:
                            print(f"Worker {worker_pid}: Error selecting leaflet batch: {str(e)}")
                    
                    print(f"Worker {worker_pid}: Leaflet reconstructed with {len(leaflet0.residues)} residues")
                except Exception as e:
                    print(f"Worker {worker_pid}: Error loading leaflet info: {str(e)}")
                    # If error, detect leaflet on the spot
                    leaflet0, _ = identify_lipid_leaflets(universe)
            else:
                # If file doesn't exist, detect leaflet
                print(f"Worker {worker_pid}: No leaflet info file found, detecting leaflets")
                leaflet0, _ = identify_lipid_leaflets(universe)
        
        # Select lipids
        lipid_sels = select_lipids(universe, leaflet0)
        
        # Select proteins
        proteins = select_proteins(universe)
        
        # Detect protein pairs
        close_pairs = calculate_protein_com_distances(universe, proteins)
        
        # Calculate protein-protein contacts
        protein_contacts = {}
        for pair_name in close_pairs:
            protein1_name, protein2_name = pair_name.split('-')
            protein1 = proteins[protein1_name]
            protein2 = proteins[protein2_name]
            
            print(f"Worker {worker_pid}: Calculating contacts between {protein1_name} and {protein2_name}")
            
            p1_contacts, p2_contacts, contact_matrix, p1_min_dist, p2_min_dist, residue_ids1, residue_ids2 = calculate_protein_protein_contacts(
                protein1, protein2, box, cutoff=protein_cutoff
            )
            
            protein_contacts[pair_name] = {
                'protein1': p1_contacts,
                'protein2': p2_contacts,
                'contact_matrix': contact_matrix,
                'residue_ids1': residue_ids1,
                'residue_ids2': residue_ids2,
                'min_distances1': p1_min_dist,
                'min_distances2': p2_min_dist
            }
        
        # Calculate lipid-protein contacts (both residue-level and unique molecule counts)
        lipid_contacts = {}
        unique_lipid_contacts = {}
        for protein_name, protein in proteins.items():
            if len(protein) == 0:
                continue
                
            print(f"Worker {worker_pid}: Calculating lipid-protein contacts for {protein_name}")
            lipid_contacts[protein_name] = calculate_lipid_protein_contacts(
                protein, lipid_sels, box, cutoff=contact_cutoff
            )
            unique_lipid_contacts[protein_name] = calculate_unique_lipid_protein_contacts(
                protein, lipid_sels, box, cutoff=contact_cutoff
            )
        
        # Extract residue-level contact information
        residue_contacts = extract_residue_contacts(universe, frame_idx, proteins, close_pairs, leaflet0)
        
        # Track target lipid binding state for each protein (CAUSAL DATA COLLECTION)
        target_lipid_binding_state = {}
        lipid_composition_with_target_lipid = {}
        lipid_composition_without_target_lipid = {}
        unique_lipid_composition_with_target_lipid = {}
        unique_lipid_composition_without_target_lipid = {}
        
        # DEBUG: Print available lipid types and target lipid
        print(f"Worker {worker_pid}: DEBUG: TARGET_LIPID = {TARGET_LIPID}")
        print(f"Worker {worker_pid}: DEBUG: Available lipid types in lipid_sels = {list(lipid_sels.keys()) if lipid_sels else 'None'}")
        if proteins:
            first_protein = list(proteins.keys())[0]
            print(f"Worker {worker_pid}: DEBUG: Available lipid types in lipid_contacts[{first_protein}] = {list(lipid_contacts[first_protein].keys()) if first_protein in lipid_contacts else 'None'}")
        
        for protein_name, protein in proteins.items():
            # Check if target lipid is in contact with this protein
            target_lipid_contact = False
            print(f"Worker {worker_pid}: DEBUG: Checking {protein_name} for {TARGET_LIPID} contacts...")
            print(f"Worker {worker_pid}: DEBUG: lipid_contacts[{protein_name}] keys = {list(lipid_contacts[protein_name].keys()) if protein_name in lipid_contacts else 'None'}")
            
            if TARGET_LIPID in lipid_contacts[protein_name]:
                target_lipid_contacts = lipid_contacts[protein_name][TARGET_LIPID]['contacts']
                contact_sum = np.sum(target_lipid_contacts)
                print(f"Worker {worker_pid}: DEBUG: {protein_name} - {TARGET_LIPID} contact sum = {contact_sum}")
                if contact_sum > 0:
                    target_lipid_contact = True
            else:
                print(f"Worker {worker_pid}: DEBUG: {TARGET_LIPID} not found in lipid_contacts[{protein_name}]")
            
            target_lipid_binding_state[protein_name] = target_lipid_contact
            
            # Track lipid composition based on target lipid binding state
            if target_lipid_contact:
                lipid_composition_with_target_lipid[protein_name] = {}
                unique_lipid_composition_with_target_lipid[protein_name] = {}
                for lipid_type in lipid_sels:
                    if lipid_type != TARGET_LIPID and lipid_type in lipid_contacts[protein_name]:
                        # Residue contact count
                        contact_count = np.sum(lipid_contacts[protein_name][lipid_type]['contacts'])
                        lipid_composition_with_target_lipid[protein_name][lipid_type] = contact_count
                        # Unique molecule count
                        unique_count = unique_lipid_contacts[protein_name][lipid_type]
                        unique_lipid_composition_with_target_lipid[protein_name][lipid_type] = unique_count
            else:
                lipid_composition_without_target_lipid[protein_name] = {}
                unique_lipid_composition_without_target_lipid[protein_name] = {}
                for lipid_type in lipid_sels:
                    if lipid_type != TARGET_LIPID and lipid_type in lipid_contacts[protein_name]:
                        # Residue contact count
                        contact_count = np.sum(lipid_contacts[protein_name][lipid_type]['contacts'])
                        lipid_composition_without_target_lipid[protein_name][lipid_type] = contact_count
                        # Unique molecule count
                        unique_count = unique_lipid_contacts[protein_name][lipid_type]
                        unique_lipid_composition_without_target_lipid[protein_name][lipid_type] = unique_count
        
        print(f"Worker {worker_pid}: Completed processing frame {frame_idx}")
        
        # Return successful result (including leaflet information and causal data)
        return {
            'success': True,
            'frame': frame_idx,
            'protein_contacts': protein_contacts,
            'lipid_contacts': lipid_contacts,
            'unique_lipid_contacts': unique_lipid_contacts,  # NEW: Unique molecule counts
            'residue_contacts': residue_contacts,
            'close_pairs': close_pairs,
            'leaflet0': leaflet0,  # Include leaflet information
            'target_lipid_binding_state': target_lipid_binding_state,  # NEW: Target lipid binding state
            'lipid_composition_with_target_lipid': lipid_composition_with_target_lipid,  # NEW: Composition when target lipid bound
            'lipid_composition_without_target_lipid': lipid_composition_without_target_lipid,  # NEW: Composition when target lipid not bound
            'unique_lipid_composition_with_target_lipid': unique_lipid_composition_with_target_lipid,  # NEW: Unique molecule composition when target bound
            'unique_lipid_composition_without_target_lipid': unique_lipid_composition_without_target_lipid  # NEW: Unique molecule composition when target not bound
        }
        
    except Exception as e:
        error_str = traceback.format_exc()
        print(f"Error in worker processing frame {args[0]}: {str(e)}")
        print(error_str)
        
        # Return failed result
        return {
            'success': False,
            'frame': args[0],
            'error': str(e),
            'traceback': error_str
        }

def process_batch(frames, top_file, traj_file, mp_context='fork', batch_size=50, min_cores=2, 
                 contact_cutoff=CONTACT_CUTOFF, protein_cutoff=PROTEIN_CONTACT_CUTOFF, 
                 force_parallel=False, selected_leaflet0=None):
    """Process frame batch in parallel
    
    Parameters
    ----------
    frames : list
        List of frame indices
    top_file : str
        Path to topology file
    traj_file : str
        Path to trajectory file
    mp_context : str
        Multiprocessing context
    batch_size : int
        Batch size
    min_cores : int
        Minimum number of cores
    contact_cutoff : float
        Lipid-protein contact cutoff
    protein_cutoff : float
        Protein-protein contact cutoff
    force_parallel : bool
        Force parallel processing
    selected_leaflet0 : MDAnalysis.AtomGroup, optional
        Selected leaflet
        
    Returns
    -------
    list
        Processing results
    """
    print("\n===== Processing Batch =====")
    print(f"Batch contains {len(frames)} frames")
    
    # Ensure temp files directory
    os.makedirs(TEMP_FILES_DIR, exist_ok=True)
    
    # Handle leaflet information
    leaflet_info_file = os.path.join(TEMP_FILES_DIR, "leaflet_info.pickle")
    leaflet_info_txt = os.path.join(TEMP_FILES_DIR, "leaflet_info.txt")
    
    regenerate_leaflet = False
    
    # Consider regeneration only if no selected leaflet
    if selected_leaflet0 is None:
        if os.path.exists(leaflet_info_file):
            print(f"Removing existing leaflet information file: {leaflet_info_file}")
            os.remove(leaflet_info_file)
            regenerate_leaflet = True
        
        if os.path.exists(leaflet_info_txt):
            print(f"Removing existing leaflet text report: {leaflet_info_txt}")
            os.remove(leaflet_info_txt)
    else:
        print("Using provided leaflet information - not regenerating")
    
    # Generate or reuse leaflet information
    leaflet0 = selected_leaflet0
    
    # Need to regenerate leaflet information
    if leaflet0 is None or regenerate_leaflet:
        try:
            # Load trajectory
            from ..core.trajectory_loader import load_universe
            universe = load_universe(top_file, traj_file)
            if universe:
                universe.trajectory[frames[0]]  # Load first frame
                # Detect and save leaflets
                print("Forcing regeneration of leaflet information...")
                leaflet0, _ = identify_lipid_leaflets(universe)
                if os.path.exists(leaflet_info_file):
                    print(f"Successfully created new leaflet information file: {leaflet_info_file}")
                else:
                    print(f"Failed to create leaflet information file")
                    return []
            else:
                print("Could not load universe to create leaflet information")
                return []
        except Exception as e:
            print(f"Error creating leaflet information: {str(e)}")
            traceback.print_exc()
            return []
    
    # Determine available CPU cores
    n_available_cores = mp.cpu_count()
    print(f"System has {n_available_cores} available CPU cores")
    
    # Determine number of cores (use at least MIN_CORES)
    num_cores = max(min_cores, min(int(n_available_cores * 0.75), len(frames)))
    
    # Check if parallel processing is possible (force_parallel=True forces parallel)
    parallel_possible = (num_cores > 1 and len(frames) > 1) or force_parallel
    
    if parallel_possible:
        print(f"Using {num_cores} cores for parallel processing")
        
        # Pass trajectory file information and calculation parameters to each worker process
        worker_args = []
        for i, frame_idx in enumerate(frames):
            # Set is_first=True for first frame
            is_first = (i == 0)
            # Pass necessary parameters to each worker process
            worker_args.append((
                frame_idx, 
                top_file, 
                traj_file, 
                contact_cutoff,  # Lipid-protein contact cutoff
                protein_cutoff,  # Protein-protein contact cutoff
                is_first,
                selected_leaflet0  # Also pass selected leaflet
            ))
        
        # Execute multiprocessing
        try:
            # Use 'fork' for macOS, 'spawn' for others
            ctx = mp.get_context(mp_context)
            
            print(f"Creating process pool with {num_cores} workers using {mp_context} context")
            print("Each worker will perform COMPLETE contact analysis independently")
            
            with ctx.Pool(processes=num_cores) as pool:
                # Execute processing for each frame in parallel (complete calculation in each process)
                results = pool.map(_frame_processor_worker, worker_args)
                
                # Keep only successful frames as results
                successful_results = [res for res in results if res['success']]
                
                print(f"Successfully processed {len(successful_results)}/{len(frames)} frames in parallel")
                
                return successful_results
                
        except Exception as e:
            print(f"Error in parallel processing: {str(e)}")
            print("Error details:")
            traceback.print_exc()
            print("Falling back to serial processing...")
            
            # Load trajectory
            from ..core.trajectory_loader import load_universe
            universe = load_universe(top_file, traj_file)
            if universe is None:
                return []
                
            # Select proteins
            proteins = select_proteins(universe)
            
            # Execute serial processing
            from ..core.frame_processor import process_frames_serially
            results, _ = process_frames_serially(universe, proteins, frames, leaflet0)
            return results
    else:
        # Serial processing if parallel not possible
        print("Using serial processing (single core or small batch)")
        
        # Load trajectory
        from ..core.trajectory_loader import load_universe
        universe = load_universe(top_file, traj_file)
        if universe is None:
            return []
            
        # Select proteins
        proteins = select_proteins(universe)
        
        # Execute serial processing (pass leaflet information)
        from ..core.frame_processor import process_frames_serially
        results, _ = process_frames_serially(universe, proteins, frames, leaflet0)
        return results