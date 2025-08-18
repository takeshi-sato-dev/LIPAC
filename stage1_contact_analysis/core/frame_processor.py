"""Frame processing functions for trajectory analysis"""

import numpy as np
import traceback
from .trajectory_loader import select_lipids, select_proteins, identify_lipid_leaflets
from .contact_calculator import (
    calculate_protein_protein_contacts,
    calculate_lipid_protein_contacts,
    calculate_protein_com_distances
)
from ..analysis.residue_contacts import extract_residue_contacts
from ..config import PROTEIN_CONTACT_CUTOFF, CONTACT_CUTOFF, TARGET_LIPID

def process_frame(universe, frame_idx, proteins, reference_leaflet0=None, force_update_leaflets=False):
    """Process one frame
    
    Parameters
    ----------
    universe : MDAnalysis.Universe
        The universe object
    frame_idx : int
        Frame index to process
    proteins : dict
        Dictionary of protein selections
    reference_leaflet0 : MDAnalysis.AtomGroup, optional
        Reference leaflet to use
    force_update_leaflets : bool
        Force leaflet update
        
    Returns
    -------
    dict
        Processing results for the frame
    """
    try:
        print(f"\n==== Processing frame {frame_idx} ====")
        
        # Load frame
        universe.trajectory[frame_idx]
        box = universe.dimensions[:3]
        print(f"Frame {frame_idx} loaded, box dimensions: {box}")
        
        # Detect leaflets (only use LeafletFinder for first frame or when forced)
        if force_update_leaflets or reference_leaflet0 is None:
            print(f"Frame {frame_idx}: Performing leaflet detection with LeafletFinder")
            leaflet0, leaflet1 = identify_lipid_leaflets(universe)
            # Update reference leaflet
            reference_leaflet0 = leaflet0
        else:
            # Use previous leaflet as is
            print(f"Frame {frame_idx}: Using reference leaflet from previous frame")
            leaflet0 = reference_leaflet0
            leaflet1 = None  # Lower leaflet not used
        
        # Select lipids - use only leaflet0 (upper leaflet)
        lipid_sels = select_lipids(universe, leaflet0)
        
        # Calculate COM distances between proteins and get close pairs
        print("Calculating protein-protein COM distances...")
        close_pairs = calculate_protein_com_distances(universe, proteins)
        print(f"Found {len(close_pairs)} close protein pairs")
        
        # Calculate protein-protein contacts
        print("Calculating protein-protein contacts...")
        protein_contacts = {}
        for pair_name in close_pairs:
            protein1_name, protein2_name = pair_name.split('-')
            protein1 = proteins[protein1_name]
            protein2 = proteins[protein2_name]
            
            print(f"  Calculating contacts between {protein1_name} and {protein2_name}...")
            p1_contacts, p2_contacts, contact_matrix, p1_min_dist, p2_min_dist, residue_ids1, residue_ids2 = calculate_protein_protein_contacts(
                protein1, protein2, box, cutoff=PROTEIN_CONTACT_CUTOFF
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
            
            p1_total = np.sum(p1_contacts)
            p2_total = np.sum(p2_contacts)
            print(f"  Found {p1_total} contacts for {protein1_name} and {p2_total} contacts for {protein2_name}")
        
        # Calculate lipid-protein contacts
        print("Calculating lipid-protein contacts (3D distance)...")
        lipid_contacts = {}
        for protein_name, protein in proteins.items():
            print(f"  Processing {protein_name}...")
            lipid_contacts[protein_name] = calculate_lipid_protein_contacts(
                protein, lipid_sels, box, cutoff=CONTACT_CUTOFF
            )
            
            # Display contact counts for each lipid type
            for lipid_type in lipid_sels:
                if protein_name in lipid_contacts and lipid_type in lipid_contacts[protein_name]:
                    contact_values = lipid_contacts[protein_name][lipid_type]['contacts']
                    contact_count = np.sum(contact_values)
                    print(f"  {protein_name} has {contact_count} contacts with {lipid_type}")
        
        # Get residue-level contact information
        residue_contacts = extract_residue_contacts(universe, frame_idx, proteins, close_pairs, leaflet0)
        
        print(f"Frame {frame_idx} processed successfully")
        
        # Track target lipid binding state for each protein
        target_lipid_binding_state = {}
        lipid_composition_with_target_lipid = {}
        lipid_composition_without_target_lipid = {}
        
        # DEBUG: Print available lipid types and target lipid
        print(f"DEBUG: TARGET_LIPID = {TARGET_LIPID}")
        print(f"DEBUG: Available lipid types in lipid_sels = {list(lipid_sels.keys()) if lipid_sels else 'None'}")
        if proteins:
            first_protein = list(proteins.keys())[0]
            print(f"DEBUG: Available lipid types in lipid_contacts[{first_protein}] = {list(lipid_contacts[first_protein].keys()) if first_protein in lipid_contacts else 'None'}")
        
        for protein_name, protein in proteins.items():
            # Check if target lipid is in contact with this protein
            target_lipid_contact = False
            print(f"DEBUG: Checking {protein_name} for {TARGET_LIPID} contacts...")
            print(f"DEBUG: lipid_contacts[{protein_name}] keys = {list(lipid_contacts[protein_name].keys()) if protein_name in lipid_contacts else 'None'}")
            
            if TARGET_LIPID in lipid_contacts[protein_name]:
                target_lipid_contacts = lipid_contacts[protein_name][TARGET_LIPID]['contacts']
                contact_sum = np.sum(target_lipid_contacts)
                print(f"DEBUG: {protein_name} - {TARGET_LIPID} contact sum = {contact_sum}")
                if contact_sum > 0:
                    target_lipid_contact = True
            else:
                print(f"DEBUG: {TARGET_LIPID} not found in lipid_contacts[{protein_name}]")
            
            target_lipid_binding_state[protein_name] = target_lipid_contact
            
            # Track lipid composition based on target lipid binding state
            if target_lipid_contact:
                lipid_composition_with_target_lipid[protein_name] = {}
                for lipid_type in lipid_sels:
                    if lipid_type != TARGET_LIPID and lipid_type in lipid_contacts[protein_name]:
                        contact_count = np.sum(lipid_contacts[protein_name][lipid_type]['contacts'])
                        lipid_composition_with_target_lipid[protein_name][lipid_type] = contact_count
            else:
                lipid_composition_without_target_lipid[protein_name] = {}
                for lipid_type in lipid_sels:
                    if lipid_type != TARGET_LIPID and lipid_type in lipid_contacts[protein_name]:
                        contact_count = np.sum(lipid_contacts[protein_name][lipid_type]['contacts'])
                        lipid_composition_without_target_lipid[protein_name][lipid_type] = contact_count
        
        # Add residue-level information to results
        return {
            'frame': frame_idx,
            'protein_contacts': protein_contacts,
            'lipid_contacts': lipid_contacts,
            'residue_contacts': residue_contacts,
            'target_lipid_binding_state': target_lipid_binding_state,  # NEW: Target lipid binding state
            'lipid_composition_with_target_lipid': lipid_composition_with_target_lipid,  # NEW: Composition when target lipid bound
            'lipid_composition_without_target_lipid': lipid_composition_without_target_lipid,  # NEW: Composition when target lipid not bound
            'leaflet0': reference_leaflet0  # Keep leaflet information
        }
    except Exception as e:
        print(f"Error processing frame {frame_idx}: {str(e)}")
        traceback.print_exc()
        return None

def process_frames_serially(universe, proteins, frames, reference_leaflet0=None):
    """Process frames serially
    
    Parameters
    ----------
    universe : MDAnalysis.Universe
        The universe object
    proteins : dict
        Dictionary of protein selections
    frames : list
        List of frame indices to process
    reference_leaflet0 : MDAnalysis.AtomGroup, optional
        Reference leaflet to use
        
    Returns
    -------
    tuple
        (results, reference_leaflet0)
    """
    results = []
    
    for i, frame_idx in enumerate(frames):
        try:
            # Force LeafletFinder use for first frame
            force_update = (i == 0 and reference_leaflet0 is None)
            
            print(f"Processing frame {frame_idx}...")
            result = process_frame(universe, frame_idx, proteins, reference_leaflet0, force_update)
            
            if result:
                print(f"Frame {frame_idx} processed successfully")
                results.append(result)
                # Update reference leaflet from first frame
                if i == 0 and reference_leaflet0 is None and 'leaflet0' in result:
                    reference_leaflet0 = result['leaflet0']
                    print(f"Updated reference leaflet from first frame result")
            else:
                print(f"Frame {frame_idx} processing returned None")
        except Exception as e:
            print(f"Error processing frame {frame_idx}: {str(e)}")
            traceback.print_exc()
    
    return results, reference_leaflet0