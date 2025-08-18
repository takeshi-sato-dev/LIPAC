"""Trajectory loading and leaflet identification functions"""

import os
import pickle
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.leaflet import LeafletFinder
from ..config import LEAFLET_CUTOFF, LIPID_TYPES, TEMP_FILES_DIR

def load_universe(top_file, traj_file):
    """Load trajectory
    
    Parameters
    ----------
    top_file : str
        Path to topology file
    traj_file : str
        Path to trajectory file
        
    Returns
    -------
    MDAnalysis.Universe or None
        Loaded universe object or None if failed
    """
    print(f"Loading topology: {top_file}")
    print(f"Loading trajectory: {traj_file}")
    
    # Check file existence
    if not os.path.exists(top_file):
        print(f"ERROR: Topology file not found: {top_file}")
        return None
        
    if not os.path.exists(traj_file):
        print(f"ERROR: Trajectory file not found: {traj_file}")
        return None
        
    try:
        universe = mda.Universe(top_file, traj_file)
        print(f"Trajectory loaded successfully with {len(universe.trajectory)} frames.")
        return universe
    except Exception as e:
        print(f"Error loading universe: {str(e)}")
        import traceback
        traceback.print_exc()
        return None

def identify_lipid_leaflets(universe, leaflet_cutoff=LEAFLET_CUTOFF):
    """Identify and return lipid leaflets - one-time detection followed by reuse
    
    Parameters
    ----------
    universe : MDAnalysis.Universe
        The universe object
    leaflet_cutoff : float
        Cutoff distance for leaflet detection
        
    Returns
    -------
    tuple
        (upper_leaflet, lower_leaflet) AtomGroups
    """
    print("Identifying lipid leaflets...")
    
    # Path for saving leaflet information
    leaflet_info_file = os.path.join(TEMP_FILES_DIR, "leaflet_info.pickle")
    
    # Load existing leaflet information if available
    if os.path.exists(leaflet_info_file):
        print(f"Loading existing leaflet information from {leaflet_info_file}")
        try:
            with open(leaflet_info_file, 'rb') as f:
                leaflet_info = pickle.load(f)
            
            # Reconstruct leaflets from saved resids
            upper_resids = leaflet_info['upper_leaflet_resids']
            lower_resids = leaflet_info.get('lower_leaflet_resids', [])
            
            # Batch selection to avoid recursion errors
            leaflet0 = mda.AtomGroup([], universe)
            leaflet1 = mda.AtomGroup([], universe)
            
            # Select only lipid atoms first
            lipid_resnames = " ".join(LIPID_TYPES)
            lipid_atoms = universe.select_atoms(f"resname {lipid_resnames}")
            
            batch_size = 100
            # Upper leaflet reconstruction
            print(f"Reconstructing upper leaflet with {len(upper_resids)} residue IDs...")
            for i in range(0, len(upper_resids), batch_size):
                batch = upper_resids[i:i+batch_size]
                batch_sel_str = " or ".join([f"resid {r}" for r in batch])
                try:
                    batch_atoms = lipid_atoms.select_atoms(batch_sel_str)
                    leaflet0 = leaflet0.union(batch_atoms)
                except Exception as e:
                    print(f"Error selecting upper leaflet batch: {str(e)}")
            
            # Lower leaflet reconstruction
            if lower_resids:
                print(f"Reconstructing lower leaflet with {len(lower_resids)} residue IDs...")
                for i in range(0, len(lower_resids), batch_size):
                    batch = lower_resids[i:i+batch_size]
                    batch_sel_str = " or ".join([f"resid {r}" for r in batch])
                    try:
                        batch_atoms = lipid_atoms.select_atoms(batch_sel_str)
                        leaflet1 = leaflet1.union(batch_atoms)
                    except Exception as e:
                        print(f"Error selecting lower leaflet batch: {str(e)}")
            else:
                leaflet1 = leaflet0
            
            print(f"Leaflet reconstruction completed")
            print(f"  Upper leaflet: {len(leaflet0.residues)} molecules")
            print(f"  Lower leaflet: {len(leaflet1.residues)} molecules")
            
            return leaflet0, leaflet1
            
        except Exception as e:
            print(f"Error loading existing leaflet information: {str(e)}")
            print("Proceeding with new leaflet detection...")
    
    # Perform new leaflet detection if no existing file
    print("Running LeafletFinder to detect lipid leaflets...")
    
    # Target all lipid head groups including CHOL
    L = LeafletFinder(universe, "name PO4 ROH GL1 GL2 AM1 AM2 GM1 GM2")
    
    # Use appropriate cutoff value
    cutoff = L.update(leaflet_cutoff)
    print(f"Cutoff distance: {cutoff}")
    print(f"Number of leaflets found: {len(L.components)}")
    
    if len(L.components) < 2:
        print("Warning: Less than 2 leaflets found. Using the same leaflet for both.")
        leaflet0 = L.groups(0)
        leaflet1 = leaflet0
        upper_resids = [res.resid for res in leaflet0.residues]
        lower_resids = []
    else:
        temp_leaflet0 = L.groups(0)
        temp_leaflet1 = L.groups(1)
        
        # Count DPSM content to determine upper leaflet
        dpsm_count0 = len(temp_leaflet0.select_atoms("resname DPSM").residues)
        dpsm_count1 = len(temp_leaflet1.select_atoms("resname DPSM").residues)
        
        print(f"DPSM content: Leaflet 0 = {dpsm_count0}, Leaflet 1 = {dpsm_count1}")
        
        # Use the one with more DPSM as upper leaflet
        if dpsm_count0 >= dpsm_count1:
            leaflet0 = temp_leaflet0
            leaflet1 = temp_leaflet1
        else:
            leaflet0 = temp_leaflet1
            leaflet1 = temp_leaflet0
        
        upper_resids = [res.resid for res in leaflet0.residues]
        lower_resids = [res.resid for res in leaflet1.residues]
    
    # Save leaflet information
    os.makedirs(TEMP_FILES_DIR, exist_ok=True)
    
    leaflet_info = {
        'upper_leaflet_resids': upper_resids,
        'lower_leaflet_resids': lower_resids,
    }
    
    try:
        with open(leaflet_info_file, 'wb') as f:
            pickle.dump(leaflet_info, f)
        print(f"Leaflet information saved to {leaflet_info_file}")
    except Exception as e:
        print(f"Warning: Could not save leaflet information: {str(e)}")
    
    print("Lipid leaflets identified successfully.")
    return leaflet0, leaflet1

def select_lipids(universe, leaflet0):
    """Select lipids - using only leaflet0 (upper leaflet)
    Note: Lipid-protein contact analysis is performed only for leaflet0
    
    Parameters
    ----------
    universe : MDAnalysis.Universe
        The universe object
    leaflet0 : MDAnalysis.AtomGroup
        Upper leaflet atoms
        
    Returns
    -------
    dict
        Dictionary of lipid selections by type
    """
    from ..config import LIPID_TYPES
    
    # Initialize lipid selections
    lipid_sels = {}
    
    # Statistics for upper leaflet
    total_molecules = len(leaflet0.residues)
    total_atoms = len(leaflet0.atoms)
    print(f"\nUpper leaflet statistics: {total_molecules} molecules ({total_atoms} atoms)")
    
    # Select each lipid type and count molecules accurately
    for lipid_type in LIPID_TYPES:
        try:
            lipid_sel = leaflet0.select_atoms(f"resname {lipid_type}")
            lipid_molecule_count = len(lipid_sel.residues)
            lipid_atom_count = len(lipid_sel.atoms)
            
            # Save selection
            lipid_sels[lipid_type] = {'sel': [lipid_sel]}
            
            # Display accurate molecule and atom counts
            percentage = (lipid_molecule_count / total_molecules) * 100 if total_molecules > 0 else 0
            print(f"{lipid_type}: {lipid_molecule_count} molecules ({percentage:.1f}%) - {lipid_atom_count} atoms")
            
        except Exception as e:
            print(f"Error selecting {lipid_type}: {str(e)}")
            lipid_sels[lipid_type] = {'sel': [mda.AtomGroup([], universe)]}
            print(f"{lipid_type}: 0 molecules (0.0%)")
    
    return lipid_sels

def select_proteins(universe, n_proteins=4):
    """Select proteins - set residue numbers correctly
    
    Parameters
    ----------
    universe : MDAnalysis.Universe
        The universe object
    n_proteins : int
        Number of proteins to select
        
    Returns
    -------
    dict
        Dictionary of protein selections
    """
    from ..config import RESIDUE_OFFSET
    
    proteins = {}
    segids = ['PROA', 'PROB', 'PROC', 'PROD'][:n_proteins]
    
    for i, segid in enumerate(segids, 1):
        protein_name = f"Protein_{i}"
        # Select transmembrane helix region (residue range can be adjusted)
        # Original residue range 65:103 converted to: 621:659
        try:
            orig_selection = universe.select_atoms(f"segid {segid} and resid 65:103")
            proteins[protein_name] = orig_selection
            
            # Display converted residue ID range (for confirmation)
            converted_range = f"{65+RESIDUE_OFFSET}:{103+RESIDUE_OFFSET}"
            print(f"{protein_name}: {len(proteins[protein_name])} atoms, "
                  f"{len(proteins[protein_name].residues)} residues")
            print(f"  Original residue range: 65:103")
            print(f"  Converted residue range: {converted_range}")
        except Exception as e:
            print(f"Warning: Could not select {protein_name} with segid {segid}: {str(e)}")
            proteins[protein_name] = mda.AtomGroup([], universe)
    
    return proteins