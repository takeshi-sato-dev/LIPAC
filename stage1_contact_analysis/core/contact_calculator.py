"""Contact calculation functions for protein-protein and lipid-protein interactions"""

import numpy as np
from ..config import RESIDUE_OFFSET, CONTACT_CUTOFF, PROTEIN_CONTACT_CUTOFF

def calculate_protein_protein_contacts(protein1, protein2, box, cutoff=PROTEIN_CONTACT_CUTOFF):
    """Calculate residue-residue contacts between two proteins with residue number conversion
    
    Parameters
    ----------
    protein1, protein2 : MDAnalysis.AtomGroup
        Protein selections
    box : array-like
        Box dimensions for PBC
    cutoff : float
        Contact cutoff distance in Angstrom
        
    Returns
    -------
    tuple
        Contact arrays and matrices with converted residue IDs
    """
    # Initialize contact matrices
    protein1_contacts = np.zeros(len(protein1.residues))
    protein2_contacts = np.zeros(len(protein2.residues))
    min_distances1 = np.ones(len(protein1.residues)) * float('inf')
    min_distances2 = np.ones(len(protein2.residues)) * float('inf')
    contact_matrix = np.zeros((len(protein1.residues), len(protein2.residues)))
    
    # Optimization: rough screening before distance calculation
    p1_com = protein1.center_of_mass()
    p2_com = protein2.center_of_mass()
    
    # COM distance calculation with PBC correction
    diff = p1_com - p2_com
    for dim in range(3):
        if diff[dim] > box[dim] * 0.5:
            diff[dim] -= box[dim]
        elif diff[dim] < -box[dim] * 0.5:
            diff[dim] += box[dim]
    
    com_dist = np.sqrt(np.sum(diff * diff))
    
    # Skip calculation if COMs are too far apart
    if com_dist > 30.0:  # If farther than 30Å, assume no contact
        print(f"Proteins too far apart: {com_dist:.2f}Å > 30Å, skipping contact calculation")
        
        # Store converted residue IDs
        residue_ids1 = [res.resid + RESIDUE_OFFSET for res in protein1.residues]
        residue_ids2 = [res.resid + RESIDUE_OFFSET for res in protein2.residues]
        
        return protein1_contacts, protein2_contacts, contact_matrix, min_distances1, min_distances2, residue_ids1, residue_ids2
    
    print(f"Calculating protein-protein contacts between {len(protein1.residues)} and {len(protein2.residues)} residues")
    
    # Calculate minimum distance between residues
    for i, res1 in enumerate(protein1.residues):
        for j, res2 in enumerate(protein2.residues):
            # Find minimum distance between atoms in residues
            min_dist = float('inf')
            
            # Optimization: use residue COMs for rough screening
            res1_com = res1.atoms.center_of_mass()
            res2_com = res2.atoms.center_of_mass()
            
            # Residue COM distance calculation with PBC
            res_diff = res1_com - res2_com
            for dim in range(3):
                if res_diff[dim] > box[dim] * 0.5:
                    res_diff[dim] -= box[dim]
                elif res_diff[dim] < -box[dim] * 0.5:
                    res_diff[dim] += box[dim]
            
            res_com_dist = np.sqrt(np.sum(res_diff * res_diff))
            
            # Skip if residue COMs are too far apart
            max_atom_dist = 10.0  # Estimated maximum distance between atoms in residue
            if res_com_dist > (cutoff + max_atom_dist):
                continue
            
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
                    
                    # Early termination if distance below cutoff found
                    if min_dist <= cutoff:
                        break
                
                if min_dist <= cutoff:
                    break
            
            # Save minimum distance
            min_distances1[i] = min(min_distances1[i], min_dist)
            min_distances2[j] = min(min_distances2[j], min_dist)
                
            # Record contact as binary (0/1) if within cutoff
            if min_dist <= cutoff:
                protein1_contacts[i] += 1.0
                protein2_contacts[j] += 1.0
                contact_matrix[i, j] = 1.0
    
    # Store converted residue IDs (add conversion offset to original resid)
    residue_ids1 = [res.resid + RESIDUE_OFFSET for res in protein1.residues]
    residue_ids2 = [res.resid + RESIDUE_OFFSET for res in protein2.residues]
    
    return protein1_contacts, protein2_contacts, contact_matrix, min_distances1, min_distances2, residue_ids1, residue_ids2

def calculate_lipid_protein_contacts(protein, lipid_sels, box, cutoff=CONTACT_CUTOFF):
    """Calculate residue-lipid contacts between protein and each lipid type (3D distance)
    
    Parameters
    ----------
    protein : MDAnalysis.AtomGroup
        Protein selection
    lipid_sels : dict
        Dictionary of lipid selections by type
    box : array-like
        Box dimensions for PBC
    cutoff : float
        Contact cutoff distance in Angstrom
        
    Returns
    -------
    dict
        Contact information for each lipid type
    """
    # Initialize results
    lipid_contacts = {}
    for lipid_type in lipid_sels:
        lipid_contacts[lipid_type] = np.zeros(len(protein.residues))
    
    # Loop through each residue
    for i, res in enumerate(protein.residues):
        # Save original residue ID (for debug output)
        original_resid = res.resid
        # Converted residue ID (for result reporting)
        converted_resid = original_resid + RESIDUE_OFFSET
        
        # Calculate average Z coordinate (height) of residue
        res_z_avg = np.mean([atom.position[2] for atom in res.atoms])
        
        # Loop through each lipid type
        for lipid_type, sel_info in lipid_sels.items():
            # Process only upper leaflet
            leaflet_sel = sel_info['sel'][0]
            if len(leaflet_sel) == 0:
                continue
            
            # Calculate average Z coordinate of leaflet (for optimization)
            leaflet_z_avg = np.mean([atom.position[2] for atom in leaflet_sel.atoms])
            
            # Skip if Z distance is too large (important optimization)
            z_diff = abs(res_z_avg - leaflet_z_avg)
            if z_diff > 15.0:  # Skip if separated by more than 15Å
                continue
            
            # Loop through each lipid residue
            lipid_count = 0
            for lipid_res in leaflet_sel.residues:
                # Calculate minimum distance between residue and lipid residue
                min_dist = float('inf')
                
                # Optimization: first check residue COM distance
                res_com = res.atoms.center_of_mass()
                lipid_com = lipid_res.atoms.center_of_mass()
                
                # COM distance calculation with PBC correction
                com_diff = res_com - lipid_com
                for dim in range(3):
                    if com_diff[dim] > box[dim] * 0.5:
                        com_diff[dim] -= box[dim]
                    elif com_diff[dim] < -box[dim] * 0.5:
                        com_diff[dim] += box[dim]
                
                com_dist = np.sqrt(np.sum(com_diff * com_diff))
                
                # Skip if COMs are too far apart (important optimization)
                max_atom_dist = 8.0  # Estimated maximum distance between atoms in residue
                if com_dist > (cutoff + max_atom_dist):
                    continue
                
                for atom in res.atoms:
                    for lipid_atom in lipid_res.atoms:
                        # 3D distance calculation with PBC correction
                        diff = atom.position - lipid_atom.position
                        for dim in range(3):  # Calculate for all dimensions (XYZ)
                            if diff[dim] > box[dim] * 0.5:
                                diff[dim] -= box[dim]
                            elif diff[dim] < -box[dim] * 0.5:
                                diff[dim] += box[dim]
                        
                        # 3D distance calculation
                        dist = np.sqrt(np.sum(diff * diff))
                        min_dist = min(min_dist, dist)
                        
                        # Early termination if distance below cutoff found
                        if min_dist <= cutoff:
                            break
                    
                    if min_dist <= cutoff:
                        break
                
                # Count as contact if within cutoff (binary method)
                if min_dist <= cutoff:
                    lipid_contacts[lipid_type][i] += 1
                    lipid_count += 1
    
    # Return results including residue IDs and contact information
    result = {}
    for lipid_type in lipid_sels:
        # Create list of converted residue IDs
        residue_ids = [res.resid + RESIDUE_OFFSET for res in protein.residues]
        result[lipid_type] = {
            'contacts': lipid_contacts[lipid_type],
            'residue_ids': residue_ids
        }
    
    return result

def calculate_unique_lipid_protein_contacts(protein, lipid_sels, box, cutoff=CONTACT_CUTOFF):
    """Calculate unique lipid molecule counts in contact with protein (no double counting)
    
    Parameters
    ----------
    protein : MDAnalysis.AtomGroup
        Protein selection
    lipid_sels : dict
        Dictionary of lipid selections by type
    box : array-like
        Box dimensions for PBC
    cutoff : float
        Contact cutoff distance in Angstrom
        
    Returns
    -------
    dict
        Unique lipid molecule counts for each lipid type
    """
    unique_contacts = {}
    
    # Pre-calculate protein COM and positions for optimization
    protein_com = protein.center_of_mass()
    protein_positions = np.array([atom.position for atom in protein.atoms])
    
    # Loop through each lipid type
    for lipid_type, sel_info in lipid_sels.items():
        # Process only upper leaflet
        leaflet_sel = sel_info['sel'][0]
        if len(leaflet_sel) == 0:
            unique_contacts[lipid_type] = 0
            continue
        
        # Set to store unique lipid molecules in contact
        contacted_lipids = set()
        
        # Group lipids by molecule (residue)
        lipid_residues = leaflet_sel.residues
        
        # Remove overly aggressive Z-coordinate filtering that was causing issues
        # Will rely on COM distance filtering instead
        
        # Check each lipid molecule with optimization
        for lipid_res in lipid_residues:
            # Optimization: first check XY-plane distance only
            lipid_com = lipid_res.atoms.center_of_mass()
            
            # XY-plane distance calculation with PBC correction
            xy_diff = lipid_com[:2] - protein_com[:2]  # Only X and Y components
            xy_diff = xy_diff - box[:2] * np.round(xy_diff / box[:2])  # PBC in XY
            xy_dist = np.sqrt(np.sum(xy_diff * xy_diff))
            
            # Skip if XY distance is too large (10Å cutoff for XY plane)
            if xy_dist > 10.0:
                continue
            
            # Fully vectorized distance calculation for all atom pairs
            lipid_positions = np.array([atom.position for atom in lipid_res.atoms])
            
            # Calculate all pairwise distances at once
            # Shape: (n_lipid_atoms, n_protein_atoms, 3)
            diff = lipid_positions[:, np.newaxis, :] - protein_positions[np.newaxis, :, :]
            diff = diff - box * np.round(diff / box)
            
            # Calculate distances: shape (n_lipid_atoms, n_protein_atoms)
            distances = np.sqrt(np.sum(diff * diff, axis=2))
            
            # Check if any distance is within cutoff
            if np.any(distances <= cutoff):
                contacted_lipids.add(lipid_res.resid)
        
        unique_contacts[lipid_type] = len(contacted_lipids)
    
    return unique_contacts

def check_tm_helix_interactions(protein1, protein2, box, protein_cutoff=6.0):
    """Check for interactions between TM helices
    More stringent dimer detection by checking if TM domains actually interact
    
    Parameters
    ----------
    protein1, protein2 : MDAnalysis.AtomGroup
        Protein selections
    box : array-like
        Box dimensions for PBC
    protein_cutoff : float
        Cutoff distance for protein contacts
        
    Returns
    -------
    bool
        True if TM domains interact
    """
    # Select central part of TM helix (middle residues)
    # Residue range is 65:103, so define 75-95 as central part
    protein1_core = protein1.select_atoms("resid 68:74")
    protein2_core = protein2.select_atoms("resid 68:74")
    
    if len(protein1_core) == 0 or len(protein2_core) == 0:
        return False
    
    # Calculate minimum distance between TM helix centers
    min_distance = float('inf')
    
    for atom1 in protein1_core.atoms:
        for atom2 in protein2_core.atoms:
            # Distance calculation with PBC correction
            diff = atom1.position - atom2.position
            for dim in range(3):
                if diff[dim] > box[dim] * 0.5:
                    diff[dim] -= box[dim]
                elif diff[dim] < -box[dim] * 0.5:
                    diff[dim] += box[dim]
            
            dist = np.sqrt(np.sum(diff * diff))
            min_distance = min(min_distance, dist)
            
            # Early termination: if any close atom pair found, consider as interaction
            if min_distance <= protein_cutoff:
                return True
    
    # Judge based on minimum distance
    return min_distance <= protein_cutoff

def calculate_protein_com_distances(universe, proteins):
    """Calculate COM distances between protein pairs and identify close pairs with TM domain interactions"""
    from ..config import DIMER_CUTOFF, TARGET_LIPID
    
    print("Calculating protein-protein COM distances and TM helix interactions...")
    box = universe.dimensions[:3]
    
    # Calculate COM for each protein
    protein_coms = {}
    for protein_name, protein in proteins.items():
        if len(protein) > 0:
            protein_coms[protein_name] = protein.center_of_mass()
    
    # Identify close protein pairs
    close_pairs = {}
    
    protein_names = list(protein_coms.keys())
    for i in range(len(protein_names)):
        for j in range(i + 1, len(protein_names)):
            protein1_name = protein_names[i]
            protein2_name = protein_names[j]
            
            if protein1_name not in protein_coms or protein2_name not in protein_coms:
                continue
            
            # Calculate COM distance (considering periodic boundary conditions)
            com1 = protein_coms[protein1_name]
            com2 = protein_coms[protein2_name]
            
            diff = com1 - com2
            # PBC correction
            for dim in range(3):
                if diff[dim] > box[dim] * 0.5:
                    diff[dim] -= box[dim]
                elif diff[dim] < -box[dim] * 0.5:
                    diff[dim] += box[dim]
            
            dist = np.sqrt(np.sum(diff * diff))
            
            # Apply basic distance cutoff
            if dist <= DIMER_CUTOFF:
                # Check for target lipid presence
                has_target_lipid = len(universe.select_atoms(f"resname {TARGET_LIPID}")) > 0
                
                # Check TM domain interactions (for logging purposes only)
                if has_target_lipid:
                    has_tm_interactions = check_tm_helix_interactions(
                        proteins[protein1_name], 
                        proteins[protein2_name], 
                        box,
                        protein_cutoff=6.0
                    )
                    
                    if not has_tm_interactions:
                        print(f"  {protein1_name}-{protein2_name}: distance {dist:.2f} Å (no TM interactions)")
                    else:
                        print(f"  {protein1_name}-{protein2_name}: distance {dist:.2f} Å (with TM interactions)")
                
                # Record the pair
                pair_name = f"{protein1_name}-{protein2_name}"
                close_pairs[pair_name] = dist
                print(f"  Found close pair: {pair_name}, distance: {dist:.2f} Å")
    
    
    print(f"Found {len(close_pairs)} close protein pairs")
    return close_pairs  