"""Contact calculation functions for protein-protein and lipid-protein interactions"""

import numpy as np
from ..config import RESIDUE_OFFSET, CONTACT_CUTOFF, PROTEIN_CONTACT_CUTOFF, TM_HELIX_RESID_RANGE, APPLY_RESIDUE_OFFSET

from MDAnalysis.lib.distances import distance_array as _distance_array


def _full_box(box):
    """[lx, ly, lz, 90, 90, 90] from a box given as [lx, ly, lz] or as the full 6 values.
    The earlier code applied the minimum image per axis, i.e. it assumed a rectangular box;
    the same assumption is kept here."""
    b = np.asarray(box, dtype=np.float32)
    return np.r_[b[:3], 90.0, 90.0, 90.0].astype(np.float32)


def _whole_com(atomgroup, box):
    """Center of mass of an AtomGroup after making it whole across the periodic box.
    Every atom is placed at its minimum image relative to the first atom, which is
    exact whenever the group is smaller than half the box (true for a residue, a lipid
    or a peptide)."""
    pos = atomgroup.positions.astype(float)
    b = np.asarray(box[:3], dtype=float)
    d = pos - pos[0]
    pos = pos[0] + d - b * np.round(d / b)
    m = atomgroup.masses
    return (pos * m[:, None]).sum(0) / m.sum()


def _residue_index(atomgroup):
    """Index (0..n_res-1) of the residue of every atom, in the order of atomgroup.residues."""
    order = {ri: k for k, ri in enumerate(atomgroup.residues.resindices)}
    return np.array([order[ri] for ri in atomgroup.resindices])


def _residue_min_distance(ag1, ag2, box):
    """Exact minimum bead-bead distance (PBC) for every residue of ag1 x every residue of ag2."""
    D = _distance_array(ag1.positions, ag2.positions, box=_full_box(box))
    r1 = _residue_index(ag1); r2 = _residue_index(ag2)
    n1 = len(ag1.residues); n2 = len(ag2.residues)
    M = np.full((n1, D.shape[1]), np.inf)
    np.minimum.at(M, r1, D)
    R = np.full((n1, n2), np.inf)
    np.minimum.at(R.T, r2, M.T)
    return R

def calculate_protein_protein_contacts(protein1, protein2, box, cutoff=PROTEIN_CONTACT_CUTOFF):
    """Residue-residue contacts between two proteins.

    A residue pair is in contact when any bead of one residue lies within `cutoff` of any
    bead of the other (minimum image). All residue pairs are evaluated; the COM
    prescreens of the earlier version, which dropped real contacts when a protein was
    split across the periodic boundary, are removed.
    Returns (protein1_contacts, protein2_contacts, contact_matrix, min_distances1,
    min_distances2, residue_ids1, residue_ids2) as before; min_distances are exact.
    """
    R = _residue_min_distance(protein1, protein2, box)
    contact_matrix = (R <= cutoff).astype(float)
    protein1_contacts = contact_matrix.sum(1)
    protein2_contacts = contact_matrix.sum(0)
    min_distances1 = R.min(1)
    min_distances2 = R.min(0)
    if APPLY_RESIDUE_OFFSET:
        residue_ids1 = [res.resid + RESIDUE_OFFSET for res in protein1.residues]
        residue_ids2 = [res.resid + RESIDUE_OFFSET for res in protein2.residues]
    else:
        residue_ids1 = [res.resid for res in protein1.residues]
        residue_ids2 = [res.resid for res in protein2.residues]
    return protein1_contacts, protein2_contacts, contact_matrix, min_distances1, min_distances2, residue_ids1, residue_ids2

def calculate_lipid_protein_contacts(protein, lipid_sels, box, cutoff=CONTACT_CUTOFF):
    """Residue-lipid contacts between a protein and each lipid type (3D, minimum image).

    For every protein residue and lipid type, the count is the number of lipid molecules
    of that type in the upper leaflet with any bead within `cutoff` of any bead of the
    residue. Every bead of every lipid molecule is evaluated. The z prescreen and the
    COM prescreen of the earlier version are removed: both dropped lipids that were in
    contact (molecules split across the periodic boundary, and lipid tails next to
    residues in the bilayer core).
    """
    lipid_contacts = {}
    n_res = len(protein.residues)
    for lipid_type, sel_info in lipid_sels.items():
        leaflet_sel = sel_info['sel'][0]
        counts = np.zeros(n_res)
        if len(leaflet_sel) > 0:
            lip_atoms = leaflet_sel.residues.atoms
            R = _residue_min_distance(protein, lip_atoms, box)
            counts = (R <= cutoff).sum(1).astype(float)
        lipid_contacts[lipid_type] = counts
    if APPLY_RESIDUE_OFFSET:
        residue_ids = [res.resid + RESIDUE_OFFSET for res in protein.residues]
    else:
        residue_ids = [res.resid for res in protein.residues]
    return {lt: {'contacts': lipid_contacts[lt], 'residue_ids': residue_ids} for lt in lipid_sels}

def calculate_unique_lipid_protein_contacts(protein, lipid_sels, box, cutoff=CONTACT_CUTOFF):
    """Number of distinct lipid molecules of each type with any bead within `cutoff` of any
    bead of the protein selection (minimum image). The 10 A in-plane COM prescreen of the
    earlier version, which dropped molecules in contact, is removed."""
    unique_contacts = {}
    for lipid_type, sel_info in lipid_sels.items():
        leaflet_sel = sel_info['sel'][0]
        if len(leaflet_sel) == 0:
            unique_contacts[lipid_type] = 0
            continue
        lip_atoms = leaflet_sel.residues.atoms
        D = _distance_array(protein.positions, lip_atoms.positions, box=_full_box(box))
        hit = (D <= cutoff).any(0)
        unique_contacts[lipid_type] = int(len(np.unique(lip_atoms.resindices[hit])))
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
    # Calculate central part of TM helix from config
    if TM_HELIX_RESID_RANGE:
        start, end = map(int, TM_HELIX_RESID_RANGE.split(':'))
        total_residues = end - start + 1
        # Use central ~1/4 of TM helix (approximately 7 residues for a 26-residue helix)
        center = (start + end) // 2
        half_width = max(3, total_residues // 8)  # At least 3 residues on each side
        core_start = center - half_width
        core_end = center + half_width
        core_range = f"{core_start}:{core_end}"
    else:
        # Fallback to original hardcoded range
        core_range = "68:74"

    protein1_core = protein1.select_atoms(f"resid {core_range}")
    protein2_core = protein2.select_atoms(f"resid {core_range}")
    
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
    # Calculate TM region COM for each protein (for screening only)
    protein_coms = {}
    for protein_name, protein in proteins.items():
        if len(protein) > 0:
            # Use TM helix region COM for dimer screening
            if TM_HELIX_RESID_RANGE:
                tm_region = protein.select_atoms(f"resid {TM_HELIX_RESID_RANGE}")
                if len(tm_region) > 0:
                    protein_coms[protein_name] = _whole_com(tm_region, box)
                else:
                    # Fallback to full protein if TM region not found
                    protein_coms[protein_name] = _whole_com(protein, box)
            else:
                protein_coms[protein_name] = _whole_com(protein, box)
    
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