"""Unit test for the exact contact functions of stage1_contact_analysis/core/contact_calculator.py.

Synthetic membranes are built with molecules deliberately split across the periodic
boundary and with lipid tails reaching residues in the bilayer core. Every count
returned by the LIPAC functions is compared with a brute-force count (any bead of the
lipid molecule within the cutoff of any bead of the residue, minimum image).
Run from the repository root:  python validation/test_exact_contacts.py
"""
import sys, os, types, importlib.util
import numpy as np, MDAnalysis as mda

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
pkg = types.ModuleType('stage1_contact_analysis'); pkg.__path__ = []
cfg = types.ModuleType('stage1_contact_analysis.config')
cfg.RESIDUE_OFFSET = 0; cfg.CONTACT_CUTOFF = 6.0; cfg.PROTEIN_CONTACT_CUTOFF = 6.0
cfg.TM_HELIX_RESID_RANGE = None; cfg.APPLY_RESIDUE_OFFSET = False; cfg.DIMER_CUTOFF = 20.0; cfg.TARGET_LIPID = 'DPG3'
core = types.ModuleType('stage1_contact_analysis.core'); core.__path__ = []
sys.modules.update({'stage1_contact_analysis': pkg, 'stage1_contact_analysis.config': cfg, 'stage1_contact_analysis.core': core})
spec = importlib.util.spec_from_file_location('stage1_contact_analysis.core.contact_calculator',
                                              os.environ.get('CC_FILE', os.path.join(ROOT, 'stage1_contact_analysis', 'core', 'contact_calculator.py')))
cc = importlib.util.module_from_spec(spec); spec.loader.exec_module(cc)


def build(seed, L=50.0, nres=10, nlip=120, nb=6):
    rng = np.random.default_rng(seed)
    n_prot = 2 * nres * 3; n_at = n_prot + nlip * nb
    resindex = np.concatenate([np.repeat(np.arange(2 * nres), 3), np.repeat(np.arange(2 * nres, 2 * nres + nlip), nb)])
    u = mda.Universe.empty(n_at, n_residues=2 * nres + nlip, n_segments=3, atom_resindex=resindex,
                           residue_segindex=np.r_[np.zeros(nres, int), np.ones(nres, int), np.full(nlip, 2)], trajectory=True)
    u.add_TopologyAttr('segids', ['PROA', 'PROB', 'MEMB'])
    types_ = ['CHOL', 'DPSM', 'DIPC', 'DPG3']
    u.add_TopologyAttr('resnames', ['LEU'] * (2 * nres) + [types_[i % 4] for i in range(nlip)])
    u.add_TopologyAttr('resids', list(range(1, nres + 1)) * 2 + list(range(1, nlip + 1)))
    u.add_TopologyAttr('names', ['BB', 'SC1', 'SC2'] * (2 * nres) + ['PO4'] + ['C'] * (nb - 1) * 0 + (['PO4'] + [f'C{i}' for i in range(1, nb)]) * nlip if False else ['BB', 'SC1', 'SC2'] * (2 * nres) + (['PO4'] + [f'C{i}' for i in range(1, nb)]) * nlip)
    u.add_TopologyAttr('masses', rng.choice([45., 72.], n_at))
    pos = np.zeros((n_at, 3)); k = 0
    # two proteins placed near the box edge so that they straddle the boundary
    for p, (x0, y0) in enumerate([(L - 1.0, 20.0), (3.0, 22.0)]):
        for r in range(nres):
            for b in range(3):
                pos[k] = [x0 + rng.normal(0, 1.5), y0 + rng.normal(0, 1.5), 20 + 3.5 * r + b]; k += 1
    for l in range(nlip):
        x0, y0 = rng.uniform(0, L, 2); up = l % 2 == 0
        for b in range(nb):                       # tails run toward the core and tilt laterally
            pos[k] = [x0 + 1.8 * b * rng.normal(), y0 + 1.8 * b * rng.normal(), (50 - 3 * b) if up else (26 + 3 * b)]; k += 1
    pos[:, :2] %= L                                # wrap: many molecules and both proteins become split
    u.atoms.positions = pos; u.dimensions = np.r_[L, L, 80.0, 90, 90, 90]
    return u, types_


def brute(prot, lip_atoms, box, cutoff=6.0):
    counts = np.zeros(len(prot.residues))
    for i, res in enumerate(prot.residues):
        hit = set()
        for lr in lip_atoms.residues:
            d = res.atoms.positions[:, None, :] - lr.atoms.positions[None, :, :]
            d -= box[:3] * np.round(d / box[:3])
            if (np.sqrt((d ** 2).sum(-1)) <= cutoff).any():
                hit.add(lr.resindex)
        counts[i] = len(hit)
    return counts


ok = True
for seed in range(6):
    u, types_ = build(seed)
    box = u.dimensions
    upper = u.select_atoms('segid MEMB and resid ' + ' '.join(str(r) for r in range(1, 121, 2)))
    sels = {t: {'sel': [upper.select_atoms(f'resname {t}')]} for t in types_}
    for s in ['PROA', 'PROB']:
        prot = u.select_atoms(f'segid {s}')
        res = cc.calculate_lipid_protein_contacts(prot, sels, box[:3])
        uni = cc.calculate_unique_lipid_protein_contacts(prot, sels, box[:3])
        for t in types_:
            b = brute(prot, sels[t]['sel'][0].residues.atoms, box)
            same = np.array_equal(res[t]['contacts'], b)
            lipall = sels[t]['sel'][0].residues.atoms
            d = prot.positions[:, None, :] - lipall.positions[None, :, :]; d -= box[:3] * np.round(d / box[:3])
            ub = len(set(lipall.resindices[(np.sqrt((d ** 2).sum(-1)) <= 6.0).any(0)]))
            ok &= same and uni[t] == ub
    p1, p2, M, *_ = cc.calculate_protein_protein_contacts(u.select_atoms('segid PROA'), u.select_atoms('segid PROB'), box[:3])
    A, B = u.select_atoms('segid PROA'), u.select_atoms('segid PROB')
    Mb = np.zeros_like(M)
    for i, r1 in enumerate(A.residues):
        for j, r2 in enumerate(B.residues):
            d = r1.atoms.positions[:, None, :] - r2.atoms.positions[None, :, :]; d -= box[:3] * np.round(d / box[:3])
            Mb[i, j] = float((np.sqrt((d ** 2).sum(-1)) <= 6.0).any())
    ok &= np.array_equal(M, Mb)
    print('seed', seed, 'protein-protein contacts', int(M.sum()), 'match' if np.array_equal(M, Mb) else 'MISMATCH')
print('ALL PASS' if ok else 'FAIL')
