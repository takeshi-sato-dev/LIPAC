"""Unit test: the exact counts of fast_stage1.py must equal the contact functions of
stage1_contact_analysis/core/contact_calculator.py on synthetic membranes, frame by frame.
Run from the repository root:  python stage1_fast/test_against_lipac.py
"""
import sys, os, types, importlib.util
import numpy as np, MDAnalysis as mda

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
pkg = types.ModuleType('stage1_contact_analysis'); pkg.__path__ = []
cfg = types.ModuleType('stage1_contact_analysis.config')
cfg.RESIDUE_OFFSET = 0; cfg.CONTACT_CUTOFF = 6.0; cfg.PROTEIN_CONTACT_CUTOFF = 6.0
cfg.TM_HELIX_RESID_RANGE = None; cfg.APPLY_RESIDUE_OFFSET = False; cfg.DIMER_CUTOFF = 20.0; cfg.TARGET_LIPID = 'DPG3'
core = types.ModuleType('stage1_contact_analysis.core'); core.__path__ = []
sys.modules.update({'stage1_contact_analysis': pkg, 'stage1_contact_analysis.config': cfg, 'stage1_contact_analysis.core': core})
spec = importlib.util.spec_from_file_location('stage1_contact_analysis.core.contact_calculator',
                                              os.path.join(ROOT, 'stage1_contact_analysis', 'core', 'contact_calculator.py'))
cc = importlib.util.module_from_spec(spec); spec.loader.exec_module(cc)
sys.path.insert(0, HERE)
import fast_stage1 as fs


def make_universe(seed, L=60.0, nres=12, nlip=160, nb=6):
    rng = np.random.default_rng(seed)
    types_ = ['CHOL', 'DPSM', 'DIPC', 'DPG3']
    n_at = 2 * nres * 3 + nlip * nb
    resindex = np.concatenate([np.repeat(np.arange(2 * nres), 3), np.repeat(np.arange(2 * nres, 2 * nres + nlip), nb)])
    u = mda.Universe.empty(n_at, n_residues=2 * nres + nlip, n_segments=3, atom_resindex=resindex,
                           residue_segindex=np.r_[np.zeros(nres, int), np.ones(nres, int), np.full(nlip, 2)], trajectory=True)
    u.add_TopologyAttr('segids', ['PROA', 'PROB', 'MEMB'])
    u.add_TopologyAttr('resnames', ['LEU'] * (2 * nres) + [types_[(i // 2) % 4] for i in range(nlip)])
    u.add_TopologyAttr('resids', list(range(54, 54 + nres)) * 2 + list(range(1, nlip + 1)))
    u.add_TopologyAttr('names', ['BB', 'SC1', 'SC2'] * (2 * nres) + (['PO4'] + [f'C{i}' for i in range(1, nb)]) * nlip)
    u.add_TopologyAttr('masses', rng.choice([45., 72.], n_at))
    pos = np.zeros((n_at, 3)); k = 0
    for x0, y0 in [(L - 1.0, 25.0), (2.0, 28.0)]:           # both copies straddle the boundary
        for r in range(nres):
            for b in range(3):
                pos[k] = [x0 + rng.normal(0, 1.5), y0 + rng.normal(0, 1.5), 22 + 3.5 * r + b]; k += 1
    for l in range(nlip):
        x0, y0 = rng.uniform(0, L, 2); up = l % 2 == 0
        for b in range(nb):
            pos[k] = [x0 + 1.8 * b * rng.normal(), y0 + 1.8 * b * rng.normal(), (55 - 3 * b) if up else (30 + 3 * b)]; k += 1
    pos[:, :2] %= L
    u.atoms.positions = pos; u.dimensions = np.r_[L, L, 90.0, 90, 90, 90]
    return u, [l + 1 for l in range(nlip) if l % 2 == 0], types_


ok = True
for seed in range(8):
    u, upper, types_ = make_universe(seed)
    args = types.SimpleNamespace(lipids=types_, leaflet_pickle=None, subset=None, start=0, out='/tmp', leaflet_z='all')
    fs.upper_leaflet_resindices = lambda u_, a_, f_: (np.array([r.resindex for r in u_.residues
                                                                if r.segid == 'MEMB' and int(r.resid) in set(upper)]), 'test')
    orig = fs.mda.Universe; fs.mda.Universe = lambda *a, **k: u
    args.psf = args.xtc = 'x'
    S = fs.attach_helpers(fs.build_state(args)); fs.mda.Universe = orig
    res_exact, _, mol_full, _, _ = fs.process_frame(S, 6.0)
    up_atoms = u.select_atoms('segid MEMB and resid ' + ' '.join(map(str, upper)))
    sels = {t: {'sel': [up_atoms.select_atoms(f'resname {t}')]} for t in types_}
    for c, s in enumerate(['PROA', 'PROB']):
        prot = u.select_atoms(f'segid {s}')
        r = cc.calculate_lipid_protein_contacts(prot, sels, u.dimensions[:3])
        m = cc.calculate_unique_lipid_protein_contacts(prot, sels, u.dimensions[:3])
        a = np.array([r[t]['contacts'] for t in types_]).T.astype(int)
        ok &= np.array_equal(res_exact[c], a) and np.array_equal(mol_full[c], [m[t] for t in types_])
    print('seed', seed, 'contacts', int(res_exact.sum()), 'OK' if ok else 'MISMATCH')
print('ALL PASS' if ok else 'FAIL')
