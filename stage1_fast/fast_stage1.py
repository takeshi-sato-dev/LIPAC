"""
LIPAC Stage 1, fast and parallel re-implementation (per-frame contacts at full time resolution).

What it computes for every frame and every protein copy
--------------------------------------------------------
Two contact definitions are computed side by side.

1. exact  : a residue-lipid contact is counted when any bead of the lipid molecule lies within
            CUTOFF (6 A) of any bead of the residue (minimum image). No heuristic prefilters.
            This is the definition of core/contact_calculator.py from LIPAC 3 on.
2. legacy : the definition of calculate_lipid_protein_contacts in LIPAC 1 and 2
            (every commit up to a010a80), reproduced for comparison with earlier results, including its
            two prefilters, which discard real contacts:
              (a) the residue is skipped for a lipid type when |mean z of the residue beads -
                  mean z of all beads of that lipid type in the upper leaflet| > 15 A;
              (b) a lipid is skipped when the distance between the residue COM and the lipid
                  COM (mass-weighted, raw coordinates, per-dimension minimum image) exceeds
                  CUTOFF + 8 A.

Per frame the script stores
  res_exact[copy, residue, lipid_type]   residue-lipid contact counts (exact definition)
  res_legacy[copy, residue, lipid_type]  residue-lipid contact counts (legacy definition)
  mol_full[copy, lipid_type]             number of distinct lipid molecules in contact with any
                                         residue of the copy (exact definition, full protein)
  mol_subset[copy, lipid_type]           same, restricted to the residue subset (if given)
  mol_legacy[copy, lipid_type]           the legacy "unique molecules" count, which keeps only
                                         lipids whose COM lies within 10 A in the xy plane of
                                         the COM of the selected protein (reproduced for checks)
Residue-level arrays cover the full protein, so any residue subset can be summed afterwards.

Leaflet membership is fixed, as in the published Stage 1: it is read from an existing
leaflet_info.pickle (upper_leaflet_resids) or detected once with LeafletFinder at the first
analysed frame, using the published headgroup selection, a 10 A cutoff, and the rule that the
leaflet with more DPSM is the upper leaflet.

Parallelisation: the frame range is split into chunks; every worker opens the trajectory once
and processes whole chunks. Each chunk is written to its own .npz, so an interrupted run
resumes where it stopped. merge_and_export.py assembles the chunks and writes the CSV files.

Usage
  python fast_stage1.py --psf step5_assembly.psf --xtc step7_production.xtc \
      --start 20000 --stop 80000 --step 1 --target DPG3 \
      --lipids CHOL DPSM DIPC DPG3 DOPS --subset 65:103 \
      --leaflet-pickle .../temp_files/leaflet_info.pickle \
      --out out_EGFR1_GM3 --workers 8 --chunk 500
"""
import os
import sys
import json
import time
import pickle
import argparse
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import capped_distance

HEADGROUP_SEL = "name PO4 ROH GL1 GL2 AM1 AM2 GM1 GM2"
SEGIDS = ['PROA', 'PROB', 'PROC', 'PROD', 'PROE', 'PROF', 'PROG', 'PROH']

_G = {}   # per-worker state


# ----------------------------------------------------------------------------- setup
def upper_leaflet_resindices(u, args, frame):
    """Resindices of upper-leaflet lipid molecules (fixed for the run)."""
    lipid_sel = "resname " + " ".join(args.lipids)
    lipids = u.select_atoms(lipid_sel)
    if args.leaflet_pickle and os.path.exists(args.leaflet_pickle):
        info = pickle.load(open(args.leaflet_pickle, 'rb'))
        up = set(int(r) for r in info['upper_leaflet_resids'])
        res = [r for r in lipids.residues if int(r.resid) in up]
        return np.array([r.resindex for r in res]), 'pickle'
    from MDAnalysis.analysis.leaflet import LeafletFinder
    u.trajectory[frame]
    L = LeafletFinder(u, HEADGROUP_SEL)
    L.update(10.0)
    if len(L.components) < 2:
        sys.exit('LeafletFinder found fewer than two leaflets')
    g0, g1 = L.groups(0), L.groups(1)
    n0 = len(g0.select_atoms('resname DPSM').residues); n1 = len(g1.select_atoms('resname DPSM').residues)
    upper = g0 if n0 >= n1 else g1
    up = set(int(r) for r in upper.residues.resids)
    res = [r for r in lipids.residues if int(r.resid) in up]
    os.makedirs(args.out, exist_ok=True)
    pickle.dump({'upper_leaflet_resids': sorted(up),
                 'lower_leaflet_resids': sorted(int(r) for r in (g1 if upper is g0 else g0).residues.resids)},
                open(os.path.join(args.out, 'leaflet_info.pickle'), 'wb'))
    return np.array([r.resindex for r in res]), 'LeafletFinder'


def build_state(args):
    u = mda.Universe(args.psf, args.xtc)
    segs = [s for s in SEGIDS if len(u.select_atoms(f'segid {s}')) > 0]
    if not segs:
        sys.exit('no protein segments PROA.. found')
    prot = u.select_atoms('segid ' + ' '.join(segs))
    copies = [u.select_atoms(f'segid {s}') for s in segs]
    nres = [len(c.residues) for c in copies]
    if len(set(nres)) != 1:
        sys.exit(f'protein copies differ in length: {nres}')
    nres = nres[0]
    resids = copies[0].residues.resids.copy()
    # per protein atom: copy index and residue index within copy
    a_copy = np.empty(len(prot), int); a_res = np.empty(len(prot), int)
    pos_in_prot = {ix: k for k, ix in enumerate(prot.indices)}
    for c, cg in enumerate(copies):
        for r, res in enumerate(cg.residues):
            for ix in res.atoms.indices:
                k = pos_in_prot[ix]; a_copy[k] = c; a_res[k] = r
    upper_ri, leaflet_source = upper_leaflet_resindices(u, args, args.start)
    upper_res = u.residues[upper_ri]
    types = list(args.lipids)
    tmap = {t: i for i, t in enumerate(types)}
    keep = np.array([rn in tmap for rn in upper_res.resnames])
    upper_res = upper_res[keep]
    lip = upper_res.atoms
    mol_of_resindex = {ri: m for m, ri in enumerate(upper_res.resindices)}
    l_mol = np.array([mol_of_resindex[ri] for ri in lip.resindices])
    mol_type = np.array([tmap[rn] for rn in upper_res.resnames])
    l_type = mol_type[l_mol]
    # beads used for the legacy leaflet mean z: the published code averages over the atoms of its
    # leaflet AtomGroup, which is headgroup beads only when the leaflet comes from LeafletFinder
    # (fresh run) and whole molecules when it is rebuilt from leaflet_info.pickle
    if getattr(args, 'leaflet_z', 'headgroup') == 'headgroup':
        hg = set(HEADGROUP_SEL.split()[1:])
        l_zmask = np.array([n in hg for n in lip.names])
    else:
        l_zmask = np.ones(len(lip), bool)
    subset = None
    if args.subset:
        a, b = map(int, args.subset.split(':'))
        subset = (resids >= a) & (resids <= b)
    return dict(u=u, prot=prot, copies=copies, nc=len(copies), nres=nres, resids=resids,
                a_copy=a_copy, a_res=a_res, lip=lip, upper_res=upper_res, l_mol=l_mol,
                l_type=l_type, l_zmask=l_zmask, mol_type=mol_type, types=types, ntype=len(types),
                subset=subset, segs=segs, leaflet_source=leaflet_source,
                prot_res=prot.residues)


# ----------------------------------------------------------------------------- per frame
def minimg(d, box):
    return d - box * np.round(d / box)


def legacy_min_image(d, box):
    """per-dimension minimum image as written in the published code (> half box: subtract)"""
    d = d.copy()
    for k in range(3):
        hi = d[:, k] > box[k] * 0.5; lo = d[:, k] < -box[k] * 0.5
        d[hi, k] -= box[k]; d[lo, k] += box[k]
    return d


def process_frame(S, cutoff):
    u = S['u']; box = u.dimensions[:3].astype(float)
    ppos = S['prot'].positions; lpos = S['lip'].positions
    pairs = capped_distance(ppos, lpos, max_cutoff=cutoff, box=u.dimensions, return_distances=False)
    nc, nres, nt = S['nc'], S['nres'], S['ntype']
    nmol = len(S['upper_res'])
    res_exact = np.zeros((nc, nres, nt), np.int16)
    res_legacy = np.zeros((nc, nres, nt), np.int16)
    mol_full = np.zeros((nc, nt), np.int16)
    mol_subset = np.zeros((nc, nt), np.int16)
    mol_legacy = np.zeros((nc, nt), np.int16)
    if len(pairs) == 0:
        return res_exact, res_legacy, mol_full, mol_subset, mol_legacy
    pc = S['a_copy'][pairs[:, 0]]; pr = S['a_res'][pairs[:, 0]]; lm = S['l_mol'][pairs[:, 1]]
    # unique (copy, residue, molecule)
    key = (pc * nres + pr).astype(np.int64) * nmol + lm
    uk = np.unique(key)
    m = uk % nmol; cr = uk // nmol; c = cr // nres; r = cr % nres
    t = S['mol_type'][m]
    np.add.at(res_exact, (c, r, t), 1)
    # distinct molecules per copy (full protein and subset)
    k2 = np.unique(c.astype(np.int64) * nmol + m)
    np.add.at(mol_full, (k2 // nmol, S['mol_type'][k2 % nmol]), 1)
    if S['subset'] is not None:
        ins = S['subset'][r]
        k3 = np.unique(c[ins].astype(np.int64) * nmol + m[ins])
        np.add.at(mol_subset, (k3 // nmol, S['mol_type'][k3 % nmol]), 1)
    # ---- legacy filters on the (copy, residue, molecule) triples
    prot_res = S['prot_res']
    res_com = prot_res.center_of_mass(compound='residues')            # (nc*nres, 3)
    gi = c * nres + r                                                    # global residue index
    mol_com = S['upper_res'].atoms.center_of_mass(compound='residues')  # (nmol, 3)
    lz = S['lip'].positions[:, 2]
    type_zavg = np.array([lz[(S['l_type'] == k) & S['l_zmask']].mean() if np.any((S['l_type'] == k) & S['l_zmask']) else np.nan
                          for k in range(nt)])
    zres = S['res_zavg_fn']()
    ok_z = np.abs(zres[gi] - type_zavg[t]) <= 15.0
    dcom = legacy_min_image(res_com[gi] - mol_com[m], box)
    ok_com = np.sqrt((dcom ** 2).sum(1)) <= cutoff + 8.0
    okL = ok_z & ok_com
    np.add.at(res_legacy, (c[okL], r[okL], t[okL]), 1)
    # legacy unique molecules: protein-selection COM, xy distance <= 10 A, any bead within cutoff
    sel_com = S['legacy_sel_com_fn']()                                   # (nc, 3)
    cm = np.unique(c.astype(np.int64) * nmol + m)
    cc = cm // nmol; mm = cm % nmol
    if S['subset'] is not None:
        # the legacy unique count used the selected residues only: keep molecules touching the subset
        ins = S['subset'][r]
        cm_s = set((c[ins].astype(np.int64) * nmol + m[ins]).tolist())
        keepm = np.array([x in cm_s for x in cm.tolist()])
        cc, mm = cc[keepm], mm[keepm]
    dxy = mol_com[mm, :2] - sel_com[cc, :2]
    dxy = dxy - box[:2] * np.round(dxy / box[:2])
    okU = np.sqrt((dxy ** 2).sum(1)) <= 10.0
    np.add.at(mol_legacy, (cc[okU], S['mol_type'][mm[okU]]), 1)
    return res_exact, res_legacy, mol_full, mol_subset, mol_legacy


def attach_helpers(S):
    prot_res = S['prot_res']
    # atom index lists per protein residue, for unweighted mean z (legacy res_z_avg)
    idx_lists = [res.atoms.ix for res in prot_res]
    S['res_atoms_ix'] = idx_lists
    u = S['u']
    def zfn():
        pos = u.atoms.positions
        return np.array([pos[ix, 2].mean() for ix in idx_lists])
    S['res_zavg_fn'] = zfn
    if S['subset'] is not None:
        sel = [cg.residues[S['subset']].atoms for cg in S['copies']]
    else:
        sel = [cg for cg in S['copies']]
    S['legacy_sel_com_fn'] = lambda: np.array([g.center_of_mass() for g in sel])
    return S


# ----------------------------------------------------------------------------- workers
def worker_init(args_dict):
    args = argparse.Namespace(**args_dict)
    S = build_state(args)
    _G['S'] = attach_helpers(S); _G['args'] = args


def run_chunk(chunk):
    S = _G['S']; args = _G['args']; u = S['u']
    ci, frames = chunk
    out = os.path.join(args.out, 'chunks', f'chunk_{ci:05d}.npz')
    if os.path.exists(out):
        return ci, 0, 0.0
    t0 = time.time()
    R = {k: [] for k in ('res_exact', 'res_legacy', 'mol_full', 'mol_subset', 'mol_legacy')}
    times = []
    for f in frames:
        ts = u.trajectory[f]; times.append(ts.time)
        a = process_frame(S, args.cutoff)
        for k, v in zip(R, a):
            R[k].append(v)
    tmp = out + '.tmp.npz'
    np.savez_compressed(tmp, frames=np.array(frames), time_ps=np.array(times),
                        **{k: np.stack(v) for k, v in R.items()})
    os.replace(tmp, out)
    return ci, len(frames), time.time() - t0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--psf', required=True); ap.add_argument('--xtc', required=True)
    ap.add_argument('--start', type=int, required=True); ap.add_argument('--stop', type=int, required=True)
    ap.add_argument('--step', type=int, default=1)
    ap.add_argument('--target', default='DPG3')
    ap.add_argument('--lipids', nargs='+', required=True)
    ap.add_argument('--subset', default=None, help='residue range (PSF numbering), e.g. 65:103')
    ap.add_argument('--leaflet-pickle', default=None)
    ap.add_argument('--cutoff', type=float, default=6.0)
    ap.add_argument('--leaflet-z', choices=['headgroup', 'all'], default='headgroup',
                    help='beads for the legacy leaflet mean z (headgroup = published fresh run)')
    ap.add_argument('--out', required=True)
    ap.add_argument('--workers', type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument('--chunk', type=int, default=500)
    ap.add_argument('--frames', type=int, nargs='*', default=None, help='explicit frame list (validation)')
    args = ap.parse_args()
    os.makedirs(os.path.join(args.out, 'chunks'), exist_ok=True)

    S = attach_helpers(build_state(args))
    nframes = len(S['u'].trajectory)
    frames = args.frames if args.frames else list(range(args.start, min(args.stop, nframes), args.step))
    meta = dict(psf=args.psf, xtc=args.xtc, start=args.start, stop=args.stop, step=args.step,
                target=args.target, lipids=S['types'], subset=args.subset, cutoff=args.cutoff, leaflet_z=args.leaflet_z,
                segids=S['segs'], resids=S['resids'].tolist(), n_upper_molecules=int(len(S['upper_res'])),
                upper_by_type={t: int((S['mol_type'] == i).sum()) for i, t in enumerate(S['types'])},
                leaflet_source=S['leaflet_source'], n_frames_total=nframes, n_frames_analysed=len(frames),
                mdanalysis=mda.__version__)
    json.dump(meta, open(os.path.join(args.out, 'meta.json'), 'w'), indent=1)
    print(json.dumps({k: v for k, v in meta.items() if k != 'resids'}, indent=1), flush=True)
    del S
    chunks = [(i, frames[a:a + args.chunk]) for i, a in enumerate(range(0, len(frames), args.chunk))]
    todo = [ch for ch in chunks if not os.path.exists(os.path.join(args.out, 'chunks', f'chunk_{ch[0]:05d}.npz'))]
    print(f'{len(frames)} frames in {len(chunks)} chunks; {len(todo)} to run on {args.workers} workers', flush=True)
    t0 = time.time(); done = 0
    if args.workers <= 1:
        worker_init(vars(args))
        for ch in todo:
            ci, n, dt = run_chunk(ch); done += n
            print(f'chunk {ci}: {n} frames in {dt:.1f} s', flush=True)
    else:
        import multiprocessing as mp
        ctx = mp.get_context('spawn')
        with ctx.Pool(args.workers, initializer=worker_init, initargs=(vars(args),)) as pool:
            for ci, n, dt in pool.imap_unordered(run_chunk, todo):
                done += n
                el = time.time() - t0
                print(f'chunk {ci}: {n} frames in {dt:.1f} s | {done}/{sum(len(c[1]) for c in todo)} frames, '
                      f'{el/60:.1f} min elapsed', flush=True)
    print('finished; run merge_and_export.py next', flush=True)


if __name__ == '__main__':
    main()
