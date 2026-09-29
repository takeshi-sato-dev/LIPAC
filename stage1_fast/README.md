# Stage 1, fast and parallel (for full time resolution)

`fast_stage1.py` computes the lipid–protein contacts of Stage 1 for every frame of a
trajectory. It uses the contact definition of `stage1_contact_analysis` (LIPAC 3
and later): a residue and a lipid molecule are in contact when any bead of one lies within 6 Å
of any bead of the other, under the minimum image. The search uses MDAnalysis
`capped_distance`, which takes about 13 ms per frame for four copies of a 50-residue
construct in an upper leaflet of 2600 molecules (reading the trajectory not included).

Frames are processed in chunks by a pool of workers. Each worker opens the trajectory
once and selects every atom group from its own Universe, and every chunk is written to
its own file, which lets an interrupted run resume where it stopped.

## Run

```bash
pip install MDAnalysis pandas networkx
python stage1_fast/fast_stage1.py --psf step5_assembly.psf --xtc step7_production.xtc \
    --start 20000 --stop 80000 --step 1 --lipids CHOL DPSM DIPC DPG3 DOPS \
    --out runs/system_with_GM3 --workers 8
python stage1_fast/merge_and_export.py --out runs/system_with_GM3
```

`--subset 65:103` restricts the tables to a residue range (the per-residue arrays keep
the full protein). `--leaflet-pickle` reads the upper leaflet from an existing
`leaflet_info.pickle`; otherwise the leaflet is detected once with LeafletFinder at the
first analyzed frame, and only the membership is kept fixed.

## Output

| file | content |
|---|---|
| `stage1_fast_all.npz` | per frame: residue × lipid-type contact counts (`res_exact`), molecules in contact per copy (`mol_full`, `mol_subset`), frames and times |
| `causal_data_exact_full.csv` | Stage 2 input table, current definition, full protein |
| `causal_data_exact_subset.csv` | the same for the residue subset |
| `causal_data_legacy_*.csv` | the definition of LIPAC 1 and 2, with its two prefilters, reproduced only for comparison with earlier results |

The tables carry the columns of the Stage 1 output (`frame`, `protein`,
`target_lipid_bound`, `<lipid>_contacts`, `<lipid>_unique_molecules`) and in addition
`time_ps`, `<target>_contacts` and `<target>_molecules`, the latter two as continuous
measures of target-lipid contact.

## Tests

`python stage1_fast/test_against_lipac.py` compares the exact counts with the contact
functions of `stage1_contact_analysis/core/contact_calculator.py` on synthetic membranes
whose molecules straddle the periodic boundary. `validate_real.py` compares two Stage 1
tables row by row.
