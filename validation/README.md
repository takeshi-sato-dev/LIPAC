# Validation of LIPAC 3

`test_exact_contacts.py` compares every contact function of
`stage1_contact_analysis/core/contact_calculator.py` (residue–lipid contacts, molecules in
contact, protein–protein contacts) with a brute-force count on synthetic membranes in
which proteins and lipids are split across the periodic boundary. LIPAC 3 passes;
the prescreened functions of LIPAC 2 fail.

    python validation/test_exact_contacts.py

Checks on Martini trajectories (29 September 2026):

1. The serial path (`--no-parallel`) and the parallel path return identical Stage 1
   output. LIPAC 2 fails this check, because its parallel path held the lipid
   coordinates at the first analyzed frame.
2. The contact numbers, the numbers of molecules in contact and the binding states of
   Stage 1 agree with a brute-force count by MDAnalysis `distance_array` (EGFR
   transmembrane–juxtamembrane construct, four copies, 1200 residue rows over six random
   frames, and 360 values over ten consecutive frames).
3. `stage1_fast` and `stage1_contact_analysis` return the same counts
   (`stage1_fast/test_against_lipac.py`, and the same trajectory frames as in 2).
4. `run_stage2_calibrated.py` on real CHOL contact numbers (3000 frames × 4 copies,
   autocorrelation time about 170 frames): 12 of 12 data sets in which the binding state
   was shifted in time are classified "no detectable effect"; a planted linear effect of
   +8 contacts is classified linear, and a planted cooperative effect of +15 contacts in
   40% of bound frames is classified cooperative. Effects of +3 contacts (linear) and of
   +8 contacts in 30% of bound frames (cooperative) are not detected at this trajectory
   length. The uncalibrated rule of the Bayesian models (delta WAIC > 2) classifies the
   same null data sets as cooperative in 95% of cases.
