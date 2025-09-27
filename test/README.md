# Test Data for GM3-EGFR Analysis

This directory contains representative test datasets for code demonstration.

## Files

- `test_system.psf`: Topology file (61.57 MB)
- `test_trajectory.xtc`: Trajectory with 199 frames (344.59 MB)

## Important Note

**This is a subset of the full trajectory for demonstration purposes.**

- Full analysis uses frames 20000-80000 with step 20 (3000 frames total)
- Test data uses frames 80000-89900 with step 50 (199 frames)
- Full trajectory size: ~110 GB (8-10 μs simulation)

## Generation

These files were extracted from the full MD trajectory:
- Original PSF: step5_assembly.psf
- Original XTC: md_wrapped.xtc
- Original trajectory: 88799 frames

## Usage

```python
import MDAnalysis as mda

# Load test data
u = mda.Universe('test_data/test_system.psf', 'test_data/test_trajectory.xtc')
print(f'Loaded {len(u.atoms)} atoms, {len(u.trajectory)} frames')

# Note: Statistical results from this subset may differ from full analysis
# For publication results, full trajectory (frames 20000-80000) was used
```

## Full Data Access

Full trajectories (110 GB) are available upon request from the authors.
Pre-computed contact matrices and analysis results from the full trajectory
are provided in the `results/` directory for validation.
