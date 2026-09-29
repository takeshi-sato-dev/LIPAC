# LIPAC: Lipid-Protein Analysis with Causal Inference

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18080231.svg)](https://doi.org/10.5281/zenodo.18080231)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue)](https://www.python.org/downloads/)

## Notice of errors in LIPAC 1 and 2 (LIPAC 3, 29 September 2026)

Stage 1 of LIPAC 1 and 2, which covers every commit up to a010a80 (9 April 2026),
contains two errors in the calculation of lipid–protein contacts. Both are corrected in
LIPAC 3.

1. **Lipid coordinates held at the first analyzed frame.** In the parallel version of
   Stage 1, the upper-leaflet AtomGroup selected at the first analyzed frame was passed
   to every worker process. An AtomGroup passed to another process is pickled together
   with the Universe it belongs to. The lipids carried the coordinates of
   the first analyzed frame in every later frame, while the protein coordinates were
   updated correctly. The binding state of the target lipid and every lipid contact
   number were affected. The serial version was not affected. The worker now receives
   the atom indices of the leaflet and selects the lipids from its own Universe.
2. **Prescreens that discarded real contacts.** Two prescreens were applied before the
   6 Å bead–bead criterion: a comparison of the mean height of a residue with the mean
   height of a lipid type (15 Å), and a center-of-mass distance (14 Å). A molecule split
   across the periodic boundary has a center of mass near the middle of the box, and
   every contact with such a molecule was discarded. The prescreens are removed, and a
   residue and a lipid molecule are in contact whenever any bead of one lies within
   6 Å of any bead of the other, under the minimum image. The unique-molecule count and
   the protein–protein contacts are calculated in the same way.

The results of *J. Chem. Inf. Model.* DOI 10.1021/acs.jcim.5c02497 were produced with the
affected parallel version, and the author has requested the retraction of that article.

**Checks.** `validation/README.md` lists the checks of this version: the serial and the
parallel paths return identical output, the contacts agree with a brute-force count on
real trajectories, and `validation/test_exact_contacts.py` compares every contact
function with a brute-force count on synthetic membranes whose molecules are split across
the periodic boundary (LIPAC 2 fails this test). `stage1_fast/` adds a fast parallel
implementation of the same contact definition for analysis at full time resolution.

**Stage 2.** The Bayesian models treat every frame as an independent observation.
Frames of a molecular dynamics trajectory are strongly correlated, and on real data a
credible interval of a per-copy effect can exclude zero when no effect is present. In
the same data, the mixture model beats the linear model when the binding state is
shifted in time against the contact numbers, because the contact numbers are not
normally distributed. `run_stage2_calibrated.py` (module
`stage2_contact_analysis/analysis/calibration.py`) gives the classification that holds
for correlated frames: the effect of binding is estimated with a block bootstrap over
time blocks chosen from the autocorrelation time, and both the effect and the
mixture-over-linear comparison are tested against binding states shifted circularly in
time within each copy (at least 100 shifts; 200 by default). A lipid type is classified
as linear or cooperative only when it passes these tests, and the uncalibrated rule
(delta WAIC > 2) of the Bayesian models should not be used on its own.



LIPAC is a Python package for comprehensive analysis of lipid-protein interactions from molecular dynamics simulations with integrated causal inference capabilities.

## Features

- **Causal Inference**: Determine causal effects of specific lipid binding on membrane organization
- **Dual Metrics**: Analyze both residue-level contacts and unique molecule counts
- **High Performance**: Parallel processing with checkpoint/restart capability
- **Comprehensive Visualization**: Publication-ready plots and statistical reports
- **Flexible**: Supports multiple MD trajectory formats (GROMACS, CHARMM, AMBER)

## Installation

### From PyPI
```bash
pip install lipac
```

### From source
```bash
git clone https://github.com/takeshi-sato-dev/LIPAC.git
cd lipac
pip install -e .
```

### Dependencies
- Python ≥ 3.10
- MDAnalysis ≥ 2.0.0
- NumPy ≥ 1.20.0
- PyMC ≥ 5.0.0
- Matplotlib ≥ 3.3.0
- Pandas ≥ 1.3.0
- Seaborn ≥ 0.11.0

## Test Data

LIPAC includes sample datasets for testing and demonstration purposes. Download the test data from Zenodo:

[![Test Data DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17211644.svg)](https://doi.org/10.5281/zenodo.17211644)

### Test Dataset Contents

The test dataset includes:
- **`test_system_with_mediator.psf/xtc`**: System containing target lipid (DPG3/GM3) that prevents protein dimerization
- **`test_system_without_mediator.psf/xtc`**: System without target lipid, allowing protein dimerization
- Sample trajectories for demonstrating causal inference analysis

### Download and Setup

```bash
# Create test directory
mkdir test

# Download test data from Zenodo
wget -O test/test_system_with_mediator.psf "https://zenodo.org/records/17211644/files/test_system_with_mediator.psf"
wget -O test/test_system_without_mediator.psf "https://zenodo.org/records/17211644/files/test_system_without_mediator.psf"
wget -O test/test_trajectory_with_mediator.xtc "https://zenodo.org/records/17211644/files/test_trajectory_with_mediator.xtc"
wget -O test/test_trajectory_without_mediator.xtc "https://zenodo.org/records/17211644/files/test_trajectory_without_mediator.xtc"
```

Or use the provided download script:
```bash
python scripts/download_test_data.py
```

## Quick Start

### Basic Contact Analysis
```python
from lipac import ContactAnalysis

# Initialize analyzer
analyzer = ContactAnalysis(
    topology='system.gro',
    trajectory='traj.xtc',
    target_lipid='GM3'  # Specify target lipid for competition analysis
)

# Run analysis
results = analyzer.run(
    start=0, 
    stop=100000, 
    step=50,
    n_jobs=8  # Number of parallel workers
)

# Save results
results.save('contact_results.pkl')
```

### Causal Analysis
```python
from lipac import CausalAnalysis

# Load contact data
causal = CausalAnalysis('contact_results.pkl')

# Compute causal effects
effects = causal.compute_causal_effects(
    n_samples=2000,
    n_chains=4,
    target_lipid='GM3'
)

# Generate visualizations
causal.plot_causal_effects(output_dir='results/')
causal.plot_competition_analysis(output_dir='results/')
```

### Advanced Usage with Custom Parameters
```python
from lipac import LIPACPipeline

# Complete pipeline with custom parameters
pipeline = LIPACPipeline(
    topology='system.gro',
    trajectory='traj.xtc',
    config={
        'target_lipid': 'GM3',
        'contact_cutoff': 6.0,  # Å
        'xy_plane_cutoff': 10.0,  # Å for optimization
        'frame_step': 50,
        'n_jobs': 16,
        'mcmc_samples': 4000,
        'mcmc_chains': 4
    }
)

# Run complete analysis
results = pipeline.run_full_analysis()

# Generate report
pipeline.generate_report('analysis_report.html')
```

## Documentation

Full documentation is available at [https://lipac.readthedocs.io](https://lipac.readthedocs.io)

### Tutorials
- [Basic Usage](docs/tutorials/basic_usage.md)
- [Causal Analysis](docs/tutorials/causal_analysis.md)
- [Visualization](docs/tutorials/visualization.md)
- [Performance Optimization](docs/tutorials/optimization.md)

## Example Data

Example trajectories and analysis scripts are available in the `examples/` directory:
- `examples/simple_membrane/`: Single protein in POPC membrane
- `examples/complex_membrane/`: Transmembrane peptide in complex lipid membrane with GM3
- `examples/benchmark/`: Performance benchmark scripts

**Test Data**: Complete test datasets including MD trajectories and expected outputs are available on Zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17211644.svg)](https://doi.org/10.5281/zenodo.17211644)

## Citation

If you use LIPAC in your research, please cite:

```bibtex
@article{lipac_2026,
  title = {A Computational Framework for Causal Inference in Molecular Dynamics Analysis of Lipid-Protein Interactions},
  author = {Sato, Takeshi},
  journal = {Journal of Chemical Information and Modelling},
  year = {2026},
  volume = {X},
  number = {XX},
  pages = {XXXX},
  doi = {10.1021/acs.jcim.5c02497}
}
```

## Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## License

LIPAC is licensed under the MIT License. See [LICENSE](LICENSE) for details.

## Support

- **Issues**: [GitHub Issues](https://github.com/takeshi-sato-dev/lipac/issues)
- **Discussions**: [GitHub Discussions](https://github.com/takeshi-sato-dev/lipac/discussions)
- **Email**: takeshi@mb.kyoto-phu.ac.jp

## Acknowledgments

We acknowledge contributions and support from Kyoto Pharmaceutical University Fund for the Promotion of Collaborative Research. This work was partially supported by JSPS KAKENHI Grant Number 21K06038.

## Related Projects

- [MDAnalysis](https://www.mdanalysis.org/): Trajectory analysis framework
- [PyMC](https://www.pymc.io/): Bayesian statistical modeling
- [FATSLiM](https://github.com/FATSLiM/fatslim): Membrane analysis tools

## v2: Mixture Causal Model

Bayesian mixture causal model for detecting cooperative lipid dynamics in MD simulations. When a linear model produces sign-inconsistent β values across replicates, this may reflect multimodal responses rather than insufficient sampling. The mixture model decomposes bound-state observations into two latent subpopulations—one with a direct binding effect only (β_direct) and one with an additional cooperative effect (β_direct + β_coop)—with mixing probability π estimated from the data. Automated model selection via WAIC classifies each lipid species as linear or cooperative.

**Reference:** Sato, T. "Nonlinear Causal Inference for Cooperative Lipid Dynamics in Molecular Dynamics Simulations" (submitted to *J. Chem. Theory Comput.*)

### Usage

```bash
python analysis/mixture_causal_analysis.py --input stage1_output.csv
```

### Recommended Workflow

1. Run the linear LIPAC analysis first to identify causal effects and assess consistency across protein copies.
2. For any lipid type where β values show inconsistent signs across copies, apply the mixture model.
3. Classify: ΔWAIC > 2 **and** 95% CI of β_coop excludes zero → **cooperative**; otherwise → **linear**.
4. Report π × β_coop alongside individual estimates of β_coop and π.
5. Interpret β_direct and β_coop quantitatively only when β_coop/σ > 2.

### Dependencies

PyMC, ArviZ, NumPy, pandas
