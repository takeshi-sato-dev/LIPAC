# LIPAC: Lipid-Protein Analysis with Causal Inference

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16891334.svg)](https://doi.org/10.5281/zenodo.16891334)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue)](https://www.python.org/downloads/)
[![JOSS](https://joss.theoj.org/papers/10.21105/joss.XXXXX/status.svg)](https://doi.org/10.21105/joss.XXXXX)

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

**Test Data**: Complete test datasets including MD trajectories and expected outputs are available on Zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16891479.svg)](https://doi.org/10.5281/zenodo.16891479)

## Citation

If you use LIPAC in your research, please cite:

```bibtex
@article{lipac_2025,
  title = {LIPAC: Lipid-Protein Analysis with Causal Inference for Molecular Dynamics Simulations},
  author = {Sato, Takeshi},
  journal = {Journal of Open Source Software},
  year = {2025},
  volume = {X},
  number = {XX},
  pages = {XXXX},
  doi = {10.21105/joss.XXXXX}
}
```

## Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## License

LIPAC is licensed under the MIT License. See [LICENSE](LICENSE) for details.

## Support

- **Issues**: [GitHub Issues](https://github.com/takeshi-sato-dev/lipac/issues)
- **Discussions**: [GitHub Discussions](https://github.com/takeshi-sato-dev/lipac/discussions)
- **Email**: takeshi.sato@example.com

## Acknowledgments

We acknowledge contributions and support from Kyoto Pharmaceutical University Fund for the Promotion of Collaborative Research. This work was partially supported by JSPS KAKENHI Grant Number 21K06038.

## Related Projects

- [MDAnalysis](https://www.mdanalysis.org/): Trajectory analysis framework
- [PyMC](https://www.pymc.io/): Bayesian statistical modeling
- [FATSLiM](https://github.com/FATSLiM/fatslim): Membrane analysis tools