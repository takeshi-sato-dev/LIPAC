# BayesianLipidAnalysis: Bayesian Analysis of Lipid-Protein Interactions

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.XXXXX/status.svg)](https://doi.org/10.21105/joss.XXXXX)

## Overview

This code is a Python package for Bayesian statistical analysis of lipid-protein interactions from molecular dynamics simulations. It provides comprehensive tools to analyze how specific lipids (such as GM3 ganglioside) affect protein-protein and protein-lipid contact patterns in membrane systems.

## Features

- **Exploratory Data Analysis**: Statistical comparison of contact patterns between systems with and without lipid modifiers
- **Bayesian Regression Models**: Quantitative assessment of lipid effects on protein interactions
- **Hierarchical Bayesian Models**: Protein-specific effect analysis accounting for inter-protein variability
- **Residue-Level Analysis**: Identification of competitive binding sites and interface regions
- **Lipid-Specific Analysis**: Differential effects on various lipid types (cholesterol, phospholipids, etc.)
- **Comprehensive Visualization**: Publication-ready plots with convergence diagnostics

## Installation

### Prerequisites

- Python 3.10 or higher
- pip package manager

### Install from source

```bash
# Create dedicated environment for Bayesian analysis
cd stage2_contact_analysis
python3.10 -m venv bayesian_env
source bayesian_env/bin/activate  # On Windows: bayesian_env\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

### Environment Setup

The Bayesian analysis stage requires specific packages that may conflict with Stage 1 dependencies. We recommend using a separate virtual environment:

```bash
# If bayesian_env doesn't exist, create it
python3.10 -m venv bayesian_env
source bayesian_env/bin/activate

# Upgrade pip and install requirements
pip install --upgrade pip
pip install -r requirements.txt

# Verify installation
python -c "import pymc; import arviz; print('Bayesian environment ready!')"
```

### Dependencies

- numpy >= 1.20.0
- pandas >= 1.3.0
- matplotlib >= 3.4.0
- seaborn >= 0.11.0
- scipy >= 1.7.0
- pymc >= 5.0.0
- arviz >= 0.14.0

## Quick Start

### Basic Usage

```bash
python main.py --with-lipid-data data/with_target_lipid.csv --without-lipid-data data/without_target_lipid.csv --output-dir results/
```

### Advanced Options

```bash
python bayesian_lipid_analysis_main.py \
    --with-gm3-data data/with_modifier.csv \
    --no-gm3-data data/without_modifier.csv \
    --output-dir results/ \
    --mcmc-samples 3000 \
    --chains 6 \
    --skip-hierarchical \
    --verbose
```

## Input Data Format

The input CSV files should contain contact complementarity data with the following columns:

### Required Columns
- `protein`: Protein identifier (e.g., "Protein_1", "Protein_2")
- `residue`: Residue identifier
- `protein_contact`: Protein-protein contact frequency
- `lipid_contact`: Total lipid-protein contact frequency
- `ratio`: Lipid-to-protein contact ratio

### Optional Columns
- `DPG3_contact`: Contact frequency with the lipid modifier (e.g., GM3)
- `CHOL_contact`: Cholesterol contact frequency
- `DOPS_contact`: DOPS lipid contact frequency
- `DIPC_contact`: DIPC lipid contact frequency
- `DPSM_contact`: DPSM lipid contact frequency
- `partner_protein`: Partner protein in dimer (for pair analysis)

## Analysis Pipeline

1. **Data Loading and Preprocessing**: Merges and aligns data from systems with and without lipid modifiers
2. **Exploratory Analysis**: Statistical tests and visualization of contact differences
3. **Bayesian Regression**: Quantifies the relationship between lipid binding and contact changes
4. **Hierarchical Modeling**: Accounts for protein-specific variations in lipid effects
5. **Residue-Level Analysis**: Identifies competitive binding sites at protein interfaces
6. **Lipid-Specific Analysis**: Assesses differential effects on various membrane lipids
7. **Report Generation**: Creates comprehensive summary of all findings

## Output Files

The analysis generates organized output in the specified directory:

```
output_dir/
├── figures/                      # Exploratory analysis plots
├── models/                        # Bayesian model results and diagnostics
├── residue_analysis/              # Residue-level analysis results
├── lipid_specific_analysis/       # Lipid-specific effect analysis
├── combined_data.csv              # Merged dataset
├── matched_data.csv               # Aligned comparison data
├── exploratory_analysis_results.txt
├── bayesian_analysis_results.txt
├── hierarchical_model_results.txt
└── analysis_summary.txt          # Comprehensive summary report
```

## Command-Line Options

| Option | Default | Description |
|--------|---------|-------------|
| `--with-gm3-data` | `bayesian_lipid_analysis/with_gm3/contact_complementarity.csv` | CSV file with lipid modifier |
| `--no-gm3-data` | `bayesian_lipid_analysis/no_gm3/contact_complementarity.csv` | CSV file without lipid modifier |
| `--output-dir` | `bayesian_lipid_analysis/bayesian_analysis` | Output directory |
| `--mcmc-samples` | 2000 | Number of MCMC samples |
| `--tune-samples` | 1000 | Number of tuning samples |
| `--chains` | 4 | Number of MCMC chains |
| `--seed` | 42 | Random seed for reproducibility |
| `--skip-exploratory` | False | Skip exploratory analysis |
| `--skip-bayesian` | False | Skip Bayesian regression |
| `--skip-hierarchical` | False | Skip hierarchical model |
| `--skip-residue-analysis` | False | Skip residue-level analysis |
| `--skip-lipid-specific` | False | Skip lipid-specific analysis |
| `--verbose` | False | Display detailed output |

## Example Analysis

### Preparing Input Data

Your CSV files should be structured as follows:

```csv
protein,residue,protein_contact,lipid_contact,ratio,DPG3_contact,CHOL_contact
Protein_1,R101,0.85,0.15,0.176,0.05,0.08
Protein_1,K102,0.72,0.28,0.389,0.12,0.10
...
```

### Running the Analysis

```python
# For programmatic use
from bayesian_lipid_analysis.data_loader import load_data
from bayesian_lipid_analysis.bayesian_models import bayesian_regression_model

# Load data
combined_df, matched_df = load_data('with_gm3.csv', 'no_gm3.csv')

# Run Bayesian analysis
results = bayesian_regression_model(matched_df, 'output_dir/')
```

## Interpreting Results

### Bayesian Regression Output
- **β < 0**: Lipid modifier decreases protein-protein contacts (competitive binding)
- **β > 0**: Lipid modifier increases protein-protein contacts (cooperative effect)
- **95% HDI**: Credible interval for the effect size

### Residue-Level Analysis
- **Competitive binding sites**: Residues where lipid and protein compete for binding
- **Enrichment > 1**: Lipid preferentially binds at protein interfaces
- **p < 0.05**: Statistically significant overlap between binding sites

## Citation

If you use BayesianLipidAnalysis in your research, please cite:

```bibtex
@article{sato2024bayesian_lipid_analysis,
  title={BayesianLipidAnalysis: Bayesian Analysis of Lipid-Protein Interactions from Molecular Dynamics Simulations},
  author={Sato, Takeshi},
  journal={Journal of Open Source Software},
  year={2024},
  volume={X},
  number={XX},
  pages={XXXXX},
  doi={10.21105/joss.XXXXX}
}
```

## Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Support

For questions and support, please open an issue on GitHub or contact the author at [your.email@example.com].

## Acknowledgments

This work was supported by Kyoto Pharmaceutical University. We thank the molecular dynamics simulation community for valuable feedback and suggestions.