# Stage 1: Contact Analysis

Protein-protein and lipid-protein contact analysis tool for molecular dynamics trajectories.

## Installation

```bash
# Clone the repository
git clone <repository-url>
cd lipac

# Install dependencies
pip install -r requirements.txt
```

## Usage

### 1. Configure your analysis settings

Edit `stage1_contact_analysis/config.py` to set your file paths and parameters:

```python
# Input file paths
DEFAULT_WITH_LIPID_PSF = 'path/to/your/with_target_lipid_system.psf'
DEFAULT_WITH_LIPID_XTC = 'path/to/your/with_target_lipid_trajectory.xtc'
DEFAULT_WITHOUT_LIPID_PSF = 'path/to/your/without_target_lipid_system.psf'
DEFAULT_WITHOUT_LIPID_XTC = 'path/to/your/without_target_lipid_trajectory.xtc'

# Output directories
DEFAULT_WITH_LIPID_OUTPUT = 'lipac_results/with_target_lipid'
DEFAULT_WITHOUT_LIPID_OUTPUT = 'lipac_results/without_target_lipid'

# Analysis parameters
START_FRAME = 0      # First frame to analyze
STOP_FRAME = 1000    # Last frame to analyze
STEP_FRAME = 1       # Frame interval
BATCH_SIZE = 50      # Frames per batch for parallel processing
MIN_CORES = 2        # Minimum CPU cores to use
```

### 2. Run the analysis

Basic execution using config.py settings:

```bash
python stage1_contact_analysis/main.py
```

Or override specific parameters via command line:

```bash
# Override input files (useful for testing with different data)
python stage1_contact_analysis/main.py \
    --with-lipid-psf test_data/test_system_with_target_lipid.psf \
    --with-lipid-xtc test_data/test_trajectory_with_target_lipid.xtc \
    --without-lipid-psf test_data/test_system_without_target_lipid.psf \
    --without-lipid-xtc test_data/test_trajectory_without_target_lipid.xtc

# Override frame range for quick testing
python stage1_contact_analysis/main.py --start 0 --stop 10 --debug

# Skip plot generation for faster testing
python stage1_contact_analysis/main.py --skip-plots --debug
```

### 3. Check the output

The analysis generates:

- `contact_complementarity.csv` - Protein-protein and lipid-protein contact statistics
- `residue_contacts.csv` - Per-residue contact information
- Various plots in the output directories (unless `--skip-plots` is used)

## Command Line Options

| Option | Description | Default (from config.py) |
|--------|-------------|---------------------------|
| `--with-lipid-psf` | PSF file for system with target lipid | `DEFAULT_WITH_LIPID_PSF` |
| `--with-lipid-xtc` | Trajectory file for system with target lipid | `DEFAULT_WITH_LIPID_XTC` |
| `--without-lipid-psf` | PSF file for system without target lipid | `DEFAULT_WITHOUT_LIPID_PSF` |
| `--without-lipid-xtc` | Trajectory file for system without target lipid | `DEFAULT_WITHOUT_LIPID_XTC` |
| `--start` | First frame to analyze | `START_FRAME` |
| `--stop` | Last frame to analyze | `STOP_FRAME` |
| `--step` | Frame interval | `STEP_FRAME` |
| `--batch-size` | Frames per batch | `BATCH_SIZE` |
| `--cores` | Minimum CPU cores to use | `MIN_CORES` |
| `--debug` | Enable debug mode | `False` |
| `--skip-plots` | Skip plot generation | `False` |
| `--no-parallel` | Disable parallel processing | `False` |

## Testing with Sample Data

Test data is provided in the `test_data/` directory. To run a quick test:

```bash
# Quick test with 10 frames
python stage1_contact_analysis/main.py \
    --with-lipid-psf test_data/test_system_with_target_lipid.psf \
    --with-lipid-xtc test_data/test_trajectory_with_target_lipid.xtc \
    --without-lipid-psf test_data/test_system_without_target_lipid.psf \
    --without-lipid-xtc test_data/test_trajectory_without_target_lipid.xtc \
    --stop 10 \
    --debug
```

## Requirements

- Python 3.10+
- MDAnalysis
- NumPy
- Pandas
- Matplotlib
- Seaborn
- tqdm

See `requirements.txt` for specific versions.

## Output Structure

```
lipac_results/
├── with_target_lipid/
│   ├── contact_complementarity.csv
│   ├── protein_contact_maps.png
│   └── lipid_contact_heatmaps.png
├── without_target_lipid/
│   ├── contact_complementarity.csv
│   ├── protein_contact_maps.png
│   └── lipid_contact_heatmaps.png
├── comparison/
│   └── contact_comparison_plots.png
└── leaflet_info.pickle
```