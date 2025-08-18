"""Configuration parameters for lipid-protein contact analysis"""

# Default file paths - WITH target lipid system (GM3 prevents dimerization)
DEFAULT_WITH_LIPID_PSF = 'test_data/with_target_lipid/system.psf'
DEFAULT_WITH_LIPID_XTC = 'test_data/with_target_lipid/trajectory.xtc'

# Default file paths - WITHOUT target lipid system (allows dimerization)
DEFAULT_WITHOUT_LIPID_PSF = 'test_data/without_target_lipid/system.psf'
DEFAULT_WITHOUT_LIPID_XTC = 'test_data/without_target_lipid/trajectory.xtc'

# Base output directory
BASE_OUTPUT_DIR = 'output/lipac_results'

# Output directories
DEFAULT_WITH_LIPID_OUTPUT = f'{BASE_OUTPUT_DIR}/with_target_lipid'
DEFAULT_WITHOUT_LIPID_OUTPUT = f'{BASE_OUTPUT_DIR}/without_target_lipid'
DEFAULT_COMPARISON_OUTPUT = f'{BASE_OUTPUT_DIR}/comparison'

# Temporary files directory (for leaflet info, etc.)
TEMP_FILES_DIR = f'{BASE_OUTPUT_DIR}/temp_files'

# Frame processing parameters
START_FRAME = 0
STOP_FRAME = 200
STEP_FRAME = 1

# Contact detection parameters
CONTACT_CUTOFF = 6.0  # Angstrom - cutoff for lipid-protein contacts (3D)
PROTEIN_CONTACT_CUTOFF = 6.0  # Angstrom - cutoff for protein-protein contacts
DIMER_CUTOFF = 20.0  # Angstrom - cutoff distance for protein pairs to be considered as dimer

# Parallel processing parameters
BATCH_SIZE = 50  # Batch size for parallel processing
MIN_CORES = 2  # Minimum number of cores to use (>1 to force parallel processing)

# Leaflet detection parameters
LEAFLET_CUTOFF = 10.0  # Cutoff value for leaflet detection
RESIDUE_OFFSET = 556  # Residue numbering offset for conversion

# Target lipid for special analysis (can be any lipid type)
TARGET_LIPID = 'DPG3'  # Default target lipid (GM3), can be changed to any lipid type

# All lipid types in the system
LIPID_TYPES = ['CHOL', 'DPSM', 'DIPC', 'DPG3', 'DOPS']