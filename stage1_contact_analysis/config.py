"""Configuration parameters for lipid-protein contact analysis"""

# Default file paths - WITH target lipid system (GM3 prevents dimerization)
<<<<<<< HEAD
DEFAULT_WITH_LIPID_PSF = '../../test/test_system_with_mediator.psf'
DEFAULT_WITH_LIPID_XTC = '../../test/test_trajectory_with_mediator.xtc'

# Default file paths - WITHOUT target lipid system (allows dimerization)
DEFAULT_WITHOUT_LIPID_PSF = '../../test/test_system_without_mediator.psf'
DEFAULT_WITHOUT_LIPID_XTC = '../../test/test_trajectory_without_mediator.xtc'
=======
DEFAULT_WITH_LIPID_PSF = 'test_data/with_target_lipid/system.psf'
DEFAULT_WITH_LIPID_XTC = 'test_data/with_target_lipid/trajectory.xtc'

# Default file paths - WITHOUT target lipid system (allows dimerization)
DEFAULT_WITHOUT_LIPID_PSF = 'test_data/without_target_lipid/system.psf'
DEFAULT_WITHOUT_LIPID_XTC = 'test_data/without_target_lipid/trajectory.xtc'
>>>>>>> 4880c9ee1958ba009c2a2f29d635256470a257d5

# Base output directory
BASE_OUTPUT_DIR = 'output/lipac_results'

# Output directories
DEFAULT_WITH_LIPID_OUTPUT = f'{BASE_OUTPUT_DIR}/with_target_lipid'
DEFAULT_WITHOUT_LIPID_OUTPUT = f'{BASE_OUTPUT_DIR}/without_target_lipid'
DEFAULT_COMPARISON_OUTPUT = f'{BASE_OUTPUT_DIR}/comparison'

# Temporary files directory (for leaflet info, etc.)
TEMP_FILES_DIR = f'{BASE_OUTPUT_DIR}/temp_files'

<<<<<<< HEAD
# Frame processing parameters (adjusted for test trajectory)
START_FRAME = 0
STOP_FRAME = 199
=======
# Frame processing parameters
START_FRAME = 0
STOP_FRAME = 200
>>>>>>> 4880c9ee1958ba009c2a2f29d635256470a257d5
STEP_FRAME = 1

# Contact detection parameters
CONTACT_CUTOFF = 6.0  # Angstrom - cutoff for lipid-protein contacts (3D)
PROTEIN_CONTACT_CUTOFF = 6.0  # Angstrom - cutoff for protein-protein contacts
<<<<<<< HEAD
DIMER_CUTOFF = 20.0  # Angstrom - cutoff distance for protein pairs to be considered as dimer (same as mirage)
=======
DIMER_CUTOFF = 20.0  # Angstrom - cutoff distance for protein pairs to be considered as dimer
>>>>>>> 4880c9ee1958ba009c2a2f29d635256470a257d5

# Parallel processing parameters
BATCH_SIZE = 50  # Batch size for parallel processing
MIN_CORES = 2  # Minimum number of cores to use (>1 to force parallel processing)

# Leaflet detection parameters
LEAFLET_CUTOFF = 10.0  # Cutoff value for leaflet detection
<<<<<<< HEAD
RESIDUE_OFFSET = 556  # Residue conversion offset (same as premirage4.py)
=======
RESIDUE_OFFSET = 556  # Residue numbering offset for conversion
>>>>>>> 4880c9ee1958ba009c2a2f29d635256470a257d5

# Target lipid for special analysis (can be any lipid type)
TARGET_LIPID = 'DPG3'  # Default target lipid (GM3), can be changed to any lipid type

# All lipid types in the system
LIPID_TYPES = ['CHOL', 'DPSM', 'DIPC', 'DPG3', 'DOPS']