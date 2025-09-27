"""
Configuration parameters for Bayesian analysis
"""

# Default file paths (from stage1 config)
<<<<<<< HEAD
DEFAULT_WITH_LIPID_DATA = '../stage1_contact_analysis/output/lipac_results/with_target_lipid/contact_complementarity.csv'
DEFAULT_WITHOUT_LIPID_DATA = '../stage1_contact_analysis/output/lipac_results/without_target_lipid/contact_complementarity.csv'
DEFAULT_OUTPUT_DIR = '../stage1_contact_analysis/output/bayesian_analysis'
DEFAULT_CAUSAL_DATA_DIR = '../stage1_contact_analysis/output/lipac_results/with_target_lipid'
=======
DEFAULT_WITH_LIPID_DATA = 'your/path/to/with_target_lipid/contact_complementarity.csv'
DEFAULT_WITHOUT_LIPID_DATA = 'your/path/to/without_target_lipid/contact_complementarity.csv'
DEFAULT_OUTPUT_DIR = 'your/path/to/bayesian_analysis'
DEFAULT_CAUSAL_DATA_DIR = 'your/path/to/with_target_lipid'
>>>>>>> 4880c9ee1958ba009c2a2f29d635256470a257d5

# MCMC parameters
MCMC_SAMPLES = 2000   # MCMC samples default 2000
TUNE_SAMPLES = 1000   # Tuning samples default 1000
CHAINS = 4            # Number of chains
RANDOM_SEED = 42      # Random seed
TARGET_ACCEPT = 0.95  # Target acceptance rate

# Target lipid (imported from stage1 config if available)
TARGET_LIPID = 'DPG3'        # Default target lipid, should match stage1 config