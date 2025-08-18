"""
Configuration parameters for Bayesian analysis
"""

# Default file paths (from stage1 config)
DEFAULT_WITH_LIPID_DATA = 'output/lipac_results/with_target_lipid/contact_complementarity.csv'
DEFAULT_WITHOUT_LIPID_DATA = 'output/lipac_results/without_target_lipid/contact_complementarity.csv'
DEFAULT_OUTPUT_DIR = 'output/bayesian_analysis'
DEFAULT_CAUSAL_DATA_DIR = 'output/lipac_results/with_target_lipid'

# MCMC parameters
MCMC_SAMPLES = 2000   # MCMC samples default 2000
TUNE_SAMPLES = 1000   # Tuning samples default 1000
CHAINS = 4            # Number of chains
RANDOM_SEED = 42      # Random seed
TARGET_ACCEPT = 0.95  # Target acceptance rate

# Target lipid (imported from stage1 config if available)
TARGET_LIPID = 'DPG3'        # Default target lipid, should match stage1 config