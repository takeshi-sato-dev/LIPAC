"""
BayesianLipidAnalysis - Bayesian analysis of lipid-protein interactions
Analysis of contact complementarity data to understand the relationship between lipid binding and protein dimer formation

Author: Takeshi Sato, PhD
Kyoto Pharmaceutical University
2024
"""

__version__ = "1.0.0"
__author__ = "Takeshi Sato"
__email__ = "your.email@example.com"  # Update with actual email

from . import core
from . import analysis
from . import visualization
from . import utils

__all__ = [
    'core',
    'analysis',
    'visualization',
    'utils'
]
