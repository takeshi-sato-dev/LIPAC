"""Stage 1: Contact Analysis from MD Trajectories

This module provides comprehensive analysis of lipid-protein and protein-protein
interactions from molecular dynamics simulations.
"""

__version__ = "1.0.0"

from .main import run_analysis, main
from .config import *

__all__ = ['run_analysis', 'main']