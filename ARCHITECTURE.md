# LIPAC Architecture Overview

## Package Structure

```
lipac/
├── stage1_contact_analysis/          # Stage 1: Contact detection and quantification
│   ├── core/
│   │   ├── trajectory_loader.py      # MD trajectory loading and selection
│   │   ├── frame_processor.py        # Frame-by-frame processing
│   │   └── contact_calculator.py     # Contact calculation algorithms
│   ├── analysis/
│   │   ├── complementarity.py        # Competition analysis
│   │   ├── residue_contacts.py       # Residue-level analysis
│   │   └── comparison.py             # System comparison
│   ├── utils/
│   │   ├── checkpoint.py             # Checkpoint/restart functionality
│   │   └── parallel.py               # Parallel processing
│   ├── visualization/
│   │   └── plotting.py               # Contact visualization
│   ├── config.py                     # Configuration parameters
│   └── main.py                       # Stage 1 entry point
│
├── stage2_contact_analysis/          # Stage 2: Bayesian statistical analysis
│   ├── core/
│   │   └── data_loader.py            # Data loading and preprocessing
│   ├── analysis/
│   │   ├── exploratory.py            # Exploratory data analysis
│   │   ├── bayesian_models.py        # Core Bayesian models
│   │   ├── hierarchical_models.py    # Hierarchical Bayesian analysis
│   │   ├── causal_analysis.py        # Causal inference methods
│   │   ├── temporal_clustering.py    # Phase transition detection
│   │   ├── domain_formation.py       # Nanodomain analysis
│   │   ├── lipid_analysis.py         # Lipid-specific effects
│   │   └── residue_analysis.py       # Residue-level statistics
│   ├── utils/
│   │   └── reports.py                # Report generation
│   ├── visualization/
│   │   └── plots.py                  # Statistical visualization
│   ├── config.py                     # Configuration parameters
│   └── main.py                       # Stage 2 entry point
│
├── test_data/                        # Test datasets
├── examples/                         # Example scripts and data
├── docs/                             # Documentation
└── paper/                            # JOSS paper
```

## Analysis Pipeline

### Stage 1: Contact Detection and Quantification

1. **Trajectory Loading** (`trajectory_loader.py`)
   - Load MD trajectories using MDAnalysis
   - Protein and lipid selection
   - Leaflet identification

2. **Frame Processing** (`frame_processor.py`)
   - Parallel frame-by-frame analysis
   - Contact calculation with cutoff distances
   - Target lipid binding state tracking

3. **Contact Calculation** (`contact_calculator.py`)
   - Protein-protein contact matrices
   - Lipid-protein contact counting
   - Distance-based analysis

4. **Data Aggregation** (`analysis/`)
   - Contact complementarity analysis
   - Residue-level aggregation
   - Competition analysis

### Stage 2: Bayesian Statistical Analysis

1. **Data Loading** (`data_loader.py`)
   - Load Stage 1 results
   - Data preprocessing and matching
   - Format conversion

2. **Exploratory Analysis** (`exploratory.py`)
   - Descriptive statistics
   - Distribution visualization
   - Correlation analysis

3. **Bayesian Modeling** (`bayesian_models.py`, `hierarchical_models.py`)
   - Hierarchical regression models
   - Uncertainty quantification
   - Effect size estimation

4. **Causal Analysis** (`causal_analysis.py`)
   - Target lipid binding effects
   - Compositional changes
   - Causal inference

5. **Temporal Analysis** (`temporal_clustering.py`)
   - Phase transition detection
   - Equilibrium/non-equilibrium identification
   - Dynamic clustering

6. **Domain Analysis** (`domain_formation.py`)
   - Local enrichment calculation
   - Spatial heterogeneity analysis
   - Nanodomain detection

## Key Algorithms

### Contact Detection Algorithm
```python
def calculate_contacts(protein, lipids, cutoff=6.0):
    """
    Calculate contacts using distance matrix approach
    - O(N*M) complexity where N=protein atoms, M=lipid atoms
    - Periodic boundary condition handling
    - Efficient memory usage through chunking
    """
```

### Bayesian Changepoint Detection
```python
def detect_phase_transition(flux_data):
    """
    Bayesian changepoint detection for phase transitions
    - PyMC implementation
    - Automatic threshold detection
    - Uncertainty quantification
    """
```

### Hierarchical Bayesian Model
```python
def hierarchical_regression(data):
    """
    Multi-level model accounting for:
    - Global effects (all proteins)
    - Protein-specific effects
    - Lipid-specific effects
    - Interaction effects
    """
```

## Performance Optimizations

1. **Parallel Processing**
   - Frame-level parallelization
   - Batch processing with checkpoints
   - Memory-efficient chunking

2. **Caching and Checkpoints**
   - Intermediate result caching
   - Automatic restart capability
   - Progress tracking

3. **Memory Management**
   - Streaming data processing
   - Efficient data structures
   - Garbage collection optimization

## Dependencies and Requirements

### Core Dependencies
- **MDAnalysis**: Trajectory analysis
- **PyMC**: Bayesian statistical modeling
- **NumPy/SciPy**: Numerical computing
- **Pandas**: Data manipulation
- **ArviZ**: Bayesian diagnostics

### Visualization
- **Matplotlib**: Core plotting
- **Seaborn**: Statistical visualization
- **Plotly**: Interactive plots

### Performance
- **Multiprocessing**: Parallel processing
- **NumPy**: Vectorized operations
- **Cython** (optional): Performance critical sections

## Design Principles

1. **Modularity**: Clear separation of concerns between stages
2. **Reproducibility**: Deterministic results with seed control
3. **Scalability**: Efficient handling of large trajectories
4. **Extensibility**: Easy addition of new analysis methods
5. **Robustness**: Error handling and data validation

## Data Flow

```
MD Trajectories → Stage 1 → Contact Data → Stage 2 → Statistical Results
     ↓              ↓            ↓            ↓
  Topology      Parallel     Causal       Bayesian
    Files      Processing     Data         Analysis
     ↓              ↓            ↓            ↓
  Selection    Checkpoints   Competition   Reports
     ↓              ↓            ↓            ↓
   Contacts     Aggregation   Temporal     Plots
```

## Configuration Management

- Centralized configuration files for each stage
- Command-line interface for parameter override
- Environment-specific settings
- Default parameter validation

This architecture provides a robust, scalable framework for comprehensive lipid-protein interaction analysis with state-of-the-art statistical methods.