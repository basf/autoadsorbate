"""Top-level package for autoadsorbate."""

__author__ = """Fakoe Edvin"""
__email__ = "edvin.fako@basf.com"
__version__ = "0.2.0"

# Core classes
from .core import Fragment, Surface, Intermediate

# Molecular operations
from .molecular import get_marked_smiles

# String utilities
from .string_utils import _example_config, construct_smiles

# Utilities
from .utils import get_drop_snapped, compute_energy

# Visualization
from .viz import docs_plot_conformers, docs_plot_sites

__all__ = [
    # Core classes
    "Fragment",
    "Surface", 
    "Intermediate",
    
    # Molecular operations
    "get_marked_smiles",
    
    # String utilities
    "construct_smiles",
    "_example_config",
    
    # Utilities
    "get_drop_snapped",
    "compute_energy",
    
    # Visualization
    "docs_plot_conformers",
    "docs_plot_sites",
]
