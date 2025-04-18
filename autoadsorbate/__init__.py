"""Top-level package for autoadsorbate."""

__author__ = """Fakoe Edvin"""
__email__ = "edvin.fako@basf.com"
__version__ = "0.2.0"

from autoadsorbate.autoadsorbate import Fragment, Surface
from autoadsorbate.Smile import get_marked_smiles
from autoadsorbate.string_utils import _example_config, construct_smiles
from autoadsorbate.utils import docs_plot_conformers, docs_plot_sites, get_drop_snapped, compute_energy

__all__ = [
    "Fragment",
    "Surface",
    "docs_plot_conformers",
    "get_marked_smiles",
    "docs_plot_sites",
    "construct_smiles",
    "_example_config",
    "get_drop_snapped",
    "compute_energy",
]
