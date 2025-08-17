"""Molecular placement utilities."""

import numpy as np
from ase import Atoms


def get_drop_snapped(atoms, surface, height=3.0):
    """
    Place atoms above a surface at a specified height.
    
    Args:
        atoms: ASE Atoms object to place
        surface: Surface atoms object
        height: Height above surface in Angstroms
        
    Returns:
        Atoms: Positioned atoms object
    """
    atoms_copy = atoms.copy()
    
    # Get surface height (assuming surface is in xy plane)
    if hasattr(surface, 'positions'):
        surface_height = np.max(surface.positions[:, 2])
    else:
        surface_height = 0.0
    
    # Move atoms to specified height above surface
    atoms_copy.positions[:, 2] += surface_height + height
    
    return atoms_copy
