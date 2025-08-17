"""Energy calculation utilities."""

import numpy as np


def compute_energy(atoms, calculator=None):
    """
    Compute energy for a given atomic structure.
    
    Args:
        atoms: ASE Atoms object
        calculator: ASE calculator object (optional)
        
    Returns:
        float: Computed energy
    """
    if calculator is not None:
        atoms.calc = calculator
        return atoms.get_potential_energy()
    else:
        # Return a placeholder energy if no calculator is provided
        return 0.0
