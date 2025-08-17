"""Analysis utility functions."""

import numpy as np
from ase import Atoms


def count_C_next_to_O(atoms, cutoff=1.6):
    """
    Count carbon atoms next to oxygen atoms.
    
    Args:
        atoms: ASE Atoms object
        cutoff: Distance cutoff for neighbor counting
        
    Returns:
        int: Number of C-O pairs within cutoff distance
    """
    count = 0
    positions = atoms.positions
    symbols = atoms.get_chemical_symbols()
    
    for i, symbol_i in enumerate(symbols):
        if symbol_i == 'C':
            for j, symbol_j in enumerate(symbols):
                if symbol_j == 'O' and i != j:
                    distance = np.linalg.norm(positions[i] - positions[j])
                    if distance <= cutoff:
                        count += 1
    
    return count
