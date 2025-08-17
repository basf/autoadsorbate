"""Sorting utility functions."""

import numpy as np
from ase import Atoms


def get_sorted_by_snap_dist(conformers):
    """
    Sorts conformers by their snap distance to ensure diverse orientations.

    Args:
        conformers (list): List of ASE Atoms objects representing conformers.

    Returns:
        list: Sorted list of conformers.
    """
    if len(conformers) <= 1:
        return conformers
    
    # Calculate snap distances between all pairs of conformers
    snap_distances = []
    for i, conf1 in enumerate(conformers):
        for j, conf2 in enumerate(conformers[i+1:], i+1):
            dist = calculate_snap_distance(conf1, conf2)
            snap_distances.append((dist, i, j))
    
    # Sort by distance (largest first to maximize diversity)
    snap_distances.sort(reverse=True)
    
    # Use a greedy algorithm to select diverse conformers
    selected = [0]  # Start with the first conformer
    available = set(range(1, len(conformers)))
    
    while available and len(selected) < len(conformers):
        best_score = -1
        best_conformer = None
        
        for conf_idx in available:
            # Calculate minimum distance to already selected conformers
            min_dist = float('inf')
            for selected_idx in selected:
                for dist, i, j in snap_distances:
                    if (i == conf_idx and j == selected_idx) or (i == selected_idx and j == conf_idx):
                        min_dist = min(min_dist, dist)
                        break
            
            if min_dist > best_score:
                best_score = min_dist
                best_conformer = conf_idx
        
        if best_conformer is not None:
            selected.append(best_conformer)
            available.remove(best_conformer)
        else:
            break
    
    return [conformers[i] for i in selected]


def calculate_snap_distance(atoms1, atoms2, bond_cutoff=1.6):
    """
    Calculates the snap distance between two atomic structures.

    Args:
        atoms1 (Atoms): First atomic structure.
        atoms2 (Atoms): Second atomic structure.
        bond_cutoff (float): Distance cutoff for bond detection.

    Returns:
        float: Snap distance between the two structures.
    """
    # Get connectivity matrices
    conn1 = get_connectivity_matrix(atoms1, bond_cutoff)
    conn2 = get_connectivity_matrix(atoms2, bond_cutoff)
    
    # Calculate difference
    diff = np.abs(conn1 - conn2)
    return np.sum(diff)


def get_connectivity_matrix(atoms, bond_cutoff=1.6):
    """
    Creates a connectivity matrix for an atomic structure.

    Args:
        atoms (Atoms): Atomic structure.
        bond_cutoff (float): Distance cutoff for bond detection.

    Returns:
        np.ndarray: Connectivity matrix.
    """
    distances = atoms.get_all_distances()
    connectivity = (distances < bond_cutoff).astype(int)
    # Remove self-connections
    np.fill_diagonal(connectivity, 0)
    return connectivity
