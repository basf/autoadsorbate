"""Molecular alignment and positioning functions."""

import numpy as np
from ase import Atoms
from ..utils.geometry import rotation_matrix_from_vectors


def _reset_position(atoms: Atoms) -> Atoms:
    """
    Resets the position of atoms by centering them at the origin.

    Args:
        atoms (Atoms): The atoms object to reposition.

    Returns:
        Atoms: The repositioned atoms object.
    """
    atoms_copy = atoms.copy()
    center_of_mass = atoms_copy.get_center_of_mass()
    atoms_copy.positions -= center_of_mass
    return atoms_copy


def _reset_rotation(
    atoms: Atoms, 
    n_vector: np.ndarray = np.array([0, 0, 1]), 
    rot_deg: float = 0
) -> Atoms:
    """
    Resets the rotation of atoms by aligning them with a specified normal vector.

    Args:
        atoms (Atoms): The atoms object to rotate.
        n_vector (np.ndarray, optional): The normal vector for alignment. Defaults to [0, 0, 1].
        rot_deg (float, optional): Additional rotation angle in degrees. Defaults to 0.

    Returns:
        Atoms: The rotated atoms object.
    """
    atoms_copy = atoms.copy()
    
    # Get the current normal vector (assuming it's along the z-axis)
    current_normal = np.array([0, 0, 1])
    
    # Calculate rotation matrix to align with the target normal vector
    rotation_matrix = rotation_matrix_from_vectors(current_normal, n_vector)
    
    # Apply rotation
    atoms_copy.positions = np.dot(atoms_copy.positions, rotation_matrix.T)
    
    # Apply additional rotation if specified
    if rot_deg != 0:
        rot_rad = np.radians(rot_deg)
        cos_rot = np.cos(rot_rad)
        sin_rot = np.sin(rot_rad)
        
        # Rotation matrix around the normal vector
        rot_matrix = np.array([
            [cos_rot, -sin_rot, 0],
            [sin_rot, cos_rot, 0],
            [0, 0, 1]
        ])
        
        atoms_copy.positions = np.dot(atoms_copy.positions, rot_matrix.T)
    
    return atoms_copy
