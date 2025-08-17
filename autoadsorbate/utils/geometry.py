"""Geometry utility functions."""

import numpy as np
from ase import Atoms


def rotation_matrix_from_vectors(vec1, vec2):
    """
    Find the rotation matrix that aligns vec1 to vec2.
    
    Args:
        vec1: A 3d "source" vector
        vec2: A 3d "destination" vector
        
    Returns:
        mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """
    a, b = (
        (vec1 / np.linalg.norm(vec1)).reshape(3),
        (vec2 / np.linalg.norm(vec2)).reshape(3),
    )
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    if all(a == b):
        rotation_matrix = np.eye(3)
    else:
        kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s**2))
    return rotation_matrix


def random_three_vector():
    """
    Generates a random 3D unit vector.

    Returns:
        tuple: A tuple containing the x, y, and z components of the random 3D unit vector.
    """
    phi = np.random.uniform(0, np.pi * 2)
    costheta = np.random.uniform(-1, 1)

    theta = np.arccos(costheta)
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)
    return (x, y, z)


def random_rotate(atoms):
    """
    Applies a random rotation to the given atoms object.

    Args:
        atoms (object): An object containing an array of atoms to be rotated.

    Returns:
        object: The atoms object after applying the random rotation.
    """
    axis = random_three_vector()
    angle = np.random.uniform(0, 2 * np.pi)
    atoms.rotate(angle, axis)
    return atoms
