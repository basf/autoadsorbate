"""Tests for utility functions."""

import pytest
import numpy as np
from ase import Atoms

from autoadsorbate.utils import (
    rotation_matrix_from_vectors,
    random_three_vector,
    random_rotate,
    get_sorted_by_snap_dist
)


class TestGeometry:
    """Test cases for geometry utilities."""
    
    def test_rotation_matrix_from_vectors_identity(self):
        """Test rotation matrix for identical vectors."""
        vec1 = np.array([0, 0, 1])
        vec2 = np.array([0, 0, 1])
        
        matrix = rotation_matrix_from_vectors(vec1, vec2)
        assert np.allclose(matrix, np.eye(3))
    
    def test_rotation_matrix_from_vectors_orthogonal(self):
        """Test rotation matrix for orthogonal vectors."""
        vec1 = np.array([1, 0, 0])
        vec2 = np.array([0, 1, 0])
        
        matrix = rotation_matrix_from_vectors(vec1, vec2)
        rotated = matrix @ vec1
        
        # Should align with vec2 (allowing for sign differences)
        assert np.allclose(np.abs(rotated), np.abs(vec2))
    
    def test_rotation_matrix_from_vectors_normalized(self):
        """Test that rotation matrix works with non-normalized vectors."""
        vec1 = np.array([2, 0, 0])
        vec2 = np.array([0, 3, 0])
        
        matrix = rotation_matrix_from_vectors(vec1, vec2)
        rotated = matrix @ vec1
        
        # Should align with vec2 direction
        assert np.allclose(rotated / np.linalg.norm(rotated), 
                          vec2 / np.linalg.norm(vec2))
    
    def test_random_three_vector(self):
        """Test random 3D vector generation."""
        vector = random_three_vector()
        
        assert len(vector) == 3
        assert isinstance(vector, tuple)
        assert np.isclose(np.linalg.norm(vector), 1.0)
    
    def test_random_three_vector_multiple(self):
        """Test multiple random vectors are different."""
        vectors = [random_three_vector() for _ in range(10)]
        vector_arrays = [np.array(v) for v in vectors]
        
        # Check that at least some vectors are different
        different_pairs = 0
        for i in range(len(vector_arrays)):
            for j in range(i+1, len(vector_arrays)):
                if not np.allclose(vector_arrays[i], vector_arrays[j]):
                    different_pairs += 1
        
        assert different_pairs > 0
    
    def test_random_rotate(self):
        """Test random rotation of atoms."""
        atoms = Atoms("CO", positions=[[0, 0, 0], [1, 0, 0]])
        original_positions = atoms.positions.copy()
        
        rotated_atoms = random_rotate(atoms)
        
        assert isinstance(rotated_atoms, Atoms)
        assert len(rotated_atoms) == len(atoms)
        # Positions should be different after rotation
        assert not np.allclose(rotated_atoms.positions, original_positions)


class TestSorting:
    """Test cases for sorting utilities."""
    
    def test_get_sorted_by_snap_dist_single(self):
        """Test sorting with single conformer."""
        atoms = Atoms("CO", positions=[[0, 0, 0], [1, 0, 0]])
        conformers = [atoms]
        
        sorted_conformers = get_sorted_by_snap_dist(conformers)
        assert len(sorted_conformers) == 1
        assert sorted_conformers[0] is conformers[0]
    
    def test_get_sorted_by_snap_dist_multiple(self):
        """Test sorting with multiple conformers."""
        atoms1 = Atoms("CO", positions=[[0, 0, 0], [1, 0, 0]])
        atoms2 = Atoms("CO", positions=[[0, 0, 0], [0, 1, 0]])
        atoms3 = Atoms("CO", positions=[[0, 0, 0], [0, 0, 1]])
        
        conformers = [atoms1, atoms2, atoms3]
        sorted_conformers = get_sorted_by_snap_dist(conformers)
        
        assert len(sorted_conformers) == 3
        assert all(isinstance(conf, Atoms) for conf in sorted_conformers)
    
    def test_get_sorted_by_snap_dist_empty(self):
        """Test sorting with empty list."""
        sorted_conformers = get_sorted_by_snap_dist([])
        assert sorted_conformers == []
    
    def test_calculate_snap_distance(self):
        """Test snap distance calculation."""
        from autoadsorbate.utils.sorting import calculate_snap_distance
        
        atoms1 = Atoms("CO", positions=[[0, 0, 0], [1, 0, 0]])
        atoms2 = Atoms("CO", positions=[[0, 0, 0], [0, 1, 0]])
        
        distance = calculate_snap_distance(atoms1, atoms2)
        assert isinstance(distance, (int, float))
        assert distance >= 0
    
    def test_get_connectivity_matrix(self):
        """Test connectivity matrix generation."""
        from autoadsorbate.utils.sorting import get_connectivity_matrix
        
        atoms = Atoms("CO", positions=[[0, 0, 0], [1, 0, 0]])
        connectivity = get_connectivity_matrix(atoms)
        
        assert isinstance(connectivity, np.ndarray)
        assert connectivity.shape == (2, 2)
        assert np.allclose(connectivity.diagonal(), 0)  # No self-connections


class TestIntegration:
    """Integration tests for utility functions."""
    
    def test_rotation_and_sorting_integration(self):
        """Test integration of rotation and sorting."""
        # Create conformers
        atoms1 = Atoms("CO", positions=[[0, 0, 0], [1, 0, 0]])
        atoms2 = Atoms("CO", positions=[[0, 0, 0], [0, 1, 0]])
        conformers = [atoms1, atoms2]
        
        # Rotate one conformer
        rotated_conformers = [random_rotate(conf.copy()) for conf in conformers]
        
        # Sort rotated conformers
        sorted_conformers = get_sorted_by_snap_dist(rotated_conformers)
        
        assert len(sorted_conformers) == 2
        assert all(isinstance(conf, Atoms) for conf in sorted_conformers)
    
    def test_geometry_operations_chain(self):
        """Test chaining of geometry operations."""
        vec1 = np.array([1, 0, 0])
        vec2 = np.array([0, 1, 0])
        
        # Create rotation matrix
        matrix = rotation_matrix_from_vectors(vec1, vec2)
        
        # Apply to atoms
        atoms = Atoms("CO", positions=[[0, 0, 0], [1, 0, 0]])
        atoms.positions = atoms.positions @ matrix.T
        
        # Randomly rotate
        rotated_atoms = random_rotate(atoms)
        
        assert isinstance(rotated_atoms, Atoms)
        assert len(rotated_atoms) == 2
