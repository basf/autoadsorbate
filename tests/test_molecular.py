"""Tests for molecular operations."""

import pytest
import numpy as np
from ase import Atoms

from autoadsorbate.molecular import (
    conformers_from_smile,
    _reset_position,
    _reset_rotation,
    get_marked_smiles,
    check_smile,
    remove_canonical_duplicates,
    sort_smiles
)


class TestConformers:
    """Test cases for conformer generation."""
    
    def test_conformers_from_smile_basic(self):
        """Test basic conformer generation."""
        conformers = conformers_from_smile("COC", conformer_count=3)
        assert len(conformers) == 3
        assert all(isinstance(conf, Atoms) for conf in conformers)
        assert all(len(conf) > 0 for conf in conformers)
    
    def test_conformers_from_smile_custom_seed(self):
        """Test conformer generation with custom seed."""
        conformers1 = conformers_from_smile("CCO", conformer_count=2, random_seed=42)
        conformers2 = conformers_from_smile("CCO", conformer_count=2, random_seed=42)
        
        # With same seed, should get same conformers
        assert len(conformers1) == len(conformers2)
    
    def test_conformers_from_smile_invalid_smiles(self):
        """Test conformer generation with invalid SMILES."""
        with pytest.raises(Exception):
            conformers_from_smile("invalid_smiles", conformer_count=1)


class TestAlignment:
    """Test cases for molecular alignment."""
    
    def test_reset_position(self):
        """Test position resetting."""
        # Create a simple molecule
        atoms = Atoms("CO", positions=[[0, 0, 0], [1, 1, 1]])
        original_com = atoms.get_center_of_mass()
        
        reset_atoms = _reset_position(atoms)
        new_com = reset_atoms.get_center_of_mass()
        
        # Center of mass should be at origin
        assert np.allclose(new_com, [0, 0, 0], atol=1e-10)
        assert not np.allclose(original_com, [0, 0, 0])
    
    def test_reset_rotation_basic(self):
        """Test basic rotation resetting."""
        atoms = Atoms("CO", positions=[[0, 0, 0], [1, 0, 0]])
        rotated_atoms = _reset_rotation(atoms)
        
        assert isinstance(rotated_atoms, Atoms)
        assert len(rotated_atoms) == len(atoms)
    
    def test_reset_rotation_with_custom_vector(self):
        """Test rotation resetting with custom normal vector."""
        atoms = Atoms("CO", positions=[[0, 0, 0], [1, 0, 0]])
        n_vector = np.array([0, 1, 0])
        
        rotated_atoms = _reset_rotation(atoms, n_vector=n_vector)
        assert isinstance(rotated_atoms, Atoms)
    
    def test_reset_rotation_with_angle(self):
        """Test rotation resetting with custom angle."""
        atoms = Atoms("CO", positions=[[0, 0, 0], [1, 0, 0]])
        
        rotated_atoms = _reset_rotation(atoms, rot_deg=45)
        assert isinstance(rotated_atoms, Atoms)


class TestSmilesProcessing:
    """Test cases for SMILES processing."""
    
    def test_get_marked_smiles(self):
        """Test SMILES marking."""
        marked = get_marked_smiles("COC")
        assert marked.startswith("*")
        assert "COC" in marked
    
    def test_check_smile_valid(self):
        """Test SMILES validation with valid SMILES."""
        assert check_smile("COC") == True
        assert check_smile("CCO") == True
        assert check_smile("c1ccccc1") == True  # Benzene
    
    def test_check_smile_invalid(self):
        """Test SMILES validation with invalid SMILES."""
        assert check_smile("invalid") == False
        assert check_smile("") == False
    
    def test_remove_canonical_duplicates(self):
        """Test removal of canonical duplicates."""
        smiles_list = ["COC", "CCO", "COC", "OCC"]  # COC and CCO are the same molecule
        unique_smiles = remove_canonical_duplicates(smiles_list)
        
        assert len(unique_smiles) < len(smiles_list)
        assert all(smiles in smiles_list for smiles in unique_smiles)
    
    def test_sort_smiles(self):
        """Test SMILES sorting by type."""
        smiles_list = ["COC", "*CCO", "CCO", "*COC"]
        top_smiles, nontop_smiles = sort_smiles(smiles_list)
        
        assert all(smiles.startswith("*") for smiles in top_smiles)
        assert all(not smiles.startswith("*") for smiles in nontop_smiles)
        assert len(top_smiles) + len(nontop_smiles) == len(smiles_list)
    
    def test_sort_smiles_empty(self):
        """Test SMILES sorting with empty list."""
        top_smiles, nontop_smiles = sort_smiles([])
        assert top_smiles == []
        assert nontop_smiles == []
