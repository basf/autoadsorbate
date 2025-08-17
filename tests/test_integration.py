"""Integration tests for autoadsorbate."""

import pytest
import numpy as np
from ase import Atoms
from ase.build import fcc111

from autoadsorbate import Fragment, Surface, Intermediate
from autoadsorbate.molecular import conformers_from_smile, get_marked_smiles
from autoadsorbate.utils import rotation_matrix_from_vectors, get_sorted_by_snap_dist
from autoadsorbate.viz import docs_plot_conformers


class TestFragmentSurfaceIntegration:
    """Integration tests for Fragment and Surface classes."""
    
    def test_fragment_surface_basic_workflow(self):
        """Test basic workflow with Fragment and Surface."""
        # Create fragment
        fragment = Fragment(smile="COC", to_initialize=3)
        assert len(fragment.conformers) == 3
        
        # Create surface
        slab = fcc111("Cu", (2, 2, 2), periodic=True, vacuum=10)
        surface = Surface(slab)
        assert isinstance(surface.site_dict, dict)
        
        # Get a conformer
        conformer = fragment.get_conformer(0)
        assert isinstance(conformer, Atoms)
        assert len(conformer) > 0
    
    def test_fragment_with_surface_normal(self):
        """Test Fragment with surface normal alignment."""
        # Create surface
        slab = fcc111("Cu", (2, 2, 2), periodic=True, vacuum=10)
        surface = Surface(slab)
        
        # Get surface normal
        normal = surface.get_surface_normal()
        assert np.isclose(np.linalg.norm(normal), 1.0)
        
        # Create fragment aligned with surface normal
        fragment = Fragment(smile="COC", to_initialize=2)
        conformer = fragment.get_conformer(0, n_vector=normal)
        assert isinstance(conformer, Atoms)
    
    def test_multiple_fragments_on_surface(self):
        """Test multiple fragments on the same surface."""
        slab = fcc111("Cu", (2, 2, 2), periodic=True, vacuum=10)
        surface = Surface(slab)
        
        fragments = [
            Fragment(smile="COC", to_initialize=2),
            Fragment(smile="CCO", to_initialize=2)
        ]
        
        for fragment in fragments:
            conformer = fragment.get_conformer(0)
            assert isinstance(conformer, Atoms)
            assert len(conformer) > 0


class TestMolecularOperationsIntegration:
    """Integration tests for molecular operations."""
    
    def test_conformer_generation_and_sorting(self):
        """Test conformer generation and sorting integration."""
        # Generate conformers
        conformers = conformers_from_smile("COC", conformer_count=5)
        assert len(conformers) == 5
        
        # Sort conformers
        sorted_conformers = get_sorted_by_snap_dist(conformers)
        assert len(sorted_conformers) == 5
        
        # All should be Atoms objects
        assert all(isinstance(conf, Atoms) for conf in sorted_conformers)
    
    def test_smiles_processing_workflow(self):
        """Test complete SMILES processing workflow."""
        # Test SMILES marking
        original_smiles = "COC"
        marked_smiles = get_marked_smiles(original_smiles)
        assert marked_smiles.startswith("*")
        assert original_smiles in marked_smiles
        
        # Generate conformers from marked SMILES
        conformers = conformers_from_smile(marked_smiles, conformer_count=2)
        assert len(conformers) == 2
    
    def test_geometry_operations_integration(self):
        """Test integration of geometry operations."""
        # Create vectors
        vec1 = np.array([1, 0, 0])
        vec2 = np.array([0, 1, 0])
        
        # Get rotation matrix
        rotation_matrix = rotation_matrix_from_vectors(vec1, vec2)
        assert rotation_matrix.shape == (3, 3)
        
        # Apply to atoms
        atoms = Atoms("CO", positions=[[0, 0, 0], [1, 0, 0]])
        original_positions = atoms.positions.copy()
        
        atoms.positions = atoms.positions @ rotation_matrix.T
        assert not np.allclose(atoms.positions, original_positions)


class TestVisualizationIntegration:
    """Integration tests for visualization."""
    
    def test_conformer_plotting(self):
        """Test conformer plotting functionality."""
        # Create fragment with multiple conformers
        fragment = Fragment(smile="COC", to_initialize=4)
        
        # Plot conformers
        fig = docs_plot_conformers(fragment.conformers)
        assert fig is not None
        
        # Test with different number of conformers
        fragment_small = Fragment(smile="COC", to_initialize=2)
        fig_small = docs_plot_conformers(fragment_small.conformers)
        assert fig_small is not None
    
    def test_surface_plotting(self):
        """Test surface plotting functionality."""
        from autoadsorbate.viz import docs_plot_sites
        
        # Create surface
        slab = fcc111("Cu", (2, 2, 2), periodic=True, vacuum=10)
        surface = Surface(slab)
        
        # Plot surface with sites
        fig = docs_plot_sites(surface)
        assert fig is not None


class TestCompleteWorkflow:
    """Test complete workflow from SMILES to surface adsorption."""
    
    def test_complete_workflow(self):
        """Test the complete workflow from SMILES to surface adsorption."""
        # 1. Create surface
        slab = fcc111("Cu", (2, 2, 2), periodic=True, vacuum=10)
        surface = Surface(slab)
        
        # 2. Create fragment
        fragment = Fragment(smile="COC", to_initialize=3)
        
        # 3. Get surface normal
        normal = surface.get_surface_normal()
        
        # 4. Get aligned conformer
        conformer = fragment.get_conformer(0, n_vector=normal)
        
        # 5. Create intermediate
        if surface.site_dict:
            site_id = list(surface.site_dict.keys())[0]
            site_info = surface.get_site(site_id)
            intermediate = Intermediate(site_info, [fragment])
            
            assert intermediate.ActiveSite == site_info
            assert len(intermediate.fragments) == 1
            assert intermediate.fragments[0] is fragment
        
        # 6. Verify all objects are properly created
        assert isinstance(surface, Surface)
        assert isinstance(fragment, Fragment)
        assert isinstance(conformer, Atoms)
        assert len(conformer) > 0
    
    def test_workflow_with_multiple_molecules(self):
        """Test workflow with multiple different molecules."""
        # Create surface
        slab = fcc111("Cu", (2, 2, 2), periodic=True, vacuum=10)
        surface = Surface(slab)
        
        # Create multiple fragments
        molecules = ["COC", "CCO", "c1ccccc1"]  # Dimethyl ether, ethanol, benzene
        fragments = []
        
        for smiles in molecules:
            fragment = Fragment(smile=smiles, to_initialize=2)
            fragments.append(fragment)
            
            # Get conformer
            conformer = fragment.get_conformer(0)
            assert isinstance(conformer, Atoms)
            assert len(conformer) > 0
        
        assert len(fragments) == 3
        assert all(isinstance(f, Fragment) for f in fragments)


class TestErrorHandling:
    """Test error handling in integration scenarios."""
    
    def test_invalid_smiles_in_fragment(self):
        """Test handling of invalid SMILES in Fragment creation."""
        with pytest.raises(Exception):
            Fragment(smile="invalid_smiles", to_initialize=1)
    
    def test_empty_surface_sites(self):
        """Test handling of surface with no detected sites."""
        # Create a very small surface that might not have sites
        slab = fcc111("Cu", (1, 1, 1), periodic=True, vacuum=5)
        surface = Surface(slab, precision=0.1, touch_sphere_size=1)
        
        # Should not raise an error even if no sites are found
        assert isinstance(surface.site_dict, dict)
    
    def test_fragment_with_zero_conformers(self):
        """Test Fragment with zero conformers."""
        with pytest.raises(Exception):
            Fragment(smile="COC", to_initialize=0)
