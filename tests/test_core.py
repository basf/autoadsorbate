"""Tests for core classes."""

import pytest
import numpy as np
from ase import Atoms
from ase.build import fcc111

from autoadsorbate.core import Fragment, Surface, Intermediate


class TestFragment:
    """Test cases for Fragment class."""
    
    def test_fragment_initialization(self):
        """Test Fragment initialization with basic parameters."""
        fragment = Fragment(smile="COC", to_initialize=5)
        assert fragment.smile == "COC"
        assert fragment.to_initialize == 5
        assert fragment.randomSeed == 2104
        assert len(fragment.conformers) == 5
        assert len(fragment.conformers_aligned) == 5
        assert all(not aligned for aligned in fragment.conformers_aligned)
    
    def test_fragment_custom_parameters(self):
        """Test Fragment initialization with custom parameters."""
        fragment = Fragment(
            smile="CCO", 
            to_initialize=3, 
            random_seed=42, 
            sort_conformers=False
        )
        assert fragment.smile == "CCO"
        assert fragment.to_initialize == 3
        assert fragment.randomSeed == 42
        assert fragment.sort_conformers == False
        assert len(fragment.conformers) == 3
    
    def test_get_conformer(self):
        """Test getting a specific conformer."""
        fragment = Fragment(smile="COC", to_initialize=3)
        conformer = fragment.get_conformer(0)
        assert isinstance(conformer, Atoms)
        assert len(conformer) > 0
        assert fragment.conformers_aligned[0] == True
    
    def test_get_conformer_invalid_index(self):
        """Test getting conformer with invalid index."""
        fragment = Fragment(smile="COC", to_initialize=3)
        with pytest.raises(ValueError):
            fragment.get_conformer(5)
    
    def test_get_conformer_with_rotation(self):
        """Test getting conformer with custom rotation."""
        fragment = Fragment(smile="COC", to_initialize=3)
        n_vector = np.array([1, 0, 0])
        conformer = fragment.get_conformer(0, n_vector=n_vector, rot_deg=45)
        assert isinstance(conformer, Atoms)


class TestSurface:
    """Test cases for Surface class."""
    
    def test_surface_initialization(self):
        """Test Surface initialization."""
        slab = fcc111("Cu", (2, 2, 2), periodic=True, vacuum=10)
        surface = Surface(slab)
        assert isinstance(surface.atoms, Atoms)
        assert isinstance(surface.site_dict, dict)
        assert surface.precision == 0.5
        assert surface.touch_sphere_size == 3
    
    def test_surface_custom_parameters(self):
        """Test Surface initialization with custom parameters."""
        slab = fcc111("Cu", (2, 2, 2), periodic=True, vacuum=10)
        surface = Surface(
            slab, 
            precision=1.0, 
            touch_sphere_size=5, 
            sort_sites=False
        )
        assert surface.precision == 1.0
        assert surface.touch_sphere_size == 5
    
    def test_get_site(self):
        """Test getting a specific site."""
        slab = fcc111("Cu", (2, 2, 2), periodic=True, vacuum=10)
        surface = Surface(slab)
        
        if surface.site_dict:
            site_id = list(surface.site_dict.keys())[0]
            site_info = surface.get_site(site_id)
            assert isinstance(site_info, dict)
    
    def test_get_site_invalid_id(self):
        """Test getting site with invalid ID."""
        slab = fcc111("Cu", (2, 2, 2), periodic=True, vacuum=10)
        surface = Surface(slab)
        
        with pytest.raises(ValueError):
            surface.get_site(999)
    
    def test_get_sites_by_type(self):
        """Test getting sites by type."""
        slab = fcc111("Cu", (2, 2, 2), periodic=True, vacuum=10)
        surface = Surface(slab)
        
        sites = surface.get_sites_by_type("top")
        assert isinstance(sites, list)
    
    def test_get_surface_area(self):
        """Test surface area calculation."""
        slab = fcc111("Cu", (2, 2, 2), periodic=True, vacuum=10)
        surface = Surface(slab)
        
        area = surface.get_surface_area()
        assert isinstance(area, float)
        assert area > 0
    
    def test_get_surface_normal(self):
        """Test surface normal calculation."""
        slab = fcc111("Cu", (2, 2, 2), periodic=True, vacuum=10)
        surface = Surface(slab)
        
        normal = surface.get_surface_normal()
        assert isinstance(normal, np.ndarray)
        assert len(normal) == 3
        assert np.isclose(np.linalg.norm(normal), 1.0)


class TestIntermediate:
    """Test cases for Intermediate class."""
    
    def test_intermediate_initialization(self):
        """Test Intermediate initialization."""
        active_site = {"type": "top", "position": [0, 0, 0]}
        intermediate = Intermediate(active_site)
        assert intermediate.ActiveSite == active_site
        assert intermediate.fragments == []
    
    def test_intermediate_with_fragments(self):
        """Test Intermediate initialization with fragments."""
        active_site = {"type": "bridge", "position": [1, 1, 0]}
        fragments = [Fragment("COC", to_initialize=2)]
        intermediate = Intermediate(active_site, fragments)
        assert intermediate.ActiveSite == active_site
        assert len(intermediate.fragments) == 1
        assert isinstance(intermediate.fragments[0], Fragment)
