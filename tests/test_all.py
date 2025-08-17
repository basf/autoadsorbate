"""Comprehensive test suite for autoadsorbate."""

import pytest
import numpy as np
from ase import Atoms
from ase.build import fcc111

from autoadsorbate import Fragment, Surface, Intermediate


def test_Surface():
    """Test Surface class initialization."""
    slab = fcc111("Cu", (2, 2, 2), periodic=True, vacuum=10)
    s = Surface(slab)
    assert type(s.site_dict) == dict
    assert hasattr(s, 'atoms')
    assert hasattr(s, 'site_dataframe')


def test_Fragment():
    """Test Fragment class initialization."""
    f = Fragment(smile="COC", to_initialize=5)
    assert f.smile == "COC"
    assert len(f.conformers) == 5
    assert all(isinstance(conf, Atoms) for conf in f.conformers)


def test_Intermediate():
    """Test Intermediate class initialization."""
    active_site = {"type": "top", "position": [0, 0, 0]}
    i = Intermediate(active_site)
    assert i.ActiveSite == active_site
    assert i.fragments == []


def test_basic_workflow():
    """Test basic workflow from SMILES to surface adsorption."""
    # Create surface
    slab = fcc111("Cu", (2, 2, 2), periodic=True, vacuum=10)
    surface = Surface(slab)
    
    # Create fragment
    fragment = Fragment(smile="COC", to_initialize=3)
    
    # Get conformer
    conformer = fragment.get_conformer(0)
    
    # Verify all objects are properly created
    assert isinstance(surface, Surface)
    assert isinstance(fragment, Fragment)
    assert isinstance(conformer, Atoms)
    assert len(conformer) > 0


def test_surface_properties():
    """Test surface property calculations."""
    slab = fcc111("Cu", (2, 2, 2), periodic=True, vacuum=10)
    surface = Surface(slab)
    
    # Test surface area
    area = surface.get_surface_area()
    assert isinstance(area, float)
    assert area > 0
    
    # Test surface normal
    normal = surface.get_surface_normal()
    assert isinstance(normal, np.ndarray)
    assert len(normal) == 3
    assert np.isclose(np.linalg.norm(normal), 1.0)


def test_fragment_conformers():
    """Test fragment conformer generation and manipulation."""
    fragment = Fragment(smile="CCO", to_initialize=4)
    
    # Test conformer retrieval
    for i in range(len(fragment.conformers)):
        conformer = fragment.get_conformer(i)
        assert isinstance(conformer, Atoms)
        assert len(conformer) > 0
    
    # Test invalid index
    with pytest.raises(ValueError):
        fragment.get_conformer(10)


def test_fragment_alignment():
    """Test fragment alignment with custom vectors."""
    fragment = Fragment(smile="COC", to_initialize=2)
    
    # Test alignment with custom normal vector
    n_vector = np.array([1, 0, 0])
    conformer = fragment.get_conformer(0, n_vector=n_vector, rot_deg=45)
    assert isinstance(conformer, Atoms)


if __name__ == "__main__":
    # Run basic tests
    test_Surface()
    test_Fragment()
    test_Intermediate()
    test_basic_workflow()
    test_surface_properties()
    test_fragment_conformers()
    test_fragment_alignment()
    print("All basic tests passed!")
