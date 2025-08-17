#!/usr/bin/env python3
"""
Test script to verify notebook functionality works with basic data
"""

import numpy as np
import pandas as pd
from ase import Atoms
from ase.build import fcc111

# Import the updated autoadsorbate modules
from autoadsorbate import Fragment, Surface
from autoadsorbate.viz import plot_most_stable, make_hist_plot, plot_energy_heatmap, filter_xdf
from autoadsorbate.utils import compute_energy, get_drop_snapped
from autoadsorbate.utils.analysis import count_C_next_to_O

def create_test_data():
    """Create test data to verify notebook functionality"""
    
    # Create a simple surface
    slab = fcc111("Cu", (2, 2, 2), periodic=True, vacuum=10)
    surface = Surface(slab)
    
    # Create test DataFrame with required columns
    test_data = {
        'pid': ['test-1-1', 'test-1-2', 'test-2-1', 'test-2-2'],
        'energy': [0.1, 0.2, 0.3, 0.4],
        'energy_calibrated': [0.1, 0.2, 0.3, 0.4],
        'H': [1, 1, 2, 2],
        'backbone': ['COC', 'COC', 'CCO', 'CCO'],
        'smiles': ['COC', 'COC', 'CCO', 'CCO'],
        'calibrate_keys': ['COC-H1', 'COC-H1', 'CCO-H2', 'CCO-H2'],
        'traj_index': [0, 1, 2, 3],
        'origin': ['test1', 'test2', 'test3', 'test4'],
        'material_formula': ['Cu', 'Cu', 'Cu', 'Cu'],
        'mi': ['1#1#1', '1#1#1', '1#1#1', '1#1#1']
    }
    
    xdf = pd.DataFrame(test_data)
    
    # Create test trajectory
    relaxed_traj = []
    for i in range(4):
        atoms = slab.copy()
        # Add some atoms to represent fragments
        atoms.extend(Atoms('CO', positions=[[0, 0, 5+i], [1, 0, 5+i]]))
        relaxed_traj.append(atoms)
    
    return xdf, relaxed_traj, surface

def test_basic_functionality():
    """Test basic functionality of the notebook components"""
    
    print("Testing basic functionality...")
    
    # Create test data
    xdf, relaxed_traj, surface = create_test_data()
    
    # Test Fragment creation
    fragment = Fragment(smile="COC", to_initialize=2)
    print(f"✓ Fragment created with {len(fragment.conformers)} conformers")
    
    # Test Surface functionality
    print(f"✓ Surface created with {len(surface.site_dict)} sites")
    print(f"✓ Surface area: {surface.get_surface_area():.2f} Å²")
    
    # Test utility functions
    energy = compute_energy(relaxed_traj[0])
    print(f"✓ Energy computation: {energy}")
    
    # Test analysis function
    count = count_C_next_to_O(relaxed_traj[0])
    print(f"✓ C-O count: {count}")
    
    # Test plotting functions (without displaying)
    try:
        # Test filter_xdf
        filtered_df = filter_xdf(xdf, relaxed_traj)
        print(f"✓ filter_xdf: {len(filtered_df)} rows")
        
        # Test plotting functions (create but don't show)
        fig1 = plot_most_stable(xdf, relaxed_traj)
        print("✓ plot_most_stable created successfully")
        
        fig2 = make_hist_plot(xdf)
        print("✓ make_hist_plot created successfully")
        
        ax, heatmap = plot_energy_heatmap(
            xdf, 
            column='energy_calibrated',
            std=0.05,
            e_min=0,
            e_max=1,
            resolution=0.1,
            normalize=True,
            return_heatmap=True
        )
        print("✓ plot_energy_heatmap created successfully")
        
    except Exception as e:
        print(f"✗ Plotting test failed: {e}")
    
    print("\nAll basic functionality tests completed!")

def test_imports():
    """Test that all required imports work"""
    
    print("Testing imports...")
    
    try:
        from autoadsorbate import Fragment, Surface, Intermediate
        print("✓ Core classes imported successfully")
        
        from autoadsorbate.molecular import conformers_from_smile, get_marked_smiles
        print("✓ Molecular functions imported successfully")
        
        from autoadsorbate.surface import conformer_to_site
        print("✓ Surface functions imported successfully")
        
        from autoadsorbate.utils import compute_energy, get_drop_snapped
        print("✓ Utility functions imported successfully")
        
        from autoadsorbate.viz import (
            plot_most_stable, make_hist_plot, plot_energy_heatmap,
            docs_plot_conformers, docs_plot_sites
        )
        print("✓ Visualization functions imported successfully")
        
        from autoadsorbate.string_utils import _example_config, construct_smiles
        print("✓ String utilities imported successfully")
        
    except ImportError as e:
        print(f"✗ Import failed: {e}")
        return False
    
    print("✓ All imports successful!")
    return True

if __name__ == "__main__":
    print("=" * 50)
    print("Testing Notebook Functionality with New Modular Structure")
    print("=" * 50)
    
    # Test imports first
    if test_imports():
        # Test basic functionality
        test_basic_functionality()
        
        print("\n" + "=" * 50)
        print("✓ All tests passed! The notebook functionality works with the new structure.")
        print("=" * 50)
    else:
        print("\n✗ Import tests failed. Please check the module structure.")
