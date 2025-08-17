"""Surface class for handling surface structures."""

import numpy as np
import pandas as pd
from ase import Atoms
from ase.build.tools import sort as sort_atoms

from .site_detection import get_shrinkwrap_ads_sites


class Surface:
    """
    Base class for initializing surface structures.

    Attributes:
        atoms (Atoms): The ASE Atoms object representing the surface.
        site_dict (dict): Dictionary containing adsorption site information.
        site_dataframe (pd.DataFrame): DataFrame containing site information.
    """

    def __init__(
        self,
        atoms: Atoms,
        precision: float = 0.5,
        touch_sphere_size: float = 3,
        sort_sites: bool = True,
    ):
        """
        Initialize surface with adsorption site detection.

        Args:
            atoms (Atoms): The ASE Atoms object representing the surface.
            precision (float, optional): Precision for site detection. Defaults to 0.5.
            touch_sphere_size (float, optional): Size of touch sphere for site detection. Defaults to 3.
            sort_sites (bool, optional): Whether to sort sites. Defaults to True.
        """
        self.atoms = atoms.copy()
        self.precision = precision
        self.touch_sphere_size = touch_sphere_size
        
        # Detect adsorption sites
        self.site_dict, self.site_dataframe = get_shrinkwrap_ads_sites(
            self.atoms, 
            precision=precision, 
            touch_sphere_size=touch_sphere_size
        )
        
        if sort_sites:
            self._sort_sites()
    
    def _sort_sites(self):
        """Sort sites by their z-coordinate (height)."""
        if self.site_dataframe is not None and len(self.site_dataframe) > 0:
            self.site_dataframe = self.site_dataframe.sort_values('z', ascending=False)
    
    def get_site(self, site_id: int) -> dict:
        """
        Get information about a specific adsorption site.

        Args:
            site_id (int): The ID of the site to retrieve.

        Returns:
            dict: Dictionary containing site information.
        """
        if site_id not in self.site_dict:
            raise ValueError(f"Site ID {site_id} not found. Available sites: {list(self.site_dict.keys())}")
        
        return self.site_dict[site_id]
    
    def get_sites_by_type(self, site_type: str) -> list:
        """
        Get all sites of a specific type.

        Args:
            site_type (str): The type of sites to retrieve (e.g., 'top', 'bridge', '3-fold').

        Returns:
            list: List of site dictionaries of the specified type.
        """
        sites = []
        for site_id, site_info in self.site_dict.items():
            if site_info.get('type') == site_type:
                sites.append(site_info)
        return sites
    
    def get_surface_area(self) -> float:
        """
        Calculate the surface area of the surface.

        Returns:
            float: The surface area in square angstroms.
        """
        # Calculate surface area based on cell dimensions
        cell = self.atoms.cell
        area = np.linalg.norm(np.cross(cell[0], cell[1]))
        return area
    
    def get_surface_normal(self) -> np.ndarray:
        """
        Get the surface normal vector.

        Returns:
            np.ndarray: The surface normal vector.
        """
        cell = self.atoms.cell
        normal = np.cross(cell[0], cell[1])
        normal /= np.linalg.norm(normal)
        return normal
