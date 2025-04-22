"""Main module."""

from typing import Dict, List

import ase
import numpy as np
import pandas as pd
from ase import Atom, Atoms
from ase.visualize import view

from .Smile import (
    _reset_position,
    _reset_rotation,
    conformers_from_smile,
)
from .Surf import conformer_to_site, get_shrinkwrap_ads_sites
from .utils import (
    get_sorted_by_snap_dist,
    make_site_info_writable,
)


class Intermediate:
    """
    Base class for initializing reaction intermediates.

    Attributes:
        ActiveSite: The active site for the intermediate.
        fragments: A list of fragments associated with the intermediate.
    """

    def __init__(self, ActiveSite, fragments=None):
        """
        Initialize attributes.

        Args:
            ActiveSite: The active site for the intermediate.
            fragments (list, optional): A list of fragments associated with the intermediate. Defaults to an empty list.
        """
        self.ActiveSite = ActiveSite
        self.fragments = fragments if fragments is not None else []


class Fragment:
    """
    Base class for initializing reaction fragments.

    Attributes:
        smile (str): The SMILES string of the fragment.
        to_initialize (int): The number of conformers to initialize.
        randomSeed (int): The random seed for conformer generation.
        conformers (List[Atoms]): A list of ASE Atoms objects representing the conformers.
        conformers_aligned (List[bool]): A list indicating whether each conformer is aligned.
    """

    def __init__(
        self,
        smile: str,
        to_initialize: int = 10,
        random_seed: int = 2104,
        sort_conformers: bool = True,
    ):
        """
        Initialize attributes.

        Args:
            smile (str): The SMILES string of the fragment.
            to_initialize (int, optional): The number of conformers to initialize. Defaults to 10.
            random_seed (int, optional): The random seed for conformer generation. Defaults to 2104.
            sort_conformers (bool, optional): Decides if the initial orientation of the fragment conformations is diverse.
        """
        self.smile = smile
        self.to_initialize = to_initialize
        self.randomSeed = random_seed

        self.conformers = conformers_from_smile(
            smile, to_initialize, random_seed=random_seed
        )
        self.conformers_aligned = [False for _ in self.conformers]

        self.sort_conformers = sort_conformers
        if self.sort_conformers:
            self.conformers = get_sorted_by_snap_dist(self.conformers)

    def get_conformer(
        self, i: int, n_vector: np.ndarray = np.array([0, 0, 1]), rot_deg: float = 0
    ) -> Atoms:
        """
        Returns a copy of the i-th conformer, aligned and rotated as specified.

        Args:
            i (int): The index of the conformer to retrieve.
            n_vector (np.ndarray, optional): The normal vector for rotation. Defaults to [0, 0, 1].
            rot_deg (float, optional): The rotation angle in degrees. Defaults to 0.

        Returns:
            Atoms: A copy of the aligned and rotated conformer.
        """
        if not self.conformers_aligned[i]:
            self.conformers[i] = _reset_position(self.conformers[i])
            self.conformers[i] = _reset_rotation(self.conformers[i])
            self.conformers_aligned[i] = True

        self.conformers[i].rotate(rot_deg, n_vector)
        self.conformers[i].info["smiles"] = self.smile
        return self.conformers[i].copy()

    def view(self, return_traj=False):
        traj = [self.get_conformer(i) for i, _ in enumerate(self.conformers)]
        if return_traj:
            return traj
        view(traj)


class Surface:
    """
    Base class for initializing a reactive surface.

    Attributes:
        atoms (Atoms): The ASE Atoms object representing the surface.
        precision (float): The precision for the grid spacing.
        touch_sphere_size (float): The size of the touch sphere.
        site_dict (Dict): Dictionary containing site information.
        site_df (pd.DataFrame): DataFrame containing site information.
    """

    def __init__(
        self,
        atoms: Atoms,
        precision: float = 0.25,
        touch_sphere_size: float = 3,
        dummy: bool = False,
    ):
        """
        Initialize attributes.

        Args:
            atoms (Atoms): The ASE Atoms object representing the surface.
            precision (float, optional): The precision for the grid spacing. Defaults to 0.25.
            touch_sphere_size (float, optional): The size of the touch sphere. Defaults to 3.
        """
        self.dummy = dummy
        self.atoms = atoms
        self.precision = precision
        self.touch_sphere_size = touch_sphere_size
        if self.dummy:
            self.site_dict = {}
            self.site_df = pd.DataFrame(self.site_dict)
        else:
            self.site_dict = get_shrinkwrap_ads_sites(
                atoms=self.atoms,
                precision=self.precision,
                touch_sphere_size=self.touch_sphere_size,
            )
            self.site_df = pd.DataFrame(self.site_dict)
            self.sort_site_df()

    def sort_site_df(self, by: str = "xyz"):
        """
        Sorts the site DataFrame by coordinates or distance.

        Args:
            by (str, optional): The sorting criterion ('xyz' or 'dist'). Defaults to 'xyz'.
        """
        if by == "xyz":
            sort = {}
            for c in [0, 1, 2]:
                sort[f"sort_{c}"] = [
                    np.round(coord[c], 1) for coord in self.site_df.coordinates
                ]
            for k, v in sort.items():
                self.site_df[k] = v
            self.site_df = self.site_df.sort_values(
                by=list(sort.keys()), ignore_index=True
            )
            for k in sort.keys():
                self.site_df.pop(k)
        elif by == "dist":
            self.site_df["sort"] = [
                np.round(np.linalg.norm(coord), 1) for coord in self.site_df.coordinates
            ]
            self.site_df = self.site_df.sort_values(by="sort", ignore_index=True)
            self.site_df.pop("sort")

    def get_site(self, index: int) -> Atoms:
        """
        Returns the atoms object for a specific site.

        Args:
            index (int): The index of the site.

        Returns:
            Atoms: The ASE Atoms object for the site.
        """
        site_atoms = self.atoms.copy()

        if "adsorbate_info" in site_atoms.info.keys():
            site_atoms.info.pop("adsorbate_info")

        info = self.site_df.loc[index].to_dict()
        site_atoms.info.update(info)
        site_atoms.append(Atom("X", position=self.site_df["coordinates"].loc[index]))
        del site_atoms[:-1]
        return site_atoms

    def view_site(self, index: int, return_atoms: bool = False) -> Atoms:
        """
        Visualizes a specific site.

        Args:
            index (int): The index of the site.
            return_atoms (bool, optional): Whether to return the atoms object. Defaults to False.

        Returns:
            Atoms: The ASE Atoms object for the site if return_atoms is True.
        """
        site_atoms = self.get_site(index)
        site_atoms += self.atoms[site_atoms.info["topology"]]
        for x in [np.round(x, 1) for x in np.arange(0.1, 2.1, 0.1)]:
            site_atoms.append(
                Atom(
                    "X",
                    position=site_atoms.info["coordinates"]
                    + site_atoms.info["n_vector"] * x,
                )
            )
        if return_atoms:
            return site_atoms
        else:
            view(site_atoms)

    def view_surface(
        self,
        return_atoms: bool = False,
        explicit_marker: str = None,
        mode: str = "normal",
    ) -> Atoms:
        """
        Visualizes the entire surface.

        Args:
            return_atoms (bool, optional): Whether to return the atoms object. Defaults to False.
            explicit_marker (str, optional): The marker symbol to use for visualization. Defaults to None.
            mode (str, optional): choose from ['normal', hedgehog]. Defaults to 'normal'.

        Returns:
            Atoms: The ASE Atoms object for the surface if return_atoms is True.
        """
        view_atoms = self.atoms.copy()
        inds = list(set([i for ind_ls in self.site_dict["topology"] for i in ind_ls]))

        if explicit_marker:
            for i in inds:
                view_atoms[i].symbol = explicit_marker
        else:
            marker_map = _get_marker_map(self.atoms)
            for i in inds:
                view_atoms[i].number = marker_map[view_atoms[i].number]

        if mode == "hedgehog":
            for i in self.site_df.index.values:
                view_atoms += self.view_site(i, return_atoms=True)[
                    [
                        atom.index
                        for atom in self.view_site(i, return_atoms=True)
                        if atom.symbol == "X"
                    ]
                ]
        if return_atoms:
            return view_atoms
        else:
            view(view_atoms)

    def compare_sites(self, site_index1: int, site_index2: int) -> bool:
        """
        Compares two sites for symmetry equivalence.

        Args:
            site_index1 (int): The index of the first site.
            site_index2 (int): The index of the second site.

        Returns:
            bool: True if the sites are equivalent, False otherwise.
        """
        from ase.utils.structure_comparator import SymmetryEquivalenceCheck

        SEC = SymmetryEquivalenceCheck()
        site1 = self.get_site(site_index1)
        site2 = self.get_site(site_index2)

        for s in [site1, site2]:
            if len(s) == 1:
                s.positions += [0, 0, 0.01]

        return SEC.compare(self.atoms + site1, self.atoms + site2)

    def get_nonequivalent_sites(self) -> List[int]:
        """
        Returns a list of indices for nonequivalent sites.

        Returns:
            List[int]: A list of indices for nonequivalent sites.
        """
        original = []
        i_s = self.site_df.index.values
        matches = np.array([False for _ in i_s])

        for i in i_s:
            if not matches[i]:
                m = [self.compare_sites(i, j) for j in i_s]
                matches += m
                matches = matches > 0
                original.append(i)
            if all(matches):
                break
        return original

    def sym_reduce(self):
        """
        Reduces the site DataFrame to nonequivalent sites.
        """
        include = self.get_nonequivalent_sites()
        include_filter = [i in include for i in self.site_df.index.values]
        self.site_df = self.site_df[include_filter]
        self.site_dict = self.site_df.to_dict(orient="list")

    def get_populated_sites(
        self,
        fragment,
        site_index="all",
        sample_rotation=True,
        mode="heuristic",
        conformers_per_site_cap=None,
        overlap_thr=1.5,
        verbose=False,
    ):
        """
        Populates the specified sites with the given fragment, optimizing the orientation to minimize overlap.

        Parameters:
        fragment (object): An object containing the fragment to be attached.
        site_index (str or int): The index of the site to be populated. Default is 'all'.
        sample_rotation (bool): Whether to sample different rotations of the fragment. Default is True.
        mode (str): The mode of operation. Can be 'heuristic' or 'all'. Default is 'heuristic'.
        conformers_per_site_cap (int or None): The maximum number of conformers per site. Default is None.
        overlap_thr (float): The overlap threshold. Default is 1.5.
        verbose (bool): Whether to print detailed information during execution. Default is False.

        Returns:
        list: A list containing the optimized atoms objects for each site.

        Raises:
        ValueError: If the mode is not implemented or if the fragment object is invalid.
        """

        all_sites = {}
        site_df = self.site_df

        if mode.lower() == "all":
            sites = [site_df.loc[i].to_dict() for i in site_df.index.values]

        elif mode.lower() == "heuristic":
            all_sites["S1S"] = [
                site_df.loc[i].to_dict()
                for i in site_df[site_df.connectivity > 1].index.values
            ]
            all_sites["Cl"] = [
                site_df.loc[i].to_dict()
                for i in site_df[site_df.connectivity == 1].index.values
            ]

            if fragment.smile[:3] == "S1S":
                sites = all_sites["S1S"]

            if fragment.smile[:2] == "Cl":
                sites = all_sites["Cl"]

        else:
            raise ValueError("argument 'mode' can be 'heuristic' or 'all'")
            return

        if sample_rotation:
            conformers = []
            for i, _ in enumerate(fragment.conformers):
                c = fragment.get_conformer(i)
                if c.info["smiles"][:3] == "S1S":
                    angles = [0, 180]
                if c.info["smiles"][:2] == "Cl":
                    angles = [a for a in range(0, 360, 45)]
                for a in angles:
                    ca = c.copy()
                    ca.rotate(a, "z")
                    conformers.append(ca)
        else:
            conformers = [c.copy() for c in fragment.conformers]

        out_trj = []

        if verbose:
            print("conformers", len(conformers))
            print("sites", len(sites))

        for site in sites:
            c_trj = []
            for conformer in conformers:
                c_trj += conformer_to_site(
                    self.atoms, site, conformer, mode="optimize", overlap_thr=0
                )  # the zero is intentional

            if conformers_per_site_cap != None:
                c_trj = [atoms for atoms in c_trj if atoms.info["mdf"] > overlap_thr]

                if len(c_trj) > 1:
                    c_trj = get_sorted_by_snap_dist(c_trj)[
                        : int(np.min([conformers_per_site_cap, len(c_trj)]))
                    ]
                elif len(c_trj) == 0:
                    c_trj = []
                else:
                    pass

                if len(c_trj) < conformers_per_site_cap and verbose:
                    print(
                        f"WARNING: Failed to find requested number of conformers with condition: ovelap_thr = {overlap_thr}. Found {len(c_trj)} / {conformers_per_site_cap}. Consider setting a higher Fragment(to_initialize = < N >)"
                    )
                if len(c_trj) == conformers_per_site_cap and verbose:
                    print(
                        f"SUCCESS! Found the requested numer of conformers with condition: ovelap_thr = {overlap_thr}. Found {len(c_trj)} / {conformers_per_site_cap}."
                    )

                for atoms in c_trj:
                    atoms.info["adsorbate_info"] = {}
                    atoms.info["adsorbate_info"]["site"] = make_site_info_writable(site)
                    atoms.info["adsorbate_info"]["smiles"] = fragment.smile
                    atoms.info["adsorbate_info"]["mdf"] = atoms.info.pop("mdf")
                    formula = atoms[
                        [
                            atom.index
                            for atom in atoms
                            if atoms.arrays["fragments"][atom.index]
                            == max(atoms.arrays["fragments"])
                        ]
                    ].get_chemical_formula()
                    atoms.info["adsorbate_info"]["adsorbate_formula"] = formula

                out_trj += c_trj

        return out_trj


def _get_marker_map(atoms: Atoms) -> Dict[int, int]:
    """
    Generates a marker map for visualizing the surface.

    Args:
        atoms (Atoms): The ASE Atoms object representing the surface.

    Returns:
        Dict[int, int]: A dictionary mapping atomic numbers to marker atomic numbers.
    """
    marker_map = {}
    at_no = set(atoms.get_atomic_numbers())
    for atomic_number in at_no:
        done = False
        swap_atomic_number = atomic_number + 1
        while not done:
            if swap_atomic_number not in at_no:
                done = True
            else:
                swap_atomic_number += 1
        marker_map[atomic_number] = swap_atomic_number

    for k, v in marker_map.items():
        print(
            f"Visualizing surface {ase.symbols.chemical_symbols[k]} atoms as {ase.symbols.chemical_symbols[v]}"
        )
        # print(f'Visualizing surface {ase.symbols.symbols([k]).get_chemical_formula()} atoms as {ase.symbols.symbols([v]).get_chemical_formula()}')
    return marker_map


class ActiveSite(Surface):
    """
    Base class for initializing reaction fragments.

    Attributes:
        atoms (Atoms): The ASE Atoms object representing the surface.
        must_include (List[int]): List of atom indices that must be included in the active site.
        must_exclude (List[int]): List of atom indices that must be excluded from the active site.
        keep_tops (bool): Whether to keep top sites.
    """

    def __init__(
        self,
        atoms: Atoms,
        must_include: List[int] = [],
        must_exclude: List[int] = [],
        keep_tops: bool = True,
    ):
        """
        Initialize attributes.

        Args:
            atoms (Atoms): The ASE Atoms object representing the surface.
            must_include (List[int], optional): List of atom indices that must be included in the active site. Defaults to [].
            must_exclude (List[int], optional): List of atom indices that must be excluded from the active site. Defaults to [].
            keep_tops (bool, optional): Whether to keep top sites. Defaults to True.
        """
        super().__init__(atoms)
        self.must_include = must_include
        self.must_exclude = must_exclude
        self.keep_tops = keep_tops
        self.pop_sites()

    def pop_sites(self):
        """
        Filters the sites based on must_include and must_exclude lists.
        """
        if not self.must_include:
            return

        include_filter = np.array(
            [
                any(i in v for i in self.must_include)
                for v in self.site_df.topology.values
            ]
        )

        for e in self.must_exclude:
            include_filter &= np.array(
                [e not in v for v in self.site_df.topology.values]
            )

        if self.keep_tops:
            temp_site_df = self.site_df[include_filter]
            tops = list(set(t for v in temp_site_df.topology.values for t in v))

            for i, v in enumerate(self.site_df.topology.values):
                if v in [[t] for t in tops]:
                    include_filter[i] = True

        self.site_df = self.site_df[include_filter]
        self.site_dict = self.site_df.to_dict(orient="list")
