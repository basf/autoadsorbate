import itertools
import math
from copy import deepcopy

import numpy as np
import pandas as pd
from ase import Atom, Atoms
from scipy.spatial.distance import cdist

from .utils import polar2cart, rotation_matrix_from_vectors


def get_shrinkwrap_site_n_vector(
    atoms: Atoms, site: np.ndarray, grid: Atoms, touch_sphere_size: float
) -> np.ndarray:
    """Takes atoms, site, and shrinkwrap grid, and returns a site vector based on
    grid points in the touch_sphere_size * 1.1 radius.

    Args:
        atoms (Atoms): Atoms slab.
        site (np.ndarray): Site (point in 3D space on the surface).
        grid (Atoms): Shrinkwrap grid of atoms.
        touch_sphere_size (float): Radius to consider for grid points.

    Returns:
        np.ndarray: Normalized vector of the site.
    """

    atoms_copy = atoms.copy()
    grid_copy = grid.copy()
    shift_vector = atoms_copy.cell[0] * 0.5 + atoms_copy.cell[1] * 0.5 - site
    atoms_copy.positions += shift_vector
    grid_copy.positions += shift_vector
    shifted_site = site + shift_vector
    atoms_copy.wrap()
    grid_copy.wrap()

    grid_distances = np.linalg.norm(grid_copy.positions - shifted_site, axis=1)
    # grid_site_indices = np.where(grid_distances < touch_sphere_size * 1.2)[0]
    grid_site_indices = np.where(grid_distances < np.min(grid_distances) * 1.2)[0]

    if len(grid_site_indices) == 0:
        raise ValueError(
            f"No grid points found within the specified radius. Closest grid point at {np.min(grid_distances)}, cutoff {touch_sphere_size * 1.2}"
        )

    center_of_mass = grid_copy[grid_site_indices].get_center_of_mass()
    n_vector = center_of_mass - shifted_site

    n_vector /= np.linalg.norm(n_vector)

    return n_vector


def get_shrinkwrap_site_h_vector(
    ads_site_atoms: Atoms, n_vector: np.ndarray
) -> np.ndarray:
    """Calculates the h_vector for a given site based on the adsorbed site atoms and n_vector.

    Args:
        ads_site_atoms (Atoms): Atoms at the adsorption site.
        n_vector (np.ndarray): Normalized vector of the site.

    Returns:
        np.ndarray: Normalized h_vector of the site.
    """
    if len(ads_site_atoms) == 0:
        raise ValueError("ads_site_atoms cannot be empty.")

    if len(ads_site_atoms) == 1:
        h_vector = np.array([1.0, 0.0, 0.0])
    else:
        virtual_point = np.mean(ads_site_atoms.positions, axis=0) + 1.5 * n_vector
        distances = np.linalg.norm(ads_site_atoms.positions - virtual_point, axis=1)
        closest_indices = np.argsort(distances)[:2]

        h_vector = (
            ads_site_atoms.positions[closest_indices[0]]
            - ads_site_atoms.positions[closest_indices[1]]
        )
        h_vector /= np.linalg.norm(h_vector)

    return h_vector


def get_shrinkwrap_ads_sites(
    atoms: Atoms,
    precision: float = 0.5,
    touch_sphere_size: float = 3,
    return_trj: bool = False,
):
    """Identifies adsorption sites on a surface using a shrinkwrap grid.

    Args:
        atoms (Atoms): Atoms slab.
        precision (float): Precision for the shrinkwrap grid.
        touch_sphere_size (float): Radius to consider for grid points.
        return_trj (bool): Whether to return the trajectory for demo mode.

    Returns:
        dict: Dictionary containing site information.
        list (optional): Trajectory list if return_trj is True.
    """

    grid = get_shrinkwrap_grid(
        atoms, precision=precision, touch_sphere_size=touch_sphere_size
    )
    surf_ind = shrinkwrap_surface(
        atoms, precision=precision, touch_sphere_size=touch_sphere_size
    )
    targets = get_list_of_touching(atoms, grid, surf_ind, epsilon=0.1)

    trj = []
    coordinates = []
    connectivity = []
    topology = []
    n_vector = []
    h_vector = []
    site_formula = []

    for target in targets:
        atoms_copy = atoms.copy()

        for index in target:
            atoms_copy.append(Atom("X", atoms_copy[index].position + [0, 0, 0]))

        extended_atoms = atoms_copy.copy() * [2, 2, 1]
        extended_grid = grid.copy() * [2, 2, 1]

        if len(target) == 1:
            site_atoms = atoms_copy[target]
            site_coord = site_atoms.positions[0]

        else:
            combs = []
            min_std_devs = []

            for c in itertools.combinations(
                [atom.index for atom in extended_atoms if atom.symbol == "X"],
                len(target),
            ):
                c = list(c)
                min_std_devs.append(max(extended_atoms.positions[c].std(axis=0)))
                combs.append(c)

            min_std_devs = np.array(min_std_devs)
            min_comb_index = np.argmin(min_std_devs)

            site_atoms = extended_atoms[combs[min_comb_index]]
            site_coord = np.mean(site_atoms.positions, axis=0)
            site_coord = get_wrapped_site(site_coord, atoms_copy)
            site_coord = np.array(site_coord)

        n_vec = get_shrinkwrap_site_n_vector(
            extended_atoms, site_coord, extended_grid, touch_sphere_size
        )
        h_vec = get_shrinkwrap_site_h_vector(site_atoms, n_vec)
        site_form = atoms[target].symbols.formula.count()

        coordinates.append(site_coord)
        n_vector.append(n_vec)
        h_vector.append(h_vec)
        topology.append(target)
        connectivity.append(len(target))
        site_formula.append(site_form)

    sites_dict = {
        "coordinates": coordinates,
        "connectivity": connectivity,
        "topology": topology,
        "n_vector": n_vector,
        "h_vector": h_vector,
        "site_formula": site_formula,
    }

    if return_trj:  # demo mode
        extended_atoms = extended_atoms[
            [
                atom.index
                for atom in extended_atoms
                if np.linalg.norm(atom.position - site_coord) < 7
            ]
        ]
        for m in range(20):
            extended_atoms.append(Atom("H", site_coord + n_vec * m * 0.5))
        trj.append(extended_atoms)
        return sites_dict, trj

    return sites_dict


def get_list_of_touching(
    atoms: Atoms,
    grid: Atoms,
    surf_ind: np.ndarray,
    touch_sphere_size: float = 3,
    epsilon: float = 0.3,
    precision: float = 0.5,
) -> list:
    """Identifies groups of atoms in the surface that are within a specified distance from grid points.

    Args:
        atoms (Atoms): Atoms slab.
        grid (Atoms): Shrinkwrap grid of atoms.
        surf_ind (np.ndarray): Indices of surface atoms.
        touch_sphere_size (float): Radius to consider for grid points.
        epsilon (float): Tolerance for distance calculations.
        precision (float): Precision for the shrinkwrap grid.

    Returns:
        list: List of groups of touching atoms.
    """

    # Extend surface indices for periodic boundary conditions
    surf_ind_extended = np.tile(surf_ind, 9)
    surf_ind_extended_plusone = (
        surf_ind_extended + 1
    )  # Avoid 0 index and account for PBC

    surface = atoms[surf_ind].copy()
    large_surface = _special_scale(surface)
    surf_positions = large_surface.positions

    list_of_touching = []

    for g in grid:
        distances = np.linalg.norm(surf_positions - g.position, axis=1)
        mask = (distances > touch_sphere_size - epsilon) & (
            distances < touch_sphere_size + epsilon
        )
        touching_plusone = surf_ind_extended_plusone[mask]
        touching = np.unique(touching_plusone) - 1  # Convert back to original indices
        touching.sort()

        if len(touching) > 0:
            touching_str = "-".join(map(str, touching))
            list_of_touching.append(touching_str)

    list_of_touching = list(set(list_of_touching))

    if "" in list_of_touching:
        list_of_touching.remove("")

    list_of_touching = [list(map(int, x.split("-"))) for x in list_of_touching]
    list_of_touching.sort(key=len)

    return list_of_touching


def _special_scale(atoms: Atoms) -> Atoms:
    """Scales the atoms object to create a larger surface by replicating it in the x and y directions.

    Args:
        atoms (Atoms): Atoms object representing the surface.

    Returns:
        Atoms: A larger Atoms object with the surface replicated in the x and y directions.
    """
    large_atoms = atoms.copy()
    del large_atoms[:]

    shifts = np.array([[i, j] for i in [0, 1, -1] for j in [0, 1, -1]])
    base_positions = atoms.positions
    cell = atoms.cell

    for shift in shifts:
        shifted_positions = base_positions + shift[0] * cell[0] + shift[1] * cell[1]
        shifted_atoms = atoms.copy()
        shifted_atoms.positions = shifted_positions
        large_atoms += shifted_atoms

    return large_atoms


def get_shrinkwrap_grid(
    slab: Atoms,
    precision: float,
    drop_increment: float = 0.1,
    touch_sphere_size: float = 3,
    marker: str = "He",
    raster_speed_boost: bool = False,
) -> Atoms:
    """Generates a shrinkwrap grid for a given slab.

    Args:
        slab (Atoms): Atoms object representing the slab.
        precision (float): Precision for the starting grid.
        drop_increment (float): Increment for dropping grid points.
        touch_sphere_size (float): Radius to consider for grid points.
        marker (str): Marker atom type for the grid.
        raster_speed_boost (bool): Whether to use raster speed boost (experimental feature).

    Returns:
        Atoms: Atoms object representing the shrinkwrap grid.
    """

    from scipy.spatial.distance import cdist

    from .raster_utilities import get_surface_from_rasterized_top_view

    if raster_speed_boost:
        raster_surf_index = get_surface_from_rasterized_top_view(
            slab, pixel_per_angstrom=10
        )
        slab = slab[raster_surf_index]

    grid_positions = _get_starting_grid(slab, precision=precision).positions
    large_slab = get_large_atoms(slab)
    slab_positions = large_slab.positions

    distances_to_grid = cdist(grid_positions, slab_positions).min(axis=1)
    drop_vectors = np.array([[0, 0, drop_increment] for _ in grid_positions])

    while (distances_to_grid > touch_sphere_size).any():
        grid_positions -= (
            drop_vectors * (distances_to_grid > touch_sphere_size)[:, np.newaxis]
        )
        distances_to_grid = cdist(grid_positions, slab_positions).min(axis=1)

        if (distances_to_grid > touch_sphere_size).all() and (
            grid_positions[:, 2] <= 0
        ).all():
            break

    grid = Atoms(
        [marker for _ in grid_positions],
        grid_positions,
        pbc=[True, True, True],
        cell=slab.cell,
    )
    grid = grid[[atom.index for atom in grid if atom.position[2] > 0]]

    return grid


def shrinkwrap_surface(
    slab: Atoms, precision: float = 0.5, touch_sphere_size: float = 3.1
) -> np.ndarray:
    """Identifies surface atoms in a slab using a shrinkwrap grid.

    Args:
        slab (Atoms): Atoms object representing the slab.
        precision (float): Precision for the shrinkwrap grid.
        touch_sphere_size (float): Radius to consider for grid points.

    Returns:
        np.ndarray: Indices of surface atoms.
    """
    grid_positions = get_shrinkwrap_grid(
        slab, precision, touch_sphere_size=touch_sphere_size - 0.2
    ).positions
    slab_positions = slab.positions

    distances_to_grid = cdist(slab_positions, grid_positions).min(axis=1)
    mask = distances_to_grid < touch_sphere_size

    surf_indices = np.where(mask)[0]
    return surf_indices


def shrinkwrap_particle(
    particle: Atoms,
    precision: float,
    touch_sphere_size: float,
    angular_resolution: float,
    return_grid: bool = False,
):
    """Finds surface atoms by rotating the particle and applying shrinkwrap grid.

    Args:
        particle (Atoms): Atoms object representing the particle.
        precision (float): Precision for the shrinkwrap grid.
        touch_sphere_size (float): Radius to consider for grid points.
        angular_resolution (float): Angular resolution for rotations.
        return_grid (bool): Whether to return the grid along with surface indices.

    Returns:
        tuple: Surface indices and optionally the grid.
    """
    surf_indices = set()
    rot_trj = []
    grid = Atoms()

    angles = np.arange(0, 360, angular_resolution)
    angle_x, angle_y = np.meshgrid(angles, angles)
    angle_x = angle_x.flatten()
    angle_y = angle_y.flatten()

    center = np.sum(particle.cell * 0.5, axis=0)

    for ax, ay in zip(angle_x, angle_y):
        ratoms = particle.copy()
        print(
            f"At angle_x: {ax}, angle_y: {ay}. Surface atoms found: {len(surf_indices)}"
        )

        ratoms.rotate("x", ax, center=center)
        ratoms.rotate("y", ay, center=center)
        rot_trj.append(ratoms)

        if return_grid:
            grid.rotate("x", ax, center=center)
            grid.rotate("y", ay, center=center)
            grid += get_shrinkwrap_grid(
                ratoms, precision=precision, touch_sphere_size=touch_sphere_size
            )

        surf_indices.update(
            shrinkwrap_surface(
                ratoms, precision=precision, touch_sphere_size=touch_sphere_size
            )
        )

    if return_grid:
        return list(surf_indices), grid
    else:
        return list(surf_indices)


def assign_fragment(atoms: Atoms, fragment: Atoms) -> tuple:
    """Assigns fragment indices to atoms and fragment.

    Args:
        atoms (Atoms): Atoms object representing the main structure.
        fragment (Atoms): Atoms object representing the fragment to be assigned.

    Returns:
        tuple: Updated atoms and fragment with assigned fragment indices.
    """
    if "fragments" not in atoms.arrays:
        fragment_index = 0
    else:
        fragment_index = max(atoms.arrays["fragments"])

    atoms.arrays["fragments"] = np.full(len(atoms), fragment_index, dtype=int)
    fragment.arrays["fragments"] = np.full(len(fragment), fragment_index + 1, dtype=int)

    return atoms, fragment


def _place_adsorbate(
    slab: Atoms,
    adsorbate: Atoms,
    surf_nvector: np.ndarray,
    site: np.ndarray,
    surf_hvector: np.ndarray = np.array([1, 0, 0]),
) -> Atoms:
    """Places an adsorbate on a slab at a specified site with given surface normal and horizontal vectors.

    Args:
        slab (Atoms): Atoms object representing the slab.
        adsorbate (Atoms): Atoms object representing the adsorbate.
        surf_nvector (np.ndarray): Surface normal vector.
        site (np.ndarray): Site coordinates on the slab.
        surf_hvector (np.ndarray): Surface horizontal vector.

    Returns:
        Atoms: Updated slab with the adsorbate placed.
    """
    from ase.build.tools import sort as ase_sort

    adsorbate_copy = adsorbate.copy()

    # Normalize the surface horizontal vector
    surf_hvector[2] = 0
    surf_hvector /= np.linalg.norm(surf_hvector)

    # Rotate the adsorbate around the z-axis
    rotation_angle = np.sign(surf_hvector[1]) * np.arccos(surf_hvector[0]) / np.pi * 180
    adsorbate_copy.rotate("z", rotation_angle)

    # Position the adsorbate based on its chemical formula
    if adsorbate_copy.get_chemical_formula() == "H":
        adsorbate_copy[0].position = [0, 0, 1.51]
    else:
        adsorbate_copy = align_to_vector(adsorbate_copy, surf_nvector)

    # Remove specific atoms based on their symbols
    if adsorbate_copy[0].symbol == "Cl":
        adsorbate_copy = adsorbate_copy[1:]
    elif adsorbate_copy[0].symbol == "S" and adsorbate_copy[1].symbol == "S":
        adsorbate_copy = adsorbate_copy[2:]

    # Sort the adsorbate atoms
    adsorbate_copy = ase_sort(adsorbate_copy)

    # Translate the adsorbate to the specified site
    adsorbate_copy.positions += site

    # Add the adsorbate to the slab
    slab += adsorbate_copy

    return slab


def check_sites_dict(sites_dict):
    """check if the lenghts of all contents in sites_dict are the same.

    Args:
        sites_dict (_type_): as taken by Surf.NEW_attach_intermediate

    Returns:
        bool: true if all lenghts are the same
    """

    result = True

    # coords_len = len(sites_dict.keys()[0])

    # for key in sites_dict.keys():
    #     if len(sites_dict.keys()) != coords_len:
    #         result = False

    return result


def get_wrapped_site(site: np.ndarray, atoms: Atoms) -> np.ndarray:
    """Wraps a site position within the unit cell of the given atoms.

    Args:
        site (np.ndarray): Site coordinates to be wrapped.
        atoms (Atoms): Atoms object representing the unit cell.

    Returns:
        np.ndarray: Wrapped site coordinates.
    """
    atoms_copy = atoms.copy()
    atoms_copy.append(Atom("X", site))
    atoms_copy.wrap()
    return atoms_copy[-1].position


def get_adsorbate_formula(atoms: Atoms) -> tuple:
    """Gets the chemical formula and count of the adsorbate fragments in the atoms object.

    Args:
        atoms (Atoms): Atoms object containing the fragments.

    Returns:
        tuple: Chemical formula and count of the adsorbate fragments.
    """
    fragment_indices = np.unique(atoms.arrays["fragments"])
    fragment_indices = fragment_indices[
        fragment_indices != 0
    ]  # Exclude the first fragment

    adsorbate_atoms = atoms[
        [
            i
            for i, atom in enumerate(atoms)
            if atoms.arrays["fragments"][i] in fragment_indices
        ]
    ]

    formula = str(adsorbate_atoms.symbols.formula)
    formula_count = adsorbate_atoms.symbols.formula.count()

    return formula, formula_count


def NEW_attach_intermediate(
    atoms=None,
    smiles_list=[],
    sites_dict={
        "coordinates": None,
        "connectivity": None,
        "topology": None,
        "n_vector": None,
        "h_vector": None,
        "site_formula": None,
    },
    configs_per_smiles=1,
    fragment_overlap_thr=1.5,
    try_perturb=10,
    perturb_scale=0.2,
    sample_rotation=False,
    verbose=False,
):
    """
    Takes an atoms object and lists of molecules/intermediates (as ase.atoms objects) and places them on
    all provided top positions (if "Cl") and "bridge" (S1S) positions.

    Args:
        atoms (Atoms, optional): Atoms object representing the main structure. Defaults to None.
        smiles_list (list, optional): List of SMILES strings representing the adsorbates. Defaults to [].
        sites_dict (dict, required): Dictionary containing site information.
        configs_per_smiles (int, optional): Number of configurations per SMILES. Defaults to 1.
        fragment_overlap_thr (float, optional): Overlap threshold for fragments. Defaults to 1.5.
        try_perturb (int, optional): Number of perturbation attempts. Defaults to 10.
        perturb_scale (float, optional): Scale for perturbation. Defaults to 0.2.
        sample_rotation (bool, optional): Whether to sample rotation. Defaults to False.
        verbose (bool, optional): Whether to print verbose output. Defaults to False.

    Returns:
        list: List of updated atoms objects with adsorbates placed.
    """

    fragment_overlap_thr_default = deepcopy(fragment_overlap_thr)
    sd = {
        "coordinates": None,
        "connectivity": None,
        "topology": None,
        "n_vector": None,
        "h_vector": None,
        "site_formula": None,
    }
    sd.update(sites_dict)

    if not check_sites_dict(sd):
        raise ValueError(
            "Site dictionary has inconsistent length. All values in sites_dict shall be lists of equal length."
        )

    coordinates = sd["coordinates"]
    connectivity = sd["connectivity"]
    topology = sd["topology"]
    n_vector = sd["n_vector"]
    h_vector = sd["h_vector"]
    site_formula = sd["site_formula"]

    out_trj = []

    for i, coord in enumerate(coordinates):
        for smile in smiles_list:
            fragment_overlap_thr = fragment_overlap_thr_default

            if smile == "H":
                new_atoms = handle_hydrogen(
                    atoms,
                    coord,
                    i,
                    smile,
                    coordinates,
                    connectivity,
                    topology,
                    n_vector,
                    h_vector,
                    site_formula,
                )
                out_trj.append(new_atoms)
            else:
                handle_other_adsorbates(
                    atoms,
                    smile,
                    i,
                    coord,
                    fragment_overlap_thr,
                    configs_per_smiles,
                    try_perturb,
                    perturb_scale,
                    sample_rotation,
                    verbose,
                    coordinates,
                    connectivity,
                    topology,
                    n_vector,
                    h_vector,
                    site_formula,
                    out_trj,
                )

    return out_trj


def handle_hydrogen(
    atoms: Atoms,
    coord,
    i: int,
    smile: str,
    coordinates,
    connectivity,
    topology,
    n_vector,
    h_vector,
    site_formula,
) -> Atoms:
    """
    Handles the placement of a hydrogen adsorbate on the slab.

    Args:
        atoms (Atoms): The slab atoms object.
        coord: The coordinates of the adsorption site.
        i (int): The index of the current site.
        smile (str): The SMILES string of the adsorbate.
        coordinates: The list of coordinates for all sites.
        connectivity: The list of connectivity for all sites.
        topology: The list of topology for all sites.
        n_vector: The list of normal vectors for all sites.
        h_vector: The list of h vectors for all sites.
        site_formula: The list of formulas for all sites.

    Returns:
        Atoms: The new atoms object with the hydrogen adsorbate.
    """
    adsorbate = Atoms("H", positions=[[0, 0, 1.1]])
    adsorbate.positions += coord
    new_atoms = atoms.copy() + adsorbate
    new_atoms, adsorbate = assign_fragment(new_atoms, adsorbate)
    new_atoms = new_atoms + adsorbate

    info = {
        "smiles": smile,
        "adsorbate_formula": "H",
        "adsorbate_formula_count": 1,
        "ads_site_coord": coordinates[i],
        "ads_site_connectivity": connectivity[i],
        "ads_site_topology": topology[i],
        "ads_site_hvector": h_vector[i],
        "ads_site_nvector": n_vector[i],
        "ads_site_formula": site_formula[i],
    }
    new_atoms.info.update(info)
    return new_atoms


def handle_other_adsorbates(
    atoms: Atoms,
    smile: str,
    i: int,
    coord,
    fragment_overlap_thr: float,
    configs_per_smiles: int,
    try_perturb: int,
    perturb_scale: float,
    sample_rotation: bool,
    verbose: bool,
    coordinates,
    connectivity,
    topology,
    n_vector,
    h_vector,
    site_formula,
    out_trj: list,
) -> None:
    """
    Handles the placement of non-hydrogen adsorbates on the slab.

    Args:
        atoms (Atoms): The slab atoms object.
        smile (str): The SMILES string of the adsorbate.
        i (int): The index of the current site.
        coord: The coordinates of the adsorption site.
        fragment_overlap_thr (float): The threshold for fragment overlap.
        configs_per_smiles (int): Number of conformations per SMILES string.
        try_perturb (int): The number of perturbation attempts.
        perturb_scale (float): The scale of perturbation.
        sample_rotation (bool): Whether to sample rotations.
        verbose (bool): Whether to print verbose output.
        coordinates: The list of coordinates for all sites.
        connectivity: The list of connectivity for all sites.
        topology: The list of topology for all sites.
        n_vector: The list of normal vectors for all sites.
        h_vector: The list of h vectors for all sites.
        site_formula: The list of formulas for all sites.
        out_trj (list): The list to append successful perturbations to.
    """
    from . import Smile

    if smile == "Cl":
        fragment_overlap_thr = 1.0

    adsorbates = Smile.NEW_create_adsorbates(
        [smile], conformations_per_smiles=configs_per_smiles
    )

    for adsorbate in adsorbates:
        if (adsorbate[0].symbol == "Cl" and connectivity[i] == 1) or (
            adsorbate[0].symbol == "S" and connectivity[i] > 1
        ):
            if verbose:
                print(f"Trying smile: {smile}")

            new_slab = atoms.copy()
            new_slab, adsorbate = assign_fragment(new_slab, adsorbate)
            new_slab = _place_adsorbate(
                new_slab, adsorbate, n_vector[i], coord, surf_hvector=h_vector[i]
            )

            mfd = minimum_fragment_distance(new_slab)

            if mfd > fragment_overlap_thr:
                _update_atoms_info(
                    new_slab,
                    smile,
                    coordinates[i],
                    connectivity[i],
                    topology[i],
                    n_vector[i],
                    h_vector[i],
                    site_formula[i],
                )
                out_trj.append(new_slab)
                if verbose:
                    print(f"SUCCESS at smile: {smile}")
            else:
                if verbose:
                    print(f"    Fragment too close: {mfd}. ")
                _try_perturb_adsorbate(
                    atoms,
                    adsorbate,
                    fragment_overlap_thr,
                    try_perturb,
                    perturb_scale,
                    sample_rotation,
                    verbose,
                    n_vector[i],
                    coord,
                    h_vector[i],
                    smile,
                    coordinates[i],
                    connectivity[i],
                    topology[i],
                    site_formula[i],
                    out_trj,
                )


def _try_perturb_adsorbate(
    atoms,
    adsorbate,
    fragment_overlap_thr: float,
    try_perturb: int,
    perturb_scale: float,
    sample_rotation: bool,
    verbose: bool,
    n_vector,
    coord,
    h_vector,
    smile: str,
    site_coord,
    site_connectivity,
    site_topology,
    site_formula,
    out_trj: list,
) -> None:
    """
    Attempts to perturb the adsorbate and place it on the slab until the minimum fragment distance is greater than the threshold.

    Args:
        atoms (Atoms): The slab atoms object.
        adsorbate (Atoms): The adsorbate atoms object.
        fragment_overlap_thr (float): The threshold for fragment overlap.
        try_perturb (int): The number of perturbation attempts.
        perturb_scale (float): The scale of perturbation.
        sample_rotation (bool): Whether to sample rotations.
        verbose (bool): Whether to print verbose output.
        n_vector: The normal vector at the adsorption site.
        coord: The coordinates of the adsorption site.
        h_vector: The h vector at the adsorption site.
        smile (str): The SMILES string of the adsorbate.
        site_coord: The coordinates of the adsorption site.
        site_connectivity: The connectivity of the adsorption site.
        site_topology: The topology of the adsorption site.
        site_formula: The formula of the adsorption site.
        out_trj (list): The list to append successful perturbations to.
    """
    counter = 0
    rot_counter = 1

    while counter < try_perturb:
        new_slab = atoms.copy()
        perturbed_adsorbate = adsorbate.copy()

        if not sample_rotation:
            perturbed_adsorbate = perturb_adsorbate(
                adsorbate, scale=perturb_scale, force_rotation=False
            )
        else:
            perturbed_adsorbate = perturb_adsorbate(
                adsorbate, scale=perturb_scale, force_rotation=rot_counter
            )
            rot_counter += 1

        new_slab, perturbed_adsorbate = assign_fragment(new_slab, perturbed_adsorbate)
        new_slab = _place_adsorbate(
            new_slab, perturbed_adsorbate, n_vector, coord, surf_hvector=h_vector
        )

        mfd = minimum_fragment_distance(new_slab)

        if verbose:
            print(f"Tried fragment perturb {counter} times. Fragment too close: {mfd}.")

        if mfd > fragment_overlap_thr:
            _update_atoms_info(
                new_slab,
                smile,
                site_coord,
                site_connectivity,
                site_topology,
                n_vector,
                h_vector,
                site_formula,
            )
            out_trj.append(new_slab)
            break

        counter += 1


def _update_atoms_info(
    atoms,
    smile: str,
    site_coord,
    site_connectivity,
    site_topology,
    n_vector,
    h_vector,
    site_formula,
) -> None:
    """
    Updates the info dictionary of an Atoms object with adsorbate and site information.

    Args:
        atoms (Atoms): The Atoms object to update.
        smile (str): The SMILES string of the adsorbate.
        site_coord: The coordinates of the adsorption site.
        site_connectivity: The connectivity of the adsorption site.
        site_topology: The topology of the adsorption site.
        n_vector: The normal vector at the adsorption site.
        h_vector: The h vector at the adsorption site.
        site_formula: The formula of the adsorption site.
    """
    adsorbate_formula, adsorbate_formula_count = get_adsorbate_formula(atoms)

    atoms_info = {
        "smiles": smile,
        "adsorbate_formula": adsorbate_formula,
        "adsorbate_formula_count": adsorbate_formula_count,
        "ads_site_coord": site_coord,
        "ads_site_connectivity": site_connectivity,
        "ads_site_topology": site_topology,
        "ads_site_hvector": h_vector,
        "ads_site_nvector": n_vector,
        "ads_site_formula": site_formula,
    }

    atoms.info.update(atoms_info)


def perturb_adsorbate(
    adsorbate: Atoms, scale: float = 0.2, force_rotation: bool = False
) -> Atoms:
    """
    Perturbs the adsorbate by applying random translation and rotation.

    Args:
        adsorbate (Atoms): The adsorbate atoms object.
        scale (float): The scale of perturbation.
        force_rotation (bool): Whether to force a specific rotation.

    Returns:
        Atoms: The perturbed adsorbate atoms object.
    """
    rand_vector = np.random.normal(loc=0.0, scale=scale, size=3)

    if adsorbate[0].symbol == "Cl":
        rand_angle = np.random.random() * 360

        if force_rotation:
            rand_angle = force_rotation * 30

        rand_vector[2] = 1
        rand_vector /= np.linalg.norm(rand_vector)

        adsorbate.rotate("z", rand_angle)
        adsorbate = align_to_vector(adsorbate, rand_vector)

    elif adsorbate[0].symbol == "S":
        rand_angle = 0
        if force_rotation:
            rand_angle = force_rotation * 180

        rand_vector[0] = 0
        rand_vector[2] = 1
        rand_vector /= np.linalg.norm(rand_vector)

        adsorbate.rotate("z", rand_angle)
        adsorbate = align_to_vector(adsorbate, rand_vector)

    else:
        raise ValueError(
            "Adsorbate must originate from a marked smiles (starting with Cl or S)"
        )

    return adsorbate


def minimum_fragment_distance(atoms: Atoms, compare_inds: list = []) -> float:
    """
    Calculates the minimum distance between fragments in the given Atoms object.

    Args:
        atoms (Atoms): The Atoms object containing the fragments.

    Returns:
        float: The minimum distance between any two fragments.
    """
    fragments = list(set(atoms.arrays["fragments"]))

    fragment_atoms_dict = {}
    for f in fragments:
        fragment_atoms = Atoms()
        for i, a in enumerate(atoms):
            if atoms.arrays["fragments"][i] == f:
                fragment_atoms.append(atoms[i])
        fragment_atoms.cell = atoms.cell
        fragment_atoms_dict[f] = fragment_atoms * [2, 2, 2]

    if len(compare_inds) == 0:
        pairs = [(a, b) for i, a in enumerate(fragments) for b in fragments[i + 1 :]]
    else:
        pairs = [x for x in itertools.combinations(compare_inds, 2)]

    min_distances = []
    for p in pairs:
        positions_a = fragment_atoms_dict[p[0]].positions
        positions_b = fragment_atoms_dict[p[1]].positions

        min_dist = cdist(positions_a, positions_b).min()
        min_distances.append(min_dist)

    return min(min_distances) if min_distances else np.inf


def _trilaterate(
    P1: np.ndarray, P2: np.ndarray, P3: np.ndarray, r1: float, r2: float, r3: float
) -> tuple:
    """
    Find the intersection of three spheres.

    Args:
        P1 (np.ndarray): The center of the first sphere.
        P2 (np.ndarray): The center of the second sphere.
        P3 (np.ndarray): The center of the third sphere.
        r1 (float): The radius of the first sphere.
        r2 (float): The radius of the second sphere.
        r3 (float): The radius of the third sphere.

    Returns:
        tuple: The two intersection points of the three spheres.

    Raises:
        ValueError: If the three spheres do not intersect.
    """
    from numpy import cross, dot, sqrt
    from numpy.linalg import norm

    temp1 = P2 - P1
    e_x = temp1 / norm(temp1)
    temp2 = P3 - P1
    i = dot(e_x, temp2)
    temp3 = temp2 - i * e_x
    e_y = temp3 / norm(temp3)
    e_z = cross(e_x, e_y)
    d = norm(P2 - P1)
    j = dot(e_y, temp2)
    x = (r1**2 - r2**2 + d**2) / (2 * d)
    y = (r1**2 - r3**2 - 2 * i * x + i**2 + j**2) / (2 * j)
    temp4 = r1**2 - x**2 - y**2

    if temp4 < 0:
        raise ValueError("The three spheres do not intersect!")

    z = sqrt(temp4)
    p_12_a = P1 + x * e_x + y * e_y + z * e_z
    p_12_b = P1 + x * e_x + y * e_y - z * e_z

    return p_12_a, p_12_b


def _get_pbc_mean(atoms: Atoms, targets: list, marker: str = "X") -> np.ndarray:
    """
    Calculates the periodic boundary condition (PBC) mean position of the target atoms.

    Args:
        atoms (Atoms): The Atoms object containing the atoms.
        targets (list): The list of target atom indices.
        marker (str): The marker symbol to use for temporary atoms.

    Returns:
        np.ndarray: The mean position of the target atoms under PBC.

    Raises:
        ValueError: If the marker is already in use in the atoms object.
    """
    if marker in atoms.get_chemical_symbols():
        raise ValueError(
            f"Marker: {marker} in use in atoms object. Please use a different marker."
        )

    catoms = atoms.copy()

    for t in targets:
        catoms.append(Atom(marker, catoms[t].position + [0, 0, 1]))

    catoms = catoms * [2, 2, 1]

    xlist = [atom.index for atom in catoms if atom.symbol == marker]

    dmin = float("inf")
    cmin = None
    for c in itertools.combinations(xlist, len(targets)):
        plist = catoms[c].positions
        vlist = [
            plist[pair[0]] - plist[pair[1]]
            for pair in itertools.combinations(range(len(targets)), 2)
        ]
        d = sum(np.linalg.norm(v) for v in vlist)
        if d < dmin:
            dmin = d
            cmin = c

    return catoms[cmin].get_center_of_mass()


def get_nvector(atoms: Atoms) -> np.ndarray:
    """
    Calculates the normal vector for the given atoms.

    Args:
        atoms (Atoms): The Atoms object containing the atoms.

    Returns:
        np.ndarray: The normal vector.
    """
    nvector = np.array([0, 0, 1])  # Default normal vector

    if len(atoms) == 1:
        return atoms[0].position
    atom1 = atoms[0]  # First surrogate atom
    atom2 = atoms[1]  # Second atom

    if atom1.symbol == "Cl":
        nvector = atom2.position - atom1.position
    elif atom1.symbol == "S" and atom2.symbol == "S":
        atom3 = atoms[
            2
        ]  # Atom directly attached to the site (atom following S1S in SMILES)
        # Compute normal vector to S-S "defect marker"
        nvector = np.cross(
            np.cross(atom1.position - atom2.position, atom2.position - atom3.position),
            atom1.position - atom2.position,
        )
        nvector = nvector / -np.linalg.norm(nvector)

    return nvector


def align_to_vector(atoms: Atoms, vector: list = [0, 0, 1]) -> Atoms:
    """
    Aligns the atoms to a given vector by applying a rotation.

    Args:
        atoms (Atoms): The Atoms object to be aligned.
        vector (list): The target vector to align the atoms to. Default is [0, 0, 1].

    Returns:
        Atoms: The aligned Atoms object.
    """
    atomsr = atoms.copy()
    vec1 = get_nvector(atomsr)
    vec2 = np.array(vector)
    mat = rotation_matrix_from_vectors(vec1, vec2)

    for atom in atomsr:
        vec = atom.position
        vec_rot = mat.dot(vec)
        atom.position = vec_rot

    return atomsr


def rotate_mol_to_hvector(atoms: Atoms, site: int, mol: Atoms) -> Atoms:
    """
    Rotates the molecule to align with the h-vector of the given site.

    Args:
        atoms (Atoms): The Atoms object containing the slab.
        site (int): The index of the site to align the molecule to.
        mol (Atoms): The molecule to be rotated.

    Returns:
        Atoms: The rotated molecule.
    """
    m = mol.copy()
    hvector = get_hvector(atoms, site)
    angle = (np.arccos(hvector[0])) / np.pi * 180  # Convert radians to degrees
    m.rotate("z", angle)
    return m


def get_hvector(atoms: Atoms, site: np.ndarray) -> np.ndarray:
    """
    Calculates the h-vector for the given site in the atoms object.

    Args:
        atoms (Atoms): The Atoms object containing the atoms.
        site (np.ndarray): The coordinates of the site.

    Returns:
        np.ndarray: The h-vector.
    """
    distances = {i: np.linalg.norm(pos - site) for i, pos in enumerate(atoms.positions)}
    sorted_distances = sorted(distances.items(), key=lambda item: item[1])
    id_list = [pair[0] for pair in sorted_distances[:2]]
    dvec = atoms.positions[id_list[0]] - atoms.positions[id_list[1]]
    dvec[2] = 0  # Project onto the xy-plane
    hvector = dvec / np.linalg.norm(dvec)
    return hvector


def get_empty_atoms(atoms: Atoms) -> Atoms:
    """
    Returns a copy of the given ASE Atoms object with all atoms removed.

    Args:
        atoms (Atoms): The ASE Atoms object to be emptied.

    Returns:
        Atoms: A copy of the ASE Atoms object with all atoms removed.
    """
    empty_atoms = atoms.copy()
    del empty_atoms[:]
    return empty_atoms


def get_large_atoms(atoms: Atoms, gutter: float = 2) -> Atoms:
    """
    Returns a copy of the given ASE Atoms object, expanded and filtered based on the gutter value.

    Args:
        atoms (Atoms): The ASE Atoms object to be expanded.
        gutter (float): The gutter value used for filtering atoms. Default is 2.

    Returns:
        Atoms: A copy of the ASE Atoms object, expanded and filtered.
    """
    gutter = min(gutter, 2)

    large_atoms = atoms.copy() * [3, 3, 1]

    del_list = []
    for i, p in enumerate(large_atoms.get_scaled_positions(wrap=True)):
        if any(
            [
                p[0] < 1 / 3 - gutter / np.linalg.norm(atoms.cell[0]),
                p[0] > 2 / 3 + gutter / np.linalg.norm(atoms.cell[0]),
                p[1] < 1 / 3 - gutter / np.linalg.norm(atoms.cell[1]),
                p[1] > 2 / 3 + gutter / np.linalg.norm(atoms.cell[1]),
            ]
        ):
            del_list.append(i)
    del large_atoms[del_list]

    for a in large_atoms:
        a.position -= atoms.cell[0] + atoms.cell[1]

    return large_atoms


def _check_a_b_vectors(atoms: Atoms) -> bool:
    """
    Checks if the cell follows the common convention (as per ASE):
    1) atoms.cell[0] vector is collinear with the Cartesian x-axis.
    2) atoms.cell[1] vector is heading towards the positive y-axis.
    3) atoms.cell[0] and atoms.cell[1] vectors belong to the Cartesian x,y plane.

    Args:
        atoms (Atoms): The ASE Atoms object to be checked.

    Returns:
        bool: True if the checks passed, False otherwise.
    """
    cell = atoms.cell
    check = all(
        [
            cell[0][1] == 0,
            cell[0][2] == 0,
            cell[1][1] > 0,
            np.isclose(np.dot(np.cross(cell[0], cell[1]), [1, 0, 0]), 0),
        ]
    )
    return check


def _get_starting_grid(atoms: Atoms, precision: float, marker: str = "He") -> Atoms:
    """
    Generates a starting grid of marker atoms above the given atoms object.

    Args:
        atoms (Atoms): The ASE Atoms object to be used as a reference.
        precision (float): The precision for the grid spacing.
        marker (str): The marker symbol to use for the grid atoms. Default is 'He'.

    Returns:
        Atoms: An ASE Atoms object containing the grid of marker atoms.
    """
    if _check_a_b_vectors(atoms):
        range_x = np.arange(0, atoms.cell[0][0] + atoms.cell[1][0], precision)
        if atoms.cell[1][0] < 0:
            range_x = np.arange(
                atoms.cell[1][0], atoms.cell[0][0] + np.abs(atoms.cell[1][0]), precision
            )
        range_y = np.arange(0, atoms.cell[1][1], precision)
    else:
        max_range_limit_a = max(
            [
                np.linalg.norm(atoms.cell[0] + atoms.cell[1]),
                np.linalg.norm(atoms.cell[0] - atoms.cell[1]),
                np.linalg.norm(-atoms.cell[0] + atoms.cell[1]),
                np.linalg.norm(-atoms.cell[0] - atoms.cell[1]),
            ]
        )

        max_range_limit_b = max(
            [
                np.linalg.norm(atoms.cell[0] + atoms.cell[1]),
                np.linalg.norm(atoms.cell[0] - atoms.cell[1]),
                np.linalg.norm(-atoms.cell[0] + atoms.cell[1]),
                np.linalg.norm(-atoms.cell[0] - atoms.cell[1]),
            ]
        )

        range_x = np.arange(0, max_range_limit_a, precision)
        range_y = np.arange(0, max_range_limit_b, precision)

    grid = get_empty_atoms(atoms)
    start_grid_z = max([a.position[2] for a in atoms]) + 3.5

    for x, y in [(x, y) for x in range_x for y in range_y]:
        grid += Atom(marker, [x, y, start_grid_z])
        if any(grid.get_scaled_positions(wrap=False)[-1] > 1.01) or any(
            grid.get_scaled_positions(wrap=False)[-1] < -0.01
        ):
            del grid[-1]
    return grid


def _drop_marker(
    grid: Atoms, atoms: Atoms, marker_index: int, touch_sphere_size: float = 3
) -> Atoms:
    """
    Helper function that lowers each grid marker atom along the Cartesian z-axis until the distance criteria is reached.

    Args:
        grid (Atoms): The initial grid of marker atoms.
        atoms (Atoms): The atoms object onto which the grid is dropped.
        marker_index (int): The index of the marker atom in the grid.
        touch_sphere_size (float): The size of the touch sphere. Default is 3.

    Returns:
        Atoms: The updated grid with the marker atom dropped.
    """
    del_list = []
    a = grid[marker_index]
    dist = 10

    while dist > 0.1:
        dist = (
            min([np.linalg.norm(a.position - b.position) for b in atoms])
            - touch_sphere_size
        )
        a.position[2] -= dist * 0.5

        if a.position[2] < 0.1:
            del_list.append(marker_index)
            break

    if del_list:
        del grid[del_list[0]]

    return grid


def get_AFM_cartoon(atoms: Atoms, precision: float = 1, show_figure=False) -> None:
    """
    Generates an AFM cartoon visualization of the given atoms object.

    Args:
        atoms (Atoms): The ASE Atoms object to be visualized.
        precision (float): The precision for the grid spacing. Default is 1.

    Returns:
        If show_figure is False, returns a DataFrame with z values.
    """
    grid = get_shrinkwrap_grid(atoms, precision)
    min_g = min([g[2] for g in grid.positions])

    grid.positions -= [0, 0, min_g]

    df = []
    for g in grid:
        df.append({"x": g.position[0], "y": g.position[1], "z": g.position[2]})
    df = pd.DataFrame(df)
    df["x"] = df["x"].round(3)
    df["y"] = df["y"].round(3)

    df = df.pivot(index="y", columns="x", values="z")
    df = df.reindex(index=df.index[::-1])

    if show_figure:
        import matplotlib.pyplot as plt
        import sns

        sns.heatmap(df)
        plt.axis("scaled")
        plt.show()
    else:
        return df


def to_bottom_of_cell(atoms: Atoms) -> Atoms:
    """
    Drops atoms so that the lowest (z) position is at z=0.

    Args:
        atoms (Atoms): The ASE Atoms object to be shifted.

    Returns:
        Atoms: The shifted ASE Atoms object.
    """
    a = atoms.copy()

    min_z = min(atom.position[2] for atom in a)
    for atom in a:
        atom.position[2] -= min_z
    return a


def freeze_atoms(atoms: Atoms, fract_from_bottom: float = 0.66) -> Atoms:
    """
    Freezes a fraction of atoms from the bottom of the structure.

    Args:
        atoms (Atoms): The ASE Atoms object to freeze.
        fract_from_bottom (float, optional): The fraction of atoms to freeze from the bottom. Defaults to 0.66.

    Returns:
        Atoms: The ASE Atoms object with constraints set.
    """
    from ase.constraints import FixAtoms

    a = atoms.copy()

    min_z = min(atom.position[2] for atom in a)
    max_z = max(atom.position[2] for atom in a)

    freeze_z = min_z + (max_z - min_z) * fract_from_bottom

    c = FixAtoms(indices=[atom.index for atom in a if atom.position[2] < freeze_z])
    a.set_constraint(c)

    return a


def attach_fragment(
    atoms: Atoms,
    site_dict: dict,
    fragment: Atoms,
    n_rotation: float,
    height: float = None,
) -> Atoms:
    """
    Places a fragment on a slab in a simple way.

    Args:
        atoms (Atoms): The slab ASE Atoms object.
        site_dict (dict): Dictionary containing site information.
            Example:
                {
                    'coordinates': array([1.29400541, 2.24128311, 22.33930596]),
                    'connectivity': 1,
                    'topology': [30],
                    'n_vector': array([0.00752117, -0.00320119, 0.99996659]),
                    'h_vector': array([1., 0., 0.]),
                    'site_formula': {'Cu': 1}
                }
        fragment (Atoms): The fragment to be attached.
        n_rotation (float): Angle to rotate around n_vector.
        height (float, optional): Height to place the fragment. Defaults to None.

    Returns:
        Atoms: The slab with the fragment attached.
    """
    fragment_copy = fragment.copy()

    drop_fragment_to_z_zero(fragment_copy)
    if height is not None:
        fragment_copy.positions += [0, 0, height]

    if fragment_copy[0].symbol == "S" and fragment_copy[1].symbol == "S":
        shift_x = fragment_copy[:2].get_center_of_mass()[0]
        fragment_copy.positions -= np.array([shift_x, 0, 0])

    fragment_copy.rotate(n_rotation, "z")
    fragment_copy = align_to_vector(fragment_copy, site_dict["n_vector"])

    if fragment_copy[0].symbol == "Cl":
        del fragment_copy[:1]
    if fragment_copy[0].symbol == "S" and fragment_copy[1].symbol == "S":
        del fragment_copy[:2]

    fragment_copy.positions += site_dict["coordinates"]

    atoms += fragment_copy

    return atoms


def drop_fragment_to_z_zero(atoms: Atoms) -> None:
    """
    Drops the fragment so that the lowest (z) position is at z=0.

    Args:
        atoms (Atoms): The ASE Atoms object to be shifted.
    """
    if atoms[0].symbol == "Cl":
        min_z = np.min(atoms.positions[1:, 2])
    elif atoms[0].symbol == "S" and atoms[1].symbol == "S":
        min_z = np.min(atoms.positions[2:, 2])
    else:
        raise ValueError(
            "Unsupported fragment type. The fragment must start with 'Cl' or 'S'."
        )

    atoms.positions -= np.array([0, 0, min_z])


def conformer_to_site(atoms, site, conformer, mode="optimize", overlap_thr=0):
    """
    Aligns and attaches a conformer to a specified site on a slab of atoms, optimizing the orientation to minimize overlap.

    Parameters:
    atoms (object): An object containing an array of atoms.
    site (dict): A dictionary containing site information, including coordinates, n_vector, and h_vector.
    conformer (object): An object containing the conformer to be attached.
    mode (str): The mode of operation. Currently, only 'optimize' is supported. Default is 'optimize'.
    overlap_thr (float): The overlap threshold. Default is 0.0.

    Returns:
    object: The combined atoms object with the conformer attached and optimized.

    Raises:
    ValueError: If the mode is not implemented or if the atoms object is invalid.
    """
    atoms = atoms.copy()
    conformer = conformer.copy()

    if "fragments" not in atoms.arrays.keys():
        atoms.arrays["fragments"] = np.array([0 for a in atoms])

    n_f = int(np.max(atoms.arrays["fragments"])) + 1
    conformer.arrays["fragments"] = np.array([n_f for a in conformer])

    if conformer.info["smiles"][:3] == "S1S":
        conformer.positions -= conformer[:2].get_center_of_mass()
        conformer = conformer[2:]
    if conformer.info["smiles"][:2] == "Cl":
        conformer = conformer[1:]

    conformer = align_to_vector(conformer, site["n_vector"])

    h_rot_angle = (
        np.sign(site["h_vector"][1]) * np.arccos(site["h_vector"][0]) / np.pi * 180
    )
    conformer.rotate(h_rot_angle, site["n_vector"])
    conformer.positions += site["coordinates"]

    out_atoms = atoms + conformer
    # out_atoms.info['site_info'] = site

    if conformer.info["smiles"][:3] == "S1S":
        out_atoms = swing_fragment(
            atoms=out_atoms,
            fragment_index=n_f,
            site=site,
            resolution=10,
            mode=mode,
            span_angle=50,
            overlap_thr=overlap_thr,
        )

    if conformer.info["smiles"][:2] == "Cl":
        out_atoms = swirl_fragment(
            atoms=out_atoms,
            fragment_index=n_f,
            site=site,
            resolution=10,
            mode="optimize",
            overlap_thr=0.0,
        )

    return out_atoms


def swing_fragment(
    atoms,
    fragment_index,
    site,
    resolution=10,
    mode="optimize",
    span_angle=50,
    overlap_thr=0,
):
    """
    Optimizes the orientation of a fragment relative to a slab to minimize overlap and distance.

    Parameters:
    atoms (object): An object containing an array of atoms.
    fragment_index (int): The index of the target fragment to be separated.
    site (dict): A dictionary containing site information, including coordinates and h_vector.
    resolution (int): The resolution for angle generation in degrees. Default is 10.
    mode (str): The mode of operation. Currently, only 'optimize' is supported. Default is 'optimize'.
    span_angle (int): The span angle for generating swings. Default is 50.
    overlap_thr (float): The overlap threshold. Default is 0.0.

    Returns:
    list: A list containing the optimized atoms object if successful.

    Raises:
    ValueError: If the mode is not implemented or if the atoms object is invalid.
    """
    slab, fragment = split_slab_from_target_fragment(atoms, fragment_index)

    swings = [r for r in range(0, span_angle, resolution)]
    swings = [sw[0] * sw[1] for sw in itertools.product(swings, [-1, 1])]

    mdfs = []
    view_atoms = []
    if mode.lower() == "optimize":
        for swing in swings:
            f = fragment.copy()
            f.rotate(
                swing, site["h_vector"], center=site["coordinates"], rotate_cell=False
            )
            xatoms = slab.copy() + f

            mdf = minimum_fragment_distance(xatoms)
            mdfs.append((xatoms, mdf))
            view_atoms.append(xatoms)

        mdfs = sorted(mdfs, key=lambda tup: tup[1])
        out_atoms, mdf = mdfs[-1]
        out_atoms.info["mdf"] = mdf

        if mdf > overlap_thr:
            return [out_atoms]
        else:
            print(
                f"Failed to generate for: overlap_thr={overlap_thr}, maximal_distance={mdfs[-1][-1]}"
            )

    if mode.lower() != "optimize":
        raise ValueError(
            f"Mode {mode} not yet implemented... Please set mode='optimize'."
        )


def swirl_fragment(
    atoms, fragment_index, site, resolution=10, mode="optimize", overlap_thr=0.0
):
    """
    Optimizes the orientation of a fragment relative to a slab to minimize overlap and distance.

    Parameters:
    atoms (object): An object containing an array of atoms.
    fragment_index (int): The index of the target fragment to be separated.
    site (dict): A dictionary containing site information, including coordinates.
    resolution (int): The resolution for angle generation in degrees. Default is 10.
    mode (str): The mode of operation. Currently, only 'optimize' is supported. Default is 'optimize'.
    overlap_thr (float): The overlap threshold. Default is 0.0.

    Returns:
    list: A list containing the optimized atoms object if successful.

    Raises:
    ValueError: If the mode is not implemented or if the atoms object is invalid.
    """
    slab, fragment = split_slab_from_target_fragment(atoms, fragment_index)

    mdfs = []
    view_atoms = []

    if mode.lower() == "optimize":
        thetas = [math.radians(a) for a in range(5, 85, resolution)]
        phis = [math.radians(a) for a in range(0, 360, resolution * 4)]

        swirls = itertools.product(thetas, phis)
        for swirl in swirls:
            f = fragment.copy()

            f.positions -= np.array(site["coordinates"])
            f = align_to_vector(f, vector=polar2cart(swirl[0], swirl[1]))
            f.positions += np.array(site["coordinates"])

            xatoms = slab.copy() + f

            mdf = minimum_fragment_distance(xatoms)
            mdfs.append((xatoms, mdf))
            view_atoms.append(xatoms)

        mdfs = sorted(mdfs, key=lambda tup: tup[1])
        out_atoms, mdf = mdfs[-1]
        out_atoms.info["mdf"] = mdf

        if mdf > overlap_thr:
            return [out_atoms]
        else:
            print(
                f"Failed to generate for: overlap_thr={overlap_thr}, maximal_distance={mdfs[-1][-1]}"
            )

    if mode.lower() != "optimize":
        raise ValueError(
            f"Mode {mode} not yet implemented... Please set mode='optimize'."
        )


def split_slab_from_target_fragment(atoms, fragment_index):
    """
    Splits the atoms object into two parts: slab and fragment based on the fragment index.

    Parameters:
    atoms (object): An object containing an array of atoms with an 'arrays' attribute that includes 'fragments'.
    fragment_index (int): The index of the target fragment to be separated.

    Returns:
    tuple: A tuple containing two elements:
        - slab (object): The atoms object excluding the target fragment.
        - fragment (object): The atoms object containing only the target fragment.

    Raises:
    ValueError: If the 'atoms' object does not have an 'arrays' attribute with 'fragments' key.
    """

    if not hasattr(atoms, "arrays") or "fragments" not in atoms.arrays:
        raise ValueError(
            "The 'atoms' object must have an 'arrays' attribute with 'fragments' key."
        )

    slab_indices = []
    fragment_indices = []

    for i, atom in enumerate(atoms):
        if atoms.arrays["fragments"][i] == fragment_index:
            fragment_indices.append(i)
        else:
            slab_indices.append(i)

    slab = atoms[slab_indices]
    fragment = atoms[fragment_indices]

    return slab, fragment
