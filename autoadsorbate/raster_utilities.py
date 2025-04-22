import numpy as np
from ase.build.tools import sort as ase_sort

# import matplotlib.image as mpimg
# from ase.io import read, write
# from ase.visualize import view
# import matplotlib.pyplot as plt


def snap_to_grid(atoms):
    """
    Returns the rounded positions of atoms in the given atomic structure.

    Parameters:
    atoms (object): An ASE atoms object representing the atomic structure.

    Returns:
    numpy.ndarray: A 2D array of rounded atomic positions.
    """
    a = atoms.copy()
    pos = np.round(a.positions, 2)
    return pos


def get_pixel_positions(atoms):
    """
    Returns the pixel positions of atoms in the given atomic structure by snapping them to a grid.

    Parameters:
    atoms (object): An ASE atoms object representing the atomic structure.

    Returns:
    numpy.ndarray: A 2D array of pixel positions of the atoms, excluding the last column.
    """
    return snap_to_grid(atoms)[:, :-1]


def get_pixel_grid(atoms, pixel_per_angstrom=100):
    """
    Generates a pixel grid based on the atomic structure and specified pixel density.

    Parameters:
    atoms (object): An ASE atoms object representing the atomic structure.
    pixel_per_angstrom (int): The number of pixels per angstrom. Default is 100.

    Returns:
    numpy.ndarray: A 2D array representing the pixel grid.
    """
    x_footprint = atoms.cell[0][0] + atoms.cell[1][0]
    y_footprint = atoms.cell[0][1] + atoms.cell[1][1]

    x_size = int(np.ceil((x_footprint) * pixel_per_angstrom))
    y_size = int(np.ceil((y_footprint) * pixel_per_angstrom))

    grid = np.zeros([y_size, x_size])

    return grid


def createKernel(radius, value):
    """
    Creates a circular kernel with the specified radius and value.

    Parameters:
    radius (int): The radius of the circular kernel.
    value (float): The value to assign to the elements within the circular area.

    Returns:
    numpy.ndarray: A 2D array representing the circular kernel.
    """
    kernel = np.zeros((2 * radius + 1, 2 * radius + 1))
    y, x = np.ogrid[-radius : radius + 1, -radius : radius + 1]
    mask = x**2 + y**2 <= radius**2
    kernel[mask] = value
    return kernel


def get_surface_from_rasterized_top_view(
    atoms_org, pixel_per_angstrom=10, return_raster=False
):
    """
    Identifies surface atoms from a rasterized top view of the atomic structure.

    Parameters:
    atoms_org (object): An ASE atoms object representing the original atomic structure.
    pixel_per_angstrom (int): The number of pixels per angstrom for the rasterization. Default is 10.
    return_raster (bool): Whether to return the rasterized grid along with the surface indices. Default is False.

    Returns:
    list or tuple: If return_raster is False, returns a list of surface atom indices.
                   If return_raster is True, returns a tuple containing the list of surface atom indices and the rasterized grid.
    """
    atoms = atoms_org.copy()
    atoms.arrays["original_index"] = np.array([atom.index for atom in atoms])
    atoms = atoms * [2, 2, 1]
    atoms = ase_sort(atoms, tags=atoms.positions[:, 2])
    pixel_grid = get_pixel_grid(atoms, pixel_per_angstrom=pixel_per_angstrom)
    h, w = pixel_grid.shape[1], pixel_grid.shape[0]  # img size
    A = pixel_grid.copy()

    mapping = {}

    for i, pos in enumerate(get_pixel_positions(atoms)):
        y, x = pos[1] * pixel_per_angstrom, pos[0] * pixel_per_angstrom
        xx, yy = np.meshgrid(np.linspace(0, h - 1, h), np.linspace(0, w - 1, w))
        radius = pixel_per_angstrom * 2
        mask = (xx - x) ** 2 + (yy - y) ** 2 < radius**2
        A[mask] = i
        mapping[i] = atoms.arrays["original_index"][i]

    B = A[
        int(len(A) * 0.2) : int(len(A) * 0.8), int(len(A.T) * 0.2) : int(len(A.T) * 0.8)
    ]
    surf_inds = list(set([mapping[val] for val in np.unique(B)]))

    if return_raster:
        return surf_inds, B
    else:
        return surf_inds
