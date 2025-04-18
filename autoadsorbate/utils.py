import itertools

import ase
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ase import Atom
from ase.build.tools import sort as sort_atoms
from ase.constraints import FixAtoms
from ase.io import read
from ase.optimize import BFGS
from ase.optimize.minimahopping import MinimaHopping
from ase.visualize.plot import plot_atoms


def rotation_matrix_from_vectors(vec1, vec2):
    """Find the rotation matrix that aligns vec1 to vec2
    :param vec1: A 3d "source" vector
    :param vec2: A 3d "destination" vector
    :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """
    a, b = (
        (vec1 / np.linalg.norm(vec1)).reshape(3),
        (vec2 / np.linalg.norm(vec2)).reshape(3),
    )
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    if all(a == b):
        rotation_matrix = np.eye(3)
    else:
        kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s**2))
    return rotation_matrix


def random_three_vector():
    """
    Generates a random 3D unit vector.

    Returns:
    tuple: A tuple containing the x, y, and z components of the random 3D unit vector.
    """
    phi = np.random.uniform(0, np.pi * 2)
    costheta = np.random.uniform(-1, 1)

    theta = np.arccos(costheta)
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)
    return (x, y, z)


def random_rotate(atoms):
    """
    Applies a random rotation to the given atoms object.

    Parameters:
    atoms (object): An object containing an array of atoms to be rotated.

    Returns:
    object: The atoms object after applying the random rotation.
    """
    axis = random_three_vector()
    angle = np.random.uniform(0, 2 * np.pi)
    atoms.rotate(angle, axis)
    return atoms


def get_backbone_bond_change(trj, bond_cutoff=1.6):
    """
    Calculates the change in backbone bonds between the first and last frames of a trajectory.

    Parameters:
    trj (list): A list of atomic structures representing the trajectory.
    bond_cutoff (float): The distance cutoff to consider a bond. Default is 1.6.

    Returns:
    numpy.ndarray: An array representing the changes in backbone bonds.
    """
    H_count = len(trj[0][[atom.index for atom in trj[0] if atom.symbol == "H"]])

    traj = [trj[0].copy(), trj[-1].copy()]
    a_dict = {}

    for i in [0, -1]:
        traj[i] = traj[i][
            [atom.index for atom in traj[i] if atom.symbol in ["C", "O", "H"]]
        ]
        traj[i] = sort_atoms(traj[i], tags=traj[i].get_atomic_numbers())
        a_dict[i] = (traj[i].get_all_distances() < bond_cutoff) * 1
        # a_dict[i] = a_dict[i][H_count:,H_count:] # ignore "H-H" bonds
        a_dict[i][:H_count, :H_count] = 0  # ignore "H-H" bonds

    # print(a_dict.keys())
    a0 = a_dict[0]
    a1 = a_dict[-1]

    a = a_dict[-1] - a_dict[0]

    a = a[np.triu_indices(len(a))]

    # return a1,a0
    return sum(abs(a))


def read_relax_traj(file, pop_site_info=True):
    """
    Processes a trajectory file to extract and update information about the final atomic structure.

    Parameters:
    file (str): The path to the trajectory file.
    pop_site_info (bool): Whether to remove site-specific information from the final structure's info. Default is True.

    Returns:
    object: The final atomic structure with updated information.
    """
    traj = read(file, index=":")
    if len(traj) == 0:
        atoms = traj
    else:
        atoms = traj[-1]
        atoms.info["bond_change"] = get_backbone_bond_change(traj)
        atoms.info["snap_pos_compare"] = snap_pos_compare(traj[0], traj[-1])

    atoms.info["backbone_formula"] = atoms[
        [atom.index for atom in atoms if atom.symbol in ["C", "O"]]
    ].get_chemical_formula()
    # atoms.info['H_count'] = atoms[[atom.index for atom in atoms if atom.symbol in ['H']]].get_chemical_formula()
    if "adsorbate_formula_count" not in atoms.info.keys():
        atoms.info["adsorbate_formula_count"] = atoms[
            [atom.index for atom in atoms if atom.symbol in ["C", "H", "O"]]
        ].symbols.formula.count()

    adsorbate_formula_count = atoms.info["adsorbate_formula_count"]
    for k, v in adsorbate_formula_count.items():
        v = int(v)
    atoms.info.update(adsorbate_formula_count)

    if pop_site_info:
        pop_keys = [k for k in atoms.info.keys() if "ads_" in k] + [
            "adsorbate_formula_count"
        ]
        for k in pop_keys:
            atoms.info.pop(k)

    return atoms


def read_relax_dir(files):
    """
    Reads and processes a list of trajectory files, extracting information and creating a DataFrame.

    Parameters:
    files (list): A list of file paths to the trajectory files.

    Returns:
    tuple: A tuple containing:
        - rdf (pd.DataFrame): A DataFrame with information extracted from the trajectory files.
        - relaxed_traj (list): A list of atomic structures from the trajectory files.
    """
    files.sort()
    relaxed_traj = []
    rdf = []

    for i, file in enumerate(files):
        atoms = read_relax_traj(file)
        relaxed_traj.append(atoms)
        _a = atoms.copy()
        info = _a.info.pop("adsorbate_info")
        info.update(_a.info)
        info["traj_index"] = i
        #        print(info)
        rdf.append(info)

    rdf = pd.DataFrame(rdf)
    rdf = rdf.fillna(0.0)
    return rdf, relaxed_traj


def compute_energy(df, ref_dict, parent_en_dict):
    """
    Computes the energy of atomic structures based on reference energies and parent energy.

    Parameters:
    df (pd.DataFrame): A DataFrame containing information about the atomic structures.
    ref_dict (dict): A dictionary containing reference energies for different elements.
    parent_en (float): The energy of the parent structure.

    Returns:
    pd.DataFrame: The updated DataFrame with computed energy values.
    """
    for symbol in ["C", "H", "O"]:
        df[f"{symbol}_en"] = [ref_dict[symbol] for i in df.index.values]
    df["parent_en"] = [parent_en_dict[pid] for pid in df.pid.values]

    df["energy"] = df["mlff_energy"] - (
        df["parent_en"]
        + df["C"] * df["C_en"]
        + df["O"] * df["O_en"]
        + df["H"] * df["H_en"]
    )

    return df


def snap_pos_compare(atoms1, atoms2, sort=True, return_float=None):
    """
    Compares the positions of atoms in two atomic structures to determine their similarity.

    Parameters:
    atoms1 (object): The first atomic structure to compare.
    atoms2 (object): The second atomic structure to compare.
    sort (bool): Whether to sort the atoms by their positions before comparison. Default is True.
    return_float (float or None): The value to return if the chemical formulas of the structures do not match. Default is None.

    Returns:
    float: The sum of positional differences between the two structures. If the chemical formulas do not match, returns return_float.
    """
    if atoms1.get_chemical_formula() != atoms2.get_chemical_formula():
        return return_float

    d = 0

    for symbol in atoms1.symbols.formula.count().keys():
        a1 = atoms1[[atom.index for atom in atoms1 if atom.symbol == symbol]]
        a2 = atoms2[[atom.index for atom in atoms2 if atom.symbol == symbol]]
        if len(a1) != len(a2):
            return return_float

        trj = [a1, a2]

        for i, _ in enumerate(trj):
            for j in [0, 1, 2]:
                trj[i] = sort_atoms(trj[i], tags=trj[i].positions[:, j])

        d += _compare_pos(trj[0].positions, trj[1].positions)

    return d


def _compare_pos(pos1, pos2):
    """
    Computes the sum of Euclidean distances between corresponding positions in two arrays.

    Parameters:
    pos1 (numpy.ndarray): A 2D array of positions (shape: [n_atoms, 3]).
    pos2 (numpy.ndarray): A 2D array of positions (shape: [n_atoms, 3]).

    Returns:
    float: The sum of Euclidean distances between corresponding positions in pos1 and pos2.
    """
    return sum(np.linalg.norm(pos1 - pos2, axis=1))


def slice_traj_by_formula(traj):
    """
    Slices a trajectory into sub-trajectories based on unique chemical formulas.

    Parameters:
    traj (list): A list of atomic structures representing the trajectory.

    Returns:
    list: A list of lists, where each sublist contains atomic structures with the same chemical formula.
    """
    unique_formulas = list(set([atoms.get_chemical_formula() for atoms in traj]))

    return_list = []

    for formula in unique_formulas:
        return_list.append(
            [atoms for atoms in traj if atoms.get_chemical_formula() == formula]
        )

    return return_list


def get_drop_snapped(check_traj, d_cut, verbose=False):
    """
    Groups atomic structures in a trajectory based on positional similarity and drops duplicates.

    Parameters:
    check_traj (list): A list of atomic structures to be checked for positional similarity.
    d_cut (float): The distance cutoff to consider two structures as similar.
    verbose (bool): Whether to print detailed information during execution. Default is False.

    Returns:
    list: A list of atomic structures with duplicates removed based on positional similarity.
    """
    snapped_traj = []
    grouping = [0 for atoms in check_traj]

    n = 1
    while (np.array(grouping) == 0).any():
        unset_inds = [g for g, group in enumerate(grouping) if group == 0]
        ref_atoms = check_traj[unset_inds[0]]
        snapped_traj.append(ref_atoms)

        for i in unset_inds:
            d = snap_pos_compare(ref_atoms, check_traj[i], return_float=100)
            if d < d_cut:
                grouping[i] = n
        n += 1
        if verbose:
            print(f"sorted {len(check_traj) - len(unset_inds)} / {len(check_traj)}")

    return snapped_traj


def count_C_next_to_O(atoms):
    """
    Counts the number of carbon (C) atoms that are within a certain distance of an oxygen (O) atom.

    Parameters:
    atoms (object): An object containing an array of atoms.

    Returns:
    int: The number of carbon atoms within 1.6 Å of an oxygen atom. Returns 0 if no oxygen atom is present.
    """
    if "O" not in atoms.symbols:
        return 0

    O_index = [atom.index for atom in atoms if atom.symbol == "O"][0]

    xs = 0
    for atom in atoms:
        dist = atoms.get_distance(O_index, atom.index)
        if atom.symbol == "C" and dist < 1.6:
            xs += 1
    return xs


def polar2cart(theta, phi, r=1):
    """
    Converts polar coordinates to Cartesian coordinates.

    Parameters:
    theta (float): The polar angle in radians.
    phi (float): The azimuthal angle in radians.
    r (float): The radius. Default is 1.

    Returns:
    list: A list containing the x, y, and z Cartesian coordinates.
    """
    return [
        r * np.sin(theta) * np.cos(phi),
        r * np.sin(theta) * np.sin(phi),
        r * np.cos(theta),
    ]


def get_sorted_by_snap_dist(traj):
    """
    Sorts a list of atomic structures based on their snap position comparison distance.

    Parameters:
    traj (list): A list of atomic structures to be sorted.

    Returns:
    list: A list of atomic structures sorted by their snap position comparison distance.
    """

    ref_atoms = traj[0].copy()
    lst = [(atoms, snap_pos_compare(ref_atoms, atoms)) for atoms in traj]
    lst = sorted(lst, key=lambda tup: tup[1])

    slice_index = int(np.floor(len(lst) * 0.5))

    if len(traj) % 2 != 0:
        lst = lst[1:]

    a = lst[:slice_index]
    b = list(reversed(lst[slice_index:]))
    lst = list(itertools.chain.from_iterable(zip(a, b)))

    if len(traj) % 2 != 0:
        lst = [(ref_atoms, 0.0)] + lst

    return list(list(zip(*lst))[0])


def make_site_info_writable(site):
    """
    Converts numpy arrays in a site dictionary to lists to make the dictionary JSON serializable.

    Parameters:
    site (dict): A dictionary containing site information.

    Returns:
    dict: The updated dictionary with numpy arrays converted to lists.
    """
    for k, v in site.items():
        if type(v) == np.ndarray:
            site[k] = list(v)
    return site


def docs_plot_conformers(conformer_trajectory, rotation="-90x,0y,0z"):
    """
    Helper function to plot a series of conformers.

    Parameters:
    conformer_trajectory (list): A list of atomic structures representing the conformers.
    rotation (str): The rotation to apply to the plot. Default is '-90x,0y,0z'.

    Returns:
    matplotlib.figure.Figure: The figure object containing the plots.
    """
    fig, ax = plt.subplots(1, 5, figsize=(10, 2), sharex=True, sharey=True)

    for i, atoms in enumerate(ax):
        plot_atoms(conformer_trajectory[i], ax=ax[i], rotation=rotation)
        ax[i].set_xlim(0, 7), ax[i].set_ylim(0, 7)
        ax[i].set_axis_off()
    fig.suptitle("Generated structures viewed from +X axis", fontsize=12)
    fig.tight_layout()
    return fig


def docs_plot_sites(surface_object, rotation="-45x,0y,0z"):
    """
    Helper function to plot a series of sites from a surface object.

    Parameters:
    surface_object (object): An object containing surface site information.
    rotation (str): The rotation to apply to the plot. Default is '-45x,0y,0z'.

    Returns:
    matplotlib.figure.Figure: The figure object containing the plots.
    """
    import random

    inds = random.sample(list(surface_object.site_df.index.values), 6)

    fig, ax = plt.subplots(1, 6, figsize=(10, 4), sharex=True, sharey=True)

    for ax_i, i in enumerate(inds):
        view_atoms = surface_object.view_site(i, return_atoms=True)
        parent_atoms = surface_object.atoms.copy()

        for atom in parent_atoms:
            if atom.index not in view_atoms.info["topology"]:
                parent_atoms[atom.index].symbol = "Zn"

        view_atoms = parent_atoms + view_atoms
        plot_atoms(view_atoms, ax=ax[ax_i], rotation=rotation)
        ax[ax_i].set_axis_off()
    fig.suptitle(
        "Generated structures viewed from +X+Z axis, site atoms shown as Cu, other as Zn",
        fontsize=12,
    )
    fig.tight_layout()
    return fig


def freeze_atoms(atoms, from_top=3):
    """
    Freezes the bottom layers of atoms in a structure to prevent them from moving during relaxation.

    Parameters:
    atoms (object): An ASE atoms object representing the atomic structure.
    from_top (int): The number of layers from the top to exclude from freezing. Default is 3.

    Returns:
    object: The ASE atoms object with the bottom layers frozen.
    """

    c = FixAtoms(
        indices=[
            atom.index
            for atom in atoms
            if atom.position[2]
            < max(
                [
                    a.position[2]
                    for a in atoms
                    if a.symbol not in ["C", "O", "S", "N", "H"]
                ]
            )
            - from_top
        ]
    )
    atoms.set_constraint(c)
    return atoms


def save_relax_config(dyn, fname):
    """
    Saves the current configuration of the atoms during a relaxation process.

    Parameters:
    dyn (object): The ASE dynamics object containing the atoms being relaxed.
    fname (str): The filename where the configuration will be saved.

    Returns:
    None
    """
    atomsi = dyn.atoms
    en = atomsi.get_potential_energy()
    frcs = atomsi.get_forces()

    atomsi.info.update(
        {
            "mlff_energy": en,
        }
    )
    atomsi.arrays.update({"mlff_forces": frcs})
    ase.io.write(fname, atomsi, append=True)


def is_clear_to_start(fname):
    """
    Checks if a file exists and creates it if it does not.

    Parameters:
    fname (str): The path to the file.

    Returns:
    bool: False if the file exists, True if the file did not exist and was created.
    """
    import os.path
    from pathlib import Path

    if os.path.isfile(fname):
        return False
    else:
        Path(fname).touch()
        return True


def relaxatoms(atoms, calc, prefix, steps=300, freeze_bottom=False, fmax=0.05):
    """
    Relaxes the atomic positions of a given structure using a specified calculator and saves the relaxation trajectory.

    Parameters:
    atoms (object): An ASE atoms object representing the atomic structure to be relaxed.
    calc (object): An ASE calculator object to be used for the relaxation.
    prefix (str): A prefix for the filename where the relaxation trajectory will be saved.
    steps (int): The maximum number of optimization steps. Default is 300.
    freeze_bottom (bool): Whether to freeze the bottom layers of atoms during relaxation. Default is False.

    Returns:
    float: The potential energy of the relaxed atomic structure.
    """
    fname = f"{prefix}_{atoms.info['uid']}.xyz"
    atoms.calc = calc

    if freeze_bottom:
        atoms = freeze_atoms(atoms, from_top=3)

    dyn = BFGS(atoms)
    dyn.attach(save_relax_config, interval=1, dyn=dyn, fname=fname)
    dyn.run(fmax=fmax, steps=steps)

    return atoms.get_potential_energy()


def minhopatoms(atoms, calc, prefix, steps=10, fmax=0.05, temperature=1000):
    fname = f"{prefix}_{atoms.info['uid']}.xyz"
    atoms.calc = calc

    dyn = BFGS(atoms)

    opt = MinimaHopping(
        atoms=atoms,
        T0=1000.0,  # K, initial MD ‘temperature’
        Ediff0=0.5,  # eV, initial energy acceptance threshold
        logfile=f"hop_{fname}.log",  # text log
        minima_threshold=0.5,  # A, threshold for identical configs
        timestep=1.0,
        minima_traj=f"{fname}.traj",
        fmax=fmax,
    )
    #    opt.attach(save_relax_config, interval=1, dyn=dyn, fname=fname)
    opt(totalsteps=steps)

    # dyn.run(fmax=fmax, steps=steps)

    return atoms.get_potential_energy()


def get_blenderized(traj, scale=[1, 1, 1], hide_spot="default"):
    """
    Prepares a trajectory for visualization in Blender by scaling and adding hidden atoms to ensure consistent atom counts.

    Parameters:
    traj (list): A list of atomic structures representing the trajectory.
    scale (list): A list of three integers representing the scaling factors for the x, y, and z dimensions. Default is [3, 3, 1].
    hide_spot (str or list): The position to place hidden atoms. If 'default', the position is set to the center of the cell. Default is 'default'.

    Returns:
    list: A list of atomic structures prepared for visualization in Blender.
    """
    blenderized_trj = []
    max_atoms_dict = get_max_atoms_dict(traj)

    for t in traj:
        a = t.copy()
        a.positions += [0, 0.4, 0]
        a.wrap()

        if hide_spot == "default":
            hide_spot = a.cell[0] * 0.5 + a.cell[1] * 0.5

        for s, no in max_atoms_dict.items():
            n_atom_to_add = no - len([atom for atom in a if atom.symbol == s])

            for n in range(n_atom_to_add):
                a.append(Atom(s, hide_spot))

        a = a * scale
        blenderized_trj.append(a)

    return blenderized_trj


def get_max_atoms_dict(traj):
    """
    Computes the maximum number of each type of atom present in any structure within a trajectory.

    Parameters:
    traj (list): A list of atomic structures representing the trajectory.

    Returns:
    dict: A dictionary where the keys are atomic symbols and the values are the maximum number of atoms of that type found in any structure within the trajectory.
    """
    max_atoms_dict = {}
    for t in traj:
        for s in t.symbols.species():
            if s not in max_atoms_dict.keys():
                max_atoms_dict[s] = []

            max_atoms_dict[s].append(
                len(t[[atom.index for atom in t if atom.symbol == s]])
            )
    for k, v in max_atoms_dict.items():
        max_atoms_dict[k] = max(v)
    return max_atoms_dict


# def conformer_to_site(atoms, site, conformer, mode='optimize', overlap_thr = 0):
#
#    atoms = atoms.copy()
#    conformer = conformer.copy()
#
#    if 'fragments' not in atoms.arrays.keys():
#        atoms.arrays['fragments'] = np.array([0 for a in atoms])
#
#    n_f = int(np.max(atoms.arrays['fragments'])) + 1
#    conformer.arrays['fragments'] = np.array([n_f for a in conformer])
#
#    if conformer[0].symbol == 'S' and conformer[1].symbol == 'S':
#        conformer.positions -= conformer[:2].get_center_of_mass()
#
#    conformer=align_to_vector(conformer, site['n_vector'])
#
#    h_rot_angle = np.sign(site['h_vector'][1]) * np.arccos(site['h_vector'][0]) / np.pi * 180
#    conformer.rotate(h_rot_angle, site['n_vector'])
#    conformer.positions += site['coordinates']
#
#    out_atoms = atoms+conformer
#    out_atoms.info['site_info'] = site
#
#    if conformer[0].symbol == 'S' and conformer[1].symbol == 'S':
#        out_atoms = swing_fragment(
#            atoms = out_atoms,
#            fragment_index=n_f,
#            site = site,
#            resolution = 10,
#            mode = mode,
#            span_angle=50,
#            overlap_thr=0
#        )
#
#    if conformer[0].symbol == 'Cl':
#        out_atoms = swirl_fragment(
#            atoms = out_atoms,
#            fragment_index = n_f,
#            site = site,
#            resolution = 10,
#            mode = 'optimize',
#            overlap_thr=0.
#        )
#
#    return out_atoms
#
# def conformer_to_site(atoms, site, conformer, mode='optimize', overlap_thr = 0):
#
#    atoms = atoms.copy()
#    conformer = conformer.copy()
#
#    if 'fragments' not in atoms.arrays.keys():
#        atoms.arrays['fragments'] = np.array([0 for a in atoms])
#
#    n_f = int(np.max(atoms.arrays['fragments'])) + 1
#    conformer.arrays['fragments'] = np.array([n_f for a in conformer])
#
#    if conformer[0].symbol == 'S' and conformer[1].symbol == 'S':
#        conformer.positions -= conformer[:2].get_center_of_mass()
#
#    conformer=align_to_vector(conformer, site['n_vector'])
#
#    h_rot_angle = np.sign(site['h_vector'][1]) * np.arccos(site['h_vector'][0]) / np.pi * 180
#    conformer.rotate(h_rot_angle, site['n_vector'])
#    conformer.positions += site['coordinates']
#
#    out_atoms = atoms+conformer
#    out_atoms.info['site_info'] = site
#
#    if conformer[0].symbol == 'S' and conformer[1].symbol == 'S':
#        out_atoms = swing_fragment(
#            atoms = out_atoms,
#            fragment_index=n_f,
#            site = site,
#            resolution = 10,
#            mode = mode,
#            span_angle=50,
#            overlap_thr=0
#        )
#
#    if conformer[0].symbol == 'Cl':
#        out_atoms = swirl_fragment(
#            atoms = out_atoms,
#            fragment_index = n_f,
#            site = site,
#            resolution = 10,
#            mode = 'optimize',
#            overlap_thr=0.
#        )
#
#    return out_atoms

# def write_chemiscope(trajectory, write_path='./', filename='chemiscope.json.gz'):
#     """Get a chemiscope with all the atoms.info of a trajectory.

#     Args:
#         trajectory (list of Atoms): _description_
#         filename (str): Defaults to 'chemiscope.json.gz'
#         write_path (str, optional): Where to write chemiscope.json.gz. Defaults to './'.
#     """
#     from chemiscope import write_input

#     properties={}
#     info_keys = set(trajectory[0].info.keys())


#     for prop in info_keys:

#         properties[prop]={"target":"structure",
#                         "values":[atoms.info[prop] for atoms in trajectory]
#                         }

#     print(f'Writing {os.path.join(write_path,filename)}')
#     write_input(os.path.join(write_path,filename),
#                 frames=trajectory,
#                 properties=properties,
#                 # environments=[(i, [atom.index for atom in a if atom.symbol == 'He'][0], 4.5) for i, a in enumerate(check_trj)]
#             )
#     return

# def do_TSNE(
#         soapdesc,
#         n_components=2,
#         perplexity='default',
#         metric="euclidean",
#         n_jobs=25,
#         random_state=42,
#         verbose=False,
#     ):

#     if perplexity == 'default':  perplexity = int(len(soapdesc)*0.05)

#     from openTSNE import TSNE

#     tsne = TSNE(
#         n_components=n_components,
#         perplexity=perplexity,
#         metric=metric,
#         n_jobs=n_jobs,
#         random_state=random_state,
#         verbose=verbose,
#     )
#     env_embedding=tsne.fit(soapdesc)

#     return env_embedding


# def get_soapdesc(trajectory:list,
#         element:str='all',
#         soap_cutoff:float=4.5
#     ):
#     """_summary_

#     Args:
#         trajectory (list): list of ase.Atoms structures for a subset to be selected.
#         element (str, optional): Element eg 'H', 'Zr', soap descriptor will be constructed for all sites of this element. Defaults to 'all'.
#         soap_cutoff (float, optional): Defaults to 4.5.

#     Returns: normalized soapdesc.
#     """
#     from soapml.descriptors.soap import SOAP
#     from sklearn.preprocessing import normalize
#     import numpy as np
#     from soapml.selection.BasicSelection import AtomSelector
#     from pymatgen.core import Element

#     traj = trajectory.copy()

#     if element == 'all':
#         print(f'    Generating SOAP with {soap_cutoff} cuttoff for all atoms in trajectory...')
#         soap=SOAP(rcut=4.5,nmax=8,lmax=4,sigma=0.5,periodic=True,rbf='gto',crossover=True)
#         soap.fit(traj,site_to_structure_method='inner')
#         soapdesc=soap.featurize_many(traj,n_jobs=20)

#     else:
#         print(f'    Selecting {element} atoms in the trajectory...')
#         selector=AtomSelector()
#         selector.fit({'numbers':[Element(element).number]}) #atomic number of He that is the marker
#         selected_atoms=selector.transform(traj)

#         print(f'    Generating SOAP with {soap_cutoff} cuttoff...')
#         soap=SOAP(rcut=soap_cutoff,nmax=8,lmax=4,sigma=0.5,periodic=True,rbf='gto',crossover=True)
#         soap.fit(traj,site_to_structure_method='off')
#         soapdesc=soap.featurize_many(traj,idx_list=selected_atoms, n_jobs=20)

#     soapdesc = normalize(soapdesc, axis=1)
#     return soapdesc

# def plot_tsne_selection(selected_filter:list, env_embedding:list, show_plot = True):
#     """_summary_

#     Args:
#         selected_filter (list): list of ints eg [1,0,0,1,0], 1 mean selected, 0 means ot selected.
#         env_embedding (list): eg tsne coordinates.
#     """
#     import plotly.express as px
#     import pandas as pd
#     import seaborn as sns

#     if len(selected_filter) != len(env_embedding):
#         print('Lists must be same lenght!')

#     else:
#         df = []
#         for i, _ in enumerate(selected_filter):
#             data = {'selected': selected_filter[i]}

#             for j, val in enumerate(env_embedding[i]):
#                 data[f'tsne{j}'] = val
#             df.append(data)
#         df = pd.DataFrame(df)
#         df = df.sort_values(by='selected', ascending=True)

#         if show_plot:
#             # px.defaults.width=800
#             # px.defaults.height=600
#             # print('Plotting first two tsne components.')
#             # return px.scatter(df, x='tsne0', y='tsne1', color='selected')
#             return sns.scatterplot(data=df, x='tsne0', y='tsne1', hue='selected')

#         return(df)


# def quick_soap_fps(
#         trajectory:list,
#         n_to_select:int,
#         element:str='all',
#         soap_cutoff:float=4.5,
#         plot_tsne=False,
#         get_chemiscope=False,
#         soapdesc='generate'
#     ):
#     """Selects subset of desired lenght based on most different soap.

#     Args:
#         trajectory (list): list of ase.Atoms structures for a subset to be selected.
#         n_to_select (int): number of structures to be selected from trajectory.
#         element (str, optional): Element eg 'H', 'Zr', soap descriptor will be constructed for all sites of this element. Defaults to 'all'.
#         soap_cutoff (float, optional): Defaults to 4.5.

#     Returns:
#         selected_trajectory (list): list of most unique ase.Atoms structures.
#     """
#     from soapml.selection.SampleSelection import FPS

#     if soapdesc == 'generate':
#         soapdesc = get_soapdesc(trajectory, element=element)
#         print(len(soapdesc))

#     print(f'    Starting FPS, selecting {n_to_select} configurations...')
#     small_fps = FPS(initialize=0, n_to_select=n_to_select)
#     small_fps.fit(soapdesc)

#     print('    Done!')
#     selected_traj = []
#     selected_filter = []

#     for i, a in enumerate(trajectory):
#         sf = 0
#         if i in small_fps.selected_idx_:
#             selected_traj.append(a.copy())
#             sf=1
#         selected_filter.append(sf)

#     if plot_tsne:
#         env_embedding = do_TSNE(soapdesc)
#         return plot_tsne_selection(selected_filter, env_embedding, show_plot = True)

#     return selected_traj


# def write_packmol_input(slab, tolerance = 2.0, output = 'out.xyz', structure = './input.xyz', number = 1):

#     xmin = 0 ; xmax = slab.cell[0][0] - 1
#     ymin = 0 ; ymax = slab.cell[1][1] - 1
#     zmin = max([a.position[2] for a in slab]) + 2 ; zmax = slab.cell[2][2] - 2

#     ############################
#     packmol_template = """
#     tolerance {}

#     filetype xyz

#     output {}

#     structure {}
#       number {}
#       inside box {} {} {} {} {} {}
#     end structure

#     """.format(
#         tolerance,
#         output,
#         structure,
#         number,
#         xmin, ymin, zmin,
#         xmax, ymax, zmax
#     )

#     path = os.path.join('./run_packmol.inp')
#     f = open(path, 'w')
#     f.write(packmol_template)
#     f.close()

#     # cmd = [
#     #     '/gpfs/nobackup/projects/qm_inorganics/notebook_edvin/bin/packmol-20.3.5/packmol ',
#     #     path
#     # ]
#     # cmd = " ".join(cmd)
#     # run_subprocess(cmd)


# def write_packmol_input_mix(slab, tolerance = 2.0, output = 'out.xyz', structure1 = './input1.xyz',
#                             structure2 = './input2.xyz', number = 1, fraction1=0.5):

#     xmin = 0 ; xmax = slab.cell[0][0] - 1
#     ymin = 0 ; ymax = slab.cell[1][1] - 1
#     zmin = max([a.position[2] for a in slab]) + 2 ; zmax = slab.cell[2][2] - 2

#     number1 = int(number*fraction1)
#     number2 = int(number*(1 - fraction1))

#     ############################
#     packmol_template = """
#     tolerance {}

#     filetype xyz

#     output {}

#     structure {}
#       number {}
#       inside box {} {} {} {} {} {}
#     end structure

#     structure {}
#       number {}
#       inside box {} {} {} {} {} {}
#     end structure

#     """.format(
#         tolerance,
#         output,
#         structure1,
#         number1,
#         xmin, ymin, zmin,
#         xmax, ymax, zmax,
#         structure2,
#         number2,
#         xmin, ymin, zmin,
#         xmax, ymax, zmax
#     )

#     path = os.path.join('./run_packmol.inp')
#     f = open(path, 'w')
#     f.write(packmol_template)
#     f.close()

#     # cmd = [
#     #     '/gpfs/nobackup/projects/qm_inorganics/notebook_edvin/bin/packmol-20.3.5/packmol ',
#     #     path
#     # ]
#     # cmd = " ".join(cmd)
#     # run_subprocess(cmd)

# def detect_unique_structures(structures, prec=.001, rcut=5.5, element='all'):
#     """
#     Author: Sandip De, Edvin Fako
#     Given a list of pymatgen structure object this routine identify the unique structure index using SOAP descriptors.

#     """

#     from pymatgen.io.ase import AseAtomsAdaptor
#     from soapml.descriptors.soap import SOAP
#     from sklearn.preprocessing import normalize
#     import numpy as np
#     from openTSNE import TSNE
#     from soapml.selection.SampleSelection import FPS
#     from soapml.selection.BasicSelection import AtomSelector
#     from pymatgen.core import Element
#     from scipy.spatial.distance import pdist, squareform

#     species = list(
#         structures[0].composition.element_composition.as_dict().keys())

#     print("Using species:{}. Selecting: '{}'".format(species, element))
#     rcut = rcut
#     nmax = 8
#     lmax = 4
#     sigma = 0.5

#     atoms = [AseAtomsAdaptor.get_atoms(s) for s in structures]


#     if element == 'all':
#         soap=SOAP(rcut=rcut,nmax=nmax,lmax=lmax,sigma=sigma,periodic=True,rbf='gto',crossover=True)
#         soap.fit(atoms,site_to_structure_method='inner')
#         soapdesc=soap.featurize_many(atoms,n_jobs=20)

#     else:
#         selector=AtomSelector()
#         selector.fit({'numbers':[Element(element).number]}) #atomic number of He that is the marker
#         selected_atoms=selector.transform(atoms)

#         soap=SOAP(rcut=rcut,nmax=nmax,lmax=lmax,sigma=sigma,periodic=True,rbf='gto',crossover=True)
#         soap.fit(atoms,site_to_structure_method='off')
#         soapdesc=soap.featurize_many(atoms,idx_list=selected_atoms, n_jobs=20)

#     soapdesc = normalize(soapdesc, axis=1)

#     dist = squareform(pdist(soapdesc))
#     skip = np.array([False for i in range(len(dist))])
#     ilist, jlist = np.where(dist < prec)

#     for i, j in zip(ilist, jlist):
#         if not skip[i] and j > i: skip[j] = True
#     (uid, ) = np.where(skip == False)

#     return (uid)
