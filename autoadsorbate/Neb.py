import itertools
from typing import Union

import numpy as np
import pandas as pd
from ase import Atoms
from ase.neb import NEB


def permute_image(atoms, fix_species=None):
    """
    Permutes the positions of atoms in the given atomic structure, optionally fixing the positions of specified species.

    Parameters:
    atoms (object): An ASE atoms object representing the atomic structure.
    fix_species (list or None): A list of atomic numbers representing species whose positions should not be permuted. Default is None.

    Returns:
    list: A list of ASE atoms objects, each representing a unique permutation of the atomic structure.
    """
    ind_dict = {}

    atoms_numbers = list(set(atoms.arrays["numbers"]))
    atoms_numbers.sort()

    perturb_numbers = atoms_numbers.copy()
    if fix_species != None:
        for n in fix_species:
            perturb_numbers.remove(n)
            ind_dict[n] = [[atom.index for atom in atoms if atom.number == n]]

    for n in perturb_numbers:
        ind_dict[n] = []

        inds = [atom.index for atom in atoms if atom.number == n]

        for perm in [p for p in itertools.permutations(inds)]:
            ind_dict[n].append(list(perm))

    arrays = [v for _, v in ind_dict.items()]

    all_perm_inds = []

    for i in itertools.product(*arrays):
        l = [item for sublist in i for item in sublist]
        all_perm_inds.append(l)

    perm_ini_traj = []
    for perm in all_perm_inds:
        a = atoms.copy()
        a = a[perm]
        perm_ini_traj.append(a)

    return perm_ini_traj


def get_connectivity(atoms, bond_range=[0.0, 2]):
    """
    Computes the connectivity matrix for the given atomic structure based on the specified bond range.

    Parameters:
    atoms (object): An ASE atoms object representing the atomic structure.
    bond_range (list): A list containing the minimum and maximum bond distances to consider. Default is [0., 2.].

    Returns:
    numpy.ndarray: A 2D array representing the connectivity matrix, where each element indicates the presence of a bond.
    """
    distm = atoms.get_all_distances()
    distm = np.triu(distm, k=0)

    bondm = np.logical_and(distm > bond_range[0], distm < bond_range[1]) * 1
    distm = distm * bondm

    return distm


def get_neb_images(ini, fin, images_no=10, method="linear"):
    """
    Generates a series of images for the Nudged Elastic Band (NEB) method by interpolating between initial and final structures.

    Parameters:
    ini (object): An ASE atoms object representing the initial structure.
    fin (object): An ASE atoms object representing the final structure.
    images_no (int): The number of intermediate images to generate. Default is 10.
    method (str): The interpolation method to use. Default is 'linear'.

    Returns:
    list: A list of ASE atoms objects representing the interpolated images.
    """
    traj = [ini.copy() for _ in range(images_no)] + [fin.copy()]
    neb = NEB(traj)
    neb.interpolate(method=method)

    return neb.images


def get_neb_norm(neb_images):
    """
    Computes the NEB (Nudged Elastic Band) norm, which is the minimum distance between atoms across all NEB images.

    Parameters:
    neb_images (list): A list of ASE atoms objects representing the NEB images.

    Returns:
    float: The minimum distance between atoms across all NEB images.
    """
    im_dist = []
    for i, im in enumerate(neb_images):
        d = im.get_all_distances()
        d = d[d > 0]
        im_dist.append(np.min(d))

    neb_norm = np.min(im_dist)
    return neb_norm


def get_neb_dists(neb_images):
    """
    Computes the minimum distances between atoms for each NEB (Nudged Elastic Band) image.

    Parameters:
    neb_images (list): A list of ASE atoms objects representing the NEB images.

    Returns:
    tuple: A tuple containing two lists:
        - inds (list): The indices of the NEB images.
        - im_dist (list): The minimum distances between atoms for each NEB image.
    """
    im_dist = []
    inds = []
    for i, im in enumerate(neb_images):
        d = im.get_all_distances()
        d = d[d > 0]
        im_dist.append(np.min(d))
        inds.append(i)

    return inds, im_dist


def get_distm(ini, fin):
    """
    Computes the absolute difference in connectivity matrices between the initial and final atomic structures.

    Parameters:
    ini (object): An ASE atoms object representing the initial structure.
    fin (object): An ASE atoms object representing the final structure.

    Returns:
    numpy.ndarray: A 2D array representing the absolute difference in connectivity matrices between the initial and final structures.
    """
    distm_fin = get_connectivity(fin)
    distm_ini = get_connectivity(ini)
    distm = distm_fin - distm_ini

    distm = np.abs(distm)

    return distm


def get_distm_sum(ini, fin, unit="A", tolerance=0.5):
    """
    Computes the sum of the absolute differences in connectivity matrices between the initial and final atomic structures.

    Parameters:
    ini (object): An ASE atoms object representing the initial structure.
    fin (object): An ASE atoms object representing the final structure.
    unit (str): The unit of measurement for the differences. Can be 'A' for angstroms or 'bonds' for covalent bonds. Default is 'A'.
    tolerance (float): The tolerance value for considering a difference as significant. Default is 0.5.

    Returns:
    float: The sum of the absolute differences in connectivity matrices, considering the specified unit and tolerance.

    Raises:
    ValueError: If the unit is not 'A' or 'bonds'.
    """
    distm = get_distm(ini, fin)
    distm = np.abs(distm)

    if unit == "bonds":
        distm = distm > tolerance
        return np.sum(distm)

    elif unit == "A":
        distm = distm * (distm > tolerance)
        return np.sum(distm)

    else:
        raise ValueError(
            "Value of variable 'unit' can be 'A' - angstrom or 'bonds' - covalent bond"
        )


def plot_distm(ini, fin):
    """
    Plots the connectivity matrices and their differences for the initial and final atomic structures, along with the NEB distances.

    Parameters:
    ini (object): An ASE atoms object representing the initial structure.
    fin (object): An ASE atoms object representing the final structure.

    Returns:
    None
    """
    import seaborn as sns

    f, (ax1, ax2, ax3, ax4) = plt.subplots(
        1, 4, figsize=(12, 2.3)
    )  # , gridspec_kw = {'width_ratios': [1, 3]})

    distm_fin = Neb.get_connectivity(fin)
    distm_ini = Neb.get_connectivity(ini)
    distm = distm_fin - distm_ini

    neb_images = get_neb_images(ini, fin)
    x, y = get_neb_dists(neb_images)

    sns.heatmap(distm_ini, linewidth=0.5, ax=ax1)
    sns.heatmap(distm_fin, linewidth=0.5, ax=ax2)
    sns.heatmap(np.abs(distm), linewidth=0.5, ax=ax3)
    sns.scatterplot(x=x, y=y, ax=ax4)
    f.tight_layout()


def arrange_backbone(ini, fin, bond_len_tolerance=0.5):
    """
    Arranges the backbone of the initial and final atomic structures by permuting the initial structure to best match the final structure.

    Parameters:
    ini (object): An ASE atoms object representing the initial structure.
    fin (object): An ASE atoms object representing the final structure.
    bond_len_tolerance (float): The tolerance value for considering bond length differences as significant. Default is 0.5.

    Returns:
    tuple: A tuple containing two ASE atoms objects:
        - best_ini (object): The permuted initial structure that best matches the final structure.
        - best_fin (object): The final structure.
    """
    ini = ini[ini.numbers.argsort()]
    fin = fin[fin.numbers.argsort()]

    h = ini.copy()
    h = h[[atom.index for atom in h if atom.symbol == "H"]]

    h_fin = fin.copy()
    h_fin = h_fin[[atom.index for atom in h_fin if atom.symbol == "H"]]

    b_ini = ini.copy()
    b_ini = b_ini[[atom.index for atom in b_ini if atom.symbol != "H"]]

    b_fin = fin.copy()
    b_fin = b_fin[[atom.index for atom in b_fin if atom.symbol != "H"]]

    b_ini_perms = permute_image(b_ini)

    best_ini_perm = get_best_perm(
        b_ini_perms, b_fin, bond_len_tolerance=bond_len_tolerance
    )

    if fin[0].symbol == "H":
        best_fin = h_fin.copy()
        best_fin += b_fin

        best_ini = h.copy()
        best_ini += best_ini_perm
    else:
        best_fin = b_fin.copy()
        best_fin += h_fin

        best_ini = best_ini_perm.copy()
        best_ini += h

    # print('best_ini: ', best_ini)
    # print('best_fin: ', best_fin)

    return best_ini, best_fin


def get_sorted(atoms):
    """
    Sorts the atoms in the given atomic structure by their atomic numbers.

    Parameters:
    atoms (object): An ASE atoms object representing the atomic structure.

    Returns:
    object: An ASE atoms object with atoms sorted by their atomic numbers.
    """
    atoms_numbers = list(set(atoms.arrays["numbers"]))
    atoms_numbers.sort()

    inds = []
    for n in atoms_numbers:
        inds += [atom.index for atom in atoms if atom.number == n]

    b = atoms.copy()
    b = b[inds]
    return b


def get_best_perm(ini_perms, fin, bond_len_tolerance, plot=False):
    """
    Finds the best permutation of the initial atomic structure that best matches the final structure based on bond length tolerance.

    Parameters:
    ini_perms (list): A list of ASE atoms objects representing the permuted initial structures.
    fin (object): An ASE atoms object representing the final structure.
    bond_len_tolerance (float): The tolerance value for considering bond length differences as significant.
    plot (bool): Whether to plot the connectivity matrices and their differences. Default is False.

    Returns:
    object: The best permuted ASE atoms object that matches the final structure.
    """
    best_ds = 1000000000000000  # very big number
    fin = get_sorted(fin)

    df = []
    for j, a in enumerate(ini_perms):
        a = get_sorted(a)
        # print('ini: ', a)
        # print('fin: ', fin)

        if plot:
            plot_distm(a, fin)

        neb_images = get_neb_images(a, fin)

        ds = get_distm_sum(a, fin, unit="bonds", tolerance=bond_len_tolerance)
        if ds <= best_ds:
            df.append({"ds": ds, "touch": get_neb_norm(neb_images), "perm_index": j})
            best_ds = ds
    df = pd.DataFrame(df)

    min_ds = np.min(df.ds.values)

    df = df.sort_values(by="touch", ascending=False)

    df = df[df.ds == min_ds]
    print(df)
    best_index = df.perm_index.values[0]
    best_touch = df.touch.values[0]
    print("best index: ", best_index)
    best_perm = ini_perms[best_index].copy()
    best_perm = get_sorted(best_perm)

    info = {"ds": min_ds, "touch": best_touch}
    best_perm.info.update(info)

    return best_perm


def arrange_images(ini, fin, bond_len_tolerance=0.5, plot=False):
    """
    Arranges the initial and final atomic structures by permuting the initial structure to best match the final structure.

    Parameters:
    ini (object): An ASE atoms object representing the initial structure.
    fin (object): An ASE atoms object representing the final structure.
    bond_len_tolerance (float): The tolerance value for considering bond length differences as significant. Default is 0.5.
    plot (bool): Whether to plot the connectivity matrices and their differences. Default is False.

    Returns:
    tuple: A tuple containing two ASE atoms objects:
        - best_ini (object): The permuted initial structure that best matches the final structure.
        - fin (object): The final structure.
    """
    t = []

    ini, fin = arrange_backbone(ini, fin, bond_len_tolerance=bond_len_tolerance)
    fix_species = [atom.number for atom in ini if atom.number != 1]
    perms = permute_image(ini, fix_species=[6, 8])
    print(len(perms))
    best_ini = get_best_perm(perms, fin, bond_len_tolerance=0.1, plot=plot)
    return best_ini, fin


# def make_best_nebs(ini, fin, rnorm_cutoff=2, bond_range = [0., 2], neb_images=10, neb_dist=0.7, verbose=False):
#
#    ini = ini[ini.numbers.argsort()]
#    fin = fin[fin.numbers.argsort()]
#
#    permuted_ini = permute_image(ini)
#
#    best_nebs = []
#
#    best_ini = select_best_neb(permuted_ini, fin, rnorm_cutoff=rnorm_cutoff, bond_range = bond_range, neb_images = neb_images, neb_dist = neb_dist, verbose=False)
#
#    if best_ini:
#
#        neb_trj = [best_ini.copy() for _ in range(10)]
#        neb_trj.append(fin.copy())
#
#        neb = NEB(neb_trj)
#
#        neb.interpolate()
#
#        best_nebs.append(neb.images)
#
#    return best_nebs

# def check_neb(iatoms, fatoms, bonds=1, bond_range = [0., 2]):
#
#    reactm = check_image(iatoms) - check_image(fatoms)
#
#    flat_reactm = reactm.flatten()
#    flat_reactm = flat_reactm[flat_reactm>0]
#    flat_reactm.sort()
#
#    return np.linalg.norm(reactm), np.linalg.norm(flat_reactm[:-bonds])
#
# def select_best_ini(ini, fin, rnorm_cutoff=3, bond_range = [0., 2], neb_images=10, neb_dist=0.8, verbose=False):
#    df = []
#
#    perm_ini_traj = permute_atoms(ini)
#
#    for i, a in enumerate(perm_ini_traj):
#
#        INI = perm_ini_traj[i].copy()
#        FIN = fin.copy()
#
#        rnorm, _ = check_neb(INI, FIN, bonds=1, bond_range = [0., 2])
#
#        neb = NEB([INI.copy() for _ in range(neb_images)]+[FIN])
#        neb.interpolate()
#
#        im_dist = []
#        for im in neb.images:
#            d = im.get_all_distances()
#            d = d[d>0]
#            im_dist.append(np.min(d))
#        neb_norm = np.min(im_dist)
#
#        info = {'neb_index':i, 'rnorm': rnorm, 'neb_norm': neb_norm}
#        df.append(info)
#    df = pd.DataFrame(df)
#    df=df[df['rnorm'] < rnorm_cutoff]
#    df=df[df['neb_norm'] > neb_dist]
#    if verbose:
#        print(df)
#    df = df.sort_values(by='rnorm', ascending=True)[:1]
#    if len(df.index.values) > 0:
#        best_index = df.neb_index.values[0]
#        best_ini = perm_ini_traj[best_index].copy()
#        best_ini.info.update(info)
#        return best_ini
#    else:
#        if verbose:
#            print(f'select_best_neb returns empty list.')
#        return False


# Credits Lars Leon Schaaf

# def permute_like_species(sp_t, select):
#    """
#    Returns new index of all posible permutations between like elements of the species
#    total array (spt_t), but only for the select configurations
#
#    Args:
#        spt_t (list): list of species
#        select (list): if only a subset of the indicies need to be permuted
#
#    Returns:
#        list: index of new permutations
#    """
#    sp_t = np.array(sp_t)
#    select = np.array(select)
#
#    # Indexes -> which will be permuted
#    indexs = np.arange(len(sp_t))  # the value we're transforming is simple its position
#
#    # Select for permutations
#    indexs_sel = indexs[select]
#    species = sp_t[select]
#
#    # Unique species
#    grouped_s = np.unique(species)
#    # ['A', 'B']
#    grouped_i = [indexs_sel[species == s] for s in grouped_s]
#    # [array([0, 2, 3]), array([5, 7])]
#
#    indiv_permutations = [list(itertools.permutations(i)) for i in grouped_i]
#    indiv_permutations
#    #   [[(0, 2, 3), (0, 3, 2), (2, 0, 3), (2, 3, 0), (3, 0, 2), (3, 2, 0)],
#    #    [(5, 7), (7, 5)]]
#
#    all_perm = list(itertools.product(*indiv_permutations))
#    #    [((0, 2, 3), (5, 7)),
#    #     ((0, 2, 3), (7, 5)),
#    #     ((0, 3, 2), (5, 7)),
#    #     ((0, 3, 2), (7, 5)),
#    #     ((2, 0, 3), (5, 7)),
#    #     ((2, 0, 3), (7, 5)),
#    #     ((2, 3, 0), (5, 7)),
#    #     ((2, 3, 0), (7, 5)),
#    #     ((3, 0, 2), (5, 7)),
#    #     ((3, 0, 2), (7, 5)),
#    #     ((3, 2, 0), (5, 7)),
#    #     ((3, 2, 0), (7, 5))]
#
#    all_perm_flattened = [list(itertools.chain.from_iterable(i)) for i in all_perm]
#
#    indexs_permutations = []
#
#    for perm in all_perm:
#        indexi = indexs.copy()
#        for i, s in enumerate(grouped_s):
#            indexi[select & (sp_t == s)] = perm[i]
#        indexs_permutations.append(indexi)
#        # del s
#        # del i
#
#    # del perm
#
#    indexs_permutations = np.array(indexs_permutations)
#
#    return indexs_permutations
#
#
# def get_values_for_all_permutations(val_t, permutations):
#    # Val_t could be an atoms object
#    return np.array([val_t[perm] for perm in permutations])
#
#
# def get_distances_between_images(imagesi):
#    """Returns distance between each image ie 2norm of d2-d1"""
#
#    spring_lengths = []
#    for j in range(len(imagesi) - 1):
#        spring_vec = imagesi[j + 1].get_positions() - imagesi[j].get_positions()
#        spring_lengths.append(np.linalg.norm(spring_vec))
#    return np.array(spring_lengths)
#
#
# def add_intermediary_images(
#    imagesi, dist_cutoff, interpolate_method="idpp", max_number=100, verbose=False,
# ):
#    """Add additional images inbetween existing ones, purely based on geometry"""
#    # create copy of images
#    imagesi = [at.copy() for at in imagesi]
#    interp_images = []
#    max_dist_images = max(get_distances_between_images(imagesi))
#    for iter in range(max_number):
#        if max_dist_images <= dist_cutoff:
#            print(f"Max distance readched after {iter} iterations")
#            break
#        distances = get_distances_between_images(imagesi)
#        jmax = np.argmax(distances)
#
#        toInterpolate = [imagesi[jmax]]
#        toInterpolate += [toInterpolate[0].copy()]
#        toInterpolate += [imagesi[jmax + 1]]
#
#        neb = NEB(toInterpolate)
#        neb.interpolate(method=interpolate_method, apply_constraint=True)
#
#        interp_images.append([jmax, toInterpolate[1].copy()])
#        # Add images
#        imagesi.insert(jmax + 1, toInterpolate[1].copy())
#        if verbose:
#            print(f"Additional image added at {jmax} with distances {max(distances)}")
#        max_dist_images = max(get_distances_between_images(imagesi))
#
#    return interp_images, imagesi
#
# def get_slice(mystring):
#    return slice(*[{True: lambda n: None, False: int}[x == ''](x) for x in (mystring.split(':') + ['', '', ''])[:3]])
#
# def sort_atom_species(at, mask_slice):
#    cs = at.get_chemical_symbols()
#    order = np.arange(len(at))
#    order[get_slice(mask_slice)] = order[get_slice(mask_slice)][np.argsort(cs[get_slice(mask_slice)])]
#    nat = ase.Atoms([at[i] for i in order], cell=at.get_cell(), pbc=at.get_pbc())
#    return nat
#
#
## permuting atoms
# def get_all_permutations(data, tags):
#    all_perms = [data]
#    for tag in np.unique(tags):
#        # print(tag)
#        # print(all_perms)
#        all_permsi = all_perms.copy()
#        all_perms = []
#        for dataj in all_permsi:
#            for perm in list(itertools.permutations(dataj[tags == tag])):
#                datai = dataj.copy()
#                datai[tags == tag] = perm
#                all_perms.append(datai)
#
#    return np.array(all_perms)
#
# def permute_atoms_to_smallest_distance(atstart, atend, mol_index):
#    indices = np.arange(len(atend))
#    species = np.array(atend.get_chemical_symbols())
#    positions = atend.get_positions()
#
#    data = indices[mol_index:]
#    tags = species[mol_index:]
#
#    permuts = get_all_permutations(data, tags)
#
#    # get distance between all permutations
#    pos_start = atstart.get_positions()[mol_index:]
#    dists = [np.linalg.norm(positions[perm] - pos_start) for perm in permuts]
#    at_new = atend.copy()
#
#    at_new.set_positions(np.concatenate([positions[:mol_index], positions[permuts[np.argmin(dists)]]]))
#    dist = np.min(dists)
#    # print(dist)
#    at_new.info['neb_dist'] = dist
#    return at_new
#
# def permute_trajectory_to_smallest_distance_between_images(traj,mol_index):
#    traj = [at.copy() for at in traj]
#    for i, at in enumerate(traj):
#        if i == 0:
#            continue
#        traj[i] = permute_atoms_to_smallest_distance(traj[i-1], traj[i], mol_index).copy()
#    return traj


def make_neb(initial: Atoms, final: Atoms, images_no: int = 10, idpp: bool = True):
    """Function that creates a trajectory of linear interpolations between two ase atoms objects.

    Args:
        initial (Atoms): Initial image of the neb chain
        final (Atoms): Final image of the neb chain
        images_no (int, optional): Number of images the neb chain. Defaults to 10.
        idpp (bool): Use preoptimization on linear interpolated images. Defaults to True.

    Returns:
        list of Atoms: Interpolated images
    """
    from ase.neb import NEB

    images = [initial]
    for i in range(images_no):
        images.append(initial.copy())
    images.append(final)

    neb = NEB(images, images_no)

    if idpp == True:
        neb.interpolate(method="idpp")
    else:
        neb.interpolate()

    return neb.images


def _swap_atoms_positions(
    atoms: Atoms, target_ind: Union[tuple, list], permutation: Union[tuple, list]
):
    """Internal routine that

    Args:
        atoms (Atoms): _description_
        target_ind (Union[tuple,list]): _description_
        permutation (Union[tuple,list]): _description_

    Returns:
        _type_: _description_
    """

    if not set(target_ind) == set(permutation):
        print("Target_ind and permutation must contain same indicies.")
        return
    else:
        positions = [atoms[i].position.copy() for i in permutation]

        # print(positions)

        for i, j in enumerate(target_ind):
            # print(f'setting position of atom {atoms[j]} from {atoms[j].position} to {positions[i]}', )
            atoms[j].position = positions[i]
    return atoms


def _get_permuted_indices(atoms):
    from collections import OrderedDict
    from itertools import permutations

    species_in_order = [atom.symbol for atom in atoms]
    species = list(OrderedDict.fromkeys(species_in_order))

    mol_atom_list = [i for i, f in enumerate(atoms.arrays["fragments"]) if f == 1]

    permuted_trj = []
    groups = []
    for symbol in species:
        group = [
            a.index for a in atoms if a.symbol == symbol and a.index in mol_atom_list
        ]
        groups.append(group)

    pgroups = []
    for group in groups:
        if len(group) > 0:
            pgroup = []
            for p in permutations(group):
                pgroup.append(p)
            pgroups.append(pgroup)

    return pgroups


def get_swapped_by_species(atoms: Atoms):
    """_summary_

    Args:
        atoms (Atoms): _description_
    """

    from itertools import permutations

    pgroups = _get_permuted_indices(atoms)

    permutations = []
    for pgroup in itertools.product(*pgroups):
        tup = ()
        for t in pgroup:
            tup += t
        permutations.append(tup)

    swapped_trj = []
    for permutation in permutations:
        a = atoms.copy()
        a = _swap_atoms_positions(a, permutations[0], permutation)
        swapped_trj.append(a)

    return swapped_trj
