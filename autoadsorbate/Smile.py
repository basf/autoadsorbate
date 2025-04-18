import re
from tempfile import NamedTemporaryFile
from typing import List, Tuple

import numpy as np
from ase import Atoms
from ase.io import read
from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom
from rdkit.Chem.rdForceFieldHelpers import (
    MMFFGetMoleculeForceField,
    MMFFGetMoleculeProperties,
    OptimizeMoleculeConfs,
)

from .utils import rotation_matrix_from_vectors


def conformers_from_smile(
    smiles: str, conformer_count: int = 10, random_seed: int = 0xF00D
) -> List[Atoms]:
    """
    Generates conformers from a SMILES string.

    Args:
        smiles (str): The SMILES string of the molecule.
        conformer_count (int, optional): The number of conformers to generate. Defaults to 10.
        random_seed (int, optional): The random seed for conformer generation. Defaults to 0xf00d.

    Returns:
        List[Atoms]: A list of ASE Atoms objects representing the conformers.
    """
    naked_mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(naked_mol)

    rdDistGeom.EmbedMultipleConfs(
        mol, conformer_count, randomSeed=random_seed
    )  # Generate conformer_count conformers

    conformer_trj = []
    for conf_id in range(conformer_count):
        conformer = mol.GetConformer(conf_id)
        positions = conformer.GetPositions()
        symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
        atoms = Atoms(symbols, positions=positions)
        conformer_trj.append(atoms)

    return conformer_trj


def NEW_create_adsorbates(
    smiles_list: List[str] = ["ClH"], conformations_per_smiles: int = 10
) -> List[Atoms]:
    """
    Generates adsorbates from a list of SMILES strings.

    Args:
        smiles_list (List[str], optional): List of SMILES strings formatted to start with a decorator (Cl, S1S). Defaults to ['ClH'].
        conformations_per_smiles (int, optional): Number of conformations to generate per SMILES string. Defaults to 10.

    Returns:
        List[Atoms]: List of ASE Atoms objects representing the adsorbates.
    """
    trj = []

    for smiles in smiles_list:
        if smiles == "ClH":
            trj.append(Atoms(["Cl", "H"], [[0, 0, 0], [0, 0, 1.5]]))
            continue

        if check_smile(smiles):
            conformer_trj = conformers_from_smile(
                smiles, conformer_count=conformations_per_smiles
            )

            for atoms in conformer_trj:
                atoms = _reset_rotation(atoms)

                if atoms[0].symbol == "S":
                    slide = atoms[1].position - atoms[0].position
                    atoms.positions = atoms.positions - slide * 0.5

                trj.append(atoms)

    return trj


def sort_smiles(smiles_list: List[str]) -> Tuple[List[str], List[str]]:
    """
    Sorts a list of SMILES strings into top and non-top categories based on their prefixes.

    Args:
        smiles_list (List[str]): List of SMILES strings to be sorted.

    Returns:
        Tuple[List[str], List[str]]: Two lists, one for top SMILES and one for non-top SMILES.
    """
    top_smiles = []
    nontop_smiles = []

    for s in smiles_list:
        if s.startswith("Cl"):
            top_smiles.append(s)
        elif s.startswith("S1S"):
            nontop_smiles.append(s)
        else:
            raise ValueError(
                f"SMILES: {s} in smiles_list must be a marked SMILES (starting with S1S or Cl)."
            )

    return top_smiles, nontop_smiles


def atoms_from_smile(smiles: str) -> Atoms:
    """
    Generates an ASE Atoms object from a SMILES string.

    Args:
        smiles (str): The SMILES string of the molecule.

    Returns:
        Atoms: The ASE Atoms object representing the molecule.
    """
    m = Chem.MolFromSmiles(smiles)
    m2 = Chem.AddHs(m)
    AllChem.EmbedMolecule(m2)
    AllChem.MMFFOptimizeMolecule(m2)

    mprop = MMFFGetMoleculeProperties(m2)
    ff = MMFFGetMoleculeForceField(m2, mprop)
    E_cff = OptimizeMoleculeConfs(m2, ff)[0][1]

    with NamedTemporaryFile(suffix=".xyz") as tmp_file:
        print(Chem.MolToXYZBlock(m2), file=open(tmp_file.name, "w"))
        atoms = read(tmp_file.name)

    atoms.pbc = [True, True, True]
    atoms.cell = [[10, 0, 0], [0, 10, 0], [0, 0, 10]]
    atoms.center()

    return atoms


def create_adsorbates(
    smiles: str = "ClH",
    dist: float = 1.1,
    check_dist: bool = True,
    conf_no: int = 10,
    max_tries: int = 20,
) -> List[Atoms]:
    """
    Generates adsorbates from a SMILES string.

    Args:
        smiles (str, optional): Special SMILES formatted to start with a decorator (Cl, S1S). Defaults to 'ClH'.
        dist (float, optional): Minimum distance from the adsorption site plane. Defaults to 1.1 Å.
        check_dist (bool, optional): Check if the generated conformation is far enough from the site plane. Defaults to True.
        conf_no (int, optional): Number of conformations to generate. Defaults to 10.
        max_tries (int, optional): Maximum number of tries for generating a 'far-enough' conformation before giving up. Defaults to 20.

    Returns:
        List[Atoms]: List of CFF optimized molecules.
    """
    trj = []

    if not check_smile(smiles):
        return trj

    for i in range(conf_no):
        min_distance = 0
        n = 1
        while min_distance < dist:
            atoms = atoms_from_smile(smiles)

            if check_dist:
                min_distance = get_closest_atom(atoms)
            else:
                min_distance, dist = 100, 99  # large number to keep all configs

            if n % 5 == 0:
                print(f"Tried {n}, min_distance {min_distance}")
            n += 1
            if n > max_tries:
                print(f"Failed to find configuration in {max_tries} steps")
                break
            if min_distance > dist:
                atoms.info["smile"] = smiles
                atoms = _reset_rotation(atoms)

                if atoms[0].symbol == "S":
                    slide = atoms[1].position - atoms[0].position
                    atoms.positions -= slide * 0.5
                trj.append(atoms)

        print(
            f"{smiles} found {len(trj)} out of {conf_no} configurations that are {round(min_distance, 3)} Å far from the defect plane"
        )

    return trj


def get_nvector(atoms: Atoms) -> np.ndarray:
    """
    Calculates the normal vector for the given atoms.

    Args:
        atoms (Atoms): The ASE Atoms object containing the atoms.

    Returns:
        np.ndarray: The normal vector.
    """
    atom1 = atoms[0]  # first atom
    atom2 = atoms[1]  # second atom

    if atom1.symbol == "Cl":
        nvector = atom2.position - atom1.position
    elif atom1.symbol == "S" and atom2.symbol == "S":
        # Compute normal vector to S-S "defect marker"
        nvector = np.cross(
            np.cross(
                atom1.position - atom2.position,
                atoms.get_center_of_mass() - (atom1.position - atom2.position),
            ),
            atom1.position - atom2.position,
        )
        nvector = nvector / -np.linalg.norm(nvector)
    else:
        nvector = np.array([0, 0, 1])  # Default normal vector

    return nvector


def get_closest_atom(atoms: Atoms) -> float:
    """
    Calculates the minimum distance from the adsorption motif to the closest atom in the molecule.

    Args:
        atoms (Atoms): The ASE Atoms object containing the atoms.

    Returns:
        float: The minimum distance from the adsorption motif to the closest atom.
    """
    from sympy import Plane, Point3D

    if atoms[0].symbol == "S":
        first_molecule_atom_index = 2
    elif atoms[0].symbol == "Cl":
        first_molecule_atom_index = 1
    else:
        raise ValueError(
            "Molecule must be formatted to start with a marker (Cl or S1S)."
        )

    nvector = get_nvector(atoms)
    point = tuple(atoms[0].position)
    plane = Plane(Point3D(point), normal_vector=nvector)

    if len(atoms) > 3:
        dist_list = [
            float(plane.distance(Point3D(tuple(atom.position))))
            for atom in atoms[first_molecule_atom_index:]
        ]
        minimum_distance = min(dist_list)
    else:
        minimum_distance = 10.0

    return minimum_distance


def align_to_z(atoms: Atoms) -> Atoms:
    """
    Aligns the atoms to the z-axis by applying a rotation.

    Args:
        atoms (Atoms): The ASE Atoms object to be aligned.

    Returns:
        Atoms: The aligned ASE Atoms object.
    """
    atomsr = atoms.copy()
    vec1 = get_nvector(atomsr)
    vec2 = np.array([0, 0, 1])
    mat = rotation_matrix_from_vectors(vec1, vec2)

    for atom in atomsr:
        vec = atom.position
        vec_rot = mat.dot(vec)
        atom.position = vec_rot

    return atomsr


def align_to_vector(atoms: Atoms, vector: List[float] = [0, 0, 1]) -> Atoms:
    """
    Aligns the atoms to a given vector by applying a rotation.

    Args:
        atoms (Atoms): The ASE Atoms object to be aligned.
        vector (List[float], optional): The target vector to align the atoms to. Defaults to [0, 0, 1].

    Returns:
        Atoms: The aligned ASE Atoms object.
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


def _reset_position(atoms: Atoms) -> Atoms:
    """
    Resets the position of the atoms so that the first atom is at the origin.

    Args:
        atoms (Atoms): The ASE Atoms object to be shifted.

    Returns:
        Atoms: The shifted ASE Atoms object with the first atom at the origin.
    """
    atomsp = atoms.copy()
    atomsp.positions -= atoms[0].position
    return atomsp


def align_to_xz(atoms: Atoms) -> Atoms:
    """
    Aligns the atoms to the xz-plane by applying rotations.

    Args:
        atoms (Atoms): The ASE Atoms object to be aligned.

    Returns:
        Atoms: The aligned ASE Atoms object.
    """
    if atoms[0].symbol == "Cl":
        return atoms

    if atoms[0].symbol == "S" and atoms[1].symbol == "S":
        angle = 0
        cms = atoms.get_center_of_mass()

        while abs(atoms[1].position[1]) > 0.01 or np.sign(atoms[1].position[1]) != 1:
            atoms.rotate(angle, [0, 0, 1], center=atoms[0].position, rotate_cell=False)
            angle += 0.01

        angle = 0
        while abs(cms[1]) > 0.01 or np.sign(cms[2]) != 1:
            atoms.rotate(angle, [1, 0, 0], center=atoms[0].position, rotate_cell=False)
            angle += 0.01
            cms = atoms.get_center_of_mass()

        if cms[0] < 0:
            # if atoms[1].position[0] < 0:
            atoms.rotate(180, [0, 0, 1], center=atoms[0].position, rotate_cell=False)

    return atoms


def _reset_rotation(atoms: Atoms) -> Atoms:
    """
    Resets the rotation of the atoms by aligning them to the z-axis and xz-plane.

    Args:
        atoms (Atoms): The ASE Atoms object to be aligned.

    Returns:
        Atoms: The aligned ASE Atoms object.
    """
    atomsr = atoms.copy()
    atomsr = align_to_z(atomsr)
    atomsr = _reset_position(atomsr)
    atomsr = align_to_xz(atomsr)
    return atomsr


def get_hvector(atoms: Atoms, site: np.ndarray) -> np.ndarray:
    """
    Calculates the h-vector for the given site in the atoms object.

    Args:
        atoms (Atoms): The ASE Atoms object containing the atoms.
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


def check_smile(smiles: str) -> bool:
    """
    Checks if a SMILES string is valid and meets specific criteria.

    Args:
        smiles (str): The SMILES string to be checked.

    Returns:
        bool: True if the SMILES string is valid and meets the criteria, False otherwise.
    """
    proceed = True
    try:
        m = Chem.MolFromSmiles(smiles)
        m2 = Chem.AddHs(m)
        val0 = m2.GetAtoms()[0].GetExplicitValence()
        val1 = m2.GetAtoms()[1].GetExplicitValence()
        atom0_symbol = m2.GetAtoms()[0].GetSymbol()
        atom1_symbol = m2.GetAtoms()[1].GetSymbol()

        if atom0_symbol == "Cl" and val0 > 1:
            proceed = False
        if atom0_symbol == "S" and val0 > 2:
            proceed = False
        if atom1_symbol == "S" and val1 > 2:
            proceed = False

    except Exception as e:
        print(f"Failed to check SMILES {smiles}: {e}")
        proceed = False

    return proceed


def rotate_mol_to_hvector(atoms: Atoms, site: np.ndarray, mol: Atoms) -> Atoms:
    """
    Rotates the molecule to align with the h-vector of the given site.

    Args:
        atoms (Atoms): The ASE Atoms object containing the slab.
        site (np.ndarray): The coordinates of the site to align the molecule to.
        mol (Atoms): The molecule to be rotated.

    Returns:
        Atoms: The rotated molecule.
    """
    m = mol.copy()
    hvector = get_hvector(atoms, site)
    angle = np.arccos(hvector[0]) / np.pi * 180  # Convert radians to degrees
    m.rotate("z", angle)
    return m


def get_bce() -> list:
    """
    Returns a list of elements that have 2-character symbols.

    Returns:
        list: A list of elements with 2-character symbols.
    """
    from ase.data import chemical_symbols

    list_of_two_character_elements = []
    for symbol in chemical_symbols[1:]:
        if len(symbol) == 2:
            list_of_two_character_elements.append(symbol)
    return list_of_two_character_elements


def rearrange_first_double_bond(smile: str) -> str:
    """
    Rearranges the first double or triple bond in a SMILES string.

    Args:
        smile (str): The SMILES string to be rearranged.

    Returns:
        str: The rearranged SMILES string.
    """
    out_smile = smile
    if smile[1] in ["=", "#"]:
        sm_list = _parse_smile(smile)
        out_smile = sm_list[2] + "(" + smile[1] + smile[0] + ")"
        for s in sm_list[3:]:
            out_smile += s
    if smile[0] in ["O"] and smile[1] not in ["=", "#"]:
        sm_list = _parse_smile(smile)
        out_smile = sm_list[1] + "(" + smile[0] + ")"
        for s in sm_list[2:]:
            out_smile += s
    return out_smile


def rearrange_last_double_bond(smile: str) -> str:
    """
    Rearranges the last double or triple bond in a SMILES string.

    Args:
        smile (str): The SMILES string to be rearranged.

    Returns:
        str: The rearranged SMILES string.
    """
    if smile[-2] in ["=", "#"]:
        smile = smile[:-2] + "(" + smile[-2] + smile[-1] + ")"
    if smile[-1] in ["O"] and smile[-2] not in ["=", "#"]:
        smile = smile[:-1] + "(" + smile[-1] + ")"

    return smile


def rearrange_edge_double_bonds(smile: str) -> str:
    """
    Rearranges the first and last double or triple bonds in a SMILES string.

    Args:
        smile (str): The SMILES string to be rearranged.

    Returns:
        str: The rearranged SMILES string.
    """
    smile = rearrange_first_double_bond(smile)
    smile = rearrange_last_double_bond(smile)
    return smile


def get_keep_together(smile: str) -> List[Tuple[int, int]]:
    """
    Returns tuples with indices of characters in the SMILES string that should remain "next to each other".

    Args:
        smile (str): The SMILES string to be parsed.

    Returns:
        List[Tuple[int, int]]: A list of tuples, each containing the start and end indices of characters that should remain together.
    """
    motifs = {
        "()": r"\(([A-Za-z0-9_=\+\-\*]+)\)",  # Matches content within parentheses
        "[]": r"\[([A-Za-z0-9_=\+\-\*]+)\]",  # Matches content within square brackets
    }

    keep_together = []

    for bracket, motif in motifs.items():
        for match in re.finditer(motif, smile):
            keep_together.append((match.start(), match.end()))

    keep_together.sort(key=lambda tup: tup[0])
    return keep_together


def get_biatomic(smile: str) -> List[Tuple[int, int]]:
    """
    Identifies biatomic motifs in a SMILES string and returns their positions.

    Args:
        smile (str): The SMILES string to be analyzed.

    Returns:
        List[Tuple[int, int]]: A list of tuples representing the start and end positions of biatomic motifs.
    """
    keep_together = []
    for motif in get_bce():
        for i in re.finditer(motif, smile):
            if smile[i.start() - 1] != "(" and smile[i.end()] != ")":
                keep_together.append((i.start(), i.end()))
    return keep_together


def get_group_logical(smile: str) -> List[Tuple[int, int]]:
    """
    Tries to group knot atoms (atoms that have one or multiple side chains) in a SMILES string.

    Args:
        smile (str): The SMILES string to be analyzed.

    Returns:
        List[Tuple[int, int]]: A list of tuples representing the start and end positions of grouped knot atoms.
    """
    keep_together = get_keep_together(smile)

    lst = [index for pair in keep_together for index in pair]

    ls = []
    previous = None
    for i in lst:
        if i != previous:
            ls.append(i)
        else:
            del ls[-1]
        previous = i

    l1 = ls[0::2]
    l2 = ls[1::2]
    ls = list(zip(l1, l2))

    for i, pair in enumerate(ls):
        if smile[pair[0]] == "(":
            ls[i] = (pair[0] - 1, pair[1])

    return ls


def get_parsed_tup(smile: str) -> List[Tuple[int, int]]:
    """
    Parses a SMILES string and returns a list of tuples representing the start and end positions of grouped atoms.

    Args:
        smile (str): The SMILES string to be parsed.

    Returns:
        List[Tuple[int, int]]: A list of tuples representing the start and end positions of grouped atoms.
    """
    keep_together = get_group_logical(smile)
    biatomic = get_biatomic(smile)

    ls = keep_together + biatomic
    ls.sort(key=lambda tup: tup[0])

    skip_ls = set()
    for i in ls:
        skip_ls.update(range(i[0], i[-1]))

    parsed_tup = []
    for i, c in enumerate(smile):
        if i not in skip_ls:
            parsed_tup.append((i, i + 1))
        else:
            for pair in ls:
                if pair[0] == i:
                    parsed_tup.append(pair)
                    break

    return parsed_tup


def parse_smile(smile: str) -> List[str]:
    """
    Parses a SMILES string and corrects for ring closures.

    Args:
        smile (str): The SMILES string to be parsed.

    Returns:
        List[str]: A list of parsed SMILES components with corrected ring closures.
    """
    smile = rearrange_edge_double_bonds(smile)
    smile_parsed = _parse_smile(smile)

    previous = None
    parsed_corrected_for_rings = []

    for s in smile_parsed:
        try:
            int(s[0])
            ss = previous + s
            del parsed_corrected_for_rings[-1]
            parsed_corrected_for_rings.append(ss)
        except ValueError:
            parsed_corrected_for_rings.append(s)
        previous = s

    return parsed_corrected_for_rings


def _parse_smile(smile: str) -> List[str]:
    """
    Parses a SMILES string into its components based on parsed tuples.

    Args:
        smile (str): The SMILES string to be parsed.

    Returns:
        List[str]: A list of parsed SMILES components.
    """
    parsed_smile = []
    parsed_tup = get_parsed_tup(smile)
    for i in parsed_tup:
        parsed_smile.append(smile[i[0] : i[-1]])

    return parsed_smile


def reformat_smile(smile: str) -> str:
    """
    Reformats a SMILES string by parsing it and then concatenating the parsed components.

    Args:
        smile (str): The SMILES string to be reformatted.

    Returns:
        str: The reformatted SMILES string.
    """
    happy_smile = parse_smile(smile)
    reformatted_smile = "".join(happy_smile)
    return reformatted_smile


def flip_smile(smile: str) -> str:
    """
    Flips a SMILES string by reversing the order of its parsed components.

    Args:
        smile (str): The SMILES string to be flipped.

    Returns:
        str: The flipped SMILES string.
    """
    happy_smile = parse_smile(smile)
    sad_smile = "".join(happy_smile[::-1])
    return sad_smile


def attack_backbone(smile: str, backbone_index: int, atom_target: str) -> List[str]:
    """
    Returns a list of sites in the proper order for marking (with Cl).

    Args:
        smile (str): The SMILES string to be parsed.
        backbone_index (int): The index of the backbone to be targeted.
        atom_target (str): The atom to be targeted in the backbone.

    Returns:
        List[str]: A list of sites in the proper order for marking.
    """
    normal_list = parse_smile(smile)
    backbone_target = normal_list[backbone_index]

    if atom_target not in backbone_target:
        print(f"Targeted atom: {atom_target}, not found in target: {backbone_target}")
        return []

    rearranged_list = [backbone_target]

    if normal_list[0:backbone_index]:
        rearranged_list.append("(")
        rearranged_list.extend(normal_list[0:backbone_index][::-1])
        rearranged_list.append(")")

    if normal_list[backbone_index + 1 :]:
        rearranged_list.append("(")
        rearranged_list.extend(normal_list[backbone_index + 1 :])
        rearranged_list.append(")")

    return rearranged_list


def get_reformated_target(site: str) -> List[str]:
    """
    Reformats the target site by replacing specific markers with new markers.

    Args:
        site (str): The target site to be reformatted.

    Returns:
        List[str]: A list of reformatted target sites.
    """
    hacked_site_list = []
    for i, j in get_keep_together(site):
        marker = site[i:j]

        special_case_offset = 0
        swap_marker = ""

        if marker == "(O)":
            swap_marker = "Cl[OH+]"
        elif marker == "(=O)":
            swap_marker = "Cl[O+]="
        elif site[i - 1 : j] == "N(C)":
            swap_marker = "Cl[N+](C)"
            special_case_offset = 1

        new_site_list = [marker] + site.split(marker)

        check_len_site = "".join(new_site_list)

        while len(check_len_site) < len(site):
            check_len_site += marker

        hacked_site = swap_marker + check_len_site[len(marker) + special_case_offset :]
        hacked_site_list.append(hacked_site)

    hacked_site_list = list(set(hacked_site_list))

    if site == "O":
        hacked_site_list.append("Cl[O+]")
    elif site == "N":
        hacked_site_list.append("Cl[N+]")
    elif site == "N1":
        hacked_site_list.append("Cl[N+]1")

    return hacked_site_list


def get_marked_smiles(
    smiles_list: List[str],
    attack_atoms: List[str] = ["O", "N"],
    print_output: bool = False,
) -> List[str]:
    """
    Given a list of SMILES strings of isolated molecules, returns reformatted SMILES strings of equivalent molecules
    with coordination sites marked with Cl. The marked SMILES are returned as a list.
    Target atoms supported are O and N.

    Args:
        smiles_list (List[str]): List of SMILES strings to be marked.
        attack_atoms (List[str], optional): List of target atoms to be marked. Defaults to ['O', 'N'].
        print_output (bool, optional): Whether to print the output for debugging. Defaults to False.

    Returns:
        List[str]: A list of marked SMILES strings.
    """
    out_smiles = []

    for s in smiles_list:
        if print_output:
            print(s)
        parsed_smile = parse_smile(s)
        for i, ST in enumerate(parsed_smile):
            if print_output:
                print(i, ST)
            for attack_atom in attack_atoms:
                if attack_atom not in ST:
                    continue

                site_ls = attack_backbone(s, i, attack_atom)
                site = site_ls[0]
                heads = get_reformated_target(site)

                if print_output:
                    print(heads)

                for head in heads:
                    out_smile = head + "".join(site_ls[1:])
                    out_smiles.append(out_smile)

        if print_output:
            print()

    return out_smiles


def insert_in_str(string: str, insert_str: str = "", pos: int = 0) -> str:
    """
    Inserts a string at a specific position in another string.

    Args:
        string (str): The original string.
        insert_str (str, optional): The string to be inserted. Defaults to ''.
        pos (int, optional): The position at which to insert the string. Defaults to 0.

    Returns:
        str: The modified string with the inserted substring.
    """
    return string[:pos] + insert_str + string[pos:]


def drop_motif(string: str, motif: str) -> List[str]:
    """
    Removes all instances of a specific motif from a string and returns a list of all possibilities.

    Args:
        string (str): The original string.
        motif (str): The motif to be removed.

    Returns:
        List[str]: A list of strings with the motif removed in all possible ways.
    """
    return_strings = []
    splits = string.split(motif)

    ssplits = []
    for s in splits:
        ssplits.append(s)
        ssplits.append(motif)
    del ssplits[-1]

    for i, s in enumerate(ssplits):
        if s == motif:
            return_strings.append("".join(ssplits[:i] + ssplits[i + 1 :]))

    return return_strings


def remove_canonical_duplicates(smile_list: List[str]) -> List[str]:
    """
    Makes a list of canonical SMILES and removes duplicates.

    Args:
        smile_list (List[str]): List of SMILES strings to be processed.

    Returns:
        List[str]: List of unique SMILES strings in canonical form.
    """
    from rdkit.Chem import rdMolHash

    c_smiles = []
    unique_list = []

    for s in smile_list:
        m = Chem.MolFromSmiles(s)
        c_smile = rdMolHash.MolHash(m, rdMolHash.HashFunction.CanonicalSmiles)

        if c_smile not in c_smiles:
            unique_list.append(s)
            c_smiles.append(c_smile)

    return unique_list
