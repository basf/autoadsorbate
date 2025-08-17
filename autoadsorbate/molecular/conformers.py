"""Conformer generation from SMILES strings."""

from typing import List
import numpy as np
from ase import Atoms
from rdkit import Chem
from rdkit.Chem import rdDistGeom


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
