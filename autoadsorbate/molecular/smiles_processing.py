"""SMILES processing and manipulation functions."""

import re
from typing import List, Tuple
from rdkit import Chem
from rdkit.Chem import AllChem


def get_marked_smiles(smiles: str) -> str:
    """
    Adds marking to a SMILES string to indicate surface attachment points.

    Args:
        smiles (str): The original SMILES string.

    Returns:
        str: The marked SMILES string.
    """
    # This is a placeholder - the actual implementation would depend on the specific marking scheme
    # used in the original codebase
    return f"*{smiles}"


def check_smile(smiles: str) -> bool:
    """
    Validates a SMILES string.

    Args:
        smiles (str): The SMILES string to validate.

    Returns:
        bool: True if the SMILES string is valid, False otherwise.
    """
    mol = Chem.MolFromSmiles(smiles)
    return mol is not None


def remove_canonical_duplicates(smiles_list: List[str]) -> List[str]:
    """
    Removes duplicate SMILES strings by converting to canonical form.

    Args:
        smiles_list (List[str]): List of SMILES strings.

    Returns:
        List[str]: List with duplicates removed.
    """
    canonical_smiles = []
    seen = set()
    
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            canonical = Chem.MolToSmiles(mol, canonical=True)
            if canonical not in seen:
                canonical_smiles.append(smiles)
                seen.add(canonical)
    
    return canonical_smiles


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
    
    for smiles in smiles_list:
        if smiles.startswith("*"):
            top_smiles.append(smiles)
        else:
            nontop_smiles.append(smiles)
    
    return top_smiles, nontop_smiles
