"""Fragment class for handling molecular fragments."""

from typing import List, Union
import numpy as np
from ase import Atoms

from ..molecular.conformers import conformers_from_smile
from ..utils.sorting import get_sorted_by_snap_dist


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
        self,
        i: Union[int, float],
        n_vector: np.ndarray = np.array([0, 0, 1]),
        rot_deg: float = 0
    ) -> Atoms:
        """
        Returns a copy of the i-th conformer, aligned and rotated as specified.

        Args:
            i (int): The index of the conformer to retrieve.
            n_vector (np.ndarray, optional): The normal vector for rotation. Defaults to [0, 0, 1].
            rot_deg (float, optional): The rotation angle in degrees. Defaults to 0.

        Returns:
            Atoms: The aligned and rotated conformer.
        """
        from ..molecular.alignment import _reset_position, _reset_rotation
        
        if i >= len(self.conformers):
            raise ValueError(f"Conformer index {i} out of range. Max index: {len(self.conformers) - 1}")
        
        conformer = self.conformers[i].copy()
        
        if not self.conformers_aligned[i]:
            conformer = _reset_position(conformer)
            conformer = _reset_rotation(conformer, n_vector, rot_deg)
            self.conformers_aligned[i] = True
        
        return conformer
