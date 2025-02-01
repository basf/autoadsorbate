import time
import os
from glob import glob
from joblib import Parallel, delayed

import pandas as pd
import numpy as np

import ase
from ase.optimize import BFGS
from ase import units
from ase.constraints import FixAtoms

from mace.calculators import mace_mp
from pathlib import Path

import sys
sys.path.insert(0, '/gpfs/projects/qm_inorganics/notebook_edvin/git/autoadsorbate/')
from autoadsorbate.utils freeze_atoms, save_relax_config, relaxatoms

rcut=40
steps = 15
prefix = 'MACE_prerelax'
freeze_bottom=False

#traj = ase.io.read('./traj.xyz', index = ':')

print('Reading generated traj . . . ')
traj = []
for file in glob('/gpfs/projects/qm_inorganics/notebook_edvin/aads_paper/RELAX/relax_Cu211_get_population/generate/generated_*.xyz'):
    traj += ase.io.read(file, index=':')
print(f'Done reading generated traj . . . found {len(traj)}')

models = '/gpfs/projects/qm_inorganics/notebook_edvin/mace_paper/revision_nature_09012025_mpa-0/mace-mpa-0-medium.model'

print(models)

macemp = mace_mp(model=models,dispersion=True,device="cuda",dispersion_cutoff=rcut* units.Bohr, default_dtype="float64", enable_cueq=True) # return ASE calculator

for atoms in traj:
    c = FixAtoms(mask=[atom.symbol == 'Cu'  for atom in atoms]) 
    atoms.set_constraint(c)

def main():
    for atoms in traj:
        relaxatoms(atoms, macemp, prefix, steps=steps, freeze_bottom=freeze_bottom)
    
if __name__ == "__main__":
    main()
