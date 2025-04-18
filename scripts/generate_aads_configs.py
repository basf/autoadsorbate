import numpy as np
from glob import glob
import uuid

from joblib import Parallel, delayed

from ase.build import fcc211
from ase.io import read, write

import sys

sys.path.insert(0, "/gpfs/projects/qm_inorganics/notebook_edvin/git/autoadsorbate/")
from autoadsorbate.autoadsorbate import Fragment, Surface

smiles_file_path = "/gpfs/projects/qm_inorganics/notebook_edvin/aads_paper/relax_Cu111_Cu2_try2/MACE_relax_*.xyz"
slab = fcc211(symbol="Cu", size=(6, 3, 3), vacuum=10)
touch_sphere_size = 2.8
to_init_dict = {"S1": 10, "Cl": 5}

print("preparing slab and Surface")
slab = fcc211(symbol="Cu", size=(3, 3, 2), vacuum=10)
slab.positions += np.array([5, 0, 0])
slab.wrap()
s = Surface(slab, touch_sphere_size=touch_sphere_size)
s.site_df = s.site_df[s.site_df.index.isin([15, 16])]

files = glob(smiles_file_path)
files.sort()


def run_generate(file, slab, touch_sphere_size, to_init_dict):
    try:
        smiles = read(file).info["smiles"]
    except:
        print("Failed at file: ", file)

    try:
        fragment = Fragment(smiles, to_initialize=to_init_dict[smiles[:2]])
        print("Success: working on SMILES: ", smiles)
        print("running get_populated_sites...")

        pop_trj = s.get_populated_sites(
            fragment,
            site_index="all",
            sample_rotation=True,
            mode="heuristic",
            conformers_per_site_cap=5,
            overlap_thr=1.6,
            verbose=True,
        )
        for atoms in pop_trj:
            atoms.info["uid"] = str(uuid.uuid4())

        write(
            f"./generated_{atoms.info['adsorbate_info']['adsorbate_formula']}_{smiles}.xyz",
            pop_trj,
        )
        print(
            f"at file {len(glob('./generated*xyz'))} / {len(files)}. Generated {len(pop_trj)} structures for smiles {smiles}"
        )

    except:
        print("Failed at SMILES: ", smiles)


# for i, file in enumerate(files):
#    run_generate(file, slab, touch_sphere_size, to_init_dict)

Parallel(n_jobs=64)(
    delayed(run_generate)(file, slab, touch_sphere_size, to_init_dict) for file in files
)
