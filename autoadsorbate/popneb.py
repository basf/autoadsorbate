import uuid
from tempfile import NamedTemporaryFile
from typing import List

import numpy as np
import pandas as pd
from ase.io import read, write
from ase.mep import NEB, NEBTools
from ase.optimize import BFGS, FIRE


class popNEB:
    """
    A class to perform Nudged Elastic Band (NEB) calculations in
    "popped" segments, by relaxing each image allong the
    initial chain and running NEB chains between local minima.

    Attributes:
    images (List[str]): A list of file paths to the initial images.
    n_ini_chain (int): The number of initial chains.
    n_chain (int): The number of chains.
    max_steps (int): The maximum number of optimization steps.
    f_max (float): The maximum force convergence criterion.
    ini_chain_steps (int): The number of steps for the initial chain.
    compare_mode (str): The mode of comparison for structures. Default is 'cartesian'.
    compare_threshold (float): The threshold for comparison.
    compare_kwargs (dict): Additional keyword arguments for comparison.
    calc (object): The calculator to be used for NEB calculations.
    calc_kwargs (dict): Additional keyword arguments for the calculator.
    store_all_structures (bool): Whether to store all structures during the calculation.

    Methods:
    __init__: Initializes the popNEB class with the provided parameters.
    """

    def __init__(
        self,
        images: List[str] = [],
        n_ini_chain=10,
        n_chain: int = 10,
        max_steps: int = 100,
        f_max: float = 0.05,
        ini_chain_steps: int = 5,
        compare_mode="cartesian",
        compare_threshold=0.3,
        compare_kwargs={},
        calc=None,
        calc_kwargs={},
        store_all_structures=True,
    ):
        self._id = str(uuid.uuid1())
        self.store_all_structures = store_all_structures
        self.images = images
        self.n_chain = n_chain
        self.n_ini_chain = n_ini_chain
        self.max_steps = max_steps
        self.f_max = f_max
        self.ini_chain_steps = ini_chain_steps
        self.compare_mode = compare_mode
        self.compare_threshold = compare_threshold
        self.compare_kwargs = compare_kwargs
        self.images_original = self.images.copy()
        self.calc = calc
        self.calc_kwargs = calc_kwargs
        self.minima_allong_chain = None
        self.neb_path = None
        self.results = {}
        self.tmp_file = NamedTemporaryFile(suffix=".xyz")

    def save(self, path="./", fname="none"):
        write(path + f"{self._id}.xyz", self.result_traj)

    def get_unique_images(self, images=None):
        """
        Saves the result trajectory to a file.

        Parameters:
        path (str): The directory path where the result file will be saved.

        Returns:
        None
        """
        if self.compare_mode.lower() == "cartesian":
            if images == None:
                return compare_cartesian(self.images, self.compare_threshold)
            else:
                return compare_cartesian(images, self.compare_threshold)

        if self.compare_mode.lower() == "default":
            print("Using default SOAP and clustering to determine unique minima!")
            if images == None:
                return compare_soap_default(self.images)
            else:
                unique_trj, self.best_clustering = compare_soap_default(images)
                print(f"Found {len(unique_trj)} unique local minima!")
                return unique_trj

        else:
            print(f"compare_mode: '{self.compare_mode}' not implemented ")

    def interpolate_pair(self, initial, final, n_chain):
        """
        Generates a series of images for the Nudged Elastic Band (NEB) method by interpolating between initial and final structures.

        Parameters:
        initial (object): An ASE atoms object representing the initial structure.
        final (object): An ASE atoms object representing the final structure.
        n_chain (int): The number of intermediate images to generate.

        Returns:
        list: A list of ASE atoms objects representing the interpolated images.
        """
        images = [initial.copy()]
        for i in range(n_chain):
            image = initial.copy()
            images.append(image)
        images.append(final)
        neb = NEB(images)
        neb.interpolate()
        return neb.images

    def run_few_steps(self, initial, final):
        """
        Runs a few steps of the NEB (Nudged Elastic Band) method between the initial and final structures.

        Parameters:
        initial (object): An ASE atoms object representing the initial structure.
        final (object): An ASE atoms object representing the final structure.

        Returns:
        list: A list of ASE atoms objects representing the trajectory after running a few NEB steps.
        """
        dyn_index = "iniNEB"
        few_step_trj = self.run_neb(
            initial,
            final,
            self.ini_chain_steps,
            self.n_ini_chain,
            dyn_index=dyn_index,
            neb_type="iniNEB",
        )
        return few_step_trj

    def run_neb(self, initial, final, max_steps, n_chain, dyn_index, neb_type="NEB"):
        """
        Runs the NEB (Nudged Elastic Band) method between the initial and final structures.

        Parameters:
        initial (object): An ASE atoms object representing the initial structure.
        final (object): An ASE atoms object representing the final structure.
        max_steps (int): The maximum number of optimization steps.
        n_chain (int): The number of intermediate images to generate.
        dyn_index (str): The index for the dynamics observer.
        neb_type (str): The type of NEB calculation. Default is 'NEB'.

        Returns:
        list: A list of ASE atoms objects representing the NEB images after running the NEB method.
        """
        images = self.interpolate_pair(initial, final, n_chain)
        images = self.set_calculator(images, self.calc)

        neb = NEB(images)
        qn = BFGS(neb)
        qn.attach(
            neb_observer,
            dyn=qn,
            dyn_index=dyn_index,
            uid=self._id,
            traj_file=self.tmp_file.name,
            interval=1,
        )
        qn.run(fmax=self.f_max, steps=max_steps)
        if "ini" not in neb_type:
            self.update_result()
        return neb.images

    def relax_images(self, images):
        """
        Relaxes a series of images using the FIRE optimizer.

        Parameters:
        images (list): A list of ASE atoms objects representing the images to be relaxed.

        Returns:
        list: A list of ASE atoms objects representing the relaxed images.
        """
        relax_chain = images.copy()
        self.set_calculator(relax_chain, self.calc)
        for i, atoms in enumerate(relax_chain):
            dyn = FIRE(atoms, trajectory=None)
            dyn.attach(
                relax_observer,
                dyn=dyn,
                dyn_index=i,
                uid=self._id,
                traj_file=self.tmp_file.name,
            )
            dyn.run(fmax=self.f_max, steps=self.max_steps)

        self.update_result()
        return relax_chain

    def set_calculator(self, images):
        """
        Sets the calculator for a series of images.

        Parameters:
        images (list): A list of ASE atoms objects representing the images.
        calculator (object): The calculator to be set for each image.

        Returns:
        list: A list of ASE atoms objects with the calculator set.
        """
        for a in images:
            calc = self.get_calculator()
            a.calc = calc
        return images

    def get_calculator(self):
        """
        Returns the calculator object based on the specified calculator name in the configuration.

        Returns:
        object: The calculator object to be used for NEB calculations.

        Raises:
        ValueError: If the specified calculator name is not supported.
        """
        if self.calc["name"].lower() == "emt":
            from ase.calculators.emt import EMT

            return EMT()

        elif self.calc["name"].lower() == "mace":
            from mace.calculators import mace_mp

            calc = mace_mp(
                model=self.calc["path"],
                dispersion=False,
                default_dtype="float64",
                device="cpu",
                **self.calc_kwargs,
            )
            return calc

        else:
            print(f"not supported: {self.calc}")

    def run(self):
        self.get_minima_allong_chain()
        self.run_neb_from_minima()

    def get_minima_allong_chain(self):
        """
        Finds and relaxes unique images along the NEB chain to identify local minima.

        Returns:
        None
        """
        unique_images = self.get_unique_images()

        c_relaxed = []
        for a1, a2 in zip(unique_images, unique_images[1:]):
            c_images = self.run_few_steps(a1, a2)
            c_relaxed += self.relax_images(c_images)

        self.minima_allong_chain = c_relaxed
        self.minima_allong_chain = self.get_unique_images(c_relaxed)

    def run_neb_from_minima(self):
        """
        Runs the NEB (Nudged Elastic Band) method between local minima identified along the NEB chain.

        Returns:
        None
        """
        if self.minima_allong_chain == None:
            print("run function: get_minima_allong_chain()")
            return

        self.neb_path = []
        dyn_index = 0
        for a1, a2 in zip(self.minima_allong_chain, self.minima_allong_chain[1:]):
            self.neb_path += self.run_neb(
                a1, a2, self.max_steps, self.n_chain, dyn_index
            )
            dyn_index += 1

    def update_result(self):
        """
        Updates the result trajectory and result DataFrame with information from the temporary file.

        Returns:
        None
        """
        result_df = []
        if type(self.tmp_file) == str:
            self.result_traj = read(self.tmp_file, index=":")
        else:
            self.result_traj = read(self.tmp_file.name, index=":")
        for a in self.result_traj:
            result_df.append(a.info["popneb_info"])
        result_df = pd.DataFrame(result_df)
        result_df["converged"] = [
            max_f < self.f_max for max_f in result_df.max_f.values
        ]
        self.result_df = result_df

    def show_neb_report(self, x_scale=1, y_scale=1.5, converged=[True]):
        """
        Generates and displays a NEB (Nudged Elastic Band) report plot.

        Parameters:
        x_scale (float): The scaling factor for the x-axis of the plot. Default is 1.
        y_scale (float): The scaling factor for the y-axis of the plot. Default is 1.5.
        converged (list): A list of boolean values indicating whether to include converged results. Default is [True].

        Returns:
        None
        """
        import matplotlib.pyplot as plt
        from ase.visualize.plot import plot_atoms
        from matplotlib.gridspec import GridSpec

        type = "NEB"
        pdf = self.result_df[
            self.result_df.type.isin([type]) & self.result_df.converged.isin(converged)
        ]

        fig_size_x = len(pdf.image_index.unique()) * x_scale
        fig_size_y = len(pdf.dyn_index.unique()) * y_scale

        fig = plt.figure(layout="constrained", figsize=[fig_size_x, fig_size_y])
        gs = GridSpec(
            len(pdf.dyn_index.unique()) + 1, len(pdf.image_index.unique()), figure=fig
        )

        ax1 = fig.add_subplot(gs[0, :])

        _x = []
        _y = []
        _c = []

        for i in pdf.dyn_index.unique():
            _df = pdf[pdf.dyn_index == i]
            _x += [x * 9 + i * 10 for x in _df.image_index.values / len(_df)]
            _y += list(_df.e.values)
            _c += [i for v in _df.dyn_index.values]

            structures = [self.result_traj[j].copy() for j in _df.index.values]

            for j in _df.image_index.values:
                j = int(j)
                ax = fig.add_subplot(gs[i + 1, j])
                _a = structures[j]
                plot_atoms(_a, ax=ax)
                ax.set_axis_off()

        ax1.scatter(_x, _y, c=_c)
        ax1.set_xticks([])

        fig.suptitle(type)
        fig.set_layout_engine(layout="tight")

        plt.show()

    def show_relax_report(self, x_scale=1, y_scale=4):
        """
        Generates and displays a relaxation report plot.

        Parameters:
        x_scale (float): The scaling factor for the x-axis of the plot. Default is 1.
        y_scale (float): The scaling factor for the y-axis of the plot. Default is 4.

        Returns:
        None
        """
        import matplotlib.pyplot as plt
        from ase.visualize.plot import plot_atoms
        from matplotlib.gridspec import GridSpec

        type = "RELAX"
        pdf = self.result_df[self.result_df.type.isin([type])]

        fig_size_x = len(pdf.dyn_index.unique()) * x_scale
        fig_size_y = 1 * y_scale

        fig = plt.figure(layout="constrained", figsize=[fig_size_x, fig_size_y])
        gs = GridSpec(2, len(pdf.dyn_index.unique()), figure=fig)

        ax1 = fig.add_subplot(gs[0, :])

        cmap = {}
        for i, cluster in enumerate(self.best_clustering):
            cmap[i] = cluster

        _x = []
        _y = []
        _c = []
        structures = []

        for i in pdf.dyn_index.unique():
            _df = pdf[pdf.dyn_index == i]
            _x += [x * 9 + i * 10 for x in _df.nsteps / len(_df)]
            _y += list(_df.e.values)
            _c += [cmap[v] for v in _df.dyn_index.values]

            structure_inds = _df[_df.converged == True].index.values
            if len(structure_inds) == 1:
                structure_i = structure_inds[0]

            structures.append(self.result_traj[structure_i].copy())

        ax1.scatter(_x, _y, c=_c)
        ax1.set_xticks([])

        axs = []
        for i in pdf.dyn_index.unique():
            ax = fig.add_subplot(gs[1, i])
            axs.append(ax)
            plot_atoms(structures[i], ax=axs[i])
            ax.set_axis_off()

        fig.suptitle(type)
        fig.set_layout_engine(layout="tight")

        plt.show()


def get_all_speacies(traj):
    """
    Extracts all unique atomic species from a trajectory.

    Parameters:
    traj (list): A list of ASE atoms objects representing the trajectory.

    Returns:
    list: A list of unique atomic species present in the trajectory.
    """
    all_syms = []
    for atoms in traj:
        all_syms += [atom.symbol for atom in atoms]
    all_syms = list(set(all_syms))
    return all_syms


def get_flattened_soap_desc(soap, traj):
    """
    Generates flattened SOAP (Smooth Overlap of Atomic Positions) descriptors for a trajectory.

    Parameters:
    soap (object): A SOAP descriptor object.
    traj (list): A list of ASE atoms objects representing the trajectory.

    Returns:
    numpy.ndarray: A 2D array where each row is a flattened SOAP descriptor for an atomic structure in the trajectory.
    """
    flat_soap = []
    for atoms in traj:
        flat_soap.append(soap.create(atoms).flatten())
    return np.array(flat_soap)


def get_default_soap(traj):
    """
    Creates a default SOAP (Smooth Overlap of Atomic Positions) descriptor object for a trajectory.

    Parameters:
    traj (list): A list of ASE atoms objects representing the trajectory.

    Returns:
    object: A SOAP descriptor object configured with the species present in the trajectory.
    """
    from dscribe.descriptors import SOAP

    species = get_all_speacies(traj)
    soap = SOAP(
        species=species,
        r_cut=4.5,
        n_max=8,
        l_max=4,
        sigma=0.5,
        periodic=False,
        rbf="gto",
        crossover=True,
    )
    return soap


def get_best_clustering(X):
    """
    Finds the best clustering for the given data using KMeans and silhouette score.

    Parameters:
    X (numpy.ndarray): A 2D array where each row represents a data point.

    Returns:
    numpy.ndarray: An array of cluster labels for the best clustering.
    """
    from sklearn.cluster import KMeans
    from sklearn.metrics import silhouette_score

    model_sizes = [k for k in range(3, len(X))]

    kmeans_per_k = [KMeans(n_clusters=k, random_state=42).fit(X) for k in model_sizes]
    # inertias = [model.inertia_ for model in kmeans_per_k]
    y = [model.fit_predict(X) for model in kmeans_per_k]

    silhouette_scores = [silhouette_score(X, model.labels_) for model in kmeans_per_k]
    best_clustering_ind = silhouette_scores.index(
        max(silhouette_scores)
    )  # index shift due to cut list
    best_clustering = y[best_clustering_ind]
    return best_clustering


def get_best_soap_clustering(traj):
    """
    Finds the best clustering for a trajectory using SOAP (Smooth Overlap of Atomic Positions) descriptors.

    Parameters:
    traj (list): A list of ASE atoms objects representing the trajectory.

    Returns:
    numpy.ndarray: An array of cluster labels for the best clustering.
    """
    soap = get_default_soap(traj)
    flat_soap = get_flattened_soap_desc(soap, traj)
    return get_best_clustering(flat_soap)


def compare_soap_default(traj):
    if len(traj) <= 2:
        return traj
    unique_trj = []
    best_soap_clustering = get_best_soap_clustering(traj)
    # print(best_soap_clustering)
    included_clusters = []
    for i, a in enumerate(traj):
        if best_soap_clustering[i] not in included_clusters:
            unique_trj.append(a.copy())
            included_clusters.append(best_soap_clustering[i])
    return unique_trj, best_soap_clustering


def compare_cartesian(self, images, compare_threshold):
    """
    Compares SOAP (Smooth Overlap of Atomic Positions) descriptors for a trajectory and returns unique structures.

    Parameters:
    traj (list): A list of ASE atoms objects representing the trajectory.

    Returns:
    tuple: A tuple containing:
        - unique_trj (list): A list of unique ASE atoms objects based on SOAP clustering.
        - best_soap_clustering (numpy.ndarray): An array of cluster labels for the best clustering.
    """
    unique_index = [0]
    for i, atoms in enumerate(images):
        deltas = [
            np.sum(np.abs(atoms.positions - images[j].positions)) for j in unique_index
        ]
        if all([d > compare_threshold for d in deltas]):
            unique_index.append(i)
    return [images[i].copy() for i in unique_index]


def neb_observer(dyn, dyn_index, neb_type="NEB", uid=None, traj_file="default"):
    """
    Observes and records the state of the NEB (Nudged Elastic Band) calculation at each step.

    Parameters:
    dyn (object): The ASE dynamics object containing the NEB calculation.
    dyn_index (int): The index for the dynamics observer.
    neb_type (str): The type of NEB calculation. Default is 'NEB'.
    uid (str or None): A unique identifier for the NEB calculation. Default is None.
    traj_file (str): The file path to save the trajectory. Default is 'default'.

    Returns:
    None
    """
    neb = dyn.atoms
    epot = []

    nebtools = NEBTools(neb.images)
    max_f_chain = nebtools.get_fmax()

    for i, a in enumerate(neb.images):
        e = a.get_potential_energy()
        forces_popneb = a.get_forces()
        max_f = np.max(np.linalg.norm(forces_popneb, axis=1))

        converged = 0
        last_step = 0
        if max_f < dyn.fmax:
            converged = 1
        if max_f < dyn.fmax or dyn.nsteps >= dyn.max_steps - 1:
            last_step = 1
        #        converged = max_f_chain < dyn.fmax

        a.info["popneb_info"] = {
            "type": neb_type,
            "e": e,
            "max_f": max_f_chain,
            #            'chain_max_f': max_f_chain,
            #            'max_f':max_f,
            "image_index": i,
            "dyn_index": dyn_index,
            "nsteps": dyn.nsteps,
            "converged": converged,
            "last_step": last_step,
        }
        if uid != None:
            a.info["popneb_info"]["_id"] = uid

        a.arrays["forces_popneb"] = forces_popneb

        if traj_file == "default":
            write(f"./{uid}.xyz", a, append=True)
        else:
            write(traj_file, a, append=True)
    # out_traj.append([a for a in neb.images])


def relax_observer(dyn, dyn_index, uid=None, traj_file="default"):
    """
    Observes and records the state of the relaxation calculation at each step.

    Parameters:
    dyn (object): The ASE dynamics object containing the relaxation calculation.
    dyn_index (int): The index for the dynamics observer.
    uid (str or None): A unique identifier for the relaxation calculation. Default is None.
    traj_file (str): The file path to save the trajectory. Default is 'default'.

    Returns:
    None
    """
    a = dyn.atoms

    e = a.get_potential_energy()
    forces_popneb = a.get_forces()
    max_f = np.max(np.linalg.norm(forces_popneb, axis=1))

    converged = 0
    last_step = 0
    if max_f < dyn.fmax:
        converged = 1
    if max_f < dyn.fmax or dyn.nsteps == dyn.max_steps - 1:
        last_step = 1

    a.info["popneb_info"] = {
        "type": "RELAX",
        "e": e,
        "max_f": max_f,
        "dyn_index": dyn_index,
        "nsteps": dyn.nsteps,
        "converged": converged,
        "last_step": last_step,
    }
    if uid != None:
        a.info["popneb_info"]["_id"] = uid

    a.arrays["forces_popneb"] = forces_popneb

    if traj_file == "default":
        write(f"./{uid}.xyz", a, append=True)
    else:
        write(traj_file, a, append=True)


#    out_traj.append(a)
