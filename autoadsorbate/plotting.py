import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from ase.visualize.plot import plot_atoms

from autoadsorbate.utils import count_C_next_to_O


def gaussian(x, mu, sig):
    return (
        1.0
        / (np.sqrt(2.0 * np.pi) * sig)
        * np.exp(-np.power((x - mu) / sig, 2.0) / 2)
        / 100
    )


def normalize_energy_values(energy_values, mode="intergal"):
    if mode.lower() == "integral":
        energy_values = energy_values / np.sum(energy_values)
    elif mode.lower() == "max":
        energy_values = energy_values / np.max(energy_values)
    else:
        raise ValueError("normalize_energy_values supports modes: 'integral', 'max'")
    return energy_values


def get_gaussian_vector(
    e,
    std=0.05,
    e_min=-0.2,
    e_max=3,
    resolution=0.01,
    normalize=False,
    normalize_mode="integral",
):
    energy_range = np.linspace(e_min, e_max, int((e_max - e_min) / resolution))
    energy_values = np.zeros(len(energy_range))
    for i, energy in enumerate(energy_range):
        energy_values[i] = gaussian(energy, e, std)

    if normalize:
        energy_values = normalize_energy_values(energy_values, mode=normalize_mode)
    return energy_range, energy_values


def get_gaussian_vectors(
    energies,
    std=0.05,
    e_min=-0.2,
    e_max=3,
    resolution=0.01,
    normalize=True,
    normalize_mode="integral",
):
    energy_range = np.linspace(e_min, e_max, int((e_max - e_min) / resolution))
    energy_values = np.zeros(len(energy_range))
    for e in energies:
        energy_values += get_gaussian_vector(
            e, std=std, e_min=e_min, e_max=e_max, resolution=resolution
        )[1]

    if normalize:
        energy_values = normalize_energy_values(energy_values, mode=normalize_mode)

    return energy_range, energy_values


def energy_descriptor_from_slice(
    df_slice,
    column="energy_calibrated",
    std=0.05,
    e_min="auto",
    e_max="auto",
    resolution="auto",
    normalize=True,
    normalize_mode="integral",
):
    if e_min == "auto":
        e_min = df_slice[column].min() - 5 * std

    if e_max == "auto":
        e_max = df_slice[column].max() + 5 * std

    if resolution == "auto":
        resolution = std / 7

    energy_range, energy_values = get_gaussian_vectors(
        df_slice[column].values,
        std=std,
        e_min=e_min,
        e_max=e_max,
        resolution=resolution,
        normalize=normalize,
        normalize_mode=normalize_mode,
    )

    return energy_range, energy_values


def filter_xdf(xdf, relaxed_traj):
    _xdf = xdf[
        (xdf["energy"] > -100)
        & (xdf["energy"] < 40)
        & (xdf.bond_change == 0)
        & xdf.backbone_formula.isin(["C", "C2", "C2O", "CO", "CO2", "O", "O2"])
        # & (xdf.origin == 'aads')
        # & (xdf.H < 7)
    ]

    if "array_from_ocp" in _xdf.columns:
        _xdf.pop("array_from_ocp")

    C_bonds_O = []
    for i in _xdf.traj_index.values:
        sx = count_C_next_to_O(relaxed_traj[i])
        C_bonds_O.append(sx)
    _xdf["C_bonds_O"] = C_bonds_O

    backbone = []
    for i in _xdf.index.values:
        f = _xdf.backbone_formula.loc[i]
        if f == "C2O" and (_xdf.loc[i]["C_bonds_O"] == 2):
            f = "COC"
        backbone.append(f)
    _xdf["backbone"] = backbone

    # set H_max for each backbone manually
    map_H_max = {"C": 3, "C2": 6, "C2O": 6, "CO": 4, "CO2": 4, "O": 2, "O2": 2}
    map_backbone = dict(
        [(v, i) for i, v in enumerate(sorted(list(_xdf.backbone.unique())))]
    )

    # make plot float for plotting energy verus multiple in 2d
    map_origin = {"aads": 0, "ocp": 1}
    _xdf["plot_float"] = [
        map_backbone[_xdf.backbone.values[i]]
        + _xdf.H.values[i] / 20
        + map_origin[_xdf.origin.values[i]] * 0.5
        for i in range(len(_xdf))
    ]
    _xdf["plot_float"] = _xdf["plot_float"] * 0.2
    _xdf["H_max"] = [
        map_H_max[_xdf.backbone_formula.values[i]] for i in range(len(_xdf))
    ]
    _xdf = _xdf[_xdf["H"] <= _xdf["H_max"]]

    # make reference energy
    _xdf["calibrate_keys"] = _xdf.backbone + "-H" + _xdf.H.astype(int).astype(str)

    set_zero_dict = {}

    for k, v in dict(_xdf.groupby(["backbone", "H"]).energy.min()).items():
        key = f"{k[0]}-H{int(k[1])}"
        set_zero_dict[key] = v

    group_ref_energy = []
    for k in _xdf["calibrate_keys"]:
        group_ref_energy.append(set_zero_dict[k])

    _xdf["group_ref_energy"] = group_ref_energy
    _xdf["energy_calibrated"] = _xdf["energy"] - _xdf["group_ref_energy"]

    print(f"remaining values in DF: {len(_xdf)}")
    return _xdf


def center_fragment_in_cell(atoms, fragment_inds):
    a = atoms.copy()
    fragment = atoms.copy()[
        [
            atom.index
            for atom in atoms
            if atoms.arrays["fragments"][atom.index] in fragment_inds
        ]
    ]
    fragment_center = fragment.get_center_of_mass()
    fragment_center[2] = 0
    a.positions += -fragment_center + (a.cell[0] + a.cell[1]) / 2
    a.wrap()
    return a


def get_fragment_center(atoms, fragment_index):
    a = atoms.copy()
    a = a[
        [
            atom.index
            for atom in atoms
            if a.arrays["fragments"][atom.index] == fragment_index
        ]
    ]
    center = []
    for i in [0, 1, 2]:
        center.append(
            (np.max(a.positions[:, i]) - np.min(a.positions[:, i])) * 0.5
            + np.min(a.positions[:, i])
        )
    return np.array(center)


def plot_most_stable(_xdf, relaxed_traj):
    fig, axs = plt.subplots(
        ncols=len(_xdf.H.unique()), nrows=len(_xdf.backbone.unique()), figsize=[10, 8]
    )

    _xdf = _xdf.sort_values(by=["H", "backbone"])

    view_atoms = []

    for i, backbone in enumerate(_xdf.backbone.unique()):
        for j, H in enumerate(_xdf.H.unique()):
            ax = axs[i, j]

            df_slice = _xdf[_xdf.H.isin([H]) & _xdf.backbone.isin([backbone])]
            # df_slice.sort_values(by=['energy', 'backbone', 'H'],ascending=True)
            df_slice = df_slice[df_slice.energy == df_slice.energy.min()]

            if len(df_slice) > 0:
                e = np.round(df_slice.iloc[0].energy, 2)
                origin = df_slice.iloc[0].origin
                traj_index = df_slice.iloc[0].traj_index

                atoms = relaxed_traj[traj_index].copy()
                atoms_center = get_fragment_center(atoms, fragment_index=1)
                atoms_center[2] = 0
                half_cell = atoms.cell[1] * 0.5 + atoms.cell[0] * 0.5
                atoms.positions += -atoms_center + half_cell
                atoms.wrap()
                # atoms.positions -= half_cell

                plot_atoms(atoms, ax, rotation=("0x,0y,0z"), show_unit_cell=0)
                # ax.set_title(atoms.info['adsorbate_info']['smiles'], size=8)
                ax.set_title(df_slice.smiles.iloc[0], size=8)

            ax.set_axis_off()
            # x = cell[0][0] + cell[1][0]
            # y = cell[0][1] + cell[1][1]
            ax.set_xlim(half_cell[0] - 3, half_cell[0] + 3)
            ax.set_ylim(half_cell[1] - 3, half_cell[1] + 3)

            view_atoms.append(atoms)

    # fig.set_layout_engine(layout='tight')
    plt.tight_layout(pad=0.01, w_pad=0.4, h_pad=0.01)


def make_hist_plot(_xdf):
    fig, axs = plt.subplots(
        ncols=len(_xdf.H.unique()),
        nrows=len(_xdf.backbone.unique()),
        figsize=[10, 10],
        sharex=True,
        sharey=True,
    )

    _xdf = _xdf.sort_values(by=["H", "backbone"], ascending=True)

    view_atoms = []

    for i, backbone in enumerate(_xdf.backbone.unique()):
        for j, H in enumerate(_xdf.H.unique()):
            ax = axs[i, j]

            df_slice = _xdf[_xdf.H.isin([H]) & _xdf.backbone.isin([backbone])]
            # df_slice.sort_values(by=['energy', 'backbone', 'H'],ascending=True)
            # df_slice=df_slice[df_slice.energy==df_slice.energy.min()]

            if len(df_slice) > 0:
                sns.histplot(df_slice, x="energy_calibrated", ax=ax, bins=3, kde=True)
                ax.set_title(df_slice.calibrate_keys.values[0], size=8)
                # ax.set_ylim(0, 150)
                ax.tick_params(axis="x", labelsize=6)
                ax.tick_params(axis="y", labelsize=6)

            else:
                ax.set_axis_off()
            # ax.set_xlim(1, 6)

            # x = cell[0][0] + cell[1][0]
            # y = cell[0][1] + cell[1][1]
            # ax.set_xlim(x+2, x+9)
            # ax.set_ylim(y-3, y+3)

            # view_atoms.append(atoms)

    # fig.set_layout_engine(layout='tight')
    plt.tight_layout(pad=0.01, w_pad=0.4, h_pad=0.01)


def plot_energy_heatmap(
    _xdf,
    column,
    std,
    e_min,
    e_max,
    resolution,
    normalize,
    return_heatmap=False,
    T=False,
    cmap="viridis",
    normalize_mode="max",
    ax=None,
):
    heat_map = []
    yticklabels = []

    for i, backbone in enumerate(_xdf.backbone.unique()):
        for j, H in enumerate(_xdf.H.unique()):
            df_slice = _xdf[_xdf.H.isin([H]) & _xdf.backbone.isin([backbone])]
            if len(df_slice) > 0:
                v = energy_descriptor_from_slice(
                    df_slice,
                    column=column,
                    std=std,
                    e_min=e_min,
                    e_max=e_max,
                    resolution=resolution,
                    normalize=normalize,
                    normalize_mode=normalize_mode,
                )
                heat_map.append(v[1])
                yticklabels.append(df_slice.calibrate_keys.values[0])
    heat_map = np.array(heat_map)

    xticklabels = []
    wanted_labels = np.arange(-10, 10, 0.4)
    for i, e in enumerate(v[0]):
        if any(np.abs(e - wanted_labels) < 1e-2):
            label = str(np.round(e, 1))
            if label not in xticklabels + ["-0.0"]:
                xticklabels.append(label)
            else:
                xticklabels.append("")

        else:
            xticklabels.append("")

    if ax == None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    if T == False:
        sns.heatmap(
            heat_map,
            xticklabels=xticklabels,
            yticklabels=yticklabels,
            cbar=False,
            ax=ax,
        )
        for i in range(heat_map.shape[0] + 1):
            ax.axhline(i, color="white", lw=2)

    else:
        ax = sns.heatmap(
            heat_map.T,
            xticklabels=yticklabels,
            yticklabels=xticklabels,
            cbar=False,
            cmap=cmap,
            ax=ax,
        )
        for i in range(heat_map.shape[1] + 1):
            ax.axvline(i, color="white", lw=0.5)

    ax.tick_params(axis="both", which="both", length=0)
    ax.tick_params(axis="x", labelsize=6)
    ax.tick_params(axis="y", labelsize=6)
    ax.invert_yaxis()

    if return_heatmap:
        return heat_map
