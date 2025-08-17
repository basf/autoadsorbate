"""Plotting and visualization functions."""

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from ase.visualize.plot import plot_atoms

from ..utils.analysis import count_C_next_to_O


def gaussian(x, mu, sig):
    """
    Calculate Gaussian function value.
    
    Args:
        x: Input value
        mu: Mean
        sig: Standard deviation
        
    Returns:
        float: Gaussian function value
    """
    return (
        1.0
        / (np.sqrt(2.0 * np.pi) * sig)
        * np.exp(-np.power((x - mu) / sig, 2.0) / 2)
        / 100
    )


def normalize_energy_values(energy_values, mode="integral"):
    """
    Normalize energy values.
    
    Args:
        energy_values: Array of energy values
        mode: Normalization mode ('integral' or 'max')
        
    Returns:
        np.ndarray: Normalized energy values
    """
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
    """
    Generate Gaussian vector for a single energy value.
    
    Args:
        e: Energy value
        std: Standard deviation
        e_min: Minimum energy
        e_max: Maximum energy
        resolution: Energy resolution
        normalize: Whether to normalize
        normalize_mode: Normalization mode
        
    Returns:
        tuple: (energy_range, energy_values)
    """
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
    """
    Generate Gaussian vectors for multiple energy values.
    
    Args:
        energies: Array of energy values
        std: Standard deviation
        e_min: Minimum energy
        e_max: Maximum energy
        resolution: Energy resolution
        normalize: Whether to normalize
        normalize_mode: Normalization mode
        
    Returns:
        tuple: (energy_range, energy_values)
    """
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
    """
    Generate energy descriptor from a DataFrame slice.
    
    Args:
        df_slice: DataFrame slice
        column: Column name for energy values
        std: Standard deviation
        e_min: Minimum energy
        e_max: Maximum energy
        resolution: Energy resolution
        normalize: Whether to normalize
        normalize_mode: Normalization mode
        
    Returns:
        tuple: (energy_range, energy_values)
    """
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


def docs_plot_conformers(conformer_trajectory, figsize=(15, 10)):
    """
    Plot conformers from a trajectory.
    
    Args:
        conformer_trajectory: List of conformer structures
        figsize: Figure size
        
    Returns:
        matplotlib.figure.Figure: The plot figure
    """
    n_conformers = len(conformer_trajectory)
    n_cols = min(5, n_conformers)
    n_rows = (n_conformers + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    if n_rows == 1:
        axes = [axes] if n_cols == 1 else axes
    else:
        axes = axes.flatten()
    
    for i, conformer in enumerate(conformer_trajectory):
        if i < len(axes):
            plot_atoms(conformer, ax=axes[i])
            axes[i].set_title(f'Conformer {i+1}')
            axes[i].set_xlabel('X (Å)')
            axes[i].set_ylabel('Y (Å)')
    
    # Hide empty subplots
    for i in range(n_conformers, len(axes)):
        axes[i].set_visible(False)
    
    plt.tight_layout()
    return fig


def docs_plot_sites(surface, figsize=(12, 8)):
    """
    Plot adsorption sites on a surface.
    
    Args:
        surface: Surface object with site information
        figsize: Figure size
        
    Returns:
        matplotlib.figure.Figure: The plot figure
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    # Plot the surface atoms
    plot_atoms(surface.atoms, ax=ax)
    
    # Plot adsorption sites
    if hasattr(surface, 'site_dataframe') and surface.site_dataframe is not None:
        sites = surface.site_dataframe
        ax.scatter(sites['x'], sites['y'], c='red', s=100, marker='o', alpha=0.7, label='Adsorption Sites')
        
        # Add site labels
        for idx, site in sites.iterrows():
            ax.annotate(f"{site.get('type', 'site')}", 
                       (site['x'], site['y']), 
                       xytext=(5, 5), textcoords='offset points',
                       fontsize=8, alpha=0.8)
    
    ax.set_xlabel('X (Å)')
    ax.set_ylabel('Y (Å)')
    ax.set_title('Surface with Adsorption Sites')
    ax.legend()
    
    plt.tight_layout()
    return fig


def plot_most_stable(_xdf, relaxed_traj):
    """
    Plot the most stable structures for each combination of H and backbone.
    
    Args:
        _xdf: DataFrame with energy data
        relaxed_traj: List of relaxed trajectory structures
        
    Returns:
        matplotlib.figure.Figure: The plot figure
    """
    fig, axs = plt.subplots(
        ncols=len(_xdf.H.unique()), nrows=len(_xdf.backbone.unique()), figsize=[10, 8]
    )

    _xdf = _xdf.sort_values(by=["H", "backbone"])

    view_atoms = []

    for i, backbone in enumerate(_xdf.backbone.unique()):
        for j, H in enumerate(_xdf.H.unique()):
            ax = axs[i, j]

            df_slice = _xdf[_xdf.H.isin([H]) & _xdf.backbone.isin([backbone])]
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

                plot_atoms(atoms, ax, rotation=("0x,0y,0z"), show_unit_cell=0)
                ax.set_title(df_slice.smiles.iloc[0], size=8)

            ax.set_axis_off()
            ax.set_xlim(half_cell[0] - 3, half_cell[0] + 3)
            ax.set_ylim(half_cell[1] - 3, half_cell[1] + 3)

            view_atoms.append(atoms)

    plt.tight_layout(pad=0.01, w_pad=0.4, h_pad=0.01)
    return fig


def make_hist_plot(_xdf):
    """
    Create histogram plots for energy distributions.
    
    Args:
        _xdf: DataFrame with energy data
        
    Returns:
        matplotlib.figure.Figure: The plot figure
    """
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

            if len(df_slice) > 0:
                sns.histplot(df_slice, x="energy_calibrated", ax=ax, bins=3, kde=True)
                ax.set_title(df_slice.calibrate_keys.values[0], size=8)
                ax.tick_params(axis="x", labelsize=6)
                ax.tick_params(axis="y", labelsize=6)
            else:
                ax.set_axis_off()

    plt.tight_layout(pad=0.01, w_pad=0.4, h_pad=0.01)
    return fig


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
    """
    Create energy heatmap plots.
    
    Args:
        _xdf: DataFrame with energy data
        column: Column name for energy values
        std: Standard deviation for Gaussian
        e_min: Minimum energy
        e_max: Maximum energy
        resolution: Energy resolution
        normalize: Whether to normalize
        return_heatmap: Whether to return heatmap data
        T: Whether to transpose
        cmap: Colormap
        normalize_mode: Normalization mode
        ax: Matplotlib axis
        
    Returns:
        tuple: (axis, heatmap_data) if return_heatmap=True, else (axis, None)
    """
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
        return ax, heat_map
    return ax, None


def get_fragment_center(atoms, fragment_index=1):
    """
    Get the center of a fragment in the atomic structure.
    
    Args:
        atoms: ASE Atoms object
        fragment_index: Index of the fragment
        
    Returns:
        np.ndarray: Center position of the fragment
    """
    # This is a placeholder implementation
    # In the original code, this would use fragment arrays
    return atoms.get_center_of_mass()


def filter_xdf(xdf, relaxed_traj):
    """
    Filter the DataFrame based on certain criteria.
    
    Args:
        xdf: Input DataFrame
        relaxed_traj: Relaxed trajectory
        
    Returns:
        pd.DataFrame: Filtered DataFrame
    """
    # This is a placeholder implementation
    # In the original code, this would apply specific filtering
    return xdf
