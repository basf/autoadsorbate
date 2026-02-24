"""Visualise adsorption sites on the Ni/Ru slab, colour-coded by SOAP cluster."""
import copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle
from ase.build import fcc111
from autoadsorbate import Surface
from autoadsorbate.utils import _compute_soap_site_vectors, _soap_similarity_matrix
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

# ── Build Ni/Ru slab ──────────────────────────────────────────────
slab = fcc111("Ni", (3, 3, 3), periodic=True, vacuum=10)
top_z = slab.positions[:, 2].max()
for atom in slab:
    if abs(atom.position[2] - top_z) < 0.1:
        atom.symbol = "Ru"
        break

s = Surface(slab)
s_ase = copy.deepcopy(s)
s_soap = copy.deepcopy(s)
s_ase.sym_reduce(method="ase")
s_soap.sym_reduce(method="soap", similarity_threshold=0.99)
ase_idx = set(s_ase.site_df.index)
soap_idx = set(s_soap.site_df.index)

# ── SOAP cluster labels ───────────────────────────────────────────
vecs = _compute_soap_site_vectors(s.atoms, s.site_df)
sim = _soap_similarity_matrix(vecs)
dist = np.clip(1.0 - sim, 0.0, 2.0)
np.fill_diagonal(dist, 0.0)
Z = linkage(squareform(dist, checks=False), method="average")
labels = fcluster(Z, t=0.01, criterion="distance")

coords = np.array(s.site_df["coordinates"].tolist())
connectivity = np.array(s.site_df["connectivity"].tolist())

# ── Colour map: one colour per SOAP cluster ──────────────────────
n_clusters = len(set(labels))
cmap = plt.colormaps.get_cmap("tab10")
cluster_ids = sorted(set(labels))
colour_map = {cl: cmap(i / max(n_clusters - 1, 1)) for i, cl in enumerate(cluster_ids)}

# Build cluster descriptions
cluster_desc = {}
for cl in cluster_ids:
    mask = labels == cl
    idx_in_cl = np.where(mask)[0]
    row = s.site_df.iloc[idx_in_cl[0]]
    conn = row["connectivity"]
    formula = str(row["site_formula"])
    site_type = {1: "atop", 2: "bridge", 3: "hollow"}.get(conn, f"conn{conn}")
    cluster_desc[cl] = f"{site_type} {formula}"

# ── Helper: draw slab atoms as circles (top-down, x-y plane) ────
atom_colours = {"Ni": "#8CBE6E", "Ru": "#D4534B"}
atom_radius = {"Ni": 0.9, "Ru": 0.95}

def draw_slab(ax, atoms):
    """Draw slab atoms as filled circles, top layer on top."""
    pos = atoms.positions
    symbols = atoms.get_chemical_symbols()
    # Sort by z so that top-layer atoms are drawn last (on top)
    order = np.argsort(pos[:, 2])
    for i in order:
        z_frac = (pos[i, 2] - pos[:, 2].min()) / (pos[:, 2].max() - pos[:, 2].min() + 1e-9)
        alpha = 0.25 + 0.75 * z_frac  # deeper atoms more transparent
        r = atom_radius.get(symbols[i], 0.8)
        c = atom_colours.get(symbols[i], "#AAAAAA")
        circle = Circle(
            (pos[i, 0], pos[i, 1]), r,
            facecolor=c, edgecolor="black", linewidth=0.4,
            alpha=alpha, zorder=2 + z_frac,
        )
        ax.add_patch(circle)

    # Draw unit cell outline
    cell = atoms.cell
    corners = np.array([
        [0, 0], [cell[0, 0], cell[0, 1]],
        [cell[0, 0] + cell[1, 0], cell[0, 1] + cell[1, 1]],
        [cell[1, 0], cell[1, 1]], [0, 0],
    ])
    ax.plot(corners[:, 0], corners[:, 1], "k-", lw=0.8, zorder=1)

# ── Figure ────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for ax_i, (title, highlight_idx) in enumerate([
    (f"ASE representatives ({len(ase_idx)} sites)", ase_idx),
    (f"SOAP representatives ({len(soap_idx)} sites)", soap_idx),
]):
    ax = axes[ax_i]

    # Draw slab atoms in the same coordinate system
    draw_slab(ax, s.atoms)

    # Plot ALL sites as small grey dots
    ax.scatter(
        coords[:, 0], coords[:, 1],
        s=15, c="lightgrey", edgecolors="grey", linewidths=0.3,
        zorder=5, label="_nolegend_",
    )

    # Overlay highlighted representative sites
    for cl in cluster_ids:
        mask = labels == cl
        for site_i in np.where(mask)[0]:
            df_idx = s.site_df.index[site_i]
            if df_idx not in highlight_idx:
                continue
            conn = connectivity[site_i]
            marker = {1: "o", 2: "s", 3: "^"}.get(conn, "D")
            ax.scatter(
                coords[site_i, 0], coords[site_i, 1],
                s=120, c=[colour_map[labels[site_i]]],
                edgecolors="black", linewidths=1.0,
                marker=marker, zorder=10,
            )
            # Label with index
            ax.annotate(
                str(df_idx), (coords[site_i, 0], coords[site_i, 1]),
                textcoords="offset points", xytext=(5, 5),
                fontsize=7, fontweight="bold", zorder=11,
            )

    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_xlabel("x (Å)")
    ax.set_ylabel("y (Å)")
    ax.set_aspect("equal")
    ax.autoscale_view()
    # Add a small margin
    pad = 1.0
    ax.set_xlim(coords[:, 0].min() - pad, coords[:, 0].max() + pad)
    ax.set_ylim(coords[:, 1].min() - pad, coords[:, 1].max() + pad)

# ── Legend ─────────────────────────────────────────────────────────
legend_handles = []
# Cluster colours
for cl in cluster_ids:
    legend_handles.append(
        Line2D([0], [0], marker="o", color="w",
               markerfacecolor=colour_map[cl], markeredgecolor="black",
               markersize=8, label=f"Cluster {cl}: {cluster_desc[cl]}")
    )
# Marker shapes for connectivity
legend_handles.append(Line2D([0], [0], marker="o", color="w", markerfacecolor="grey",
                             markersize=8, label="atop (conn=1)"))
legend_handles.append(Line2D([0], [0], marker="s", color="w", markerfacecolor="grey",
                             markersize=8, label="bridge (conn=2)"))
legend_handles.append(Line2D([0], [0], marker="^", color="w", markerfacecolor="grey",
                             markersize=8, label="hollow (conn=3)"))
# Slab atoms
legend_handles.append(Line2D([0], [0], marker="o", color="w", markerfacecolor="#8CBE6E",
                             markeredgecolor="black", markersize=10, label="Ni atom"))
legend_handles.append(Line2D([0], [0], marker="o", color="w", markerfacecolor="#D4534B",
                             markeredgecolor="black", markersize=10, label="Ru atom"))

fig.legend(handles=legend_handles, loc="lower center", ncol=4, fontsize=8,
           frameon=True, bbox_to_anchor=(0.5, -0.02))

fig.suptitle("Ni(111) 3×3 + 1 Ru dopant — Adsorption sites (top view)", fontsize=14, fontweight="bold")
fig.tight_layout(rect=[0, 0.10, 1, 0.95])
fig.savefig("README_files/site_comparison.png", dpi=150, bbox_inches="tight")
plt.show()
print("Saved to README_files/site_comparison.png")
