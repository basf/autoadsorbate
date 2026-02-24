"""Detailed comparison of ASE vs SOAP site reduction on Ni/Ru slab."""
import copy
import numpy as np
from ase.build import fcc111
from autoadsorbate import Surface
from autoadsorbate.utils import _compute_soap_site_vectors, _soap_similarity_matrix
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

# Build Ni/Ru slab
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

# Compute SOAP cluster labels for all sites
df = s.site_df.copy()
df["in_ASE"] = [i in ase_idx for i in df.index]
df["in_SOAP"] = [i in soap_idx for i in df.index]
df["formula_str"] = [str(f) for f in df["site_formula"]]

vecs = _compute_soap_site_vectors(s.atoms, s.site_df)
sim = _soap_similarity_matrix(vecs)
dist = np.clip(1.0 - sim, 0.0, 2.0)
np.fill_diagonal(dist, 0.0)
Z = linkage(squareform(dist, checks=False), method="average")
labels = fcluster(Z, t=0.01, criterion="distance")
df["soap_cluster"] = labels

print("=" * 90)
print("ALL SITES — grouped by SOAP cluster")
print("=" * 90)
for cl in sorted(df["soap_cluster"].unique()):
    members = df[df["soap_cluster"] == cl]
    rep_ase = [i for i in members.index if i in ase_idx]
    rep_soap = [i for i in members.index if i in soap_idx]
    print(f"\nSOAP Cluster {cl}  ({len(members)} sites)")
    print(f"  ASE representatives : {rep_ase}")
    print(f"  SOAP representative : {rep_soap}")
    print(f"  Members:")
    for idx, row in members.iterrows():
        conn = row["connectivity"]
        formula = str(row["site_formula"])
        marker = ""
        if idx in soap_idx:
            marker = " <-- SOAP rep"
        elif idx in ase_idx:
            marker = " <-- ASE only"
        print(
            f"    idx={idx:3d}  conn={conn:2d}  formula={formula:20s}"
            f"  ASE={row['in_ASE']}  SOAP={row['in_SOAP']}{marker}"
        )

print()
print("=" * 90)
print("SUMMARY BY SITE TYPE")
print("=" * 90)

type_groups = df.groupby(["connectivity", "formula_str"])
header = f"{'Type':30s} {'Total':>6s} {'ASE reps':>8s} {'SOAP reps':>9s} {'SOAP clusters':>14s}"
print(header)
print("-" * 70)
for (conn, formula), grp in type_groups:
    n_total = len(grp)
    n_ase = int(grp["in_ASE"].sum())
    n_soap = int(grp["in_SOAP"].sum())
    clusters = sorted(grp["soap_cluster"].unique())
    print(f"conn={conn} {formula:22s} {n_total:6d} {n_ase:8d} {n_soap:9d}   {clusters}")

print()
print(f"Total sites: {len(df)}")
print(f"ASE unique:  {len(ase_idx)}")
print(f"SOAP unique: {len(soap_idx)}")
print(f"SOAP clusters: {len(df['soap_cluster'].unique())}")
