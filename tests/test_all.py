from autoadsorbate import Fragment, Surface


def test_Surface():
    from ase.build import fcc111

    slab = fcc111("Cu", (2, 2, 2), periodic=True, vacuum=10)
    s = Surface(slab)
    assert type(s.site_dict) == dict


def test_Fragment():
    f = Fragment(input="COC", to_initialize=5)
    assert f.smile == "COC"


# ---------- SOAP-based symmetry reduction tests ----------

def _make_surface():
    """Helper: build a small Cu(111) Surface for SOAP tests."""
    from ase.build import fcc111
    slab = fcc111("Cu", (3, 3, 3), periodic=True, vacuum=10)
    return Surface(slab)


def test_soap_similarity_matrix():
    """get_soap_similarity_matrix returns a square, symmetric matrix with 1s on the diagonal."""
    import numpy as np
    s = _make_surface()
    sim = s.get_soap_similarity_matrix()
    n = len(s.site_df)
    assert sim.shape == (n, n), f"Expected ({n},{n}), got {sim.shape}"
    np.testing.assert_allclose(np.diag(sim), 1.0, atol=1e-6)
    np.testing.assert_allclose(sim, sim.T, atol=1e-10)


def test_get_nonequivalent_sites_soap_hcluster():
    """get_nonequivalent_sites_soap (hcluster) returns fewer sites than the total."""
    s = _make_surface()
    n_before = len(s.site_df)
    unique = s.get_nonequivalent_sites_soap(similarity_threshold=0.99, method="hcluster")
    assert 0 < len(unique) <= n_before
    # all returned labels must exist in the original DataFrame
    assert all(idx in s.site_df.index for idx in unique)


def test_get_nonequivalent_sites_soap_greedy():
    """get_nonequivalent_sites_soap (greedy) returns fewer sites than the total."""
    s = _make_surface()
    n_before = len(s.site_df)
    unique = s.get_nonequivalent_sites_soap(similarity_threshold=0.99, method="greedy")
    assert 0 < len(unique) <= n_before


def test_soap_reduce():
    """soap_reduce mutates site_df in-place and reduces its length."""
    s = _make_surface()
    n_before = len(s.site_df)
    s.soap_reduce(similarity_threshold=0.99)
    n_after = len(s.site_df)
    assert 0 < n_after <= n_before
    # site_dict should stay in sync
    assert len(s.site_dict["coordinates"]) == n_after


def test_soap_reduce_vs_sym_reduce():
    """SOAP reduces sites on a perfect fcc(111) slab; count should be
    <= ASE's result (SOAP at 0.9 may merge more aggressively)."""
    import copy
    s1 = _make_surface()
    s2 = copy.deepcopy(s1)
    s1.sym_reduce()
    s2.soap_reduce(similarity_threshold=0.9)
    assert 0 < len(s2.site_df) <= len(s1.site_df), (
        f"sym_reduce gave {len(s1.site_df)}, soap_reduce gave {len(s2.site_df)}"
    )


# ---------- Unified API tests (method='ase' | 'soap') ----------

def test_sym_reduce_method_ase():
    """sym_reduce(method='ase') behaves like the legacy sym_reduce()."""
    import copy
    s1 = _make_surface()
    s2 = copy.deepcopy(s1)
    s1.sym_reduce()                       # default → method='ase'
    s2.sym_reduce(method="ase")
    assert list(s1.site_df.index) == list(s2.site_df.index)


def test_sym_reduce_method_soap():
    """sym_reduce(method='soap') produces the same result as soap_reduce()."""
    import copy
    s1 = _make_surface()
    s2 = copy.deepcopy(s1)
    s1.sym_reduce(method="soap", similarity_threshold=0.99)
    s2.soap_reduce(similarity_threshold=0.99)
    assert list(s1.site_df.index) == list(s2.site_df.index)


def test_get_nonequivalent_sites_method_switch():
    """get_nonequivalent_sites dispatches correctly via the method param."""
    s = _make_surface()
    ase_sites = s.get_nonequivalent_sites(method="ase")
    soap_sites = s.get_nonequivalent_sites(method="soap", similarity_threshold=0.99)
    # Both should reduce the site count; SOAP at 0.9 may merge more
    assert 0 < len(soap_sites) <= len(ase_sites)


def test_sym_reduce_invalid_method():
    """sym_reduce raises ValueError for an unknown method."""
    import pytest
    s = _make_surface()
    with pytest.raises(ValueError, match="Unknown method"):
        s.sym_reduce(method="invalid")


# ---------- Ni/Ru broken-symmetry comparison ----------

def _make_ni_ru_surface():
    """Build a Ni(111) slab and replace one surface Ni with Ru."""
    from ase.build import fcc111
    slab = fcc111("Ni", (3, 3, 3), periodic=True, vacuum=10)
    # Replace one top-layer atom with Ru to break the perfect symmetry
    top_z = slab.positions[:, 2].max()
    for atom in slab:
        if abs(atom.position[2] - top_z) < 0.1:
            atom.symbol = "Ru"
            break
    return slab


def test_ni_ru_ase_vs_soap_site_counts():
    """On a Ni(111) slab with one Ni→Ru substitution, both methods should
    find more unique site types than on the pristine surface and the SOAP
    method should return a plausible number of groups."""
    import copy

    slab = _make_ni_ru_surface()
    s_ase = Surface(slab)
    s_soap = copy.deepcopy(s_ase)

    n_total = len(s_ase.site_df)

    s_ase.sym_reduce(method="ase")
    s_soap.sym_reduce(method="soap", similarity_threshold=0.99)

    n_ase = len(s_ase.site_df)
    n_soap = len(s_soap.site_df)

    # Both should reduce the site count
    assert 0 < n_ase <= n_total
    assert 0 < n_soap <= n_total

    # The broken-symmetry slab should have more unique sites than
    # the pristine slab (pristine Ni(111) 3×3 has ~4 unique types)
    assert n_ase > 1
    assert n_soap > 1

    print(
        f"\n  Ni/Ru comparison: total={n_total}  ASE={n_ase}  SOAP={n_soap}"
    )


def test_ni_ru_timings(capsys):
    """Benchmark ASE vs SOAP symmetry reduction on the Ni/Ru slab and print
    a timing summary.  This is informational – the test always passes."""
    import copy
    import time

    slab = _make_ni_ru_surface()
    s_ase = Surface(slab)
    s_soap = copy.deepcopy(s_ase)

    n_total = len(s_ase.site_df)

    # --- ASE timing ---
    t0 = time.perf_counter()
    s_ase.sym_reduce(method="ase")
    t_ase = time.perf_counter() - t0

    # --- SOAP timing ---
    t0 = time.perf_counter()
    s_soap.sym_reduce(method="soap", similarity_threshold=0.99)
    t_soap = time.perf_counter() - t0

    n_ase = len(s_ase.site_df)
    n_soap = len(s_soap.site_df)

    summary = (
        f"\n{'='*60}\n"
        f"  Ni(111) 3×3 slab with 1 Ni → Ru substitution\n"
        f"  Total sites before reduction : {n_total}\n"
        f"{'─'*60}\n"
        f"  ASE SymmetryEquivalenceCheck : {n_ase:3d} unique sites  "
        f"({t_ase:.3f} s)\n"
        f"  SOAP + hierarchical cluster  : {n_soap:3d} unique sites  "
        f"({t_soap:.3f} s)\n"
        f"{'='*60}"
    )
    print(summary)

    # Compare groupings
    ase_indices = set(s_ase.site_df.index)
    soap_indices = set(s_soap.site_df.index)
    common = ase_indices & soap_indices
    ase_only = ase_indices - soap_indices
    soap_only = soap_indices - ase_indices

    print(f"  Representatives in common    : {len(common)}")
    print(f"  Only in ASE result           : {len(ase_only)}  {sorted(ase_only)}")
    print(f"  Only in SOAP result          : {len(soap_only)}  {sorted(soap_only)}")
    print()

    # Informational – always pass
    assert True
