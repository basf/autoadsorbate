from autoadsorbate import Fragment, Surface


def test_Surface():
    from ase.build import fcc111

    slab = fcc111("Cu", (2, 2, 2), periodic=True, vacuum=10)
    s = Surface(slab)
    assert type(s.site_dict) == dict


def test_Fragment():
    f = Fragment(smile="COC", to_initialize=5)
    assert f.smile == "COC"
