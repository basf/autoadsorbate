# Table of Contents

- [Installation](#installation)
- [AutoAdsorbate](#AutoAdsorbate)
  - [Fragment](#fragment)
    - [Molecules](#molecules)
    - [Reactive species](#reactive-species)
  - [Surface](#surface)
- [Making surrogate SMILES automatically](#making-surrogate-smiles-automatically)
- [Fully automatic - populate Surface with Fragment](#fully-automatic---populate-surface-with-fragment)


## Installation

The package is designed to be as lightweight as possible, to implement seamlessly into existing environments with complex dependecies. If you `git clone <autoadsorbate>` and just `sys.path.insert(0, <path/to/autoadsorbate>)`, most likely it will work.

- Built on only:
  - `ase`
  - `rdkit`
  - Basic Python packages: `pandas`, `numpy`

The package is available on PyPi:

```python
pip install autoadsorbate
```

Installation from source:
 ```python
git clone <autoadsorbate>
cd autoadsorbate
pip install .
```

## AutoAdsorbate

AutoAdsorbate is a lightweight and easy-to-use Python package for generating chemically meaningful configurations of molecules and fragments on surfaces. Built with minimal dependencies and a low barrier to entry, it enables rapid setup of surface-adsorbate systems using the Surrogate-SMILES (*SMILES) representation. Ideal for researchers in catalysis, nanotech, and materials science, AutoAdsorbate streamlines dataset generation for simulations and machine learning workflows.

The challenge of generating initial structures for heterogeneous catalysis has traditionally been addressed through manual effort. This package offers an alternative, automated approach.

To effectively simulate reactive behavior at surfaces, it is crucial to establish clear definitions within our framework. The following definitions are essential for accurately characterizing the structures of interest:

- __Fragment__:
    - <font color='red'>Molecules</font> – species that retain their corresponding geometries __even when isolated from the surface__.
    - <font color='red'>Reactive species</font> – species that adopt their corresponding geometries __only when attached to the surface__.
    
- __Surface__:
    - The surface is defined simply – <font color='red'>every atom of the slab that can be in contact with an intermediate is considered a surface atom</font>. The surface is the collection of such atoms.
    - Every atom of the surface is a "top" site.
    - When two "top" sites are close (close in its literal meaning), they form a "bridge" site.
    - When three "top" sites are close (close in its literal meaning), they form a "3-fold" site.
    - etc.
    
- __Active Site__:
    - A collection of one or more sites that can facilitate a chemical transformation is called an active site.
    - A "top" site can be an active site only for Eley-Rideal transformations.
    - All other transformations require that at least one intermediate binds through at least two sites. All involved sites compose an active site.
    
- __Intermediate__:
    - Intermediates are fragments bound to an active site.

<!-- ### basic imports -->

### Fragment

Molecules and reactive species are both initialized as the Fragment object (based on ase.Atoms). Some examples are given bellow.

#### Molecules

Let us initialize a molecule of dimethyl ether (DME):

```python
from autoadsorbate import Fragment

f = Fragment(smile = 'COC', to_initialize = 5)
```


```python
from autoadsorbate import docs_plot_conformers

conformer_trajectory = f.conformers
fig = docs_plot_conformers(conformer_trajectory)
```


    
![png](README_files/README_10_0.png)
    

Notice that the orientation of the fragment is arbitrary. While we could simply place these structures onto the surface of a material, it would be difficult to evaluate the quality of these initial random configurations. This uncertainty would force us to sample a large number of structures and run dynamic simulations to explore local minima and determine which configurations are the most stable.

However, in the case of DME, we can leverage chemical intuition to simplify the problem. The oxygen atom bridging the two methyl groups has two lone electron pairs. By using a simple trick—replacing one of these lone pairs with a marker atom (such as chlorine, Cl)—we can guide the placement more effectively.


Notice that we had to make two adjustments to the SMILES string. To replace the lone pair with a marker atom, we must "trick" the valence of the oxygen atom and rearrange the SMILES formula so that the marker atom appears first (for easier bookkeeping).
    - ```COC``` original
    - ```CO(Cl)C``` add Cl instead of the O lone pair (this is an invalid SMILES)
    - ```C[O+](Cl)C``` trick to make the valence work
    - ```Cl[O+](C)C``` rearrange so taht the SMILES string starts with the marker first (for easy book keeping)

This can be also done with a function:


```python
from autoadsorbate import get_marked_smiles
marked_smile = get_marked_smiles(['COC'])[0]
marked_smile
```
    'Cl[O+](C)(C)'

These surrogate smilles san now be used to initialize a Fragment object (we can set the number of randoms conformers to be initialized):

```python
f = Fragment(smile = 'Cl[O+](C)(C)', to_initialize = 5)
len(f.conformers)
```
    5

We can visualize these structures:
```python
from autoadsorbate import docs_plot_conformers

conformer_trajectory = f.conformers
fig = docs_plot_conformers(conformer_trajectory)
```


    
![png](README_files/README_16_0.png)
    


Now we can use the marker atom to orient our molecule:


```python
from autoadsorbate import docs_plot_sites

oriented_conformer_trajectory = [f.get_conformer(i) for i, _ in enumerate(f.conformers)]
fig = docs_plot_conformers(oriented_conformer_trajectory)
```


    
![png](README_files/README_18_0.png)
    


We can also easily remove the marker:


```python
clean_conformer_trajectory = [atoms[1:] for atoms in oriented_conformer_trajectory]
fig = docs_plot_conformers(clean_conformer_trajectory)
```


    
![png](README_files/README_20_0.png)
    


#### Reactive species 

Methoxy


```python
f = Fragment(smile = 'ClOC', to_initialize = 5)
oriented_conformer_trajectory = [f.get_conformer(i) for i, _ in enumerate(f.conformers)]
fig = docs_plot_conformers(oriented_conformer_trajectory)
```


    
![png](README_files/README_23_0.png)
    


##### Methyl


```python
f = Fragment(smile = 'ClC', to_initialize = 5)
oriented_conformer_trajectory = [f.get_conformer(i) for i, _ in enumerate(f.conformers)]
fig = docs_plot_conformers(oriented_conformer_trajectory)
```


    
![png](README_files/README_25_0.png)
    


##### Frangments with more than one binding mode (e.g. 1,2-PDO)

bound through single site:


```python
f = Fragment(smile = 'Cl[OH+]CC(O)C', to_initialize = 5)
oriented_conformer_trajectory = [f.get_conformer(i) for i, _ in enumerate(f.conformers)]
fig = docs_plot_conformers(oriented_conformer_trajectory)
```


    
![png](README_files/README_28_0.png)
    


Coordinated withboth hydroxil:


```python
f = Fragment(smile = 'S1S[OH+]CC([OH+]1)C', to_initialize = 5)
oriented_conformer_trajectory = [f.get_conformer(i) for i, _ in enumerate(f.conformers)]
fig = docs_plot_conformers(oriented_conformer_trajectory)
```


    
![png](README_files/README_30_0.png)
    


### Surface

Defining the surface of a slab may seem like a simple task, but different approaches can yield varying results depending on the context. When considering catalytic sites, we can define these as surface regions capable of binding a fragment. By using reasonable steric criteria—essentially asking, "Is there enough space for a molecule to bind to that site?"—we can identify all possible binding sites on the slab's surface. These sites can be classified as top, bridge, or multi-fold, depending on how many atoms surround the site.

As an example: First, we need to define a slab (any ```ase.Atoms``` object). A slab is an arrangement of atoms that represents the boundary between a material and another phase, such as gas, fluid, or another material. We can either read an existing slab, or a new slab:


```python
from ase.build import fcc111
slab = fcc111('Cu', (4,4,4), periodic=True, vacuum=10)
```

Now we can initalize the Surface object which associates the constructed slab (ase.Atoms) with additional information required for placing Fragments.
We can view which atoms are in the surface:


```python
s = Surface(slab)
plot_atoms(s.view_surface(return_atoms=True))
```

    Visualizing surface Cu atoms as Zn

    
![png](README_files/README_35_2.png)
    


We have access to all the sites info as a pandas dataframe:


```python
s.site_df.head()
```

<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>coordinates</th>
      <th>connectivity</th>
      <th>topology</th>
      <th>n_vector</th>
      <th>h_vector</th>
      <th>site_formula</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>[0.0, 0.0, 16.252703415323644]</td>
      <td>1</td>
      <td>[48]</td>
      <td>[-0.004670396521231514, -0.0031449903964026822...</td>
      <td>[1.0, 0.0, 0.0]</td>
      <td>{'Cu': 1}</td>
    </tr>
    <tr>
      <th>1</th>
      <td>[0.6381638700208592, 1.105332246430909, 16.252...</td>
      <td>2</td>
      <td>[48, 52]</td>
      <td>[0.0006776311857337964, -0.010516809475472271,...</td>
      <td>[-0.5000000000000001, -0.8660254037844387, 0.0]</td>
      <td>{'Cu': 2}</td>
    </tr>
    <tr>
      <th>2</th>
      <td>[1.2763277400417168, 5.162938145598479e-16, 16...</td>
      <td>2</td>
      <td>[48, 49]</td>
      <td>[-0.011576660085263627, -0.017987208564805915,...</td>
      <td>[-1.0, 0.0, 0.0]</td>
      <td>{'Cu': 2}</td>
    </tr>
    <tr>
      <th>3</th>
      <td>[1.2763277400417183, 0.7368881642872727, 16.25...</td>
      <td>3</td>
      <td>[48, 49, 52]</td>
      <td>[-0.01272989568588465, 0.0042077202541598024, ...</td>
      <td>[-0.5000000000000001, -0.8660254037844387, 0.0]</td>
      <td>{'Cu': 3}</td>
    </tr>
    <tr>
      <th>4</th>
      <td>[1.2763277400417183, 2.210664492861818, 16.252...</td>
      <td>1</td>
      <td>[52]</td>
      <td>[0.0013334161774154326, -0.007734740595549886,...</td>
      <td>[1.0, 0.0, 0.0]</td>
      <td>{'Cu': 1}</td>
    </tr>
  </tbody>
</table>
</div>



or in dict form:


```python
s.site_dict.keys()
```




    dict_keys(['coordinates', 'connectivity', 'topology', 'n_vector', 'h_vector', 'site_formula'])



One can easily get access to sites as ```ase.Atoms``` as well, and find useful information in the ```ase.Atoms.info```:


```python
site_atoms = s.view_site(0, return_atoms=True)
site_atoms.info
```




    {'coordinates': array([ 0.        ,  0.        , 16.25270342]),
     'connectivity': 1,
     'topology': [48],
     'n_vector': array([-0.0046704 , -0.00314499,  0.99998415]),
     'h_vector': array([1., 0., 0.]),
     'site_formula': {'Cu': 1}}


We can visualize a few surface sites:

```python
from autoadsorbate import docs_plot_sites
fig = docs_plot_sites(s)
```


    
![png](README_files/README_42_0.png)
    


We can reduce the complete list of sites based on symmetry (```ase.utils.structure_comparator.SymmetryEquivalenceCheck```):


```python
s.sym_reduce()
s.site_df
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>coordinates</th>
      <th>connectivity</th>
      <th>topology</th>
      <th>n_vector</th>
      <th>h_vector</th>
      <th>site_formula</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>[0.0, 0.0, 16.252703415323644]</td>
      <td>1</td>
      <td>[48]</td>
      <td>[-0.004670396521231514, -0.0031449903964026822...</td>
      <td>[1.0, 0.0, 0.0]</td>
      <td>{'Cu': 1}</td>
    </tr>
    <tr>
      <th>1</th>
      <td>[0.6381638700208592, 1.105332246430909, 16.252...</td>
      <td>2</td>
      <td>[48, 52]</td>
      <td>[0.0006776311857337964, -0.010516809475472271,...</td>
      <td>[-0.5000000000000001, -0.8660254037844387, 0.0]</td>
      <td>{'Cu': 2}</td>
    </tr>
    <tr>
      <th>3</th>
      <td>[1.2763277400417183, 0.7368881642872727, 16.25...</td>
      <td>3</td>
      <td>[48, 49, 52]</td>
      <td>[-0.01272989568588465, 0.0042077202541598024, ...</td>
      <td>[-0.5000000000000001, -0.8660254037844387, 0.0]</td>
      <td>{'Cu': 3}</td>
    </tr>
    <tr>
      <th>8</th>
      <td>[2.552655480083436, 1.4737763285745453, 16.252...</td>
      <td>3</td>
      <td>[49, 52, 53]</td>
      <td>[-0.0011596349368944389, -0.001445905668587753...</td>
      <td>[0.5000000000000002, -0.8660254037844385, 0.0]</td>
      <td>{'Cu': 3}</td>
    </tr>
  </tbody>
</table>
</div>


We can again visualize the sites:

```python
plot_atoms(s.view_surface(return_atoms=True))
```

    Visualizing surface Cu atoms as Zn





![png](README_files/README_45_2.png)
    


## Making surogate SMILES automatically

Simple methods of brute force SMILES enumeration are implemented as well. For example, only using a few lines of code we can initialize multiple conformers of all reaction intermediaries in the nitrogen hydrogenation reaction. A template of the required information can be found here:

```python
from autoadsorbate import _example_config
_example_config
```




    {'backbone_info': {'C': 1, 'N': 0, 'O': 2},
     'allow_intramolec_rings': True,
     'ring_marker': 2,
     'side_chain': ['(', ')'],
     'brackets': ['[', ']', 'H2]', 'H3]', 'H-]', 'H+]'],
     'make_labeled': True}


Now we can use (or edit) this information as we see fit:

```python
from autoadsorbate import construct_smiles
 
config = {
'backbone_info': {'C': 0, 'O': 0, 'N':2},
'allow_intramolec_rings': True,
'ring_marker': 2,
'side_chain': ['(', ')'],
'brackets': ['[', ']', 'H+]', 'H2+]', 'H3+]'],
'make_labeled': True
}

smiles = construct_smiles(config)
```
We now have a list of surrgate SMILES that can be used to initalize Fragment objects.


```python
smiles
```




    ['ClNN',
     'Cl[N]N',
     'Cl[NH+]N',
     'Cl[NH2+]N',
     'ClN[N]',
     'ClN[NH+]',
     'ClN[NH2+]',
     'ClN[NH3+]',
     'Cl[N][N]',
     'Cl[N][NH+]',
     'Cl[N][NH2+]',
     'Cl[N][NH3+]',
     'Cl[NH+][N]',
     'Cl[NH+][NH+]',
     'Cl[NH+][NH2+]',
     'Cl[NH+][NH3+]',
     'Cl[NH2+][N]',
     'Cl[NH2+][NH+]',
     'Cl[NH2+][NH2+]',
     'Cl[NH2+][NH3+]',
     'S1SN1N',
     'S1SNN1',
     'S1S[N]N1',
     'S1S[NH+]1N',
     'S1S[NH+]N1',
     'S1S[NH2+]N1',
     'S1SN1[N]',
     'S1SN1[NH+]',
     'S1SN1[NH2+]',
     'S1SN1[NH3+]',
     'S1S[N][N]1',
     'S1S[N][NH+]1',
     'S1S[N][NH2+]1',
     'S1S[NH+]1[N]',
     'S1S[NH+]1[NH+]',
     'S1S[NH+][NH+]1',
     'S1S[NH+]1[NH2+]',
     'S1S[NH+][NH2+]1',
     'S1S[NH+]1[NH3+]',
     'S1S[NH2+][NH2+]1',
     'ClN=N',
     'Cl[NH+]=N',
     'ClN=[N]',
     'ClN=[NH+]',
     'ClN=[NH2+]',
     'Cl[NH+]=[N]',
     'Cl[NH+]=[NH+]',
     'Cl[NH+]=[NH2+]',
     'S1SN=N1',
     'S1S[NH+]=N1',
     'S1S[NH+]=[NH+]1',
     'S1SN1#N']




```python
from autoadsorbate import Fragment
 
trj = []
for s in smiles:
    try:
        f = Fragment(s, to_initialize=1)
        a = f.get_conformer(0)
        trj.append(a)
    except:
        pass
 
lst = [z for z in zip([a.get_chemical_formula() for a in trj],trj)]
lst.sort(key=lambda tup: tup[0])
trj =  [a[1] for a in lst]
len(trj)
```




    52

From the list of initialized conformers we can remove the ones that are effectively identical:


```python
from autoadsorbate import get_drop_snapped
 
xtrj = get_drop_snapped(trj, d_cut=1.5)
len(xtrj)
```

    33

We can visualize these structures:


```python
import matplotlib.pyplot as plt
from ase.visualize.plot import plot_atoms
from ase import Atoms
 
fig, axs = plt.subplots(3,11, figsize=[10,5], dpi=100)
 
for i, ax in enumerate(axs.flatten()):
    try:
        platoms = xtrj[i].copy()
         
    except:
        platoms = Atoms('X', positions = [[0,0,0]])
 
    for atom in platoms:
        if atom.symbol in ['Cl', 'S']:
            atom.symbol = 'Ga'
    plot_atoms(platoms, rotation=('-90x,0y,0z'), ax=ax)
    ax.set_axis_off()
    ax.set_xlim(-1, 5)
    ax.set_ylim(-0.5, 5.5)
 
fig.set_layout_engine(layout='tight')
```


    
![png](README_files/README_52_0.png)
    


## Fully automatic - populate Surface with Fragment

A autonomous mode of Fragment placement on Surface is also implemented. The method tries to minimze the overlap of the Fragment and Surface while keeping the requested connectivity to the surface.

```python
from ase.build import fcc211
from autoadsorbate import Surface, Fragment

slab = fcc211(symbol = 'Cu', size=(6,3,3), vacuum=10)  # any ase.Atoms object
s=Surface(slab, touch_sphere_size=2.7)                 # finding all surface atoms
s.sym_reduce()                                         # keeping only non-identical sites

fragments = [
    Fragment('S1S[OH+]CC(N)[OH+]1', to_initialize=20), # For each *SMILES we can request a differnet number of conformers 
    Fragment('Cl[OH+]CC(=O)[OH+]', to_initialize=5)    # based on how much conformational complexity we expect.
]

out_trj = []
for  fragment in fragments:
    out_trj += s.get_populated_sites(
      fragment,                    # Fragment object
      site_index='all',            # a single site can be provided here
      sample_rotation=True,        # rotate the Fragment around the surface-fragment bond?
      mode='heuristic',            # 'all' or 'heuristic', if heuristic surrogate smiles with 'Cl...' will be matched with top sites, etc. 
      conformers_per_site_cap=5,   # max number of conformers to sample
      overlap_thr=1.6,             # tolerated bond overlap betwen the surface and fragment      
      verbose=True
      )
    print('out_trj ', len(out_trj))
```

    conformers 40
    sites 9
    SUCCESS! Found the requested numer of conformers with condition: ovelap_thr = 1.6. Found 5 / 5.
    WARNING: Failed to find requested number of conformers with condition: ovelap_thr = 1.6. Found 0 / 5. Consider setting a higher Fragment(to_initialize = < N >)
    WARNING: Failed to find requested number of conformers with condition: ovelap_thr = 1.6. Found 1 / 5. Consider setting a higher Fragment(to_initialize = < N >)
    SUCCESS! Found the requested numer of conformers with condition: ovelap_thr = 1.6. Found 5 / 5.
    SUCCESS! Found the requested numer of conformers with condition: ovelap_thr = 1.6. Found 5 / 5.
    SUCCESS! Found the requested numer of conformers with condition: ovelap_thr = 1.6. Found 5 / 5.
    SUCCESS! Found the requested numer of conformers with condition: ovelap_thr = 1.6. Found 5 / 5.
    WARNING: Failed to find requested number of conformers with condition: ovelap_thr = 1.6. Found 3 / 5. Consider setting a higher Fragment(to_initialize = < N >)
    WARNING: Failed to find requested number of conformers with condition: ovelap_thr = 1.6. Found 0 / 5. Consider setting a higher Fragment(to_initialize = < N >)
    out_trj  29
    conformers 40
    sites 3
    SUCCESS! Found the requested numer of conformers with condition: ovelap_thr = 1.6. Found 5 / 5.
    SUCCESS! Found the requested numer of conformers with condition: ovelap_thr = 1.6. Found 5 / 5.
    SUCCESS! Found the requested numer of conformers with condition: ovelap_thr = 1.6. Found 5 / 5.
    out_trj  44



```python

```
