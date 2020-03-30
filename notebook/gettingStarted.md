---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.2'
      jupytext_version: 1.4.0
  kernelspec:
    display_name: Julia 1.3.1
    language: julia
    name: julia-1.3
---

# Getting started

This tutorial includes following fundamental operations of molecular mining.

- Fetch test molecule data from public resources (PubChem)
- Draw chemical structure
- Preprocess molecule data
- Calculate atom/bond descriptors
- Calculate molecular properties
- Find functional groups


## Loading MolecularGraph package

If you cloned the notebook tutorial package according to [Quickstart](https://github.com/mojaie/MolecularGraph.jl_notebook), this notebook file (.ipynb) is located in `notebook` under the package root directory. So `Pkg.activate("..")` to activate the package,  and then load `MolecularGraph`. 

```julia
import Pkg
Pkg.activate("..")
using MolecularGraph
```

## Fetching test molecule data from public resources

Chemical structure data for tutorials can be downloaded from PubChem via HTTP. `_data` is created as a temporary data folder, and all test data in this tutorial will be stored in it.

```julia
cid = "6437877"  # PubChem CID
name = "Cefditoren Pivoxil"  # FIle name

# Create data directory
data_dir = "_data"
isdir(data_dir) || mkdir(data_dir)

# Fetch
url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/$(cid)/SDF"
dest = joinpath(data_dir, "$(name).mol")
isfile(dest) || download(url, dest);
```

## Chemical structure drawing

From SDFIle

1. Create molecular object from SDFile (.mol or .sdf) by `sdftomol`
1. Generate a SVG format image by `drawsvg`. Image size should be specified (width=300, height=300).
1. `display` the molecule on the notebook.

```julia
mol = sdftomol(joinpath(data_dir, "Cefditoren Pivoxil.mol"))
molsvg = drawsvg(mol, 300, 300)
display("image/svg+xml",  molsvg)
```

From SMILES

1. Create molecular object from SMILES string by `smilestomol`.
1. Apply `kekulize!` and `setdiastereo!`.
1. Generate a SVG format image by `drawsvg`. Image size should be specified (width=300, height=300).
1. `display` the molecule on the notebook.

As SMILES does not have coordinates of atoms, so the 2D coordinates will be generated when SMILES are passed to `drawsvg`. Internally `MolecularGraph` uses Schrodinger's [coordgenlibs](https://github.com/schrodinger/coordgenlibs) for 2D coords generation.

```julia
mol2 = smilestomol("O=C3N2/C(=C(/C=C\\c1scnc1C)CS[C@@H]2[C@@H]3NC(=O)C(=N\\OC)/c4nc(sc4)N)C(=O)O")
kekulize!(mol2)
setdiastereo!(mol2)
molsvg2 = drawsvg(mol2, 300, 300)
display("image/svg+xml",  molsvg2)
```

## Dealing with hydrogens

SDFiles downloaded from PubChem have hydrogen nodes. In practice, hydrogens which is not important are removed from molecular graphs for simplicity. `removehydrogens`will work for this purpose.

- If you do not want to remove hydrogens involved in stereochemistry, run `setstereocenter!` to calculate stereocenter flags. If hydrogens are removed by `removehydrogens` after that, the stereocenter atoms still have stereochemistry information in `Atom.stereo`.

- `removehydrogens(mol, all=true)`(default) will return a molecule with all hydrogen nodes removed, whereas `all=false` will return a molecule with only trivial hydrogens removed (no charge, no unpaired electron, no isotope information and no stereochemistry).

`addhydrogens` will return a molecule with all hydrogen nodes are explicitly attached. Note that newly attached hydrogens have no coordinate so no longer available for drawing (`coords2d(mol, recalculate=true)` worth trying).

```julia
setstereocenter!(mol)
mol = graphmol(removehydrogens(mol, all=false))
molsvg = drawsvg(mol, 300, 300)
display("image/svg+xml",  molsvg)
```

If hydrogens are removed from the molecular graph, we can easily know the number of hydrogens attached to non-hydrogen atoms by `hcount`. `hcount` returns a vector of hydrogen count of each atom nodes in atom index order.

```julia
println(hcount(mol))
```

## Basic descriptors

Some other practically useful molecular descriptor methods are listed below.

- `atomsymbol(mol)`: returns a vector of atom symbols
- `charge(mol)`: returns a vector of atom charges
- `nodedegree(mol)`: returns a vector of molecular graph degree of each atoms (ignore implicit Hs)
- `connectivity(mol)`: returns a vector of neighbor atoms count (consider implicit Hs)
- `heavyatoms(mol)`: returns a vector of neighbor heavy atoms count (ignore all Hs)
- `pielectron(mol)`: returns an vector of pi electron count of each atoms (-> orbital hybridization)
- `isrotatable(mol)`: returns a vector that indicates if each bonds are rotatable or not
- `isaromatic(mol)`: returns a vector that indicates if each atoms belong to aromatic rings or not
- `sssr(mol)`: returns smallest set of smallest rings(SSSR)

Some of these molecular descriptor vector is calculated when it is called for the first time and is stored in `cache` field of the molecular graph object. Cached results will be returned when it is called for the next time. Internally, `@cache` macro is used for this mechanism (see methods in src/properties.jl file). If you altered the molecular graph structure by adding/deleting graph components, call `clearcache!` to initialize the cache field.


```julia
println("atomsymbol: \n", atomsymbol(mol), "\n")
println("pielectron: \n", pielectron(mol), "\n")
println("sssr: \n", sssr(mol), "\n")
```

## Display atom indices

`SDFileAtom` node indices and `SDFileBond` edge indices are same as the order of SDFile atom/bond records. 

Similary, `SmilesAtom` and `SmilesBond` indices are same as the order of characters found in SMILES string. Implicit single bonds are also considered as there are `-` characters. Positions of ring bonds are where the ring is closed.

To know which atom is labeled with which index number, `drawatomindex!` can be used.

Internally, `drawsvg` creates `SvgCanvas` object, calls `draw2d!` to set molecule components and then calls `tosvg` to finalize SVG image. `drawatomindex!` adds atom indices to the molecular graph drawing before finalizing.


```julia
canvas = SvgCanvas()
draw2d!(canvas, mol)
drawatomindex!(canvas, mol)
mol_svg = tosvg(canvas, 400, 400)
display("image/svg+xml",  mol_svg)
```

## Calculate molecular properties

- `molweight(mol)`: returns standard molecular weight.
- `hacceptorcount(mol)`: returns number of Hydrogen bond acceptors (F, O or N).
- `hdonorcount(mol)`: returns number of Hydrogen bond donors (O or N attached to at least one hydrogen).
- `wclogp(mol)` returns predicted LogP value by Wildman and Crippen method. 
- `rotatablecount(mol)` returns number of rotatable bonds.

```julia
println("Molecular weight: ", molweight(mol))
println("Hydrogen acceptors: ", hacceptorcount(mol))
println("Hydrogen donors: ", hdonorcount(mol))
println("LogP: ", wclogp(mol))
println("Rotatable bonds: ", rotatablecount(mol))
```

## Functional group analysis

`MolecularGraph.jl` provides SMARTS and terminology graph-based functional group analysis.

SMARTS queries listed in YAML files in assets/funcgroup folder are applied to find functional groups from the molecule. The SMARTS query list in the files have graph structure like Gene Ontology or ChEMBL CHEMINF. Each SMARTS records may have isa or has fields.

- `isa`: e.g. primary alcohol `isa` alcohol but not all alcohols are primary alcohol.
- `has`: e. g. carboxyl group `has` a hydroxy group and a carbonyl group.

The substructure terminology graph indicates which substructures does the molecule have, whereas GO/CHEMINF indicates which biological function does the gene/molecule have.

```julia
fg = functionalgroupgraph(mol)

for (term, components) in fg.componentmap
    nodes = [sort(collect(comp)) for comp in components]
    println("Group: $(string(term)): ", nodes...)
end
```

Some of these terms are subset of others (for example, thiazole `isa` five membered ring). Largest in size and most specified terms often give us important information. `largest components` enriches functional group information by collecting largest and most specified terms from the functional group tree.

```julia
display("image/svg+xml",  mol_svg)
for (term, components) in largestcomponents(fg)
    nodes = [sort(collect(comp)) for comp in components]
    if !isempty(nodes)
        println("Group: $(string(term)): ", nodes...)
    end
end
```
