---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.2'
      jupytext_version: 1.4.2
  kernelspec:
    display_name: Julia 1.4.1
    language: julia
    name: julia-1.4
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
using Pkg
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
1. Generate a SVG format image by `drawsvg`. Image size should be specified (width=300, height=300).
1. `display` the molecule on the notebook.

As SMILES does not have coordinates of atoms, so the 2D coordinates will be generated when `smilestomol` is called. Internally `MolecularGraph` uses Schrodinger's [coordgenlibs](https://github.com/schrodinger/coordgenlibs) for 2D coords generation.

```julia
mol2 = smilestomol("O=C3N2/C(=C(/C=C\\c1scnc1C)CS[C@@H]2[C@@H]3NC(=O)C(=N\\OC)/c4nc(sc4)N)C(=O)O")
molsvg2 = drawsvg(mol2, 300, 300)
display("image/svg+xml",  molsvg2)
```

## Dealing with hydrogens

SDFiles downloaded from PubChem have hydrogen nodes. In practice, hydrogens which is not important are removed from molecular graphs for simplicity.

`removehydrogens(mol, all=false)` removes hyhdrogen nodes that are not important  (no charge, no unpaired electron, no isotope information and no stereochemistry).  If the `all` option is set to `true`,  all hydrogen nodes will be removed.

```julia
mol = graphmol(removehydrogens(mol, all=false))
molsvg = drawsvg(mol, 300, 300)
display("image/svg+xml",  molsvg)
```

Even if these hydrogen "nodes" are removed from the graph,  we can infer the actual number of hydrogens from the chemical structure (so called "implicit hydrogens"). `implicithconnected` returns a vector of the number of implicit hydrogens attached to each atom nodes (how to know each atom indices are described later).

```julia
println(implicithconnected(mol))
```

## Basic descriptors

Some other practically useful molecular descriptor methods are listed below.

- `atomsymbol(mol)`: returns atom symbols of each atoms in `Symbol` type
- `charge(mol)`: returns atom charges of each atoms
- `hybridization(mol)`: returns orbital hybridization modes of each atoms
- `isrotatable(mol)`: returns if each bonds are rotatable or not
- `isaromatic(mol)`: returns if each atoms belong to aromatic rings or not
- `isaromaticbond(mol)`: returns if each bonds belong to aromatic rings or not
- `sssr(mol)`: returns smallest set of smallest rings(SSSR)

Running `precalculate!(mol)` before these functions is recommended. This will cache some calculation results of costful and frequently used functions for fast computation.

```julia
precalculate!(mol)
println("atomsymbol: \n", atomsymbol(mol), "\n")
println("hybridization: \n", hybridization(mol), "\n")
println("sssr: \n", sssr(mol), "\n")
```

## Display atom indices

If the molecule is parsed from SDFile, node indices and edge indices are same as the order of SDFile atom/bond records. 

Similary, if the molecule is parsed from SMILES, node indices are same as the order of characters found in SMILES string. Implicit single bonds are also considered as there are `-` characters. Positions of ring bonds are where the ring is closed.

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

- `standardweight(Float64, mol)`: returns standard molecular weight.
- `hacceptorcount(mol)`: returns number of Hydrogen bond acceptors (F, O or N).
- `hdonorcount(mol)`: returns number of Hydrogen bond donors (O or N attached to at least one hydrogen).
- `wclogp(mol)` returns predicted LogP value by Wildman and Crippen method. 
- `rotatablecount(mol)` returns number of rotatable bonds.

```julia
println("Molecular weight: ", standardweight(Float64, mol))
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
