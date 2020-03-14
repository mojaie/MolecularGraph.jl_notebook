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

- Activate the project and import MolecularGraph.

```julia
import Pkg
Pkg.activate("..")
using MolecularGraph
```

- Chemical structure data for tutorials can be downloaded from PubChem via HTTP.

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

- `sdftomol` converts a SDFile (.mol or .sdf) into a molecular object.
- `drawsvg` generates SVG format image in the given size (width=300, height=300).
- Then, `display` the molecule on the notebook.

```julia
pcmol = sdftomol(joinpath(data_dir, "Cefditoren Pivoxil.mol"))
mol_svg = drawsvg(pcmol, 300, 300)
display("image/svg+xml",  mol_svg)
```

- Oops, there are so many trivial hydrogens.
- Hydrogens can be removed by `makehydrogensimplicit`.
- Technically, `makehydrogensimplicit` returns "subgraph view" of the original molecule. In this case, the original molecule would not be used further, so generate new molecule from the view by using `graphmol`.

```julia
mol = graphmol(makehydrogensimplicit(pcmol))
mol_svg = drawsvg(mol, 300, 300)
display("image/svg+xml",  mol_svg)
```

## Calculate molecular properties

- `molweight` shows standard molecular weight.
- `hacceptorcount` shows number of Hydrogen bond acceptors (F, O or N).
- `hdonorcount` shows number of Hydrogen bond donors (O or N attached to at least one hydrogen).
- `wclogp` shows predicted LogP value by Wildman and Crippen method. 
- `rotatablecount` shows number of rotatable bonds.

```julia
println("Molecular weight: ", molweight(mol))
println("Hydrogen acceptors: ", hacceptorcount(mol))
println("Hydrogen donors: ", hdonorcount(mol))
println("LogP: ", wclogp(mol))
println("Rotatable bonds: ", rotatablecount(mol))
```

 ## Show atom indices

- `drawsvg` is a convenient method that generates `SvgCanvas` object, calls `draw2d!` to set molecule components and then calls `tosvg` to finalize SVG image.
- `drawatomindex!` adds atom index to the molecular image canvas.
-  '!' of `draw2d!` and `drawatomindex!` is a Julia language convention which means the function is destructive. In this case, `draw2d!` modifies `SvgCanvas` by adding molecule drawing settings and components.
- The default atom indices are the same order as the atom record in SDFile.

```julia
canvas = SvgCanvas()
draw2d!(canvas, mol)
drawatomindex!(canvas, mol)
mol_svg = tosvg(canvas, 400, 400)
display("image/svg+xml",  mol_svg)
```

## Functional group detection

- MolecularGraph.jl provides a functional group detection procedure based on terminology graph.
- It is like a Gene Ontology term graph, but whose nodes are functional group terms and edges are relationship.

```julia
fg = functionalgroupgraph(mol)

display("image/svg+xml",  mol_svg)

for (term, components) in fg.componentmap
    nodes = [sort(collect(comp)) for comp in components]
    println("Group: $(string(term)): ", nodes...)
end
```

- Some of these terms are subset of others (for example, thiazole "is a" five membered ring). Largest in size and most specified terms often give us important information.
- `largest components` enriches functional group information by collecting largest and specified terms from the functional group graph.

```julia
display("image/svg+xml",  mol_svg)
for (term, components) in largestcomponents(fg)
    nodes = [sort(collect(comp)) for comp in components]
    if !isempty(nodes)
        println("Group: $(string(term)): ", nodes...)
    end
end
```
