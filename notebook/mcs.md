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

# Maximum common substructure (MCS)

MolecularGraph.js implements essential MCS methods for cheminformatics applications.

- Node induced (MCIS) / Edge induced (MCES)
- Connected / disconnected
- Topological constraint (tdMCS)
- Graph-based local similarity (GLS)


```julia
import Pkg
Pkg.activate("..")
using MolecularGraph
using MolecularGraph.MolecularGraphModel

# Download test data from PubChem

# PubChem ID->Name
compounds = Dict(
    "6047" => "Levodopa",
    "74217" => "3-Aminocoumarin",
    "6437877" => "Cefditoren Pivoxil",
    "5481173" => "Ceftazidime"
)

# Create data directory
data_dir = "_data"
isdir(data_dir) || mkdir(data_dir)

# Fetch
for (cid, name) in compounds
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/$(cid)/SDF"
    dest = joinpath(data_dir, "$(name).mol")
    isfile(dest) || download(url, dest);
end

# Convenient function to display a pair of mol images
function displayimgpair(img1, img2)
    """<div>
        <div style="float:left">$(img1)</div>
        <div style="float:left">$(img2)</div>
    </div>"""
end
```

- Let's try MCS functionalities with very small molecules.
- Drop hydrogens from downloaded molecule data by `makehydrogensimplicit`

```julia
ldopa = sdftomol(joinpath(data_dir, "Levodopa.mol"))
ldopa = graphmol(makehydrogensimplicit(ldopa))
svg1 = drawsvg(ldopa, 200, 200)

aminocoumarin = sdftomol(joinpath(data_dir, "3-Aminocoumarin.mol"))
aminocoumarin = graphmol(makehydrogensimplicit(aminocoumarin))
svg2 = drawsvg(aminocoumarin, 200, 200)

display("text/html", displayimgpair(svg1, svg2))
```

## Maximum common induced substructure (MCIS)

- `mcismol` method calculates maximum common induced substructure (MCIS).
- `mcismol` returns a tuple `(mapping, status)`.
- `mapping` (type `Dict{Int,Int}`) is a isomorphisim mapping of node(atom) indices between two molecules.
- `nodesubgraph` can be used to retrieve the matched subgraph and then display the substructure.
- `status`(type `Symbol`) shows the calculation status. `:done` means that the calculation was successfully finished.

```julia
(mapping, status) = mcismol(ldopa, aminocoumarin)

subg1 = nodesubgraph(ldopa, Set(keys(mapping)))
svg1 = drawsvg(subg1, 200, 200)

subg2 = nodesubgraph(aminocoumarin, Set(values(mapping)))
svg2 = drawsvg(subg2, 200, 200)

display("text/html", displayimgpair(svg1, svg2))
println("Status: $(status)")
println("MCIS size: $(length(mapping))")
```

- Use `drawsvg` with `highlight` keyword to display the matched substructure.

```julia
svg1 = drawsvg(ldopa, 200, 200, highlight=subg1)
svg2 = drawsvg(aminocoumarin, 200, 200, highlight=subg2)

display("text/html", displayimgpair(svg1, svg2))
println("Status: $(status)")
println("MCIS size: $(length(mapping))")
```

- Default node matcher function `nodematch` compares two atoms by their atomic symbol (`:C`, `:O`, `:N`, ...) and the number of $\pi$ electrons.
- The number of $\pi$ electrons gives as more intrinsic information of 3D geometry of the molecule than aparent bond order that cannot deal with conjugated system of the molecular orbital.
- Other atomic properties (charges, multiplicity, isotope, ...) of atoms are ignored by default.
- You can customize `nodematch` function and pass it to `mcismol` as an optional argument `nodematcher`.


## Maximum common edge-induced substructure (MCES)

- `mcesmol` calculates maximum common edge-induced substructure (MCES).
- Most chemists may feel that MCES is more intuitive than MCIS.
- `mcesmol` also returns a tuple `(mapping, status)`.
- Note that the result is an mapping of edge indices between two molecules. So use `edgesubgraph` to show the matched subgraph.

```julia
(mapping, status) = mcesmol(ldopa, aminocoumarin)

subg1 = edgesubgraph(ldopa, Set(keys(mapping)))
svg1 = drawsvg(ldopa, 200, 200, highlight=subg1)

subg2 = edgesubgraph(aminocoumarin, Set(values(mapping)))
svg2 = drawsvg(aminocoumarin, 200, 200, highlight=subg2)

display("text/html", displayimgpair(svg1, svg2))
println("Status: $(status)")
println("MCES size: $(length(mapping))")
```

## Connected MCS

- By default, `mcesmol` and `mcismol` calculate disconnected MCS.
- Pass keyword argument `connected=True` to them to enable connected MCS mode.
- Connected MCS is far faster than disconnected MCS but tend to be smaller in size.
- As disconnected MCS does not care the spatial relationship among matched fragments, connected MCS can be more appropriate option for some kind of applications.

```julia
(mapping, status) = mcismol(ldopa, aminocoumarin, connected=true)

subg1 = nodesubgraph(ldopa, Set(keys(mapping)))
svg1 = drawsvg(ldopa, 200, 200, highlight=subg1)

subg2 = nodesubgraph(aminocoumarin, Set(values(mapping)))
svg2 = drawsvg(aminocoumarin, 200, 200, highlight=subg2)

display("text/html", displayimgpair(svg1, svg2))
println("Status: $(status)")
println("Connected MCIS size: $(length(mapping))")
```

```julia
(mapping, status) = mcesmol(ldopa, aminocoumarin, connected=true)

subg1 = edgesubgraph(ldopa, Set(keys(mapping)))
svg1 = drawsvg(ldopa, 200, 200, highlight=subg1)

subg2 = edgesubgraph(aminocoumarin, Set(values(mapping)))
svg2 = drawsvg(aminocoumarin, 200, 200, highlight=subg2)

display("text/html", displayimgpair(svg1, svg2))
println("Status: $(status)")
println("Connected MCES size: $(length(mapping))")
```

## Working with larger molecules

- Then, let's try it with a little larger molecules.

```julia
cefditoren = sdftomol(joinpath(data_dir, "Cefditoren Pivoxil.mol"))
cefditoren = graphmol(makehydrogensimplicit(cefditoren))
svg1 = drawsvg(cefditoren, 250, 250)

ceftazidime = sdftomol(joinpath(data_dir, "Ceftazidime.mol"))
ceftazidime = graphmol(makehydrogensimplicit(ceftazidime))
svg2 = drawsvg(ceftazidime, 250, 250)

display("text/html", displayimgpair(svg1, svg2))
```

```julia
(mapping, status) = mcesmol(cefditoren, ceftazidime, timeout=10)

subg1 = edgesubgraph(cefditoren, Set(keys(mapping)))
svg1 = drawsvg(cefditoren, 250, 250, highlight=subg1)

subg2 = edgesubgraph(ceftazidime, Set(values(mapping)))
svg2 = drawsvg(ceftazidime, 250, 250, highlight=subg2)

display("text/html", displayimgpair(svg1, svg2))
println("Status: $(status)")
println("MCES size: $(length(mapping))")
```

- As disconnected MCS takes very long time,  keyword argument`timeout` was set to 10 second to abort MCS calculation.
- If it timed out, the status will be `timedout` and the result will be suboptimal common substructure size calculated so far.

```julia
using Profile
@profile mcesmol(cefditoren, ceftazidime, timeout=10)
Profile.print(mincount=1000)
```

- The most costful call was `maximalcliques` which yields maximal cliques. Maximum clique detection problem is known to be NP-hard.

```julia
(mapping, status) = mcesmol(cefditoren, ceftazidime, connected=true)

subg1 = edgesubgraph(cefditoren, Set(keys(mapping)))
svg1 = drawsvg(cefditoren, 250, 250, highlight=subg1)

subg2 = edgesubgraph(ceftazidime, Set(values(mapping)))
svg2 = drawsvg(ceftazidime, 250, 250, highlight=subg2)

display("text/html", displayimgpair(svg1, svg2))
println("Status: $(status)")
println("Connected MCES size: $(length(mapping))")
```

- Connected MCS is far faster than disconnected MCS.

```julia
(mapping, status) = mcesmol(cefditoren, ceftazidime, targetsize=15)

subg1 = edgesubgraph(cefditoren, Set(keys(mapping)))
svg1 = drawsvg(cefditoren, 250, 250, highlight=subg1)

subg2 = edgesubgraph(ceftazidime, Set(values(mapping)))
svg2 = drawsvg(ceftazidime, 250, 250, highlight=subg2)

display("text/html", displayimgpair(svg1, svg2))
println("Status: $(status)")
println("Connected MCES size: $(length(mapping))")
```

- Or you can use `targetsize` keyword argument to know that these compounds have at least the given size of the common substructure. This feature is quite useful in library screening.


## Topological constraint

- Disconnected MCS methods can detect as many matched fragments as possible, but it does not reflect spatial relationship of each fragments.
- On the other hand, connected MCS is too strict not to allow distant matches.
- Graph distance function can be used as a node attribute matcher to constrain the spatial arrangement of matched substructure fragments (topological constraint). This feature seems to be preferable in pharmacophore matching and structural similarity screening.
- Topological constraint also effectively decreases the calculation cost.
- Topological constraint can be used by passing keyword argument `topological=True`.
- Distance mismatch tolerance parameter $\theta$ is also available as the keyword argument `tolerance`.

```julia
(mapping, status) = mcesmol(cefditoren, ceftazidime, topological=true)

subg1 = edgesubgraph(cefditoren, Set(keys(mapping)))
svg1 = drawsvg(cefditoren, 250, 250, highlight=subg1)

subg2 = edgesubgraph(ceftazidime, Set(values(mapping)))
svg2 = drawsvg(ceftazidime, 250, 250, highlight=subg2)

display("text/html", displayimgpair(svg1, svg2))
println("Status: $(status)")
println("tcMCES size: $(length(mapping))")
```
