### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ afe097a1-6ae9-4c66-90f5-38707ae0c29c
begin
    import Pkg
	Pkg.activate(Base.current_project())
	Pkg.instantiate()
	using MolecularGraph
	using Graphs
end

# ╔═╡ 2f76d93c-221d-4c8c-a8ce-5f17b9f15cd1
md"
# Molecular graph basics

MolecularGraph.jl version: 0.14.0

This tutorial includes following molecular graph basics.

- Scope of MolecularGraph.jl
- Considerations in molecular graph implementation
- Basic operations provided by Graphs.jl interface
- MolGraph type and atom/bond properties
"

# ╔═╡ 51623203-9096-4125-a25b-1d66de219f8a
md"
### Scope of MolecularGraph.jl

MolecularGraph.jl mainly targets on small organic molecules that we deal with in medicinal chemistry. In particular, MolecularGraph.jl deals with a subset of graphs with the following characteristics:

- Undirected graph
- Simple graph, that means no self-loops and multi-edges
- Typically <100 vertices and <100 edges
- Typically degree <=4
- Almost all are planar and many of them are outerplanar
- Vertices, edges and the graph itself have multiple properties (attributes) e.g. atom symbol, bond order and metadata of the molecule
- Properties have interactions. That means changes in graph topology and vertex/edge properties can affect other vertex/edge properties.

Therefore, following molecules are not supported and may be better to use other appropriate models.
- Inorganic molecules
- Bio/synthetic polymers

"

# ╔═╡ 46afbdfc-8e40-4ef7-8693-8269819d4972
md"
### Considerations in molecular graph implementation

According to the characteristics shown above,

- Some of the graph algorithm functions in this library do not support some of general graph structure. For example, self-loops and multi-edges are not considered.
- Graph algorithm implementations for general graphs are not always optimal for molecular graphs as well. Even well-known implementations such as deposited to Graph.jl may be reimplemented depending on benchmark results.
- The properties of nodes and edges in molecular graphs are very different from the attributes like distances and labels in general graphs. At this time, MolecularGraph.jl does not depend on MetaGraphs.jl or MetaGraphsNext.jl and has its own data structure to deal with propagation and interaction of molecular graph properties.
"

# ╔═╡ b043d8ae-91b3-402f-a787-bc03d2430610
md"
### Basic operations provided by Graphs.jl interface

Main molecular graph model abstruct type `SimpleMolGraph` has `Graphs.SimpleGraph` inside to store graph structure of the molecule. Common interface of SimpleGraph can be used to access the graph structure.
"

# ╔═╡ b6252990-58a2-47c6-b8f3-bbd227e2ef86
mol = smilestomol("CC(=O)OC1=CC=CC=C1C(=O)O")

# ╔═╡ ee7d572f-fe52-4209-9757-fa62bc25f149
nv(mol)

# ╔═╡ f55360ac-cd04-44cc-8860-f1522d9bec00
ne(mol)

# ╔═╡ c5f12679-233e-416f-8883-d661cfbabee8
vertices(mol)

# ╔═╡ b35f71e3-0779-4c84-93d9-cc96785e1164
edges(mol)

# ╔═╡ db922ff9-7405-4bff-8c82-bd031ea421e2
md"
- `drawsvg` has `atomindex=true` option to display vertex indices

"

# ╔═╡ b318a824-68c0-4150-bc4a-fb9894669537
HTML(drawsvg(mol, 250, 250, atomindex=true))

# ╔═╡ 4423c072-0d94-489f-b887-3f4b0413ffa6
neighbors(mol, 2)

# ╔═╡ 190fefb6-f8f8-4924-96c9-a62d1d556bbf
has_edge(mol, 4, 5)

# ╔═╡ 44024a9f-0a57-4477-ac16-e5a4f26b34a5
md"
- `rem_vertex!(mol, v)` removes the vertex `v` and re-index the last vertex by `v` as `Graphs.rem_vertex!(g, v)` do so.
- `add_vertex!` and `add_edge!` take atom and bond properties as arguments described later. As newly added vertices and edges may not have coordinate infomation, it will be recalculated automatically.
"

# ╔═╡ 56d3c35b-d026-45d5-b3c9-bb893171069e
let
	mol = copy(mol)
	rem_vertex!(mol, 1)
	HTML(drawsvg(mol, 250, 250, atomindex=true))
end

# ╔═╡ eac02e96-e076-4c7b-9747-5eb644111ae8
let
	mol = copy(mol)
	vmap = rem_vertices!(mol, [1, 2, 3])  # remove acetyl group
	println(vmap)
	HTML(drawsvg(mol, 250, 250, atomindex=true))
end

# ╔═╡ 3a101d91-6ded-48fe-8153-9f6a0775f300
let
	mol = copy(mol)
	rem_edge!(mol, 2, 4)
	HTML(drawsvg(mol, 250, 250, atomindex=true))
end

# ╔═╡ 73e9a08c-2a18-4b83-96e8-45db09ae9dec
disconn = let
	mol = copy(mol)
	rem_edge!(mol, 2, 4)
	add_vertex!(mol, SMILESAtom(:Cl))
	add_edge!(mol, 2, 14, SMILESBond())
	mol
end

# ╔═╡ 353b1a6b-0e3c-4d58-b57c-c20b557b4412
HTML(drawsvg(disconn, 250, 250, atomindex=true))

# ╔═╡ 9a87a719-c922-4f79-848f-b3b692523d3e
md"
- 'graph' field of the molecule is used to access the graph structure (SimpleGraph) directly.
"

# ╔═╡ bcfc987e-60bb-44ad-a987-6bc5a0e4d09a
connected_components(disconn.graph)  # Graphs.connected_components

# ╔═╡ b000bf07-ff8e-4c82-9b88-fac767e4f0ba
induced_subgraph(disconn.graph, [1,2,3,14])  # Graphs.induced_subgraph

# ╔═╡ c3de9832-4175-4107-94ac-ee224e24a398
md"
- Some methods in Graphs.jl would work correctly with molecular graph objects as arguments, (e.g. connected_components(mol)), but this is not guaranteed.
- Some methods works differently depending on whether the argument is `mol` or `mol.graph`. For example, `induced_subgraph(mol.graph, nodes_or_edges)` returns a SimpleGraph object whereas `induced_subgraph(mol, nodes_or_edges)` returns a molecule object.
"

# ╔═╡ e08f56cb-7342-4f9c-b93d-577314c29b34
induced_subgraph(disconn, [1,2,3,14]) 

# ╔═╡ 9ad9e39b-e5ae-451e-be9e-3a8e1c675541
md"
### MolGraph type and atom/bond properties

- The default data type of molecules generated from SMILES, SMARTS or SDFiles is a parametric type `MolGraph{T,V,E}`
- `T` comes from SimpleGraph{T}, so this may be `Int64`
- `V` is an atom (vertex) property type. By default, this depends on the input data format (e.g. SMILESAtom, SDFAtom).
- Similarly, `E` is a bond (edge) property type.
- These type parameters can be accessed using `eltype`, `vproptype` and `eproptype`, respectively.


"

# ╔═╡ c85bf839-e205-48e9-adff-30a042551d53
typeof(mol)

# ╔═╡ ced5952a-4f99-4e48-add4-60d43e98cad1
eltype(mol)

# ╔═╡ 0c9414f5-ed5b-4310-9c62-02749c38c396
vproptype(mol)

# ╔═╡ 8f91c327-c6ca-4648-bc95-e9ebf33a6959
eproptype(mol)

# ╔═╡ 656618ef-bea2-4b35-9c8a-99771c10b379
md"
- MolecularGraph.jl has functions that can access vertex and edge properties similar to MetaGraphs.jl
"

# ╔═╡ 1308a00c-9630-4fdd-90c3-c7a81689e28a
HTML(drawsvg(mol, 250, 250, atomindex=true))

# ╔═╡ 77c4fe70-97d4-449c-bf96-eb4146af4bdb
props(mol, 5)

# ╔═╡ fd818784-3c69-4649-b7d6-3009e0f5ae22
md"
`isaromatic` here means explicit aromaticity represented as small leters in SMILES. Actual aromaticity would be calculated automatically as a descriptor array (see properties and descriptors tutorial).
"

# ╔═╡ 649e4497-7d4b-4304-a958-0e2a1e934e2b
let
	mol = smilestomol("CC(=O)Oc1ccccc1C(=O)O")
	props(mol, 5)
end

# ╔═╡ 9beb0f91-1566-4360-ba3f-1d5dd0071aa5
get_prop(mol, 12, :symbol)

# ╔═╡ 40616886-2ff7-4207-9cac-a12d00eed2f2
props(mol, 2, 3)

# ╔═╡ a8497297-b807-4630-83ea-a2a8d35e62b3
get_prop(mol, 2, 3, :order)

# ╔═╡ 0f56eff2-821f-4f4c-b0ed-e31ca19573a4
md"
- The default atom/bond types are immutable. This is intended to keep referential transparency and to obtain reproducible results.
- If you need mutable atoms and bonds, give property types explicitly as shown below. Required interface for the property type is only `getindex`, so even Dict{Symbol,Any} can be used.
"

# ╔═╡ 65798267-827d-404c-9878-57ce76f3492d
dictmol = smilestomol(MolGraph{Int64, Dict{Symbol,Any}, Dict{Symbol,Any}}, "CC(=O)OC1=CC=CC=C1C(=O)O")

# ╔═╡ 511439fb-9ae2-42ab-b13b-450ae67d2ca6
HTML(drawsvg(dictmol, 250, 250, atomindex=true))

# ╔═╡ e7ebde94-c38c-46e8-99c5-47944b844832
props(dictmol, 5)

# ╔═╡ e73c55d7-4af3-4242-91cd-2520f3aba237
md"
Let's see a case created from SDFile
"

# ╔═╡ de592d8c-7e8a-4bb5-8f95-5a978045ffa6
data_dir = let
	# Create data directory
	data_dir = "_data"
	isdir(data_dir) || mkdir(data_dir)
	data_dir
end

# ╔═╡ f6c2b832-db03-43fd-859f-5e95161f1d5e
function fetch_mol!(cid, name, datadir)
	url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/$(cid)/SDF"
	dest = joinpath(data_dir, "$(name).mol")
	isfile(dest) || download(url, dest);
	return dest
end

# ╔═╡ e635df3c-bc91-46de-b994-957df3fa54ba
molfile = fetch_mol!("6437877", "Cefditoren Pivoxil", data_dir)

# ╔═╡ a40b1869-b958-4cec-938b-b156a3dc2f07
sdfmol = let
	mol = sdftomol(molfile)
	remove_hydrogens!(mol)
	mol
end

# ╔═╡ a40609e8-2da4-4731-b133-efd75df52c3d
HTML(drawsvg(sdfmol, 350, 350, atomindex=true))

# ╔═╡ a8bab02a-c863-4a8f-8b57-0398c8ac090b
typeof(sdfmol)

# ╔═╡ a462e5b5-2978-429b-b3f9-1efa04f2b097
props(sdfmol, 14)

# ╔═╡ 8a96d876-9914-476d-821c-643788789de6
props(sdfmol, 15, 40)

# ╔═╡ 202563a2-a0f6-493e-baa4-48b531196daa
md"
- SDFile can have Metadata (called 'options'). This is stored at `MolGraph.gprops[:metadata]` and can be accessed by `metadata(mol)`.
"

# ╔═╡ 0899e5b3-bf52-4063-85d4-07253a86a09c
metadata(sdfmol)

# ╔═╡ 82984c24-4295-48bf-8358-2d6edaa2a859
metadata(sdfmol)[:PUBCHEM_COMPOUND_CID]

# ╔═╡ Cell order:
# ╟─2f76d93c-221d-4c8c-a8ce-5f17b9f15cd1
# ╟─51623203-9096-4125-a25b-1d66de219f8a
# ╟─46afbdfc-8e40-4ef7-8693-8269819d4972
# ╟─b043d8ae-91b3-402f-a787-bc03d2430610
# ╠═afe097a1-6ae9-4c66-90f5-38707ae0c29c
# ╠═b6252990-58a2-47c6-b8f3-bbd227e2ef86
# ╠═ee7d572f-fe52-4209-9757-fa62bc25f149
# ╠═f55360ac-cd04-44cc-8860-f1522d9bec00
# ╠═c5f12679-233e-416f-8883-d661cfbabee8
# ╠═b35f71e3-0779-4c84-93d9-cc96785e1164
# ╟─db922ff9-7405-4bff-8c82-bd031ea421e2
# ╠═b318a824-68c0-4150-bc4a-fb9894669537
# ╠═4423c072-0d94-489f-b887-3f4b0413ffa6
# ╠═190fefb6-f8f8-4924-96c9-a62d1d556bbf
# ╟─44024a9f-0a57-4477-ac16-e5a4f26b34a5
# ╠═56d3c35b-d026-45d5-b3c9-bb893171069e
# ╠═eac02e96-e076-4c7b-9747-5eb644111ae8
# ╠═3a101d91-6ded-48fe-8153-9f6a0775f300
# ╠═73e9a08c-2a18-4b83-96e8-45db09ae9dec
# ╠═353b1a6b-0e3c-4d58-b57c-c20b557b4412
# ╟─9a87a719-c922-4f79-848f-b3b692523d3e
# ╠═bcfc987e-60bb-44ad-a987-6bc5a0e4d09a
# ╠═b000bf07-ff8e-4c82-9b88-fac767e4f0ba
# ╟─c3de9832-4175-4107-94ac-ee224e24a398
# ╠═e08f56cb-7342-4f9c-b93d-577314c29b34
# ╟─9ad9e39b-e5ae-451e-be9e-3a8e1c675541
# ╠═c85bf839-e205-48e9-adff-30a042551d53
# ╠═ced5952a-4f99-4e48-add4-60d43e98cad1
# ╠═0c9414f5-ed5b-4310-9c62-02749c38c396
# ╠═8f91c327-c6ca-4648-bc95-e9ebf33a6959
# ╟─656618ef-bea2-4b35-9c8a-99771c10b379
# ╠═1308a00c-9630-4fdd-90c3-c7a81689e28a
# ╠═77c4fe70-97d4-449c-bf96-eb4146af4bdb
# ╟─fd818784-3c69-4649-b7d6-3009e0f5ae22
# ╠═649e4497-7d4b-4304-a958-0e2a1e934e2b
# ╠═9beb0f91-1566-4360-ba3f-1d5dd0071aa5
# ╠═40616886-2ff7-4207-9cac-a12d00eed2f2
# ╠═a8497297-b807-4630-83ea-a2a8d35e62b3
# ╟─0f56eff2-821f-4f4c-b0ed-e31ca19573a4
# ╠═65798267-827d-404c-9878-57ce76f3492d
# ╠═511439fb-9ae2-42ab-b13b-450ae67d2ca6
# ╠═e7ebde94-c38c-46e8-99c5-47944b844832
# ╟─e73c55d7-4af3-4242-91cd-2520f3aba237
# ╠═de592d8c-7e8a-4bb5-8f95-5a978045ffa6
# ╠═f6c2b832-db03-43fd-859f-5e95161f1d5e
# ╠═e635df3c-bc91-46de-b994-957df3fa54ba
# ╠═a40b1869-b958-4cec-938b-b156a3dc2f07
# ╠═a40609e8-2da4-4731-b133-efd75df52c3d
# ╠═a8bab02a-c863-4a8f-8b57-0398c8ac090b
# ╠═a462e5b5-2978-429b-b3f9-1efa04f2b097
# ╠═8a96d876-9914-476d-821c-643788789de6
# ╟─202563a2-a0f6-493e-baa4-48b531196daa
# ╠═0899e5b3-bf52-4063-85d4-07253a86a09c
# ╠═82984c24-4295-48bf-8358-2d6edaa2a859
