### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# ╔═╡ afe097a1-6ae9-4c66-90f5-38707ae0c29c
using Graphs, MolecularGraph

# ╔═╡ 2f76d93c-221d-4c8c-a8ce-5f17b9f15cd1
md"
# Molecular graph basics

MolecularGraph.jl version: 0.17.1

This tutorial includes following molecular graph basics.

- Scope of MolecularGraph.jl
- Considerations in molecular graph implementation
- Basic operations provided by Graphs.jl interface
- MolGraph type and atom/bond properties
- Molecule metadata (e.g. SDFile option block)
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
html_fixed_size(mol, 250, 250, atomindex=true)

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
	html_fixed_size(mol, 250, 250, atomindex=true)
end

# ╔═╡ eac02e96-e076-4c7b-9747-5eb644111ae8
let
	mol = copy(mol)
	vmap = rem_vertices!(mol, [1, 2, 3])  # remove acetyl group
	println(vmap)
	html_fixed_size(mol, 250, 250, atomindex=true)
end

# ╔═╡ 3a101d91-6ded-48fe-8153-9f6a0775f300
let
	mol = copy(mol)
	rem_edge!(mol, 2, 4)
	html_fixed_size(mol, 250, 250, atomindex=true)
end

# ╔═╡ 73e9a08c-2a18-4b83-96e8-45db09ae9dec
disconn = let
	mol = copy(mol)
	rem_edge!(mol, 2, 4)
	add_vertex!(mol, SMILESAtom(:Cl))
	add_edge!(mol, 2, 14, SMILESBond())
	mol
end; nothing

# ╔═╡ 353b1a6b-0e3c-4d58-b57c-c20b557b4412
html_fixed_size(disconn, 250, 250, atomindex=true)

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
html_fixed_size(mol, 250, 250, atomindex=true)

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
dictmol = smilestomol(MolGraph{Int64, Dict{Symbol,Any}, Dict{Symbol,Any}}, "CC(=O)OC1=CC=CC=C1C(=O)O"); nothing

# ╔═╡ 511439fb-9ae2-42ab-b13b-450ae67d2ca6
html_fixed_size(dictmol, 250, 250, atomindex=true)

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
end; nothing

# ╔═╡ a40609e8-2da4-4731-b133-efd75df52c3d
html_fixed_size(sdfmol, 350, 350, atomindex=true)

# ╔═╡ a8bab02a-c863-4a8f-8b57-0398c8ac090b
typeof(sdfmol)

# ╔═╡ a462e5b5-2978-429b-b3f9-1efa04f2b097
props(sdfmol, 14)

# ╔═╡ 8a96d876-9914-476d-821c-643788789de6
props(sdfmol, 15, 40)

# ╔═╡ 202563a2-a0f6-493e-baa4-48b531196daa
md"
### Molecule metadata (e.g. SDFile option block)

SDFile can have Metadata (called 'options') which is stored on `MolGraph.gprops[:metadata]`. Acccessors `has_prop`, `get_prop` and `set_prop!` methods are available for the metadata.
"

# ╔═╡ 0899e5b3-bf52-4063-85d4-07253a86a09c
has_prop(sdfmol, "PUBCHEM_COMPOUND_CID")

# ╔═╡ a0683c69-9287-4a3e-8f83-f0a8ab93dda0
get_prop(sdfmol, "PUBCHEM_COMPOUND_CID")

# ╔═╡ a7ee647f-9c34-4274-8f9d-bd010f50015b
set_prop!(sdfmol, "valid", 1)

# ╔═╡ dcbc9c40-204d-4d76-80fa-a69339527536
get_prop(sdfmol, "valid")

# ╔═╡ 97a8f1af-a1c9-47c6-96b7-874d40bf646a
md"
There are also convenient dict-like accessors.
"

# ╔═╡ 82984c24-4295-48bf-8358-2d6edaa2a859
sdfmol["PUBCHEM_COMPOUND_CID"]

# ╔═╡ bdd58925-7e6b-4de8-976c-d15ef230ce8e
sdfmol["tag"] = "test"

# ╔═╡ af0acf3a-6b53-4200-83d5-5fdecd01053f
sdfmol["tag"]

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Graphs = "86223c79-3864-5bf0-83f7-82e725a168b6"
MolecularGraph = "6c89ec66-9cd8-5372-9f91-fabc50dd27fd"

[compat]
Graphs = "~1.11.2"
MolecularGraph = "~0.17.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.0-rc1"
manifest_format = "2.0"
project_hash = "dcd5467a519a22832b1041ef7c1a94624909005c"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "d57bd3762d308bded22c3b82d033bff85f6195c6"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.4.0"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9e2a6b69137e6969bab0152632dcb3bc108c8bdd"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+1"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "d0b3f8b4ad16cb0a2988c6788646a5e6a17b6b1b"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.0.5"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a2f1c8c668c8e3cb4cca4e57a8efdb09067bb3fd"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.0+2"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "362a287c3aa50601b0bc359053d5c2468f0e7ce0"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.11"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "b1c55339b7c6c350ee89f2c1604299660525b248"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.15.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "260fd2400ed2dab602a7c15cf10c1933c59930a2"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.5"
weakdeps = ["IntervalSets", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c6317308b9dc757616f0b5cb379db10494443a7"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.2+0"

[[deps.Extents]]
git-tree-sha1 = "94997910aca72897524d2237c41eb852153b0f65"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.3"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "db16beca600632c95fc8aca29890d83788dd8b23"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.96+0"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "5c1d8ae0efc6c2e7b1fc502cbe25def8f661b7bc"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.2+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1ed150b39aebcc805c26b93a8d0122c940f64ce2"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.14+0"

[[deps.GeoInterface]]
deps = ["Extents"]
git-tree-sha1 = "9fff8990361d5127b770e3454488360443019bb3"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.3.5"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "Extents", "GeoInterface", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "b62f2b2d76cee0d61a2ef2b3118cd2a3215d3134"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.11"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "7c82e6a6cd34e9d935e9aa4051b66c6ff3af59ba"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.80.2+0"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "d61890399bc535850c4bf08e4e0d3a7ad0f21cbd"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "ebd18c326fa6cee1efb7da9a3b45cf69da2ed4d9"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.11.2"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.IntervalSets]]
git-tree-sha1 = "dba9ddf07f77f60450fe5d2e2beb9854d9a49bd0"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.10"

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

    [deps.IntervalSets.weakdeps]
    Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d986ce2d884d49126836ea94ed5bfb0f12679713"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.7+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "70c5da094887fd2cae843b8db33920bac4b6f07d"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.2+0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll"]
git-tree-sha1 = "9fd170c4bbfd8b935fdc5f8b7aa33532c991a673"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.11+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fbb1f2bef882392312feb1ede3615ddc1e9b99ed"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.49.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0c4f9c4f1a50d8f35048fa0532dabbadf702f81e"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.40.1+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "5ee6203157c120d79034c748a2acba45b82b8807"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.40.1+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.MakieCore]]
deps = ["ColorTypes", "GeometryBasics", "IntervalSets", "Observables"]
git-tree-sha1 = "c1c950560397ee68ad7302ee0e3efa1b07466a2f"
uuid = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
version = "0.8.4"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MolecularGraph]]
deps = ["Base64", "Cairo", "Colors", "Dates", "DelimitedFiles", "GeometryBasics", "Graphs", "JSON", "LinearAlgebra", "MakieCore", "Printf", "Statistics", "YAML", "coordgenlibs_jll", "libinchi_jll"]
git-tree-sha1 = "c71603f9752362e8a112f7f46ae89bb19e078069"
uuid = "6c89ec66-9cd8-5372-9f91-fabc50dd27fd"
version = "0.17.1"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "cb5a2ab6763464ae0f19c86c56c63d4a2b0f5bda"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.52.2+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "35621f10a7531bc8fa58f74610b1bfb70a3cfc6b"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.43.4+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"

    [deps.Pkg.extensions]
    REPLExt = "REPL"

    [deps.Pkg.weakdeps]
    REPL = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "eeafab08ae20c62c44c8399ccb9354a04b80db50"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.7"

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

    [deps.StaticArrays.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StringEncodings]]
deps = ["Libiconv_jll"]
git-tree-sha1 = "b765e46ba27ecf6b44faf70df40c57aa3a547dcb"
uuid = "69024149-9ee7-55f6-a4c4-859efe599b68"
version = "0.3.7"

[[deps.StructArrays]]
deps = ["ConstructionBase", "DataAPI", "Tables"]
git-tree-sha1 = "f4dc295e983502292c4c3f951dbb4e985e35b3be"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.18"

    [deps.StructArrays.extensions]
    StructArraysAdaptExt = "Adapt"
    StructArraysGPUArraysCoreExt = "GPUArraysCore"
    StructArraysSparseArraysExt = "SparseArrays"
    StructArraysStaticArraysExt = "StaticArrays"

    [deps.StructArrays.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.6.0+0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "598cd7c1f68d1e205689b1c2fe65a9f85846f297"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "d9717ce3518dc68a99e6b96300813760d887a01d"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.1+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "a54ee957f4c86b526460a720dbc882fa5edcbefc"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.41+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "d2d1a5c49fae4ba39983f63de6afcbea47194e85"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.6+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "47e45cd78224c53109495b3e324df0c37bb61fbe"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.11+0"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "bcd466676fef0878338c61e655629fa7bbc69d8e"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.YAML]]
deps = ["Base64", "Dates", "Printf", "StringEncodings"]
git-tree-sha1 = "80c3218f29cbc47111ac87e7be5e69cc05c6dd36"
uuid = "ddb6d928-2868-570f-bddf-ab3f9cf99eb6"
version = "0.4.11"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.coordgenlibs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "93fce743c1b36cde0efde11b1867cb7d16b13bf8"
uuid = "f6050b86-aaaf-512f-8549-0afff1b4d57f"
version = "3.0.2+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

[[deps.libinchi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e81a593eb6a1a57d4262a3ed18a6efbc5d8ba83c"
uuid = "172afb32-8f1c-513b-968f-184fcd77af72"
version = "1.6.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d7015d2e18a5fd9a4f47de711837e980519781a4"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.43+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

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
# ╠═a0683c69-9287-4a3e-8f83-f0a8ab93dda0
# ╠═a7ee647f-9c34-4274-8f9d-bd010f50015b
# ╠═dcbc9c40-204d-4d76-80fa-a69339527536
# ╟─97a8f1af-a1c9-47c6-96b7-874d40bf646a
# ╠═82984c24-4295-48bf-8358-2d6edaa2a859
# ╠═bdd58925-7e6b-4de8-976c-d15ef230ce8e
# ╠═af0acf3a-6b53-4200-83d5-5fdecd01053f
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
