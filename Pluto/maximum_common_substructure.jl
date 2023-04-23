### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 8b4b7b1a-8f77-4917-b55c-3f32ddd3bb4f
using Graphs, MolecularGraph


# ╔═╡ 4954986d-31e0-40f1-aa2a-33cca57273b4
using Profile

# ╔═╡ c324c158-d7fa-11ed-0300-93b6975d9deb
md"
# Maximum common substructure

MolecularGraph.jl version: 0.14.0

MolecularGraph.js implements following essential maximum common substructure (MCS) methods for cheminformatics applications.

- Maximum common induced substructure (MCIS)
- Maximum common edge-induced substructure (MCES)
- Connected or disconnected MCS
- Working with larger molecules
- Topological constraint (known as TD-MCS)
"

# ╔═╡ 62d3c057-1c65-44fd-ac25-266bdb26e536
data_dir = let
	# Create data directory
	data_dir = "_data"
	isdir(data_dir) || mkdir(data_dir)
	data_dir
end

# ╔═╡ 29a84294-8ab8-4d93-9d05-adcf67f6c57d
function fetch_mol!(cid, name, datadir)
	url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/$(cid)/SDF"
	dest = joinpath(data_dir, "$(name).mol")
	isfile(dest) || download(url, dest);
	return dest
end

# ╔═╡ 65b40864-4cbc-4122-9b58-cfbecdee34ee
md"
### Maximum common induced substructure (MCIS)

- `disconnected_mcis(mol1, mol2)` calculates maximum common induced (node-induced) substructure (MCIS).
- 'induced' (or specifically 'node-induced') means that the subgraph is defined by a subset of the vertices of the original graph, and the presence or absence of edges between all vertices of the subgraph must match that of the original graph. That also means there is at least one vertex mapping (injection) from the subgraph to the original graph.
- MCS methods in this tutorial returns a `MCSResult` type object.
  - `MCSResult.mapping` is a mapping of vertex indices between two molecules.
  - `MCSResult.status` shows the calculation status. `:done` means that the calculation was successfully finished.
- Unlike substructure match, `MCSResult` have only one of the largest vertex mapping. This is because of the high computational cost of MCS calculation as described below.
- Like substructure match, MCS calculation considers the match between of the atom symbol and the number of pi electron on the atom (not the bond order). This is why the nitrogen atom in the following example does not match (N adjacent to multiple or aromatic bond is treated as sp2 in the default descriptor function).
"

# ╔═╡ 18e3d457-bd9e-418f-a08e-25eefcd02cd0
ldopa = let
	mol = sdftomol(fetch_mol!("6047", "Levodopa", data_dir))
	remove_all_hydrogens!(mol)
	mol
end; nothing

# ╔═╡ 0c2ddc86-71ac-43d1-ae86-71510b216d92
aminocoumarin = let
	mol = sdftomol(fetch_mol!("74217", "3-Aminocoumarin", data_dir))
	remove_all_hydrogens!(mol)
	mol
end; nothing

# ╔═╡ 8b645955-b2ce-4b2a-9793-18a1612e2a75
mcis_result = disconnected_mcis(ldopa, aminocoumarin)

# ╔═╡ 7b38cdb4-7afd-4957-b6d3-116b4bd90d96
let
	nodes1 = collect(keys(mcis_result.mapping))
	# induce MCS edges from the vertex mapping
	edges1 = induced_subgraph_edges(ldopa.graph, nodes1)
	svg1 = drawsvg(ldopa,
		atomhighlight=nodes1, bondhighlight=edges1, atomindex=true)

	nodes2 = collect(values(mcis_result.mapping))
	edges2 = induced_subgraph_edges(aminocoumarin.graph, nodes2)
	svg2 = drawsvg(aminocoumarin,
		atomhighlight=nodes2, bondhighlight=edges2, atomindex=true)

	html_grid([svg1, svg2], 3, 250)
end

# ╔═╡ 154d878f-230d-4e1a-a8b1-f261e226594a
md"
### Maximum common edge-induced substructure (MCES)

- `disconnected_mces(mol1, mol2)` calculates maximum common edge-induced substructure (MCES).
- Similar to node-induced subgraph, 'edge-induced' means that the subgraph is defined by a subset of the 'edges' of the original graph. Obviously, the vertex set of edge-induced subgraph is a subset of the vertices of the original graph. Also, there is at least one 'edge' mapping (injection) from the subgraph to the original graph.
- In many cases the MCIS size is sufficient as an indicator of the degree of structural match. On the other hand, many chemists may feel MCES is more intuitive than MCIS, so it may worth doing MCES visualization with a slightly higher computational cost.
"

# ╔═╡ 8ac8565c-1f9e-42e7-a90d-705971864955
mces_result = disconnected_mces(ldopa, aminocoumarin)

# ╔═╡ 0fbaefff-78f8-4b27-aaf4-b951fe975df2
let
	edges1 = collect(keys(mces_result.mapping))
	# induce MCS nodes from the edge mapping
	nodes1, vmap1 = induced_subgraph(ldopa.graph, edges1)
	svg1 = drawsvg(ldopa,
		atomhighlight=vmap1[vertices(nodes1)], bondhighlight=edges1, atomindex=true)

	edges2 = collect(values(mces_result.mapping))
	nodes2, vmap2 = induced_subgraph(aminocoumarin.graph, edges2)
	svg2 = drawsvg(aminocoumarin,
		atomhighlight=vmap2[vertices(nodes2)], bondhighlight=edges2, atomindex=true)

	html_grid([svg1, svg2], 3, 250)
end

# ╔═╡ 3c7be77e-6116-4c93-b3c0-a058cdd20a52
md"
### Connected or disconnected MCS

- `connected_mcis(mol1, mol2)` and `connected_mces(mol1, mol2)` also calculates MCIS and MCES but with connection constraint (obtained MCS must be a connected component).
- Connected MCS is much less computationally expensive than disconnected MCS.
- As disconnected MCS does not care the spatial relationship (distance) among matched fragments, connected MCS can be more appropriate option for some kind of applications.
"

# ╔═╡ 0fb4d717-6422-4e1b-a345-a5841d164d9f
conn_mcis_result = connected_mcis(ldopa, aminocoumarin)

# ╔═╡ a0ae8cb9-a5b7-4e8d-834b-3bea0bbcc9de
let
	nodes1 = collect(keys(conn_mcis_result.mapping))
	edges1 = induced_subgraph_edges(ldopa.graph, nodes1)
	svg1 = drawsvg(ldopa,
		atomhighlight=nodes1, bondhighlight=edges1, atomindex=true)

	nodes2 = collect(values(conn_mcis_result.mapping))
	edges2 = induced_subgraph_edges(aminocoumarin.graph, nodes2)
	svg2 = drawsvg(aminocoumarin,
		atomhighlight=nodes2, bondhighlight=edges2, atomindex=true)

	html_grid([svg1, svg2], 3, 250)
end

# ╔═╡ 8b4b896d-f6e7-427c-84f5-0bcc44eb5fc9
conn_mces_result = connected_mces(ldopa, aminocoumarin)

# ╔═╡ 2cdd67cb-35d0-4bab-b8b9-5ce37552398e
let
	edges1 = collect(keys(conn_mces_result.mapping))
	nodes1, vmap1 = induced_subgraph(ldopa.graph, edges1)
	svg1 = drawsvg(ldopa,
		atomhighlight=vmap1[vertices(nodes1)], bondhighlight=edges1, atomindex=true)

	edges2 = collect(values(conn_mces_result.mapping))
	nodes2, vmap2 = induced_subgraph(aminocoumarin.graph, edges2)
	svg2 = drawsvg(aminocoumarin,
		atomhighlight=vmap2[vertices(nodes2)], bondhighlight=edges2, atomindex=true)

	html_grid([svg1, svg2], 3, 250)
end

# ╔═╡ b7a007c9-8b18-428b-a48a-a123dfcd2c8a
md"
### Working with larger molecules

Then, let's try it with a little larger molecules.
"

# ╔═╡ f9f50f92-d79e-42e2-a6cf-8107cd715f46
cefditoren = let
	mol = sdftomol(fetch_mol!("6437877", "Cefditoren Pivoxil", data_dir))
	remove_all_hydrogens!(mol)
	mol
end; nothing

# ╔═╡ c2ca480e-da0a-45b2-9ea2-4567d6655769
ceftazidime = let
	mol = sdftomol(fetch_mol!("5481173", "Ceftazidime", data_dir))
	remove_all_hydrogens!(mol)
	mol
end; nothing

# ╔═╡ 0689eeb0-ab08-47e4-b163-68d907eb8490
larger_mces = disconnected_mces(cefditoren, ceftazidime, timeout=10)

# ╔═╡ 268f7f7c-ae48-4b4a-898e-56f9aa893c46
md"
- As disconnected MCS takes very long time, keyword argument `timeout=10` was set to interrupt MCS calculation in 10 seconds.
- If it timed out, `MCSResult.status` will be `:timedout` and the result will be the mapping of suboptimal common substructure calculated so far.
"

# ╔═╡ fd1152c8-cbb4-4f92-aa5a-0d01d483df75
let
	edges1 = collect(keys(larger_mces.mapping))
	nodes1, vmap1 = induced_subgraph(cefditoren.graph, edges1)
	svg1 = drawsvg(cefditoren,
		atomhighlight=vmap1[vertices(nodes1)], bondhighlight=edges1)

	edges2 = collect(values(larger_mces.mapping))
	nodes2, vmap2 = induced_subgraph(ceftazidime.graph, edges2)
	svg2 = drawsvg(ceftazidime,
		atomhighlight=vmap2[vertices(nodes2)], bondhighlight=edges2)

	html_grid([svg1, svg2], 3, 250)
end

# ╔═╡ 1f3c9952-28bb-4e69-b84b-e6cd68ce2263
let
	@profile disconnected_mces(cefditoren, ceftazidime, timeout=10)
	Profile.print(mincount=1000)
end

# ╔═╡ 43b2240e-5c85-4b4e-8a3c-9c361b188f8d
md"
- The most costful call was `find_cliques` which yields maximal cliques. Maximum clique detection problem is known to be NP-hard. 
- Connected MCS is far faster than disconnected MCS.
"

# ╔═╡ 22766969-41f0-44c2-974f-1650e91e4fb1
larger_conn_mces = connected_mces(cefditoren, ceftazidime, timeout=10)

# ╔═╡ f0d09a83-1e9b-40bf-b83f-fc6855929f73
let
	edges1 = collect(keys(larger_conn_mces.mapping))
	nodes1, vmap1 = induced_subgraph(cefditoren.graph, edges1)
	svg1 = drawsvg(cefditoren,
		atomhighlight=vmap1[vertices(nodes1)], bondhighlight=edges1)

	edges2 = collect(values(larger_conn_mces.mapping))
	nodes2, vmap2 = induced_subgraph(ceftazidime.graph, edges2)
	svg2 = drawsvg(ceftazidime,
		atomhighlight=vmap2[vertices(nodes2)], bondhighlight=edges2)

	html_grid([svg1, svg2], 3, 250)
end

# ╔═╡ a2ca38e3-858e-4b7d-93f8-5662e4aeb220
md"
- Or you can use `targetsize` keyword argument to terminate the calculation with `:targetreached` status when the MCS reaches the target size. This feature is quite useful in library screening.
"

# ╔═╡ 8ac75243-2817-4df4-aad6-3cbf29dc098b
targeted_mces = disconnected_mces(cefditoren, ceftazidime, targetsize=15)

# ╔═╡ 3a4689f2-78fe-4f87-b31c-724a903984e8
let
	edges1 = collect(keys(targeted_mces.mapping))
	nodes1, vmap1 = induced_subgraph(cefditoren.graph, edges1)
	svg1 = drawsvg(cefditoren,
		atomhighlight=vmap1[vertices(nodes1)], bondhighlight=edges1)

	edges2 = collect(values(targeted_mces.mapping))
	nodes2, vmap2 = induced_subgraph(ceftazidime.graph, edges2)
	svg2 = drawsvg(ceftazidime,
		atomhighlight=vmap2[vertices(nodes2)], bondhighlight=edges2)

	html_grid([svg1, svg2], 3, 250)
end

# ╔═╡ b2133510-d935-41db-aee5-162509621766
md"
### Topological constraint (known as TD-MCS)

- disconnected MCS methods can detect as many matched fragments as possible, but it does not reflect spatial relationship of each fragments.
- on the other hand, connected MCS is too strict not to allow distant matches.
- graph distance function can be used as a node attribute matcher to constrain the spatial arrangement of matched substructure fragments (topological constraint). This feature seems to be preferable in pharmacophore matching and structural similarity screening.
- topological constraint also effectively decreases the calculation cost.
- topological constraint can be used by passing keyword argument `topological=True`.
- distance mismatch tolerance parameter ``θ`` is also available as the keyword argument `tolerance`
"

# ╔═╡ 53441110-c544-403e-8bc7-dd1302470ec8
tcmces_result = tcmces(cefditoren, ceftazidime, tolerance=1)

# ╔═╡ e5ddd0a1-bf60-49ad-909f-1b91cf09cc46
let
	edges1 = collect(keys(tcmces_result.mapping))
	nodes1, vmap1 = induced_subgraph(cefditoren.graph, edges1)
	svg1 = drawsvg(cefditoren,
		atomhighlight=vmap1[vertices(nodes1)], bondhighlight=edges1)

	edges2 = collect(values(tcmces_result.mapping))
	nodes2, vmap2 = induced_subgraph(ceftazidime.graph, edges2)
	svg2 = drawsvg(ceftazidime,
		atomhighlight=vmap2[vertices(nodes2)], bondhighlight=edges2)

	html_grid([svg1, svg2], 3, 250)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Graphs = "86223c79-3864-5bf0-83f7-82e725a168b6"
MolecularGraph = "6c89ec66-9cd8-5372-9f91-fabc50dd27fd"
Profile = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[compat]
Graphs = "~1.8.0"
MolecularGraph = "~0.14.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.2"
manifest_format = "2.0"
project_hash = "54652fa11b1985b020edd024aea6a367fffcb361"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "cc37d689f599e8df4f464b2fa3870ff7db7492ef"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.6.1"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "7a60c856b9fa189eb34f5f8a6f6b5529b7942957"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.6.1"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.5.2+0"

[[deps.DataAPI]]
git-tree-sha1 = "e8119c1a33d267e16108be441a287a6981ba1630"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.14.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.Extents]]
git-tree-sha1 = "5e1e4c53fa39afe63a7d356e30452249365fba99"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.1"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "1cd7f0af1aa58abc02ea1d872953a97359cb87fa"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.4"

[[deps.GeoInterface]]
deps = ["Extents"]
git-tree-sha1 = "0eb6de0b312688f852f347171aba888658e29f20"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.3.0"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "GeoInterface", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "303202358e38d2b01ba46844b92e48a3c238fd9e"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.6"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "1cf1d7dcb4bc32d7b4a5add4232db3750c27ecb4"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.8.0"

[[deps.Inflate]]
git-tree-sha1 = "5cd07aab533df5170988219191dfad0519391428"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.3"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LazyJSON]]
deps = ["JSON", "OrderedCollections", "PropertyDicts"]
git-tree-sha1 = "ce08411caa70e0c9e780f142f59debd89a971738"
uuid = "fc18253b-5e1b-504c-a4a2-9ece4944c004"
version = "0.2.2"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c7cb1f5d892775ba13767a87c7ada0b980ea0a71"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+2"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.MakieCore]]
deps = ["Observables"]
git-tree-sha1 = "9926529455a331ed73c19ff06d16906737a876ed"
uuid = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
version = "0.6.3"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MolecularGraph]]
deps = ["Colors", "Dates", "DelimitedFiles", "GeometryBasics", "Graphs", "JSON", "LinearAlgebra", "MakieCore", "Printf", "Statistics", "Unmarshal", "YAML", "coordgenlibs_jll", "libinchi_jll"]
git-tree-sha1 = "18fa7feaca1d4af33eed1e044357327bdba7f1ce"
uuid = "6c89ec66-9cd8-5372-9f91-fabc50dd27fd"
version = "0.14.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Nullables]]
git-tree-sha1 = "8f87854cc8f3685a60689d8edecaa29d2251979b"
uuid = "4d1e1d77-625e-5b40-9113-a560ec7a8ecd"
version = "1.0.0"

[[deps.Observables]]
git-tree-sha1 = "6862738f9796b3edc1c09d0890afce4eca9e7e93"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.4"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "d321bf2de576bf25ec4d3e4360faca399afca282"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.0"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "478ac6c952fddd4399e71d4779797c538d0ff2bf"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.8"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[[deps.PropertyDicts]]
git-tree-sha1 = "8cf3b5cea994cfa9f238e19c3946a39cf051896c"
uuid = "f8a19df8-e894-5f55-a973-672c1158cbca"
version = "0.1.2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "63e84b7fdf5021026d0f17f76af7c57772313d99"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.21"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StringEncodings]]
deps = ["Libiconv_jll"]
git-tree-sha1 = "33c0da881af3248dafefb939a21694b97cfece76"
uuid = "69024149-9ee7-55f6-a4c4-859efe599b68"
version = "0.3.6"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "GPUArraysCore", "StaticArraysCore", "Tables"]
git-tree-sha1 = "521a0e828e98bb69042fec1809c1b5a680eb7389"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.15"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "1544b926975372da01227b382066ab70e574a3ec"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Unmarshal]]
deps = ["JSON", "LazyJSON", "Missings", "Nullables", "Requires"]
git-tree-sha1 = "ee46863309f8f942249e1df1b74ba3088ff0f151"
uuid = "cbff2730-442d-58d7-89d1-8e530c41eb02"
version = "0.4.4"

[[deps.YAML]]
deps = ["Base64", "Dates", "Printf", "StringEncodings"]
git-tree-sha1 = "dbc7f1c0012a69486af79c8bcdb31be820670ba2"
uuid = "ddb6d928-2868-570f-bddf-ab3f9cf99eb6"
version = "0.4.8"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.coordgenlibs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "8a0fdb746dfc75758d0abea3196f5edfcbbebd79"
uuid = "f6050b86-aaaf-512f-8549-0afff1b4d57f"
version = "3.0.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.libinchi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "034ee07d3b387a4ca1a153a43a0c46549b6749ba"
uuid = "172afb32-8f1c-513b-968f-184fcd77af72"
version = "1.5.1+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╟─c324c158-d7fa-11ed-0300-93b6975d9deb
# ╠═8b4b7b1a-8f77-4917-b55c-3f32ddd3bb4f
# ╠═62d3c057-1c65-44fd-ac25-266bdb26e536
# ╠═29a84294-8ab8-4d93-9d05-adcf67f6c57d
# ╟─65b40864-4cbc-4122-9b58-cfbecdee34ee
# ╠═18e3d457-bd9e-418f-a08e-25eefcd02cd0
# ╠═0c2ddc86-71ac-43d1-ae86-71510b216d92
# ╠═8b645955-b2ce-4b2a-9793-18a1612e2a75
# ╠═7b38cdb4-7afd-4957-b6d3-116b4bd90d96
# ╟─154d878f-230d-4e1a-a8b1-f261e226594a
# ╠═8ac8565c-1f9e-42e7-a90d-705971864955
# ╠═0fbaefff-78f8-4b27-aaf4-b951fe975df2
# ╟─3c7be77e-6116-4c93-b3c0-a058cdd20a52
# ╠═0fb4d717-6422-4e1b-a345-a5841d164d9f
# ╠═a0ae8cb9-a5b7-4e8d-834b-3bea0bbcc9de
# ╠═8b4b896d-f6e7-427c-84f5-0bcc44eb5fc9
# ╠═2cdd67cb-35d0-4bab-b8b9-5ce37552398e
# ╟─b7a007c9-8b18-428b-a48a-a123dfcd2c8a
# ╠═f9f50f92-d79e-42e2-a6cf-8107cd715f46
# ╠═c2ca480e-da0a-45b2-9ea2-4567d6655769
# ╠═0689eeb0-ab08-47e4-b163-68d907eb8490
# ╟─268f7f7c-ae48-4b4a-898e-56f9aa893c46
# ╠═fd1152c8-cbb4-4f92-aa5a-0d01d483df75
# ╠═4954986d-31e0-40f1-aa2a-33cca57273b4
# ╠═1f3c9952-28bb-4e69-b84b-e6cd68ce2263
# ╟─43b2240e-5c85-4b4e-8a3c-9c361b188f8d
# ╠═22766969-41f0-44c2-974f-1650e91e4fb1
# ╠═f0d09a83-1e9b-40bf-b83f-fc6855929f73
# ╟─a2ca38e3-858e-4b7d-93f8-5662e4aeb220
# ╠═8ac75243-2817-4df4-aad6-3cbf29dc098b
# ╠═3a4689f2-78fe-4f87-b31c-724a903984e8
# ╟─b2133510-d935-41db-aee5-162509621766
# ╠═53441110-c544-403e-8bc7-dd1302470ec8
# ╠═e5ddd0a1-bf60-49ad-909f-1b91cf09cc46
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
