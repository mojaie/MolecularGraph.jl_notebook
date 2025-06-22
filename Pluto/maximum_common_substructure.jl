### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# ╔═╡ 8b4b7b1a-8f77-4917-b55c-3f32ddd3bb4f
using Graphs, MolecularGraph

# ╔═╡ c324c158-d7fa-11ed-0300-93b6975d9deb
md"
# Maximum common substructure

MolecularGraph.jl version: 0.20.0

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
- MCS methods in this library returns a result in `Tuple{Dict,Symbol}`.
  - the first element is a mapping of vertex indices between two molecules.
  - the second element is the calculation status. `:done` means that the calculation was successfully finished.
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
mcis_result, mcis_status = disconnected_mcis(ldopa, aminocoumarin)

# ╔═╡ 7b38cdb4-7afd-4957-b6d3-116b4bd90d96
let
	nodes1 = collect(keys(mcis_result))
	# induce MCS edges from the vertex mapping
	edges1 = induced_subgraph_edges(ldopa.graph, nodes1)
	svg1 = html_fixed_size(ldopa, 250, 250,
		atomhighlight=nodes1, bondhighlight=edges1, atomindex=true)

	nodes2 = collect(values(mcis_result))
	edges2 = induced_subgraph_edges(aminocoumarin.graph, nodes2)
	svg2 = html_fixed_size(aminocoumarin, 250, 250,
		atomhighlight=nodes2, bondhighlight=edges2, atomindex=true)

	[svg1, svg2]
end

# ╔═╡ 154d878f-230d-4e1a-a8b1-f261e226594a
md"
### Maximum common edge-induced substructure (MCES)

- `disconnected_mces(mol1, mol2)` calculates maximum common edge-induced substructure (MCES).
- Similar to node-induced subgraph, 'edge-induced' means that the subgraph is defined by a subset of the 'edges' of the original graph. Obviously, the vertex set of edge-induced subgraph is a subset of the vertices of the original graph. Also, there is at least one 'edge' mapping (injection) from the subgraph to the original graph.
- In many cases the MCIS size is sufficient as an indicator of the degree of structural match. On the other hand, many chemists may feel MCES is more intuitive than MCIS, so it may worth doing MCES visualization with a slightly higher computational cost.
"

# ╔═╡ 8ac8565c-1f9e-42e7-a90d-705971864955
mces_result, mces_status = disconnected_mces(ldopa, aminocoumarin)

# ╔═╡ 0fbaefff-78f8-4b27-aaf4-b951fe975df2
let
	edges1 = collect(keys(mces_result))
	# induce MCS nodes from the edge mapping
	nodes1, vmap1 = induced_subgraph(ldopa.graph, edges1)
	svg1 = html_fixed_size(ldopa, 250, 250,
		atomhighlight=vmap1[vertices(nodes1)], bondhighlight=edges1, atomindex=true)

	edges2 = collect(values(mces_result))
	nodes2, vmap2 = induced_subgraph(aminocoumarin.graph, edges2)
	svg2 = html_fixed_size(aminocoumarin, 250, 250,
		atomhighlight=vmap2[vertices(nodes2)], bondhighlight=edges2, atomindex=true)

	[svg1, svg2]
end

# ╔═╡ 3c7be77e-6116-4c93-b3c0-a058cdd20a52
md"
### Connected or disconnected MCS

- `connected_mcis(mol1, mol2)` and `connected_mces(mol1, mol2)` also calculates MCIS and MCES but with connection constraint (obtained MCS must be a connected component). RDKit FMCS also adopts this approach.
- Connected MCS is much less computationally expensive than disconnected MCS.
- As disconnected MCS does not care the spatial relationship (distance) among matched fragments, connected MCS can be more appropriate option for some kind of applications.
"

# ╔═╡ 0fb4d717-6422-4e1b-a345-a5841d164d9f
conn_mcis_result, conn_mcis_status = connected_mcis(ldopa, aminocoumarin)

# ╔═╡ a0ae8cb9-a5b7-4e8d-834b-3bea0bbcc9de
let
	nodes1 = collect(keys(conn_mcis_result))
	edges1 = induced_subgraph_edges(ldopa.graph, nodes1)
	svg1 = html_fixed_size(ldopa, 250, 250,
		atomhighlight=nodes1, bondhighlight=edges1, atomindex=true)

	nodes2 = collect(values(conn_mcis_result))
	edges2 = induced_subgraph_edges(aminocoumarin.graph, nodes2)
	svg2 = html_fixed_size(aminocoumarin, 250, 250,
		atomhighlight=nodes2, bondhighlight=edges2, atomindex=true)

	[svg1, svg2]
end

# ╔═╡ 8b4b896d-f6e7-427c-84f5-0bcc44eb5fc9
conn_mces_result, conn_mces_status = connected_mces(ldopa, aminocoumarin)

# ╔═╡ 2cdd67cb-35d0-4bab-b8b9-5ce37552398e
let
	edges1 = collect(keys(conn_mces_result))
	nodes1, vmap1 = induced_subgraph(ldopa.graph, edges1)
	svg1 = html_fixed_size(ldopa, 250, 250,
		atomhighlight=vmap1[vertices(nodes1)], bondhighlight=edges1, atomindex=true)

	edges2 = collect(values(conn_mces_result))
	nodes2, vmap2 = induced_subgraph(aminocoumarin.graph, edges2)
	svg2 = html_fixed_size(aminocoumarin, 250, 250,
		atomhighlight=vmap2[vertices(nodes2)], bondhighlight=edges2, atomindex=true)

	[svg1, svg2]
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
larger_mces, larger_mces_status = disconnected_mces(cefditoren, ceftazidime, timeout=10)

# ╔═╡ 268f7f7c-ae48-4b4a-898e-56f9aa893c46
md"
- As disconnected MCS takes very long time, keyword argument `timeout=10` was set to interrupt MCS calculation in 10 seconds.
- If it timed out, `MCSResult.status` will be `:timedout` and the result will be the mapping of suboptimal common substructure calculated so far.
"

# ╔═╡ fd1152c8-cbb4-4f92-aa5a-0d01d483df75
let
	edges1 = collect(keys(larger_mces))
	nodes1, vmap1 = induced_subgraph(cefditoren.graph, edges1)
	svg1 = html_fixed_size(cefditoren, 250, 250,
		atomhighlight=vmap1[vertices(nodes1)], bondhighlight=edges1)

	edges2 = collect(values(larger_mces))
	nodes2, vmap2 = induced_subgraph(ceftazidime.graph, edges2)
	svg2 = html_fixed_size(ceftazidime, 250, 250,
		atomhighlight=vmap2[vertices(nodes2)], bondhighlight=edges2)

	[svg1, svg2]
end

# ╔═╡ 43b2240e-5c85-4b4e-8a3c-9c361b188f8d
md"
- By using `Profile.@profile`, we can conclude that the most costful task is `maximum_clique`.
- `maximum_clique` is used to solve maximum common subgraph problem. This problem is known to be NP-hard, means that computation time may increase exponentially with the number of nodes in the graph.
- As described above, connected MCS is far faster than disconnected MCS (but still there is a possibility of hang up).
"

# ╔═╡ 22766969-41f0-44c2-974f-1650e91e4fb1
larger_conn_mces, larger_conn_mces_status = connected_mces(cefditoren, ceftazidime, timeout=10)

# ╔═╡ f0d09a83-1e9b-40bf-b83f-fc6855929f73
let
	edges1 = collect(keys(larger_conn_mces))
	nodes1, vmap1 = induced_subgraph(cefditoren.graph, edges1)
	svg1 = html_fixed_size(cefditoren, 250, 250,
		atomhighlight=vmap1[vertices(nodes1)], bondhighlight=edges1)

	edges2 = collect(values(larger_conn_mces))
	nodes2, vmap2 = induced_subgraph(ceftazidime.graph, edges2)
	svg2 = html_fixed_size(ceftazidime, 250, 250,
		atomhighlight=vmap2[vertices(nodes2)], bondhighlight=edges2)

	[svg1, svg2]
end

# ╔═╡ a2ca38e3-858e-4b7d-93f8-5662e4aeb220
md"
- Or you can use `targetsize` keyword argument to terminate the calculation with `:targetreached` status when the MCS reaches the target size. This feature is quite useful in library screening.
"

# ╔═╡ 8ac75243-2817-4df4-aad6-3cbf29dc098b
targeted_mces, targeted_mces_status = disconnected_mces(cefditoren, ceftazidime, targetsize=15)

# ╔═╡ 3a4689f2-78fe-4f87-b31c-724a903984e8
let
	edges1 = collect(keys(targeted_mces))
	nodes1, vmap1 = induced_subgraph(cefditoren.graph, edges1)
	svg1 = html_fixed_size(cefditoren, 250, 250,
		atomhighlight=vmap1[vertices(nodes1)], bondhighlight=edges1)

	edges2 = collect(values(targeted_mces))
	nodes2, vmap2 = induced_subgraph(ceftazidime.graph, edges2)
	svg2 = html_fixed_size(ceftazidime, 250, 250,
		atomhighlight=vmap2[vertices(nodes2)], bondhighlight=edges2)

	[svg1, svg2]
end

# ╔═╡ b2133510-d935-41db-aee5-162509621766
md"
### Topological constraint (known as TD-MCS)

- disconnected MCS methods can detect as many matched fragments as possible, but it does not reflect spatial relationship of each fragments.
- on the other hand, connected MCS is too strict not to allow distant matches.
- graph distance function can be used as a node attribute matcher to constrain the spatial arrangement of matched substructure fragments (topological constraint). This feature seems to be preferable in pharmacophore matching and structural similarity screening.
- topological constraint also effectively decreases the calculation cost (but note that mathmatically this approach still have a possibility of hang up, so the default timeout is set to 10 seconds).
- topological constraint can be used by passing keyword argument `topological=True`.
- distance mismatch tolerance parameter ``θ`` is also available as the keyword argument `tolerance`
"

# ╔═╡ 53441110-c544-403e-8bc7-dd1302470ec8
tdmces_result, tdmces_status = tdmces(cefditoren, ceftazidime, tolerance=1)

# ╔═╡ e5ddd0a1-bf60-49ad-909f-1b91cf09cc46
let
	edges1 = collect(keys(tdmces_result))
	nodes1, vmap1 = induced_subgraph(cefditoren.graph, edges1)
	svg1 = html_fixed_size(cefditoren, 250, 250,
		atomhighlight=vmap1[vertices(nodes1)], bondhighlight=edges1)

	edges2 = collect(values(tdmces_result))
	nodes2, vmap2 = induced_subgraph(ceftazidime.graph, edges2)
	svg2 = html_fixed_size(ceftazidime, 250, 250,
		atomhighlight=vmap2[vertices(nodes2)], bondhighlight=edges2)

	[svg1, svg2]
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Graphs = "86223c79-3864-5bf0-83f7-82e725a168b6"
MolecularGraph = "6c89ec66-9cd8-5372-9f91-fabc50dd27fd"

[compat]
Graphs = "~1.13.0"
MolecularGraph = "~0.20.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.0-beta3"
manifest_format = "2.0"
project_hash = "2893ac0a6f9efcc7049fa731a60ca36cadf84c28"

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

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.0+1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "4e1fe97fdaed23e9dc21d4d664bea76b65fc50a0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.22"

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

[[deps.Extents]]
git-tree-sha1 = "b309b36a9e02fe7be71270dd8c0fd873625332b4"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.6"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.GeoFormatTypes]]
git-tree-sha1 = "8e233d5167e63d708d41f87597433f59a0f213fe"
uuid = "68eda718-8dee-11e9-39e7-89f7f65f511f"
version = "0.4.4"

[[deps.GeoInterface]]
deps = ["DataAPI", "Extents", "GeoFormatTypes"]
git-tree-sha1 = "294e99f19869d0b0cb71aef92f19d03649d028d5"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.4.1"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "Extents", "GeoInterface", "IterTools", "LinearAlgebra", "PrecompileTools", "Random", "StaticArrays"]
git-tree-sha1 = "2670cf32dcf0229c9893b895a9afe725edb23545"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.5.9"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "c5abfa0ae0aaee162a3fbb053c13ecda39be545b"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.13.0"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.IntervalSets]]
git-tree-sha1 = "5fbb102dcb8b1a858111ae81d56682376130517d"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.11"

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

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "a007feb38b422fbdab534406aeca1b86823cb4d6"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.11.1+1"

[[deps.LibGit2]]
deps = ["LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.9.0+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "OpenSSL_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.3+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.MakieCore]]
deps = ["ColorTypes", "GeometryBasics", "IntervalSets", "Observables"]
git-tree-sha1 = "6ec4aba99ce1a8c247822ce346b273d0d686e141"
uuid = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
version = "0.10.0"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MolecularGraph]]
deps = ["Base64", "Colors", "Dates", "DelimitedFiles", "GeometryBasics", "Graphs", "JSON", "LinearAlgebra", "MakieCore", "OrderedCollections", "Printf", "Statistics", "YAML", "coordgenlibs_jll", "libinchi_jll"]
git-tree-sha1 = "06c9b999d0476cfe3b6393272239556ac65ddf23"
uuid = "6c89ec66-9cd8-5372-9f91-fabc50dd27fd"
version = "0.20.0"

    [deps.MolecularGraph.extensions]
    CairoExt = "Cairo"
    RDKitExt = "RDKitMinimalLib"

    [deps.MolecularGraph.weakdeps]
    Cairo = "159f3aea-2a34-519c-b102-8c37f9878175"
    RDKitMinimalLib = "44044271-7623-48dc-8250-42433c44e4b7"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2025.2.25"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.0+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.12.0"

    [deps.Pkg.extensions]
    REPLExt = "REPL"

    [deps.Pkg.weakdeps]
    REPL = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "516f18f048a195409d6e072acf879a9f017d3900"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.2"

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
version = "1.12.0"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "0feb6b9031bd5c51f9072393eb5ab3efd31bf9e4"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.13"

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

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.8.3+2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

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

[[deps.YAML]]
deps = ["Base64", "Dates", "Printf", "StringEncodings"]
git-tree-sha1 = "2f58ac39f64b41fb812340347525be3b590cce3b"
uuid = "ddb6d928-2868-570f-bddf-ab3f9cf99eb6"
version = "0.4.14"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.coordgenlibs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "93fce743c1b36cde0efde11b1867cb7d16b13bf8"
uuid = "f6050b86-aaaf-512f-8549-0afff1b4d57f"
version = "3.0.2+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.12.0+0"

[[deps.libinchi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e81a593eb6a1a57d4262a3ed18a6efbc5d8ba83c"
uuid = "172afb32-8f1c-513b-968f-184fcd77af72"
version = "1.6.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.64.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.5.0+2"
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
