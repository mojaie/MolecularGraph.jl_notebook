### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ e24f1a5d-40d3-497b-849f-066fb5c4fdbb
using Graphs, MolecularGraph

# ╔═╡ 7c1264b4-d720-11ed-065f-ebac300f8216
md"
# Preprocessing

MolecularGraph.jl version: 0.14.0

This tutorial includes following preprocessing strategies.

- Remove hydrogen vertices
- Extract molecules of interest
- Standardize charges
- Dealing with resonance structure
- Customize property updater
"

# ╔═╡ bc1ec157-9874-420d-b2c5-55a802a7e172
md"
Public databases (e.g. PubChem, ChEMBL) and flat file databases (e.g. SDFile) have different formats and may not always be used for your analysis as is. For example, 

- whether hydrogens are explicitly written or omitted
- whether salt and water molecules are included in the molecular graph or provided as metadata
- representation of resonance structure (e.g. diazo group; [C-]-[N+]#N <-> C-[N+]=[N-])
- charges depend on the condition - powder, dissolved or in physiological condition
"

# ╔═╡ 68e70311-ee46-488b-a075-969425c91e89
data_dir = let
	# Create data directory
	data_dir = "_data"
	isdir(data_dir) || mkdir(data_dir)
	data_dir
end

# ╔═╡ 3a65fd63-bfcd-437c-816f-a785fde8ea66
function fetch_mol!(cid, name, datadir)
	url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/$(cid)/SDF"
	dest = joinpath(data_dir, "$(name).mol")
	isfile(dest) || download(url, dest);
	return dest
end

# ╔═╡ c3f5d668-2bee-484d-aed8-27c506bf26b7
molfile = fetch_mol!("6437877", "Cefditoren Pivoxil", data_dir)

# ╔═╡ 2a21a50f-7def-4f0e-8016-ffe9210186f0
md"
### Remove hydrogen vertices

- SDFiles downloaded from PubChem have hydrogen nodes. In practice, hydrogens which is not important are removed from molecular graphs for simplicity.
- `remove_hydrogens!(mol)` removes hydrogen vertices that are not important (no charge, no unpaired electron, no specific isotope composition and not involved in stereochemistry).
- `remove_all_hydrogens!(mol)` removes all hydrogen vertices.

"

# ╔═╡ ab1f83ed-11d3-4115-95e0-d12bfed2ba4f
mol = sdftomol(molfile)

# ╔═╡ 2c4ab89b-7ddc-444f-bb14-9ac4178e0a20
let
	remove_hydrogens!(mol)
	mol
end

# ╔═╡ cfdef631-2d01-4149-b3ab-ac12d7e48952
md"
### Extract molecules of interest

`connected_components(mol)` returns connected components that are sets of vertices of the individual molecules in the molecular graph object.
"

# ╔═╡ b0e885a3-d4d0-4ebb-80fc-15ca29ecb0da
mol2 = smilestomol("CC(=O)OC1=CC=CC=C1C(=O)O.O.CCO"); nothing

# ╔═╡ 45e113f4-17d3-4a43-aa99-68e7a3a76466
html_fixed_size(drawsvg(mol2, atomindex=true), 250, 250)

# ╔═╡ dd7ac9f5-a265-49bc-a4b4-36ec5e44cbad
connected_components(mol2)

# ╔═╡ af8cfae7-4e59-4cfb-bd6d-afac7b118249
md"
- To extract the molecule of interest, you can iterate over the connected components and apply `induced_subgraph(mol, vertices)` to extract the molecules and filter them one by one.
- Or simply `extract_largest_component!(mol)` can be used. This removes vertices not belong to the largest component (connected component which has the largest number of vertices) from the graph.
"

# ╔═╡ 379a507d-798d-4bdf-8198-52654b783110
let
	mol = copy(mol2)
	extract_largest_component!(mol)
	mol
end

# ╔═╡ 33265b05-212a-4427-a7b1-fe96148f1486
md"
### Standardize charges

- `protonate_acids!(mol)` removes charges on oxo/thio acid anions
- `deprotonate_oniums!(mol)` removes charges on ammonium/oxonium cations
"

# ╔═╡ 6ed23d57-e615-4a86-a490-4911fe856528
charged = smilestomol("CCCC(=O)[O-].[N+]CCCC[N+]")

# ╔═╡ 5ebcd1ba-9242-4a18-93ae-b4efabcaf8cb
let
	mol = copy(charged)
	protonate_acids!(mol)
	deprotonate_oniums!(mol)
	mol
end

# ╔═╡ f5d8bfb0-0dc1-4aed-bccc-cbd92f64f04e
md"
### Dealing with resonance structure

- Substructure match methods in this library compares atom symbols and the number of ``\pi`` electrons, so in many cases you don't have to care about fluctuations in resonance structure.
"

# ╔═╡ 87c7cbf0-b31a-4c08-90b9-8606551c782f
quinoline1 = smilestomol("N1=CC=CC2=C1C=CC=C2")

# ╔═╡ 9ab00ca9-dde0-40a0-afbf-89da54eec0cf
quinoline2 = smilestomol("N=1C=CC=C2C1C=CC=C2")

# ╔═╡ 92c654ea-1399-49d0-b99a-e35f01145030
has_exact_match(quinoline1, quinoline2)

# ╔═╡ fbdf8bd4-144f-488f-9aa4-a4b504577567
pi_electron(quinoline1)

# ╔═╡ 54c0e4f8-978a-4662-9621-c7b2532f65a0
pi_electron(quinoline2)

# ╔═╡ 74a9168a-30c2-41bc-80a1-8f44a2896880
inchikey(quinoline1)

# ╔═╡ 3fda8505-fd5f-4f97-8f47-432465c6544b
inchikey(quinoline2)

# ╔═╡ 76639ed0-f887-40e2-bcec-db176019ed9a
md"
- Property auto-update mechanisms described later assigns double bonds and single bonds to aromatic rings which consist of SMILES small letter atoms (so called 'kekulization'). We can see how kekulization works by the script as shown below.
"

# ╔═╡ 70d15699-ce3e-4cb3-9d06-c285d7f475b2
function custom_updater!(mol)
    update_edge_rank!(mol)
    clear_caches!(mol)
    set_state!(mol, :has_updates, false)
	# skip kekulize!
end

# ╔═╡ 114518d1-8374-4ef7-9ec9-ea9326a9ab20
pyrene = smilestomol(
	"c1cc2cccc3c2c4c1cccc4cc3", config=Dict(:updater => custom_updater!))

# ╔═╡ 3b3b3da7-059b-4c20-aae6-70364d822559
let
	kekulize!(pyrene)
	pyrene
end

# ╔═╡ a1c2a97b-d6cd-4571-b50f-fb3e204eb470
md"
### Custom property updater

The following code is the default property updater for SMILES molecule

```julia
function smiles_on_update!(mol)
    update_edge_rank!(mol)
    clear_caches!(mol)
    set_state!(mol, :has_updates, false)
    # preprocessing
    kekulize!(mol)
    # recalculate bottleneck descriptors
    sssr!(mol)
    lone_pair!(mol)
    apparent_valence!(mol)
    valence!(mol)
    is_ring_aromatic!(mol)
end
```

- The first three lines in the function are required. `update_edge_rank!` recalculate the edge rank cache if the graph topology was changed (e.g. by `rem_vertex!` method). Then, clear caches and set the state `:has_updates` back to false.
- The next line (kekulize!) also should be kept. Add additional preprocessing methods after this line if you need (e.g. standardization of electric charges).  
- The lines after that are to cache bottleneck descriptors. These can be removed if these descriptors would not be used.

"

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Graphs = "86223c79-3864-5bf0-83f7-82e725a168b6"
MolecularGraph = "6c89ec66-9cd8-5372-9f91-fabc50dd27fd"

[compat]
Graphs = "~1.8.0"
MolecularGraph = "~0.14.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.2"
manifest_format = "2.0"
project_hash = "39cbfe9a3b79115b4c7c15351f8731e6f70a5385"

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
# ╟─7c1264b4-d720-11ed-065f-ebac300f8216
# ╟─bc1ec157-9874-420d-b2c5-55a802a7e172
# ╠═e24f1a5d-40d3-497b-849f-066fb5c4fdbb
# ╠═68e70311-ee46-488b-a075-969425c91e89
# ╠═3a65fd63-bfcd-437c-816f-a785fde8ea66
# ╠═c3f5d668-2bee-484d-aed8-27c506bf26b7
# ╟─2a21a50f-7def-4f0e-8016-ffe9210186f0
# ╠═ab1f83ed-11d3-4115-95e0-d12bfed2ba4f
# ╠═2c4ab89b-7ddc-444f-bb14-9ac4178e0a20
# ╟─cfdef631-2d01-4149-b3ab-ac12d7e48952
# ╠═b0e885a3-d4d0-4ebb-80fc-15ca29ecb0da
# ╠═45e113f4-17d3-4a43-aa99-68e7a3a76466
# ╠═dd7ac9f5-a265-49bc-a4b4-36ec5e44cbad
# ╟─af8cfae7-4e59-4cfb-bd6d-afac7b118249
# ╠═379a507d-798d-4bdf-8198-52654b783110
# ╟─33265b05-212a-4427-a7b1-fe96148f1486
# ╠═6ed23d57-e615-4a86-a490-4911fe856528
# ╠═5ebcd1ba-9242-4a18-93ae-b4efabcaf8cb
# ╟─f5d8bfb0-0dc1-4aed-bccc-cbd92f64f04e
# ╠═87c7cbf0-b31a-4c08-90b9-8606551c782f
# ╠═9ab00ca9-dde0-40a0-afbf-89da54eec0cf
# ╠═92c654ea-1399-49d0-b99a-e35f01145030
# ╠═fbdf8bd4-144f-488f-9aa4-a4b504577567
# ╠═54c0e4f8-978a-4662-9621-c7b2532f65a0
# ╠═74a9168a-30c2-41bc-80a1-8f44a2896880
# ╠═3fda8505-fd5f-4f97-8f47-432465c6544b
# ╟─76639ed0-f887-40e2-bcec-db176019ed9a
# ╠═70d15699-ce3e-4cb3-9d06-c285d7f475b2
# ╠═114518d1-8374-4ef7-9ec9-ea9326a9ab20
# ╠═3b3b3da7-059b-4c20-aae6-70364d822559
# ╟─a1c2a97b-d6cd-4571-b50f-fb3e204eb470
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
