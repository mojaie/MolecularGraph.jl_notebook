### A Pluto.jl notebook ###
# v0.20.9

using Markdown
using InteractiveUtils

# ╔═╡ 7f5ac271-a399-4ca2-a6e6-4ce08277ce41
using Graphs, MolecularGraph

# ╔═╡ 3d9ec0c8-d7f8-11ed-1e08-0983ce040e4e
md"
# Substructure and query

MolecularGraph.jl version: 0.19.0

This tutorial includes following database search and substructure mining tasks

- Substructure match
- InChI and InChIKey
- SMARTS query
- Structural alerts (e.g. PAINS)
- Functional group analysis
- Query containment
"

# ╔═╡ 22c15ca6-784a-428f-b851-890f7146ff52
data_dir = let
	# Create data directory
	data_dir = "_data"
	isdir(data_dir) || mkdir(data_dir)
	data_dir
end

# ╔═╡ 606ef609-5411-467f-8536-0d93ef2c409c
function fetch_mol!(cid, name, datadir)
	url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/$(cid)/SDF"
	dest = joinpath(data_dir, "$(name).mol")
	isfile(dest) || download(url, dest);
	return dest
end

# ╔═╡ c371014c-2cc6-4542-9061-29bf6fade656
function fetch_sdf_by_sid!(sid, name, datadir)
	# Thiostrepton cid:16154490 does not have 2D coords
	# 2D coords generation from SMILES takes long time with such macrocycle compounds
	url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/$(sid)/SDF"
	dest = joinpath(data_dir, "$(name).mol")
	isfile(dest) || download(url, dest);
	return dest
end

# ╔═╡ 43e8d598-3889-4954-a6ab-be936fca159a
md"
### Substructure match

- `has_exact_match(mol1, mol2)` returns whether the structures of `mol1` and `mol2` match exactly.
- `has_substruct_match(mol1, mol2)` returns whether `mol2` is a substructure of `mol1`
"

# ╔═╡ b175ed91-1cfd-48f4-a3e7-892be19c5661
thiostrepton = sdftomol(fetch_sdf_by_sid!("175266353", "Thiostrepton", data_dir))

# ╔═╡ d081bb98-0b70-446f-b544-870ac305fdf8
has_exact_match(thiostrepton, thiostrepton)

# ╔═╡ 82d0ca4c-96b3-4792-8cc9-d9e3947fa4b4
thiazole = smilestomol("s1cncc1")

# ╔═╡ 0c72e2ea-b6cb-46b1-bbcc-48e0bfa0a618
has_exact_match(thiostrepton, thiazole)

# ╔═╡ fe178836-f0d8-409b-baa0-61b34c257cbc
has_substruct_match(thiostrepton, thiazole)

# ╔═╡ b632e57f-6ee8-424c-a6c0-d75d0a285769
md"
- `exact_matches(mol1, mol2)` and `substruct_matches(mol1, mol2)` returns an iterator that yields vertex mappings between the molecules.
"

# ╔═╡ 69f7a580-bf68-41c6-ae96-6ecc587e058f
let
	matched_nodes = Set{Int}()
	for mapping in substruct_matches(thiostrepton, thiazole)
		push!(matched_nodes, keys(mapping)...)
	end
	nodes = collect(matched_nodes)
	edges = induced_subgraph_edges(thiostrepton.graph, nodes)
	html_fixed_size(
		thiostrepton, 300, 300, atomhighlight=nodes, bondhighlight=edges)
end

# ╔═╡ 526dd590-2858-4eaa-a9a3-e77375e9ff7f
md"
### InChI and InChIKey

- `inchi(mol)` and `inchikey(mol)` are available to generate InChI and InChIKey by using InChI library (libinchi) published by InChI Trust, IUPAC.
- InChI and InChIKey provides structure-based identifier to chemical database and enables efficient exact structure match. See InChI Trust technical FAQ for usage and limitation of InChI and InChIKey [https://www.inchi-trust.org/technical-faq/](https://www.inchi-trust.org/technical-faq/)
"

# ╔═╡ 3187aaf2-0c17-40e0-909c-da5d5bf13a36
inchi(thiazole)

# ╔═╡ f563ce55-fcc6-4b18-9050-1617121d1c52
inchikey(thiazole)

# ╔═╡ f110144e-8685-4ac7-bde4-6f31f293981d
md"
### SMARTS query

SMARTS queries also can be used as an argument of substructure match methods.
"

# ╔═╡ e8391325-b7bc-4f4f-b4d6-5637adb75775
let
	matched_nodes = Set{Int}()
	# nitrogens with one or two hydrogens connected
	for mapping in substruct_matches(
			thiostrepton, smartstomol(raw"[#7;H1,H2]"))
		push!(matched_nodes, keys(mapping)...)
	end
	nodes = collect(matched_nodes)
	edges = induced_subgraph_edges(thiostrepton.graph, nodes)
	html_fixed_size(
		thiostrepton, 300, 300, atomhighlight=nodes, bondhighlight=edges)
end

# ╔═╡ 3cd09b95-75d7-463d-a3d2-9af10f46eb78
let
	matched_nodes = Set{Int}()
	# amides next to thiazoles
	for mapping in substruct_matches(
			thiostrepton, smartstomol(raw"[$(Cc1ncsc1)](=O)N"))
		push!(matched_nodes, keys(mapping)...)
	end
	nodes = collect(matched_nodes)
	edges = induced_subgraph_edges(thiostrepton.graph, nodes)
	html_fixed_size(
		thiostrepton, 300, 300, atomhighlight=nodes, bondhighlight=edges)
end

# ╔═╡ cd47e23d-636a-453e-82eb-bf5bfa12c56d
has_substruct_match(thiostrepton, smartstomol(raw"[$(Cc1ncsc1)](=O)N"))

# ╔═╡ 19feb1b5-f810-4d84-92c3-4425ca157150
md"
### Structural alerts

- MolecularGraph.jl has default query sets for functional group detection and structural alerts (Substructures often avoided in medicinal chemistry. These are often called 'PAINS' but PAINS is one of those).
- Note that the structural alerts are not a universal structural filter. Some of them are based solely on chemists' intuition, others are compounds that interfere with only some assay systems.
- `query_containment_diagram(; sources=[])` generates query containment diagram (described later) from the default query sets which is used for the analysis.
- The structural alerts dataset is from PubChem database. Which datasets to use can be specified by `sources` keyword arguments with the following keywords that indicate the origin: `BMS`, `Dundee`, `Glaxo`, `MLSMR` and `PAINS`. If `sources` is not specified or is an empty vector, all default datasets will be used.

"

# ╔═╡ 13090468-209e-43a0-97c0-6ffb9861f88f
sa_diagram = query_containment_diagram(); nothing

# ╔═╡ 28f5b1f7-a7f3-4278-b02d-fbb0603f52e7
curcumin = smilestomol(raw"O=C(\C=C\c1ccc(O)c(OC)c1)CC(=O)\C=C\c2cc(OC)c(O)cc2");nothing

# ╔═╡ 12899136-feef-43bd-867d-42f18a0486ad
html_fixed_size(curcumin, 500, 300, atomindex=true)

# ╔═╡ d6881fe7-9f91-4ea1-b7b3-00a4121e8445
md"
- `find_queries(mol, query_diagram; sources=[])` returns matched queries in the default datasets as a subgraph of the query containment diagram.
"

# ╔═╡ 43080b6d-d0a6-49d3-b8fc-c02daae5bab6
curcumin_sa = find_queries(curcumin, sa_diagram; sources=["BMS", "Dundee", "Glaxo", "MLSMR", "PAINS"]); nothing

# ╔═╡ e84d1e3b-4617-49bd-938a-593c77bd24a1
let
	results = []
	for rcd in values(curcumin_sa[2])
		push!(results, Dict(
			"key" => rcd["key"],
			"name" => rcd["info"][1]["name"],
			"query" => rcd["info"][1]["query"],
			"source" => rcd["info"][1]["source"],
			"matched" => rcd["matched"]  # matched vertex indices
		))
	end
	results
end

# ╔═╡ 1f5cf7ec-110e-40be-b5ef-a5b6ae33383c
md"
### Functional group analysis

- Manually curated default functional group datasets are available for efficient caluclation and systematic annotaton of functional groups. These datasets can be specified in `sources` keyword of the `query_containment_diagram` method with the following keywords: `default`, `scaffold`, `biomolecule`
"

# ╔═╡ 50720c92-2082-45b5-af16-4c6b04b2ed67
fg_diagram = query_containment_diagram(sources=["default", "scaffold", "biomolecule"]); nothing

# ╔═╡ f9097463-464a-47ee-aaa3-5501513c044b
cefditoren = let
	molfile = fetch_mol!("6437877", "Cefditoren Pivoxil", data_dir)
	mol = sdftomol(molfile)
	remove_hydrogens!(mol)
	mol
end; nothing

# ╔═╡ 61062c66-e4e0-4177-a2cd-21b84ac6c572
html_fixed_size(cefditoren, 500, 300, atomindex=true)

# ╔═╡ c151cb13-b17f-46c9-be32-328c2f794da0
cefditoren_fg = find_queries(cefditoren, fg_diagram); nothing

# ╔═╡ a928a761-91eb-4416-87a6-ca4c52d9636e
let
	results = Dict[]
	for n in vertices(cefditoren_fg[1])
	    indegree(cefditoren_fg[1], n) == 0 || continue
		rcd = cefditoren_fg[2][n]
	    push!(results, Dict(
			"key" => rcd["key"],
			"name" => rcd["info"][1]["name"],
			"query" => rcd["info"][1]["query"],
			"source" => rcd["info"][1]["source"],
			"matched" => rcd["matched"]
		))
	end
	results
end

# ╔═╡ 87b379c1-0740-4d63-b9d5-13db607d4e61
md"
### Query containment

- Carbonyl group query (`[#6]=O`) always matches the molecule which matches carboxylate query (`O=CO`), similarly carboxylate query always matches the molecule which matches ester (`[$(C[#6])](=O)[$(O(C)[#6])]`). Thus, if clearly a particular query will return a subset of the search results for another query, the calculation can be omitted. MolecularGraph.jl implements a mechanism to resolve such query containment.
"

# ╔═╡ f3ec83c3-7034-4062-a41d-30c46ec5f309
has_substruct_match(
	smartstomol(raw"[$(C[#6])](=O)[$(O(C)[#6])]"),
	smartstomol("O=CO")
)

# ╔═╡ 1d0bf723-3762-43f3-bda7-369ff1b8a326
md"
##### Limitation

- In the current implementation, truthtable size and the query containment evaluation is ``O(2^n)`` where ``n`` is the number of query literals of an atom. To avoid hangup, default `maxsize` parameter of `generate_truthtable` method is set to 14.

"

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Graphs = "86223c79-3864-5bf0-83f7-82e725a168b6"
MolecularGraph = "6c89ec66-9cd8-5372-9f91-fabc50dd27fd"

[compat]
Graphs = "~1.12.1"
MolecularGraph = "~0.19.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.5"
manifest_format = "2.0"
project_hash = "0d907938e8b33d82664195bb0235a852e72bb0ca"

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
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "71aa551c5c33f1a4415867fe06b7844faadb0ae9"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.1.1"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "fde3bf89aead2e723284a8ff9cdf5b551ed700e8"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.5+0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

    [deps.ColorTypes.weakdeps]
    StyledStrings = "f489334b-da3d-4c2e-b8f0-e476e12c162b"

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
version = "1.1.1+0"

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

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d55dffd9ae73ff72f1c0482454dcf2ec6c6c4a63"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.5+0"

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

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "301b5d5d731a0654825f1f2e906990f7141a106b"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.16.0+0"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "2c5512e11c791d1baed2049c5652441b28fc6a31"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7a214fdac5ed5f59a22c2d9a885a16da1c74bbc7"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.17+0"

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

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "fee60557e4f19d0fe5cd169211fdda80e494f4e8"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.84.0+0"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "a641238db938fff9b2f60d08ed9030387daf428c"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.3"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a6dbda1fd736d60cc477d99f2e7a042acfa46e8"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.15+0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "3169fd3440a02f35e549728b0890904cfd4ae58a"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.12.1"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "f923f9a774fcf3f5cb761bfa43aeadd689714813"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.5.1+0"

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

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eb62a3deb62fc6d8822c0c4bef73e4412419c5d8"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.8+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c602b1127f4751facb671441ca72715cc95938a"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.3+0"

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
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c8da7e6a91781c41a863611c7e966098d783c57a"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.4.7+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a31572773ac1b745e0343fe5e2c8ddda7a37e997"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.41.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "321ccef73a96ba828cd51f2ab5b9f917fa73945a"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.41.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.MakieCore]]
deps = ["ColorTypes", "GeometryBasics", "IntervalSets", "Observables"]
git-tree-sha1 = "733d910c70805e7114c82508bae99c6cdf004466"
uuid = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
version = "0.9.3"

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
deps = ["Base64", "Cairo", "Colors", "Dates", "DelimitedFiles", "GeometryBasics", "Graphs", "JSON", "LinearAlgebra", "MakieCore", "OrderedCollections", "Printf", "RDKitMinimalLib", "Statistics", "YAML", "coordgenlibs_jll", "libinchi_jll"]
git-tree-sha1 = "70876530a29764becf6ae83f987de78145eaf3ef"
uuid = "6c89ec66-9cd8-5372-9f91-fabc50dd27fd"
version = "0.19.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

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
version = "0.8.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "275a9a6d85dc86c24d03d1837a0010226a96f540"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.56.3+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "db76b1ecd5e9715f3d043cec13b2ec93ce015d53"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.44.2+0"

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

[[deps.RDKitMinimalLib]]
deps = ["JSON", "RDKit_jll"]
git-tree-sha1 = "56837668e23c773b2537aceae7f3588ad4227077"
uuid = "44044271-7623-48dc-8250-42433c44e4b7"
version = "1.2.0"

[[deps.RDKit_jll]]
deps = ["Artifacts", "FreeType2_jll", "JLLWrappers", "Libdl", "Zlib_jll", "boost_jll"]
git-tree-sha1 = "37ebe7296ae1e018be4cc1abb53518bbd58a3c7a"
uuid = "03d1d220-30e6-562a-9e1a-3062d7b56d75"
version = "2022.9.5+0"

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

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

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

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "b8b243e47228b4a3877f1dd6aee0c5d56db7fcf4"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.6+1"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "b5899b25d17bf1889d25906fb9deed5da0c15b3b"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.12+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "a4c0ee07ad36bf8bbce1c3bb52d21fb1e0b987fb"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.7+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "7ed9347888fac59a618302ee38216dd0379c480d"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.12+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.YAML]]
deps = ["Base64", "Dates", "Printf", "StringEncodings"]
git-tree-sha1 = "2f58ac39f64b41fb812340347525be3b590cce3b"
uuid = "ddb6d928-2868-570f-bddf-ab3f9cf99eb6"
version = "0.4.14"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.boost_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "7a89efe0137720ca82f99e8daa526d23120d0d37"
uuid = "28df3c45-c428-5900-9ff8-a3135698ca75"
version = "1.76.0+1"

[[deps.coordgenlibs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "93fce743c1b36cde0efde11b1867cb7d16b13bf8"
uuid = "f6050b86-aaaf-512f-8549-0afff1b4d57f"
version = "3.0.2+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.libinchi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e81a593eb6a1a57d4262a3ed18a6efbc5d8ba83c"
uuid = "172afb32-8f1c-513b-968f-184fcd77af72"
version = "1.6.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "002748401f7b520273e2b506f61cab95d4701ccf"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.48+0"

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
# ╟─3d9ec0c8-d7f8-11ed-1e08-0983ce040e4e
# ╠═7f5ac271-a399-4ca2-a6e6-4ce08277ce41
# ╠═22c15ca6-784a-428f-b851-890f7146ff52
# ╠═606ef609-5411-467f-8536-0d93ef2c409c
# ╠═c371014c-2cc6-4542-9061-29bf6fade656
# ╟─43e8d598-3889-4954-a6ab-be936fca159a
# ╠═b175ed91-1cfd-48f4-a3e7-892be19c5661
# ╠═d081bb98-0b70-446f-b544-870ac305fdf8
# ╠═82d0ca4c-96b3-4792-8cc9-d9e3947fa4b4
# ╠═0c72e2ea-b6cb-46b1-bbcc-48e0bfa0a618
# ╠═fe178836-f0d8-409b-baa0-61b34c257cbc
# ╟─b632e57f-6ee8-424c-a6c0-d75d0a285769
# ╠═69f7a580-bf68-41c6-ae96-6ecc587e058f
# ╟─526dd590-2858-4eaa-a9a3-e77375e9ff7f
# ╠═3187aaf2-0c17-40e0-909c-da5d5bf13a36
# ╠═f563ce55-fcc6-4b18-9050-1617121d1c52
# ╟─f110144e-8685-4ac7-bde4-6f31f293981d
# ╠═e8391325-b7bc-4f4f-b4d6-5637adb75775
# ╠═3cd09b95-75d7-463d-a3d2-9af10f46eb78
# ╠═cd47e23d-636a-453e-82eb-bf5bfa12c56d
# ╟─19feb1b5-f810-4d84-92c3-4425ca157150
# ╠═13090468-209e-43a0-97c0-6ffb9861f88f
# ╠═28f5b1f7-a7f3-4278-b02d-fbb0603f52e7
# ╠═12899136-feef-43bd-867d-42f18a0486ad
# ╟─d6881fe7-9f91-4ea1-b7b3-00a4121e8445
# ╠═43080b6d-d0a6-49d3-b8fc-c02daae5bab6
# ╠═e84d1e3b-4617-49bd-938a-593c77bd24a1
# ╟─1f5cf7ec-110e-40be-b5ef-a5b6ae33383c
# ╠═50720c92-2082-45b5-af16-4c6b04b2ed67
# ╠═f9097463-464a-47ee-aaa3-5501513c044b
# ╠═61062c66-e4e0-4177-a2cd-21b84ac6c572
# ╠═c151cb13-b17f-46c9-be32-328c2f794da0
# ╠═a928a761-91eb-4416-87a6-ca4c52d9636e
# ╟─87b379c1-0740-4d63-b9d5-13db607d4e61
# ╠═f3ec83c3-7034-4062-a41d-30c46ec5f309
# ╟─1d0bf723-3762-43f3-bda7-369ff1b8a326
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
