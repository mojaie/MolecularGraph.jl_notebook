### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 7f5ac271-a399-4ca2-a6e6-4ce08277ce41
begin
    import Pkg
	Pkg.activate(Base.current_project())
	Pkg.instantiate()
	using MolecularGraph
	using Graphs
end

# ╔═╡ 3d9ec0c8-d7f8-11ed-1e08-0983ce040e4e
md"
# Substructure and query

MolecularGraph.jl version: 0.14.0

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
	svg = drawsvg(thiostrepton, atomhighlight=nodes, bondhighlight=edges)
	html_fixed_size(svg, 300, 300)
end

# ╔═╡ 526dd590-2858-4eaa-a9a3-e77375e9ff7f
md"
### InChI and InChIKey

- `inchi(mol)` and `inchikey(mol)` are available to generate InChI and InChIKey by using InChI library (libinchi) published by InChI Trust, IUPAC.
- InChI and InChIKey provides structure-based identifier to chemical database and enables efficient exact structure match. See InChI Trust technical FAQ for usage and limitation of InChI and InChIKey [https://www.inchi-trust.org/technical-faq/](https://www.inchi-trust.org/technical-faq/)
"

# ╔═╡ 3187aaf2-0c17-40e0-909c-da5d5bf13a36
inchi(thiostrepton)

# ╔═╡ f563ce55-fcc6-4b18-9050-1617121d1c52
inchikey(thiostrepton)

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
	svg = drawsvg(thiostrepton, atomhighlight=nodes, bondhighlight=edges)
	html_fixed_size(svg, 300, 300)
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
	svg = drawsvg(thiostrepton, atomhighlight=nodes, bondhighlight=edges)
	html_fixed_size(svg, 300, 300)
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
html_fixed_size(drawsvg(curcumin, atomindex=true), 500, 300)

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
html_fixed_size(drawsvg(cefditoren, atomindex=true), 500, 300)

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
