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
- Functional group analysis
- Structural alert (e.g. PAINS)
- Query containment (experimental)
"

# ╔═╡ 22c15ca6-784a-428f-b851-890f7146ff52
data_dir = let
	# Create data directory
	data_dir = "_data"
	isdir(data_dir) || mkdir(data_dir)
	data_dir
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

# ╔═╡ 98d11d83-01c7-4e5c-8a44-6f3d06224985
md"
# Substructure match

- `has_exact_match(mol1, mol2)` returns whether the structures of `mol1` and `mol2` match exactly.
- `has_substruct_match(mol1, mol2)` returns whether `mol2` is a substructure of `mol1`

"

# ╔═╡ dbb371af-fd50-48dd-b653-87095b6c2acc
thiostrepton = sdftomol(fetch_sdf_by_sid!("175266353", "Thiostrepton", data_dir))

# ╔═╡ 13dea6cc-6626-4a75-aeca-20b08908e8d9
has_exact_match(thiostrepton, thiostrepton)

# ╔═╡ a3b9fc71-3a64-44e3-8ba3-42ed860e3d8c
thiazole = smilestomol("s1cncc1")

# ╔═╡ e24b5427-5f69-488c-ab76-7cd7c510f4b9
has_exact_match(thiostrepton, thiazole)

# ╔═╡ d3ad2ec6-1f02-4664-a9ed-ea1e79de7e40
has_substruct_match(thiostrepton, thiazole)

# ╔═╡ b632e57f-6ee8-424c-a6c0-d75d0a285769
md"
- `exact_matches(mol1, mol2)` and `substruct_matches(mol1, mol2)` returns an iterator that yields vertex mappings between the molecules.
"

# ╔═╡ 58584903-f737-4654-8e48-0d3b28e19fb9
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

# ╔═╡ e216758d-c84e-4cbb-9358-73e31fa53bdb
md"
# InChI and InChIKey

- `inchi(mol)` and `inchikey(mol)` are available to generate InChI and InChIKey by using InChI library (libinchi) published by InChI Trust, IUPAC.
- InChI and InChIKey provides structure-based identifier to chemical database and enables efficient exact structure match. See InChI Trust technical FAQ for usage and limitation of InChI and InChIKey [https://www.inchi-trust.org/technical-faq/](https://www.inchi-trust.org/technical-faq/)
"

# ╔═╡ d03046ef-5cdc-4bae-a9a1-ab3f282d305b
inchi(thiostrepton)

# ╔═╡ 50e3e3e2-aa10-432b-a8db-e9c3fe8298d1
inchikey(thiostrepton)

# ╔═╡ 76f3acc5-3209-422c-9d67-8db8b24dd0c3
md"
# SMARTS query

A SMARTS query also can be used as an argument of substructure match methods.
"

# ╔═╡ ab5213fb-da1f-453e-9db1-1972b829da9c
md"
# Functional group analysis


"

# ╔═╡ 87b379c1-0740-4d63-b9d5-13db607d4e61
md"
# Query containment (experimental)

"

# ╔═╡ Cell order:
# ╠═3d9ec0c8-d7f8-11ed-1e08-0983ce040e4e
# ╠═7f5ac271-a399-4ca2-a6e6-4ce08277ce41
# ╠═22c15ca6-784a-428f-b851-890f7146ff52
# ╠═c371014c-2cc6-4542-9061-29bf6fade656
# ╟─98d11d83-01c7-4e5c-8a44-6f3d06224985
# ╠═dbb371af-fd50-48dd-b653-87095b6c2acc
# ╠═13dea6cc-6626-4a75-aeca-20b08908e8d9
# ╠═a3b9fc71-3a64-44e3-8ba3-42ed860e3d8c
# ╠═e24b5427-5f69-488c-ab76-7cd7c510f4b9
# ╠═d3ad2ec6-1f02-4664-a9ed-ea1e79de7e40
# ╟─b632e57f-6ee8-424c-a6c0-d75d0a285769
# ╠═58584903-f737-4654-8e48-0d3b28e19fb9
# ╟─e216758d-c84e-4cbb-9358-73e31fa53bdb
# ╠═d03046ef-5cdc-4bae-a9a1-ab3f282d305b
# ╠═50e3e3e2-aa10-432b-a8db-e9c3fe8298d1
# ╠═76f3acc5-3209-422c-9d67-8db8b24dd0c3
# ╠═ab5213fb-da1f-453e-9db1-1972b829da9c
# ╠═87b379c1-0740-4d63-b9d5-13db607d4e61
