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
### Substructure match

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
### InChI and InChIKey

- `inchi(mol)` and `inchikey(mol)` are available to generate InChI and InChIKey by using InChI library (libinchi) published by InChI Trust, IUPAC.
- InChI and InChIKey provides structure-based identifier to chemical database and enables efficient exact structure match. See InChI Trust technical FAQ for usage and limitation of InChI and InChIKey [https://www.inchi-trust.org/technical-faq/](https://www.inchi-trust.org/technical-faq/)
"

# ╔═╡ d03046ef-5cdc-4bae-a9a1-ab3f282d305b
inchi(thiostrepton)

# ╔═╡ 50e3e3e2-aa10-432b-a8db-e9c3fe8298d1
inchikey(thiostrepton)

# ╔═╡ 76f3acc5-3209-422c-9d67-8db8b24dd0c3
md"
### SMARTS query

SMARTS queries also can be used as an argument of substructure match methods.
"

# ╔═╡ 4e6c750d-2bf9-4f2b-ac28-8b2aeb1cb0d6
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

# ╔═╡ a649fb33-53e2-4e93-8267-be2a595a640f
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

# ╔═╡ 760dca56-97a5-4074-aa1f-b7bdd7a38e2c
has_substruct_match(thiostrepton, smartstomol(raw"[$(Cc1ncsc1)](=O)N"))

# ╔═╡ ab5213fb-da1f-453e-9db1-1972b829da9c
md"
### Functional group analysis


"

# ╔═╡ 19feb1b5-f810-4d84-92c3-4425ca157150
md"
### Structural alerts


"

# ╔═╡ 87b379c1-0740-4d63-b9d5-13db607d4e61
md"
### Query containment (experimental)

"

# ╔═╡ c2e69ed7-4f54-447b-b592-b65f7a2c567e
mol = smartstomol("CCC(C)C1C(=O)NC(C(=O)NC(=C)C(=O)NC(C(=O)NC23CCC(=NC2C4=CSC(=N4)C(C(OC(=O)C5=NC6=C(C=CC(C6O)N1)C(=C5)C(C)O)C)NC(=O)C7=CSC(=N7)C(NC(=O)C8CSC(=N8)C(=CC)NC(=O)C(NC(=O)C9=CSC3=N9)C(C)O)C(C)(C(C)O)O)C1=NC(=CS1)C(=O)NC(=C)C(=O)NC(=C)C(=O)N)C)C");nothing

# ╔═╡ 5c5673fc-14b6-4f25-85d2-bb4e7aaa434b
sssr(mol)

# ╔═╡ d67a7898-4574-422c-8af6-4e148951001b
q1 = smartstomol("C")

# ╔═╡ dc93d51e-cb6c-44e5-9f40-4cdb689146f1
q2 = smartstomol(raw"[$([Ru]),$([Rh]),$([Se]),$([Pd]),$([Sc]),$([Bi]),$([Sb]),$([Ag]),$([Ti]),$([Al]),$([Cd]),$([V]),$([In]),$([Cr]),$([Sn]),$([Mn]),$([La]),$([Fe]),$([Er]),$([Tm]),$([Yb]),$([Lu]),$([Hf]),$([Ta]),$([W]),$([Re]),$([Co]),$([Os]),$([Ni]),$([Ir]),$([Cu]),$([Zn]),$([Ga]),$([Ge]),$([As]),$([Y]),$([Zr]),$([Nb]),$([Ce]),$([Pr]),$([Nd]),$([Sm]),$([Eu]),$([Gd]),$([Tb]),$([Dy]),$([Ho]),$([Pt]),$([Au]),$([Hg]),$([Tl]),$([Pb]),$([Ac]),$([Th]),$([Pa]),$([Mo]),$([U]),$([Tc]),$([Te]),$([Po]),$([At])]")

# ╔═╡ 601f71ad-0d74-4708-9200-0834ce04aa8f
MolecularGraph.optimize_query(q1.vprops[1])

# ╔═╡ bcd4013f-e748-4154-8491-db66e8a8e0b9
MolecularGraph.optimize_query(q2.vprops[1])

# ╔═╡ 98b38ccb-6a09-4aec-abd2-284c1793fc42
MolecularGraph.resolve_recursive(
	q1.vprops[1].tree,
	MolecularGraph.querypropmap(q1.vprops[1].tree)
)

# ╔═╡ 27839cba-d644-4180-bfbb-ba60c72038b7
MolecularGraph.resolve_recursive(
	q2.vprops[1].tree,
	MolecularGraph.querypropmap(q2.vprops[1].tree)
)

# ╔═╡ Cell order:
# ╟─3d9ec0c8-d7f8-11ed-1e08-0983ce040e4e
# ╠═7f5ac271-a399-4ca2-a6e6-4ce08277ce41
# ╠═22c15ca6-784a-428f-b851-890f7146ff52
# ╠═c371014c-2cc6-4542-9061-29bf6fade656
# ╠═98d11d83-01c7-4e5c-8a44-6f3d06224985
# ╠═dbb371af-fd50-48dd-b653-87095b6c2acc
# ╠═13dea6cc-6626-4a75-aeca-20b08908e8d9
# ╠═a3b9fc71-3a64-44e3-8ba3-42ed860e3d8c
# ╠═e24b5427-5f69-488c-ab76-7cd7c510f4b9
# ╠═d3ad2ec6-1f02-4664-a9ed-ea1e79de7e40
# ╟─b632e57f-6ee8-424c-a6c0-d75d0a285769
# ╠═58584903-f737-4654-8e48-0d3b28e19fb9
# ╠═e216758d-c84e-4cbb-9358-73e31fa53bdb
# ╠═d03046ef-5cdc-4bae-a9a1-ab3f282d305b
# ╠═50e3e3e2-aa10-432b-a8db-e9c3fe8298d1
# ╠═76f3acc5-3209-422c-9d67-8db8b24dd0c3
# ╠═4e6c750d-2bf9-4f2b-ac28-8b2aeb1cb0d6
# ╠═a649fb33-53e2-4e93-8267-be2a595a640f
# ╠═760dca56-97a5-4074-aa1f-b7bdd7a38e2c
# ╠═ab5213fb-da1f-453e-9db1-1972b829da9c
# ╠═19feb1b5-f810-4d84-92c3-4425ca157150
# ╠═87b379c1-0740-4d63-b9d5-13db607d4e61
# ╠═c2e69ed7-4f54-447b-b592-b65f7a2c567e
# ╠═5c5673fc-14b6-4f25-85d2-bb4e7aaa434b
# ╠═d67a7898-4574-422c-8af6-4e148951001b
# ╠═dc93d51e-cb6c-44e5-9f40-4cdb689146f1
# ╠═601f71ad-0d74-4708-9200-0834ce04aa8f
# ╠═bcd4013f-e748-4154-8491-db66e8a8e0b9
# ╠═98b38ccb-6a09-4aec-abd2-284c1793fc42
# ╠═27839cba-d644-4180-bfbb-ba60c72038b7
