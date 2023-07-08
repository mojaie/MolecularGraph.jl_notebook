### A Pluto.jl notebook ###
# v0.19.25

using Markdown
using InteractiveUtils

# ╔═╡ 21cea84e-126c-11ee-21b7-992b40f8d4ee
begin
	using Pkg
	Pkg.develop(path="../../MolecularGraph.jl")
	Pkg.activate(Base.current_project())
	Pkg.instantiate()
	using Profile, Graphs, MolecularGraph
	Pkg.status()
end

# ╔═╡ 3e067d65-930b-42b1-a399-7312660b4db1
let
	# 
	cefditoren = smilestomol("CC1=C(SC=N1)C=CC2=C(N3C(C(C3=O)NC(=O)C(=NOC)C4=CSC(=N4)N)SC2)C(=O)OCOC(=O)C(C)(C)C"
		)  # Cefditoren pivoxil
    ceftazidime = smilestomol("CC(C)(C(=O)O)ON=C(C1=CSC(=N1)N)C(=O)NC2C3N(C2=O)C(=C(CS3)C[N+]4=CC=CC=C4)C(=O)[O-]") # Ceftazidime
	cefditoren_ = tdmces_constraints(cefditoren)
	ceftazidime_ = tdmces_constraints(ceftazidime)
	mapping, _ = maximum_common_subgraph(cefditoren_, ceftazidime_, tolerance=1)
	#mapping, _ = tdmces(cefditoren, ceftazidime, tolerance=1)

	# atoms and bonds highlighting
	edges1 = collect(keys(mapping))
	nodes1, vmap1 = induced_subgraph(cefditoren.graph, edges1)
	svg1 = drawsvg(cefditoren, atomindex=true,
		atomhighlight=vmap1[vertices(nodes1)], bondhighlight=edges1)

	edges2 = collect(values(mapping))
	nodes2, vmap2 = induced_subgraph(ceftazidime.graph, edges2)
	svg2 = drawsvg(ceftazidime, atomindex=true,
		atomhighlight=vmap2[vertices(nodes2)], bondhighlight=edges2)
	println(length(mapping))
	html_grid([svg1, svg2], 3, 250)
end

# ╔═╡ 61294436-762a-4f02-bd9b-1ae8ecbd0c3c
let
	cefditoren = smilestomol(
		"CC1=C(SC=N1)C=CC2=C(N3C(C(C3=O)NC(=O)C(=NOC)C4=CSC(=N4)N)SC2)C(=O)OCOC(=O)C(C)(C)C")  # Cefditoren pivoxil
    ceftazidime = smilestomol(
		"CC(C)(C(=O)O)ON=C(C1=CSC(=N1)N)C(=O)NC2C3N(C2=O)C(=C(CS3)C[N+]4=CC=CC=C4)C(=O)[O-]") # Ceftazidime
	cefditoren_ = tdmcis_constraints(cefditoren)
	ceftazidime_ = tdmcis_constraints(ceftazidime)
	mapping, _ = maximum_common_subgraph(cefditoren_, ceftazidime_, tolerance=1)

	println(length(mapping))
	# atoms and bonds highlighting
	nodes1 = collect(keys(mapping))
	edges1 = induced_subgraph_edges(cefditoren.graph, nodes1)
	svg1 = drawsvg(cefditoren, atomhighlight=nodes1, bondhighlight=edges1)

	nodes2 = collect(values(mapping))
	edges2 = induced_subgraph_edges(ceftazidime.graph, nodes2)
	svg2 = drawsvg(ceftazidime, atomhighlight=nodes2, bondhighlight=edges2)

	html_grid([svg1, svg2], 3, 250)
end

# ╔═╡ 5ba7d2bd-bd30-4b97-a425-e1f67df75199
smiles = [
    "c1c2ccccc2ccc1",
    "C1=CC=C2C(=C1)C=CC3=CC=CC=C32",
    "c1c2cc3cc4ccccc4cc3cc2ccc1",
    "c1ccc2c(c1)ccc3c2ccc4c3cccc4",
    "c1cc2cccc3c2c4c1cccc4cc3",
    "c1cc2ccc3ccc4ccc5ccc6ccc1c7c2c3c4c5c67",
    "c16ccc2ccc3ccc5c4c(c1c2c34)c(cc5)cc6",
    raw"c1ccc3c6c1cccc6\C=C/2\c4cccc5cccc(\C=C\23)c45",
    "c1ccc2c(c1)c3cccc4c3c5c2cccc5c6c4cccc6"
];

# ╔═╡ 2f543a99-f251-465a-99c0-d2682826fc5c
let
    for (i, j) in combinations(1:length(smiles))
        score = tdmces(smiles[i], smiles[j])
        println(f'{i} {j} {score}')

# ╔═╡ Cell order:
# ╠═21cea84e-126c-11ee-21b7-992b40f8d4ee
# ╠═3e067d65-930b-42b1-a399-7312660b4db1
# ╠═61294436-762a-4f02-bd9b-1ae8ecbd0c3c
# ╠═5ba7d2bd-bd30-4b97-a425-e1f67df75199
# ╠═2f543a99-f251-465a-99c0-d2682826fc5c
