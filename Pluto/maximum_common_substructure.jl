### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 8b4b7b1a-8f77-4917-b55c-3f32ddd3bb4f
begin
    import Pkg
	Pkg.activate(Base.current_project())
	Pkg.instantiate()
	using Graphs, MolecularGraph
end

# ╔═╡ c324c158-d7fa-11ed-0300-93b6975d9deb
md"
# Maximum common substructure

MolecularGraph.jl version: 0.14.0

MolecularGraph.js implements following essential maximum common substructure (MCS) methods for cheminformatics applications.

- Maximum common induced substructure (MCIS)
- Maximum common edge-induced substructure (MCES)
- Connected or disconnected MCS
- Topological constraint (tdMCS)
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

# ╔═╡ 18e3d457-bd9e-418f-a08e-25eefcd02cd0
ldopa = let
	mol = sdftomol(fetch_mol!("6047", "Levodopa", data_dir))
	remove_hydrogens!(mol)
	mol
end; nothing

# ╔═╡ 0c2ddc86-71ac-43d1-ae86-71510b216d92
aminocoumarin = let
	mol = sdftomol(fetch_mol!("74217", "3-Aminocoumarin", data_dir))
	remove_hydrogens!(mol)
	mol
end; nothing

# ╔═╡ 65b40864-4cbc-4122-9b58-cfbecdee34ee
md"
### Maximum common induced substructure (MCIS)

- `disconnected_mcis(mol1, mol2)` calculates maximum common induced (node-induced) substructure (MCIS).
- 'induced' (or specifically 'node-induced') means that the subgraph is defined by a subset of the vertices of the original graph, and the presence or absence of edges between all vertices of the subgraph must match that of the original graph. That also means there is at least one vertex mapping (injection) from the subgraph to the original graph.

- MCS methods in this tutorial returns a `MCSResult` type object.
  - `MCSResult.mapping` is a mapping of vertex indices between two molecules.
  - `MCSResult.status` shows the calculation status. `:done` means that the calculation was successfully finished.
"

# ╔═╡ 8b645955-b2ce-4b2a-9793-18a1612e2a75
mcis_result = disconnected_mcis(ldopa, aminocoumarin)

# ╔═╡ 7b38cdb4-7afd-4957-b6d3-116b4bd90d96
let
	nodes1 = collect(keys(mcis_result.mapping))
	edges1 = induced_subgraph_edges(ldopa.graph, nodes1)
	svg1 = drawsvg(ldopa,
		bondhighlight=edges1, atomindex=true)

	nodes2 = collect(values(mcis_result.mapping))
	edges2 = induced_subgraph_edges(aminocoumarin.graph, nodes2)
	svg2 = drawsvg(aminocoumarin,
		atomhighlight=nodes2, bondhighlight=edges2, atomindex=true)

	html_grid([svg1, svg2], 3, 250)
end

# ╔═╡ 3c7be77e-6116-4c93-b3c0-a058cdd20a52
md"
### Connected MCS

- `connected_mcis(mol1, mol2)` also calculates MCIS but with connection constraint (obtained MCIS must be a connected component).
"

# ╔═╡ f9f50f92-d79e-42e2-a6cf-8107cd715f46
cefditoren = let
	mol = sdftomol(fetch_mol!("6437877", "Cefditoren Pivoxil", data_dir))
	remove_hydrogens!(mol)
	mol
end

# ╔═╡ c2ca480e-da0a-45b2-9ea2-4567d6655769
ceftazidime = let
	mol = sdftomol(fetch_mol!("5481173", "Ceftazidime", data_dir))
	remove_hydrogens!(mol)
	mol
end

# ╔═╡ fd1152c8-cbb4-4f92-aa5a-0d01d483df75


# ╔═╡ Cell order:
# ╟─c324c158-d7fa-11ed-0300-93b6975d9deb
# ╠═8b4b7b1a-8f77-4917-b55c-3f32ddd3bb4f
# ╠═62d3c057-1c65-44fd-ac25-266bdb26e536
# ╠═29a84294-8ab8-4d93-9d05-adcf67f6c57d
# ╠═18e3d457-bd9e-418f-a08e-25eefcd02cd0
# ╠═0c2ddc86-71ac-43d1-ae86-71510b216d92
# ╟─65b40864-4cbc-4122-9b58-cfbecdee34ee
# ╠═8b645955-b2ce-4b2a-9793-18a1612e2a75
# ╠═7b38cdb4-7afd-4957-b6d3-116b4bd90d96
# ╠═3c7be77e-6116-4c93-b3c0-a058cdd20a52
# ╠═f9f50f92-d79e-42e2-a6cf-8107cd715f46
# ╠═c2ca480e-da0a-45b2-9ea2-4567d6655769
# ╠═fd1152c8-cbb4-4f92-aa5a-0d01d483df75
