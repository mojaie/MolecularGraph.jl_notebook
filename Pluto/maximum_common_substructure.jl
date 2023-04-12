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
	remove_hydrogens!(mol)
	mol
end; nothing

# ╔═╡ 0c2ddc86-71ac-43d1-ae86-71510b216d92
aminocoumarin = let
	mol = sdftomol(fetch_mol!("74217", "3-Aminocoumarin", data_dir))
	remove_hydrogens!(mol)
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
	remove_hydrogens!(mol)
	mol
end; nothing

# ╔═╡ c2ca480e-da0a-45b2-9ea2-4567d6655769
ceftazidime = let
	mol = sdftomol(fetch_mol!("5481173", "Ceftazidime", data_dir))
	remove_hydrogens!(mol)
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
### Topological constraint (tcMCS)

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
