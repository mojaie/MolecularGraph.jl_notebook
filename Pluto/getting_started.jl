### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ afe097a1-6ae9-4c66-90f5-38707ae0c29c
begin
    import Pkg
	Pkg.activate(Base.current_project())
	Pkg.instantiate()
	using MolecularGraph
	using Graphs
end

# ╔═╡ 2f76d93c-221d-4c8c-a8ce-5f17b9f15cd1
md"
# Getting started

MolecularGraph.jl version: 0.14.0

This tutorial includes following fundamental operations of molecular mining.

- Chemical structure drawing
- Remove hydrogens
- Calculate molecular properties
- Substructure match and molecule query
- Maximum common substructure (MCS)
"

# ╔═╡ 51623203-9096-4125-a25b-1d66de219f8a
md"
#### Loading packages

Add `MolecularGraph` and `Graphs` packages to the workspace.
"

# ╔═╡ 632276aa-90c8-4e99-a4c0-4472dfa60afc
md"
#### Fetching SDFile from PubChem

Chemical structure data for tutorials can be downloaded from PubChem via HTTP. `_data` is created as a temporary data folder, and all test data in this tutorial will be stored in it."

# ╔═╡ 337bad91-f2c6-4d14-b8ae-df179fae46d5
data_dir = let
	# Create data directory
	data_dir = "_data"
	isdir(data_dir) || mkdir(data_dir)
	data_dir
end

# ╔═╡ 402d63c7-f1d4-47d8-9bd7-680f7fa0eef2
function fetch_mol!(cid, name, datadir)
	url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/$(cid)/SDF"
	dest = joinpath(data_dir, "$(name).mol")
	isfile(dest) || download(url, dest);
	return dest
end

# ╔═╡ bb14b0ff-7caa-4379-8524-56eefa04d24a
molfile = fetch_mol!("6437877", "Cefditoren Pivoxil", data_dir)

# ╔═╡ 5c9ce573-809f-47d4-901c-db744c054f76
molfile2 = fetch_mol!("5481173", "Ceftazidime", data_dir)

# ╔═╡ c15a6692-08ee-4d6d-94e6-09800b78ab36
md"
### Chemical structure drawing

#### From SDFIle

- Create molecular object from SDFile (.mol or .sdf) by using `sdftomol(filepath)`
- The default image size in pixel is set to (250, 250)
"

# ╔═╡ a368d692-f5d2-4392-b86d-ce431fac6b70
mol = sdftomol(molfile)

# ╔═╡ 3a762797-f60c-44cf-848b-4b8386aad162
md"
#### From SMILES

- Create molecular object from SMILES string by smilestomol(string)
- The default image size in pixel is set to (250, 250)
- The coordinates of atoms in SMILES will be automatically calculated. This library internally uses Schrodinger's coordgenlibs ([https://github.com/schrodinger/coordgenlibs](https://github.com/schrodinger/coordgenlibs)) for 2D coordinates generation.
"

# ╔═╡ b6252990-58a2-47c6-b8f3-bbd227e2ef86
smilestomol("O=C3N2/C(=C(/C=C\\c1scnc1C)CS[C@@H]2[C@@H]3NC(=O)C(=N\\OC)/c4nc(sc4)N)C(=O)O")

# ╔═╡ e0e12cc7-38dd-41b3-a92b-a167dd7af641
md"
### Remove hydrogens

- SDFiles downloaded from PubChem have hydrogen nodes. In practice, hydrogens which is not important are removed from molecular graphs for simplicity.
- `remove_hydrogens!(mol)` removes hyhdrogen vertices that are not important (no charge, no unpaired electron, no specific isotope composition and not involved in stereochemistry).
"

# ╔═╡ 090e782d-4892-480b-81f7-264c0f57327d
let
	remove_hydrogens!(mol)
	mol
end

# ╔═╡ 5b5c92af-afb5-4ef2-970d-dfdec0a1daa3
md"
- Even if these hydrogen 'nodes' are removed from the graph, we can infer the actual number of hydrogens from the atom valence (so called 'implicit hydrogens').
- Most of the molecule properties and descriptor calculations in the library consider such apparently removed hydrogens.
"

# ╔═╡ b4b46b10-dbf7-4c75-82ea-52d5ce05ea1c
nv(sdftomol(molfile))  # number of vertices in the graph

# ╔═╡ 190a645d-1117-4d39-a10c-3b2ce1c4f96c
standard_weight(sdftomol(molfile), 2)

# ╔═╡ e673c9e8-8183-46c6-a4ec-83666dcaf31b
nv(mol)

# ╔═╡ 86951ff0-ce9a-494c-bdd1-c12c474c9a90
standard_weight(mol, 2)

# ╔═╡ a24bc4d5-2120-4555-b6af-77bb376adfde
md"
### Calculate molecular properties

- As already shown in above, `standard_weight(mol[, digits])` returns the standard molecular weight.
- Other fundamental molecular property calculation methods are as folows
"

# ╔═╡ 15abe066-7b95-4878-93a3-385b46d8940d
hydrogen_acceptor_count(mol)  # F, O and N

# ╔═╡ 5214ffbd-2db1-4179-8553-b9b00d54fd55
hydrogen_donor_count(mol)  # O or N attached to at least one hydrogen

# ╔═╡ 4d271dd5-973f-46bc-a1d7-0342a9937e3d
rotatable_count(mol)  # rotatable bonds

# ╔═╡ 0dfb3c21-467e-42e0-aaa3-43b06d085f39
wclogp(mol, 3)  # wclogp(mol[, digits]), predicted logP by Wildman and Crippen method

# ╔═╡ c9dd2955-e873-410d-a7d6-d8ccef1bad40
molecular_formula(mol)  # formula in Hill system

# ╔═╡ 8f71c5b0-cd46-42be-9652-b455733b06dd
md"
### Substructure match and molecule query

- `has_substruct_match(mol, subst)` returns whether `subst` is the substructure of `mol`.
"

# ╔═╡ 9f021239-fda2-4cdf-9fd4-deba94d885c5
let
	substr = smilestomol("C1=CCSC2N1C(=O)C2")
	has_substruct_match(mol, substr)
end

# ╔═╡ 920ac17e-1195-4c6b-9e84-084591006c79
md"
- Similar to SMILES, SMARTS query can be parsed by using `smartstomol(string)` method
"

# ╔═╡ 4671fcc9-7293-42fc-b636-eb9059903bfa
let
	# non-aromatic ring nitrogen adjacent to double bond
	substr = smartstomol(raw"[$([N;r]A=A)]")
	has_substruct_match(mol, substr)
end

# ╔═╡ 5c9ad0b9-7d9a-4d12-83c9-aa0a867ba799
md"
- See 'Substructure match and molecule query' tutorial for more detailed substructure mathcing and visualization tasks.
"

# ╔═╡ 928c875b-e786-49b3-8395-7a5bcc81cdb6
let
	substr = smilestomol("s1cncc1")
	# matched vertex indices
	matched = vcat([collect(keys(m)) for m in substruct_matches(mol, substr)]...)
	# matched edges can be obtained from the node-induced subgraph
	matched_edges = induced_subgraph_edges(mol.graph, matched)
	# highlight matched vertices and edges 
	img = drawsvg(mol, atomhighlight=matched, bondhighlight=matched_edges)
	html_fixed_size(img, 300, 300)
end

# ╔═╡ 72f4aecb-edee-43fe-943c-c9eb2951d1ab
md"
### Maximum common substructure (MCS)

`MolecularGraph.jl` provides following fundamental MCS calculation methods for cheminformatics applications.

- Node-induced (MCIS) or edge-induced (MCES)
- Connected or disconnected
- Topological constraint (tdMCS)
"

# ╔═╡ 1e4102f4-e104-4acf-9409-260a37c2146a
cefditoren = let
	mol = sdftomol(molfile)
	remove_hydrogens!(mol)
	mol
end; nothing

# ╔═╡ 5aae8dc5-a22b-43bb-9439-352928a5a5d8
ceftazidime = let
	mol = sdftomol(molfile2)
	remove_hydrogens!(mol)
	mol
end; nothing

# ╔═╡ 980a93c5-8d8b-484d-9d63-02abaf53d41f
let
	# Connected maximum common edge-induced subgraph
	result = connected_mces(cefditoren, ceftazidime, tolerance=1)

	# atoms and bonds highlighting
	edges1 = collect(keys(result.mapping))
	nodes1, vmap1 = induced_subgraph(cefditoren.graph, edges1)
	svg1 = drawsvg(cefditoren,
		atomhighlight=vmap1[vertices(nodes1)], bondhighlight=edges1)

	edges2 = collect(values(result.mapping))
	nodes2, vmap2 = induced_subgraph(ceftazidime.graph, edges2)
	svg2 = drawsvg(ceftazidime,
		atomhighlight=vmap2[vertices(nodes2)], bondhighlight=edges2)

	html_grid([svg1, svg2], 3, 250)
end

# ╔═╡ a67e3a37-30a8-4ccf-8d1e-31a9d8cec359
let
	# Maximum common node-induced subgraph with topological constraint
	result = tcmcis(cefditoren, ceftazidime, tolerance=1)

	# atoms and bonds highlighting
	nodes1 = collect(keys(result.mapping))
	edges1 = induced_subgraph_edges(cefditoren.graph, nodes1)
	svg1 = drawsvg(cefditoren, atomhighlight=nodes1, bondhighlight=edges1)

	nodes2 = collect(values(result.mapping))
	edges2 = induced_subgraph_edges(ceftazidime.graph, nodes2)
	svg2 = drawsvg(ceftazidime, atomhighlight=nodes2, bondhighlight=edges2)

	html_grid([svg1, svg2], 3, 250)
end

# ╔═╡ Cell order:
# ╟─2f76d93c-221d-4c8c-a8ce-5f17b9f15cd1
# ╟─51623203-9096-4125-a25b-1d66de219f8a
# ╠═afe097a1-6ae9-4c66-90f5-38707ae0c29c
# ╟─632276aa-90c8-4e99-a4c0-4472dfa60afc
# ╠═337bad91-f2c6-4d14-b8ae-df179fae46d5
# ╠═402d63c7-f1d4-47d8-9bd7-680f7fa0eef2
# ╠═bb14b0ff-7caa-4379-8524-56eefa04d24a
# ╠═5c9ce573-809f-47d4-901c-db744c054f76
# ╟─c15a6692-08ee-4d6d-94e6-09800b78ab36
# ╠═a368d692-f5d2-4392-b86d-ce431fac6b70
# ╟─3a762797-f60c-44cf-848b-4b8386aad162
# ╠═b6252990-58a2-47c6-b8f3-bbd227e2ef86
# ╟─e0e12cc7-38dd-41b3-a92b-a167dd7af641
# ╠═090e782d-4892-480b-81f7-264c0f57327d
# ╟─5b5c92af-afb5-4ef2-970d-dfdec0a1daa3
# ╠═b4b46b10-dbf7-4c75-82ea-52d5ce05ea1c
# ╠═190a645d-1117-4d39-a10c-3b2ce1c4f96c
# ╠═e673c9e8-8183-46c6-a4ec-83666dcaf31b
# ╠═86951ff0-ce9a-494c-bdd1-c12c474c9a90
# ╟─a24bc4d5-2120-4555-b6af-77bb376adfde
# ╠═15abe066-7b95-4878-93a3-385b46d8940d
# ╠═5214ffbd-2db1-4179-8553-b9b00d54fd55
# ╠═4d271dd5-973f-46bc-a1d7-0342a9937e3d
# ╠═0dfb3c21-467e-42e0-aaa3-43b06d085f39
# ╠═c9dd2955-e873-410d-a7d6-d8ccef1bad40
# ╟─8f71c5b0-cd46-42be-9652-b455733b06dd
# ╠═9f021239-fda2-4cdf-9fd4-deba94d885c5
# ╟─920ac17e-1195-4c6b-9e84-084591006c79
# ╠═4671fcc9-7293-42fc-b636-eb9059903bfa
# ╟─5c9ad0b9-7d9a-4d12-83c9-aa0a867ba799
# ╠═928c875b-e786-49b3-8395-7a5bcc81cdb6
# ╟─72f4aecb-edee-43fe-943c-c9eb2951d1ab
# ╠═1e4102f4-e104-4acf-9409-260a37c2146a
# ╠═5aae8dc5-a22b-43bb-9439-352928a5a5d8
# ╠═980a93c5-8d8b-484d-9d63-02abaf53d41f
# ╠═a67e3a37-30a8-4ccf-8d1e-31a9d8cec359
