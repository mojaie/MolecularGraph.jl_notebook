### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ e24f1a5d-40d3-497b-849f-066fb5c4fdbb
begin
    import Pkg
	Pkg.activate(Base.current_project())
	Pkg.instantiate()
	using MolecularGraph
	using Graphs
end

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
- `remove_hydrogens!(mol)` removes hyhdrogen vertices that are not important (no charge, no unpaired electron, no specific isotope composition and not involved in stereochemistry).

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

# ╔═╡ 74a9168a-30c2-41bc-80a1-8f44a2896880
pi_electron(quinoline1)

# ╔═╡ 3fda8505-fd5f-4f97-8f47-432465c6544b
pi_electron(quinoline2)

# ╔═╡ 76639ed0-f887-40e2-bcec-db176019ed9a
md"
- Property auto-update mechanisms described later assigns double bonds and single bonds to aromatic rings which consist of SMILES small letter atoms (so called 'kekulization'). We can see how kekulization works by the script as shown below.
"

# ╔═╡ 70d15699-ce3e-4cb3-9d06-c285d7f475b2
function custom_updater(mol)
	MolecularGraph.update_edge_rank!(mol)
    mol.state[:has_updates] = false
    stereocenter_from_smiles!(mol)
    stereobond_from_smiles!(mol)
	# skip kekulize!
end

# ╔═╡ 114518d1-8374-4ef7-9ec9-ea9326a9ab20
pyrene = smilestomol("c1cc2cccc3c2c4c1cccc4cc3", updater=custom_updater)

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
    mol.state[:has_updates] = false
    stereocenter_from_smiles!(mol)
    stereobond_from_smiles!(mol)
    kekulize!(mol)
    coordgen!(mol)
    # recalculate bottleneck descriptors
    sssr!(mol)
    lone_pair!(mol)
    apparent_valence!(mol)
    valence!(mol)
    is_ring_aromatic!(mol)
end
```

- The first two lines in the function are required. `update_edge_rank!` recalculate the edge rank cache if the graph topology was changed (e.g. by `rem_vertex!` method). Then, set the state `:has_updates` back to false.
- Next four lines specific to SMILES(stereochemistry resolvers, kekulization and 2D coords generation) also should be kept.
- The lines after that are to cache bottleneck descriptors. These can be removed if these descriptors would not be used.

"

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
# ╠═74a9168a-30c2-41bc-80a1-8f44a2896880
# ╠═3fda8505-fd5f-4f97-8f47-432465c6544b
# ╟─76639ed0-f887-40e2-bcec-db176019ed9a
# ╠═70d15699-ce3e-4cb3-9d06-c285d7f475b2
# ╠═114518d1-8374-4ef7-9ec9-ea9326a9ab20
# ╠═3b3b3da7-059b-4c20-aae6-70364d822559
# ╟─a1c2a97b-d6cd-4571-b50f-fb3e204eb470
