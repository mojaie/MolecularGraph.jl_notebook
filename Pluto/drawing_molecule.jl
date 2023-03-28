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
# Drawing molecules

MolecularGraph.jl version: 0.14.0

This tutorial deals with 2D/3D rendering of molecules.

- Image output settings
- Molecule generation
- 3D molecule rendering using Makie.jl
"

# ╔═╡ 51623203-9096-4125-a25b-1d66de219f8a
md"
#### Loading packages

Add `MolecularGraph` and `Graphs` packages to the workspace.
"

# ╔═╡ 632276aa-90c8-4e99-a4c0-4472dfa60afc
md"
#### Fetching test molecule data from public resources

Chemical structure data for tutorials can be downloaded from PubChem via HTTP. `_data` is created as a temporary data folder, and all test data in this tutorial will be stored in it."

# ╔═╡ bb14b0ff-7caa-4379-8524-56eefa04d24a
molfile = let
	cid = "6437877"  # PubChem CID
	name = "Cefditoren Pivoxil"  # FIle name

	# Create data directory
	data_dir = "_data"
	isdir(data_dir) || mkdir(data_dir)

	# Fetch
	url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/$(cid)/SDF"
	dest = joinpath(data_dir, "$(name).mol")
	isfile(dest) || download(url, dest);
	dest
end

# ╔═╡ c15a6692-08ee-4d6d-94e6-09800b78ab36
md"
### Structure drawing options

#### Change image size

1. Generate a SVG image by `drawsvg(mol, width, height)`
1. Use `HTML` to display it
"

# ╔═╡ a368d692-f5d2-4392-b86d-ce431fac6b70
mol = let
	mol = sdftomol(joinpath(molfile))
	remove_hydrogens!(mol, all=false)
	HTML(drawsvg(mol, 400, 400))
end

# ╔═╡ 3a762797-f60c-44cf-848b-4b8386aad162
md"
#### Generate 2D coordinate of atoms

As SMILES does not have coordinates of atoms, so the 2D coordinates will be generated when smilestomol is called. Internally `MolecularGraph` uses Schrodinger's coordgenlibs ([https://github.com/schrodinger/coordgenlibs](https://github.com/schrodinger/coordgenlibs)) for 2D coords generation.
"

# ╔═╡ 69d8f6fa-d518-49a8-b802-466ce2750062
smol = let 
	mol = smilestomol("O=C3N2/C(=C(/C=C\\c1scnc1C)CS[C@@H]2[C@@H]3NC(=O)C(=N\\OC)/c4nc(sc4)N)C(=O)O")
	coordgen!(mol)
	mol
end

# ╔═╡ ec4cfaff-0cb2-4e9a-bf8f-a487b9fcc0ac
md"
- `coordgen!` also works with molecules from SDFile. This means you can recalculate 2D coordinate of that after molecular graph modification.
"

# ╔═╡ ab3d6fef-79f0-4458-a5b0-68f774733e55
let
	canvas = SvgCanvas()
	draw2d!(canvas, smol)
	drawatomindex!(canvas, smol)
	molsvg = tosvg(canvas, 250, 250)
	HTML(molsvg)
end

# ╔═╡ Cell order:
# ╟─2f76d93c-221d-4c8c-a8ce-5f17b9f15cd1
# ╟─51623203-9096-4125-a25b-1d66de219f8a
# ╠═afe097a1-6ae9-4c66-90f5-38707ae0c29c
# ╟─632276aa-90c8-4e99-a4c0-4472dfa60afc
# ╠═bb14b0ff-7caa-4379-8524-56eefa04d24a
# ╟─c15a6692-08ee-4d6d-94e6-09800b78ab36
# ╠═a368d692-f5d2-4392-b86d-ce431fac6b70
# ╠═3a762797-f60c-44cf-848b-4b8386aad162
# ╠═69d8f6fa-d518-49a8-b802-466ce2750062
# ╠═ec4cfaff-0cb2-4e9a-bf8f-a487b9fcc0ac
# ╠═ab3d6fef-79f0-4458-a5b0-68f774733e55
