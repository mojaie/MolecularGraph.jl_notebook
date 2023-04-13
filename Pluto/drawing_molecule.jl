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

- Settings of 2D structure images
- Layout for web and Pluto notebook
- Regenerate 2D coordinates
- 3D molecule rendering using Makie.jl
"

# ╔═╡ d761d4ef-6d37-4382-ba4e-3dcbae4feced
data_dir = let
	# Create data directory
	data_dir = "_data"
	isdir(data_dir) || mkdir(data_dir)
	data_dir
end

# ╔═╡ 9d12e251-852d-431d-be4e-5a5314375dda
function fetch_mol!(cid, name, datadir)
	url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/$(cid)/SDF"
	dest = joinpath(data_dir, "$(name).mol")
	isfile(dest) || download(url, dest);
	return dest
end

# ╔═╡ 4feb7c09-20c6-4828-8e0e-cd819e1bacf6
molfile = fetch_mol!("6437877", "Cefditoren Pivoxil", data_dir)

# ╔═╡ c15a6692-08ee-4d6d-94e6-09800b78ab36
md"
### Settings of 2D structure images

#### Change image size

1. Generate a SVG image by `drawsvg(mol)`
1. Use `html_fixed_size(mol, width, height)` to display it
"

# ╔═╡ a368d692-f5d2-4392-b86d-ce431fac6b70
mol = let
	mol = sdftomol(joinpath(molfile))
	remove_hydrogens!(mol, all=false)
	html_fixed_size(drawsvg(mol), 400, 400)
end

# ╔═╡ eb399948-2ff1-47a3-8fc6-b44a440876e7
md"
### Layout for web and Pluto notebook
"

# ╔═╡ 3a762797-f60c-44cf-848b-4b8386aad162
md"
### Regenerate 2D coordinates

As SMILES does not have coordinates of atoms, so the 2D coordinates will be generated when smilestomol is called. Internally `MolecularGraph` uses Schrodinger's coordgenlibs ([https://github.com/schrodinger/coordgenlibs](https://github.com/schrodinger/coordgenlibs)) for 2D coords generation.
"

# ╔═╡ 69d8f6fa-d518-49a8-b802-466ce2750062
smol = smilestomol("O=C3N2/C(=C(/C=C\\c1scnc1C)CS[C@@H]2[C@@H]3NC(=O)C(=N\\OC)/c4nc(sc4)N)C(=O)O"); nothing

# ╔═╡ ec4cfaff-0cb2-4e9a-bf8f-a487b9fcc0ac
md"
- `coordgen!` also works with molecules from SDFile. This means you can recalculate 2D coordinate of that after molecular graph modification.
"

# ╔═╡ ab3d6fef-79f0-4458-a5b0-68f774733e55
html_fixed_size(drawsvg(smol, atomindex=true), 250, 250)

# ╔═╡ Cell order:
# ╟─2f76d93c-221d-4c8c-a8ce-5f17b9f15cd1
# ╠═afe097a1-6ae9-4c66-90f5-38707ae0c29c
# ╠═d761d4ef-6d37-4382-ba4e-3dcbae4feced
# ╠═9d12e251-852d-431d-be4e-5a5314375dda
# ╠═4feb7c09-20c6-4828-8e0e-cd819e1bacf6
# ╟─c15a6692-08ee-4d6d-94e6-09800b78ab36
# ╠═a368d692-f5d2-4392-b86d-ce431fac6b70
# ╠═eb399948-2ff1-47a3-8fc6-b44a440876e7
# ╠═3a762797-f60c-44cf-848b-4b8386aad162
# ╠═69d8f6fa-d518-49a8-b802-466ce2750062
# ╠═ab3d6fef-79f0-4458-a5b0-68f774733e55
# ╠═ec4cfaff-0cb2-4e9a-bf8f-a487b9fcc0ac
