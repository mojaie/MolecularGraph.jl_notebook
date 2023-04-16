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
  - Change image size
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

# ╔═╡ c15a6692-08ee-4d6d-94e6-09800b78ab36
md"
### Settings of 2D structure images

##### Change image size

1. Generate a SVG image by `drawsvg(mol)`
1. Use `html_fixed_size(mol, width, height)` to display it

- The default output of the molecule object in Pluto notebook is `html_fixed_size(drawsvg(mol), 250, 250)`, which displays the SVG image of the molecule with a width of 250px and a height of 250px.
"

# ╔═╡ a368d692-f5d2-4392-b86d-ce431fac6b70
mol = let
	molfile = fetch_mol!("6437877", "Cefditoren Pivoxil", data_dir)
	mol = sdftomol(molfile)
	remove_hydrogens!(mol, all=false)
	html_fixed_size(drawsvg(mol), 400, 400)
end

# ╔═╡ eb399948-2ff1-47a3-8fc6-b44a440876e7
md"
##### Layout for web and Pluto notebook

- The default output of the vector of molecule objects in Pluto notebook is `html_grid(drawsvg.(mols), 3, 250))`, which displays the SVG images of the molecules in a 3-column grid with a row height of 250px. The width is determined accorting to the HTML parent objects.
"

# ╔═╡ 384d0cd4-7d12-49f5-a432-2e472664ee55
mols = [
	smilestomol(raw"O=S(=O)(N)c1c(Cl)cc(c(C(=O)O)c1)NCc2occc2"),
	smilestomol(raw"CCC(CC)OC1C=C(CC(C1NC(=O)C)N)C(=O)OCC"),
	smilestomol(raw"O=C([C@H](CCC(O)=O)NC(C1=CC=C(N(CC2=CN=C(N=C(N)N=C3N)C3=N2)C)C=C1)=O)O"),
	smilestomol(raw"O=C2\N=C(/OC2c1ccccc1)N"),
	smilestomol(raw"O(CCN(C)C)C(c1ccccc1)c2ccccc2"),
	smilestomol(raw"C[C@@H]1C[C@H]2[C@@H](C[C@H]3[C@H](O2)[C@H]([C@@H]([C@H]4[C@H](O3)[C@H]([C@@H]([C@]5(O4)C[C@@H](CO5)O)C)C)O)C)O[C@H]6C[C@@H]7[C@]([C@@H](C[C@@H]8[C@@H](O7)C/C=C\C[C@@H]9[C@@H](O8)C=C[C@@H]2[C@@H](O9)C=C[C@@H]3[C@@H](O2)C[C@@H]2[C@@H](O3)[C@@H]([C@@H]3[C@@H](O2)CC=C[C@@H](O3)/C=C/[C@@H](CO)O)O)O)(O[C@@H]6C1)C")
]

# ╔═╡ b2cd3c97-7a5a-48eb-a412-83ed1c8080e3
html_grid(drawsvg.(mols), 4, 200)

# ╔═╡ 3a762797-f60c-44cf-848b-4b8386aad162
md"
### Regenerate 2D coordinates

As SMILES does not have coordinates of atoms, so the 2D coordinates will be generated when smilestomol is called. Internally `MolecularGraph` uses Schrodinger's coordgenlibs ([https://github.com/schrodinger/coordgenlibs](https://github.com/schrodinger/coordgenlibs)) for 2D coords generation.
"

# ╔═╡ 69d8f6fa-d518-49a8-b802-466ce2750062
smol = smilestomol(raw"O=C3N2/C(=C(/C=C\c1scnc1C)CS[C@@H]2[C@@H]3NC(=O)C(=N\OC)/c4nc(sc4)N)C(=O)O"); nothing

# ╔═╡ ab3d6fef-79f0-4458-a5b0-68f774733e55
html_fixed_size(drawsvg(smol, atomindex=true), 250, 250)

# ╔═╡ ec4cfaff-0cb2-4e9a-bf8f-a487b9fcc0ac
md"
- `coordgen!` also works with molecules from SDFile. This means you can recalculate 2D coordinate of that after molecular graph modification.
"

# ╔═╡ 4708cf8b-15b6-4769-92ec-a3b165b42f3b
md"
### 3D molecule rendering using Makie.jl

- To generate 3D molecular image, Makie.jl recipes `spacefilling` and `ballstick` are implemented.
- Makie recipes require one of the Makie backends to be installed. If GLMakie is installed, the code would be like the following.
- (Still installation and call of Makie backends takes very long time, and the appearance of 3D rendering is not good, help wanted)

```julia
using MolecularGraph
using GLMakie
mol = sdfmol(SDF_3D_FILE_PATH)
spacefilling(mol)
ballstick(mol)
```
"

# ╔═╡ Cell order:
# ╟─2f76d93c-221d-4c8c-a8ce-5f17b9f15cd1
# ╠═afe097a1-6ae9-4c66-90f5-38707ae0c29c
# ╠═d761d4ef-6d37-4382-ba4e-3dcbae4feced
# ╠═9d12e251-852d-431d-be4e-5a5314375dda
# ╟─c15a6692-08ee-4d6d-94e6-09800b78ab36
# ╠═a368d692-f5d2-4392-b86d-ce431fac6b70
# ╟─eb399948-2ff1-47a3-8fc6-b44a440876e7
# ╠═384d0cd4-7d12-49f5-a432-2e472664ee55
# ╠═b2cd3c97-7a5a-48eb-a412-83ed1c8080e3
# ╟─3a762797-f60c-44cf-848b-4b8386aad162
# ╠═69d8f6fa-d518-49a8-b802-466ce2750062
# ╠═ab3d6fef-79f0-4458-a5b0-68f774733e55
# ╟─ec4cfaff-0cb2-4e9a-bf8f-a487b9fcc0ac
# ╟─4708cf8b-15b6-4769-92ec-a3b165b42f3b
