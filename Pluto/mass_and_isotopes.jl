### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 7fd1a4a3-7cf5-447f-a8c8-0513c90deda0
begin
    import Pkg
	Pkg.activate(Base.current_project())
	Pkg.instantiate()
	using MolecularGraph
	using Graphs
	using Plots
end

# ╔═╡ 690be766-d7e9-11ed-38b1-cd534cfdf86b
md"
# Mass and isotopes

MolecularGraph.jl version: 0.14.0

This tutorial includes following fundamental operations related to molecular mass and isotopes.

- Molecular weight and exact mass
- Uncertainty
- Isotopic composition
- Simulate mass spectrum

These are based on datasource of experimental isotopic composition provided by NIST [https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses](https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses)
and not intended to perform structure-based mass spec analysis (e.g. predict fragmentation pattern)
"

# ╔═╡ 8574b8bf-76f1-449f-95aa-7b9a1c84ab05
md"
### Molecular weight and exact mass

`standard_weight(mol[, digits])` returns the standard molecular weight (relative molecular mass).
"

# ╔═╡ f92dd039-4794-4337-a696-349d6dfae3a8
standard_weight(smilestomol("CCO"), 2)

# ╔═╡ 77eb309f-7130-436e-937f-a189eb914ecd
md"
Both `exact_mass(mol[, digits])` and `monoiso_mass(mol[, digits])` return exact mass of the given molecule. The difference is that `exact_mass` considers atomic mass specified in `mass` field of atom property type whereas `monoiso_mass` always returns exact mass of the molecule consists of the most abundant isotopes in terrestrial sources (monoisotopic mass).
"

# ╔═╡ 239a9f2a-8aef-4549-87f7-188a1fcacfbf
etohd6 = smilestomol("[2H]C([2H])([2H])C([2H])([2H])O[2H]")

# ╔═╡ 94320020-d5dd-4e2d-a365-916f9c8f6975
standard_weight(etohd6, 2)

# ╔═╡ 5de6d99e-1bec-40b8-be6b-a810af1a62b7
monoiso_mass(etohd6, 6)

# ╔═╡ dc467d05-c9ef-42ef-a187-0dee54c6e87c
exact_mass(etohd6, 6)

# ╔═╡ 9e90b521-9323-4e09-a4f3-81b97e186738
md"
`nominal_mass(mol)` is an aliase of `round(monoiso_mass(mol), digits=digits)`
"

# ╔═╡ 33fe8972-2567-4f49-9efb-45907b7a5629
nominal_mass(etohd6)

# ╔═╡ 920ff62f-cd3d-43cb-a6d7-bb6e85e0d9ca
md"
### Uncertainty

The above methods have `_unc` versions that returns a tuple of the mass value and its uncertainty. The following result means the molecular weight of C2H5O is 46.069±0.005.
"

# ╔═╡ 8e05f314-50b8-434f-9844-da96e0f4c4d7
standard_weight_unc(smilestomol("CCO"))

# ╔═╡ 3061517d-234a-491c-a9ee-27c2aec95578
md"
### Isotopic composition

- `isotopic_composition(atomsymbol, num)` returns total masses of the given number of atoms and their isotopic compositions.
- Isotopes that have lower abundance than the threshold parameter will be filtered out (default 0.001 = 0.1%)
- `isotopic_composition(mol)` returns isotopic composition of the molecule.
"

# ╔═╡ 2cb2e039-0003-48fe-b492-442f6e0aaa15
isotopic_composition(:C, 100; threshold=0.01)

# ╔═╡ fe73f4aa-6534-4049-b700-62bbb7351f1a
isotopic_composition(:C, 1000; threshold=0.01)

# ╔═╡ 6afef280-9279-49bd-9c63-610346bf36c2
isotopic_composition(smilestomol("CCO"))

# ╔═╡ 7acae275-5885-4136-95a1-41782f00616f
md"
### Simulate mass spectrum

According to the isotopic composition, simulated mass spectrum can be plotted by using `simulate_massspec` method and available plot libraries. Following example uses Plots.jl to draw the simulated spectrum.
"

# ╔═╡ e5760dc0-3c6d-44e4-9612-611f65de87e8
let
	mol = smilestomol("c1cc(ccc1[C@H]([C@@H](CO)NC(=O)C(Cl)Cl)O)[N+](=O)[O-]")
	structure = drawsvg(mol)
	data = simulate_massspec(mol)
	p = plot(data[:, 1], data[:, 2], size=(320, 240), leg=false, xlabel = "Mass", ylabel = "Intensity")
	buf = IOBuffer()
	show(buf, MIME("image/svg+xml"), p)
	spectrum = String(take!(buf))
	close(buf)
	html_grid([structure, spectrum], 2, 250)
end

# ╔═╡ Cell order:
# ╟─690be766-d7e9-11ed-38b1-cd534cfdf86b
# ╠═7fd1a4a3-7cf5-447f-a8c8-0513c90deda0
# ╟─8574b8bf-76f1-449f-95aa-7b9a1c84ab05
# ╠═f92dd039-4794-4337-a696-349d6dfae3a8
# ╟─77eb309f-7130-436e-937f-a189eb914ecd
# ╠═239a9f2a-8aef-4549-87f7-188a1fcacfbf
# ╠═94320020-d5dd-4e2d-a365-916f9c8f6975
# ╠═5de6d99e-1bec-40b8-be6b-a810af1a62b7
# ╠═dc467d05-c9ef-42ef-a187-0dee54c6e87c
# ╟─9e90b521-9323-4e09-a4f3-81b97e186738
# ╠═33fe8972-2567-4f49-9efb-45907b7a5629
# ╟─920ff62f-cd3d-43cb-a6d7-bb6e85e0d9ca
# ╠═8e05f314-50b8-434f-9844-da96e0f4c4d7
# ╟─3061517d-234a-491c-a9ee-27c2aec95578
# ╠═2cb2e039-0003-48fe-b492-442f6e0aaa15
# ╠═fe73f4aa-6534-4049-b700-62bbb7351f1a
# ╠═6afef280-9279-49bd-9c63-610346bf36c2
# ╟─7acae275-5885-4136-95a1-41782f00616f
# ╠═e5760dc0-3c6d-44e4-9612-611f65de87e8
