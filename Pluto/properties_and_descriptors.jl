### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 10155c71-76c5-4e13-bce0-18cbdee039cd
begin
    import Pkg
	Pkg.activate(Base.current_project())
	Pkg.instantiate()
	using MolecularGraph
	using Graphs
end

# ╔═╡ 7df34918-d68f-11ed-0191-6b281f9c8a5f
md"
# Properties and descriptors

MolecularGraph.jl version: 0.14.0

This tutorial includes following concepts of molecular graph properties.

- Built-in molecule properties and descriptors
  - Lipinski's Rule of five (RO5)
  - Molecular formula
  - Atom and bond properties
  - Graph topology (ring and fused ring)
- Auto-update mechanism of properties
"

# ╔═╡ b1d563c7-a3ee-4563-93cb-cf1cfd4fca14
md"
### Built-in molecule properties

Here are examples of built-in properties often used in molecular mining.

"

# ╔═╡ f6e1c81f-04f4-4250-b980-b2a730586dbf
md"
##### Lipinski's Rule of five (RO5)

`digits` (after the decimal point) should be specfied in some methods return Float values

- `standard_weight(mol[, digits])`: molecular weight
- `hydrogen_acceptor_count(mol)`: number of hydrogen acceptor
- `hydrogen_donor_count(mol)`: number of hydrogen donor
- `wclogp(mol[, digits])`: predicted logP by Wildman and Clippen (1999), https://doi.org/10.1021/ci990307l
- `rotatable_count(mol)`: number of rotatable bonds (in some RO5 variants)
"

# ╔═╡ da8f8753-db36-48f3-89d4-63b46d50bc4f
mol = smilestomol("CC(=O)OC1=CC=CC=C1C(=O)O")

# ╔═╡ 7dd3cda6-341d-4dd3-8419-8ab44c739fa0
standard_weight(mol, 2)

# ╔═╡ 76ec3420-36ac-4782-a0d7-4f982312cd2f
hydrogen_acceptor_count(mol)

# ╔═╡ 39bf3f2e-d9af-4282-91bc-61c0cd16d183
hydrogen_donor_count(mol)

# ╔═╡ 0a5692bf-59a6-4edd-a6b9-b1b5ca31c19e
wclogp(mol, 2)

# ╔═╡ 11c57c9d-bf4b-4f40-bb5a-1d3357b69d0f
rotatable_count(mol)

# ╔═╡ e5281cc4-1e41-4697-b76c-4ce429284ba3
md"
##### Molecular formula

- `molecular_formula(mol)`: molecular formula in Hill system
- `empirical_formula(mol)`: empirical formula in Hill system

"

# ╔═╡ d2d96e2c-1a50-453c-9629-3d061e1961fa
glc = smilestomol("OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O")

# ╔═╡ f4004c63-168a-4a82-a4a6-2a75a14eb332
molecular_formula(glc)

# ╔═╡ 38a72c37-c180-4a2d-a19b-103dda1ce80e
empirical_formula(glc)

# ╔═╡ 6261f0a7-997b-4444-8286-b59c0325e6bc
md"
##### Atom and bond properties

The following methods return a vector of atom properties in index order

- `atom_symbol(mol)`: atom letters as a symbol e.g. :C, :O and :N
- `charge(mol)`: electric charge of the atom. only integer charge is allowed in the model
- `multiplicity(mol)`: 1: no unpaired electron(default), 2: radical, 3: biradical

- `lone_pair(mol)`: number of lone pair on the atom
- `implicit_hydrogens(mol)`: number of implicit hydrogens that are not appear as graph vertices but automatically calculated, drawn in image and used for calculation of other descriptors.
- `valence(mol)`: number of atom valence, specific to each atom species and considering electric charge. Implicit number of hydrogens is obtained by subtracting the degree of the vertex from the valence.
- `is_aromatic(mol)`: whether the atom is aromatic or not. only binary aromaticity is allowed in the model.
- `pi_electron(mol)`: number of pi electrons
- `hybridization(mol)`: orbital hybridization e.g. sp, sp2 and sp3
"

# ╔═╡ 419e9f72-c528-4f57-afb2-ff079d5ac494
atom_symbol(mol)

# ╔═╡ bc4319e7-a247-4b46-bfee-baf8f1431d27
charge(mol)

# ╔═╡ 0ffeba60-457c-437b-8160-a5987ec384a6
multiplicity(mol)

# ╔═╡ 354eebcc-703c-45e9-9f12-8713fd0c3885
lone_pair(mol)

# ╔═╡ 8eef6837-08f2-4687-b030-c64d3cca3975
implicit_hydrogens(mol)

# ╔═╡ 336054f9-e727-40a2-93c6-d5156458433d
valence(mol)

# ╔═╡ b50eee82-7c09-4ea6-942d-7252906e21b4
is_aromatic(mol)

# ╔═╡ 52258140-e642-4d8a-9ab0-ebdf6898fdd8
pi_electron(mol)

# ╔═╡ eae0ed0f-9723-4304-bfe5-87efca6b7ec5
hybridization(mol)

# ╔═╡ 110a2830-af9a-4bb7-b60f-c7eaae1e3461
md"
The following methods return a vector of bond properties in lexicographic order of edges (called edge rank, the same order as `edges(mol)` method iterate). MolGraph caches the order and `edge_rank(mol, e)` or `edge_rank(mol, src(e), dst(e))` can be used to access the property vector.

- `bond_order(mol)`: only integer bond order is allowed in the model
- `is_edge_aromatic(mol)`: whether the bond is aromatic or not
- `is_rotatable(mol)`: whether the bond is rotatable or not
"

# ╔═╡ b53afd5e-d5b8-4f51-b52f-f07cb82ceb6e
bond_order(mol)

# ╔═╡ 6a4c75f9-6441-4477-bb70-f1f431701166
is_edge_aromatic(mol)

# ╔═╡ bc130afe-a9af-40ba-bc12-d64295a30350
is_rotatable(mol)

# ╔═╡ 1223eb64-69d2-4d42-9fa8-43ca9123993d
md"
#### Graph topology

- `sssr(mol)`: smallest set of smallest rings
- `degree(mol)`: the degree of the vertex
- `connected_components(mol)`
- `fused_rings(mol)`: 2-edge connected components

"

# ╔═╡ 31286a95-9df9-4604-9aa1-692f6643cdb2
sssr(mol)

# ╔═╡ 49194431-58a9-4579-ab99-2e737319fb5d
degree(mol)

# ╔═╡ e5ec8f5d-6c76-4e67-8618-749f171772b7
connected_components(mol)

# ╔═╡ 3c7c2152-5739-4e11-8106-983928b32599
fused_rings(mol)

# ╔═╡ a14b158a-5cee-44a1-bedc-4fbce52a67a3
md"
# Auto-update mechanism of properties

The following are examples of properties often considered and implemented in molecular graph models

- atom (vertex) properties
  - atom symbol e.g. C, O, N, ...
  - atomic charge
  - whether aromatic or not
- bond (edge) properties
  - bond order
  - whether aromatic or not
- graph topology
  - connected components (a molecule object with multiple molecules)
  - smallest set of smallest rings (SSSR)
- graph properties
  - metadata e.g. name, compound ID

What is important is that these properties can have interactions. In other words, many of the vertex/edge properties depend on other properties in the vertex/edge or adjacent vertexes/edges. This also means that removal or addition of vertices/edges or changes in graph topology can affect individual properties.

To reduce this kind of side-effects, keep consistency and obtain reproducible results, molecule properties in this library is implemented based on following categories.

- primary property, or simply 'property': Information on atom symbol, charges, bond orders, etc. that define the molecule obtained from SMILES, SDFile or databases.
- secondary property, or 'descriptor': Properties that depend on primary properties and graph topology. Methods to generate secondary properties are required to be 'pure function'. That means these methods should take only graph object (SimpleGraph) and primary property vector as arguments and should not alter objects outside the method.

(Note that these terms are not common technical terms. Just for implementation convenience, they are categorized as above.)

Primary properties are stored in atom/bond property objects (e.g. SMILESAtom). On the other hand, Secondary properties are calculated ad-hoc, or called from caches (stored in MolGraph.state field) if it is considered to be expensive and frequently called.

If primary properties or graph topology were changed (e.g. remove vertices), the molecule object is marked as `:has_updates`, and recalculated with a pre-defined routine (see preprocessing tutorial for details) when the next time the secondary property is called.
"

# ╔═╡ 684c7188-6e78-44d1-819c-acaf599b229b
HTML(drawsvg(mol, 250, 250, atomindex=true))

# ╔═╡ 37916b47-b446-4658-bf28-9f91e540e39e
is_aromatic(mol)

# ╔═╡ 27b1df0e-59f4-4132-8105-13e1fa154a84
let
	mol = copy(mol)
	rem_vertex!(mol, 6)
	is_aromatic(mol)
end

# ╔═╡ Cell order:
# ╠═7df34918-d68f-11ed-0191-6b281f9c8a5f
# ╠═10155c71-76c5-4e13-bce0-18cbdee039cd
# ╟─b1d563c7-a3ee-4563-93cb-cf1cfd4fca14
# ╟─f6e1c81f-04f4-4250-b980-b2a730586dbf
# ╠═da8f8753-db36-48f3-89d4-63b46d50bc4f
# ╠═7dd3cda6-341d-4dd3-8419-8ab44c739fa0
# ╠═76ec3420-36ac-4782-a0d7-4f982312cd2f
# ╠═39bf3f2e-d9af-4282-91bc-61c0cd16d183
# ╠═0a5692bf-59a6-4edd-a6b9-b1b5ca31c19e
# ╠═11c57c9d-bf4b-4f40-bb5a-1d3357b69d0f
# ╟─e5281cc4-1e41-4697-b76c-4ce429284ba3
# ╠═d2d96e2c-1a50-453c-9629-3d061e1961fa
# ╠═f4004c63-168a-4a82-a4a6-2a75a14eb332
# ╠═38a72c37-c180-4a2d-a19b-103dda1ce80e
# ╟─6261f0a7-997b-4444-8286-b59c0325e6bc
# ╠═419e9f72-c528-4f57-afb2-ff079d5ac494
# ╠═bc4319e7-a247-4b46-bfee-baf8f1431d27
# ╠═0ffeba60-457c-437b-8160-a5987ec384a6
# ╠═354eebcc-703c-45e9-9f12-8713fd0c3885
# ╠═8eef6837-08f2-4687-b030-c64d3cca3975
# ╠═336054f9-e727-40a2-93c6-d5156458433d
# ╠═b50eee82-7c09-4ea6-942d-7252906e21b4
# ╠═52258140-e642-4d8a-9ab0-ebdf6898fdd8
# ╠═eae0ed0f-9723-4304-bfe5-87efca6b7ec5
# ╟─110a2830-af9a-4bb7-b60f-c7eaae1e3461
# ╠═b53afd5e-d5b8-4f51-b52f-f07cb82ceb6e
# ╠═6a4c75f9-6441-4477-bb70-f1f431701166
# ╠═bc130afe-a9af-40ba-bc12-d64295a30350
# ╠═1223eb64-69d2-4d42-9fa8-43ca9123993d
# ╠═31286a95-9df9-4604-9aa1-692f6643cdb2
# ╠═49194431-58a9-4579-ab99-2e737319fb5d
# ╠═e5ec8f5d-6c76-4e67-8618-749f171772b7
# ╠═3c7c2152-5739-4e11-8106-983928b32599
# ╟─a14b158a-5cee-44a1-bedc-4fbce52a67a3
# ╠═684c7188-6e78-44d1-819c-acaf599b229b
# ╠═37916b47-b446-4658-bf28-9f91e540e39e
# ╠═27b1df0e-59f4-4132-8105-13e1fa154a84
