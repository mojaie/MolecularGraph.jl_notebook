### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 10155c71-76c5-4e13-bce0-18cbdee039cd
using Graphs, MolecularGraph

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

# ╔═╡ d99b6367-407f-416f-8291-fa5c5ed2fe63
paraquat = smilestomol("C[n+]1ccc(cc1)c2cc[n+](cc2)C.[Cl-].[Cl-]"); nothing

# ╔═╡ 3a7ebebe-317a-4a63-b90c-05d0dcb040b3
html_fixed_size(drawsvg(paraquat, atomindex=true), 250, 250)

# ╔═╡ 419e9f72-c528-4f57-afb2-ff079d5ac494
atom_symbol(paraquat)

# ╔═╡ bc4319e7-a247-4b46-bfee-baf8f1431d27
charge(paraquat)

# ╔═╡ 0ffeba60-457c-437b-8160-a5987ec384a6
multiplicity(paraquat)

# ╔═╡ 354eebcc-703c-45e9-9f12-8713fd0c3885
lone_pair(paraquat)

# ╔═╡ 8eef6837-08f2-4687-b030-c64d3cca3975
implicit_hydrogens(paraquat)

# ╔═╡ 336054f9-e727-40a2-93c6-d5156458433d
valence(paraquat)

# ╔═╡ b50eee82-7c09-4ea6-942d-7252906e21b4
is_aromatic(paraquat)

# ╔═╡ 52258140-e642-4d8a-9ab0-ebdf6898fdd8
pi_electron(paraquat)

# ╔═╡ eae0ed0f-9723-4304-bfe5-87efca6b7ec5
hybridization(paraquat)

# ╔═╡ 110a2830-af9a-4bb7-b60f-c7eaae1e3461
md"
The following methods return a vector of bond properties in lexicographic order of edges (called edge rank, the same order as `edges(mol)` method iterate). MolGraph caches the order and `edge_rank(mol, e)` or `edge_rank(mol, src(e), dst(e))` can be used to access the property vector.

- `bond_order(mol)`: only integer bond order is allowed in the model
- `is_edge_aromatic(mol)`: whether the bond is aromatic or not
- `is_rotatable(mol)`: whether the bond is rotatable or not
"

# ╔═╡ 4c1f78cc-7857-4064-b06a-19bffc9cb998
fluconazole = smilestomol("OC(Cn1cncn1)(Cn1cncn1)c1ccc(F)cc1F"); nothing

# ╔═╡ f7db2201-04a7-40a9-92a4-c41971c7888c
html_fixed_size(drawsvg(fluconazole, atomindex=true), 250, 250)

# ╔═╡ b53afd5e-d5b8-4f51-b52f-f07cb82ceb6e
bond_order(fluconazole)

# ╔═╡ 6a4c75f9-6441-4477-bb70-f1f431701166
is_edge_aromatic(fluconazole)

# ╔═╡ bc130afe-a9af-40ba-bc12-d64295a30350
is_rotatable(fluconazole)

# ╔═╡ c0e41581-7d3e-448d-99cc-6df8141b586a
is_rotatable(fluconazole)[edge_rank(fluconazole, 2, 9)]

# ╔═╡ 1223eb64-69d2-4d42-9fa8-43ca9123993d
md"
#### Graph topology

- `sssr(mol)`: smallest set of smallest rings
- `degree(mol)`: the degree of the vertex
- `connected_components(mol)`
- `fused_rings(mol)`: 2-edge connected components

"

# ╔═╡ 892b2441-0cd4-4640-886c-b9ab73381496
stannsoporfin = smilestomol(raw"[Cl-].[Cl-].[Sn+4].CCC1=C2[N-]C(C=C3N=C(C=C4[N-]C(=CC5=NC(=C2)C(C)=C5CC)C(C)=C4CCC(O)=O)C(CCC(O)=O)=C3C)=C1C"); nothing

# ╔═╡ ebc5f457-782f-41b2-b4a3-900c745c9937
html_fixed_size(drawsvg(stannsoporfin, atomindex=true), 300, 300)

# ╔═╡ 31286a95-9df9-4604-9aa1-692f6643cdb2
sssr(stannsoporfin)

# ╔═╡ 49194431-58a9-4579-ab99-2e737319fb5d
degree(stannsoporfin)

# ╔═╡ e5ec8f5d-6c76-4e67-8618-749f171772b7
connected_components(stannsoporfin)

# ╔═╡ 3c7c2152-5739-4e11-8106-983928b32599
fused_rings(stannsoporfin)

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
html_fixed_size(drawsvg(mol, atomindex=true), 250, 250)

# ╔═╡ 37916b47-b446-4658-bf28-9f91e540e39e
is_aromatic(mol)

# ╔═╡ 27b1df0e-59f4-4132-8105-13e1fa154a84
let
	mol = copy(mol)
	rem_vertex!(mol, 6)
	is_aromatic(mol)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Graphs = "86223c79-3864-5bf0-83f7-82e725a168b6"
MolecularGraph = "6c89ec66-9cd8-5372-9f91-fabc50dd27fd"

[compat]
Graphs = "~1.8.0"
MolecularGraph = "~0.14.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.2"
manifest_format = "2.0"
project_hash = "39cbfe9a3b79115b4c7c15351f8731e6f70a5385"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "cc37d689f599e8df4f464b2fa3870ff7db7492ef"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.6.1"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "7a60c856b9fa189eb34f5f8a6f6b5529b7942957"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.6.1"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.5.2+0"

[[deps.DataAPI]]
git-tree-sha1 = "e8119c1a33d267e16108be441a287a6981ba1630"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.14.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.Extents]]
git-tree-sha1 = "5e1e4c53fa39afe63a7d356e30452249365fba99"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.1"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "1cd7f0af1aa58abc02ea1d872953a97359cb87fa"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.4"

[[deps.GeoInterface]]
deps = ["Extents"]
git-tree-sha1 = "0eb6de0b312688f852f347171aba888658e29f20"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.3.0"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "GeoInterface", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "303202358e38d2b01ba46844b92e48a3c238fd9e"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.6"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "1cf1d7dcb4bc32d7b4a5add4232db3750c27ecb4"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.8.0"

[[deps.Inflate]]
git-tree-sha1 = "5cd07aab533df5170988219191dfad0519391428"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.3"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LazyJSON]]
deps = ["JSON", "OrderedCollections", "PropertyDicts"]
git-tree-sha1 = "ce08411caa70e0c9e780f142f59debd89a971738"
uuid = "fc18253b-5e1b-504c-a4a2-9ece4944c004"
version = "0.2.2"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c7cb1f5d892775ba13767a87c7ada0b980ea0a71"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+2"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.MakieCore]]
deps = ["Observables"]
git-tree-sha1 = "9926529455a331ed73c19ff06d16906737a876ed"
uuid = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
version = "0.6.3"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MolecularGraph]]
deps = ["Colors", "Dates", "DelimitedFiles", "GeometryBasics", "Graphs", "JSON", "LinearAlgebra", "MakieCore", "Printf", "Statistics", "Unmarshal", "YAML", "coordgenlibs_jll", "libinchi_jll"]
git-tree-sha1 = "18fa7feaca1d4af33eed1e044357327bdba7f1ce"
uuid = "6c89ec66-9cd8-5372-9f91-fabc50dd27fd"
version = "0.14.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Nullables]]
git-tree-sha1 = "8f87854cc8f3685a60689d8edecaa29d2251979b"
uuid = "4d1e1d77-625e-5b40-9113-a560ec7a8ecd"
version = "1.0.0"

[[deps.Observables]]
git-tree-sha1 = "6862738f9796b3edc1c09d0890afce4eca9e7e93"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.4"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "d321bf2de576bf25ec4d3e4360faca399afca282"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.0"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "478ac6c952fddd4399e71d4779797c538d0ff2bf"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.8"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.PropertyDicts]]
git-tree-sha1 = "8cf3b5cea994cfa9f238e19c3946a39cf051896c"
uuid = "f8a19df8-e894-5f55-a973-672c1158cbca"
version = "0.1.2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "63e84b7fdf5021026d0f17f76af7c57772313d99"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.21"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StringEncodings]]
deps = ["Libiconv_jll"]
git-tree-sha1 = "33c0da881af3248dafefb939a21694b97cfece76"
uuid = "69024149-9ee7-55f6-a4c4-859efe599b68"
version = "0.3.6"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "GPUArraysCore", "StaticArraysCore", "Tables"]
git-tree-sha1 = "521a0e828e98bb69042fec1809c1b5a680eb7389"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.15"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "1544b926975372da01227b382066ab70e574a3ec"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Unmarshal]]
deps = ["JSON", "LazyJSON", "Missings", "Nullables", "Requires"]
git-tree-sha1 = "ee46863309f8f942249e1df1b74ba3088ff0f151"
uuid = "cbff2730-442d-58d7-89d1-8e530c41eb02"
version = "0.4.4"

[[deps.YAML]]
deps = ["Base64", "Dates", "Printf", "StringEncodings"]
git-tree-sha1 = "dbc7f1c0012a69486af79c8bcdb31be820670ba2"
uuid = "ddb6d928-2868-570f-bddf-ab3f9cf99eb6"
version = "0.4.8"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.coordgenlibs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "8a0fdb746dfc75758d0abea3196f5edfcbbebd79"
uuid = "f6050b86-aaaf-512f-8549-0afff1b4d57f"
version = "3.0.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.libinchi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "034ee07d3b387a4ca1a153a43a0c46549b6749ba"
uuid = "172afb32-8f1c-513b-968f-184fcd77af72"
version = "1.5.1+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╟─7df34918-d68f-11ed-0191-6b281f9c8a5f
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
# ╠═d99b6367-407f-416f-8291-fa5c5ed2fe63
# ╠═3a7ebebe-317a-4a63-b90c-05d0dcb040b3
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
# ╠═4c1f78cc-7857-4064-b06a-19bffc9cb998
# ╠═f7db2201-04a7-40a9-92a4-c41971c7888c
# ╠═b53afd5e-d5b8-4f51-b52f-f07cb82ceb6e
# ╠═6a4c75f9-6441-4477-bb70-f1f431701166
# ╠═bc130afe-a9af-40ba-bc12-d64295a30350
# ╠═c0e41581-7d3e-448d-99cc-6df8141b586a
# ╟─1223eb64-69d2-4d42-9fa8-43ca9123993d
# ╠═892b2441-0cd4-4640-886c-b9ab73381496
# ╠═ebc5f457-782f-41b2-b4a3-900c745c9937
# ╠═31286a95-9df9-4604-9aa1-692f6643cdb2
# ╠═49194431-58a9-4579-ab99-2e737319fb5d
# ╠═e5ec8f5d-6c76-4e67-8618-749f171772b7
# ╠═3c7c2152-5739-4e11-8106-983928b32599
# ╟─a14b158a-5cee-44a1-bedc-4fbce52a67a3
# ╠═684c7188-6e78-44d1-819c-acaf599b229b
# ╠═37916b47-b446-4658-bf28-9f91e540e39e
# ╠═27b1df0e-59f4-4132-8105-13e1fa154a84
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
