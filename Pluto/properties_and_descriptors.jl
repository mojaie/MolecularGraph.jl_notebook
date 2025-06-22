### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# ╔═╡ 10155c71-76c5-4e13-bce0-18cbdee039cd
using Graphs, MolecularGraph

# ╔═╡ 7df34918-d68f-11ed-0191-6b281f9c8a5f
md"
# Properties and descriptors

MolecularGraph.jl version: 0.20.0

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
oleic_acid = smilestomol(raw"CCCCCCCC\C=C/CCCCCCCC(O)=O")

# ╔═╡ f4004c63-168a-4a82-a4a6-2a75a14eb332
molecular_formula(oleic_acid)

# ╔═╡ 38a72c37-c180-4a2d-a19b-103dda1ce80e
empirical_formula(oleic_acid)

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
html_fixed_size(paraquat, 250, 250, atomindex=true)

# ╔═╡ 419e9f72-c528-4f57-afb2-ff079d5ac494
atom_symbol(paraquat)

# ╔═╡ bc4319e7-a247-4b46-bfee-baf8f1431d27
atom_charge(paraquat)

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
html_fixed_size(fluconazole, 250, 250, atomindex=true)

# ╔═╡ b53afd5e-d5b8-4f51-b52f-f07cb82ceb6e
bond_order(fluconazole)

# ╔═╡ 6a4c75f9-6441-4477-bb70-f1f431701166
is_edge_aromatic(fluconazole)

# ╔═╡ bc130afe-a9af-40ba-bc12-d64295a30350
is_rotatable(fluconazole)

# ╔═╡ c0e41581-7d3e-448d-99cc-6df8141b586a
is_rotatable(fluconazole)[MolecularGraph.edge_rank(fluconazole, 2, 9)]

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
html_fixed_size(stannsoporfin, 300, 300, atomindex=true)

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
html_fixed_size(mol, 250, 250, atomindex=true)

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
Graphs = "~1.13.0"
MolecularGraph = "~0.20.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.0-beta3"
manifest_format = "2.0"
project_hash = "2893ac0a6f9efcc7049fa731a60ca36cadf84c28"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "d57bd3762d308bded22c3b82d033bff85f6195c6"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.4.0"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.0+1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "4e1fe97fdaed23e9dc21d4d664bea76b65fc50a0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.22"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

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
git-tree-sha1 = "b309b36a9e02fe7be71270dd8c0fd873625332b4"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.6"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.GeoFormatTypes]]
git-tree-sha1 = "8e233d5167e63d708d41f87597433f59a0f213fe"
uuid = "68eda718-8dee-11e9-39e7-89f7f65f511f"
version = "0.4.4"

[[deps.GeoInterface]]
deps = ["DataAPI", "Extents", "GeoFormatTypes"]
git-tree-sha1 = "294e99f19869d0b0cb71aef92f19d03649d028d5"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.4.1"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "Extents", "GeoInterface", "IterTools", "LinearAlgebra", "PrecompileTools", "Random", "StaticArrays"]
git-tree-sha1 = "2670cf32dcf0229c9893b895a9afe725edb23545"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.5.9"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "c5abfa0ae0aaee162a3fbb053c13ecda39be545b"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.13.0"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.IntervalSets]]
git-tree-sha1 = "5fbb102dcb8b1a858111ae81d56682376130517d"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.11"

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

    [deps.IntervalSets.weakdeps]
    Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "a007feb38b422fbdab534406aeca1b86823cb4d6"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.11.1+1"

[[deps.LibGit2]]
deps = ["LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.9.0+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "OpenSSL_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.3+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.MakieCore]]
deps = ["ColorTypes", "GeometryBasics", "IntervalSets", "Observables"]
git-tree-sha1 = "6ec4aba99ce1a8c247822ce346b273d0d686e141"
uuid = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
version = "0.10.0"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MolecularGraph]]
deps = ["Base64", "Colors", "Dates", "DelimitedFiles", "GeometryBasics", "Graphs", "JSON", "LinearAlgebra", "MakieCore", "OrderedCollections", "Printf", "Statistics", "YAML", "coordgenlibs_jll", "libinchi_jll"]
git-tree-sha1 = "06c9b999d0476cfe3b6393272239556ac65ddf23"
uuid = "6c89ec66-9cd8-5372-9f91-fabc50dd27fd"
version = "0.20.0"

    [deps.MolecularGraph.extensions]
    CairoExt = "Cairo"
    RDKitExt = "RDKitMinimalLib"

    [deps.MolecularGraph.weakdeps]
    Cairo = "159f3aea-2a34-519c-b102-8c37f9878175"
    RDKitMinimalLib = "44044271-7623-48dc-8250-42433c44e4b7"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2025.2.25"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.0+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.12.0"

    [deps.Pkg.extensions]
    REPLExt = "REPL"

    [deps.Pkg.weakdeps]
    REPL = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "516f18f048a195409d6e072acf879a9f017d3900"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.2"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.12.0"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "0feb6b9031bd5c51f9072393eb5ab3efd31bf9e4"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.13"

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

    [deps.StaticArrays.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StringEncodings]]
deps = ["Libiconv_jll"]
git-tree-sha1 = "b765e46ba27ecf6b44faf70df40c57aa3a547dcb"
uuid = "69024149-9ee7-55f6-a4c4-859efe599b68"
version = "0.3.7"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.8.3+2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.YAML]]
deps = ["Base64", "Dates", "Printf", "StringEncodings"]
git-tree-sha1 = "2f58ac39f64b41fb812340347525be3b590cce3b"
uuid = "ddb6d928-2868-570f-bddf-ab3f9cf99eb6"
version = "0.4.14"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.coordgenlibs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "93fce743c1b36cde0efde11b1867cb7d16b13bf8"
uuid = "f6050b86-aaaf-512f-8549-0afff1b4d57f"
version = "3.0.2+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.12.0+0"

[[deps.libinchi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e81a593eb6a1a57d4262a3ed18a6efbc5d8ba83c"
uuid = "172afb32-8f1c-513b-968f-184fcd77af72"
version = "1.6.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.64.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.5.0+2"
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
