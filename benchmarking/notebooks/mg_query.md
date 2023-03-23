---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.11.4
  kernelspec:
    display_name: Julia 1.6.0
    language: julia
    name: julia-1.6
---

```julia
using DotEnv
DotEnv.config(joinpath(ENV["PWD"], ".env"))

using Pkg
Pkg.develop(PackageSpec(path=ENV["MOLGJ_DEV_PATH"]))

using Profile
using MolecularGraph
```

```julia
qr = query_relationship();
dox = smilestomol(raw"C[C@H]1[C@H]([C@H](C[C@@H](O1)O[C@H]2C[C@@](Cc3c2c(c4c(c3O)C(=O)c5cccc(c5C4=O)OC)O)(C(=O)CO)O)N)O");
doxh = removestereohydrogens(dox);
precalculate!(doxh);
precalculate!(dox);
qn = smartstomol(raw"[$([o,n]=c1ccc(=[o,n])cc1),$([O,N]=C1C=CC(=[O,N])C=C1),$([O,N]=C1[#6]:,=[#6]C(=[O,N])[#6]:,=[#6]1)]");


Profile.clear();
@profile filter_queries(qr, doxh, filtering=true);
Profile.print(mincount=20);
```

```julia
qr = query_relationship();
furo = smilestomol(raw"O=S(=O)(N)c1c(Cl)cc(c(C(=O)O)c1)NCc2occc2");
precalculate!(furo);

Profile.clear();
@profile filter_queries(qr, furo, filtering=true);
Profile.print(mincount=50);
```

```julia
@time filter_queries(qr, furo, filtering=true);
```

```julia
qr = query_relationship()
# mol = smilestomol(raw"O=S(=O)(N)c1c(Cl)cc(c(C(=O)O)c1)NCc2occc2")
mol = smilestomol(raw"C[C@H]1[C@H]([C@H](C[C@@H](O1)O[C@H]2C[C@@](Cc3c2c(c4c(c3O)C(=O)c5cccc(c5C4=O)OC)O)(C(=O)CO)O)N)O")

precalculate!(mol)
@time filter_queries(qr, mol, filtering=true);
```
