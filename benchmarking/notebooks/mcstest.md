---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.11.2
  kernelspec:
    display_name: Julia 1.6.0
    language: julia
    name: julia-1.6
---

```julia
using DotEnv
DotEnv.config("../.env")

using Pkg
Pkg.develop(PackageSpec(path=ENV["MOLGJ_DEV_PATH"]))

using Profile
using MolecularGraph
using MolecularGraph.Graph

```

```julia
minocycline = smilestomol("C[C@H]1[C@H]([C@H](C[C@@H](O1)O[C@H]2C[C@@](Cc3c2c(c4c(c3O)C(=O)c5cccc(c5C4=O)OC)O)(C(=O)CO)O)N)O")
doxorubicin = smilestomol("CN(C)c1ccc(c2c1C[C@H]3C[C@H]4[C@@H](C(=C(C(=O)[C@]4(C(=C3C2=O)O)O)C(=O)N)O)N(C)C)O")
precalculate!(minocycline)
precalculate!(doxorubicin)

Profile.clear()
@profile connectedmces(minocycline, doxorubicin)
Profile.print(mincount=20)
```

```julia

```
