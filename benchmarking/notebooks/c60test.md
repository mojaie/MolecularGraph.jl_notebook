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

using MolecularGraph
using MolecularGraph.Graph

```

```julia
c60 = smilestomol("c12c3c4c5c1c6c7c8c2c9c1c3c2c3c4c4c%10c5c5c6c6c7c7c%11c8c9c8c9c1c2c1c2c3c4c3c4c%10c5c5c6c6c7c7c%11c8c8c9c1c1c2c3c2c4c5c6c3c7c8c1c23")
precalculate!(c60)
molsvg = drawsvg(c60, 300, 300)
display("image/svg+xml",  molsvg)
```
