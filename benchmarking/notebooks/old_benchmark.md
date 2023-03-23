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
import Pkg
Pkg.activate(".")
using Profile
using MolecularGraph

path =  joinpath("../../Workspace/resource/chemlib", "DrugBank_FDA_Approved.sdf")
# path =  joinpath(dirname(@__FILE__), "..", "_resources", "DrugBank_FDA_Approved.sdf")
```

```julia

paclitaxel = smilestomol(
        "CC1=C2[C@@]([C@]([C@H]([C@@H]3[C@]4([C@H](OC4)C[C@@H]([C@]3(C(=O)[C@@H]2OC(=O)C)C)O)OC(=O)C)OC(=O)c5ccccc5)(C[C@@H]1OC(=O)[C@H](O)[C@@H](NC(=O)c6ccccc6)c7ccccc7)O)(C)C"
    )
paclitaxel = removehydrogens(paclitaxel)
Profile.clear()
@profile functionalgroupgraph(paclitaxel)
Profile.print(mincount=400)
```

```julia
paclitaxel = smilestomol(
        "CC1=C2[C@@]([C@]([C@H]([C@@H]3[C@]4([C@H](OC4)C[C@@H]([C@]3(C(=O)[C@@H]2OC(=O)C)C)O)OC(=O)C)OC(=O)c5ccccc5)(C[C@@H]1OC(=O)[C@H](O)[C@@H](NC(=O)c6ccccc6)c7ccccc7)O)(C)C"
    )
paclitaxel = removehydrogens(paclitaxel)
Profile.clear()
@profile functionalgroupgraph(paclitaxel)
Profile.print(mincount=400)
```

```julia
Profile.clear()
@profile loadsdf()
Profile.print(mincount=400)
```

```julia
function substrsearch(smarts)
    query = parse(SMARTS, smarts)
    for (i, m) in enumerate(Iterators.take(sdfilereader(path), 200))
        isquerymatch(m, query)
    end
end
```

```julia
# Profile.clear()
@profile substrsearch("[OX2][OX2]")
Profile.print(mincount=400)
```
