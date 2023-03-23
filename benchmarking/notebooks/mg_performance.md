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
using Pkg
Pkg.activate(".")
using DelimitedFiles
using Plots
using StatsPlots
```

```julia
data = readdlm("./output/data.csv", ',')
data[:, 2:4] = data[:, 2:4] * 0.05  # 1 exec / mol
groupedbar([data[:, 2] data[:, 3] data[:, 4]], size =(600,600),
    bar_position = :dodge, bar_width=0.7, orientation = :h, fillrange=true,
    label=["No cache" "SSSR cached" "valence cached"], xlabel = "Time (ns)", ylabel = "Functions",
    xscale=:log10,  xlims= (1e2, 1e8), yticks = (1:24, data[:, 1]), yflip = true
)
```

```julia
data = readdlm("./output/data.csv", ',')
data[:, 2:4] = data[:, 2:4] * 0.05  # 1 exec / mol
groupedbar([data[:, 2] data[:, 3] data[:, 4]], size =(600,600),
    bar_position = :dodge, bar_width=0.7, orientation = :h, fillrange=true,
    label=["No cache" "SSSR cached" "valence cached"], xlabel = "Time (ns)", ylabel = "Functions",
    xlims= (0, 1e5), yticks = (1:24, data[:, 1]), yflip = true
)
```
