---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.2'
      jupytext_version: 1.4.0
  kernelspec:
    display_name: Julia 1.3.1
    language: julia
    name: julia-1.3
---

# Mass and isotopes

This tutorial includes following fundamental operations related to molecular mass and isotopes.

- Molecular weight (Relative molecular mass)
- Calculated exact mass and monoisotopic mass
- Isotopic composition
- Simulate mass spectrum


## Loading packages

This tutorial uses `PlotlyJS.jl` to display simulated mass spectrum.

```julia
import Pkg
Pkg.activate("..")
using MolecularGraph
using PlotlyJS
```

## Molecular weight (Relative molecular mass)

Molecular weight (Relative molecular mass) can be calculated by `standardweight` method.

```julia
standardweight(smilestomol("CCO"))
```

- `standardweight` returns the tuple of standard weight and its uncertainty, respectively. The result shown above means the molecular weight of C2H5O is 46.069±0.005.
- If you dont need uncertainty, give Float64 type to the method. Then the value will be rounded at 2 digits after the decimal point.

```julia
standardweight(Float64, smilestomol("CCO"))
```

## Calculated exact mass and monoisotopic mass

Both `exactmass` and `monoisotopicmass` returns exact mass of the given molecule. The difference is that `exactmass` considers atomic mass specified in `Atom.mass` fields of atoms whereas `monoisotopicmass` always returns exact mass of the molecule consists of the most abundant isotopes in terrestrial sources.

```julia
etohd6 = smilestomol("[2H]C([2H])([2H])C([2H])([2H])O[2H]")

println("Monoisotopic: ", monoisotopicmass(etohd6))
println("Exact: ", exactmass(etohd6))
implicithcount(etohd6)
```

```julia
mol = smilestomol("c1cc(ccc1[C@H]([C@@H](CO)NC(=O)C(Cl)Cl)O)[N+](=O)[O-]")
rcds = simulatemassspec(mol)
arr = zeros(Float64, length(rcds), 2)
for (i, rcd) in enumerate(rcds)
    arr[i, :] = [rcd[1], rcd[2]]
end
```

```julia
display("image/svg+xml", drawsvg(mol, 200, 200))

datapoints= scatter(
    x=arr[:, 1], y=arr[:, 2], name="Intensity",
    mode="markers", marker_size=2, marker_color="#2ca02c"
)
lines = map(eachrow(arr)) do r
    scatter(
        x=[r[1], r[1]], y=[0, r[2]],
        mode="lines", line_width=2, line_color="#2ca02c", hoverinfo="skip"
    )
end
data = Base.typed_vcat(GenericTrace, datapoints, lines)
layout = Layout(
    width=400, height=300,
    margin=Dict(:l => 50, :b => 50),
    showlegend=false,
    yaxis=Dict(:title => "Relative intensity"), xaxis=Dict(:title => "Mass")
)
plot(data, layout)
```
