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

- Standard molecular weight
- Relative molecular mass (Calculated exact mass)
- Monoisotopic mass
- Isotopic composition
- Simulate mass spectrum


## Loading packages

This tutorial uses `Plots.jl` to display simulated mass spectrum.

```julia
import Pkg
Pkg.activate("..")
using MolecularGraph
using Plots
gr()
Plots.GRBackend()

# Convenient function to display a pair of mol images
function displayimgpair(img1, img2)
    """<div>
        <div style="float:left">$(img1)</div>
        <div style="float:left">$(img2)</div>
    </div>"""
end
;
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

println("Molecular weight: ", standardweight(Float64, etohd6))
println("Monoisotopic: ", monoisotopicmass(Float64, etohd6))
println("Exact: ", exactmass(Float64, etohd6))
```

## Isotopic composition

`isotopiccomposition(atomsymbol, num)` returns total masses of the given number of atoms and their isotopic compositions.

Isotopes that have lower abundance than the `threshold` parameter will be filtered out (default 0.001 = 0.1%)

```julia
isotopiccomposition(:C, 100; threshold=0.01)
```

```julia
isotopiccomposition(:C, 1000; threshold=0.01)
```

`isotopiccomposition(mol)` returns isotopic composition of the molecule.

```julia
isotopiccomposition(smilestomol("CCO"))
```

## Simulate mass spectrum

According to the isotopic composition, simulated mass spectrum can be plotted by using `simulatemassspec` method and available plot libraries.
Follwing example uses Plots.jl to draw the simulated spectrum.

```julia
mol = smilestomol("c1cc(ccc1[C@H]([C@@H](CO)NC(=O)C(Cl)Cl)O)[N+](=O)[O-]")
structure = drawsvg(mol, 200, 200)
data = simulatemassspec(mol)
p = plot(data[:, 1], data[:, 2], size=(320, 240), leg=false, xlabel = "Mass", ylabel = "Intensity")
buf = IOBuffer()
show(buf, MIME("image/svg+xml"), p)
spectrum = String(take!(buf))
close(buf)
display("text/html", displayimgpair(structure, spectrum))
```
