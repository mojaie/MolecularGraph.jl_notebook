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



```julia
import Pkg
Pkg.activate("..")
using MolecularGraph
using Gadfly

# Convenient function to display a pair of mol images
function displayimgpair(img1, img2)
    """<div>
        <div style="float:left">$(img1)</div>
        <div style="float:left">$(img2)</div>
    </div>"""
end
```

```julia
mol = smilestomol("c1cc(ccc1[C@H]([C@@H](CO)NC(=O)C(Cl)Cl)O)[N+](=O)[O-]")
rcds = simulatemassspec(mol)
structure = drawsvg(mol, 200, 200)
buf = IOBuffer()
svg = SVG(buf, 8cm, 6cm)
svg.emit_on_finish = false
draw(svg, plot(
        x=rcds[:, 1], y=rcds[:, 2],
        Geom.hair, Guide.xlabel("Mass"), Guide.ylabel("Intensity")
))
spectrum = String(take!(buf))
display("text/html", displayimgpair(structure, spectrum))
```

```julia

```
