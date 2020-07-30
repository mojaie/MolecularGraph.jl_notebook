---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.2'
      jupytext_version: 1.4.2
  kernelspec:
    display_name: Julia 1.4.1
    language: julia
    name: julia-1.4
---

<!-- #region -->
# Calculate descriptors


This tutorial includes following fundamental operations for molecular descriptors.

- Concept of descriptor function/array
- Cache mechanism
- Frequently used descriptors
<!-- #endregion -->

```julia
using Pkg
Pkg.activate("..")
using MolecularGraph
```

```julia
# Demo molecule with atom indices
mol = smilestomol("Clc4cc2c(C(/c1ncccc1CC2)=C3/CCNCC3)cc4")
canvas = SvgCanvas()
draw2d!(canvas, mol)
drawatomindex!(canvas, mol)
molsvg = tosvg(canvas, 300, 300)
display("image/svg+xml",  molsvg)
```

## Concept of descriptor function/array

Descriptor array is typically a vector calculated from a molecular object by a descriptor function. Most of fundamental descriptor functions are coded in `src/properties.jl`.

```julia
println(hydrogenconnected(mol))
```

<!-- #region -->
`hydrogenconnected(mol)` returns a vector of the total number of hydrogens connected to each atom nodes. For example, atom \#13 is a secondary carbon so we can infer there is two hydrogens on it, therefore hydrogenconnected(mol)\[13\] is 2.

`hydrogenconnected` is defined as below in `src/properties.jl`

```julia
@cachefirst hydrogenconnected(mol::GraphMol
    ) = explicithconnected(mol) + implicithconnected(mol)
```

`@cachefirst` macro is to cache the results for fast calculation (described later). We can know that `hydrogenconnected` is derived from `explicithconnected` and `implicithconnected`.

`implicithconnected` is defined as below in `src/properties.jl`

```julia
@cachefirst function implicithconnected(mol::GraphMol)
    hcnt = (v, av) -> v === nothing ? 0 : max(0, v - av)
    return hcnt.(valence(mol), apparentvalence(mol))
end
```

so `implicithconnected` depends on `valence` and `apparentvalence`.

We can trace back the dependencies, and then can find these descriptors are derived from graph topologies (e.g. neighbors) and attributes (e.g. Atom.symbol and Bond.order). This descriptor function/array system can also deal with more complicated descriptors such as `isrotatable` and `isaromatic` as like a simple function that takes molecule object and returns an array without any adverse effects.

Note that manipulation of the molecular graph topology and attributes will break consistency in descriptor arrays. It is important to define molecular preprocessing workflow first, and then calculate descriptors.
<!-- #endregion -->

## Cache mechanism

`Graph.@cachefirst` macro enables caching descriptor arrays already calculated. This macro is set to most of default descriptor functions. We can simply call functions with `Graph.@cache` macro to return the result and also store the result to GraphMol.cache field.


```julia
using MolecularGraph.Graph

mol2 = clone(mol)
v = @cache valence(mol2)
println(v)
```

```julia
mol2.cache
```

```julia
using BenchmarkTools

@benchmark valence(mol)
```

```julia
@benchmark valence(mol2)
```

### Tips

- Some descriptor function allows keyword arguments (e.g. `coordgen(mol; forcecoordgen=true)`). If keyword arguments are given explicitly, cache will not be used.
- If the molecule object is modified (e.g. by preprocessing methods), you can recalculate and cache the result by calling the function with `@cache` macro again.
- `clearcache!(mol)` empties the cache dict.
- `precalculate!(mol)` is a convenient method to store caches of cost-effective 5 or 6 descriptors. It is recommended to apply this method after you defined your molecular model and run preprocessing methods.

```julia
clearcache!(mol2)
mol2.cache
```

```julia
precalculate!(mol2)
mol2.cache
```

## Convention of descriptor functions

For reproducibility, there is some important conventions for descriptor functions.

- All descriptor functions take a GraphMol or SubgraphView as the only argument (optional keyword arguments are acceptable but caching will not work).
- All descriptor functions never change the molecule object given as an argument.
- All descriptor functions only refer the molecular graph and attributes, other descriptor arrays and global constants. Other external mutable objects are never refered from inside the function.
- Return value of all descriptor functions must be in Julia built-in type (for serialization to JSON)

<!-- #region -->
## Frequently used descriptors


`sssr(mol)`: smallest set of smallest rings (List of rings by node set).
<!-- #endregion -->

```julia
sssr(mol)
```

`isringatom(mol)`: whether the atom is a member of ring

```julia
print(isringatom(mol))
```

`isringbond(mol)`: whether the bond is a member of ring

```julia
print(isringbond(mol))
```

`fusedrings(mol)`: list of fused ring systems by node set

```julia
fusedrings(mol)
```

`Graph.connectedcomponents(mol)`: list of connected components (molecules) by node set.

```julia
Graph.connectedcomponents(mol)
```

`ishacceptor(mol)`: whether the atom is a hydrogen acceptor or not (N, O or F with free lone pairs)

```julia
print(ishacceptor(mol))
```

`ishdonor(mol)`: whether the atom is a hydrogen donor or not (N or O with hydrogens)

```julia
print(ishdonor(mol))
```

`isrotatable(mol)`: whether the bond is rotatable or not

```julia
print(isrotatable(mol))
```

`hybridization(mol)`: hybridization of molecular orbitals

```julia
print(hybridization(mol))
```

`isaromatic(mol)`: whether the atom is in aromatic ring

```julia
print(isaromatic(mol))
```

`isaromaticbond(mol)`: whether the bond is in aromatic ring

```julia
print(isaromaticbond(mol))
```
