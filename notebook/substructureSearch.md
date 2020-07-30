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

# Search molecules from database

- Activate the project and import MolecularGraph.

```julia
using Pkg
Pkg.activate("..")
using MolecularGraph
```

- Download public domain drug dataset provided by [DrugBank](https://drugbank.ca).
- **!!!Caution!!!** the data size is a bit large (25.2 MB).

```julia
# Create data directory
data_dir = "_data"
isdir(data_dir) || mkdir(data_dir)

# Fetch
dest = joinpath(data_dir, "zippedstructures.zip")
url = "https://www.drugbank.ca/releases/5-1-1/downloads/all-open-structures"
isfile(dest) || download(url, dest)

# Unzip
run(`unzip -n -d $data_dir　$dest`);
```

`sdfilereader` loads SDFile text data from the file, and generates array of molecule objects. In this tutorial, the first 2000 molecules were extracted from the file for the test.

Calling `precalculate!` prior to substructure search is highly recommended.

```julia
path =  joinpath(data_dir, "open structures.sdf")
mols = collect(Iterators.take(sdfilereader(path), 2000))
for mol in mols
    precalculate!(mol)
end
```

Then, define convenient function for substructure search. `parse(SMARTS, string)` converts SMARTS strings into query molecule objects. `isquerymatch(mol, query)` compares a molecule-query pair and returns true if they match.

```julia
function substrsearch(smarts)
    matched = []
    query = parse(SMARTS, smarts)
    for (i, m) in enumerate(mols)
        if isquerymatch(m, query)
            push!(matched, m)
            print("@")
        else
             print("+")
        end
        if i % 50 == 0
            println(i)
        end
    end
    println()
    return matched
end
```

## Substructure search


```julia
# Peroxide -O-O-
results1 = substrsearch("[OX2][OX2]")

println("$(length(results1)) records matched")
```

```julia
## Display one of the hit compounds
mol_svg = drawsvg(results1[1], 300, 300)
display("image/svg+xml",  mol_svg)
```

```julia
# Cephem core fused rings
results2 = substrsearch("C1C(=O)N2C=CCSC12")

println("$(length(results2)) records matched")
```

```julia
# Display one of the hit compounds
mol_svg = drawsvg(results2[5], 300, 300)
display("image/svg+xml",  mol_svg)
```
