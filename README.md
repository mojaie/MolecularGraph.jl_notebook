
MolecularGraph.jl tutorial notebooks (Pluto.jl)
===========================================

To run codes in your environment, see `Edit or run this notebook` instruction shown in the top-right of the tutorial pages below.

- [Getting started](https://mojaie.github.io/MolecularGraph.jl_notebook/getting_started.jl.html)
- [Molecular graph basics](https://mojaie.github.io/MolecularGraph.jl_notebook/molecular_graph_basics.jl.html)
- [Properties and descriptors](https://mojaie.github.io/MolecularGraph.jl_notebook/properties_and_descriptors.jl.html)
- [Preprocessing](https://mojaie.github.io/MolecularGraph.jl_notebook/preprocessing.jl.html)
- [Mass and isotopes](https://mojaie.github.io/MolecularGraph.jl_notebook/mass_and_isotopes.jl.html)



Legacy tutorials (MolecularGraph.jl v0.13 and below)
-----------------------------------------


1. Install Julia kernel to Jupyter Notebook

If you are new to Jupyter Notebook with Julia kernel, please set up Julia kernel according to IJulia instruction.  
https://github.com/JuliaLang/IJulia.jl

1. Clone the repository

```
% cd your_workspace
% git clone https://github.com/mojaie/MolecularGraph.jl_notebook.git
Cloning into 'MolecularGraph.jl_notebook'...
```

2. Launch Julia REPL

```
% cd MolecularGraph_notebook
% julia

               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.6.0 (2021-03-24)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

julia>
```

3. Type `]` to enter Pkg mode and create project

```
(@1.6) pkg> activate .
Activating new environment at `your_workspace/MolecularGraph_notebook/Project.toml`
```

4. Install packages

`instantiate` the notebook tutorial project. If you do not have 'Plot.jl' yet, it may take several minites to install.

```
(MolecularGraph) pkg> instantiate
  Progress [========================================>]  1/1
1 dependency successfully precompiled in 8 seconds (16 already precompiled)
```


- [Getting started](https://nbviewer.jupyter.org/github/mojaie/MolecularGraph.jl_notebook/blob/master/notebook_v0_13/gettingStarted.ipynb)
- [Basics of molecular graph](https://nbviewer.jupyter.org/github/mojaie/MolecularGraph.jl_notebook/blob/master/notebook_v0_13/molecularGraphBasics.ipynb)
- [Preprocessing](https://nbviewer.jupyter.org/github/mojaie/MolecularGraph.jl_notebook/blob/master/notebook_v0_13/preprocess.ipynb)
- [Calculation of descriptors](https://nbviewer.jupyter.org/github/mojaie/MolecularGraph.jl_notebook/blob/master/notebook_v0_13/calculateDescriptors.ipynb)
- [Molecular/atomic mass and isotopes](https://nbviewer.jupyter.org/github/mojaie/MolecularGraph.jl_notebook/blob/master/notebook_v0_13/massAndIsotopes.ipynb)
- [Search molecules from database](https://nbviewer.jupyter.org/github/mojaie/MolecularGraph.jl_notebook/blob/master/notebook_v0_13/substructureSearch.ipynb)
- [Maximum common substructure (MCS)](https://nbviewer.jupyter.org/github/mojaie/MolecularGraph.jl_notebook/blob/master/notebook_v0_13/mcs.ipynb)



License
------------------------

[CC-by-4.0](https://creativecommons.org/licenses/by/4.0/)



Copyright
------------------------

(C) 2020-2023 Seiji Matsuoka
