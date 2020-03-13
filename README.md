
MolecularGraph.jl - Notebooks for tutorial
===========================================


Environment requirement
------------------------

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
  | | |_| | | | (_| |  |  Version 1.3.0 (2019-11-26)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

julia>
```

3. Type `]` to enter Pkg mode and create project

```
(v1.3) pkg> activate .
Activating new environment at `your_workspace/MolecularGraph_notebook/Project.toml`
```

4. Install packages

```
(v1.3) pkg> add MolecularGraph
Updating registry at `~/.julia/registries/General`
Updating git-repo `https://github.com/JuliaRegistries/General.git`
Resolving package versions...
Updating `~/Workspace/MolecularGraph.jl_notebook/Project.toml`
[6c89ec66] + MolecularGraph v0.3.2
Updating `~/Workspace/MolecularGraph.jl_notebook/Manifest.toml`
[19ecbf4d] + Codecs v0.5.0
[34da2185] + Compat v2.2.0
[ffbed154] + DocStringExtensions v0.8.1
[e30172f5] + Documenter v0.24.6
...
```


Getting started
------------------------

- [Getting started](https://nbviewer.jupyter.org/github/mojaie/graphmol.jl/blob/master/notebook/gettingStarted.ipynb)
- [Search molecules from database](https://nbviewer.jupyter.org/github/mojaie/graphmol.jl/blob/master/notebook/substructureSearch.ipynb)
- [Maximum common substructure (MCS)](https://nbviewer.jupyter.org/github/mojaie/graphmol.jl/blob/master/notebook/mcs.ipynb)



License
------------------------

[CC-by-4.0](https://creativecommons.org/licenses/by/4.0/)



Copyright
------------------------

(C) 2020 Seiji Matsuoka
