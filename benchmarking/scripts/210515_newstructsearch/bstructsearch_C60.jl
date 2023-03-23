using DotEnv
DotEnv.config()

using Pkg
Pkg.develop(PackageSpec(path=ENV["MOLGJ_DEV_PATH"]))

using MolecularGraph
using MolecularGraph.Graph
using BenchmarkTools
using DataFrames
using CSV

mol = smilestomol("c12c3c4c5c1c6c7c8c2c9c1c3c2c3c4c4c%10c5c5c6c6c7c7c%11c8c9c8c9c1c2c1c2c3c4c3c4c%10c5c5c6c6c7c7c%11c8c8c9c1c1c2c3c2c4c5c6c3c7c8c1c23")
precalculate!(mol)
query = smilestomol("C=1C=CC=CC1")


println()
display(@benchmark !isempty(structmatches($mol, $query, :exact)))
println()
display(@benchmark !isempty(structmatches($mol, $query, :substruct)))
println()
