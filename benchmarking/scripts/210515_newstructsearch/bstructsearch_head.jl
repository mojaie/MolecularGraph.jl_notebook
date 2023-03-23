using DotEnv
DotEnv.config()

using Pkg
Pkg.develop(PackageSpec(path=ENV["MOLGJ_DEV_PATH"]))

using MolecularGraph
using MolecularGraph.Graph
using BenchmarkTools
using DataFrames
using CSV

srcpath = "DrugBank 5.1.8/open structures.sdf"
fpath = joinpath(ENV["CHEMICAL_LIBRARY_PATH"], srcpath)

mols = collect(Iterators.take(sdfilereader(fpath), 1000))
for mol in mols
    precalculate!(mol)
end

queries = [
    "C=1C=CC=CC1" => smilestomol("C=1C=CC=CC1"),
    "taxol" => smilestomol("CC1=C2[C@@]([C@]([C@H]([C@@H]3[C@]4([C@H](OC4)C[C@@H]([C@]3(C(=O)[C@@H]2OC(=O)C)C)O)OC(=O)C)OC(=O)c5ccccc5)(C[C@@H]1OC(=O)[C@H](O)[C@@H](NC(=O)c6ccccc6)c7ccccc7)O)(C)C"),
    "[#6][OD2][#6]" => smartstomol("[#6][OD2][#6]"),
    "O.O.OS(O)(=O)=O" => smartstomol("O.O.OS(O)(=O)=O")
]

data = DataFrame(query=[])
for (name, query) in queries
    suite = BenchmarkGroup()
    suite["exact"] = @benchmarkable [!isempty(structmatches(m, $query, :exact)) for m in $mols]
    suite["substruct"] = @benchmarkable [!isempty(structmatches(m, $query, :substruct)) for m in $mols]
    suite["nodeinduced"] = @benchmarkable [!isempty(structmatches(m, $query, :nodeinduced)) for m in $mols]
    suite["edgeinduced"] = @benchmarkable [!isempty(structmatches(m, $query, :edgeinduced)) for m in $mols]
    tune!(suite)
    res = run(suite, verbose = true, seconds = 1)
    row = time(median(res)).data
    row["query"] = name
    push!(data, row, cols=:union)
end

display(data)
outpath = joinpath(ENV["OUTPUT_DIR"], "bstructsearch_head.csv")
CSV.write(outpath, data)

# connected MCS (MCIS)
# connected MCS (MCES)
# disconnected MCS (MCIS)
# disconnected MCS (MCES)
# tcMCS (MCIS)
# tcMCS (MCES)