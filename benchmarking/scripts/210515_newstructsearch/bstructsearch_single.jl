using DotEnv
DotEnv.config()

using Pkg
Pkg.develop(PackageSpec(path=ENV["MOLGJ_DEV_PATH"]))

using MolecularGraph
using MolecularGraph.Graph
using BenchmarkTools
using DataFrames
using CSV


queries = [
    "[Si]" => smilestomol("[Si]"),
    "[#16;X2;!R]" => smartstomol("[#16;X2;!R]"),
    raw"[NX3;H2,H1;!$(NC=O)]" => smartstomol(raw"[NX3;H2,H1;!$(NC=O)]"),
    "SN" => smilestomol("SN"),
    "[CX3]=[OX1]" => smartstomol("[CX3]=[OX1]"),
    raw"[$([CX3]([#6])[#6]),$([CX3H][#6])]=[$([NX2][#6]),$([NX2H])]" => smartstomol(raw"[$([CX3]([#6])[#6]),$([CX3H][#6])]=[$([NX2][#6]),$([NX2H])]")
]


srcpath = "DrugBank 5.1.8/open structures.sdf"
fpath = joinpath(ENV["CHEMICAL_LIBRARY_PATH"], srcpath)

mols = collect(Iterators.take(sdfilereader(fpath), 1000))
for mol in mols
    precalculate!(mol)
end

data = DataFrame(query=[])
for (name, query) in queries
    suite = BenchmarkGroup()
    suite["substruct"] = @benchmarkable [!isempty(structmatches(m, $query, :substruct, prefilter=false)) for m in $mols]
    suite["fastsingle"] = @benchmarkable [!isempty(structmatches(m, $query, :substruct, prefilter=false, fastsingleton=true)) for m in $mols]
    suite["esubstruct"] = @benchmarkable [!isempty(structmatches(m, $query, :exact, prefilter=false)) for m in $mols]
    suite["efastsingle"] = @benchmarkable [!isempty(structmatches(m, $query, :exact, prefilter=false, fastsingleton=true)) for m in $mols]
    tune!(suite)
    res = run(suite, verbose = true, seconds = 1)
    row = time(median(res)).data
    row["query"] = name
    push!(data, row, cols=:union)
end

display(data)
outpath = joinpath(ENV["OUTPUT_DIR"], "bstructsearch_single.csv")
CSV.write(outpath, data)

# connected MCS (MCIS)
# connected MCS (MCES)
# disconnected MCS (MCIS)
# disconnected MCS (MCES)
# tcMCS (MCIS)
# tcMCS (MCES)