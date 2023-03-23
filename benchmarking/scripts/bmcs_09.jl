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
    "L-Dopa, 3-Aminocoumarin" => (
        smilestomol("O=C(O)[C@@H](N)Cc1cc(O)c(O)cc1"),
        smilestomol("NC1=Cc2ccccc2OC1=O")
    ),
    "ceftazidime,Cefditoren" => (
        smilestomol(raw"O=C2N1/C(=C(\CS[C@@H]1[C@@H]2NC(=O)C(=NOC(C(=O)O)(C)C)c3nc(sc3)N)C[n+]4ccccc4)C([O-])=O"),
        smilestomol(raw"O=C3N2/C(=C(/C=C\c1scnc1C)CS[C@@H]2[C@@H]3NC(=O)C(=N\OC)/c4nc(sc4)N)C(=O)O")
    ),
    "Minocycline, Doxorubicin" => (
        smilestomol("C[C@H]1[C@H]([C@H](C[C@@H](O1)O[C@H]2C[C@@](Cc3c2c(c4c(c3O)C(=O)c5cccc(c5C4=O)OC)O)(C(=O)CO)O)N)O"),
        smilestomol("CN(C)c1ccc(c2c1C[C@H]3C[C@H]4[C@@H](C(=C(C(=O)[C@]4(C(=C3C2=O)O)O)C(=O)N)O)N(C)C)O")
    )
]

data = DataFrame(query=[])
for (name, query) in queries
    (a, b) = query
    suite = BenchmarkGroup()
    suite["substruct"] = @benchmarkable !isempty(structmatches($a, $b, :substruct))
    suite["Connected MCIS"] = @benchmarkable !isempty(mcismol($a, $b, connected=true, timeout=1))
    suite["Connected MCES"] = @benchmarkable !isempty(mcesmol($a, $b, connected=true, timeout=1))
    suite["Disconnected MCIS"] = @benchmarkable !isempty(mcismol($a, $b, timeout=1))
    suite["Disconnected MCES"] = @benchmarkable !isempty(mcesmol($a, $b, timeout=1))
    suite["tcMCIS"] = @benchmarkable !isempty(mcismol($a, $b, topological=true, timeout=1))
    suite["tcMCES"] = @benchmarkable !isempty(mcesmol($a, $b, topological=true, timeout=1))
    tune!(suite)
    res = run(suite, verbose = true, seconds = 5)
    row = time(median(res)).data
    row["query"] = name
    push!(data, row, cols=:union)
end

display(data)
outpath = joinpath(ENV["OUTPUT_DIR"], "bmcs.csv")
CSV.write(outpath, data)
