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
    "Doxorubicin, Minocycline" => (
        smilestomol("C[C@H]1[C@H]([C@H](C[C@@H](O1)O[C@H]2C[C@@](Cc3c2c(c4c(c3O)C(=O)c5cccc(c5C4=O)OC)O)(C(=O)CO)O)N)O"),
        smilestomol("CN(C)c1ccc(c2c1C[C@H]3C[C@H]4[C@@H](C(=C(C(=O)[C@]4(C(=C3C2=O)O)O)C(=O)N)O)N(C)C)O")
    )
]

data = DataFrame(query=[])
for (name, query) in queries
    (a, b) = query
    suite = Dict()
    # suite["substruct"] = @btime !isempty(structmatches($a, $b, :substruct))
    # suite["Connected MCIS"] = @btime !isempty(connectedmcis($a, $b, timeout=5))
    # suite["Connected MCES"] = @btime !isempty(connectedmces($a, $b, timeout=5))
    # suite["Disconnected MCIS"] = @btime !isempty(disconnectedmcis($a, $b, timeout=5))
    # suite["Disconnected MCES"] = @btime !isempty(disconnectedmces($a, $b, timeout=5))
    suite["tcMCIS"] = @btime !isempty(tcmcis($a, $b))
    suite["tcMCES"] = @btime !isempty(tcmces($a, $b))
    suite["tcMCISd8"] = @btime !isempty(tcmcis($a, $b, diameter=8))
    suite["tcMCESd8"] = @btime !isempty(tcmces($a, $b, diameter=8))
    println()
    # suite["query"] = name
    # push!(data, suite, cols=:union)
end

# display(data)
# outpath = joinpath(ENV["OUTPUT_DIR"], "head.csv")
# CSV.write(outpath, data)
