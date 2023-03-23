cd(@__DIR__)
using Pkg
Pkg.activate(".")
Pkg.develop(PackageSpec(path="../MolecularGraph.jl"))

using MolecularGraph
using BenchmarkTools
using DelimitedFiles

datafile = "DrugBank_FDA_Approved.sdf"
# datafile = "open structures.sdf"

# load molecules
fpath = joinpath(
    dirname(@__FILE__), "../../Workspace/resource/chemlib/", datafile)


function benchmarksuite(mols)
    suite = BenchmarkGroup()
    # suite["sdfilereader"] = @benchmarkable collect(sdfilereader(fpath))

    # Tier 1
    suite["nodedegree"] = @benchmarkable [nodedegree(mol) for mol in $mols]
    suite["atomsymbol"] = @benchmarkable atomsymbol.($mols)
    suite["charge"] = @benchmarkable charge.($mols)
    suite["multiplicity"] = @benchmarkable multiplicity.($mols)
    suite["bondorder"] = @benchmarkable bondorder.($mols)
    suite["connectedcomponents"
        ] = @benchmarkable [connectedcomponents(mol) for mol in $mols]
    suite["fusedrings"] = @benchmarkable fusedrings.($mols)
    suite["inchi"] = @benchmarkable inchi.($mols)
    suite["sdfcoords2d"] = @benchmarkable sdfcoords2d.($mols)
    suite["coordgen"] = @benchmarkable coordgen.($mols)
    suite["isplanar"] = @benchmarkable Graph.isplanar.($mols)
    suite["isouterplanar"] = @benchmarkable Graph.isouterplanar.($mols)
    suite["issimplegraph"] = @benchmarkable Graph.issimplegraph.($mols)

    # Tier 2
    # Graph.edgemincycles ->
    suite["sssr"] = @benchmarkable sssr.($mols)
    # Graph.adjacencies, atomsymbol ->
    suite["explicithconnected"] = @benchmarkable explicithconnected.($mols)
    # atomsymbol, charge ->
    suite["lonepair"] = @benchmarkable lonepair.($mols)
    # bondorder ->
    suite["apparentvalence"] = @benchmarkable MolecularGraph.apparentvalence.($mols)
    # connectedcomponents ->
    suite["connectedmembership"] = @benchmarkable connectedmembership.($mols)
    # fusedrings ->
    suite["fusedringmembership"] = @benchmarkable fusedringmembership.($mols)
    # inchi ->
    suite["inchikey"] = @benchmarkable inchikey.($mols)
    # atomsymbol ->
    suite["atomcolor"] = @benchmarkable atomcolor.($mols)
    # atomsymbol,  nodedegree ->
    suite["isatomvisible"] = @benchmarkable isatomvisible.($mols)

    # Tier 3
    # sssr ->
    suite["sssrmembership"] = @benchmarkable sssrmembership.($mols)
    # lonepair, atomsymbol ->
    suite["ishacceptor"] = @benchmarkable ishacceptor.($mols)
    # lonepair, atomsymbol ->
    suite["valence"] = @benchmarkable valence.($mols)
    # explicithconnected, nodedegree ->
    suite["heavyatomconnected"] = @benchmarkable heavyatomconnected.($mols)
    # apparentvalence, atomsymbol, charge, Graph.adjacencies ->
    suite["pielectron"] = @benchmarkable pielectron.($mols)
    # sssr, apparentvalence, lonepair, nodedegree, atomsymbol, Graph.adjacencies ->
    suite["isaromaticring"] = @benchmarkable isaromaticring.($mols)

    # Tier 4
    # valence, apparentvalence ->
    suite["implicithconnected"] = @benchmarkable implicithconnected.($mols)
    # sssrmembership ->
    suite["sssrsizes"] = @benchmarkable sssrsizes.($mols)
    # sssrmembership ->
    suite["sssrcount"] = @benchmarkable sssrcount.($mols)
    # sssrmembership ->
    suite["isringatom"] = @benchmarkable isringatom.($mols)
    # Graph.edgemincyclemembership ->
    suite["isringbond"] = @benchmarkable isringbond.($mols)
    # isaromaticring, sssr ->
    suite["isaromatic"] = @benchmarkable isaromatic.($mols)
    # isaromaticring, sssr ->
    suite["isaromaticbond"] = @benchmarkable isaromaticbond.($mols)

    # Tier 5
    # implicithconnected, explicithcount ->
    suite["hydrogenconnected"] = @benchmarkable hydrogenconnected.($mols)
    # implicithconnected, nodedegree ->
    suite["connectivity"] = @benchmarkable connectivity.($mols)
    # implicithconnected, atomsymbol ->
    suite["atomcounter"] = @benchmarkable atomcounter.($mols)
    # isringbond, nodedegree, bondorder ->
    suite["isrotatable"] = @benchmarkable isrotatable.($mols)
    # implicithconnected ->
    suite["standardweight"] = @benchmarkable standardweight.($mols)
    # implicithconnected ->
    suite["monoisotopicmass"] = @benchmarkable monoisotopicmass.($mols)
    # isringbond, sssr, (coordgen,) bondorder ->
    suite["coords2d"] = @benchmarkable coords2d.($mols)

    # Tier 6
    # hydrogenconnected, atomsymbol ->
    suite["ishdonor"] = @benchmarkable ishdonor.($mols)
    # connectivity, pielectron, lonepair, atomsymbol ->
    suite["hybridization"] = @benchmarkable hybridization.($mols)
    # atomcounter ->
    suite["molecularformula"] = @benchmarkable molecularformula.($mols)
    # atomcounter ->
    suite["empiricalformula"] = @benchmarkable empiricalformula.($mols)

    # Tier 7
    # hybridization, heavyatomconnected, pielectron, isaromatic ->
    suite["wclogptype"] = @benchmarkable wclogptype.($mols)
    # hybridization, isaromatic, heavyatomconnected ->
    suite["wclogphydrogentype"] = @benchmarkable wclogphydrogentype.($mols)
    # wclogptype, wclogphydrogentype ->
    suite["wclogpcontrib"] = @benchmarkable wclogpcontrib.($mols)

    return suite
end


# Naive
mols = collect(sdfilereader(fpath))[1:100:2000]
res1 = run(benchmarksuite(mols), verbose = true, samples = 100)


# SSSR cache
mols = collect(sdfilereader(fpath))[1:100:2000]
for mol in mols
    Graph.@cache Graph.edgemincycles(mol)
    Graph.@cache sssr(mol)
end
res2 = run(benchmarksuite(mols), verbose = true, samples = 100)


# SMARTS parameter caches
mols = collect(sdfilereader(fpath))[1:100:2000]
for mol in mols
    Graph.@cache Graph.edgemincycles(mol)
    Graph.@cache sssr(mol)
    Graph.@cache lonepair(mol)
    Graph.@cache valence(mol)
    Graph.@cache MolecularGraph.apparentvalence(mol)
end
res3 = run(benchmarksuite(mols), verbose = true, samples = 100)


# Output
fkeys = sort([k for (k, v) in res1])
ktoi = Dict(k => i for (i, k) in enumerate(fkeys))
data = hcat(fkeys, zeros(Float64, length(fkeys), 3))
for (i, r) in enumerate((res1, res2, res3))
    for (k, v) in r
        data[ktoi[k], i + 1] = time(median(v))
    end
end

writedlm("./output/data.csv",  data, ',')
