
using Pkg
Pkg.develop(PackageSpec(path=ENV["MOLGJ_DEV_PATH"]))

using DotEnv
DotEnv.config(joinpath(ENV["PWD"], ".env"))

using Revise
using Profile
using MolecularGraph

qr = query_relationship();
dox = smilestomol(raw"C[C@H]1[C@H]([C@H](C[C@@H](O1)O[C@H]2C[C@@](Cc3c2c(c4c(c3O)C(=O)c5cccc(c5C4=O)OC)O)(C(=O)CO)O)N)O");
doxh = removestereohydrogens(dox);
precalculate!(doxh);
precalculate!(dox);
qn = smartstomol(raw"[$([o,n]=c1ccc(=[o,n])cc1),$([O,N]=C1C=CC(=[O,N])C=C1),$([O,N]=C1[#6]:,=[#6]C(=[O,N])[#6]:,=[#6]1)]");

@time filter_queries(qr, doxh);

Profile.clear();
@profile filter_queries(qr, doxh, filtering=false);
Profile.print(mincount=50);

flt = filter_queries(qr, doxh);
for na in Graph.nodeattrs(flt)
    delete!(na, "parsed")
    println(na)
end

println([Graph.nodeattr(flt2, i)["key"] for i in 1:Graph.nodecount(flt2)])


hassubstructmatch(smartstomol(raw"O=CN=[N+]=[N-]"), smartstomol(raw"N=[N+]=[N-]"))
