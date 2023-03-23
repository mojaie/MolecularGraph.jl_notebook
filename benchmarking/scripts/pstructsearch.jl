using DotEnv
DotEnv.config()

using Pkg
Pkg.develop(PackageSpec(path=ENV["MOLGJ_DEV_PATH"]))

using MolecularGraph
using MolecularGraph.Graph
using Profile

srcpath = "DrugBank 5.1.8/open structures.sdf"
fpath = joinpath(ENV["CHEMICAL_LIBRARY_PATH"], srcpath)

mols = collect(Iterators.take(sdfilereader(fpath), 200))
for mol in mols
    precalculate!(mol)
end

# query = smartstomol("[#6][OD2][#6]")
# query = smilestomol("CC1=C2[C@@]([C@]([C@H]([C@@H]3[C@]4([C@H](OC4)C[C@@H]([C@]3(C(=O)[C@@H]2OC(=O)C)C)O)OC(=O)C)OC(=O)c5ccccc5)(C[C@@H]1OC(=O)[C@H](O)[C@@H](NC(=O)c6ccccc6)c7ccccc7)O)(C)C")
#
# query = smilestomol("C=1C=CC=CC1")
query = smartstomol("O.O.O.O.O.O.OS(O)(=O)=O")

# Profile.clear()
@profile [!isempty(structmatches(mol, query, :substruct)) for mol in mols]
Profile.print(mincount=100)
