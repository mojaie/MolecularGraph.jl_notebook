### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 53757a12-0d9e-11ee-1f71-27644c09c451
begin
	using Pkg
	Pkg.develop(path="../../MolecularGraph.jl")
	Pkg.activate(Base.current_project())
	Pkg.instantiate()
	using Profile, Graphs, MolecularGraph
	Pkg.status()
end


# ╔═╡ 454148c4-d796-4591-9dbe-b8851cf29bee
let
	camphor = raw"CC1(C)C2CCC1(C)C(=O)C2"
	adamantanone = raw"C1C2CC3CC1CC(C2)C3=O"
	gefitinib = raw"C1COCCN1CCCOc2c(OC)cc3ncnc(c3c2)Nc4cc(Cl)c(F)cc4"
	erlotinib = raw"COCCOc1cc2c(cc1OCCOC)ncnc2Nc3cccc(c3)C#C"
	minocycline = raw"C[C@H]1[C@H]([C@H](C[C@@H](O1)O[C@H]2C[C@@](Cc3c2c(c4c(c3O)C(=O)c5cccc(c5C4=O)OC)O)(C(=O)CO)O)N)O"
	doxorubicin = raw"CN(C)c1ccc(c2c1C[C@H]3C[C@H]4[C@@H](C(=C(C(=O)[C@]4(C(=C3C2=O)O)O)C(=O)N)O)N(C)C)O"

	function dmcs(mol1, mol2)
	    @time mcs = tcmces(mol1, mol2, timeout=2)
	    print(size(mcs))
	end

	for (s1, s2) in [(camphor, adamantanone), (gefitinib, erlotinib), (minocycline, doxorubicin)]
	    m1 = smilestomol(s1)
	    m2 = smilestomol(s2)
		dmcs(m1, m2)
	end
end

# ╔═╡ d9538764-29ba-435e-a848-8ef823eef7dc
let
	camphor = raw"CC1(C)C2CCC1(C)C(=O)C2"
	adamantanone = raw"C1C2CC3CC1CC(C2)C3=O"
	gefitinib = raw"C1COCCN1CCCOc2c(OC)cc3ncnc(c3c2)Nc4cc(Cl)c(F)cc4"
	erlotinib = raw"COCCOc1cc2c(cc1OCCOC)ncnc2Nc3cccc(c3)C#C"
	minocycline = raw"C[C@H]1[C@H]([C@H](C[C@@H](O1)O[C@H]2C[C@@](Cc3c2c(c4c(c3O)C(=O)c5cccc(c5C4=O)OC)O)(C(=O)CO)O)N)O"
	doxorubicin = raw"CN(C)c1ccc(c2c1C[C@H]3C[C@H]4[C@@H](C(=C(C(=O)[C@]4(C(=C3C2=O)O)O)C(=O)N)O)N(C)C)O"

	function dmcs(mol1, mol2)
	    mcs = tcmces(mol1, mol2, timeout=2)
	    print(size(mcs))
	end

	Profile.clear()
	for (s1, s2) in [(camphor, adamantanone), (gefitinib, erlotinib), (minocycline, doxorubicin)]
	    m1 = smilestomol(s1)
	    m2 = smilestomol(s2)
	    @profile dmcs(m1, m2)
		#@time dmcs(m1, m2)
	end

	Profile.print(mincount=20)
end

# ╔═╡ Cell order:
# ╠═53757a12-0d9e-11ee-1f71-27644c09c451
# ╠═454148c4-d796-4591-9dbe-b8851cf29bee
# ╠═d9538764-29ba-435e-a848-8ef823eef7dc
