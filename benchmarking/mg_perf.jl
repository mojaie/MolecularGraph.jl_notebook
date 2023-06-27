### A Pluto.jl notebook ###
# v0.19.26

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

	function findmcs(mol1, mol2)
	    mol1_ = tdmces_constraints(mol1, :any)
	    mol2_ = tdmces_constraints(mol2, :any)
		#prod = modular_product(mol1_, mol2_)
		# q, _ = maximum_clique(prod, timeout=2)
		# println(length(q))
		#println(nv(prod))
		#println(ne(prod))
	end

	
	for (s1, s2) in [(camphor, adamantanone), (gefitinib, erlotinib), (minocycline, doxorubicin)]
		m1 = smilestomol(s1)
		m2 = smilestomol(s2)
		@time findmcs(m1, m2)
	end
	

	"""
	Profile.clear()
	for i in 1:10
		for (s1, s2) in [(camphor, adamantanone), (gefitinib, erlotinib), (minocycline, doxorubicin)]
			m1 = smilestomol(s1)
			m2 = smilestomol(s2)
			@profile findmcs(m1, m2)
		end
	end
	Profile.print(mincount=20)
	"""
	Nothing
end

# ╔═╡ Cell order:
# ╠═53757a12-0d9e-11ee-1f71-27644c09c451
# ╠═454148c4-d796-4591-9dbe-b8851cf29bee
