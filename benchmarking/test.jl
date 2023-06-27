### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 5b4c1aac-ee43-11ed-1381-8d31643817a2
begin
	import Pkg
	Pkg.activate(Base.current_project())
	Pkg.instantiate()
	using Graphs, MolecularGraph
	Pkg.status()
end

# ╔═╡ 6959b26c-48bc-41ed-8daf-80351443bb94
c60 = smilestomol("c12c3c4c5c1c6c7c8c2c9c1c3c2c3c4c4c%10c5c5c6c6c7c7c%11c8c9c8c9c1c2c1c2c3c4c3c4c%10c5c5c6c6c7c7c%11c8c8c9c1c1c2c3c2c4c5c6c3c7c8c1c23")

# ╔═╡ f13d274d-a450-4fb4-9514-a1f92a33b3f8
kekulize(c60)

# ╔═╡ Cell order:
# ╠═5b4c1aac-ee43-11ed-1381-8d31643817a2
# ╠═6959b26c-48bc-41ed-8daf-80351443bb94
# ╠═f13d274d-a450-4fb4-9514-a1f92a33b3f8
