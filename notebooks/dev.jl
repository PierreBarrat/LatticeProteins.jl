### A Pluto.jl notebook ###
# v0.20.25

using Markdown
using InteractiveUtils

# ╔═╡ 51db6637-b4c8-44e8-a6d9-6a3be937ebb3
begin
	using Revise
	using Pkg
	Pkg.activate("..")
end

# ╔═╡ 5f6b0015-170d-41e7-bb1a-44f4a8bdfaa2
using LatticeProteins

# ╔═╡ 0db7daf2-d750-4d01-80cd-f33a8fd53614
begin
	using Plots
end

# ╔═╡ b32d4c2a-43d1-47a4-bd3d-77479088bc5d
LPS = LatticeProteins.Structures

# ╔═╡ eb7c05f6-0982-433c-9756-49e447d12033
# ╠═╡ disabled = true
#=╠═╡
chains = let
	LatticeProteins.Structures.generate_all_paths(3)
end
  ╠═╡ =#

# ╔═╡ 15863399-8360-4823-bd90-7670d380ad03
paths = begin 
	N = 2
	adj = LPS.build_adj(N)
	start = LPS.STARTS_2[1]
	visited = zeros(Bool, N, N, N)
	visited[start.x, start.y, start.z] = true
	rotation_table = LPS.generate_rotation_table(N)
	paths = Set{Vector{LPS.Site}}()
	LPS.grow!([start], visited, adj, rotation_table, paths)
	paths
end |> collect

# ╔═╡ a3277627-b691-4929-95b8-3f178c628072
length(paths)

# ╔═╡ 0b8eb74f-afac-42d5-a2ee-9c863180ad38


# ╔═╡ 427236ee-a570-48e9-b727-3bbc35f81e5a
let
	plts = map(plot, collect(paths))
	plot(plts..., layout=grid(3,1), size=(900, 900))
end

# ╔═╡ 5fc55395-1514-41f8-b26a-712d6d76cbce
plot(p)

# ╔═╡ dbc2dce5-a9e8-48d5-aeea-ec3c12584c75
let
	table = LPS.generate_rotation_table(N)
	plts = []
	for rot in LPS.rotations(p, table)
		push!(plts, plot(rot))
	end
	plts
	plot(
		plts...;
		layout=grid(8, 6),
		size=(1000, 1000),
		margin=-2*Plots.PlotMeasures.mm
	)
end

# ╔═╡ 3a83a489-bb38-41fd-b9ca-20378d8fa6b0
p

# ╔═╡ 06dfa1c8-6fff-41fd-8c28-1aa39f68a34d
table = LPS.generate_rotation_table(N)

# ╔═╡ 24eb88c4-1869-4b6b-a56d-69604f7d3dfb
plot(LPS.canonical(p, table))

# ╔═╡ 535145f4-7034-4656-b511-0eb4053c7dc1
table[4]

# ╔═╡ d72817ca-ec11-41b1-a81d-af9d77cec7a4
rot = table[1]

# ╔═╡ 4ea397ca-856e-4972-9694-7a2559cdbc32
let
	p2 = copy(p)
	for (i, site) in enumerate(p)
		p2[i] = rot[site]
	end
	p2 == p
end

# ╔═╡ b438f2b9-3564-47cb-8b8a-bf020f611f31


# ╔═╡ 9ad5a02b-2b8e-4d19-835e-0b39cf948562
function plot_chain(chain)
	dot_size = 8
	start_color = :green
	end_color = :red
	
	p = plot()
	for x in 1:3, y in 1:3, z in 1:3
		scatter!(
			p, [x], [y], [z]; 
			label="", marker=(dot_size, :circle, :blue, stroke(0))
		)
	end
	for i in 2:length(chain)
		x0, y0, z0 = chain[i-1]
		x1, y1, z1 = chain[i]
		plot!(p, [x0, x1], [y0, y1], [z0, z1], line=(:arrow, :black, 2), label="")
	end

	# mark start and end
	x, y, z = chain[1]
	scatter!(
		p, [x], [y], [z]; 
		label="", marker=(dot_size, :circle, start_color, stroke(0))
	)
	x, y, z = chain[end]
	scatter!(
		p, [x], [y], [z]; 
		label="", marker=(dot_size, :circle, end_color, stroke(0))
	)
	p
end

# ╔═╡ 34aeca8a-eb44-4880-84f0-d530057c0231
plot_chain(first(paths))

# ╔═╡ Cell order:
# ╠═51db6637-b4c8-44e8-a6d9-6a3be937ebb3
# ╠═5f6b0015-170d-41e7-bb1a-44f4a8bdfaa2
# ╠═0db7daf2-d750-4d01-80cd-f33a8fd53614
# ╠═b32d4c2a-43d1-47a4-bd3d-77479088bc5d
# ╠═eb7c05f6-0982-433c-9756-49e447d12033
# ╠═15863399-8360-4823-bd90-7670d380ad03
# ╠═a3277627-b691-4929-95b8-3f178c628072
# ╠═0b8eb74f-afac-42d5-a2ee-9c863180ad38
# ╠═427236ee-a570-48e9-b727-3bbc35f81e5a
# ╠═5fc55395-1514-41f8-b26a-712d6d76cbce
# ╠═24eb88c4-1869-4b6b-a56d-69604f7d3dfb
# ╠═535145f4-7034-4656-b511-0eb4053c7dc1
# ╠═dbc2dce5-a9e8-48d5-aeea-ec3c12584c75
# ╠═3a83a489-bb38-41fd-b9ca-20378d8fa6b0
# ╠═06dfa1c8-6fff-41fd-8c28-1aa39f68a34d
# ╠═d72817ca-ec11-41b1-a81d-af9d77cec7a4
# ╠═4ea397ca-856e-4972-9694-7a2559cdbc32
# ╠═34aeca8a-eb44-4880-84f0-d530057c0231
# ╠═b438f2b9-3564-47cb-8b8a-bf020f611f31
# ╠═9ad5a02b-2b8e-4d19-835e-0b39cf948562
