### A Pluto.jl notebook ###
# v0.20.25

using Markdown
using InteractiveUtils

# ╔═╡ 9e2d5cf8-58ff-11f1-9238-3b6f6e109168
begin
	using Revise
	using Pkg; Pkg.activate("../")
	using JLD2
	using LatticeProteins
	using Random
	using StatsPlots
end

# ╔═╡ f3e89f95-742e-4006-b3b1-6eec05bc73b9
LP = LatticeProteins

# ╔═╡ bbba7831-f8f1-4375-99e6-3c25888ca5e3
@load "../data/structures_N3_filtering0.01.jld2" structures

# ╔═╡ 985367b9-aceb-409d-aea9-85d233abd09d
target = 1

# ╔═╡ ee7980ed-3f0f-4138-a7e0-e2df50b39a03


# ╔═╡ 87ca1d39-346f-4b35-be8b-83b53ee3c1aa


# ╔═╡ d4d49708-3bee-4734-9e93-1f82dc52f0f5
parameters = MCMCParameters(;
	n_steps=1, n_sequences=5, structures, target, burnin=0,
)

# ╔═╡ 5c274d45-edfa-462d-8dd1-6304ab7a919f
init = mod.(0:26, 20) .+ 1

# ╔═╡ 075aea9c-3490-4ed0-ae80-44287b91219d
result = let
	rng = Random.seed!(1)
	sample_mcmc_chain(init, parameters; rng)
end

# ╔═╡ a081403c-ad5e-442b-b7e4-ab04425e6369
result[2]

# ╔═╡ 39def697-b623-49a4-a790-5769359fab01
chain = result[1]

# ╔═╡ d7f7feee-d053-446e-8b43-bdc6f57fdc9f
LP.fold_prob(init, structures, target)

# ╔═╡ c136f04e-5201-4a1a-bf97-ba739a0fe1bf
LP.fold_prob(chain[end], structures, target)

# ╔═╡ 98b3fca0-0ac3-4e11-8347-99bc8960effc
count(x -> x.accepted, result[2])

# ╔═╡ 7ddb69ea-1bd4-4709-8f61-233ddda8673a
let
	E = map(S -> LP.energy(init, S), structures)
	P = LP._softmin!(copy(E))
	scatter(P, E)
	plot(1:length(P), cumsum(P))
end

# ╔═╡ c56f8a82-e398-4079-85cf-aba8b94252ad
let
	E = map(S -> LP.energy(chain[end], S), structures)
	P = vcat([0.], LP._softmin!(copy(E)))
	plot(1:length(P), cumsum(P))
end

# ╔═╡ a769c42a-1147-4eae-82b9-96f9b0c9f059
unique(chain) |> length

# ╔═╡ 98f442e3-933e-4b92-8300-dccf4c837212
map(x -> LP.energy(x, structures[target]), chain) |> plot

# ╔═╡ a8bbdc02-1dd1-4e19-97fe-2475e13c34aa
plot(structures[target])

# ╔═╡ Cell order:
# ╠═9e2d5cf8-58ff-11f1-9238-3b6f6e109168
# ╠═f3e89f95-742e-4006-b3b1-6eec05bc73b9
# ╠═bbba7831-f8f1-4375-99e6-3c25888ca5e3
# ╠═985367b9-aceb-409d-aea9-85d233abd09d
# ╠═ee7980ed-3f0f-4138-a7e0-e2df50b39a03
# ╠═87ca1d39-346f-4b35-be8b-83b53ee3c1aa
# ╠═d4d49708-3bee-4734-9e93-1f82dc52f0f5
# ╠═5c274d45-edfa-462d-8dd1-6304ab7a919f
# ╠═075aea9c-3490-4ed0-ae80-44287b91219d
# ╠═a081403c-ad5e-442b-b7e4-ab04425e6369
# ╠═39def697-b623-49a4-a790-5769359fab01
# ╠═d7f7feee-d053-446e-8b43-bdc6f57fdc9f
# ╠═c136f04e-5201-4a1a-bf97-ba739a0fe1bf
# ╠═98b3fca0-0ac3-4e11-8347-99bc8960effc
# ╠═7ddb69ea-1bd4-4709-8f61-233ddda8673a
# ╠═c56f8a82-e398-4079-85cf-aba8b94252ad
# ╠═a769c42a-1147-4eae-82b9-96f9b0c9f059
# ╠═98f442e3-933e-4b92-8300-dccf4c837212
# ╠═a8bbdc02-1dd1-4e19-97fe-2475e13c34aa
