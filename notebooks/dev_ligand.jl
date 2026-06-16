### A Pluto.jl notebook ###
# v1.0.1

using Markdown
using InteractiveUtils

# ╔═╡ 5cd13410-64b0-11f1-81f3-713fcda94f0a
begin
	using Revise
	using Pkg; Pkg.activate("..")
	using ArgCheck
	using LatticeProteins
	using Plots
	using StatsBase
	# using SimplePlutoInclude
end

# ╔═╡ 62f1f8f2-a51b-4d3a-b82a-c1b1b6a6a710
LP = LatticeProteins

# ╔═╡ 53ec092a-fb32-421c-aa7f-550acec1e6c2
letters = split("ABCDE", "")

# ╔═╡ ccce05dd-77ad-40d0-b1ff-256b03ac1fce
ligand_cross = Ligand(;beads=1:5, shape=:cross)

# ╔═╡ 195dfa45-af1d-4b42-911a-7d3d2f50835a
md"""
# Symmetries 
"""

# ╔═╡ 1ca78bc1-6281-4bb1-bc1c-7cac9efd0171
function plot_cells(cells) 
    p = plot(xlim=(-1.5, 1.5), ylim=(-1.5, 1.5))
    cells = [(c[1], c[2], letters[i]) for (i,c) in enumerate(cells)]
    for c in cells 
        annotate!(p, c) 
    end
    p
end

# ╔═╡ 6eeaaea4-e20d-4e67-a144-b304dfac607f
cross_plots = map(LP.symmetries(ligand_cross)) do (r, f)
	cells = [LP.rotate(LP.flip(c, f), r) for c in ligand_cross.cells]
	plot_cells(cells)
end

# ╔═╡ 38bb847d-14eb-4d77-a428-fae8071b0005
plot(
	cross_plots..., 
	layout=grid(2, 4), size=(1200, 600)
)

# ╔═╡ 430d9438-ad7f-44fd-87a7-2d99e725a98e
md"""
# Binding patches
"""

# ╔═╡ aabb7efe-aed6-46b2-9bfb-9a15b3527349
md"""
## Relative to a face
"""

# ╔═╡ 66384e93-a183-4064-986a-a716a26d6fb4
function plot_patch(patch, N; kwargs...)
	p = plot(;
		xticks=1:N, 
		yticks=1:N,
		gridalpha=1,
		gridwidth=2,
		xlim=(0.5, N+0.5),
		ylim=(0.5, N+0.5),
		kwargs...
	)
	pal = palette(:tab10)
	for (i, r) in enumerate(patch)
		scatter!([r[1]], [r[2]], color=pal[i], label=i, marker=(12))
	end
	p
end

# ╔═╡ 93434855-a17b-48f7-a9c2-7264a55da687
LP.binding_patches_2D(ligand_cross, 3) |> length # should return 8 patches

# ╔═╡ 01b8d561-363b-46ff-b4f6-ce022fd03172
LP.binding_patches_2D(ligand_cross, 4) |> length # should return 32 patches (8x4)

# ╔═╡ 26a17f9a-8330-4a54-9b07-d2554d846ab5
patches_cross_3 = map(LP.binding_patches_2D(ligand_cross, 3)) do patch
	plot_patch(patch, 3)
end;

# ╔═╡ 9ed363e9-72ea-463f-bb64-d9f26b46826b
plot(patches_cross_3..., layout=grid(2, 4), size=(1200, 600))

# ╔═╡ abbc653d-1ee1-4eb2-8ad6-26d11dadb6ad
patches_cross_4 = map(LP.binding_patches_2D(ligand_cross, 4)) do patch
	plot_patch(patch, 4)
end;

# ╔═╡ 30d9d15f-58e4-4bd0-ac3a-61eb8c87a253
plot(patches_cross_4..., layout=grid(4, 8), size=(1200, 600), legend=false)

# ╔═╡ e29ccd3f-3b47-4c10-becb-960f1f60d996
md"""
## All patches
"""

# ╔═╡ 3c8de6d6-21c4-4fe7-935b-d6411704d6e2
LP.binding_patches(ligand_cross, 3) # 8 x 6 = 48 (symmetries x faces)

# ╔═╡ 7f264679-e8a2-4bed-8b7e-7af933471a8f
md"""
# Binding probabilities
"""

# ╔═╡ 5eec4f20-b69f-443d-af1e-a3f2711b8111
sequence = rand(1:20, 27)

# ╔═╡ 5a45302d-fa42-48a8-81f5-c7ec614f586d
LP.MJ_1996[4, 9]

# ╔═╡ 2723c3b5-0425-4830-aae0-e89151fda49e
LP.MJ_1996[3, 4]

# ╔═╡ 33dbed1f-9ea1-4234-8279-4ae68da8772c
sort(vec(LP.MJ_1996))

# ╔═╡ ac0353d0-3efc-4d71-aa31-362215821995
LP.MJ_1996[3, :] |> sortperm

# ╔═╡ 69e781be-16bd-42f5-8ee2-526c6dd918a0
scatter(LP.MJ_1996[4, :] , LP.MJ_1996[3, :] )

# ╔═╡ 104c0512-c26f-43cc-8724-e5a27c2a8f4b
LP.MJ_1996[3, :] |> sortperm

# ╔═╡ 1d77d89a-97ec-4ab2-9add-13672dd0729a
LP.MJ_1996[9, :] |> sortperm

# ╔═╡ 4858a6b2-dd8b-49d1-a5ef-9458e8e9478f
LP.MJ_1996[10, :] |> mean

# ╔═╡ 5f969c85-e8ee-4355-b63c-c7a68b2f2abe
md"""
# MCMC
"""

# ╔═╡ 93c8f931-eaa7-4486-a6b6-f8bab43466d0
f = 0.1

# ╔═╡ 8689f611-15e3-4b55-b980-877cccf3ec0d
structures =  LP.Structures.generate_structures(3, f);

# ╔═╡ 9a727f50-efd4-48bf-bc26-d3b3fd8f8878
target = structures[1]

# ╔═╡ 221be9ab-7edf-44bc-943f-f643a477aed3
bmodel = BindingModel{Int8(3)}(
	ligands=[
		Ligand(;beads=[3, 3, 3, 3, 3], shape=:cross), 
		Ligand(;beads=[9, 9, 9, 9, 9], shape=:cross) 
	],
	target=target,
	specificities=[0., 0., 0.],
	μ=15,
)

# ╔═╡ 2c7c19a5-3aea-4175-a39f-c5fdc86f54e4
LP.binding_energies(bmodel, sequence)

# ╔═╡ f85d831e-fb85-4300-80aa-59fac0160901
LP.binding_probabilities(bmodel, sequence)

# ╔═╡ 0b4f5317-2dca-4303-9e71-ffa58fb2af6e
default_target = 1

# ╔═╡ 15d39597-8a55-4175-93e6-ec3491640836
n_sequences_per_chain = 50

# ╔═╡ 3bbff889-6177-4953-8e4f-470d5ec1da2d
n_steps_between_sequences = 25

# ╔═╡ c53ae14b-f55d-4b5b-8f96-29d303e24220
β = 100

# ╔═╡ 05d7252f-27cf-4b86-86bc-b51c9cfd4852
binding_model = BindingModel{Int8(3)}(
	ligands=[
		Ligand(;beads=[3, 3, 3, 3, 3], shape=:cross), 
		# Ligand(;beads=[9, 9, 9, 9, 9], shape=:cross) 
	],
	target=target,
	specificities=[1., 0.],
	μ=20,
)

# ╔═╡ ac216d9b-fbe5-4891-86b3-6157d325a948
chains, metrics = let
	parameters = LP.MCMCParameters(;
		n_steps=10,
		n_sequences=100,
		structures,
		target=default_target,
		bmodel=binding_model,
		burnin=0.,
		β_sampling=β,
		JTT_bias=false,
	)
	chain, metrics = LP.sample_mcmc_chain(rand(1:20, 27), parameters)
end

# ╔═╡ 42906e7a-e224-4339-8d4a-6928458115ca
let
	plot([LP.fold_prob(s, structures, default_target) for s in chains])
end

# ╔═╡ de6d7750-891b-4d0e-a880-ee017a9613bd


# ╔═╡ a9ef78d5-48e6-496f-827b-b693764332ba
let
	[LP.binding_probabilities(binding_model, s) for s in chains]
end

# ╔═╡ 56913cde-e107-4914-ba8d-40fa35efd6aa
[LP.compute_ligand_phi(binding_model, s) for s in chains]

# ╔═╡ Cell order:
# ╠═5cd13410-64b0-11f1-81f3-713fcda94f0a
# ╠═62f1f8f2-a51b-4d3a-b82a-c1b1b6a6a710
# ╠═53ec092a-fb32-421c-aa7f-550acec1e6c2
# ╠═ccce05dd-77ad-40d0-b1ff-256b03ac1fce
# ╟─195dfa45-af1d-4b42-911a-7d3d2f50835a
# ╠═6eeaaea4-e20d-4e67-a144-b304dfac607f
# ╠═38bb847d-14eb-4d77-a428-fae8071b0005
# ╠═1ca78bc1-6281-4bb1-bc1c-7cac9efd0171
# ╟─430d9438-ad7f-44fd-87a7-2d99e725a98e
# ╟─aabb7efe-aed6-46b2-9bfb-9a15b3527349
# ╠═66384e93-a183-4064-986a-a716a26d6fb4
# ╠═93434855-a17b-48f7-a9c2-7264a55da687
# ╠═01b8d561-363b-46ff-b4f6-ce022fd03172
# ╠═26a17f9a-8330-4a54-9b07-d2554d846ab5
# ╠═9ed363e9-72ea-463f-bb64-d9f26b46826b
# ╠═abbc653d-1ee1-4eb2-8ad6-26d11dadb6ad
# ╠═30d9d15f-58e4-4bd0-ac3a-61eb8c87a253
# ╟─e29ccd3f-3b47-4c10-becb-960f1f60d996
# ╠═3c8de6d6-21c4-4fe7-935b-d6411704d6e2
# ╟─7f264679-e8a2-4bed-8b7e-7af933471a8f
# ╠═9a727f50-efd4-48bf-bc26-d3b3fd8f8878
# ╠═221be9ab-7edf-44bc-943f-f643a477aed3
# ╠═5eec4f20-b69f-443d-af1e-a3f2711b8111
# ╠═2c7c19a5-3aea-4175-a39f-c5fdc86f54e4
# ╠═f85d831e-fb85-4300-80aa-59fac0160901
# ╠═5a45302d-fa42-48a8-81f5-c7ec614f586d
# ╠═2723c3b5-0425-4830-aae0-e89151fda49e
# ╠═33dbed1f-9ea1-4234-8279-4ae68da8772c
# ╠═ac0353d0-3efc-4d71-aa31-362215821995
# ╠═69e781be-16bd-42f5-8ee2-526c6dd918a0
# ╠═104c0512-c26f-43cc-8724-e5a27c2a8f4b
# ╠═1d77d89a-97ec-4ab2-9add-13672dd0729a
# ╠═4858a6b2-dd8b-49d1-a5ef-9458e8e9478f
# ╠═5f969c85-e8ee-4355-b63c-c7a68b2f2abe
# ╠═93c8f931-eaa7-4486-a6b6-f8bab43466d0
# ╠═8689f611-15e3-4b55-b980-877cccf3ec0d
# ╠═0b4f5317-2dca-4303-9e71-ffa58fb2af6e
# ╠═15d39597-8a55-4175-93e6-ec3491640836
# ╠═3bbff889-6177-4953-8e4f-470d5ec1da2d
# ╠═c53ae14b-f55d-4b5b-8f96-29d303e24220
# ╠═ac216d9b-fbe5-4891-86b3-6157d325a948
# ╠═42906e7a-e224-4339-8d4a-6928458115ca
# ╠═05d7252f-27cf-4b86-86bc-b51c9cfd4852
# ╠═de6d7750-891b-4d0e-a880-ee017a9613bd
# ╠═a9ef78d5-48e6-496f-827b-b693764332ba
# ╠═56913cde-e107-4914-ba8d-40fa35efd6aa
