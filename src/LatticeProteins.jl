module LatticeProteins

using ArgCheck
using DelimitedFiles
using LinearAlgebra
using ProgressMeter
using Random
using StatsBase
using TreeTools
using UnPack

include("Structures.jl")
using .Structures
export Structure

include("recipes.jl")
include("utils.jl")

include("ligand.jl")
export Ligand, BindingModel

include("metropolis.jl")
export MCMCParameters
export sample_mcmc_chain
export FIELDS

include("tree_evolution.jl")
export evolve_on_tree

end # module LatticeProteins
