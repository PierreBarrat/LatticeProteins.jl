# General imports
using Revise
using StaticArrays
using DataStructures
using CSV
using DataFrames
using Random
using ProgressLogging
using ProgressMeter
using Plots
using Statistics
using TreeTools
import PlotlyJS #not "using" avoiding keeping PlotlyJS plot function


println("Imports done, importing hamiltonian path")
includet("hamiltonian_path.jl")

println("hamiltonian path imported, importing metropolis")
includet("metropolis.jl")

println("metropolis imported, importing benchmark")
includet("benchmark.jl")

println("benchmark imported, importing visualization")
includet("visualization.jl")

println("visualization imported, importing tree")
includet("tree.jl")

println("tree imported, importing lattice_protein")