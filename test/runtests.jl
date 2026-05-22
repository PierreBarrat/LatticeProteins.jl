using Test
using LatticeProteins

@testset "LatticeProteins.jl" begin
    include("Structures.jl")
end

@testset "Aqua.jl" begin
    include("Aqua.jl")
end
