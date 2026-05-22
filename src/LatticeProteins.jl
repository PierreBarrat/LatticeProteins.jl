module LatticeProteins

using ArgCheck

include("Structures.jl")
using .Structures
export Structure

include("recipes.jl")

end # module LatticeProteins
