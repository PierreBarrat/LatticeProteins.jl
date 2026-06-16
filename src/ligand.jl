using ArgCheck # repeated here for standalone include()

@kwdef struct Ligand
    beads::Vector{Int}
    shape::Symbol
    # cells: [(0, 0), (1, 0), (0, 1)] would be an `L` ligand.
    cells::Vector{Tuple{Int,Int}} = cells_from_shape(shape)
    reflection::Bool = true
    rotation::Bool = true
    function Ligand(beads::AbstractVector, shape::Symbol, cells, reflection, rotation)
        @argcheck length(beads) == length(cells)
        @argcheck cells[1] == (0, 0)
        return new(beads, shape, cells, reflection, rotation)
    end
end

function Base.:(==)(l1::Ligand, l2::Ligand)
    return l1.beads == l2.beads &&
           l1.shape == l2.shape &&
           l1.reflection == l2.reflection &&
           l1.rotation == l2.rotation
end

cross_cells() = [(0, 0), (0, 1), (1, 0), (0, -1), (-1, 0)]
rod3_cells() = [(0, 0), (-1, 0), (1, 0)]

function cells_from_shape(shape)
    return if shape == :cross
        cross_cells()
    elseif shape == :rod3
        rod3_cells()
    else
        allowed_shapes = [:cross]
        error("Shape $shape not recognized. Allowed shapes are: $(allowed_shapes).")
    end
end

############################################################################################
#= BINDING MODES
For a given ligand, generate all binding patches as vectors [(x,y,z)] of length `cells`.
(x,y,z) refers to a position on the cube: this is INDEPENDENT of structure.
In a simple case, there are 6 x 4 x 2 patches (cube faces, rotations, reflection). There could be more as we will also consider translation of the ligand (if N=4 for instance).
This set has to be generated ONCE per ligand
=#
############################################################################################

struct Face
    axis::Int # axis orthogonal to the face (1 -> x, 2 -> y, 3 -> z)
    side::Int # side of the cube: 1 or N
end
faces(N) = Face[Face(1, 1), Face(1, N), Face(2, 1), Face(2, N), Face(3, 1), Face(3, N)]

"""
    embed(r2::Tuple{Int, Int}, face::Face)

Embed 2D coordinates `r2` of the form `(u, v)` into a 3D site that lives on `face`.

### Example
`face=Face(1, 1)`, `r2=(1, 2)`: output `(1, 1, 2)`

or `face=Face(2, 4)`, `r2=(1, 2)`: output `(1, 4, 2)` (far face axis to y for N=4)
"""
function embed(r2::Tuple{Int,Int}, face::Face)
    r3 = Int[0, 0, 0]
    r3[face.axis] = face.side
    u, v = filter(!=(face.axis), [1, 2, 3]) # inplace axes
    r3[u] = r2[1]
    r3[v] = r2[2]
    return Structures.site(r3...)
end
embed(patch::Vector{Tuple{Int,Int}}, face::Face) = [embed(r2, face) for r2 in patch]

function symmetries(ligand::Ligand)
    # Always rotate around first cell
    rotations = ligand.rotation ? 4 : 1
    reflections = ligand.reflection ? 2 : 1
    return [(r, f) for r in 1:rotations for f in 1:reflections]
end
function rotate(cell::Tuple{Int,Int}, r)
    @argcheck r in 1:4
    u, v = cell
    return if r == 1
        (u, v)
    elseif r == 2
        (-v, u)
    elseif r == 3
        (-u, -v)
    elseif r == 4
        (v, -u)
    end
end
# flip across y-axis (arbitrary choice)
flip(cell::Tuple{Int,Int}, f) = f == 1 ? cell : (-cell[1], cell[2])

"""

Return a vector of binding patches `bps` and a vector of faces.
A binding patch is a vector `bp = [(x1, y1, z1), ...]`.
`bp[i]=(x,y,z)` means that bead `i` of the ligand touches site `(x,y,z)` of the cube.

The vector `F` of faces is such that `bps[i]` is a patch corresponding to `F[i]`.
"""
function binding_patches(ligand::Ligand, N::Integer)
    patches_2D = binding_patches_2D(ligand, N) # for an arbitrary face
    patches = Vector{Vector{Structures.Site}}(undef, 0)
    F = Face[] # storing binding faces for potential analysis
    for f in faces(N)
        # embed each (u, v) 2D patch onto the face, returning (x, y, z)
        append!(patches, [embed(p, f) for p in patches_2D])
        append!(F, repeat([f], length(patches_2D)))
    end
    return patches, F
end
function binding_patches_2D(ligand, N)
    # patches for an arbitrary face, as a 2D grid
    patches = []
    for (r, f) in symmetries(ligand)
        cells = [rotate(flip(c, f), r) for c in ligand.cells]
        # Try to put first cell (0, 0) at any position (u0, v0) in (1:N, 1:N)
        for u0 in 1:N, v0 in 1:N
            patched_cells = [(c[1] + u0, c[2] + v0) for c in cells]
            # if the ligand fits, valid patch
            fits = all(patched_cells) do (u, v)
                u in 1:N && v in 1:N
            end
            if fits
                push!(patches, patched_cells)
            end
        end
    end
    return patches
end

############################################################################################
### BINDING
############################################################################################

struct BindingModel{N}
    ligands::Vector{Ligand}
    target::Structure{N}
    specificities::Vector{Float64} # one per ligand + unbound
    μ::Float64 # chemical potential
    sim::Vector{Int} # site_index_map
    patches::Vector{Vector{Vector{Structures.Site}}} # nested N_ligands x N_patches x patch

    # μ = Inf --> complete competition between ligands (each has infinite concentration)
    # μ = -Inf (numerically bad but ...) binding probabilities are independent
    # as the unbound state gets almost all the mass.

    function BindingModel{N}(
        ligands, target::Structure, specificities, μ, sim, patches
    ) where {N}
        @argcheck length(specificities) == length(ligands) + 1 """Expected one specificity
        for each ligand plus one for unbound state.
        Instead $(length(ligands)) ligands and $(length(specificities)) specificities.
        """

        return new{N}(ligands, target, specificities, Float64(μ), sim, patches)
    end
end

"""
    BindingModel(; ligands, target, specificities, μ)
"""
function BindingModel(;
    ligands::Vector{Ligand}, target::Structure{N}, specificities::Vector{Float64}, μ::Real
) where {N}
    patches = [binding_patches(ligand, N)[1] for ligand in ligands]
    sim = Structures.site_index_map(target)
    return BindingModel{N}(ligands, target, specificities, μ, sim, patches)
end

"""
    binding_energies(ligand::Ligand, structure::Structure, sequence::Vector{Int})

Return all binding energies for this ligand (all faces & orientations), as well as vector with the faces the binding occurs on.
"""
function binding_energies(
    ligand::Ligand,
    patches::Vector{Vector{Structures.Site}},
    structure::Structure{N},
    sequence::Vector{Int};
    site_index_map=Structures.site_index_map(structure),
) where {N}
    site_id = LatticeProteins.Structures._site_id
    energies = zeros(Float64, length(patches))
    for (p, patch) in enumerate(patches)
        E = 0.0
        for (i_ligand, site) in enumerate(patch)
            i_chain = site_index_map[site_id(site, N)]
            aa = sequence[i_chain]
            bead = ligand.beads[i_ligand]
            E += MJ_1996[aa, bead]
        end
        energies[p] = E
    end
    return energies
end
# Version that computes binding patches
function binding_energies(
    ligand::Ligand, structure::Structure{N}, sequence::Vector{Int}; kwargs...
) where {N}
    patches, faces = binding_patches(ligand, N)
    energies = binding_energies(ligand, patches, structure, sequence; kwargs...)
    return energies
end
"""
    binding_energies(model::BindingModel, sequence::Vector{Int})

For each ligand in the model, return (i) binding energies of all configurations (ii) corresponding faces.
"""
function binding_energies(model::BindingModel, sequence::Vector{Int})
    return map(zip(model.ligands, model.patches)) do (ligand, patches)
        binding_energies(ligand, patches, model.target, sequence; site_index_map=model.sim)
    end
end
"""
    binding_probabilities(model::BindingModel, sequence::Vector{Int})

Return binding probabilities of all ligands in `model` to `sequence`. Ligands are in competition for binding. Last element of return array is the probability of the unbound state.
"""
function binding_probabilities(model::BindingModel, sequence::Vector{Int})
    energies = [x for x in binding_energies(model, sequence)]
    push!(energies, [model.μ]) # unbound state
    # μ ~ energy cost of being in the solution

    # Softmin - Zs contains ligand-specific partition functions
    E_min = minimum(Iterators.flatten(energies))
    Zs = map(energies) do E # E ~ energies of given ligand for all faces / orientations
        sum(x -> exp(-x + E_min), E)
    end
    Z_tot = sum(Zs)
    return Zs / Z_tot
end

function compute_ligand_phi(model::BindingModel, sequence::Vector{Int})
    ϕs = map(log, binding_probabilities(model, sequence))
    return sum(λ * ϕ for (λ, ϕ) in zip(model.specificities, ϕs))
end
compute_ligand_phi(::Nothing, sequence::Vector{Int}) = 0.0
