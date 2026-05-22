"""
Module for generating and managing Hamiltonian paths (structures) on NxNxN cubic lattices.
"""
module Structures

import Base.isless

using ArgCheck
using ProgressMeter

export Structure, Site, generate_structures, is_neighbour

############################################################################################
# Types
############################################################################################
"""
    Site

A 1-based coordinate (x, y, z) in the lattice, where 1 ≤ x,y,z ≤ N.
Stored as a NamedTuple for clarity: `(x=1, y=2, z=3)`.
"""
const Site = NamedTuple{(:x, :y, :z),NTuple{3,Int8}}
"""
    site(x, y, z) -> Site

Convenience constructor for Site.
"""
site(x, y, z) = (x=Int8(x), y=Int8(y), z=Int8(z))

"""
    Structure{N <: Integer}

A Hamiltonian path through an NxNxN cubic lattice.

# Fields
- `path::Vector{Site}`: Sequence of all N³ sites, each visited exactly once
- `contacts::Vector{Tuple{Int, Int}}`: Path index pairs (i, j) where sites are
  lattice-adjacent but not consecutive in the path
"""
struct Structure{N}
    path::Vector{Site}
    contacts::Vector{Tuple{Int,Int}}

    function Structure{N}(path::Vector{Site}) where {N<:Int8}
        # ---- Validation ----
        @argcheck length(path) == N^3 "Path must contain N³ sites"
        @argcheck allunique(path) "All sites in path must be unique"

        # Check all coordinates are in valid range
        for site in path
            @argcheck 1 ≤ site.x ≤ N && 1 ≤ site.y ≤ N && 1 ≤ site.z ≤ N
        end

        # Check consecutive sites are neighbours
        for i in 1:(length(path) - 1)
            @argcheck is_neighbour(path[i], path[i + 1]) "Consecutive sites must be neighbours"
        end

        return new{N}(path, compute_contacts(path))
    end
end

# Cast any integer type parameter to Int8 before calling the inner constructor
Structure{N}(path::Vector{Site}) where {N<:Integer} = Structure{Int8(N)}(path)
# Convenience constructor (infers N from path)
Structure(path::Vector{Site}) = Structure{_infer_N(path)}(path)

"""
    _infer_N(path) -> Int8

Infer N from a path by finding the maximum coordinate.
"""
_infer_N(path::Vector{Site}) = maximum(max(site.x, site.y, site.z) for site in path)

Base.isless(chain::Structure, other::Structure) = chain.path < other.path

############################################################################################
# Methods
############################################################################################

"""
    is_neighbour(s1::Site, s2::Site) -> Bool

Check if two sites are adjacent (differ by exactly 1 in exactly one coordinate).
"""
function is_neighbour(s1::Site, s2::Site)
    dx = abs(s1.x - s2.x)
    dy = abs(s1.y - s2.y)
    dz = abs(s1.z - s2.z)
    return dx + dy + dz == 1
end

"""
    neighbours(s::Site, N::Int) -> Iterator{Site}

Return an iterator over all lattice-adjacent neighbours of site s.
"""
function neighbours(site::NamedTuple, N::Int)
    x, y, z = site
    Δs = (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)
    neighbours = Iterators.map(Δs) do (Δx, Δy, Δz)
        (x=x + Δx, y=y + Δy, z=z + Δz)
    end
    neighbours = Iterators.filter(neighbours) do (x, y, z)
        1 ≤ x ≤ N && 1 ≤ y ≤ N && 1 ≤ z ≤ N
    end
    return neighbours
end

"""
    compute_contacts(path) -> Vector{Tuple{Int, Int}}

Compute all contact pairs (i, j) where i < j, |i-j| > 1, and path[i] and path[j] are neighbours.
"""
function compute_contacts(path::Vector{Site})
    n = length(path)
    contacts = Tuple{Int,Int}[]
    for i in 1:n, j in (i + 2):n  # Skip consecutive (j = i+1)
        if is_neighbour(path[i], path[j])
            push!(contacts, (i, j))
        end
    end
    return contacts
end

############################################################################################
# Symetries
############################################################################################

function generate_rotation_table(N::Integer) # this hardcodes dimension 3
    table = Vector{Dict{Site,Site}}(undef, 48)

    # permutations of axes
    perms = [[1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1]]

    # sign flips of axes
    signs = [[s1, s2, s3] for s1 in [1, -1], s2 in [1, -1], s3 in [1, -1]]

    # cube center for easy rotation
    x0, y0, z0 = (1 + N) / 2, (1 + N) / 2, (1 + N) / 2

    id = 1
    for p in perms, s in signs
        table[id] = Dict{Site,Site}()
        for x in 1:N, y in 1:N, z in 1:N
            r = [x - x0, y - y0, z - z0]
            r = [r[p[i]] * s[i] for i in 1:3]
            r += [x0, y0, z0]
            table[id][site(x, y, z)] = site(r[1], r[2], r[3])
        end
        id += 1
    end

    return table
end

function rotations(path::Vector{Site}, table)
    return map(table) do rot
        rotated_path = copy(path)
        for (i, site) in enumerate(path)
            rotated_path[i] = rot[site]
        end
        rotated_path
    end
end

function canonical(path::Vector{Site}, rotation_table)
    best_rotated_path = copy(path)
    rotated_path = copy(path)
    # find best symetry
    idx = 1
    for (i, rot) in enumerate(rotation_table)
        for (j, site) in enumerate(path)
            rotated_path[j] = rot[site]
        end
        if rotated_path < best_rotated_path
            best_rotated_path = copy(rotated_path)
            idx = i
        end
    end

    # apply best symetry
    for (i, site) in enumerate(path)
        rotated_path[i] = rotation_table[idx][site]
    end

    return rotated_path
end

############################################################################################
# Generation algorithm
############################################################################################

"""
    build_adj(N::Int) -> Array{Vector{Site}, 3}

Precompute neighbour lists for every site in the NxNxN lattice.
Indexed as `adj[x, y, z]` for zero-allocation lookup during path generation.
"""
function build_adj(N::Int)
    adj = Array{Vector{Site}}(undef, N, N, N)
    for x in Int8(1):Int8(N), y in Int8(1):Int8(N), z in Int8(1):Int8(N)
        adj[x, y, z] = collect(neighbours(site(x, y, z), N))
    end
    return adj
end

function add_canonical!(paths::Set{Vector{Site}}, candidate::Vector{Site}, rotation_table)
    return push!(paths, canonical(candidate, rotation_table))
end
"""
    grow!(path, visited, adj, results)

In-place recursive backtracking: extend `path` with unvisited neighbours until
all sites are visited, then push a copy into `results`.
`adj[x,y,z]` gives the precomputed neighbour list for site (x,y,z).
`visited[x,y,z]` is a plain `Bool` array for O(1) lookup with no hashing.
"""
function grow!(
    path::Vector{Site},
    visited::Array{Bool,3}, # (x, y, z) -> bool
    adj::Array{Vector{Site},3}, # (x, y, z) -> neighbours
    rotation_table::Vector{Dict{Site,Site}}, # rot_id -> Dict(site => rotated_site)
    results::Set{Vector{Site}}, # Set since we'll filter for symetries
)
    if length(path) == length(visited)
        # push!(results, copy(path))
        if length(results) == 4
            # Main.@infiltrate
        end
        add_canonical!(results, path, rotation_table) # path copied internally
        return nothing
    end
    s = path[end]
    for nb in adj[s.x, s.y, s.z]
        visited[nb.x, nb.y, nb.z] && continue
        visited[nb.x, nb.y, nb.z] = true
        push!(path, nb)
        grow!(path, visited, adj, rotation_table, results)
        pop!(path)
        visited[nb.x, nb.y, nb.z] = false
    end
end

const STARTS_3 = [
    # Center of a face or corner (by symetry one gives all)
    site(1, 1, 1), # corner
    site(2, 2, 1), # face center
]
const STARTS_2 = [site(1, 1, 1)] # by symetry it's useless to look at other corners
const STARTS_1 = [site(1, 1, 1)]

const STARTS = (STARTS_1, STARTS_2, STARTS_3)

"""
    generate_all_paths(N::Int) -> Vector{Vector{Site}}

Generate all Hamiltonian paths for an NxNxN lattice using recursive backtracking.
"""
function generate_all_paths(N::Int)
    starting_sites = STARTS[N]
    adj = build_adj(N)
    rotation_table = generate_rotation_table(N)
    all_paths = Set{Vector{Site}}()
    for start in starting_sites
        path = [start]
        visited = zeros(Bool, N, N, N)
        visited[start.x, start.y, start.z] = true
        grow!(path, visited, adj, rotation_table, all_paths)
    end
    return all_paths
end

"""
    generate_structures(N::Int) -> Vector{Structure{N}}

Generate all unique Hamiltonian paths for an NxNxN lattice.
Uses hybrid approach: grow() recursion + rotation canonicalization.
"""
function generate_structures(N::Int)
    # 1. Generate all raw paths using grow()
    return raw_paths = generate_all_paths(N)

    # 2. Canonicalize each path (find minimal rotation)

    # 3. Convert to Structure
    # return [Structure{N}(collect(path)) for path in canonical_paths]
end

# Placeholder for canonicalization
function canonicalize(path::Vector{Site}, N::Int)
    @warn "Not implemented yet"
    return path
end

end # module Structure
