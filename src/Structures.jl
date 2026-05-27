"""
Module for generating and managing Hamiltonian paths (structures) on NxNxN cubic lattices.
"""
module Structures

import Base.isless

using ArgCheck
using JLD2
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
- `path::Vector{Site}`: Sequence of all N^3 sites, each visited exactly once
- `contacts::Vector{Tuple{Int, Int}}`: Path index pairs (i, j) where sites are
  lattice-adjacent but not consecutive in the path
- `contact_partners::Vector{Vector{Int}}`: Adjacency list; `contact_partners[i]`
  lists all j in contact with position i. For O(degree) delta-energy computation.
"""
struct Structure{N}
    path::Vector{Site}
    contacts::Vector{Tuple{Int,Int}}
    contact_partners::Vector{Vector{Int}}

    function Structure{N}(path::Vector{Site}) where {N}
        # ---- Validation ----
        @argcheck N isa Integer && N > 0 "N must be a positive integer, instead $N"
        @argcheck length(path) == N^3 "Path must contain N^3 sites, instead $(length(path))"
        @argcheck allunique(path) "All sites in path must be unique"

        # Check all coordinates are in valid range
        for site in path
            @argcheck 1 ≤ site.x ≤ N && 1 ≤ site.y ≤ N && 1 ≤ site.z ≤ N
        end

        # Check consecutive sites are neighbours
        for i in 1:(length(path) - 1)
            @argcheck is_neighbour(path[i], path[i + 1]) "Consecutive sites must be neighbours"
        end

        contacts = compute_contacts(path)
        return new{Int8(N)}(
            path, contacts, compute_contact_partners(contacts, length(path))
        )
    end
end

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

"""
    compute_contact_partners(contacts, n) -> Vector{Vector{Int}}

Build an adjacency list from a contact list: `result[i]` contains all j such
that (i,j) or (j,i) is a contact. For O(degree) delta-energy updates.
"""
function compute_contact_partners(contacts::Vector{Tuple{Int,Int}}, n::Int)
    partners = [Int[] for _ in 1:n]
    for (i, j) in contacts
        push!(partners[i], j)
        push!(partners[j], i)
    end
    return partners
end

############################################################################################
# Symetries
############################################################################################

# Linear index for a site in the NxNxN lattice (1-based)
_site_id(s::Site, N::Int) = Int((s.x - 1) * N^2 + (s.y - 1) * N + s.z)
function _id_to_site(id::Integer, N::Int)
    return site((id - 1) ÷ N^2 + 1, (id - 1) % N^2 ÷ N + 1, (id - 1) % N + 1)
end

"""
    generate_rotation_table(N) -> Matrix{Int8}  (N^3 x 48)

Each column is one of the 48 cube symmetries. Entry [i, r] is the site ID that
site i maps to under rotation r. Integer indexing — no hashing.
"""
function generate_rotation_table(N::Integer)
    table = Matrix{Int8}(undef, N^3, 48)

    # permutations of axes
    perms = [[1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1]]

    # sign flips of axes
    signs = [[s1, s2, s3] for s1 in [1, -1], s2 in [1, -1], s3 in [1, -1]]

    # cube center for easy rotation
    x0, y0, z0 = (1 + N) / 2, (1 + N) / 2, (1 + N) / 2

    rot_id = 1
    for p in perms, s in signs
        for x in 1:N, y in 1:N, z in 1:N
            r = [x - x0, y - y0, z - z0]
            r = [r[p[i]] * s[i] for i in 1:3]
            r .+= [x0, y0, z0]
            table[_site_id(site(x, y, z), N), rot_id] = _site_id(site(r[1], r[2], r[3]), N)
        end
        rot_id += 1
    end

    return table
end

"""
    canonical(path, rotation_table, N) -> Vector{Site}

Return the lexicographically minimal rotation of `path`.
Works on integer site IDs internally — no hashing.
"""
function canonical(path::Vector{Site}, rotation_table::Matrix{Int8}, N::Int)
    path_ids = Int8[_site_id(s, N) for s in path]
    best = copy(path_ids)
    current = similar(path_ids)

    for r in axes(rotation_table, 2)
        for i in eachindex(path_ids)
            current[i] = rotation_table[path_ids[i], r]
        end
        if current < best
            best .= current
        end
    end

    return [_id_to_site(id, N) for id in best]
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

function add_canonical!(
    paths::Set{Vector{Site}}, candidate::Vector{Site}, rotation_table::Matrix{Int8}, N::Int
)
    return push!(paths, canonical(candidate, rotation_table, N))
end

"""
    grow!(path, visited, adj, rotation_table, N, results)

In-place recursive backtracking: extend `path` with unvisited neighbours until
all sites are visited, then push a copy into `results`.
`adj[x,y,z]` gives the precomputed neighbour list for site (x,y,z).
`visited[x,y,z]` is a plain `Bool` array for O(1) lookup with no hashing.
"""
function grow!(
    path::Vector{Site},
    visited::Array{Bool,3},
    adj::Array{Vector{Site},3},
    rotation_table::Matrix{Int8},
    N::Int,
    results::Set{Vector{Site}},
)
    if length(path) == length(visited)
        add_canonical!(results, path, rotation_table, N)
        return nothing
    end
    s = path[end]
    for nb in adj[s.x, s.y, s.z]
        visited[nb.x, nb.y, nb.z] && continue
        visited[nb.x, nb.y, nb.z] = true
        push!(path, nb)
        grow!(path, visited, adj, rotation_table, N, results)
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
        grow!(path, visited, adj, rotation_table, N, all_paths)
    end
    return all_paths
end

"""
    generate_structures(N::Int) -> Vector{Structure{N}}

Generate all unique Hamiltonian paths for an NxNxN lattice.
Uses grow() recursion + rotation canonicalization to deduplicate.
"""
function generate_structures(N::Int)
    paths = generate_all_paths(N)
    return [Structure{N}(path) for path in paths]
end

function generate_and_save_structures(
    N::Int;
    filtering=[0.01, 0.1, 1.0],
    filename=x -> "data/structures_N$(N)_filtering$(x).jld2",
)
    all_structures = generate_structures(N)
    println("Generated $(length(all_structures)) unique structures for N=$N")

    for f in filtering
        n_save = Int(round(Int, length(all_structures) * f))
        idx = map(x -> round(Int, x), range(1, length(all_structures); length=n_save))
        structures = all_structures[idx]
        JLD2.@save filename(f) structures
    end
end

end # module Structure
