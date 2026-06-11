const FIELDS = Dict{Symbol,Matrix{Float64}}()

function _resolve_field(field::Symbol, L::Int)
    field == :none && return nothing
    h = FIELDS[field]
    return size(h, 1) == 1 ? repeat(h, L, 1) : h
end

# const MJ_1996 = rand(Float64, 20, 20)
"""
    MCMCParameters

Parameters for controlling the Metropolis-Hastings MCMC simulation.

# Fields
- `n_steps::Int`: Number of MCMC steps between each saved sequence
- `n_sequences::Int`: Number of sequences to sample
- `burnin::Int`: Number of burn-in steps to discard before sampling (default: `10 * n_steps`)
- `β_sampling::Float64`: Inverse temperature parameter for sampling (default: `100.0`)
- `β_thermo::Float64`: Thermodynamic inverse temperature (default: `1.0`)
- `structures::Vector{Structure}`: Collection of protein structures used for energy calculations
- `target::Int`: Index of the target structure in `structures` (default: `1`)
"""
Base.@kwdef struct MCMCParameters
    n_steps::Int # number of steps between sequences
    n_sequences::Int # number of sequences the chain should sample
    burnin::Int = 10 * n_steps
    β_sampling::Float64 = 100.0 # sampling temperature
    β_thermo::Float64 = 1.0 # thermodynamic temperature # to remove, won't be used in MCMC

    JTT_bias::Bool = false # whether to pick mutation according to the JTT matrix

    structures::Vector{<:Structure}
    target::Int = 1 # id of target structure
    field::Symbol = :none

    bmodel::Union{Nothing,BindingModel} = nothing
end

"""
    MCMCState

Mutable runtime state of the Metropolis chain. Holds the current sequence, per-structure
energies, scratch buffers reused each step, and the cached fold probability `ϕ_old`.
"""
mutable struct MCMCState
    sequence::Vector{Int}
    energies::Vector{Float64}
    energy_buffer::Vector{Float64} # compute energies for proposed flip
    delta_mj::Vector{Float64} # store relevant rows of MJ for a given sequence
    ϕ_old::Float64 # old log fitness value to avoid one softmin for structures
end

"""
    copy(state::MCMCState) -> MCMCState

Deep-copy the real state (`sequence`, `energies`, `ϕ_old`); allocate fresh scratch buffers.
Used to branch the chain at each tree node.
"""
function Base.copy(state::MCMCState)
    return MCMCState(
        copy(state.sequence),
        copy(state.energies),
        similar(state.energy_buffer),
        similar(state.delta_mj),
        state.ϕ_old,
    )
end

# K x N x L representation of contacts for quick access in core mcmc loop
# K ~ max number of neighbours for a given site
# N ~ number of structures
function build_contact_tensor(structures::Vector{<:Structure{N}}) where {N}
    L = Int(N)^3
    n_structures = length(structures)
    K = maximum(maximum(length(s.contact_partners[i]) for i in 1:L) for s in structures)
    tensor = zeros(Int8, K, n_structures, L)
    for (n, s) in enumerate(structures)
        for i in 1:L
            for (k, j) in enumerate(s.contact_partners[i])
                tensor[k, n, i] = Int8(j)
            end
        end
    end
    return tensor
end

"""
    _mcmc_init(init, parameters) -> (MCMCState, contact_tensor, field)

Build the precomputed structures needed for MCMC: the contact tensor (K×N_structures×L),
the resolved field matrix (or `nothing`), and an initialised `MCMCState` for `init`.
Intended to be called once and shared across all steps or tree branches.
"""
function _mcmc_init(init::Vector{Int}, parameters::MCMCParameters)
    contact_tensor = build_contact_tensor(parameters.structures)
    field = _resolve_field(parameters.field, length(init))
    energies = map(S -> energy(init, S), parameters.structures)
    state = MCMCState(
        copy(init),
        energies,
        similar(energies),
        zeros(length(init)),
        log(_softmin(energies, parameters.target)),
    )
    return state, contact_tensor, field
end

"""
    _advance_mcmc!(state, contact_tensor, field, target, β, n_steps; JTT_bias, rng)

Run `n_steps` Metropolis steps in-place on `state`. No output — mutates `state` only.
"""
function _advance_mcmc!(
    state::MCMCState,
    contact_tensor::Array{Int8,3},
    field,
    target::Int,
    β::Float64,
    n_steps::Int;
    JTT_bias=false,
    rng=Random.GLOBAL_RNG,
)
    for _ in 1:n_steps
        metropolis_swap!(state, contact_tensor, field, target, β; JTT_bias, rng)
    end
    return nothing
end

"""
    sample_mcmc_chain(init::Vector{Int}, parameters::MCMCParameters; rng, progress) -> (chains, metrics)

Run a Metropolis-Hastings chain starting from `init` and return sampled sequences.

`parameters.burnin` steps are discarded first. Then `parameters.n_sequences` sequences are
collected, with `parameters.n_steps` steps between each saved sequence.

Returns:
- `chains`: vector of `n_sequences` sequences (each a `Vector{Int}` of length L=N³)
- `metrics`: vector of `(; accepted, i, a)` named tuples, one per step after burnin

Sequences are integers in 1:20 using the alphabet in `ref_alphabet` (alphabetical order).
Set `progress=true` to display a progress bar over the `n_sequences` sequences.
"""
function sample_mcmc_chain(
    init::Vector{Int}, parameters::MCMCParameters; rng=Random.GLOBAL_RNG, progress=false
)
    @unpack target, β_sampling, JTT_bias = parameters
    state, contact_tensor, field = _mcmc_init(init, parameters)

    _advance_mcmc!(
        state, contact_tensor, field, target, β_sampling, parameters.burnin; JTT_bias, rng
    )

    chains = Vector{Vector{Int}}(undef, parameters.n_sequences)
    metrics = []
    chains[1] = copy(state.sequence)
    prog = Progress(parameters.n_sequences; enabled=progress, desc="Sampling sequences")
    next!(prog)
    for s in 2:(parameters.n_sequences)
        for _ in 1:(parameters.n_steps)
            accepted, i, a = metropolis_swap!(
                state, contact_tensor, field, target, β_sampling; JTT_bias, rng
            )
            push!(metrics, (; accepted, i, a))
        end
        chains[s] = copy(state.sequence)
        next!(prog)
    end

    return chains, metrics
end

"""
    metropolis_swap!(state, contact_tensor, field, target, β; JTT_bias, rng)
        -> (accepted::Bool, position::Int, new_aa::Int)

Propose and conditionally accept one single-site mutation. Mutates `state` in-place.
Returns whether the move was accepted, the mutated position, and the resulting amino acid.
"""
function metropolis_swap!(
    state::MCMCState,
    contact_tensor::Array{Int8,3},
    field::Union{Nothing,Matrix{Float64}},
    target::Int,
    β::Float64;
    JTT_bias=false,
    rng=Random.GLOBAL_RNG,
)
    # choose mutation
    i, a, b = pick_mutation(rng, state.sequence, JTT_bias)

    # change in log Pfold
    ϕ_fold_new = compute_log_pfold!(state, contact_tensor, target, (i, a, b))

    # Combine fitnesses
    ϕ_new = ϕ_fold_new

    # potentially add bias
    r = if isnothing(field)
        β * (ϕ_new - state.ϕ_old)
    else
        β * (ϕ_new - state.ϕ_old) - β * (field[i, b] - field[i, a])
    end

    if r > 0 || rand(rng) < exp(r)
        state.ϕ_old = ϕ_new
        state.energies .= state.energy_buffer
        state.sequence[i] = b
        return true, i, b
    else
        return false, i, a
    end
end

function compute_log_pfold!( # ! because buffers in `state` are changed
    state::MCMCState,
    contact_tensor::Array{Int8,3},
    target::Int,
    mutation::Tuple{Int,Int,Int}, # i a b
)
    # unpack state variables, buffer, and mutation
    (; sequence, energies, energy_buffer, delta_mj) = state
    i, a, b = mutation

    # Precompute ΔE_ij if i: a --> b
    @inbounds for j in eachindex(sequence)
        delta_mj[j] = MJ_1996[b, sequence[j]] - MJ_1996[a, sequence[j]]
    end

    # structure loop
    K, n_structures, _ = size(contact_tensor)
    partners = view(contact_tensor, :, :, i)  # (K, n_structures), contiguous
    for n in 1:n_structures
        ΔE = 0.0
        for k in 1:K
            j = partners[k, n] # k-th contact of pos i in structure n
            j == 0 && break
            ΔE += delta_mj[j]
        end
        energy_buffer[n] = ΔE
    end
    energy_buffer .+= energies

    P_fold = _softmin(energy_buffer, target)
    return log(P_fold)
end

"""
    pick_mutation(rng, sequence, JTT_bias) -> (position, old_aa, new_aa)

Sample a random position and a new amino acid. If `JTT_bias`, draw the new state from
the JTT row of the current amino acid; otherwise draw uniformly from the other 19 states.
"""
function pick_mutation(rng, sequence, JTT_bias)
    i = rand(rng, 1:length(sequence))
    a = sequence[i]
    b = if !JTT_bias
        mod(rand(rng, 1:19) + a - 1, 20) + 1
    else
        wsample(rng, 1:20, JTT[a, :])
    end
    return i, a, b
end

function _softmin!(X::AbstractVector)
    min_X = minimum(X)
    for (i, x) in enumerate(X)
        X[i] = exp(-x + min_X)
    end
    X ./= sum(X)
    return X
end
_softmin(X) = _softmin!(copy(X))

function _softmin(X, idx::Integer)
    min_X = minimum(X)
    Z = sum(x -> exp(-x + min_X), X)
    return exp(-X[idx] + min_X) / Z
end

"""
    energy(sequence, structure; field=nothing) -> Float64
    energy(sequences, structure; field=nothing) -> Vector{Float64}

MJ_1996 contact energy of `sequence` (or each sequence in `sequences`) folded into
`structure`. If `field` is provided (L×20 matrix), adds the per-site bias term.
"""
function energy(
    sequence::Vector{Int},
    structure::Structure;
    field::Union{Nothing,Matrix{Float64}}=nothing,
)
    e = 0.0
    for (i, j) in structure.contacts
        e += MJ_1996[sequence[i], sequence[j]]
    end
    if !isnothing(field)
        for i in eachindex(sequence)
            e += field[i, sequence[i]]
        end
    end
    return e
end
function energy(
    sequences::Vector{Vector{Int}},
    structure::Structure;
    field::Union{Nothing,Matrix{Float64}}=nothing,
)
    out = zeros(length(sequences))
    for (i, j) in structure.contacts
        for (n, seq) in enumerate(sequences)
            out[n] += MJ_1996[seq[i], seq[j]]
        end
    end
    if !isnothing(field)
        L = size(field, 1)
        for i in 1:L
            for (n, seq) in enumerate(sequences)
                out[n] += field[i, seq[i]]
            end
        end
    end
    return out
end

"""
    thermostability_curve(sequence, structures, target, βvals; field) -> Vector{Float64}
    thermostability_curve(sequences, structures, target, βvals; field) -> Matrix{Float64}

Fraction-folded curve P(target | β) for each β in `βvals`.
For multiple sequences returns an (n_β × n_sequences) matrix.
"""
function thermostability_curve(
    sequence::Vector{Int},
    structures,
    target::Integer,
    βvals::AbstractVector;
    field::Union{Nothing,Matrix{Float64}}=nothing,
)
    E = map(S -> energy(sequence, S; field), structures)
    return map(β -> _softmin!(β * copy(E))[target], βvals)
end
"""
    thermostability(sequence, structures, target, βvals, threshold; field) -> Float64 or nothing
    thermostability(sequences, structures, target, βvals, threshold; field) -> Vector

Estimate the melting temperature as the β⁻¹ at which the fraction-folded curve crosses
`threshold`. Returns `nothing` if the threshold is never crossed within `βvals`.
"""
function thermostability(
    sequence::Vector{Int},
    structures,
    target::Integer,
    βvals::AbstractVector,
    threshold::AbstractFloat;
    field::Union{Nothing,Matrix{Float64}}=nothing,
)
    frac_folded = thermostability_curve(sequence, structures, target, βvals; field)
    i1 = findlast(>=(threshold), frac_folded)
    i2 = findfirst(<(threshold), frac_folded)
    return isnothing(i1) || isnothing(i2) ? nothing : (1 / βvals[i1] + 1 / βvals[i2]) / 2
end

function thermostability_curve(
    sequences::Vector{Vector{Int}},
    structures,
    target::Integer,
    βvals::AbstractVector;
    field::Union{Nothing,Matrix{Float64}}=nothing,
)
    n_str, n_seq = length(structures), length(sequences)
    E = zeros(n_str, n_seq)  # E[m, n] = energy of sequence n in structure m
    for (m, S) in enumerate(structures)
        E[m, :] .= energy(sequences, S; field)
    end

    return mapreduce(hcat, 1:n_seq) do n # nβ x n: curve of sequences as columns
        map(β -> _softmin!(β * copy(E[:, n]))[target], βvals)
    end
end
function thermostability(
    sequences::Vector{Vector{Int}},
    structures,
    target::Integer,
    βvals::AbstractVector,
    threshold;
    field::Union{Nothing,Matrix{Float64}}=nothing,
)
    frac_folded = thermostability_curve(sequences, structures, target, βvals; field)
    return map(eachcol(frac_folded)) do frac
        i1 = findlast(>=(threshold), frac)
        i2 = findfirst(<(threshold), frac)
        isnothing(i1) || isnothing(i2) ? nothing : (1 / βvals[i1] + 1 / βvals[i2]) / 2
    end
end

"""
    fold_prob(sequence, structures; field=nothing, β=1.0) -> Vector{Float64}
    fold_prob(sequence, structures, target; field=nothing, β=1.0) -> Float64
    fold_prob(sequences, structures; field=nothing, β=1.0) -> Matrix{Float64}
    fold_prob(sequences, structures, target; field=nothing, β=1.0) -> Vector{Float64}

Boltzmann folding probabilities at inverse temperature `β`.
Without `target`: returns the full probability vector over all structures.
With `target`: returns the scalar probability of folding into that structure.
For multiple sequences the matrix form is (n_structures × n_sequences).
"""
function fold_prob(
    sequence::Vector{Int}, structures; field::Union{Nothing,Matrix{Float64}}=nothing, β=1.0
)
    E = map(S -> β * energy(sequence, S; field), structures)
    _softmin!(E)
    return E
end
function fold_prob(
    sequence::Vector{Int},
    structures,
    target::Integer;
    field::Union{Nothing,Matrix{Float64}}=nothing,
    β=1.0,
)
    return fold_prob(sequence, structures; field, β)[target]
end
function fold_prob(
    sequences::Vector{Vector{Int}},
    structures;
    field::Union{Nothing,Matrix{Float64}}=nothing,
    β=1.0,
)
    n_str, n_seq = length(structures), length(sequences)
    E = zeros(n_str, n_seq)  # E[m, n] = energy of sequence n in structure m
    for (m, S) in enumerate(structures)
        E[m, :] .= β * energy(sequences, S; field)
    end
    for n in 1:n_seq
        _softmin!(view(E, :, n))
    end
    return E  # E[m, n] = P(structure m | sequence n)
end
function fold_prob(
    sequences::Vector{Vector{Int}},
    structures,
    target::Integer;
    field::Union{Nothing,Matrix{Float64}}=nothing,
    β=1.0,
)
    return fold_prob(sequences, structures; field, β)[target, :]
end
