const MJ_1996 = readdlm(joinpath(pkgdir(LatticeProteins), "data/MJ_1996.csv"), Float64) # 20x20 matrix

const MJ_alphabet = [
    "C",
    "M",
    "F",
    "I",
    "L",
    "V",
    "W",
    "Y",
    "A",
    "G",
    "T",
    "S",
    "N",
    "Q",
    "D",
    "E",
    "H",
    "R",
    "K",
    "P",
]
const MJ_alphabet_3 = [
    "Cys",
    "Met",
    "Phe",
    "Ile",
    "Leu",
    "Val",
    "Trp",
    "Tyr",
    "Ala",
    "Gly",
    "Thr",
    "Ser",
    "Asn",
    "Gln",
    "Asp",
    "Glu",
    "His",
    "Arg",
    "Lys",
    "Pro",
]

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

    structures::Vector{<:Structure}
    target::Int = 1 # id of target structure
end

"""
Buffers and mutable states of the MCMC chain
"""
mutable struct MCMCState
    sequence::Vector{Int}
    energies::Vector{Float64}
    energy_buffer::Vector{Float64} # compute energies for proposed flip
    delta_mj::Vector{Float64} # store relevant rows of MJ for a given sequence
    ϕ_old::Float64 # old log_fold_prob value to avoid one softmin
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
    sample_mcmc_chain(
        init::Vector{Int}, parameters::MCMCParameters; rng=Random.GLOBAL_RNG
    )
"""
function sample_mcmc_chain(
    init::Vector{Int}, parameters::MCMCParameters; rng=Random.GLOBAL_RNG
)
    contact_tensor = build_contact_tensor(parameters.structures)
    energies = map(S -> energy(init, S), parameters.structures)
    state = MCMCState(
        copy(init),
        energies,
        similar(energies),
        zeros(length(init)),
        _softmin(energies, parameters.target),
    )

    for _ in 1:(parameters.burnin)
        metropolis_swap!(
            state, contact_tensor, parameters.target, parameters.β_sampling; rng
        )
    end
    chains = Vector{Vector{Int}}(undef, parameters.n_sequences)
    metrics = []
    chains[1] = copy(state.sequence)
    for s in 2:(parameters.n_sequences)
        for _ in 1:(parameters.n_steps)
            accepted, i, a = metropolis_swap!(
                state, contact_tensor, parameters.target, parameters.β_sampling; rng
            )
            push!(metrics, (; accepted, i, a))
        end
        chains[s] = copy(state.sequence)
    end

    return chains, metrics
end

function metropolis_swap!(
    state::MCMCState,
    contact_tensor::Array{Int8,3},
    target::Int,
    β::Float64;
    rng=Random.GLOBAL_RNG,
)
    (; sequence, energies, energy_buffer, delta_mj) = state
    i = rand(rng, 1:length(sequence))
    a = sequence[i]
    b = mod(rand(rng, 1:19) + a, 20) + 1
    @inbounds for j in 1:length(sequence)
        delta_mj[j] = MJ_1996[b, sequence[j]] - MJ_1996[a, sequence[j]]
    end
    K, n_structures, _ = size(contact_tensor)
    partners = view(contact_tensor, :, :, i)  # (K, n_structures), contiguous
    for n in 1:n_structures
        ΔE = 0.0
        for k in 1:K
            j = partners[k, n]
            j == 0 && break
            ΔE += delta_mj[j]
        end
        energy_buffer[n] = ΔE
    end
    energy_buffer .+= energies

    ϕ_new = _softmin(energy_buffer, target)

    if ϕ_new > state.ϕ_old || rand(rng) < (ϕ_new / state.ϕ_old)^β
        state.ϕ_old = ϕ_new
        energies .= energy_buffer
        sequence[i] = b
        return true, i, b
    else
        return false, i, a
    end
end

function _softmin!(X::Vector{Float64})
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

function energy(sequence::Vector{Int}, structure::Structure)
    energy = 0.0
    for (i, j) in structure.contacts
        energy += MJ_1996[sequence[i], sequence[j]]
    end
    return energy
end
function energy(sequences::Vector{Vector{Int}}, structure::Structure)
    out = zeros(length(sequences))
    for (i, j) in structure.contacts
        for (n, seq) in enumerate(sequences)
            out[n] += MJ_1996[seq[i], seq[j]]
        end
    end
    return out
end

function fold_prob(sequence::Vector{Int}, structures)
    E = map(S -> energy(sequence, S), structures)
    _softmin!(E)
    return E
end
function fold_prob(sequence::Vector{Int}, structures, target::Integer)
    return fold_prob(sequence, structures)[target]
end
function fold_prob(sequences::Vector{Vector{Int}}, structures)
    n_str, n_seq = length(structures), length(sequences)
    E = zeros(n_str, n_seq)  # E[m, n] = energy of sequence n in structure m
    for (m, S) in enumerate(structures)
        E[m, :] .= energy(sequences, S)
    end
    for n in 1:n_seq
        _softmin!(view(E, :, n))  # contiguous column view
    end
    return E  # E[m, n] = P(structure m | sequence n)
end
function fold_prob(sequences::Vector{Vector{Int}}, structures, target::Integer)
    return fold_prob(sequences, structures)[target, :]
end

### FOR TESTING PURPOSES

function metropolis_swap_naive!(
    sequence::Vector{Int},
    energies::Vector{Float64},
    structures::Vector{Structure},
    target::Int,
    β::Float64;
    rng=Random.GLOBAL_RNG,
)
    i = rand(rng, 1:length(sequence))
    a = sequence[i]
    b = mod(rand(rng, 1:19) + a, 20) + 1

    ϕ_old = log(fold_prob(sequence, structures, target))
    S = copy(sequence)
    S[i] = b
    ϕ_new = log(fold_prob(S, structures, target))

    if ϕ_new > ϕ_old || rand(rng) < exp(β * (ϕ_new - ϕ_old))
        # energies .= new_energies
        sequence[i] = b
        return true, i, b
    else
        return false, i, a
    end
end

function _metropolis_swap!(
    state::MCMCState,
    structures::Vector{Structure{N}},
    target::Int,
    β::Float64;
    rng=Random.GLOBAL_RNG,
) where {N}
    (; sequence, energies, energy_buffer, delta_mj) = state
    i = rand(rng, 1:length(sequence))
    a = sequence[i]
    b = mod(rand(rng, 1:19) + a, 20) + 1
    @inbounds for j in 1:length(sequence)
        delta_mj[j] = MJ_1996[b, sequence[j]] - MJ_1996[a, sequence[j]]
    end
    for (n, S) in enumerate(structures)
        ΔE = 0.0
        for j in S.contact_partners[i]
            ΔE += delta_mj[j]
        end
        energy_buffer[n] = ΔE
    end
    energy_buffer .+= energies

    ϕ_new = _softmin(energy_buffer, target)

    if ϕ_new > state.ϕ_old || rand(rng) < (ϕ_new / state.ϕ_old)^β
        state.ϕ_old = ϕ_new
        energies .= energy_buffer
        sequence[i] = b
        return true, i, b
    else
        return false, i, a
    end
end
