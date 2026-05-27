const MJ_1996 = readdlm(joinpath(pkgdir(LatticeProteins), "data/MJ_1996.csv"), Float64) # 20x20 matrix

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
    β_thermo::Float64 = 1.0 # thermodynamic temperature

    structures::Vector{<:Structure}
    target::Int = 1 # id of target structure
end

function energy(sequence::Vector{Int}, structure::Structure)
    energy = 0.0
    for (i, j) in structure.contacts
        energy += MJ_1996[sequence[i], sequence[j]]
    end
    return energy
end
function fold_prob(sequence::Vector{Int}, structures)
    E = map(S -> energy(sequence, S), structures)
    _softmin!(E)
    return E
end
function fold_prob(sequence::Vector{Int}, structures, target::Integer)
    return fold_prob(sequence, structures)[target]
end

function sample_mcmc_chain(
    init::Vector{Int}, parameters::MCMCParameters; rng=Random.GLOBAL_RNG
)
    sequence = copy(init)
    energies = map(S -> energy(sequence, S), parameters.structures)
    energy_buffer = similar(energies)

    for _ in 1:(parameters.burnin)
        metropolis_swap!(
            sequence,
            energies,
            energy_buffer,
            parameters.structures,
            parameters.target,
            parameters.β_sampling;
            rng,
        )
    end
    chains = Vector{Vector{Float64}}(undef, parameters.n_sequences)
    metrics = []
    chains[1] = copy(sequence)
    for s in 2:(parameters.n_sequences)
        for _ in 1:(parameters.n_steps)
            accepted, i, a = metropolis_swap!(
                sequence,
                energies,
                energy_buffer,
                parameters.structures,
                parameters.target,
                parameters.β_sampling;
                rng,
            )
            push!(metrics, (; accepted, i, a))
        end
        chains[s] = copy(sequence)
    end

    return chains, metrics
end

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

function metropolis_swap!(
    sequence::Vector{Int},
    energies::Vector{Float64},
    energy_buffer::Vector{Float64},
    structures::Vector{Structure{N}},
    target::Int,
    β::Float64;
    rng=Random.GLOBAL_RNG,
) where {N}
    i = rand(rng, 1:length(sequence))
    a = sequence[i]
    b = mod(rand(rng, 1:19) + a, 20) + 1
    # delta_mj = view(MJ_1996, b, :) .- view(MJ_1996, a, :)
    delta_mj = [MJ_1996[b, sequence[j]] - MJ_1996[a, sequence[j]] for j in 1:(Int(N)^3)]
    for (n, S) in enumerate(structures)
        ΔE = 0.0
        for j in S.contact_partners[i]
            ΔE += delta_mj[j]
        end
        energy_buffer[n] = ΔE
    end
    energy_buffer .+= energies

    ϕ_new = _softmin(energy_buffer, target) # probability
    ϕ_old = _softmin(energies, target) # probability

    if ϕ_new > ϕ_old || rand(rng) < (ϕ_new / ϕ_old)^β
        energies .= energy_buffer
        sequence[i] = b
        return true, i, b
    else
        return false, i, a
    end
end

# Version that does not update energies unless needed (avoids one alloc)
function _metropolis_swap!(
    sequence::Vector{Int},
    energies::Vector{Float64},
    structures::Vector{Structure},
    target::Int;
    rng=Random.GLOBAL_RNG,
)
    i = rand(rng, 1:length(sequence))
    a = sequence[i]
    b = mod(rand(rng, 1:19) + a, 20) + 1

    min_old_E = minimum(energies)
    Z = sum(enumerate(structures)) do (k, S)
        new_E = energies[k]
        for j in S.contact_partners[i]
            new_E += MJ_1996[b, sequence[j]] - MJ_1996[a, sequence[j]]
        end
        exp(-(new_E - min_old_E))
    end

    E_target = energies[target]
    for j in structures[target].contact_partners[i]
        E_target += MJ_1996[b, sequence[j]] - MJ_1996[a, sequence[j]]
    end

    p_target = exp(-(E_target - min_old_E)) / Z
    if rand() < p_target
        # NEED TO REPLACE OLD ENERGIES HERE
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
