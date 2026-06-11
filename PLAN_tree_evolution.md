# Plan: tree evolution in LatticeProteins.jl

## Dependencies
- [x] Add `TreeTools` to `LatticeProteins.jl` via `Pkg.add("TreeTools")` (registered version)
  - ended up using local dev v0.8.0 via `Pkg.develop`; compat set to `"0.7, 0.8"`

## File changes

### `src/metropolis.jl` — refactor, no behavior change
- [x] Extract `_mcmc_init(init, parameters)` → `(MCMCState, contact_tensor, field)`
- [x] Extract `_advance_mcmc!(state, contact_tensor, field, target, β, n_steps::Int; JTT_bias, rng)` → runs N steps in-place
- [x] Add `Base.copy(state::MCMCState)`: deep-copy `sequence` and `energies`, fresh `similar` for scratch buffers (`energy_buffer`, `delta_mj`), copy `ϕ_old`
- [x] Refactor `sample_mcmc_chain` to use `_mcmc_init` and `_advance_mcmc!`

### `src/tree_evolution.jl` — new file
- [x] `evolve_on_tree(root_sequence::Vector{Int}, tree::Tree, parameters::MCMCParameters; rng=Random.GLOBAL_RNG)`
  - Build `contact_tensor` and `field` once via `_mcmc_init`
  - Root state: `MCMCState` initialised from `root_sequence` directly (no burnin — caller's responsibility)
  - Traversal: `preorder_traversal(tree; root=false)`, parent lookup via `ancestor(node)`
  - Per node: `copy(parent_state)`, `n_steps = round(Int, branch_length(node))`, `_advance_mcmc!(..., n_steps)`
  - Returns `(; leaf_sequences::Dict{String,Vector{Int}}, node_sequences::Dict{String,Vector{Int}})`

### `src/LatticeProteins.jl`
- [x] `include("tree_evolution.jl")`
- [x] `export evolve_on_tree`

## Contracts / invariants
- Branch lengths must be (approximately) integer-valued floats; `round(Int, t)` is the conversion; errors propagate naturally for `NaN`/`missing`
- Zero branch length is valid: child receives an exact copy of the parent sequence
- `root_sequence` is expected to already be at equilibrium; `parameters.burnin` is ignored by `evolve_on_tree`
- `parameters.n_steps` and `parameters.n_sequences` are ignored by `evolve_on_tree`
