using TreeTools:
    Tree, TreeNode, preorder_traversal, ancestor, branch_length, label, isleaf, root

"""
    evolve_on_tree(
        root_sequence::Vector{Int},
        tree::Tree,
        parameters::MCMCParameters;
        rng = Random.GLOBAL_RNG,
        progress = false,
    )

Evolve `root_sequence` along `tree` using Metropolis MCMC.
Branch lengths must be (approximately) integer-valued: each branch of length `t` runs
`round(Int, t)` MCMC steps.

`root_sequence` is used as-is — no burnin is applied. The caller is responsible for
passing an equilibrated sequence.

Returns a named tuple `(; leaf_sequences, node_sequences)` where both values are
`Dict{String, Vector{Int}}` keyed by node label.

Set `progress=true` to display a progress bar over the tree nodes.
"""
function evolve_on_tree(
    root_sequence::Vector{Int},
    tree::Tree,
    parameters::MCMCParameters;
    rng=Random.GLOBAL_RNG,
    progress=false,
)
    @unpack target, β_sampling, JTT_bias = parameters
    root_state, contact_tensor, field = _mcmc_init(root_sequence, parameters)

    states = Dict{String,MCMCState}()
    states[label(root(tree))] = root_state

    leaf_sequences = Dict{String,Vector{Int}}()
    node_sequences = Dict{String,Vector{Int}}()
    node_sequences[label(root(tree))] = copy(root_state.sequence)

    prog = Progress(length(nodes(tree)) - 1; enabled=progress, desc="Evolving on tree")
    for node in preorder_traversal(tree; root=false)
        parent_state = states[label(ancestor(node))]
        child_state = copy(parent_state)
        n_steps = round(Int, branch_length(node))
        _advance_mcmc!(
            child_state, contact_tensor, field, target, β_sampling, n_steps; JTT_bias, rng
        )
        lbl = label(node)
        states[lbl] = child_state
        if isinternal(node)
            node_sequences[lbl] = copy(child_state.sequence)
        end
        if isleaf(node)
            leaf_sequences[lbl] = copy(child_state.sequence)
        end
        next!(prog)
    end

    return (; leaf_sequences, node_sequences)
end
