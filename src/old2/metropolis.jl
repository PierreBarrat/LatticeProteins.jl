# ------------------------------ Functions ------------------------------
function all_edges_3x3x3()
    """
    We retreive all the edges from a square : we have 54 tuples. 
    """
    edges = Set{Tuple{Int,Int}}()

    for id in 1:27
        for neighbor in ADJ[id]
            edge = (min(id, neighbor), max(id, neighbor)) #to avoid duplicates (1,2) and (2,1)
            push!(edges, edge)
        end
    end

    return edges, length(edges)
end


function contact_map(path)
    """
    Retreive all the contact for a given path
    """
    in_path = Set{Tuple{Int,Int}}() # Set with all the edges used in the path

    for k in 1:length(path)-1
        a, b = path[k], path[k+1]
        a, b = a < b ? (a, b) : (b, a) # avoid duplicates (such as (1,2) and (2,1))
        push!(in_path, (a, b))
    end

    return setdiff(ALL_EDGES, in_path) # We only keep "non-used" eedges

end


#Generate a random 27 a.a sequence
function generate_seq()
    return rand(1:20, 27)
end

#Compute the total energy of one sequence for one structure
function energy_one_struct(ctc_map, seq)
    score = 0
    for contact in ctc_map
        pos1, pos2 = contact
        aa1 = seq[pos1]
        aa2 = seq[pos2]
        score += MJ[aa1, aa2]
    end
    return score
end

function all_contact_maps(paths)
    ctc_maps = []
    for path in paths
        push!(ctc_maps, contact_map(path))
    end
    return ctc_maps
end

#compute the denominator = sum of all the exponential of minus the energy 
function total_energy(ctc_maps, seq)
    total = 0
    for ctc_map in ctc_maps
        total += exp(-energy_one_struct(ctc_map, seq))
    end
    return total
end


#Compute the score/probability of one struture for one sequence among all the other structures for all the structures
function probability_each_structure(ctc_maps, seq, total_energy)
    proba = []
    for ctc_map in ctc_maps
        score = energy_one_struct(ctc_map, seq)
        p = exp(-score) / total_energy
        push!(proba, p)
    end
    return proba
end

function probability_target_structure(ctc_maps, seq, total_energy, target_number)
    s = energy_one_struct(ctc_maps[target_number], seq)
    return exp(-s) / total_energy
end


function mutation(seq)
    pos = rand(1:27)
    current_aa = seq[pos]

    #new aa different from actual
    new_aa = rand(1:19)
    if new_aa >= current_aa #we skip if the same aa
        new_aa += 1
    end

    return (pos, new_aa)
end


#pre-calculation (only one time before the algo)
function precompute_pos_contacts(all_ctc_maps)
    #pos_contacts[p]=list of (structure_idx, partner) for position p
    #for each position, which structures have a contact here and with whom?
    #pos_contacts never change
    pos_contacts = [Vector{Tuple{Int,Int}}() for _ in 1:27]
    for (s, ctc_map) in enumerate(all_ctc_maps)
        for (a, b) in ctc_map
            push!(pos_contacts[a], (s, b))
            push!(pos_contacts[b], (s, a))
        end
    end
    return pos_contacts #pos_contacts[pos 3] = [(structure 12, pos 7), etc...]
end


#We calculate weights[s] once for the initial sequence
#Then, instead of recalculating everything when the sequence changes, 
#we will only update the weights affected by the mutation. (cf. pos_contacts)
function init_weights(all_ctc_maps, seq)
    # weights[s] = exp(−E_s(seq)) for every structure s
    return [exp(-energy_one_struct(ctc_map, seq)) for ctc_map in all_ctc_maps]
end


#metropolis function
function metropolis(mut, seq, target_number, weights, pos_contacts, ttl_energy, log_proba_target, affected, beta=1000)
    pos, new_aa = mut
    old_aa = seq[pos]

    #only structures with contact in `pos` change 
    #affected = Dict{Int, Float64}() #indices, delta

    for (s, j) in pos_contacts[pos]
        affected[s] = get(affected, s, 0.0) + MJ[new_aa, seq[j]] - MJ[old_aa, seq[j]] #sum of delta(s)
    end

    #updating ttl_energy without recalculating everything
    new_ttl_energy = ttl_energy
    for (s, delta) in affected #we iterate only on the affected structures (else delta==0, useless)
        new_ttl_energy += weights[s] * exp(-delta) - weights[s]
    end

    #equivalent to probability_target_structure
    dE_target = get(affected, target_number, 0.0)
    new_log_proba = log(weights[target_number] * exp(-dE_target)) - log(new_ttl_energy)

    #metropolis criterion
    accept = false
    if new_log_proba >= log_proba_target
        accept = true
    end
    if rand() <= exp(beta * (new_log_proba - log_proba_target))
        accept = true
    end

    "=
    if accept
        new_seq=copy(seq) #copy only if we accept
        new_seq[pos]=new_aa
        for (s, delta) in affected #update of the weights vector in place
            weights[s] *= exp(-delta)
        end
        return new_seq, new_ttl_energy, new_log_proba
    else
        return seq, ttl_energy, log_proba_target
    end
    ="

    if accept
        seq[pos] = new_aa  # modifie en place
        for (s, delta) in affected
            weights[s] *= exp(-delta)
        end
        return new_ttl_energy, new_log_proba, true
    else
        return ttl_energy, log_proba_target, false
    end

end



#algo function
function algo(seq, target_number, all_ctc_maps, epochs=1000, beta=100, threshold=2000, step=200, show_progress=true)
    #pre-calculations (only once)
    pos_contacts = precompute_pos_contacts(all_ctc_maps)
    weights = init_weights(all_ctc_maps, seq)
    ttl_energy = sum(weights)
    actual_seq = copy(seq)
    log_proba_target = log(weights[target_number]) - log(ttl_energy)
    p_struct_time_t = Float64[]
    #saved_structures = save_struct ? Dict{Int, Vector{Int}}() : nothing  # alloué seulement si nécessaire
    saved_structures = Dict{Int,Vector{Int}}()


    affected = Dict{Int,Float64}()  #alloué une seule fois

    @showprogress enabled = show_progress desc = "algo (β=$beta)" for i in 1:epochs
        #recalibration périodique, pour éviter la dérive a cause des arrondis et donc log(négatif)
        ttl_energy = sum(weights)  #recalcul exact

        push!(p_struct_time_t, exp(log_proba_target))

        empty!(affected)  # vider sans réallouer

        mut = mutation(actual_seq)
        #actual_seq, ttl_energy, log_proba_target = metropolis(mut, actual_seq, target_number, weights, pos_contacts, ttl_energy, log_proba_target, affected, beta)
        ttl_energy, log_proba_target, _ = metropolis(mut, actual_seq, target_number, weights, pos_contacts, ttl_energy, log_proba_target, affected, beta)
        if i >= threshold && i % step == 0
            saved_structures[i] = copy(actual_seq)
        end

    end
    return actual_seq, p_struct_time_t, saved_structures
end


function proba_of_seq(seq, all_ctc_maps, target_number)
    weights = init_weights(all_ctc_maps, seq)
    ttl_energy = sum(weights)
    log_proba_target = exp(log(weights[target_number]) - log(ttl_energy))

    return log_proba_target
end


function algo_count_only_if_mut_accepted(seq, target_number, all_ctc_maps, epochs=1000, beta=100, threshold=2000, step=200, show_progress=true)
    #pre-calculations (only once)
    pos_contacts = precompute_pos_contacts(all_ctc_maps)
    weights = init_weights(all_ctc_maps, seq)
    ttl_energy = sum(weights)
    actual_seq = copy(seq)
    log_proba_target = log(weights[target_number]) - log(ttl_energy)
    p_struct_time_t = Float64[]
    #saved_structures = save_struct ? Dict{Int, Vector{Int}}() : nothing  # alloué seulement si nécessaire
    saved_structures = Dict{Int,Vector{Int}}()


    affected = Dict{Int,Float64}()  #alloué une seule fois
    cpt_mut_accepted = 0

    while cpt_mut_accepted < epochs
        #recalibration périodique, pour éviter la dérive a cause des arrondis et donc log(négatif)
        ttl_energy = sum(weights)  #recalcul exact

        push!(p_struct_time_t, exp(log_proba_target))

        empty!(affected)  # vider sans réallouer

        mut = mutation(actual_seq)
        #old_seq = copy(actual_seq)

        #actual_seq, ttl_energy, log_proba_target = metropolis(mut, actual_seq, target_number, weights, pos_contacts, ttl_energy, log_proba_target, affected, beta)
        ttl_energy, log_proba_target, accepted = metropolis(mut, actual_seq, target_number, weights, pos_contacts, ttl_energy, log_proba_target, affected, beta)

        #if seq didn't change with metropolis, we don't count this iteration
        if accepted 
            cpt_mut_accepted += 1
        end

        if cpt_mut_accepted >= threshold && cpt_mut_accepted % step == 0
            saved_structures[cpt_mut_accepted] = copy(actual_seq)
        end

    end
    return actual_seq, p_struct_time_t, saved_structures
end



function algo_show_probability_distribution(seq, target_number, all_ctc_maps, epochs=100, beta=100, save_struct=false, threshold=2000)
    #pre-calculations (only once)
    pos_contacts = precompute_pos_contacts(all_ctc_maps)
    weights = init_weights(all_ctc_maps, seq)
    ttl_energy = sum(weights)
    actual_seq = copy(seq)
    log_proba_target = log(weights[target_number]) - log(ttl_energy)
    p_struct_time_t = Float64[]
    saved_structures = Dict{Int,Vector{Int}}()
    affected = Dict{Int,Float64}()

    @showprogress desc = "algo (β=$beta)" for i in 1:epochs
        ttl_energy = sum(weights)

        push!(p_struct_time_t, exp(log_proba_target))

        empty!(affected)
        mut = mutation(actual_seq)
        ttl_energy, log_proba_target, _ = metropolis(mut, actual_seq, target_number, weights, pos_contacts, ttl_energy, log_proba_target, affected, beta)

        if save_struct && i >= threshold && i % 200 == 0
            saved_structures[i] = copy(actual_seq)
        end
    end

    p_final = weights ./ ttl_energy  # distribution finale, calculée une seule fois

    return actual_seq, p_struct_time_t, p_final, saved_structures
end



function show_sorted_indices(p_final)
    sorted_indices = sortperm(p_final, rev=true)
    for i in sorted_indices[1:20]
        println("structure $i : $(p_final[i])")
    end
end



function plot_probability(seq, target_number, all_ctc_maps, epochs=1000, beta=100, threshold=2000, step=200, show_progress=true, matrix_name=matrix_name)
    _, p_struct_time_t, _ = algo(seq, target_number, all_ctc_maps, epochs, beta, threshold, step, show_progress)
    t = 1:length(p_struct_time_t)
    p = plot(t, p_struct_time_t, xlabel="Time", ylabel="Probability", title="Evolution of the probability of the sequence over time")
    savefig(p, joinpath(@__DIR__, "../figures/probability_$(length(all_ctc_maps))_$(matrix_name)_beta_$(beta).svg"))
    return p
end


function choose_MJ(matrix="MJ_1996")
    if matrix == "MJ_1996"
        dir = joinpath(@__DIR__, "../MJ_1996.csv")
        matrix_name = "MJ_1996"
        df = CSV.read(dir, DataFrame; delim='\t', header=false)
        final_matrix = Matrix{Float64}(df)
    elseif matrix == "MJ"
        dir = joinpath(@__DIR__, "../MJ.csv")
        matrix_name = "MJ*5"
        df = CSV.read(dir, DataFrame; delim='\t', header=false)
        final_matrix = Matrix{Float64}(df)
        final_matrix .*= 5.0 # Multiplie tous les éléments de la matrice par 5
    else
        error("Matrix not recognized. Please choose either 'MJ_1996' or 'MJ'.")
    end

    return final_matrix, matrix_name
end


# ------------------------------ Variables ------------------------------
edges, count = all_edges_3x3x3()
const ALL_EDGES = edges
const aa = ["CYS", "MET", "PHE", "ILE", "LEU", "VAL", "TRP", "TYR", "ALA", "GLY", "THR", "SER", "GLN", "ASN", "GLU", "ASP", "HIS", "ARG", "LYS", "PRO"]
const clas_aa = ["-", "+", "pol_sans_char", "apol_aliph", "arom"]
const aa_idx = Dict(i => aa[i] for i in eachindex(aa)) # change the logic of the dictionnary, now key=1:20 and value=aa
const aa_class = Dict(1 => "pol_sans_char", 2 => "apol_aliph", 3 => "arom", 4 => "apol_aliph", 5 => "apol_aliph", 6 => "apol_aliph", 7 => "arom", 8 => "arom", 9 => "pol_sans_char", 10 => "pol_sans_char", 11 => "pol_sans_char", 12 => "pol_sans_char", 13 => "pol_sans_char", 14 => "pol_sans_char", 15 => "pol_sans_char", 16 => "+", 17 => "+", 18 => "-", 19 => "-", 20 => "pol_sans_char")

# Reading the file and converting into a matrix 
const MJ, matrix_name = choose_MJ("MJ") # Carefull, its a const for efficiency, if we want to change it, we need to restart the kernel.