module Benchmark 
    include("metropolis.jl")
    using .Metropolis
    using ProgressLogging
    using ProgressMeter
    using Plots
    using Statistics


# ------------------------------ FUNCTIONS ------------------------------
    function benchmark_beta()
        betas = [10, 50, 100, 1000]
        results = Dict{Float64, Tuple{Vector{Int}, Vector{Float64}}}()
        plots = []

        for beta in betas
            println("Testing beta = $beta")
            final_seq, p_struct_time_t, _ = algo(seq, target_number, all_ctc_maps, 4000, beta)
            results[beta] = (final_seq, p_struct_time_t)
            
            t = 1:length(p_struct_time_t)
            push!(plots, plot(t, p_struct_time_t, xlabel="Time", ylabel="Probability", 
                            title="β = $beta"))
        end

        fig = plot(plots..., layout=(2, 2), size=(1200, 600))
        display(fig)

        return results
    end



    function hamming_distance(s1::Vector{Int}, s2::Vector{Int})
        return sum(s1 .!= s2)
    end



    function find_clusters(saved_structures::Dict{Int, Vector{Int}}, distance_threshold=3)
        epochs = sort(collect(keys(saved_structures))) # sort of epochs (increasing order)
        structures = [saved_structures[e] for e in epochs]
        n = length(structures)
        
        # calcul of all Hamming distances between saved structures
        distances = zeros(Int, n, n)
        for i in 1:n
            for j in i+1:n
                d = hamming_distance(structures[i], structures[j])
                distances[i,j] = d
                distances[j,i] = d
            end
        end

        # Clustering based on distance threshold
        cluster_id = collect(1:n)  # each structure has its own cluster at the beginning

        function find_root(x) 
            """ 
            Traceback function to find the cluster to wich a structure belongs to.
            """
            while cluster_id[x] != x 
                x = cluster_id[x]
            end
            return x
        end

        for i in 1:n
            for j in i+1:n
                if distances[i,j] <= distance_threshold # belong to the same cluster
                    ri, rj = find_root(i), find_root(j)
                    if ri != rj # if they are not already in the same cluster
                        cluster_id[ri] = rj  # merge clusters
                    end
                end
            end
        end
        # We end up with a vector where each index points to the representative of its cluster (if cluster_id[5] = 3, then structure 5 and 3 are in the same cluster)

        # Gathering the results in a dict : key = index of the cluster, value = list of structures in the cluster
        clusters = Dict{Int, Vector{Vector{Int}}}()  # cluster_root => liste d'epochs
        for i in 1:n
            root = find_root(i)
            push!(get!(clusters, root, Vector{Int}[]), structures[i])
        end

        # Retrieving the mutated positions (loci) for each cluster
        cluster_loci = Dict{Int, Vector{Int}}()
        for (root, cluster_structs) in clusters
            if length(cluster_structs) > 1 # Verifing that there are more than 1 structure in the cluster

                ref = cluster_structs[1]
                differing_positions = Int[]
                for pos in 1:length(ref)
                    if any(s[pos] != ref[pos] for s in cluster_structs) # we do not store which structure has this mutation, just if at least one is mutated
                        push!(differing_positions, pos)
                    end
                end
                cluster_loci[root] = differing_positions
            end
        end

        return clusters, cluster_loci
    end



    function show_clusters(saved_structures, distance_threshold)
        clusters, cluster_loci = find_clusters(saved_structures, distance_threshold)

        println("Nombre de clusters $(length(clusters)) pour $(length(saved_structures)) structures sauvegardées.")

        for (root, epochs) in clusters
            println("Cluster $root : avec $(length(epochs)) structure(s) :\n epochs = $epochs")
            if haskey(cluster_loci, root)
                println("  Loci touchés : $(cluster_loci[root])")
            end
            println()
        end
    end


    function multi_init_final_seq(
        target_number::Int,
        all_ctc_maps::Vector;
        beta::Float64=100.0,
        runs::Int=10000
    )

        results = Dict{Int, Vector{Int}}()
        print("ici")

        @showprogress desc="Multi-init" for i in 1:runs

            seq = generate_seq()
            final_seq, _, _ = algo(seq, target_number, all_ctc_maps, runs, beta)        
            results[i] = final_seq
        end
        return results
    end

    function multi_init_seq_over_time(
        target_number::Int,
        all_ctc_maps::Vector;
        epochs::Int=4000,
        beta::Float64=100.0,
        threshold::Int=2000,
        step::Int=200,
        runs::Int=10
    )

        results = Dict{Int, Dict{Int, Vector{Int}}}()
        print("ici")

        @showprogress desc="Multi-init" for i in 1:runs

            seq = generate_seq()
            _, _, saved_struct = algo(seq, target_number, all_ctc_maps, epochs, beta, threshold, step, false)        
            results[i] = saved_struct
        end
        return results
    end


    function calculate_variances(sequences_over_time::Dict{Int, Dict{Int, Vector{Int}}})
        runs = sort(collect(keys(sequences_over_time))) # sort of runs indexes
        timestamps = sort(collect(intersect([keys(sequences_over_time[r]) for r in runs]...))) # making sure we only keep the shared timestamps for all runs (in case of missing data)
        
        n_runs = length(runs)
        n_times = length(timestamps)
        
        intra_var = zeros(n_times)
        intra_curves = [zeros(n_times) for _ in 1:n_runs]  # une courbe par run
        inter_var = zeros(n_times)

        for (t_idx, t) in enumerate(timestamps)
            seqs_at_t = [sequences_over_time[r][t] for r in runs]

            inter_distances = [hamming_distance(seqs_at_t[i], seqs_at_t[j])
                            for i in 1:n_runs for j in i+1:n_runs]
            inter_var[t_idx] = var(inter_distances)

            if t_idx > 1
                t_prev = timestamps[t_idx - 1]
                intra_distances = [hamming_distance(sequences_over_time[r][t], sequences_over_time[r][t_prev]) for r in runs]
                intra_var[t_idx] = var(intra_distances)
                for (d_idx, d) in enumerate(intra_distances)
                    intra_curves[d_idx][t_idx] = d  # retrieving the curves for each run
                end
            end
        end

        return timestamps, intra_var, intra_curves, inter_var
    end


    
    function plot_variances(timestamps, intra_curves, inter_var, n_curves::Int=5)
        p1 = plot(title="Intra-séquence", xlabel="Epoch", ylabel="Distance de Hamming")
        for r in 1:min(n_curves, length(intra_curves))
            plot!(p1, timestamps[2:end], intra_curves[r][2:end], label="Run $r", alpha=0.7)
        end
        p2 = plot(timestamps, inter_var, title="Inter-séquences", xlabel="Epoch", ylabel="Variance", label=false)
        display(plot(p1, p2, layout=(1,2), size=(1200, 400)))
    end


    function plot_distances_inter_intra(results::Dict{Int, Dict{Int, Vector{Int}}})
    runs = sort(collect(keys(results)))
    timestamps = sort(collect(intersect([keys(results[r]) for r in runs]...)))

    n_runs = length(runs)

    # distance seq(t) vs seq(t0) for each run
    intra_dist = [zeros(length(timestamps)) for _ in 1:n_runs]
    for (r_idx, r) in enumerate(runs)
        ref = results[r][timestamps[1]]  # seq(t0)
        for (t_idx, t) in enumerate(timestamps)
            intra_dist[r_idx][t_idx] = hamming_distance(ref, results[r][t])
        end
    end

    # mean dist at each timestamp 
    inter_mean = zeros(length(timestamps))
    for (t_idx, t) in enumerate(timestamps)
        seqs_at_t = [results[r][t] for r in runs]
        distances = [hamming_distance(seqs_at_t[i], seqs_at_t[j])
                    for i in 1:n_runs for j in i+1:n_runs]
        inter_mean[t_idx] = mean(distances)
    end

    # affichage sur le même graphique
    p = plot(title="Évolution des distances de séquences", xlabel="Epoch", ylabel="Distance de Hamming")
    for (r_idx, evo_dist) in enumerate(intra_dist) 
        plot!(p, timestamps, evo_dist, label="Run $r_idx", alpha=0.5, color=:blue)
    end
    plot!(p, timestamps, inter_mean, label="Distance inter moyenne", color=:red, linewidth=2)

    display(p)
    end

    function calcul_dist_final(results::Dict{Int, Vector{Int}})
        runs = sort(collect(keys(results)))
        final_seqs = [results[r] for r in runs]
        n_runs = length(runs)

        inter_distances = [hamming_distance(final_seqs[i], final_seqs[j])
                            for i in 1:n_runs for j in i+1:n_runs]
        return inter_distances
        
    end


# ------------------------------ Exports ------------------------------
    #Functions
    export benchmark_beta, hamming_distance, find_clusters, show_clusters, multi_init_seq_over_time, multi_init_final_seq , calculate_variances, plot_variances, plot_distances_inter_intra, calcul_dist_final

end