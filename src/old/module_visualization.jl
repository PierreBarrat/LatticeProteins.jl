module Visualization
    import Main.HamiltonianPath: ADJ, id_to_coords, coords_to_id
    using Plots
    plotlyjs()

    # ------------------------------ Functions ------------------------------
    function all_edges_3x3x3()
        edges = Set{Tuple{Int,Int}}()
        for id in 1:27
            for neighbor in ADJ[id]
                edge = (min(id, neighbor), max(id, neighbor))
                push!(edges, edge)
            end
        end
        return edges, length(edges)
    end

    function get_nodes(range = 0:2)
        points = vec(collect(Iterators.product(range, range, range)))
        x = [p[1] for p in points]
        y = [p[2] for p in points]
        z = [p[3] for p in points]
        return x, y, z
    end

    function get_edges(range = 0:2)
        start, stop = first(range), last(range)
        n = length(range)
        total_elements = 3 * 3 * n^2
        lx = Vector{Float64}(undef, total_elements)
        ly = Vector{Float64}(undef, total_elements)
        lz = Vector{Float64}(undef, total_elements)

        idx = 1
        for i in range, j in range
            lx[idx:idx+2] .= [start, stop, NaN]; ly[idx:idx+2] .= [i, i, NaN]; lz[idx:idx+2] .= [j, j, NaN]
            idx += 3
            lx[idx:idx+2] .= [i, i, NaN]; ly[idx:idx+2] .= [start, stop, NaN]; lz[idx:idx+2] .= [j, j, NaN]
            idx += 3
            lx[idx:idx+2] .= [i, i, NaN]; ly[idx:idx+2] .= [j, j, NaN]; lz[idx:idx+2] .= [start, stop, NaN]
            idx += 3
        end
        return lx, ly, lz
    end

    function contact_map(path)
        in_path = Set{Tuple{Int,Int}}()
        for k in 1:length(path)-1
            a, b = path[k], path[k+1]
            a, b = a < b ? (a,b) : (b,a)
            push!(in_path, (a, b))
        end
        return setdiff(ALL_EDGES, in_path)
    end

    function plot_cube_ctc_map(x, y, z, lx, ly, lz, edges_to_highlight)
        p = plot(lx, ly, lz, color=:black, alpha=0.5, lw=1, label="")
        
        hx, hy, hz = Float64[], Float64[], Float64[]
        for (n1, n2) in edges_to_highlight
            c1 = id_to_coords(n1)
            c2 = id_to_coords(n2)
            append!(hx, [c1[1], c2[1], NaN])
            append!(hy, [c1[2], c2[2], NaN])
            append!(hz, [c1[3], c2[3], NaN])
        end
        
        plot!(p, hx, hy, hz, color=:red, lw=3, label="Contact")

        scatter!(p, x, y, z, 
                markersize=4, 
                markercolor=:blue,
                xlims = (0, 2),
                ylims = (0, 2),
                zlims = (0, 2),
                hover = indices_hover,
                showaxis=false,
                grid=false,
                ticks=false,
                aspect_ratio=:equal)

        return p
    end

    function plot_cube_path(x, y, z, lx, ly, lz, path)
        p = plot(lx, ly, lz, color=:black, alpha=0.5, lw=0.5, label="")

        n_steps = length(path) - 1
        colors = cgrad(:plasma, n_steps, categorical=true)  

        for i in 1:n_steps
            n1, n2 = path[i], path[i+1]
            c1 = id_to_coords(n1)
            c2 = id_to_coords(n2)

            plot!(p,
                [c1[1], c2[1]],
                [c1[2], c2[2]],
                [c1[3], c2[3]],
                color = colors[i],
                lw = 3,
                label = ""
            )
        end

        # highlight start and end points
        c_start = id_to_coords(path[1])
        c_end   = id_to_coords(path[end])

        scatter!(p,
            [c_start[1]], [c_start[2]], [c_start[3]],
            markersize=8, markercolor=:blue, label="Start"
        )
        scatter!(p,
            [c_end[1]], [c_end[2]], [c_end[3]],
            markersize=8, markercolor=:yellow, label="End"
        )

        scatter!(p, x, y, z,
            markersize=4,
            markercolor=:blue,
            xlims=(0, 2), ylims=(0, 2), zlims=(0, 2),
            hover=indices_hover,
            showaxis=false,
            grid=false,
            ticks=false,
            aspect_ratio=:equal
        )

        return p
    end

    # ------------------------------ Constants ------------------------------
    const ALL_EDGES, _ = all_edges_3x3x3()
    const x, y, z = invokelatest(get_nodes)
    const lx, ly, lz = invokelatest(get_edges)
    const indices_hover = invokelatest() do
        [coords_to_id(xi, yi, zi) for (xi, yi, zi) in zip(x, y, z)]
    end

    # ------------------------------ Exports ------------------------------
    #Functions
    export plot_cube_ctc_map, plot_cube_path, contact_map, x, y, z

end