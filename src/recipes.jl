using RecipesBase

@recipe function f(s::Structure)
    N = Int(LatticeProteins.Structures._infer_N(s.path))
    chain = s.path

    legend --> false
    aspect_ratio --> :equal
    showaxis --> false
    grid --> true
    xticks --> (1:0.5:N, fill("", 2N - 1))
    yticks --> (1:0.5:N, fill("", 2N - 1))
    zticks --> (1:0.5:N, fill("", 2N - 1))
    xlim --> (1, N)
    ylim --> (1, N)
    zlim --> (1, N)

    # All lattice sites
    @series begin
        seriestype := :scatter3d
        markercolor := :blue
        markersize := 7
        markerstrokewidth := 0
        label := ""
        xs = vec([x for x in 1:N, y in 1:N, z in 1:N])
        ys = vec([y for x in 1:N, y in 1:N, z in 1:N])
        zs = vec([z for x in 1:N, y in 1:N, z in 1:N])
        xs, ys, zs
    end

    # Path
    @series begin
        seriestype := :path3d
        linecolor := :black
        linewidth := 2
        label := ""
        [s.x for s in chain], [s.y for s in chain], [s.z for s in chain]
    end

    # Start (green)
    @series begin
        seriestype := :scatter3d
        markercolor := :green
        markersize := 8
        markerstrokewidth := 0
        label := ""
        [chain[1].x], [chain[1].y], [chain[1].z]
    end

    # End (red)
    @series begin
        seriestype := :scatter3d
        markercolor := :red
        markersize := 8
        markerstrokewidth := 0
        label := ""
        [chain[end].x], [chain[end].y], [chain[end].z]
    end
end
