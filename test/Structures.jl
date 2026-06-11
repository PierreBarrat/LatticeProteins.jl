using Test
using LatticeProteins

const S = LatticeProteins.Structures

@testset "Site" begin
    s = S.site(1, 2, 3)
    @test s isa S.Site
    @test s == (x=Int8(1), y=Int8(2), z=Int8(3))
end

@testset "is_neighbour" begin
    s = S.site(1, 1, 1)
    @test S.is_neighbour(s, S.site(2, 1, 1))   # +x
    @test S.is_neighbour(s, S.site(1, 2, 1))   # +y
    @test S.is_neighbour(s, S.site(1, 1, 2))   # +z
    @test !S.is_neighbour(s, s)                  # same site
    @test !S.is_neighbour(s, S.site(2, 2, 1))  # diagonal
    @test !S.is_neighbour(s, S.site(3, 1, 1))  # distance 2
end

# A valid Hamiltonian path through the 2x2x2 lattice
const PATH2 = [
    S.site(1, 1, 1),
    S.site(1, 1, 2),
    S.site(1, 2, 2),
    S.site(1, 2, 1),
    S.site(2, 2, 1),
    S.site(2, 2, 2),
    S.site(2, 1, 2),
    S.site(2, 1, 1),
]

@testset "Structure{N} constructor" begin
    @testset "N=1" begin
        st = S.Structure{1}([S.site(1, 1, 1)])
        @test length(st.path) == 1
        @test isempty(st.contacts)
    end

    @testset "N=2" begin
        st = S.Structure{2}(PATH2)
        @test length(st.path) == 8
        # site 1 (1,1,1) and site 8 (2,1,1) are lattice neighbours but not consecutive
        @test (1, 8) in st.contacts
        # consecutive pairs must not appear as contacts
        @test !any(((i, j),) -> j == i + 1, st.contacts)
        # contact_partners is consistent with contacts
        @test length(st.contact_partners) == 8
        for (i, j) in st.contacts
            @test j in st.contact_partners[i]
            @test i in st.contact_partners[j]
        end
    end

    @testset "infer N" begin
        st = S.Structure(PATH2)
        @test st isa S.Structure{Int8(2)}
        @test typeof(S.Structure{2}(PATH2)) == typeof(S.Structure(PATH2))
    end
end

@testset "site ID roundtrip" begin
    for N in (1, 2, 3)
        for x in 1:N, y in 1:N, z in 1:N
            s = S.site(x, y, z)
            id = S._site_id(s, N)
            @test 1 ≤ id ≤ N^3
            @test S._id_to_site(id, N) == s
        end
    end
    # boundary checks
    @test S._site_id(S.site(1, 1, 1), 3) == 1
    @test S._site_id(S.site(3, 3, 3), 3) == 27
end

@testset "site_index_map" begin
    @testset "N=1" begin
        st = S.Structure{1}([S.site(1, 1, 1)])
        @test S.site_index_map(st) == [1]
    end

    @testset "N=2 (explicit)" begin
        st = S.Structure{2}(PATH2)
        M = S.site_index_map(st)
        @test length(M) == 8
        # inverse of the known path: chain position sitting at each site
        @test M == [1, 2, 4, 3, 8, 7, 5, 6]
    end

    @testset "is the inverse of path" begin
        for st in (S.Structure{2}(PATH2), S.generate_structures(2)...)
            N = 2
            M = S.site_index_map(st)
            @test length(M) == N^3
            # M is a permutation of chain positions (path visits every site once)
            @test sort(M) == collect(1:(N^3))
            # documented property: path[M[site_id(s)]] == s
            for s in st.path
                @test st.path[M[S._site_id(s, N)]] == s
            end
            # and the other direction: M[site_id(path[j])] == j
            for (j, s) in enumerate(st.path)
                @test M[S._site_id(s, N)] == j
            end
        end
    end

    @testset "N=3 round-trip" begin
        st = first(S.generate_structures(3))
        N = 3
        M = S.site_index_map(st)
        @test length(M) == 27
        @test sort(M) == collect(1:27)
        for (j, s) in enumerate(st.path)
            @test M[S._site_id(s, N)] == j
            @test st.path[M[S._site_id(s, N)]] == s
        end
    end
end

@testset "rotation table" begin
    for N in (1, 2, 3)
        table = S.generate_rotation_table(N)
        @test size(table) == (N^3, 48)
        # each rotation is a bijection on site IDs
        for r in axes(table, 2)
            @test sort(table[:, r]) == Int8.(1:N^3)
        end
        # identity rotation must appear
        @test any(r -> table[:, r] == Int8.(1:N^3), axes(table, 2))
    end
end

@testset "canonical" begin
    table2 = S.generate_rotation_table(2)
    # canonical is idempotent
    c = S.canonical(PATH2, table2, 2)
    @test S.canonical(c, table2, 2) == c
    # all 48 rotations of a path share the same canonical form
    for r in axes(table2, 2)
        rotated = [S._id_to_site(table2[S._site_id(s, 2), r], 2) for s in PATH2]
        @test S.canonical(rotated, table2, 2) == c
    end
end

@testset "build_adj" begin
    adj3 = S.build_adj(3)
    # corner (1,1,1): 3 neighbours
    @test length(adj3[1, 1, 1]) == 3
    # face centre (2,2,1): 5 neighbours
    @test length(adj3[2, 2, 1]) == 5
    # interior (2,2,2): 6 neighbours
    @test length(adj3[2, 2, 2]) == 6
    # adjacency is symmetric
    for x in 1:3, y in 1:3, z in 1:3
        s = S.site(x, y, z)
        for nb in adj3[x, y, z]
            @test s in adj3[nb.x, nb.y, nb.z]
        end
    end
end

@testset "generate_structures N=2" begin
    structs = S.generate_structures(2)
    @test !isempty(structs)
    @test all(st -> st isa S.Structure{Int8(2)}, structs)
    # all paths are distinct
    @test allunique(st.path for st in structs)
end
