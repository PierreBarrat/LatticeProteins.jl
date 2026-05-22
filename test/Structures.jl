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
    end

    @testset "infer N" begin
        st = S.Structure(PATH2)
        @test st isa S.Structure{2}
    end
end
