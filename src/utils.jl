function parse_Q_matrix(file)
    # a file like data/JTT_decomposition.csv is organized as such
    # the first 19 lines are the lower triangular part of the symetric exchange matrix
    # the 20th line contains the equilibrium probabilities
    dat = readdlm(file; comments=true) # contains "" at empty positions

    # parse S
    S = zeros(Float64, 20, 20)
    for a in 2:20, b in 1:(a - 1)
        S[a, b] = dat[a - 1, b]
        S[b, a] = S[a, b]
    end

    # parse Π
    π = Float64.(dat[end, :])
    @argcheck isapprox(sum(π), 1; rtol=1e-3)
    Π = diagm(π)

    # compute Q
    Q = [S[a, b] * π[b] for a in 1:20, b in 1:20]
    for a in 1:20
        Q[a, a] = -sum(Q[a, :])
    end
    @argcheck isapprox(maximum(abs, π' * Q), 0.0, atol=1e-6)

    R = sum(a -> -π[a] * Q[a, a], 1:20)
    !isapprox(R, 1.0; rtol=1e-3) && @warn "Average rage is not one. Got $R"

    return Q
end

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

const MJ_1996::Matrix{Float64} = let
    MJ = readdlm(joinpath(pkgdir(LatticeProteins), "data/MJ_1996.csv"), Float64) # 20x20
    alphabet_permutation = sortperm(MJ_alphabet)
    MJ[alphabet_permutation, alphabet_permutation]
end

const JTT_alphabet = split("ARNDCQEGHILKMFPSTWYV", "")
const JTT = let
    Q = parse_Q_matrix(joinpath(pkgdir(LatticeProteins), "data/JTT_decomposition.csv"))
    for a in 1:20
        Q[a, a] = 0.0
    end
    alphabet_permutation = sortperm(JTT_alphabet)
    Q[alphabet_permutation, alphabet_permutation]
end

const ref_alphabet = split("ACDEFGHIKLMNPQRSTVWY", "") # sorted alphabetically
