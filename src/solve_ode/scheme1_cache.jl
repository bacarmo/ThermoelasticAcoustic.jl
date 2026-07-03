function build_cache(::Scheme1, matrices::SystemMatrices{T, I}, m₃::I) where {T, I}
    Scheme1Cache(matrices, m₃)
end

struct Scheme1Cache{T, I, S1, S2}
    # --- scratch vectors ---
    vec_m₁::Vector{T}
    vec_m₂::Vector{T}
    vec_m₃_1::Vector{T}
    vec_m₃_2::Vector{T}
    vec_m₁_m₂::Vector{T}
    # --- rhs vectors ---
    τFf₁::Vector{T}
    τFf₂::Vector{T}
    τFf₃::Vector{T}
    L₁::Vector{T}
    L₂::Vector{T}
    # --- nonlinearity vectors ---
    G::Vector{T}
    F::Vector{T}
    Λ::Vector{T}
    # --- midpoint unknowns and auxiliary quantities ---
    X::Vector{T}
    v̂ⁿ::SubArray{T, 1, Vector{T}, Tuple{UnitRange{Int}}, true}
    ĉⁿ::SubArray{T, 1, Vector{T}, Tuple{UnitRange{Int}}, true}
    d̂ⁿ::Vector{T}
    K_m₂xm₂_ĉⁿ::Vector{T}
    r̂ⁿ::Vector{T}
    # --- nonlinear system ---
    minusH::Vector{T}
    Q₁₁::Symmetric{T, SparseMatrixCSC{T, I}}
    JH₁₁::Symmetric{T, SparseMatrixCSC{T, I}}
    JH₂₂_sparse::Symmetric{T, SparseMatrixCSC{T, I}}
    JH_sparse::SparseMatrixCSC{T, I}
    # --- vectors for JH_sparse synchronization ---
    map11_direct::Vector{I}
    map11_mirror::Vector{I}
    map22_direct::Vector{I}
    map22_mirror::Vector{I}
    # --- LinearSolve.LinearCache ---
    linsolve::S1
    linsolve_m₃::S2
end

function Scheme1Cache(matrices::SystemMatrices{T, I}, m₃::I) where {T, I}
    m₁ = size(matrices.K_m₁xm₁, 1)
    m₂ = size(matrices.K_m₂xm₂, 1)

    X = zeros(T, m₁ + m₂)
    v̂ⁿ = view(X, 1:m₁)
    ĉⁿ = view(X, (m₁ + 1):(m₁ + m₂))

    τFf₃ = zeros(T, m₃)

    minusH = zeros(T, m₁ + m₂)
    Q₁₁ = copy(matrices.K_m₁xm₁)
    JH₁₁ = copy(matrices.K_m₁xm₁)
    JH₂₂_sparse = copy(matrices.K_m₂xm₂)

    JH_sparse = [matrices._2M_m₁xm₁ matrices.τA_m₁xm₂;
                 matrices.τA_m₂xm₁ matrices._2M_m₂xm₂]::SparseMatrixCSC{T, I}

    map11_direct, map11_mirror = build_upper_to_full_maps11(JH_sparse, JH₁₁)
    map22_direct, map22_mirror = build_upper_to_full_maps22(JH_sparse, JH₂₂_sparse)

    # --- Linear system 1 ---
    A1 = JH_sparse
    b1 = zeros(T, m₁ + m₂)
    prob1 = LS.LinearProblem(A1, b1)
    linsolve = LS.init(prob1, LS.KLUFactorization(; reuse_symbolic = true, check_pattern = false))

    # --- Linear system 2 ---
    A2 = matrices.M_m₃xm₃::Symmetric{T, SparseMatrixCSC{T, I}}
    b2 = τFf₃
    prob2 = LS.LinearProblem(A2, b2)
    linsolve_m₃ = LS.init(prob2, LS.CHOLMODFactorization())

    return Scheme1Cache(
        # scratch vectors
        zeros(T, m₁), zeros(T, m₂), zeros(T, m₃), zeros(T, m₃), zeros(T, m₁+m₂),
        # rhs vectors
        zeros(T, m₁), zeros(T, m₂), τFf₃, zeros(T, m₁), zeros(T, m₂),
        # nonlinearity vectors
        zeros(T, m₃), zeros(T, m₁), zeros(T, m₂),
        # midpoint unknowns
        X, v̂ⁿ, ĉⁿ, zeros(T, m₁), zeros(T, m₂), zeros(T, m₃),
        # nonlinear system
        minusH, Q₁₁, JH₁₁, JH₂₂_sparse, JH_sparse,
        # vectors for JH_sparse synchronization
        map11_direct, map11_mirror,
        map22_direct, map22_mirror,
        # LinearSolve.LinearCache
        linsolve, linsolve_m₃
    )
end