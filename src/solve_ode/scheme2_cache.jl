function build_cache(::Scheme2, matrices::SystemMatrices{T, I}, m₃::I) where {T, I}
    Scheme2Cache(matrices, m₃)
end

struct Scheme2Cache{T, I, S1, S2, S3}
    # --- scratch vectors ---
    vec_m₁_1::Vector{T}
    vec_m₁_2::Vector{T}
    vec_m₂::Vector{T}
    vec_m₃_1::Vector{T}
    vec_m₃_2::Vector{T}
    # --- rhs vectors ---
    τFf₁::Vector{T}
    τFf₂::Vector{T}
    τFf₃::Vector{T}
    L₁::Vector{T}
    L₂::Vector{T}
    # --- nonlinearity vectors ---
    G::Vector{T}
    F::Vector{T}
    # --- midpoint unknowns and auxiliary quantities ---
    v̂ⁿ::Vector{T}
    dᵃⁿ::Vector{T}
    dⁿ⁻²::Vector{T}
    cᵃⁿ::Vector{T}
    cⁿ⁻²::Vector{T}
    r̂ⁿ::Vector{T}
    # --- nonlinear system ---
    minusH::Vector{T}
    Q₁::Symmetric{T, SparseMatrixCSC{T, I}}
    Q₂_upper::Symmetric{T, SparseMatrixCSC{T, I}}
    Q₂::SparseMatrixCSC{T, I}
    JH_upper::Symmetric{T, SparseMatrixCSC{T, I}}
    JH::SparseMatrixCSC{T, I}
    # --- vectors for JH synchronization ---
    map11_direct::Vector{I}
    map11_mirror::Vector{I}
    map22_direct::Vector{I}
    map22_mirror::Vector{I}
    # --- LinearSolve.LinearCache ---
    linsolve_m₁::S1
    linsolve_m₂::S2
    linsolve_m₃::S3
end

function Scheme2Cache(matrices::SystemMatrices{T, I}, m₃::I) where {T, I}
    m₁ = size(matrices.K_m₁xm₁, 1)
    m₂ = size(matrices.K_m₂xm₂, 1)

    L₁ = zeros(T, m₁)
    L₂ = zeros(T, m₂)
    τFf₃ = zeros(T, m₃)

    Q₁ = copy(matrices.K_m₁xm₁)
    Q₂_upper = copy(matrices.K_m₂xm₂)
    Q₂ = sparse(matrices.K_m₂xm₂)::SparseMatrixCSC{T, I}

    minusH = zeros(T, m₁)
    JH_upper = copy(matrices.K_m₁xm₁)
    JH = sparse(matrices.K_m₁xm₁)::SparseMatrixCSC{T, I}

    map11_direct, map11_mirror = build_upper_to_full_maps11(JH, JH_upper)
    map22_direct, map22_mirror = build_upper_to_full_maps11(Q₂, Q₂_upper)

    # --- Linear system 1 ---
    A1 = JH
    b1 = minusH
    prob1 = LS.LinearProblem(A1, b1)
    linsolve_m₁ = LS.init(prob1, LS.KLUFactorization(; reuse_symbolic = true, check_pattern = false))

    # --- Linear system 2
    A2 = Q₂
    b2 = L₂
    prob2 = LS.LinearProblem(A2, b2)
    linsolve_m₂ = LS.init(prob2, LS.KLUFactorization(; reuse_symbolic = true, check_pattern = false))

    # --- Linear system 3 ---
    A3 = matrices.M_m₃xm₃::Symmetric{T, SparseMatrixCSC{T, I}}
    b3 = τFf₃
    prob3 = LS.LinearProblem(A3, b3)
    linsolve_m₃ = LS.init(prob3, LS.CHOLMODFactorization())

    return Scheme2Cache(
        # scratch vectors
        zeros(T, m₁), zeros(T, m₁), zeros(T, m₂), zeros(T, m₃), zeros(T, m₃),
        # rhs vectors
        zeros(T, m₁), zeros(T, m₂), τFf₃, L₁, L₂,
        # nonlinearity vectors
        zeros(T, m₃), zeros(T, m₁),
        # midpoint unknowns
        zeros(T, m₁), zeros(T, m₁), zeros(T, m₁),
        zeros(T, m₂), zeros(T, m₂), zeros(T, m₃),
        # nonlinear system
        minusH, Q₁, Q₂_upper, Q₂, JH_upper, JH,
        # vectors for JH synchronization
        map11_direct, map11_mirror, map22_direct, map22_mirror,
        # LinearSolve.LinearCache 
        linsolve_m₁, linsolve_m₂, linsolve_m₃
    )
end