function build_cache(::CrankNicolson, matrices, dof_map_m₁, dof_map_m₂, dof_map_m₃, τ)
    CrankNicolsonCache(matrices, dof_map_m₁, dof_map_m₂, dof_map_m₃, τ)
end

"""
    CrankNicolsonCache{T, I, TF, TC}

Pre-allocated workspace for the Crank-Nicolson Strategy 2 solver (joint Newton iteration in ``[\\hat{v}^n;\\,\\hat{c}^n]``).

# Fields
- `vec_m₁_1`, `vec_m₁_2`, `vec_m₁_3`: scratch vectors of length `m₁`.
- `vec_m₂_1`, `vec_m₂_2`: scratch vectors of length `m₂`.
- `vec_m₃_1`, `vec_m₃_2`: scratch vectors of length `m₃`.
- `L₁`, `L₂`: right-hand sides of the ``\\hat{v}^n`` and ``\\hat{c}^n`` systems, lengths `m₁` and `m₂`.
- `Xⁿ`: joint Newton iterate ``[\\hat{v}^n;\\,\\hat{c}^n]``, length `m₁+m₂`.
- `v̂ⁿ`, `ĉⁿ`: views into `Xⁿ[1:m₁]` and `Xⁿ[m₁+1:end]`.
- `r̂ⁿ`: midpoint acoustic unknown ``\\hat{r}^n``, length `m₃`.
- `d̂ⁿ`: midpoint wave displacement ``\\hat{d}^n = \\tfrac{d^n + d^{n-1}}{2} = \\tfrac{\\tau}{2}\\hat{v}^n + d^{n-1}``, length `m₁`.
- `K_m₂xm₂_vs_ĉⁿ`: product ``K^{m_2\\times m_2}\\hat{c}^n``, length `m₂`.
- `QXⁿ`: product ``Q X^n``, length `m₁+m₂`.
- `minusH`: residual ``-H(X^n)``, length `m₁+m₂`.
- `JH₁₁`: (1,1) block of the Jacobian ``JH``, size `m₁×m₁`.
- `JH₂₂_sparse`: (2,2) block of ``JH`` excluding its rank-1 matrix, size `m₂×m₂`.
- `JH_sparse`: global sparse Jacobian (without rank-1 matrix), size `(m₁+m₂)×(m₁+m₂)`; off-diagonal blocks are time-invariant and set at construction.
- `ΔXⁿ`: Newton update, length `m₁+m₂`.
- `Q₁₁`: (1,1) block of ``Q``, updated each time step; size `m₁×m₁`.
- `τA_m₁xm₂`, `τA_m₂xm₁`: time-invariant off-diagonal blocks ``\\tau A^{m_1\\times m_2}`` and ``\\tau A^{m_2\\times m_1}``.
- `M_m₁xm₁_vs2`, `M_m₂xm₂_vs2`: scaled mass matrices ``2M^{m_1\\times m_1}`` and ``2M^{m_2\\times m_2}``.
- `M_m₃xm₃_chol`: Cholesky factorization of ``M^{m_3\\times m_3}``.
- `linsolve`: cached linear solver handle for ``JH\\,\\Delta X^n = -H``.
- `map_direct₁₁`, `map_mirror₁₁`: index maps from `JH₁₁.data.nzval` into `JH_sparse.nzval`.
- `map_direct₂₂`, `map_mirror₂₂`: index maps from `JH₂₂_sparse.data.nzval` into `JH_sparse.nzval`.
"""
struct CrankNicolsonCache{T <: AbstractFloat, I <: Integer, TF, TC <: LS.LinearCache}
    # --- Scratch vectors — m₁ ---
    vec_m₁_1::Vector{T}
    vec_m₁_2::Vector{T}
    vec_m₁_3::Vector{T}
    # --- Scratch vectors — m₂ ---
    vec_m₂_1::Vector{T}
    vec_m₂_2::Vector{T}
    # --- Scratch vectors — m₃ ---
    vec_m₃_1::Vector{T}
    vec_m₃_2::Vector{T}
    # --- RHS vectors ---
    L₁::Vector{T}
    L₂::Vector{T}
    # --- Midpoint unknowns and auxiliary quantities ---
    Xⁿ::Vector{T}
    v̂ⁿ::SubArray{T, 1, Vector{T}, Tuple{UnitRange{Int}}, true}
    ĉⁿ::SubArray{T, 1, Vector{T}, Tuple{UnitRange{Int}}, true}
    r̂ⁿ::Vector{T}
    d̂ⁿ::Vector{T}
    K_m₂xm₂_vs_ĉⁿ::Vector{T}
    QXⁿ::Vector{T}
    # --- Newton system ---
    minusH::Vector{T}
    JH₁₁::Symmetric{T, SparseMatrixCSC{T, I}}
    JH₂₂_sparse::Symmetric{T, SparseMatrixCSC{T, I}}
    JH_sparse::SparseMatrixCSC{T, I}
    ΔXⁿ::Vector{T}
    # --- Matrices and factorizations ---
    Q₁₁::Symmetric{T, SparseMatrixCSC{T, I}}
    τA_m₁xm₂::SparseMatrixCSC{T, I}
    τA_m₂xm₁::SparseMatrixCSC{T, I}
    M_m₁xm₁_vs2::Symmetric{T, SparseMatrixCSC{T, I}}
    M_m₂xm₂_vs2::Symmetric{T, SparseMatrixCSC{T, I}}
    M_m₃xm₃_chol::TF
    linsolve::TC
    # --- Vectors for JH_sparse synchronization ---
    map_direct₁₁::Vector{I}
    map_mirror₁₁::Vector{I}
    map_direct₂₂::Vector{I}
    map_mirror₂₂::Vector{I}
end

"""
    CrankNicolsonCache(matrices, dof_map_m₁, dof_map_m₂, dof_map_m₃, τ)

Constructor for [`CrankNicolsonCache`](@ref). Allocates all fields from the
global [`SystemMatrices`](@ref), the three DOF maps, and the time-step size `τ`.

# Arguments
- `matrices::SystemMatrices{T,I}`: global FEM matrices.
- `dof_map_m₁::DOFMap`: DOF map for the `m₁`-space (wave field).
- `dof_map_m₂::DOFMap`: DOF map for the `m₂`-space (heat field).
- `dof_map_m₃::DOFMap`: DOF map for the `m₃`-space (acoustic field).
- `τ::T`: time-step size; used to set the off-diagonal blocks of `Q` and `JH`.
"""
function CrankNicolsonCache(
        matrices::SystemMatrices{T, I},
        dof_map_m₁::DOFMap,
        dof_map_m₂::DOFMap,
        dof_map_m₃::DOFMap,
        τ::T
) where {T <: AbstractFloat, I <: Integer}
    m₁ = dof_map_m₁.m
    m₂ = dof_map_m₂.m
    m₃ = dof_map_m₃.m

    Xⁿ = zeros(T, m₁ + m₂)
    v̂ⁿ = view(Xⁿ, 1:m₁)
    ĉⁿ = view(Xⁿ, (m₁ + 1):(m₁ + m₂))

    M_m₁xm₁_vs2 = 2 * matrices.M_m₁xm₁
    M_m₂xm₂_vs2 = 2 * matrices.M_m₂xm₂

    τA_m₁xm₂ = τ * matrices.A_m₁xm₂
    τA_m₂xm₁ = τ * matrices.A_m₂xm₁

    Q₁₁ = similar(M_m₁xm₁_vs2)
    JH₁₁ = similar(M_m₁xm₁_vs2)
    JH₂₂_sparse = similar(M_m₂xm₂_vs2)

    JH_sparse = [M_m₁xm₁_vs2 τA_m₁xm₂;
                 τA_m₂xm₁ M_m₂xm₂_vs2]::SparseMatrixCSC{T, I}

    prob = LS.LinearProblem(JH_sparse, zeros(T, m₁ + m₂))
    linsolve = LS.init(
        prob, LS.KLUFactorization(; reuse_symbolic = true, check_pattern = false))

    map_direct₁₁, map_mirror₁₁ = build_maps_11(JH_sparse, JH₁₁.data)
    map_direct₂₂, map_mirror₂₂ = build_maps_22(JH_sparse, JH₂₂_sparse.data, τA_m₁xm₂)

    return CrankNicolsonCache(
        # Scratch vectors — m₁
        zeros(T, m₁), zeros(T, m₁), zeros(T, m₁),
        # Scratch vectors — m₂
        zeros(T, m₂), zeros(T, m₂),
        # Scratch vectors — m₃
        zeros(T, m₃), zeros(T, m₃),
        # RHS vectors
        zeros(T, m₁), zeros(T, m₂),
        # Midpoint unknowns and auxiliary quantities
        Xⁿ, v̂ⁿ, ĉⁿ, zeros(T, m₃), zeros(T, m₁), zeros(T, m₂), zeros(T, m₁ + m₂),
        # Newton system
        zeros(T, m₁ + m₂), JH₁₁, JH₂₂_sparse, JH_sparse, zeros(T, m₁ + m₂),
        # Matrices and factorizations
        Q₁₁,
        τA_m₁xm₂, τA_m₂xm₁,
        M_m₁xm₁_vs2, M_m₂xm₂_vs2,
        cholesky(matrices.M_m₃xm₃),
        linsolve,
        # map synchronization        
        map_direct₁₁, map_mirror₁₁, map_direct₂₂, map_mirror₂₂
    )
end

"""
    build_maps_11(J, M₁₁) -> (map_direct, map_mirror)

Build index maps from `M₁₁.nzval` to `J.nzval` positions for the (1,1) block.

See the test item "Synchronization problem - case study 2" for a step-by-step analysis.
"""
function build_maps_11(J, M₁₁)
    nnz_M₁₁ = nnz(M₁₁)
    map_direct = Vector{Int}(undef, nnz_M₁₁)
    map_mirror = Vector{Int}(undef, nnz_M₁₁)
    for j in 1:size(M₁₁, 2)        # Iterate over columns of M₁₁
        for kM in nzrange(M₁₁, j)  # Iterate over nonzeros of M₁₁[:,j]
            i = M₁₁.rowval[kM]
            for kJ in nzrange(J, j)# Iterate over nonzeros of J[:,j] and find position of J[i,j]
                if J.rowval[kJ] == i
                    map_direct[kM] = kJ
                    break
                end
            end
            for kᵀJ in nzrange(J, i)# Iterate over nonzeros of J[:,i] and find position of J[j,i]
                if J.rowval[kᵀJ] == j
                    map_mirror[kM] = kᵀJ
                    break
                end
            end
        end
    end
    return map_direct, map_mirror
end

"""
    build_maps_22(J, M₂₂, A₁₂) -> (map_direct, map_mirror)

Build index maps from `M₂₂.nzval` to `J.nzval` positions for the (2,2) block.

See the test item "Synchronization problem - case study 2" for a step-by-step analysis.
"""
function build_maps_22(J, M₂₂, A₁₂)
    m₁ = size(A₁₂, 1)

    # start[j]: index within nzrange(J, j+m₁) where the M̄₂₂ entries begin,
    # i.e. right after the A₁₂[:,j] entries.
    start = [length(nzrange(A₁₂, j)) + 1 for j in 1:size(A₁₂, 2)]

    nnz_M₂₂ = nnz(M₂₂)
    map_direct = Vector{Int}(undef, nnz_M₂₂)
    map_mirror = Vector{Int}(undef, nnz_M₂₂)
    for j in 1:size(M₂₂, 2)                         # Iterate over columns of M₂₂
        jJ = j + m₁
        for kM in nzrange(M₂₂, j)                   # Iterate over nonzeros of M₂₂[:,j]
            i = M₂₂.rowval[kM]
            iJ = i + m₁
            for kJ in nzrange(J, jJ)[start[j]:end]  # find J[i+m₁, j+m₁]
                if J.rowval[kJ] == iJ
                    map_direct[kM] = kJ
                    break
                end
            end
            for kᵀJ in nzrange(J, iJ)[start[i]:end] # find J[j+m₁, i+m₁]
                if J.rowval[kᵀJ] == jJ
                    map_mirror[kM] = kᵀJ
                    break
                end
            end
        end
    end
    return map_direct, map_mirror
end