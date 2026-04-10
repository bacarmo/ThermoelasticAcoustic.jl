function build_cache(::CrankNicolson, matrices, dof_map_m₁, dof_map_m₂, dof_map_m₃, τ)
    CrankNicolsonCache(matrices, dof_map_m₁, dof_map_m₂, dof_map_m₃, τ)
end

"""
    CrankNicolsonCache{T, I, TF, TC}

Pre-allocated workspace for the Crank–Nicolson time integrator (Strategy 2:
joint Newton iteration in ``[\\hat{v}^n;\\,\\hat{c}^n]``).
All fields are initialized once before the time loop. Intended to make
`ode_solve` and `perform_step!` allocation-free, but known allocations remain
(see `FIXME` annotations in the solver source).

# Type parameters
- `T <: AbstractFloat`: floating-point precision (e.g. `Float64`).
- `I <: Integer`: integer type used by the sparse matrices (e.g. `Int64`).
- `TF`: concrete type of the Cholesky factorization of ``M^{m_3\\times m_3}``
  (inferred automatically by the constructor).
- `TC`: concrete type of linear solve cache.

# Scratch vectors
| Field         | Length | Purpose                                   |
|:------------- |:------:|:------------------------------------------|
| `vec_m₁_1–3`  | `m₁`   | General-purpose workspace for `m₁`-space. |
| `vec_m₂_1-2`  | `m₂`   | General-purpose workspace for `m₂`-space. |
| `vec_m₃_1–2`  | `m₃`   | General-purpose workspace for `m₃`-space. |

# RHS vectors
| Field | Length | Description                                                  |
|:----- |:------:|:-------------------------------------------------------------|
| `L₁`  | `m₁`   | RHS of the ``\\hat{v}^n`` system; see [`compute_L₁!`](@ref). |
| `L₂`  | `m₂`   | RHS of the ``\\hat{c}^n`` system; see [`compute_L₂!`](@ref). |

# Midpoint unknowns and auxiliary quantities
| Field           | Length | Description                                                                                                                                                  |
|:--------------- |:------:|:-------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `Xⁿ`            | `m₁+m₂`| Joint iterate ``[\\hat{v}^n;\\hat{c}^n]``.                                                                                                                   |
| `v̂ⁿ`            | `m₁`   | View into `Xⁿ[1:m₁]`; midpoint wave velocity ``\\hat{v}^n``.                                                                                                 |
| `ĉⁿ`            | `m₂`   | View into `Xⁿ[(m₁+1):end]`; midpoint temperature ``\\hat{c}^n``.                                                                                             |
| `r̂ⁿ`            | `m₃`   | Midpoint acoustic unknown ``\\hat{r}^n``; see [`solve_r̂ⁿ!`](@ref).                                                                                           |
| `d̂ⁿ`            | `m₁`   | Midpoint wave displacement ``\\hat{d}^n = \\tfrac{\\tau}{2}\\hat{v}^n + d^{n-1}``; computed in [`compute_minusH!`](@ref) and reused in [`compute_JH!`](@ref).|
| `K_m₂xm₂_vs_ĉⁿ` | `m₂`   | Product ``K^{m_2\\times m_2}\\hat{c}^n``; computed in [`compute_minusH!`](@ref) and reused in [`compute_JH!`](@ref).  |
| `QXⁿ`           | `m₁+m₂`| Product ``Q Xⁿ``; computed in [`compute_QXⁿ!`](@ref).|

# Newton system 
| Field        | Size             | Description                                            |
|:------------ |:----------------:|:-------------------------------------------------------|
| `minusH`     | `m₁+m₂`          | Residual ``-H(X^n)``; see [`compute_minusH!`](@ref).   |
| `JH₁₁`       | `m₁×m₁`          | First diagonal block of the Jacobian of ``H``.         |
| `JH₂₂_sparse`| `m₂×m₂`          | Second diagonal block of the Jacobian of ``H``, excluding the rank-1 matrix.|
| `JH_sparse`  | `(m₁+m₂)×(m₁+m₂)`| Global sparse Jacobian without the rank-1 matrix; see [`sync_JH_sparse!`](@ref). The off-diagonal blocks (``\\tau A^{m_1\\times m_2}``, ``\\tau A^{m_2\\times m_1}``) are time-invariant and set at construction; diagonal blocks are synced each Newton iteration via the pre-computed index maps. |
| `ΔXⁿ`        | `m₁+m₂`          | Newton linear system solution.|

# Matrices and factorizations
| Field          | Size      | Description                                                                       |
|:-------------- |:---------:|:----------------------------------------------------------------------------------|
| `Q₁₁`          | `m₁×m₁`   | Top-left block of ``Q``; updated each time step by [`compute_Q₁₁!`](@ref). The remaining three blocks of ``Q`` (``\\tau A^{m_1\\times m_2}``, ``\\tau A^{m_2\\times m_1}``, ``2M^{m_2\\times m_2}``) are time-invariant and available as separate cache fields.|
| `τA_m₁xm₂`     | `m₁×m₂`   | ``\\tau A^{m_1\\times m_2}`` |
| `τA_m₂xm₁`     | `m₂×m₁`   | ``\\tau A^{m_2\\times m_1}`` |
| `M_m₁xm₁_vs2`  | `m₁×m₁`   | ``2M^{m_1\\times m_1}``; used in [`compute_Q₁₁!`](@ref) and [`compute_L₁!`](@ref). |
| `M_m₂xm₂_vs2`  | `m₂×m₂`   | ``2M^{m_2\\times m_2}``; used in `Q` and [`compute_L₂!`](@ref). |
| `M_m₃xm₃_chol` | —         | Cholesky factorization of ``M^{m_3\\times m_3}``; see [`solve_r̂ⁿ!`](@ref).      |
| `linsolve`     | —         | LinearSolve.jl caching interface. |
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
        linsolve
    )
end
