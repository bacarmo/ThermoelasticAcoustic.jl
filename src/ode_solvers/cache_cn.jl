function build_cache(::CrankNicolson, matrices, dof_map_m₁, dof_map_m₂, dof_map_m₃, τ)
    CrankNicolsonCache(matrices, dof_map_m₁, dof_map_m₂, dof_map_m₃, τ)
end

"""
    CrankNicolsonCache{T, I, F, QT, JHT}

Pre-allocated workspace for the Crank–Nicolson time integrator (Strategy 2:
joint Newton iteration in ``[\\hat{v}^n;\\,\\hat{c}^n]``).
All fields are initialized once before the time loop. Intended to make
`ode_solve` and `perform_step!` allocation-free, but known allocations remain
(see `FIXME` annotations in the solver source).

# Type parameters
- `T <: AbstractFloat`: floating-point precision (e.g. `Float64`).
- `I <: Integer`: integer type used by the sparse matrices (e.g. `Int64`).
- `F`: concrete type of the Cholesky factorization of ``M^{m_3\\times m_3}``
  (inferred automatically by the constructor).
- `QT`, `JHT`: concrete `BlockMatrix` types (inferred automatically by the constructor).

# Scratch vectors
| Field         | Length | Purpose                                   |
|:------------- |:------:|:------------------------------------------|
| `vec_m₁_1–3`  | `m₁`   | General-purpose workspace for `m₁`-space. |
| `vec_m₂_1`    | `m₂`   | General-purpose workspace for `m₂`-space. |
| `vec_m₃_1–2`  | `m₃`   | General-purpose workspace for `m₃`-space. |

# Midpoint unknowns
| Field  | Length | Description                                                                                                                                                  |
|:------ |:------:|:-------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `v̂ⁿ`  | `m₁`    | View into `Xⁿ[1:m₁]`; midpoint wave velocity ``\\hat{v}^n``.                                                                                                  |
| `ĉⁿ`  | `m₂`    | View into `Xⁿ[(m₁+1):end]`; midpoint temperature ``\\hat{c}^n``.                                                                                              |
| `r̂ⁿ`  | `m₃`    | Midpoint acoustic unknown ``\\hat{r}^n``; see [`solve_r̂ⁿ!`](@ref).                                                                                            |
| `d̂ⁿ`  | `m₁`    | Midpoint wave displacement ``\\hat{d}^n = \\tfrac{\\tau}{2}\\hat{v}^n + d^{n-1}``; computed in [`compute_minusH!`](@ref) and reused in [`compute_JH!`](@ref). |

# Newton auxiliary quantities
| Field             | Size  | Description                                                                                                           |
|:----------------- |:-----:|:----------------------------------------------------------------------------------------------------------------------|
| `K_m₂xm₂_vs_ĉⁿ`   | `m₂`  | Product ``K^{m_2\\times m_2}\\hat{c}^n``; computed in [`compute_minusH!`](@ref) and reused in [`compute_JH!`](@ref).  |

# Newton system 
| Field    | Size               | Description                                                            |
|:-------- |:------------------:|:-----------------------------------------------------------------------|
| `Xⁿ`     | `m₁+m₂`            | `BlockedVector` (blocks `[m₁,m₂]`); joint iterate ``[\\hat{v}^n;\\hat{c}^n]``; updated in-place by Newton. |
| `minusH` | `m₁+m₂`            | `BlockedVector` (blocks `[m₁,m₂]`); residual ``-H(X^n)``; see [`compute_minusH!`](@ref).                   |
| `JH`     | `(m₁+m₂)×(m₁+m₂)`  | 2×2 `BlockMatrix`; jacobian of ``H``; see [`compute_JH!`](@ref). The two off-diagonal time-invariant blocks (``\\tau A^{m_1\\times m_2}``, ``\\tau A^{m_2\\times m_1}``) are set at construction and never modified; only diagonal blocks are updated each step. |

# RHS vectors
| Field | Length | Description                                                  |
|:----- |:------:|:-------------------------------------------------------------|
| `L₁`  | `m₁`   | RHS of the ``\\hat{v}^n`` system; see [`compute_L₁!`](@ref). |
| `L₂`  | `m₂`   | RHS of the ``\\hat{c}^n`` system; see [`compute_L₂!`](@ref). |

# Matrices and factorizations
| Field          | Size               | Description                                                                       |
|:-------------- |:------------------:|:----------------------------------------------------------------------------------|
| `Q`            | `(m₁+m₂)×(m₁+m₂)`  | 2×2 `BlockMatrix`; the three time-invariant blocks (``\\tau A^{m_1\\times m_2}``, ``\\tau A^{m_2\\times m_1}``, ``2M^{m_2\\times m_2}``) are set at construction and never modified; only the top-left ``m_1\\times m_1`` block is updated each step by [`compute_Q!`](@ref). |
| `M_m₁xm₁_vs2`  | `m₁×m₁`            | ``2M^{m_1\\times m_1}``; used in [`compute_Q!`](@ref) and [`compute_L₁!`](@ref). |
| `M_m₂xm₂_vs2`  | `m₂×m₂`            | ``2M^{m_2\\times m_2}``; used in [`compute_Q!`](@ref) and [`compute_L₂!`](@ref). |
| `M_m₃xm₃_chol` | —                  | Cholesky factorization of ``M^{m_3\\times m_3}``; see [`solve_r̂ⁿ!`](@ref).      |
"""
struct CrankNicolsonCache{T <: AbstractFloat, I <: Integer, F, QT, JHT}
    # --- Scratch vectors — m₁ ---
    vec_m₁_1::Vector{T}
    vec_m₁_2::Vector{T}
    vec_m₁_3::Vector{T}
    # --- Scratch vectors — m₂ ---
    vec_m₂_1::Vector{T}
    # --- Scratch vectors — m₃ ---
    vec_m₃_1::Vector{T}
    vec_m₃_2::Vector{T}
    # --- Midpoint unknowns ---
    v̂ⁿ::SubArray{T, 1, Vector{T}, Tuple{UnitRange{Int}}, true}
    ĉⁿ::SubArray{T, 1, Vector{T}, Tuple{UnitRange{Int}}, true}
    r̂ⁿ::Vector{T}
    d̂ⁿ::Vector{T}
    # --- Newton auxiliary quantities (computed in compute_minusH!, reused in compute_JH!) ---
    K_m₂xm₂_vs_ĉⁿ::Vector{T}
    # --- RHS vectors ---
    L₁::Vector{T}
    L₂::Vector{T}
    # --- Joint Newton iterate and residual ---
    Xⁿ::BlockedVector{T, Vector{T}, Tuple{BlockedOneTo{Int, Vector{Int}}}}
    minusH::BlockedVector{T, Vector{T}, Tuple{BlockedOneTo{Int, Vector{Int}}}}
    # --- Matrices ---
    Q::QT
    JH::JHT
    M_m₁xm₁_vs2::Symmetric{T, SparseMatrixCSC{T, I}}
    M_m₂xm₂_vs2::Symmetric{T, SparseMatrixCSC{T, I}}
    M_m₃xm₃_chol::F
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

    Xⁿ = BlockedVector{T}(undef, [m₁, m₂])
    v̂ⁿ = view(Xⁿ, Block(1))
    ĉⁿ = view(Xⁿ, Block(2))

    minusH = BlockedVector{T}(undef, [m₁, m₂])

    M_m₁xm₁_vs2 = 2 * matrices.M_m₁xm₁
    M_m₂xm₂_vs2 = 2 * matrices.M_m₂xm₂

    τA_m₁xm₂ = τ * matrices.A_m₁xm₂
    τA_m₂xm₁ = τ * matrices.A_m₂xm₁

    Q = BlockArray(undef_blocks, AbstractMatrix{T}, [m₁, m₂], [m₁, m₂])
    Q[Block(1, 1)] = deepcopy(matrices.M_m₁xm₁)
    Q[Block(1, 2)] = τA_m₁xm₂                    # Warning: This is not a copy
    Q[Block(2, 1)] = τA_m₂xm₁                    # Warning: This is not a copy
    Q[Block(2, 2)] = M_m₂xm₂_vs2                 # Warning: This is not a copy

    JH = BlockArray(undef_blocks, AbstractMatrix{T}, [m₁, m₂], [m₁, m₂])
    JH[Block(1, 1)] = deepcopy(matrices.M_m₁xm₁)
    JH[Block(1, 2)] = τA_m₁xm₂                   # Warning: This is not a copy
    JH[Block(2, 1)] = τA_m₂xm₁                   # Warning: This is not a copy
    JH[Block(2, 2)] = zeros(T, m₂, m₂)

    return CrankNicolsonCache(
        # Scratch vectors — m₁
        zeros(T, m₁), zeros(T, m₁), zeros(T, m₁),
        # Scratch vectors — m₂
        zeros(T, m₂),
        # Scratch vectors — m₃
        zeros(T, m₃), zeros(T, m₃),
        # Midpoint unknowns (v̂ⁿ and ĉⁿ are views into Xⁿ)
        v̂ⁿ, ĉⁿ, zeros(T, m₃), zeros(T, m₁),
        # Newton auxiliary quantities
        zeros(T, m₂),
        # RHS vectors
        zeros(T, m₁), zeros(T, m₂),
        # Joint Newton iterate and residual
        Xⁿ, minusH,
        # Matrices
        Q, JH, M_m₁xm₁_vs2, M_m₂xm₂_vs2,
        cholesky(matrices.M_m₃xm₃)
    )
end
