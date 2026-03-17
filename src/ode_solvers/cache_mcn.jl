function build_cache(::ModifiedCN, matrices, dof_map_m₁, dof_map_m₂, dof_map_m₃, τ)
    ModifiedCNCache(matrices, dof_map_m₁, dof_map_m₂, dof_map_m₃)
end

"""
    ModifiedCNCache{T, I, F}

Pre-allocated workspace for the Modified Crank–Nicolson time integrator.
All fields are initialized once before the time loop. Intended to make
`ode_solve`, `perform_first_step!`, and `perform_step!` allocation-free,
but known allocations remain (see `FIXME` annotations in the solver source).

# Type parameters
- `T <: AbstractFloat`: floating-point precision (e.g. `Float64`).
- `I <: Integer`: integer type used by the sparse matrices (e.g. `Int64`).
- `F`: concrete type of the Cholesky factorization of `M_m₃×m₃`
  (inferred automatically by the constructor).

# Scratch vectors
| Field         | Length | Purpose                                    |
|:------------- |:------:|:-------------------------------------------|
| `vec_m₁_1–5`  | `m₁`   | General-purpose workspace for `m₁`-space.  |
| `vec_m₂_1–2`  | `m₂`   | General-purpose workspace for `m₂`-space.  |
| `vec_m₃_1–2`  | `m₃`   | General-purpose workspace for `m₃`-space.  |

# Extrapolated and midpoint values
| Field      | Length | Description                                                     |
|:---------- |:------:|:----------------------------------------------------------------|
| `d_star` | `m₁`   | Extrapolated wave displacement ``d^{\\ast n}``.                   |
| `c_star` | `m₂`   | Extrapolated temperature ``c^{\\ast n}``.                         |
| `v̂ⁿ`     | `m₁`   | Midpoint wave velocity ``\\hat{v}^n``; see [`solve_v̂ⁿ!`](@ref).   |
| `ĉⁿ`     | `m₂`   | Midpoint temperature ``\\hat{c}^n``; see [`solve_ĉⁿ!`](@ref).     |
| `r̂ⁿ`     | `m₃`   | Midpoint acoustic unknown ``\\hat{r}^n``; see [`solve_r̂ⁿ!`](@ref).|
| `dⁿ⁻²`   | `m₁`   | Wave displacement two steps back (for ``n \\geq 2``).             |
| `cⁿ⁻²`   | `m₂`   | Temperature two steps back (for ``n \\geq 2``).                   |

# RHS vectors
| Field | Length | Description                                                  |
|:----- |:------:|:-------------------------------------------------------------|
| `L₁`  | `m₁`   | RHS of the ``\\hat{v}^n`` system; see [`compute_L₁!`](@ref). |
| `L₂`  | `m₂`   | RHS of the ``\\hat{c}^n`` system; see [`compute_L₂!`](@ref). |

# Matrices and factorizations
| Field          | Size    | Description                                                                                                                                                                     |
|:---------------|:-------:|:----------------------------------------------------------------------------------|
| `Q₁`           | `m₁×m₁` | LHS matrix for the ``\\hat{v}^n`` system; see [`compute_Q₁!`](@ref).              |
| `Q₂`           | `m₂×m₂` | LHS matrix for the ``\\hat{c}^n`` system; see [`compute_Q₂!`](@ref).              |
| `M_m₁xm₁_vs2`  | `m₁×m₁` | ``2M^{m_1 \\times m_1}``; used in [`compute_Q₁!`](@ref) and [`compute_L₁!`](@ref).|
| `M_m₂xm₂_vs2`  | `m₂×m₂` | ``2M^{m_2 \\times m_2}``; used in [`compute_L₂!`](@ref).                          |
| `M_m₃xm₃_chol` | —       | Cholesky factorization of ``M^{m_3 \\times m_3}``; used in [`solve_r̂ⁿ!`](@ref).   |

# Newton system
| Field      | Size    | Description                                                             |
|:---------- |:-------:|:------------------------------------------------------------------------|
| `minusH` | `m₁`    | Residual ``-H(\\hat{v}^n)``; see [`compute_minusH!`](@ref).               |
| `JH`     | `m₁×m₁` | Jacobian of ``H``; assembled each Newton step; see [`compute_JH!`](@ref). |

"""
struct ModifiedCNCache{T <: AbstractFloat, I <: Integer, F}
    # --- Scratch vectors — m₁ ---
    vec_m₁_1::Vector{T}
    vec_m₁_2::Vector{T}
    vec_m₁_3::Vector{T}
    vec_m₁_4::Vector{T}
    vec_m₁_5::Vector{T}
    # --- Scratch vectors — m₂ ---
    vec_m₂_1::Vector{T}
    vec_m₂_2::Vector{T}
    # --- Scratch vectors — m₃ ---
    vec_m₃_1::Vector{T}
    vec_m₃_2::Vector{T}
    # --- Extrapolated values ---
    d_star::Vector{T}
    c_star::Vector{T}
    # --- Midpoint unknowns ---
    v̂ⁿ::Vector{T}
    ĉⁿ::Vector{T}
    r̂ⁿ::Vector{T}
    # --- Two-step-back storage (n ≥ 2 extrapolation) ---
    dⁿ⁻²::Vector{T}
    cⁿ⁻²::Vector{T}
    # --- RHS vectors ---
    L₁::Vector{T}
    L₂::Vector{T}
    # --- Matrices ---
    Q₁::Symmetric{T, SparseMatrixCSC{T, I}}
    Q₂::Symmetric{T, SparseMatrixCSC{T, I}}
    M_m₁xm₁_vs2::Symmetric{T, SparseMatrixCSC{T, I}}
    M_m₂xm₂_vs2::Symmetric{T, SparseMatrixCSC{T, I}}
    M_m₃xm₃_chol::F
    # --- Newton system ---
    minusH::Vector{T}
    JH::Symmetric{T, SparseMatrixCSC{T, I}}
end

"""
    ModifiedCNCache(matrices, dof_map_m₁, dof_map_m₂, dof_map_m₃)

Constructor for [`ModifiedCNCache`](@ref). Allocates all fields from the
global [`SystemMatrices`](@ref) and the three DOF maps.

# Arguments
- `matrices::SystemMatrices{T,I}`: global FEM matrices.
- `dof_map_m₁::DOFMap`: DOF map for the `m₁`-space (wave field).
- `dof_map_m₂::DOFMap`: DOF map for the `m₂`-space (heat field).
- `dof_map_m₃::DOFMap`: DOF map for the `m₃`-space (acoustic field).

# Notes
- `Q₁` and `JH` are initialized via `similar(matrices.M_m₁xm₁)`.
- `Q₂` is initialized via `similar(matrices.M_m₂xm₂)`.
- `M_m₁xm₁_vs2` and `M_m₂xm₂_vs2` are computed as `2M` at construction time.
- `M_m₃xm₃_chol` is the Cholesky factorization of `matrices.M_m₃xm₃`.
"""
function ModifiedCNCache(
        matrices::SystemMatrices{T, I},
        dof_map_m₁::DOFMap,
        dof_map_m₂::DOFMap,
        dof_map_m₃::DOFMap
) where {T <: AbstractFloat, I <: Integer}
    m₁ = dof_map_m₁.m
    m₂ = dof_map_m₂.m
    m₃ = dof_map_m₃.m

    return ModifiedCNCache(
        # Scratch vectors — m₁
        zeros(T, m₁), zeros(T, m₁), zeros(T, m₁), zeros(T, m₁), zeros(T, m₁),
        # Scratch vectors — m₂
        zeros(T, m₂), zeros(T, m₂),
        # Scratch vectors — m₃
        zeros(T, m₃), zeros(T, m₃),
        # Extrapolated values
        zeros(T, m₁), zeros(T, m₂),
        # Midpoint unknowns
        zeros(T, m₁), zeros(T, m₂), zeros(T, m₃),
        # Two-step-back storage
        zeros(T, m₁), zeros(T, m₂),
        # RHS vectors
        zeros(T, m₁), zeros(T, m₂),
        # Matrices
        similar(matrices.M_m₁xm₁),
        similar(matrices.M_m₂xm₂),
        2 * matrices.M_m₁xm₁,
        2 * matrices.M_m₂xm₂,
        cholesky(matrices.M_m₃xm₃),
        # Newton system
        zeros(T, m₁),
        similar(matrices.M_m₁xm₁)
    )
end