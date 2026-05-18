function build_cache(::ModifiedCN, matrices, dof_map_m₁, dof_map_m₂, dof_map_m₃, τ)
    ModifiedCNCache(matrices, dof_map_m₁, dof_map_m₂, dof_map_m₃)
end

"""
    ModifiedCNCache{T, I, F}

Pre-allocated workspace for the Modified Crank–Nicolson time integrator.

# Fields
- `vec_m₁_1`–`vec_m₁_5`: scratch vectors of length `m₁`.
- `vec_m₂_1`, `vec_m₂_2`: scratch vectors of length `m₂`.
- `vec_m₃_1`, `vec_m₃_2`: scratch vectors of length `m₃`.
- `d_star`, `c_star`: extrapolated displacement ``d^{\\ast n}`` and temperature ``c^{\\ast n}``, lengths `m₁` and `m₂`.
- `v̂ⁿ`, `ĉⁿ`, `r̂ⁿ`: midpoint wave velocity, temperature, and acoustic unknown, lengths `m₁`, `m₂`, `m₃`.
- `dⁿ⁻²`, `cⁿ⁻²`: displacement and temperature two steps back (used for ``n \\geq 2`` extrapolation), lengths `m₁` and `m₂`.
- `L₁`, `L₂`: right-hand sides of the ``\\hat{v}^n`` and ``\\hat{c}^n`` systems, lengths `m₁` and `m₂`.
- `Q₁`, `Q₂`: LHS matrices for the ``\\hat{v}^n`` and ``\\hat{c}^n`` systems, sizes `m₁×m₁` and `m₂×m₂`; updated each time step.
- `M_m₁xm₁_vs2`, `M_m₂xm₂_vs2`: scaled mass matrices ``2M^{m_1\\times m_1}`` and ``2M^{m_2\\times m_2}``.
- `M_m₃xm₃_chol`: Cholesky factorization of ``M^{m_3\\times m_3}``.
- `minusH`: residual ``-H(\\hat{v}^n)``, length `m₁`.
- `JH`: Jacobian of ``H``, size `m₁×m₁`; assembled each Newton step.
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