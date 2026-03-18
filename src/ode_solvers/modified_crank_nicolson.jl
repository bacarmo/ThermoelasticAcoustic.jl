function ode_solve(
        ::ModifiedCN,
        cache::ModifiedCNCache,
        state::FEMState{T},
        matrices::SystemMatrices,
        dof_map_m₁::DOFMap,
        dof_map_m₂::DOFMap,
        dof_map_m₃::DOFMap,
        mesh1D::CartesianMesh{1, I},
        mesh2D::CartesianMesh{2, I},
        quad::QuadratureSetup,
        tspan::StepRangeLen{T},
        input_data::PDEInputData,
        callback::AbstractCallback
) where {T <: AbstractFloat, I <: Integer}
    # ========================================
    # Pre-compute constants
    # ========================================
    τ = T(step(tspan))
    q₅ = 2 * input_data.q₁ + τ * input_data.q₂ + (τ^2 / 2) * input_data.q₃ # q₅ = 2q₁ + τq₂ + (τ²/2)q₃

    # ==========================================
    # Step "1,0" (predictor) + step n=1 (corrector)
    # ==========================================
    perform_first_step!(ModifiedCN(), cache, state, matrices,
        dof_map_m₁, dof_map_m₂, dof_map_m₃,
        mesh1D, mesh2D, quad, τ, q₅, input_data)
    apply!(callback, state,
        mesh1D, mesh2D, dof_map_m₁, dof_map_m₂, dof_map_m₃, quad, input_data)

    # ========================================
    # Steps n ≥ 2
    # ========================================
    for n in 2:(length(tspan) - 1)
        perform_step!(ModifiedCN(), cache, state, matrices,
            dof_map_m₁, dof_map_m₂, dof_map_m₃,
            mesh1D, mesh2D, quad, n, τ, q₅, input_data)
        apply!(callback, state,
            mesh1D, mesh2D, dof_map_m₁, dof_map_m₂, dof_map_m₃, quad, input_data)
    end

    return nothing
end

# ==============================================================================
# perform_first_step!
# ==============================================================================
"""
    perform_first_step!(::ModifiedCN, cache, state, matrices, dof_map_m₁, dof_map_m₂,
                        dof_map_m₃, mesh1D, mesh2D, quad, τ, q₅, input_data)

Advance `state` to the first time level ``t_1`` using a predictor–corrector
strategy specific to the Modified Crank–Nicolson scheme:

- **Predictor** (step ``n =`` "1,0"): uses ``w^{\\ast n} = w^0`` and solves for
  ``\\hat{v}^n`` and ``\\hat{c}^n`` only (``\\hat{r}^n`` is not needed in this sub-step).
- **Corrector** (step ``n = 1``): uses ``w^{\\ast n} = (w^{\\text{"1,0"}} + w^0)/2``
  and solves the full system, including ``\\hat{r}^n``.

On entry, `state` holds the initial condition ``(v^0, d^0, c^0, r^0, z^0)``.
On exit, `state` holds the solution at ``t_1``, and `cache.dⁿ⁻²`, `cache.cⁿ⁻²`
hold ``d^0``, ``c^0`` for use in subsequent two-step extrapolations (``n \\geq 2``).

# Arguments
- `::ModifiedCN`: dispatch token for the Modified Crank–Nicolson scheme.
- `cache::ModifiedCNCache`: pre-allocated workspace.
- `state::FEMState`: solution state at time level ``0``.
- `matrices::SystemMatrices`: global FEM matrices.
- `dof_map_m₁::DOFMap`: DOF map for the ``m_1``-space.
- `dof_map_m₂::DOFMap`: DOF map for the ``m_2``-space.
- `dof_map_m₃::DOFMap`: DOF map for the ``m_3``-space.
- `mesh1D::CartesianMesh{1}`: 1D Cartesian mesh (boundary ``\\Gamma_1``).
- `mesh2D::CartesianMesh{2}`: 2D Cartesian mesh (domain ``\\Omega``).
- `quad::QuadratureSetup`: precomputed quadrature data.
- `τ::T`: time-step size.
- `q₅::T`: scalar ``q_5 = 2q_1 + \\tau q_2 + (\\tau^2/2)q_3``.
- `input_data::PDEInputData`: problem parameters.
"""
function perform_first_step!(
        ::ModifiedCN,
        cache::ModifiedCNCache,
        state::FEMState{T},
        matrices::SystemMatrices,
        dof_map_m₁::DOFMap,
        dof_map_m₂::DOFMap,
        dof_map_m₃::DOFMap,
        mesh1D::CartesianMesh{1},
        mesh2D::CartesianMesh{2},
        quad::QuadratureSetup,
        τ::T,
        q₅::T,
        input_data::PDEInputData
) where {T <: AbstractFloat}
    t_half = τ / 2
    τ_2 = τ / 2

    # --- Predictor: d*¹ = d⁰,  c*¹ = c⁰  ---
    cache.d_star .= state.d
    cache.c_star .= state.c
    solve_v̂ⁿ!(cache, state, matrices, dof_map_m₁, dof_map_m₃,
        mesh1D, mesh2D, quad, τ, t_half, q₅, input_data)
    solve_ĉⁿ!(cache, state, matrices, dof_map_m₂, mesh2D, quad, τ, t_half, input_data)

    # --- Corrector: d*¹ = (τ/2)v̂^{"1,0"} + d⁰,  c*¹ = ĉ^{"1,0"} ---
    cache.d_star .= τ_2 .* cache.v̂ⁿ .+ state.d
    cache.c_star .= cache.ĉⁿ
    solve_v̂ⁿ!(cache, state, matrices, dof_map_m₁, dof_map_m₃,
        mesh1D, mesh2D, quad, τ, t_half, q₅, input_data)
    solve_ĉⁿ!(cache, state, matrices, dof_map_m₂, mesh2D, quad, τ, t_half, input_data)
    solve_r̂ⁿ!(cache, state, mesh1D, dof_map_m₃, quad, τ, t_half, q₅, input_data)

    # Store w⁰ before state is overwritten (needed for n ≥ 2 extrapolation)
    cache.dⁿ⁻² .= state.d
    cache.cⁿ⁻² .= state.c

    update_state!(state, cache, τ)

    return nothing
end

# ==============================================================================
# perform_step!
# ==============================================================================

"""
    perform_step!(::ModifiedCN, cache, state, matrices, dof_map_m₁, dof_map_m₂,
                  dof_map_m₃, mesh1D, mesh2D, quad, n, τ, q₅, input_data)

Advance `state` from time level ``n-1`` to ``n`` (``n \\geq 2``) using the
Modified Crank–Nicolson scheme with the two-step extrapolation
``w^{\\ast n} = (3w^{n-1} - w^{n-2})/2``.

On entry, `state` holds ``(v^{n-1}, d^{n-1}, c^{n-1}, r^{n-1}, z^{n-1})`` and
`cache.dⁿ⁻²`, `cache.cⁿ⁻²` hold ``d^{n-2}``, ``c^{n-2}``.
On exit, `state` holds the solution at ``t_n``, and `cache.dⁿ⁻²`, `cache.cⁿ⁻²`
are updated to ``d^{n-1}``, ``c^{n-1}``.

# Arguments
- `::ModifiedCN`: dispatch token for the Modified Crank–Nicolson scheme.
- `cache::ModifiedCNCache`: pre-allocated workspace.
- `state::FEMState`: solution state at time level ``n-1``.
- `matrices::SystemMatrices`: global FEM matrices.
- `dof_map_m₁::DOFMap`: DOF map for the ``m_1``-space.
- `dof_map_m₂::DOFMap`: DOF map for the ``m_2``-space.
- `dof_map_m₃::DOFMap`: DOF map for the ``m_3``-space.
- `mesh1D::CartesianMesh{1}`: 1D Cartesian mesh (boundary ``\\Gamma_1``).
- `mesh2D::CartesianMesh{2}`: 2D Cartesian mesh (domain ``\\Omega``).
- `quad::QuadratureSetup`: precomputed quadrature data.
- `n::Int`: current time index (``n \\geq 2``).
- `τ::T`: time-step size.
- `q₅::T`: scalar ``q_5 = 2q_1 + \\tau q_2 + (\\tau^2/2)q_3``.
- `input_data::PDEInputData`: problem parameters.
"""
function perform_step!(
        ::ModifiedCN,
        cache::ModifiedCNCache,
        state::FEMState{T},
        matrices::SystemMatrices,
        dof_map_m₁::DOFMap,
        dof_map_m₂::DOFMap,
        dof_map_m₃::DOFMap,
        mesh1D::CartesianMesh{1},
        mesh2D::CartesianMesh{2},
        quad::QuadratureSetup,
        n::Int,
        τ::T,
        q₅::T,
        input_data::PDEInputData
) where {T <: AbstractFloat}
    t_half = (n - T(0.5)) * τ

    cache.d_star .= T(1.5) .* state.d .- T(0.5) .* cache.dⁿ⁻²
    cache.c_star .= T(1.5) .* state.c .- T(0.5) .* cache.cⁿ⁻²

    solve_v̂ⁿ!(cache, state, matrices, dof_map_m₁, dof_map_m₃,
        mesh1D, mesh2D, quad, τ, t_half, q₅, input_data)
    solve_ĉⁿ!(cache, state, matrices, dof_map_m₂, mesh2D, quad, τ, t_half, input_data)
    solve_r̂ⁿ!(cache, state, mesh1D, dof_map_m₃, quad, τ, t_half, q₅, input_data)

    # Rotate: dⁿ⁻² ← dⁿ⁻¹,  cⁿ⁻² ← cⁿ⁻¹  (before state is overwritten)
    cache.dⁿ⁻² .= state.d
    cache.cⁿ⁻² .= state.c

    update_state!(state, cache, τ)

    return nothing
end

# ==============================================================================
# solve_v̂ⁿ!
# ==============================================================================

"""
    solve_v̂ⁿ!(cache, state, matrices, dof_map_m₁, dof_map_m₃, mesh1D, mesh2D,
               quad, τ, t_half, q₅, input_data; abstol, maxiter)

Solve the nonlinear system for the midpoint wave velocity ``\\hat{v}^n``,
storing the result in `cache.v̂ⁿ`. Newton's method is warm-started from `state.v`.

Sequentially:
1. Computes `Q₁` via [`compute_Q₁!`](@ref).
2. Computes ``L_1`` via [`compute_L₁!`](@ref).
3. Solves the nonlinear system ``Q_1 \\hat{v}^n + \\tau\\alpha G^{m_1}(\\hat{v}^n) = L_1``
   via [`newton_solve!`](@ref).

Assumes `cache.d_star`, `cache.c_star` have already been populated.

# Arguments
- `cache::ModifiedCNCache`: pre-allocated workspace.
- `state::FEMState`: solution state at time level ``n-1``.
- `matrices::SystemMatrices`: global FEM matrices.
- `dof_map_m₁::DOFMap`: DOF map for the ``m_1``-space.
- `dof_map_m₃::DOFMap`: DOF map for the ``m_3``-space.
- `mesh1D::CartesianMesh{1}`: 1D Cartesian mesh (boundary ``\\Gamma_1``).
- `mesh2D::CartesianMesh{2}`: 2D Cartesian mesh (domain ``\\Omega``).
- `quad::QuadratureSetup`: precomputed quadrature data.
- `τ::T`: time-step size.
- `t_half::T`: midpoint time ``t_{n-1/2}``.
- `q₅::T`: scalar ``q_5 = 2q_1 + \\tau q_2 + (\\tau^2/2)q_3``.
- `input_data::PDEInputData`: problem parameters and forcing data.

# Keyword Arguments
- `abstol::T`: absolute tolerance passed to [`newton_solve!`](@ref) (default: `T(1e-10)`).
- `maxiter::Int`: maximum Newton iterations (default: `10`).
"""
function solve_v̂ⁿ!(
        cache::ModifiedCNCache{T},
        state::FEMState{T},
        matrices::SystemMatrices,
        dof_map_m₁::DOFMap,
        dof_map_m₃::DOFMap,
        mesh1D::CartesianMesh{1},
        mesh2D::CartesianMesh{2},
        quad::QuadratureSetup,
        τ::T,
        t_half::T,
        q₅::T,
        input_data::PDEInputData;
        abstol::T = T(1e-10),
        maxiter::Int = 10
) where {T}
    α = input_data.α(t_half)
    τα = τ * α

    compute_Q₁!(cache, matrices, τ, α, input_data.q₄, q₅)
    compute_L₁!(cache, state, matrices, dof_map_m₁, dof_map_m₃,
        mesh1D, mesh2D, quad, τ, t_half, α, q₅, input_data)

    newton_solve!(cache, state, dof_map_m₃, mesh1D, quad, τα, input_data; abstol, maxiter)

    return nothing
end

# ==============================================================================
# compute_Q₁!
# ==============================================================================

"""
    compute_Q₁!(cache, matrices, τ, α, q₄, q₅)

Update `cache.Q₁` in-place with the matrix

```math
Q_1 = 2M^{m_1 \\times m_1}
  + \\frac{\\tau^2}{2}\\alpha K^{m_1 \\times m_1}
  + \\frac{\\tau^2 q_4}{q_5}\\alpha
    \\begin{bmatrix}
    M^{m_3\\times m_3}       & 0^{m_3\\times(m_1-m_3)}\\\\[5pt]
    0^{(m_1-m_3)\\times m_3} & 0^{(m_1-m_3)\\times(m_1-m_3)}
    \\end{bmatrix}.
```

Operates only on the stored nonzeros; 
no allocation occurs except for the ``m_3 \\times m_3`` top-left block update (see FIXME in source).

# Arguments

- `cache::ModifiedCNCache`: pre-allocated workspace; `cache.Q₁` is updated in-place.
- `matrices::SystemMatrices`: global FEM matrices.
- `τ::T`: time-step size.
- `α::T`: value of ``\\alpha(t_{n-1/2})``.
- `q₄::T`: problem parameter ``q_4``.
- `q₅::T`: scalar ``q_5 = 2q_1 + \\tau q_2 + (\\tau^2/2)q_3``.
"""
function compute_Q₁!(
        cache::ModifiedCNCache{T},
        matrices::SystemMatrices,
        τ::T,
        α::T,
        q₄::T,
        q₅::T
) where {T}
    m₃ = size(matrices.M_m₃xm₃, 1)
    cst1 = (τ^2 / 2) * α
    cst2 = (τ^2 * q₄ / q₅) * α

    @. cache.Q₁.data.nzval = cache.M_m₁xm₁_vs2.data.nzval +
                             cst1 * matrices.K_m₁xm₁.data.nzval
    @. cache.Q₁.data[1:m₃, 1:m₃] += cst2 * matrices.M_m₃xm₃.data  # FIXME: allocates

    return nothing
end

# ==============================================================================
# compute_L₁!
# ==============================================================================

"""
    compute_L₁!(cache, state, matrices, dof_map_m₁, dof_map_m₃, mesh1D, mesh2D,
                quad, τ, t_half, α, q₅, input_data)

Compute ``L_1`` and store it in `cache.L₁`:
```math
L_1 =
- \\tau F^{m_1}(d^{\\ast n})
- \\tau A^{m_1\\times m_2} c^{\\ast n}
+ 2M^{m_1\\times m_1}v^{n-1}
- \\tau\\alpha K^{m_1 \\times m_1}d^{n-1}
+ \\tau\\mathcal{F}^{m_1}(f_1^{n-1/2})
+ \\frac{\\tau\\alpha}{q_5}
\\begin{bmatrix}
    M^{m_3 \\times m_3}(2q_1 r^{n-1} - \\tau q_3 z^{n-1})
    + \\tau\\mathcal{F}^{m_3}(f_3^{n-1/2})
    \\\\[5pt]
    0^{(m_1-m_3)}
\\end{bmatrix}.
```
Assumes `cache.d_star` and `cache.c_star` have already been populated.

# Arguments

- `cache::ModifiedCNCache`: pre-allocated workspace; `cache.L₁` is updated in-place.
- `state::FEMState`: solution state at time level ``n-1``.
- `matrices::SystemMatrices`: global FEM matrices.
- `dof_map_m₁::DOFMap`: DOF map for the ``m_1``-space.
- `dof_map_m₃::DOFMap`: DOF map for the ``m_3``-space.
- `mesh1D::CartesianMesh{1}`: 1D Cartesian mesh (boundary ``\\Gamma_1``).
- `mesh2D::CartesianMesh{2}`: 2D Cartesian mesh (domain ``\\Omega``).
- `quad::QuadratureSetup`: precomputed quadrature data.
- `τ::T`: time-step size.
- `t_half::T`: midpoint time ``t_{n-1/2}``.
- `α::T`: value of ``\\alpha(t_{n-1/2})``.
- `q₅::T`: scalar ``q_5 = 2q_1 + \\tau q_2 + (\\tau^2/2)q_3``.
- `input_data::PDEInputData`: problem parameters.
"""
function compute_L₁!(
        cache::ModifiedCNCache{T},
        state::FEMState{T},
        matrices::SystemMatrices,
        dof_map_m₁::DOFMap,
        dof_map_m₃::DOFMap,
        mesh1D::CartesianMesh{1, I},
        mesh2D::CartesianMesh{2, I},
        quad::QuadratureSetup,
        τ::T,
        t_half::T,
        α::T,
        q₅::T,
        input_data::PDEInputData
) where {T <: AbstractFloat, I <: Integer}
    m₁ = dof_map_m₁.m
    m₃ = dof_map_m₃.m
    τα = τ * α
    τα_q₅ = τα / q₅
    τ²α_q₅ = τ * τα_q₅
    q₁2 = 2 * input_data.q₁
    τq₃ = τ * input_data.q₃

    # vec_m₁_1 ← τF(d_star)
    assembly_nonlinearity_F!(
        cache.vec_m₁_1, τ, input_data.f, cache.d_star, mesh2D, dof_map_m₁, quad)

    # vec_m₁_2 ← A_m₁xm₂ · c_star
    mul!(cache.vec_m₁_2, matrices.A_m₁xm₂, cache.c_star)

    # vec_m₁_3 ← 2M_m₁×m₁ · vⁿ⁻¹
    mul!(cache.vec_m₁_3, cache.M_m₁xm₁_vs2, state.v)

    # vec_m₁_4 ← K_m₁×m₁ · dⁿ⁻¹
    mul!(cache.vec_m₁_4, matrices.K_m₁xm₁, state.d)

    # vec_m₁_5 ← τF(f₁(t_half))
    scale_2d = τ * mesh2D.Δx[1] * mesh2D.Δx[2] / 4
    assembly_rhs_2d!(cache.vec_m₁_5, (x, y) -> input_data.f₁(x, y, t_half),
        scale_2d, quad.W_φP, mesh2D, dof_map_m₁, quad.xP, quad.yP)

    # vec_m₃_1 ← 2q₁rⁿ⁻¹ - τq₃zⁿ⁻¹
    @. cache.vec_m₃_1 = q₁2 * state.r - τq₃ * state.z

    # vec_m₃_2 ← M_m₃×m₃ · vec_m₃_1
    mul!(cache.vec_m₃_2, matrices.M_m₃xm₃, cache.vec_m₃_1)

    # vec_m₃_1 ← (τ²α/q₅)F(f₃(t_half))
    scale_1d = τ²α_q₅ * mesh1D.Δx[1] / 2
    assembly_rhs_1d!(cache.vec_m₃_1, x -> input_data.f₃(x, t_half),
        scale_1d, quad.W_ϕP, mesh1D, dof_map_m₃, quad.xP)

    # Final assembly
    @inbounds for i in 1:m₃
        cache.L₁[i] = -cache.vec_m₁_1[i] - τ * cache.vec_m₁_2[i] +
                      cache.vec_m₁_3[i] - τα * cache.vec_m₁_4[i] +
                      cache.vec_m₁_5[i] +
                      τα_q₅ * cache.vec_m₃_2[i] + cache.vec_m₃_1[i]
    end
    @inbounds for i in (m₃ + 1):m₁
        cache.L₁[i] = -cache.vec_m₁_1[i] - τ * cache.vec_m₁_2[i] +
                      cache.vec_m₁_3[i] - τα * cache.vec_m₁_4[i] +
                      cache.vec_m₁_5[i]
    end

    return nothing
end

# ==============================================================================
# newton_solve!
# ==============================================================================

"""
    newton_solve!(cache, state, dof_map_m₃, mesh1D, quad, τα, input_data; abstol, maxiter)

Solve the nonlinear system
```math
H(\\hat{v}^n) = Q₁\\hat{v}^n
+ \\tau\\alpha G^{m_1}(\\hat{v}^n)
- L_1 = 0
```
via Newton's method, updating `cache.v̂ⁿ` in-place.
At each iteration, `cache.minusH` stores ``-H(\\hat{v}^n)``, so the Newton
step reads
```math
\\hat{v}^n \\leftarrow \\hat{v}^n + JH^{-1}(-H).
```
Convergence is declared when ``\\max_i |H_i(\\hat{v}^n)| \\leq \\texttt{abstol}``.

Assumes the following cache fields have already been populated:
- `cache.Q₁`: via [`compute_Q₁!`](@ref).
- `cache.L₁`: via [`compute_L₁!`](@ref).

# Arguments
- `cache::ModifiedCNCache`: pre-allocated workspace; `cache.v̂ⁿ` is updated in-place.
- `state::FEMState`: solution state at time level ``n-1``; provides ``d^{n-1}``
- `dof_map_m₃::DOFMap`: DOF map for the ``m_3``-space.
- `mesh1D::CartesianMesh{1}`: 1D Cartesian mesh (boundary ``\\Gamma_1``).
- `quad::QuadratureSetup`: precomputed quadrature data.
- `τα::T`: product ``\\tau\\alpha(t_{n-1/2})``.
- `input_data::PDEInputData`: problem parameters.

# Keyword Arguments
- `abstol::T`: absolute tolerance on ``\\max_i |H_i|`` (default: `T(1e-10)`).
- `maxiter::Int`: maximum number of Newton iterations (default: `10`).
"""
function newton_solve!(
        cache::ModifiedCNCache{T},
        state::FEMState{T},
        dof_map_m₃::DOFMap,
        mesh1D::CartesianMesh{1},
        quad::QuadratureSetup,
        τα::T,
        input_data::PDEInputData;
        abstol::T = T(1e-10),
        maxiter::Int = 10
) where {T}
    # Warm start: v̂ⁿ ← vⁿ⁻¹
    copyto!(cache.v̂ⁿ, state.v)

    for _ in 1:maxiter
        # Compute -H(v̂ⁿ)  →  cache.minusH
        compute_minusH!(cache, dof_map_m₃, mesh1D, quad, τα, input_data)

        maximum(abs, cache.minusH) ≤ abstol && return nothing

        # Assemble JH = Q₁ + τα [JG 0; 0 0] → cache.JH
        compute_JH!(cache, dof_map_m₃, mesh1D, quad, τα, input_data)

        # Newton step: v̂ⁿ ← v̂ⁿ + JH⁻¹ (-H)
        cache.v̂ⁿ .+= cache.JH \ cache.minusH                        # FIXME: allocates
    end

    @warn "newton_solve! did not converge within $maxiter iterations " *
          "(max|H| = $(maximum(abs, cache.minusH)), abstol = $abstol)"
    return nothing
end

# ==============================================================================
# compute_minusH!
# ==============================================================================

"""
    compute_minusH!(cache, dof_map_m₃, mesh1D, quad, τα, input_data)

Compute ``-H(\\hat{v}^n)`` and store the result in `cache.minusH`:
```math
-H = - Q₁\\hat{v}^n
   - \\tau\\alpha G^{m_1}(\\hat{v}^n)
   + L_1.
```
Assumes `cache.Q₁` holds ``Q_1`` (via [`compute_Q₁!`](@ref)) and
`cache.L₁` holds ``L_1`` (via [`compute_L₁!`](@ref)).

# Arguments

- `cache::ModifiedCNCache`: pre-allocated workspace; `cache.minusH` is updated in-place.
- `dof_map_m₃::DOFMap`: DOF map for the ``m_3``-space.
- `mesh1D::CartesianMesh{1}`: 1D Cartesian mesh (boundary ``\\Gamma_1``).
- `quad::QuadratureSetup`: precomputed quadrature data.
- `τα::T`: product ``\\tau\\alpha(t_{n-1/2})``.
- `input_data::PDEInputData`: problem parameters.
"""
function compute_minusH!(
        cache::ModifiedCNCache{T},
        dof_map_m₃::DOFMap,
        mesh1D::CartesianMesh{1},
        quad::QuadratureSetup,
        τα::T,
        input_data::PDEInputData
) where {T}
    m₃ = dof_map_m₃.m
    m₁ = length(cache.v̂ⁿ)

    # Q₁ v̂ⁿ  →  cache.minusH
    mul!(cache.minusH, cache.Q₁, cache.v̂ⁿ)

    # τα G(v̂ⁿ_{1:m₃}) →  cache.vec_m₃_1
    assembly_nonlinearity_G!(cache.vec_m₃_1, τα, input_data.g,
        @view(cache.v̂ⁿ[1:m₃]), mesh1D, dof_map_m₃, quad)

    # Final assembly: -H = L₁ - Q₁v̂ⁿ - τα G(v̂ⁿ)
    @inbounds for i in 1:m₃
        cache.minusH[i] = cache.L₁[i] - cache.minusH[i] - cache.vec_m₃_1[i]
    end
    @inbounds for i in (m₃ + 1):m₁
        cache.minusH[i] = cache.L₁[i] - cache.minusH[i]
    end

    return nothing
end

# ==============================================================================
# compute_JH!
# ==============================================================================

"""
    compute_JH!(cache, dof_map_m₃, mesh1D, quad, τα, input_data)
 
Assemble the Jacobian of ``H`` and store the result in `cache.JH`:
```math
JH = Q_1 + \\tau\\alpha
\\begin{bmatrix}
\\big[JG(\\hat{v}^n_{1:m_3})\\big]^{m_3\\times m_3} & 0^{m_3 \\times (m_1 - m_3)} \\\\
0^{(m_1 - m_3) \\times m_3}  & 0^{(m_1-m_3)\\times(m_1-m_3)}
\\end{bmatrix},
```
where ``JG`` is the ``m_3 \\times m_3`` Jacobian of the boundary
nonlinearity ``G^{m_1}``.
 
Assumes `cache.Q₁` and `cache.v̂ⁿ` have already been populated.
 
# Arguments
- `cache::ModifiedCNCache`: pre-allocated workspace; `cache.JH` is updated in-place.
- `dof_map_m₃::DOFMap`: DOF map for the ``m_3``-space.
- `mesh1D::CartesianMesh{1}`: 1D Cartesian mesh (boundary ``\\Gamma_1``).
- `quad::QuadratureSetup`: precomputed quadrature data.
- `τα::T`: product ``\\tau\\alpha(t_{n-1/2})``.
- `input_data::PDEInputData`: problem parameters; `input_data.∂ₛg` is the
  derivative of the boundary nonlinearity.
"""
function compute_JH!(
        cache::ModifiedCNCache{T},
        dof_map_m₃::DOFMap,
        mesh1D::CartesianMesh{1},
        quad::QuadratureSetup,
        τα::T,
        input_data::PDEInputData
) where {T}
    m₃ = dof_map_m₃.m

    # JG (m₃×m₃), scaled by τα                          # FIXME: allocates
    JG = assembly_global_matrix_DG(
        τα, input_data.∂ₛg, @view(cache.v̂ⁿ[1:m₃]), mesh1D, dof_map_m₃, quad)

    # JH ← Q₁, then embed τα JG into the top-left m₃×m₃ block
    cache.JH.data.nzval .= cache.Q₁.data.nzval
    @. cache.JH.data[1:m₃, 1:m₃] += JG.data             # FIXME: allocates

    return nothing
end

# ==============================================================================
# solve_ĉⁿ!
# ==============================================================================

"""
    solve_ĉⁿ!(cache, state, matrices, dof_map_m₂, mesh2D, quad, τ, t_half, input_data)

Solve the linear system for the midpoint temperature ``\\hat{c}^n``,
storing the result in `cache.ĉⁿ`:
```math
\\left[2M^{m_2 \\times m_2}
+ \\tau\\beta\\!\\left(\\mathbf{b}\\cdot c^{\\ast n}\\right)
  K^{m_2 \\times m_2}\\right]\\hat{c}^n = L_2,
```

where
```math
L_2 = - \\tau A^{m_2 \\times m_1}\\hat{v}^n
      + 2M^{m_2 \\times m_2}c^{n-1}
      + \\tau\\mathcal{F}^{m_2}(f_2^{n-1/2}).
```

Assumes `cache.c_star` and `cache.v̂ⁿ` have already been populated.

# Arguments

- `cache::ModifiedCNCache`: pre-allocated workspace; `cache.ĉⁿ` is updated in-place.
- `state::FEMState`: solution state at time level ``n-1``.
- `matrices::SystemMatrices`: global FEM matrices.
- `dof_map_m₂::DOFMap`: DOF map for the ``m_2``-space.
- `mesh2D::CartesianMesh{2}`: 2D Cartesian mesh (domain ``\\Omega``).
- `quad::QuadratureSetup`: precomputed quadrature data.
- `τ::T`: time-step size.
- `t_half::T`: midpoint time ``t_{n-1/2}``.
- `input_data::PDEInputData`: problem parameters and forcing data.
"""
function solve_ĉⁿ!(
        cache::ModifiedCNCache{T},
        state::FEMState{T},
        matrices::SystemMatrices,
        dof_map_m₂::DOFMap,
        mesh2D::CartesianMesh{2},
        quad::QuadratureSetup,
        τ::T,
        t_half::T,
        input_data::PDEInputData
) where {T}
    compute_Q₂!(cache, matrices, τ, input_data)
    compute_L₂!(cache, state, matrices, dof_map_m₂, mesh2D, quad, τ, t_half, input_data)
    cache.ĉⁿ .= cache.Q₂ \ cache.L₂  # FIXME: allocates

    return nothing
end

# ==============================================================================
# compute_Q₂!
# ==============================================================================

"""
    compute_Q₂!(cache, matrices, τ, input_data)

Update `cache.Q₂` in-place with the matrix
```math
Q_2 = 2M^{m_2 \\times m_2}
    + \\tau\\beta\\!\\left(\\mathbf{b}\\cdot c^{\\ast n}\\right)
      K^{m_2 \\times m_2},
```

where
```math
b_i = \\int_\\Omega \\psi_i \\, d\\Omega
```

with ``\\psi_i`` the ``i``-th basis function of the ``m_2``-space.

Assumes `cache.c_star` has already been populated.

# Arguments

- `cache::ModifiedCNCache`: pre-allocated workspace; `cache.Q₂` is updated in-place.
- `matrices::SystemMatrices`: global FEM matrices.
- `τ::T`: time-step size.
- `input_data::PDEInputData`: problem parameters; `input_data.β` is evaluated at
    ``\\mathbf{b} \\cdot c^{\\ast n}``.
"""
function compute_Q₂!(
        cache::ModifiedCNCache{T},
        matrices::SystemMatrices,
        τ::T,
        input_data::PDEInputData
) where {T <: AbstractFloat}
    cst1 = dot(matrices.b, cache.c_star)
    cst2 = input_data.β(cst1)
    cst3 = τ * cst2
    @. cache.Q₂.data.nzval = cache.M_m₂xm₂_vs2.data.nzval +
                             cst3 * matrices.K_m₂xm₂.data.nzval

    return nothing
end

# ==============================================================================
# compute_L₂!
# ==============================================================================

"""
    compute_L₂!(cache, state, matrices, dof_map_m₂, mesh2D, quad, τ, t_half, input_data)

Assemble ``L_2`` and store the result in `cache.L₂`:
```math
L_2 = - \\tau A^{m_2 \\times m_1}\\hat{v}^n
      + 2M^{m_2 \\times m_2}c^{n-1}
      + \\tau\\mathcal{F}^{m_2}(f_2^{n-1/2}).
```

Assumes `cache.v̂ⁿ` has already been populated.

# Arguments

- `cache::ModifiedCNCache`: pre-allocated workspace; `cache.L₂` is updated in-place.
- `state::FEMState`: solution state at time level ``n-1``.
- `matrices::SystemMatrices`: global FEM matrices.
- `dof_map_m₂::DOFMap`: DOF map for the ``m_2``-space.
- `mesh2D::CartesianMesh{2}`: 2D Cartesian mesh (domain ``\\Omega``).
- `quad::QuadratureSetup`: precomputed quadrature data.
- `τ::T`: time-step size.
- `t_half::T`: midpoint time ``t_{n-1/2}``.
- `input_data::PDEInputData`: problem parameters and forcing data.
"""
function compute_L₂!(
        cache::ModifiedCNCache{T},
        state::FEMState{T},
        matrices::SystemMatrices,
        dof_map_m₂::DOFMap,
        mesh2D::CartesianMesh{2},
        quad::QuadratureSetup,
        τ::T,
        t_half::T,
        input_data::PDEInputData
) where {T <: AbstractFloat}
    # A_m₂xm₁ · v̂ⁿ  →  cache.vec_m₂_1
    mul!(cache.vec_m₂_1, matrices.A_m₂xm₁, cache.v̂ⁿ)

    # 2M_m₂×m₂ · cⁿ⁻¹  →  cache.vec_m₂_2
    mul!(cache.vec_m₂_2, cache.M_m₂xm₂_vs2, state.c)

    # τ F(f₂(t_half))  →  cache.L₂
    scale = τ * mesh2D.Δx[1] * mesh2D.Δx[2] / 4
    assembly_rhs_2d!(cache.L₂, (x, y) -> input_data.f₂(x, y, t_half),
        scale, quad.W_φP, mesh2D, dof_map_m₂, quad.xP, quad.yP)

    # Final assembly
    @. cache.L₂ += cache.vec_m₂_2 - τ * cache.vec_m₂_1

    return nothing
end

# ==============================================================================
# solve_r̂ⁿ
# ==============================================================================
"""
    solve_r̂ⁿ!(cache, state, mesh1D, dof_map_m₃, quad, τ, t_half, q₅, input_data)

Recover the midpoint acoustic unknown ``\\hat{r}^n`` via the closed-form expression
```math
\\hat{r}^n =
- \\frac{\\tau q_4}{q_5}\\hat{v}^n_{1:m_3}
+ \\frac{2q_1}{q_5}r^{n-1}
- \\frac{\\tau q_3}{q_5}z^{n-1}
+ \\frac{\\tau}{q_5}\\left(M^{m_3 \\times m_3}\\right)^{-1}
  \\mathcal{F}^{m_3}(f_3^{n-\\frac{1}{2}}),
```
storing the result in `cache.r̂ⁿ`. The solve against ``M^{m_3 \\times m_3}``
uses the pre-computed Cholesky factorization `cache.M_m₃xm₃_chol`.

Assumes `cache.v̂ⁿ` have already been populated.

# Arguments
- `cache::ModifiedCNCache`: pre-allocated workspace; `cache.r̂ⁿ` is updated in-place.
- `state::FEMState`: solution state at time level ``n-1``.
- `mesh1D::CartesianMesh{1}`: 1D Cartesian mesh (boundary ``\\Gamma_1``).
- `dof_map_m₃::DOFMap`: DOF map for the ``m_3``-space.
- `quad::QuadratureSetup`: precomputed quadrature data.
- `τ::T`: time-step size.
- `t_half::T`: midpoint time ``t_{n-1/2}``.
- `q₅::T`: scalar ``q_5 = 2q_1 + \\tau q_2 + (\\tau^2/2)q_3``.
- `input_data::PDEInputData`: problem parameters.
"""
function solve_r̂ⁿ!(
        cache::ModifiedCNCache{T},
        state::FEMState{T},
        mesh1D::CartesianMesh{1},
        dof_map_m₃::DOFMap,
        quad::QuadratureSetup,
        τ::T,
        t_half::T,
        q₅::T,
        input_data::PDEInputData
) where {T}
    q₁ = input_data.q₁
    q₃ = input_data.q₃
    q₄ = input_data.q₄

    cst1 = -(τ * q₄ / q₅)
    cst2 = 2 * q₁ / q₅
    cst3 = τ * q₃ / q₅
    cst4 = τ / q₅

    m₃ = dof_map_m₃.m

    # F(f₃(t_half)) → cache.vec_m₃_1
    scale = mesh1D.Δx[1] / 2
    assembly_rhs_1d!(
        cache.vec_m₃_1, x -> input_data.f₃(x, t_half), scale, quad.W_ϕP, mesh1D, dof_map_m₃, quad.xP)

    # M_m₃⁻¹ F_m₃  →  cache.vec_m₃_2
    ldiv!(cache.vec_m₃_2, cache.M_m₃xm₃_chol, cache.vec_m₃_1)

    @inbounds for i in 1:m₃
        cache.r̂ⁿ[i] = cst1 * cache.v̂ⁿ[i] +
                       cst2 * state.r[i] -
                       cst3 * state.z[i] +
                       cst4 * cache.vec_m₃_2[i]
    end

    return nothing
end

# ==============================================================================
# update_state!
# ==============================================================================

"""
    update_state!(state, cache, τ)

Advance `state` from time level ``n-1`` to ``n`` using the midpoint unknowns stored in `cache`:
```math
\\begin{aligned}
v^n &= 2\\hat{v}^n - v^{n-1}, \\\\
c^n &= 2\\hat{c}^n - c^{n-1}, \\\\
r^n &= 2\\hat{r}^n - r^{n-1}, \\\\
d^n &= \\tau\\hat{v}^n + d^{n-1}, \\\\
z^n &= \\tau\\hat{r}^n + z^{n-1}.
\\end{aligned}
```
The time index `state.n` and current time `state.t` are also incremented by ``1`` and ``\\tau``, respectively.

Assumes `cache.v̂ⁿ`, `cache.ĉⁿ`, and `cache.r̂ⁿ` have already been populated
by [`solve_v̂ⁿ!`](@ref), [`solve_ĉⁿ!`](@ref), and [`solve_r̂ⁿ!`](@ref),
respectively.

# Arguments
- `state::FEMState`: solution state at time level ``n-1``; updated in-place to level ``n``.
- `cache::ModifiedCNCache`: pre-allocated workspace.
- `τ::T`: time-step size.
"""
function update_state!(
        state::FEMState{T, V},
        cache::ModifiedCNCache{T},
        τ::T
) where {T <: AbstractFloat, V <: AbstractVector{T}}
    state.n += 1
    state.t += τ
    @. state.v = 2 * cache.v̂ⁿ - state.v
    @. state.c = 2 * cache.ĉⁿ - state.c
    @. state.r = 2 * cache.r̂ⁿ - state.r
    @. state.d = τ * cache.v̂ⁿ + state.d
    @. state.z = τ * cache.r̂ⁿ + state.z

    return nothing
end