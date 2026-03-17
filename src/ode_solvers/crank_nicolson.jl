function ode_solve(
        ::CrankNicolson,
        cache::CrankNicolsonCache,
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
    q₅ = 2 * input_data.q₁ + τ * input_data.q₂ + (τ^2 / 2) * input_data.q₃

    # ========================================
    # Steps n ≥ 1
    # ========================================
    for n in 1:(length(tspan) - 1)
        perform_step!(CrankNicolson(), cache, state, matrices,
            dof_map_m₁, dof_map_m₂, dof_map_m₃,
            mesh1D, mesh2D, quad, n, τ, q₅, input_data)
        apply!(callback, state,
            mesh1D, mesh2D, dof_map_m₁, dof_map_m₂, dof_map_m₃, quad, input_data)
    end

    return nothing
end

# ==============================================================================
# perform_step!
# ==============================================================================

"""
    perform_step!(::CrankNicolson, cache, state, matrices, dof_map_m₁,
                  dof_map_m₂, dof_map_m₃, mesh1D, mesh2D, quad, n, τ, q₅,
                  input_data)
 
Advance `state` from time level ``n-1`` to ``n`` using the Crank–Nicolson
scheme (Strategy 2: joint Newton iteration in ``[\\hat{v}^n;\\,\\hat{c}^n]``).
 
On entry, `state` holds ``(v^{n-1}, d^{n-1}, c^{n-1}, r^{n-1}, z^{n-1})``.
On exit, `state` holds the solution at ``t_n``.
 
# Arguments
- `::CrankNicolson`: dispatch token.
- `cache::CrankNicolsonCache`: pre-allocated workspace.
- `state::FEMState`: solution state at time level ``n-1``.
- `matrices::SystemMatrices`: global FEM matrices.
- `dof_map_m₁::DOFMap`: DOF map for the ``m_1``-space.
- `dof_map_m₂::DOFMap`: DOF map for the ``m_2``-space.
- `dof_map_m₃::DOFMap`: DOF map for the ``m_3``-space.
- `mesh1D::CartesianMesh{1}`: 1D Cartesian mesh (boundary ``\\Gamma_1``).
- `mesh2D::CartesianMesh{2}`: 2D Cartesian mesh (domain ``\\Omega``).
- `quad::QuadratureSetup`: precomputed quadrature data.
- `n::Int`: current time index (``n \\geq 1``).
- `τ::T`: time-step size.
- `q₅::T`: scalar ``q_5 = 2q_1 + \\tau q_2 + (\\tau^2/2)q_3``.
- `input_data::PDEInputData`: problem parameters.
"""
function perform_step!(
        ::CrankNicolson,
        cache::CrankNicolsonCache,
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
    t_half = (n - T(0.5)) * τ   # t_{n-1/2}
    α = input_data.α(t_half)
    τα = τ * α

    compute_Q!(cache, matrices, τ, α, input_data.q₄, q₅)
    compute_L₁!(cache, state, matrices, dof_map_m₁, dof_map_m₃,
        mesh1D, mesh2D, quad, τ, t_half, α, q₅, input_data)
    compute_L₂!(cache, state, dof_map_m₂, mesh2D, quad, τ, t_half, input_data)

    newton_solve!(cache, state, matrices, dof_map_m₁, dof_map_m₂, dof_map_m₃,
        mesh1D, mesh2D, quad, τ, τα, input_data)

    solve_r̂ⁿ!(cache, state, mesh1D, dof_map_m₃, quad, τ, t_half, q₅, input_data)

    update_state!(state, cache, τ)

    return nothing
end

# ==============================================================================
# compute_Q!
# ==============================================================================

"""
    compute_Q!(cache, matrices, τ, α, q₄, q₅)
 
Update the top-left ``m_1 \\times m_1`` block of `cache.Q` in-place:
```math
Q_{1:m_1,1:m_1} =
2M^{m_1 \\times m_1}
+ \\frac{\\tau^2}{2}\\alpha K^{m_1 \\times m_1}
+ \\frac{\\tau^2 q_4}{q_5}\\alpha
  \\begin{bmatrix}
  M^{m_3\\times m_3}       & 0^{m_3\\times(m_1-m_3)}\\\\[5pt]
  0^{(m_1-m_3)\\times m_3} & 0^{(m_1-m_3)\\times(m_1-m_3)}
  \\end{bmatrix}.
```
 
The remaining three blocks of `cache.Q` are time-invariant and were
pre-filled at construction time (see [`CrankNicolsonCache`](@ref)):
- top-right: ``\\tau A^{m_1 \\times m_2}``,
- bottom-left: ``\\tau A^{m_2 \\times m_1}``,
- bottom-right: ``2M^{m_2 \\times m_2}``.
  
# Arguments
- `cache::CrankNicolsonCache`: pre-allocated workspace; `cache.Q` is updated in-place.
- `matrices::SystemMatrices`: global FEM matrices.
- `τ::T`: time-step size.
- `α::T`: value of ``\\alpha(t_{n-1/2})``.
- `q₄::T`: problem parameter ``q_4``.
- `q₅::T`: scalar ``q_5 = 2q_1 + \\tau q_2 + (\\tau^2/2)q_3``.
"""
function compute_Q!(
        cache::CrankNicolsonCache{T},
        matrices::SystemMatrices,
        τ::T,
        α::T,
        q₄::T,
        q₅::T
) where {T}
    m₁ = size(matrices.M_m₁xm₁, 1)
    m₃ = size(matrices.M_m₃xm₃, 1)
    cst1 = (τ^2 / 2) * α
    cst2 = (τ^2 * q₄ / q₅) * α

    # Top-left m₁×m₁ block: 2M_m₁xm₁ + cst1 * K_m₁xm₁
    @. cache.Q[1:m₁, 1:m₁] = cache.M_m₁xm₁_vs2 +
                             cst1 * matrices.K_m₁xm₁       # FIXME: allocates

    # Top-left m₃×m₃ sub-block: += cst2 * M_m₃xm₃
    @. cache.Q[1:m₃, 1:m₃] += cst2 * matrices.M_m₃xm₃      # FIXME: allocates

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
  2M^{m_1 \\times m_1}v^{n-1}
- \\tau\\alpha K^{m_1 \\times m_1}d^{n-1}
+ \\tau\\mathcal{F}^{m_1}(f_1^{n-1/2})
+ \\frac{\\tau\\alpha}{q_5}
\\begin{bmatrix}
    M^{m_3 \\times m_3}\\bigl(2q_1 r^{n-1} - \\tau q_3 z^{n-1}\\bigr)
    + \\tau\\mathcal{F}^{m_3}(f_3^{n-1/2})
    \\\\[5pt]
    0^{(m_1-m_3)}
\\end{bmatrix}.
```
 
# Arguments
- `cache::CrankNicolsonCache`: pre-allocated workspace; `cache.L₁` is updated in-place.
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
        cache::CrankNicolsonCache{T},
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

    # vec_m₁_1 ← 2M_m₁×m₁ · vⁿ⁻¹
    mul!(cache.vec_m₁_1, cache.M_m₁xm₁_vs2, state.v)

    # vec_m₁_2 ← K_m₁×m₁ · dⁿ⁻¹
    mul!(cache.vec_m₁_2, matrices.K_m₁xm₁, state.d)

    # vec_m₁_3 ← τ F(f₁(t_half))
    scale_2d = τ * mesh2D.Δx[1] * mesh2D.Δx[2] / 4
    assembly_rhs_2d!(cache.vec_m₁_3, (x, y) -> input_data.f₁(x, y, t_half),
        scale_2d, quad.W_φP, mesh2D, dof_map_m₁, quad.xP, quad.yP)

    # vec_m₃_1 ← 2q₁rⁿ⁻¹ - τq₃zⁿ⁻¹
    @. cache.vec_m₃_1 = q₁2 * state.r - τq₃ * state.z

    # vec_m₃_2 ← M_m₃×m₃ · vec_m₃_1
    mul!(cache.vec_m₃_2, matrices.M_m₃xm₃, cache.vec_m₃_1)

    # vec_m₃_1 ← (τ²α/q₅) F(f₃(t_half))
    scale_1d = τ²α_q₅ * mesh1D.Δx[1] / 2
    assembly_rhs_1d!(cache.vec_m₃_1, x -> input_data.f₃(x, t_half),
        scale_1d, quad.W_ϕP, mesh1D, dof_map_m₃, quad.xP)

    # Final assembly
    @inbounds for i in 1:m₃
        cache.L₁[i] = cache.vec_m₁_1[i] - τα * cache.vec_m₁_2[i] +
                      cache.vec_m₁_3[i] +
                      τα_q₅ * cache.vec_m₃_2[i] + cache.vec_m₃_1[i]
    end
    @inbounds for i in (m₃ + 1):m₁
        cache.L₁[i] = cache.vec_m₁_1[i] - τα * cache.vec_m₁_2[i] +
                      cache.vec_m₁_3[i]
    end

    return nothing
end

# ==============================================================================
# compute_L₂!
# ==============================================================================

"""
    compute_L₂!(cache, state, matrices, dof_map_m₂, mesh2D, quad, τ, t_half, input_data)
 
Assemble ``L_2`` and store the result in `cache.L₂`:
```math
L_2 = 2M^{m_2 \\times m_2}c^{n-1} + \\tau\\mathcal{F}^{m_2}(f_2^{n-1/2}).
```
 
# Arguments
- `cache::CrankNicolsonCache`: pre-allocated workspace; `cache.L₂` is updated in-place.
- `state::FEMState`: solution state at time level ``n-1``.
- `dof_map_m₂::DOFMap`: DOF map for the ``m_2``-space.
- `mesh2D::CartesianMesh{2}`: 2D Cartesian mesh (domain ``\\Omega``).
- `quad::QuadratureSetup`: precomputed quadrature data.
- `τ::T`: time-step size.
- `t_half::T`: midpoint time ``t_{n-1/2}``.
- `input_data::PDEInputData`: problem parameters.
"""
function compute_L₂!(
        cache::CrankNicolsonCache{T},
        state::FEMState{T},
        dof_map_m₂::DOFMap,
        mesh2D::CartesianMesh{2, I},
        quad::QuadratureSetup,
        τ::T,
        t_half::T,
        input_data::PDEInputData
) where {T <: AbstractFloat, I <: Integer}
    # vec_m₂_1 ← 2M_m₂×m₂ · cⁿ⁻¹
    mul!(cache.vec_m₂_1, cache.M_m₂xm₂_vs2, state.c)

    # L₂ ← τ F(f₂(t_half))
    scale = τ * mesh2D.Δx[1] * mesh2D.Δx[2] / 4
    assembly_rhs_2d!(cache.L₂, (x, y) -> input_data.f₂(x, y, t_half),
        scale, quad.W_φP, mesh2D, dof_map_m₂, quad.xP, quad.yP)

    # Final assembly
    @. cache.L₂ += cache.vec_m₂_1

    return nothing
end

# ==============================================================================
# newton_solve!
# ==============================================================================

"""
    newton_solve!(cache, state, matrices, dof_map_m₁, dof_map_m₂, dof_map_m₃,
                  mesh1D, mesh2D, quad, τ, τα, input_data; abstol, maxiter)
 
Solve the joint nonlinear system
```math
H(X) = Q X
+ \\begin{bmatrix}
    \\tau\\alpha G^{m_1}(\\hat{v}^n)
    + \\tau F^{m_1}\\bigl(\\tfrac{\\tau}{2}\\hat{v}^n + d^{n-1}\\bigr)
    \\\\[5pt]
    \\tau\\beta(\\mathbf{b}\\cdot\\hat{c}^n)K^{m_2\\times m_2}\\hat{c}^n
  \\end{bmatrix}
- \\begin{bmatrix} L_1 \\\\ L_2 \\end{bmatrix}
= 0
```
via Newton's method, updating `cache.Xⁿ` (and its aliases `cache.v̂ⁿ`, `cache.ĉⁿ`) in-place.
 
Convergence is declared when ``\\max_i |H_i(X)| \\leq \\texttt{abstol}``.
 
Assumes `cache.Q`, `cache.L₁`, and `cache.L₂` have already been populated.
`cache.Xⁿ` is warm-started from ``[v^{n-1};\\,c^{n-1}]`` at the beginning of each call.

# Arguments
- `cache::CrankNicolsonCache`: pre-allocated workspace.
- `state::FEMState`: solution state at time level ``n-1``; provides ``d^{n-1}``
  to evaluate ``F^{m_1}`` at the midpoint displacement.
- `matrices::SystemMatrices`: global FEM matrices; provides `matrices.b` and
  `matrices.K_m₂xm₂`.
- `dof_map_m₁::DOFMap`: DOF map for the ``m_1``-space.
- `dof_map_m₂::DOFMap`: DOF map for the ``m_2``-space.
- `dof_map_m₃::DOFMap`: DOF map for the ``m_3``-space.
- `mesh1D::CartesianMesh{1}`: 1D Cartesian mesh.
- `mesh2D::CartesianMesh{2}`: 2D Cartesian mesh.
- `quad::QuadratureSetup`: precomputed quadrature data.
- `τ::T`: time-step size.
- `τα::T`: product ``\\tau\\alpha(t_{n-1/2})``.
- `input_data::PDEInputData`: problem parameters.
 
# Keyword Arguments
- `abstol::T`: absolute tolerance on ``\\max_i|H_i|`` (default: `T(1e-10)`).
- `maxiter::Int`: maximum number of Newton iterations (default: `10`).
"""
function newton_solve!(
        cache::CrankNicolsonCache{T},
        state::FEMState{T},
        matrices::SystemMatrices,
        dof_map_m₁::DOFMap,
        dof_map_m₂::DOFMap,
        dof_map_m₃::DOFMap,
        mesh1D::CartesianMesh{1},
        mesh2D::CartesianMesh{2},
        quad::QuadratureSetup,
        τ::T,
        τα::T,
        input_data::PDEInputData;
        abstol::T = T(1e-10),
        maxiter::Int = 10
) where {T}
    # Warm start: Xⁿ ← [vⁿ⁻¹; cⁿ⁻¹]
    m₁ = dof_map_m₁.m
    m₂ = dof_map_m₂.m
    cache.Xⁿ[1:m₁] .= state.v
    cache.Xⁿ[(m₁ + 1):(m₁ + m₂)] .= state.c

    for _ in 1:maxiter
        # Compute -H(v̂ⁿ)  →  cache.minusH
        compute_minusH!(cache, state, matrices, dof_map_m₁, dof_map_m₂, dof_map_m₃,
            mesh1D, mesh2D, quad, τ, τα, input_data)

        maximum(abs, cache.minusH) ≤ abstol && return nothing

        # Assemble JH
        compute_JH!(cache, matrices, dof_map_m₁, dof_map_m₂, dof_map_m₃,
            mesh1D, mesh2D, quad, τ, τα, input_data)

        # Newton step: Xⁿ ← Xⁿ + JH⁻¹ (-H)
        cache.Xⁿ .+= cache.JH \ cache.minusH   # FIXME: allocates
    end

    @warn "newton_solve! did not converge within $maxiter iterations " *
          "(max|H| = $(maximum(abs, cache.minusH)), abstol = $abstol)"
    return nothing
end

# ==============================================================================
# compute_minusH!
# ==============================================================================

"""
    compute_minusH!(cache, state, matrices, dof_map_m₁, dof_map_m₂, dof_map_m₃,
                    mesh1D, mesh2D, quad, τ, τα, input_data)
 
Compute ``-H(X^n)`` and store the result in `cache.minusH`:
```math
-H = - Q X^n
   - \\begin{bmatrix}
       \\tau\\alpha G^{m_1}(\\hat{v}^n)
       + \\tau F^{m_1}\\bigl(\\tfrac{\\tau}{2}\\hat{v}^n + d^{n-1}\\bigr)
       \\\\[5pt]
       \\tau\\beta(\\mathbf{b}\\cdot\\hat{c}^n)K^{m_2\\times m_2}\\hat{c}^n
     \\end{bmatrix}
   + \\begin{bmatrix} L_1 \\\\ L_2 \\end{bmatrix}.
```

As a side-effect, `cache.d̂ⁿ` is populated with ``\\hat{d}^n = \\tfrac{\\tau}{2}\\hat{v}^n + d^{n-1}`` and
`cache.K_m₂xm₂_vs_ĉⁿ` with ``K^{m_2\\times m_2}\\hat{c}^n``;
both are reused by the immediately following [`compute_JH!`](@ref) call within each Newton iteration.

Assumes `cache.Q`, `cache.L₁`, `cache.L₂`, and `cache.Xⁿ` have already been populated.
 
# Arguments
- `cache::CrankNicolsonCache`: pre-allocated workspace; `cache.minusH` is updated in-place.
- `state::FEMState`: solution state at time level ``n-1``; provides ``d^{n-1}``.
- `matrices::SystemMatrices`: global FEM matrices; provides `matrices.b` and
  `matrices.K_m₂xm₂`.
- `dof_map_m₁::DOFMap`: DOF map for the ``m_1``-space.
- `dof_map_m₂::DOFMap`: DOF map for the ``m_2``-space.
- `dof_map_m₃::DOFMap`: DOF map for the ``m_3``-space.
- `mesh1D::CartesianMesh{1}`: 1D Cartesian mesh.
- `mesh2D::CartesianMesh{2}`: 2D Cartesian mesh.
- `quad::QuadratureSetup`: precomputed quadrature data.
- `τ::T`: time-step size.
- `τα::T`: product ``\\tau\\alpha(t_{n-1/2})``.
- `input_data::PDEInputData`: problem parameters.
"""
function compute_minusH!(
        cache::CrankNicolsonCache{T},
        state::FEMState{T},
        matrices::SystemMatrices,
        dof_map_m₁::DOFMap,
        dof_map_m₂::DOFMap,
        dof_map_m₃::DOFMap,
        mesh1D::CartesianMesh{1},
        mesh2D::CartesianMesh{2},
        quad::QuadratureSetup,
        τ::T,
        τα::T,
        input_data::PDEInputData
) where {T}
    m₁ = dof_map_m₁.m
    m₂ = dof_map_m₂.m
    m₃ = dof_map_m₃.m
    τ_half = τ / 2

    # Q Xⁿ → cache.minusH  (temporary; sign is flipped in the final assembly)
    mul!(cache.minusH, cache.Q, cache.Xⁿ)

    # τα G(v̂ⁿ_{1:m₃}) → vec_m₃_1
    assembly_nonlinearity_G!(cache.vec_m₃_1, τα, input_data.g,
        @view(cache.v̂ⁿ[1:m₃]), mesh1D, dof_map_m₃, quad)

    # d̂ⁿ = (τ/2) v̂ⁿ + dⁿ⁻¹  →  cache.d̂ⁿ  (reused by compute_JH!)
    @. cache.d̂ⁿ = τ_half * cache.v̂ⁿ + state.d

    # τ F(d̂ⁿ) → vec_m₁_1
    assembly_nonlinearity_F!(
        cache.vec_m₁_1, τ, input_data.f, cache.d̂ⁿ, mesh2D, dof_map_m₁, quad)

    # K_m₂xm₂ ĉⁿ → cache.K_m₂xm₂_vs_ĉⁿ  (reused by compute_JH!)
    mul!(cache.K_m₂xm₂_vs_ĉⁿ, matrices.K_m₂xm₂, cache.ĉⁿ)

    # τ β(b · ĉⁿ) → τβ
    b_dot_ĉⁿ = dot(matrices.b, cache.ĉⁿ)
    β_val = input_data.β(b_dot_ĉⁿ)
    τβ = τ * β_val

    # Final assembly: -H = L - Q Xⁿ - nonlinear contributions
    @inbounds for i in 1:m₃
        cache.minusH[i] = cache.L₁[i] - cache.minusH[i] -
                          cache.vec_m₃_1[i] - cache.vec_m₁_1[i]
    end
    @inbounds for i in (m₃ + 1):m₁
        cache.minusH[i] = cache.L₁[i] - cache.minusH[i] - cache.vec_m₁_1[i]
    end
    @inbounds for i in 1:m₂
        cache.minusH[m₁ + i] = cache.L₂[i] - cache.minusH[m₁ + i] -
                               τβ * cache.K_m₂xm₂_vs_ĉⁿ[i]
    end

    return nothing
end

# ==============================================================================
# compute_JH!
# ==============================================================================

"""
    compute_JH!(cache, state, matrices, dof_map_m₁, dof_map_m₂, dof_map_m₃,
                mesh1D, mesh2D, quad, τ, τα, input_data)
 
Assemble the Jacobian of ``H`` and store the result in `cache.JH`:
```math
JH = Q +
\\begin{bmatrix}
  \\Bigl[
    \\tau\\alpha\\, JG\\bigl(\\hat{v}^n_{1:m_3}\\bigr)
  + \\tfrac{\\tau^2}{2}\\, JF\\bigl(\\hat{d}^n\\bigr)
  \\Bigr]^{m_1 \\times m_1}
  & 0^{m_1 \\times m_2}
  \\\\[10pt]
  0^{m_2 \\times m_1}
  &
  \\tau\\,\\beta'(\\mathbf{b}\\cdot\\hat{c}^n)\\,
  \\bigl[K^{m_2\\times m_2}\\hat{c}^n\\bigr]\\,\\mathbf{b}^T
  + \\tau\\,\\beta(\\mathbf{b}\\cdot\\hat{c}^n)\\,K^{m_2 \\times m_2}
\\end{bmatrix},
```
where ``\\hat{d}^n = \\tfrac{\\tau}{2}\\hat{v}^n + d^{n-1}``.
 
Relies on the following cache fields computed by the immediately preceding [`compute_minusH!`](@ref) call:
- `cache.d̂ⁿ`: midpoint displacement ``\\hat{d}^n``.
- `cache.K_m₂xm₂_vs_ĉⁿ`: product ``K^{m_2\\times m_2}\\hat{c}^n``.
 
The ``\\beta``-nonlinearity Jacobian is assembled entry-by-entry directly into
the bottom-right ``m_2\\times m_2`` block of `cache.JH`, avoiding materialization of another full dense matrix. 
 
Assumes `cache.Q` and `cache.Xⁿ` have already been populated.
 
# Arguments
- `cache::CrankNicolsonCache`: pre-allocated workspace; `cache.JH` is updated in-place.
- `matrices::SystemMatrices`: global FEM matrices; provides `matrices.b` and  `matrices.K_m₂xm₂`.
- `dof_map_m₁::DOFMap`: DOF map for the ``m_1``-space.
- `dof_map_m₂::DOFMap`: DOF map for the ``m_2``-space.
- `dof_map_m₃::DOFMap`: DOF map for the ``m_3``-space.
- `mesh1D::CartesianMesh{1}`: 1D Cartesian mesh.
- `mesh2D::CartesianMesh{2}`: 2D Cartesian mesh.
- `quad::QuadratureSetup`: precomputed quadrature data.
- `τ::T`: time-step size.
- `τα::T`: product ``\\tau\\alpha(t_{n-1/2})``.
- `input_data::PDEInputData`: problem parameters; `input_data.∂ₛg` is the
  derivative of the boundary nonlinearity, `input_data.df` is the derivative
  of the interior nonlinearity, and `input_data.dβ` is the derivative of  ``\\beta``.
"""
function compute_JH!(
        cache::CrankNicolsonCache{T},
        matrices::SystemMatrices,
        dof_map_m₁::DOFMap,
        dof_map_m₂::DOFMap,
        dof_map_m₃::DOFMap,
        mesh1D::CartesianMesh{1},
        mesh2D::CartesianMesh{2},
        quad::QuadratureSetup,
        τ::T,
        τα::T,
        input_data::PDEInputData
) where {T}
    m₁ = dof_map_m₁.m
    m₂ = dof_map_m₂.m
    m₃ = dof_map_m₃.m
    τ²_half = τ^2 / 2

    # JG (m₃×m₃), scaled by τα                                  # FIXME: allocates
    JG = assembly_global_matrix_DG(
        τα, input_data.∂ₛg, @view(cache.v̂ⁿ[1:m₃]), mesh1D, dof_map_m₃, quad)

    # JF (m₁×m₁), evaluated at d̂ⁿ (from cache), scaled by τ²/2  # FIXME: allocates
    JF = assembly_global_matrix_DF(
        τ²_half, input_data.df, cache.d̂ⁿ, mesh2D, dof_map_m₁, quad)

    # JH ← Q
    cache.JH .= cache.Q

    # Embed τα JG + τ²/2 JF into the top-left m₁×m₁ block
    @. cache.JH[1:m₃, 1:m₃] += JG                               # FIXME: allocates
    @. cache.JH[1:m₁, 1:m₁] += JF                               # FIXME: allocates

    # Embed Jβ directly into the bottom-right m₂×m₂ block of JH
    # The rank-1 term (K ĉⁿ)(τ dβ b)ᵀ is evaluated entry-by-entry, avoiding materialization of the full dense matrix
    b_dot_ĉⁿ = dot(matrices.b, cache.ĉⁿ)
    τβ = τ * input_data.β(b_dot_ĉⁿ)
    τdβ = τ * input_data.dβ(b_dot_ĉⁿ)
    @. cache.vec_m₂_1 = τdβ * matrices.b
    @inbounds for j in 1:m₂, i in 1:m₂
        cache.JH[m₁ + i, m₁ + j] = cache.K_m₂xm₂_vs_ĉⁿ[i] * cache.vec_m₂_1[j] +
                                   τβ * matrices.K_m₂xm₂[i, j]
    end

    return nothing
end

# ==============================================================================
# solve_r̂ⁿ! — CrankNicolson
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
 
Assumes `cache.v̂ⁿ` has already been populated by [`newton_solve!`](@ref).
 
# Arguments
- `cache::CrankNicolsonCache`: pre-allocated workspace; `cache.r̂ⁿ` is updated in-place.
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
        cache::CrankNicolsonCache{T},
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
    cst2 = 2q₁ / q₅
    cst3 = τ * q₃ / q₅
    cst4 = τ / q₅

    m₃ = dof_map_m₃.m

    # F(f₃(t_half)) → cache.vec_m₃_1
    scale = mesh1D.Δx[1] / 2
    assembly_rhs_1d!(cache.vec_m₃_1, x -> input_data.f₃(x, t_half),
        scale, quad.W_ϕP, mesh1D, dof_map_m₃, quad.xP)

    # M_m₃⁻¹ F_m₃ → cache.vec_m₃_2
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
# update_state! — CrankNicolson
# ==============================================================================

"""
    update_state!(state, cache::CrankNicolsonCache, τ)
 
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
by [`newton_solve!`](@ref) and [`solve_r̂ⁿ!`](@ref), respectively.
 
# Arguments
- `state::FEMState`: solution state at time level ``n-1``; updated in-place to level ``n``.
- `cache::CrankNicolsonCache`: pre-allocated workspace.
- `τ::T`: time-step size.
"""
function update_state!(
        state::FEMState{T, V},
        cache::CrankNicolsonCache{T},
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