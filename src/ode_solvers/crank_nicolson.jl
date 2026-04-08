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

    compute_Q₁₁!(cache, matrices, τ, α, input_data.q₄, q₅)
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
# compute_Q₁₁!
# ==============================================================================

"""
    compute_Q₁₁!(cache, matrices, τ, α, q₄, q₅)
 
Update `cache.Q₁₁` in-place with the top-left ``m_1 \\times m_1`` block of ``Q``:
```math
Q₁₁ =
2M^{m_1 \\times m_1}
+ \\frac{\\tau^2}{2}\\alpha K^{m_1 \\times m_1}
+ \\frac{\\tau^2 q_4}{q_5}\\alpha
  \\begin{bmatrix}
  M^{m_3\\times m_3}       & 0^{m_3\\times(m_1-m_3)}\\\\[5pt]
  0^{(m_1-m_3)\\times m_3} & 0^{(m_1-m_3)\\times(m_1-m_3)}
  \\end{bmatrix}.
```
 
The remaining three blocks of ``Q`` are time-invariant and stored as separate
cache fields (see [`CrankNicolsonCache`](@ref)):
```math
Q =
\\begin{bmatrix}
  Q_{11}                      & \\tau A^{m_1 \\times m_2} \\\\[4pt]
  \\tau A^{m_2 \\times m_1}   & 2M^{m_2 \\times m_2}
\\end{bmatrix}.
```

# Arguments
- `cache::CrankNicolsonCache`: pre-allocated workspace; `cache.Q₁₁` is updated in-place.
- `matrices::SystemMatrices`: global FEM matrices.
- `τ::T`: time-step size.
- `α::T`: value of ``\\alpha(t_{n-1/2})``.
- `q₄::T`: problem parameter ``q_4``.
- `q₅::T`: scalar ``q_5 = 2q_1 + \\tau q_2 + (\\tau^2/2)q_3``.
"""
function compute_Q₁₁!(
        cache::CrankNicolsonCache{T, I, TF, TC},
        matrices::SystemMatrices{T, I},
        τ::T,
        α::T,
        q₄::T,
        q₅::T
) where {T, I, TF, TC}
    cst1 = (τ^2 / 2) * α
    cst2 = (τ^2 * q₄ / q₅) * α

    nzval_Q₁₁ = cache.Q₁₁.data.nzval
    nzval_2M₁₁ = cache.M_m₁xm₁_vs2.data.nzval
    nzval_K₁₁ = matrices.K_m₁xm₁.data.nzval
    nzval_M₃₃ = matrices.M_m₃xm₃.data.nzval

    nnz_m₃ = length(nzval_M₃₃)
    nnz_m₁ = length(nzval_Q₁₁)

    # Sub-block m₃×m₃: Q₁₁ ← 2M^{m₁×m₁} + cst1·K^{m₁×m₁} + cst2·M^{m₃×m₃}
    @inbounds @simd for i in 1:nnz_m₃
        nzval_Q₁₁[i] = nzval_2M₁₁[i] + cst1 * nzval_K₁₁[i] + cst2 * nzval_M₃₃[i]
    end

    # Remaining m₁×m₁ entries: Q₁₁ ← 2M^{m₁×m₁} + cst1·K^{m₁×m₁}
    @inbounds @simd for i in (nnz_m₃ + 1):nnz_m₁
        nzval_Q₁₁[i] = nzval_2M₁₁[i] + cst1 * nzval_K₁₁[i]
    end

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
H(X^n) = Q X^n
+ \\begin{bmatrix}
    \\tau\\alpha G^{m_1}(\\hat{v}^n)
    + \\tau F^{m_1}\\bigl(\\tfrac{\\tau}{2}\\hat{v}^n + d^{n-1}\\bigr)
    \\\\[5pt]
    \\tau\\beta(\\mathbf{b}\\cdot\\hat{c}^n)K^{m_2\\times m_2}\\hat{c}^n
  \\end{bmatrix}
- \\begin{bmatrix} L_1 \\\\ L_2 \\end{bmatrix}
= 0
```
for the unknown ``X^n = [\\hat{v}^n;\\,\\hat{c}^n]`` via Newton's method, updating `cache.Xⁿ` (and its aliases `cache.v̂ⁿ`, `cache.ĉⁿ`) in-place. 
 
Convergence is declared when ``\\max_i |H_i(X^n)| \\leq \\texttt{abstol}``.
 
Assumes `cache.Q₁₁`, `cache.L₁`, and `cache.L₂` have already been populated.
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
    cache.v̂ⁿ .= state.v
    cache.ĉⁿ .= state.c

    for _ in 1:maxiter
        # Compute -H(v̂ⁿ)  →  cache.minusH
        compute_minusH!(cache, state, matrices, dof_map_m₁, dof_map_m₂, dof_map_m₃,
            mesh1D, mesh2D, quad, τ, τα, input_data)

        maximum(abs, cache.minusH) ≤ abstol && return nothing

        # Assemble JH
        compute_JH_sparse!(cache, matrices, dof_map_m₁, dof_map_m₂, dof_map_m₃,
            mesh1D, mesh2D, quad, τ, τα, input_data)

        # Solve linear system: JH ΔXⁿ = minusH
        solve_newton_linear_system!(cache, matrices, τ, input_data)

        # Newton step: Xⁿ ← Xⁿ + JH⁻¹ (-H)
        cache.Xⁿ .+= cache.ΔXⁿ
    end

    @warn "newton_solve! did not converge within $maxiter iterations " *
          "(max|H| = $(maximum(abs, cache.minusH)), abstol = $abstol)"
    return nothing
end

"""
    solve_newton_linear_system!(cache, matrices, τ, input_data)

Solve the Newton linear system ``JH \\cdot \\Delta X^n = -H`` via the
Sherman–Morrison formula, exploiting the rank-1 structure of ``JH_{22}``:
```math
JH = JH_{\\text{sparse}} + u v^T,
```
where
```math
u =
\\begin{bmatrix} 0^{m_1} \\\\ \\mathbf{b} \\end{bmatrix},
\\qquad
v =
\\begin{bmatrix} 0^{m_1} \\\\ v_2 \\end{bmatrix},
\\qquad
v_2 = \\tau\\beta'(\\mathbf{b}\\cdot\\hat{c}^n)\\,K^{m_2\\times m_2}\\hat{c}^n
\\in \\mathbb{R}^{m_2}.
```
Since ``v`` vanishes on the first ``m_1`` components, all inner products reduce
to the ``m_2`` block; ``v`` is never assembled as a full vector.

Algorithm:
1. Assemble ``v_2 = \\tau\\beta'(\\mathbf{b}\\cdot\\hat{c}^n)\\,K^{m_2\\times m_2}\\hat{c}^n``.
2. Set `isfresh = true` and solve ``w_1 = JH_{\\text{sparse}}^{-1}(-H)``, triggering a fresh KLU numeric factorization.
3. Set `isfresh = false` and solve ``w_2 = JH_{\\text{sparse}}^{-1} u``, reusing the factorization from step 2.
4. Compute the Sherman–Morrison scalar
   ``\\sigma = (v_2^T w_{1,2}) / (1 + v_2^T w_{2,2})``,
   where the subscript ``2`` denotes the ``m_2`` block of each vector.
5. Correct: ``\\Delta X^n = w_1 - \\sigma\\, w_2``.

The result is stored in `cache.ΔXⁿ`.
Before calling this function, the following cache fields must be populated:
- `cache.minusH` — by [`compute_minusH!`](@ref);
- `cache.K_m₂xm₂_vs_ĉⁿ` — by [`compute_minusH!`](@ref);
- `cache.JH_sparse` — by [`compute_JH_sparse!`](@ref).

# Arguments
- `cache::CrankNicolsonCache`: pre-allocated workspace; `cache.ΔXⁿ` is updated in-place.
- `matrices::SystemMatrices`: provides `matrices.b`.
- `τ::T`: time-step size.
- `input_data::PDEInputData`: provides `input_data.dβ`.
"""
function solve_newton_linear_system!(
        cache::CrankNicolsonCache{T, I, TF, TC},
        matrices::SystemMatrices{T, I},
        τ::T,
        input_data
) where {T, I, TF, TC}
    m₁ = size(matrices.M_m₁xm₁, 1)
    m₂ = size(matrices.M_m₂xm₂, 1)

    # Step 1: assemble v₂ = τβ'(b·ĉⁿ) K^{m₂×m₂} ĉⁿ
    b_dot_ĉⁿ = dot(matrices.b, cache.ĉⁿ)
    τdβ = τ * input_data.dβ(b_dot_ĉⁿ)
    v₂ = cache.vec_m₂_1
    @inbounds @simd for i in eachindex(v₂)
        v₂[i] = τdβ * cache.K_m₂xm₂_vs_ĉⁿ[i]
    end

    # Step 2: fresh factorization of JH_sparse, solve w₁ = JH_sparse \ (−H)
    # w₁ is copied into ΔXⁿ before step 3 overwrites linsolve.u
    cache.linsolve.b .= cache.minusH
    cache.linsolve.isfresh = true
    LS.solve!(cache.linsolve)
    cache.ΔXⁿ .= cache.linsolve.u
    w₁ = cache.ΔXⁿ

    # Step 3: reuse factorization, solve w₂ = JH_sparse \ u
    @inbounds @simd for i in 1:m₁
        cache.linsolve.b[i] = zero(T)
    end
    @inbounds @simd for i in 1:m₂
        cache.linsolve.b[m₁ + i] = matrices.b[i]
    end
    cache.linsolve.isfresh = false
    LS.solve!(cache.linsolve)
    w₂ = cache.linsolve.u

    # Step 4: σ = (v₂ᵀ w₁₂) / (1 + v₂ᵀ w₂₂),  subscript 2 = m₂ block
    @views σ_num = dot(v₂, w₁[(m₁ + 1):end])
    @views σ_den = one(T) + dot(v₂, w₂[(m₁ + 1):end])
    abs(σ_den) < eps(T) && error("Sherman–Morrison breakdown: denominator ≈ 0")
    σ = σ_num / σ_den

    # Step 5: ΔXⁿ ← w₁ − σ w₂
    @inbounds @simd for i in eachindex(cache.ΔXⁿ)
        cache.ΔXⁿ[i] -= σ * w₂[i]
    end

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

Assumes `cache.Q₁₁`, `cache.L₁`, `cache.L₂`, and `cache.Xⁿ` have already been populated.
 
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

    # Q Xⁿ → cache.QXⁿ
    compute_QXⁿ!(cache)

    # τα G(v̂ⁿ_{1:m₃}) → vec_m₃_1
    assembly_nonlinearity_G!(cache.vec_m₃_1, τα, input_data.g,
        @view(cache.v̂ⁿ[1:m₃]), mesh1D, dof_map_m₃, quad)

    # d̂ⁿ = (τ/2) v̂ⁿ + dⁿ⁻¹  →  cache.d̂ⁿ  (reused by JH)
    @. cache.d̂ⁿ = τ_half * cache.v̂ⁿ + state.d

    # τ F(d̂ⁿ) → vec_m₁_1
    assembly_nonlinearity_F!(
        cache.vec_m₁_1, τ, input_data.f, cache.d̂ⁿ, mesh2D, dof_map_m₁, quad)

    # K_m₂xm₂ ĉⁿ → cache.K_m₂xm₂_vs_ĉⁿ  (reused by JH)
    mul!(cache.K_m₂xm₂_vs_ĉⁿ, matrices.K_m₂xm₂, cache.ĉⁿ)

    # τ β(b · ĉⁿ) → τβ
    b_dot_ĉⁿ = dot(matrices.b, cache.ĉⁿ)
    β_val = input_data.β(b_dot_ĉⁿ)
    τβ = τ * β_val

    # Final assembly: -H = L - Q Xⁿ - nonlinear contributions
    @inbounds for i in 1:m₃
        cache.minusH[i] = cache.L₁[i] - cache.QXⁿ[i] -
                          cache.vec_m₃_1[i] - cache.vec_m₁_1[i]
    end
    @inbounds for i in (m₃ + 1):m₁
        cache.minusH[i] = cache.L₁[i] - cache.QXⁿ[i] - cache.vec_m₁_1[i]
    end
    @inbounds for i in 1:m₂
        cache.minusH[m₁ + i] = cache.L₂[i] - cache.QXⁿ[m₁ + i] -
                               τβ * cache.K_m₂xm₂_vs_ĉⁿ[i]
    end

    return nothing
end

"""
    compute_QXⁿ!(cache)

Compute ``QX^n`` in-place into `cache.QXⁿ`, exploiting the block structure of ``Q``:
```math
QX^n =
\\begin{bmatrix}
  Q_{11}                     & \\tau A^{m_1 \\times m_2} \\\\[4pt]
  \\tau A^{m_2 \\times m_1}  & 2M^{m_2 \\times m_2}
\\end{bmatrix}
\\begin{bmatrix} \\hat{v}^n \\\\ \\hat{c}^n \\end{bmatrix}.
```

# Arguments
- `cache::CrankNicolsonCache`: pre-allocated workspace; `cache.QXⁿ` is updated in-place.
"""
function compute_QXⁿ!(cache::CrankNicolsonCache{T, I, TF, TC}) where {T, I, TF, TC}
    # (QXⁿ)₁ = Q₁₁ v̂ⁿ + τA_m₁xm₂ ĉⁿ
    mul!(cache.vec_m₁_1, cache.Q₁₁, cache.v̂ⁿ)
    mul!(cache.vec_m₁_2, cache.τA_m₁xm₂, cache.ĉⁿ)
    @inbounds @simd for i in eachindex(cache.vec_m₁_1)
        cache.QXⁿ[i] = cache.vec_m₁_1[i] + cache.vec_m₁_2[i]
    end

    # (QXⁿ)₂ = τA_m₂xm₁ v̂ⁿ + 2M_m₂xm₂ ĉⁿ
    m₁ = length(cache.v̂ⁿ)
    mul!(cache.vec_m₂_1, cache.τA_m₂xm₁, cache.v̂ⁿ)
    mul!(cache.vec_m₂_2, cache.M_m₂xm₂_vs2, cache.ĉⁿ)
    @inbounds @simd for i in eachindex(cache.vec_m₂_1)
        cache.QXⁿ[m₁ + i] = cache.vec_m₂_1[i] + cache.vec_m₂_2[i]
    end

    return nothing
end

# ==============================================================================
# compute_JH_sparse!
# ==============================================================================

"""
    compute_JH_sparse!(cache, state, matrices, dof_map_m₁, dof_map_m₂, dof_map_m₃,
                mesh1D, mesh2D, quad, τ, τα, input_data)
 
Assemble the sparse part of the Jacobian of ``H`` into `cache.JH_sparse`:
```math
JH_{\\text{sparse}} =
\\begin{bmatrix}
  Q_{11}
  + \\tau\\alpha\\, JG\\bigl(\\hat{v}^n_{1:m_3}\\bigr)
  + \\tfrac{\\tau^2}{2}\\, JF\\bigl(\\hat{d}^n\\bigr)
  & \\tau A^{m_1 \\times m_2}
  \\\\[6pt]
  \\tau A^{m_2 \\times m_1}
  &
  2M^{m_2 \\times m_2} + \\tau\\,\\beta(\\mathbf{b}\\cdot\\hat{c}^n)\\,K^{m_2 \\times m_2}
\\end{bmatrix},
```
where ``\\hat{d}^n = \\tfrac{\\tau}{2}\\hat{v}^n + d^{n-1}``. 
The off-diagonal blocks ``\\tau A^{m_1 \\times m_2}`` and ``\\tau A^{m_2 \\times m_1}`` are
time-invariant and set once at cache construction. 
The full second diagonal block of ``JH`` is
```math
JH_{22} = 2M^{m_2 \\times m_2}
  + \\tau\\,\\beta(\\mathbf{b}\\cdot\\hat{c}^n)\\,K^{m_2 \\times m_2}
  + \\underbrace{
        \\tau\\,\\beta'(\\mathbf{b}\\cdot\\hat{c}^n)
        \\bigl[K^{m_2\\times m_2}\\hat{c}^n\\bigr]\\mathbf{b}^T
    }_{\\text{rank-1, excluded from }JH_{\\text{sparse}}};
```
the rank-1 term is handled separately via Sherman-Morrison formula in the Newton solve.
 
Requires `cache.d̂ⁿ` and `cache.K_m₂xm₂_vs_ĉⁿ` populated by the preceding
[`compute_minusH!`](@ref) call. Internally calls [`compute_JH₁₁!`](@ref),
[`compute_JH₂₂_sparse!`](@ref), and [`sync_JH_sparse!`](@ref).

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
function compute_JH_sparse!(
        cache::CrankNicolsonCache{T, I, TF, TC},
        matrices::SystemMatrices{T, I},
        dof_map_m₁::DOFMap,
        dof_map_m₂::DOFMap,
        dof_map_m₃::DOFMap,
        mesh1D::CartesianMesh{1, I},
        mesh2D::CartesianMesh{2, I},
        quad::QuadratureSetup,
        τ::T,
        τα::T,
        input_data::PDEInputData
) where {T, I, TF, TC}
    compute_JH₁₁!(cache, dof_map_m₁, dof_map_m₃,
        mesh1D, mesh2D, quad, τ, τα, input_data)
    compute_JH₂₂_sparse!(cache, matrices, τ, input_data)
    sync_JH_sparse!(cache)
    return nothing
end

"""
    compute_JH₁₁!(cache, dof_map_m₁, dof_map_m₃, mesh1D, mesh2D, quad, τ, τα, input_data)

Update `cache.JH₁₁` in-place with the top-left ``m_1 \\times m_1`` block of ``JH``:
```math
JH_{11} = Q_{11}
  + \\tau\\alpha\\, JG\\bigl(\\hat{v}^n_{1:m_3}\\bigr)
  + \\tfrac{\\tau^2}{2}\\, JF\\bigl(\\hat{d}^n\\bigr),
```
where ``JG \\in \\mathbb{R}^{m_3 \\times m_3}`` is embedded in the leading sub-block
and ``JF \\in \\mathbb{R}^{m_1 \\times m_1}`` spans the full block.

`cache.JH₁₁` is stored as `Symmetric{SparseMatrixCSC}`, enabling more efficient
computations than a plain sparse matrix. Its values are then copied into the
corresponding block of `cache.JH_sparse` by [`sync_JH_sparse!`](@ref).

Requires `cache.d̂ⁿ` from the preceding [`compute_minusH!`](@ref) call.
"""
function compute_JH₁₁!(
        cache::CrankNicolsonCache{T, I, TF, TC},
        dof_map_m₁::DOFMap,
        dof_map_m₃::DOFMap,
        mesh1D::CartesianMesh{1, I},
        mesh2D::CartesianMesh{2, I},
        quad::QuadratureSetup,
        τ::T,
        τα::T,
        input_data::PDEInputData
) where {T, I, TF, TC}
    m₃ = dof_map_m₃.m

    # FIXME: allocates
    JG = assembly_global_matrix_DG(
        τα, input_data.∂ₛg, @view(cache.v̂ⁿ[1:m₃]), mesh1D, dof_map_m₃, quad)
    # FIXME: allocates
    JF = assembly_global_matrix_DF(
        τ^2 / 2, input_data.df, cache.d̂ⁿ, mesh2D, dof_map_m₁, quad)

    nzval_JH₁₁ = cache.JH₁₁.data.nzval
    nzval_Q₁₁ = cache.Q₁₁.data.nzval
    nzval_JF = JF.data.nzval
    nzval_JG = JG.data.nzval
    nnz_m₃ = length(nzval_JG)
    nnz_m₁ = length(nzval_JH₁₁)

    # Leading m₃×m₃ sub-block: Q₁₁ + JF + JG
    @inbounds @simd for i in 1:nnz_m₃
        nzval_JH₁₁[i] = nzval_Q₁₁[i] + nzval_JF[i] + nzval_JG[i]
    end

    # Remaining m₁×m₁ entries: Q₁₁ + JF
    @inbounds @simd for i in (nnz_m₃ + 1):nnz_m₁
        nzval_JH₁₁[i] = nzval_Q₁₁[i] + nzval_JF[i]
    end

    return nothing
end

"""
    compute_JH₂₂_sparse!(cache, matrices, τ, input_data)

Update `cache.JH₂₂_sparse` in-place with the sparse part of the bottom-right
``m_2 \\times m_2`` block of ``JH``:
```math
JH_{22,\\text{sparse}} = 2M^{m_2 \\times m_2}
  + \\tau\\,\\beta(\\mathbf{b}\\cdot\\hat{c}^n)\\,K^{m_2 \\times m_2}.
```
The rank-1 correction is excluded; see [`compute_JH_sparse!`](@ref).

Analogously to `cache.JH₁₁`, `cache.JH₂₂_sparse` is stored as
`Symmetric{SparseMatrixCSC}` and its values are subsequently copied into
`cache.JH_sparse` by [`sync_JH_sparse!`](@ref).
"""
function compute_JH₂₂_sparse!(
        cache::CrankNicolsonCache{T, I, TF, TC},
        matrices::SystemMatrices{T, I},
        τ::T,
        input_data::PDEInputData
) where {T, I, TF, TC}
    b_dot_ĉⁿ = dot(matrices.b, cache.ĉⁿ)
    τβ = τ * input_data.β(b_dot_ĉⁿ)

    nzval_JH₂₂ = cache.JH₂₂_sparse.data.nzval
    nzval_2M₂₂ = cache.M_m₂xm₂_vs2.data.nzval
    nzval_K₂₂ = matrices.K_m₂xm₂.data.nzval

    @inbounds @simd for i in eachindex(nzval_JH₂₂)
        nzval_JH₂₂[i] = nzval_2M₂₂[i] + τβ * nzval_K₂₂[i]
    end

    return nothing
end

"""
    sync_JH_sparse!(cache)

Copy the values of `cache.JH₁₁` and `cache.JH₂₂_sparse` into the corresponding diagonal blocks of `cache.JH_sparse`. 
The off-diagonal blocks are time-invariant and remain unchanged.
"""
function sync_JH_sparse!(cache::CrankNicolsonCache{T, I, TF, TC}) where {T, I, TF, TC}
    m₁ = length(cache.v̂ⁿ)
    m₂ = length(cache.ĉⁿ)

    # Top-left m₁×m₁ block: JH_sparse ← JH₁₁
    JH₁₁ = sparse(cache.JH₁₁)
    cache.JH_sparse[1:m₁, 1:m₁] .= JH₁₁

    # Bottom-right m₂×m₂ block: JH_sparse ← JH₂₂_sparse
    JH₂₂_sparse = sparse(cache.JH₂₂_sparse)
    cache.JH_sparse[(m₁ + 1):(m₁ + m₂), (m₁ + 1):(m₁ + m₂)] .= JH₂₂_sparse

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
    cst2 = 2 * q₁ / q₅
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