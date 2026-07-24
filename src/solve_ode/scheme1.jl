function solve_ode(
        cache::Scheme1Cache,
        tspan::StepRangeLen{T},
        state::State{T},
        matrices::SystemMatrices{T},
        nel_per_dim::NTuple{2, I},
        element_side_lengths::NTuple{2, T},
        dof_map_m₁::DOFMap,
        dof_map_m₂::DOFMap,
        dof_map_m₃::DOFMap,
        quad::QuadratureSetup,
        input_data::PDEInputData,
        callback::AbstractCallback
) where {T, I}
    τ = step(tspan)
    τ_2 = τ / 2
    τ²_2 = τ * τ / 2

    q₁, q₂, q₃, q₄ = input_data.q₁, input_data.q₂, input_data.q₃, input_data.q₄
    q₅ = 2 * q₁ + τ * q₂ + τ²_2 * q₃

    _2q₁ = 2 * q₁
    τq₃ = τ * q₃
    τ²q₄_q₅ = τ*τ*q₄/q₅

    _τq₄_q₅ = -τ * q₄ / q₅
    _2q₁_q₅ = 2 * q₁ / q₅
    _τq₃_q₅ = -τ * q₃ / q₅
    _q₅ = 1 / q₅

    Δx, Δy = element_side_lengths
    τ_jac_2d = τ * Δx * Δy / 4
    τ_jac_1d = τ * Δx / 2

    csts = (τ, τ_2, τ²_2, _2q₁, τq₃, q₅, τ²q₄_q₅,
        _τq₄_q₅, _2q₁_q₅, _τq₃_q₅, _q₅, τ_jac_1d, τ_jac_2d)

    for n in 1:(length(tspan) - 1)
        perform_step!(
            cache, state, matrices, nel_per_dim, element_side_lengths,
            dof_map_m₁, dof_map_m₂, dof_map_m₃, quad, input_data, n, csts)
        apply!(
            callback, state, nel_per_dim, element_side_lengths,
            dof_map_m₁, dof_map_m₂, dof_map_m₃, quad, input_data)
    end
    return nothing
end

# ==============================================================================
# perform_step!
# ==============================================================================
function perform_step!(
        cache::Scheme1Cache{T},
        state::State{T},
        matrices::SystemMatrices{T},
        nel_per_dim::NTuple{2, I},
        element_side_lengths::NTuple{2, T},
        dof_map_m₁::DOFMap,
        dof_map_m₂::DOFMap,
        dof_map_m₃::DOFMap,
        quad::QuadratureSetup,
        input_data::PDEInputData,
        n::Int,
        csts::NTuple{13, T}
) where {T, I}
    τ, τ_2, τ²_2, _2q₁, τq₃, q₅, τ²q₄_q₅,
    _τq₄_q₅, _2q₁_q₅, _τq₃_q₅, _q₅, τ_jac_1d, τ_jac_2d = csts
    Δx = element_side_lengths[1]

    t_half = (n - T(0.5)) * τ
    α = input_data.α(t_half)

    τα = τ * α
    τα_q₅ = τα / q₅
    τ²α_2 = τ²_2 * α
    τ²q₄α_q₅ = τ²q₄_q₅ * α

    compute_Q₁₁!(
        cache.Q₁₁, matrices.K_m₁xm₁, matrices._2M_m₁xm₁, matrices.M_m₃xm₃, τ²α_2, τ²q₄α_q₅)

    assembly_rhs_2d!(
        cache.τFf₁, (x, y) -> input_data.f₁(x, y, t_half),
        τ_jac_2d, quad.W_φP, dof_map_m₁, nel_per_dim, element_side_lengths, quad.xP, quad.yP)
    assembly_rhs_2d!(
        cache.τFf₂, (x, y) -> input_data.f₂(x, y, t_half),
        τ_jac_2d, quad.W_φP, dof_map_m₂, nel_per_dim, element_side_lengths, quad.xP, quad.yP)
    assembly_rhs_1d!(
        cache.τFf₃, x -> input_data.f₃(x, t_half),
        τ_jac_1d, quad.W_ϕP, dof_map_m₃, Δx, quad.xP)

    compute_L₁!(
        cache.L₁, state, matrices.K_m₁xm₁, matrices._2M_m₁xm₁,
        matrices.M_m₃xm₃, cache.τFf₁, cache.τFf₃,
        τα, τα_q₅, _2q₁, τq₃,
        cache.vec_m₁, cache.vec_m₃_1, cache.vec_m₃_2)
    compute_L₂!(
        cache.L₂, matrices._2M_m₂xm₂, state.c, cache.τFf₂)

    solve_nonlinear_system!(
        cache, state, matrices, element_side_lengths,
        dof_map_m₁, dof_map_m₃, quad, input_data,
        τ, τ_2, τ²_2, τα)

    LS.solve!(cache.linsolve_m₃)
    compute_r̂ⁿ!(cache.r̂ⁿ, cache.v̂ⁿ, state.r, state.z, cache.linsolve_m₃.u,
        _τq₄_q₅, _2q₁_q₅, _τq₃_q₅, _q₅)

    update_state!(state, cache.v̂ⁿ, cache.ĉⁿ, cache.r̂ⁿ, τ)

    return nothing
end

# ==============================================================================
# compute_Q₁₁!
# ==============================================================================
"""
    compute_Q₁₁!(Q₁₁, K_m₁xm₁, _2M_m₁xm₁, M_m₃xm₃, cst₁, cst₂)

Compute in-place

    Q₁₁ = cst₁*K_m₁xm₁ + _2M_m₁xm₁ + cst₂*[M_m₃xm₃ 0_m₃x(m₁-m₃); 0_(m₁-m₃)xm₃ 0_m₃xm₃]

where `cst₁ = τ²⋅α(tₙ₋₁/₂)/2` and `cst₂ = τ²⋅q₄⋅α(tₙ₋₁/₂)/q₅`
"""
function compute_Q₁₁!(
        Q₁₁::Symmetric{T, SparseMatrixCSC{T, I}},
        K_m₁xm₁::Symmetric{T, SparseMatrixCSC{T, I}},
        _2M_m₁xm₁::Symmetric{T, SparseMatrixCSC{T, I}},
        M_m₃xm₃::Symmetric{T, SparseMatrixCSC{T, I}},
        cst₁::T, cst₂::T) where {T, I}
    vec1 = Q₁₁.data.nzval
    vec2 = K_m₁xm₁.data.nzval
    vec3 = _2M_m₁xm₁.data.nzval
    vec4 = M_m₃xm₃.data.nzval

    nnz_m₁ = length(vec1)
    nnz_m₃ = length(vec4)

    for i in 1:nnz_m₃
        tmp = muladd(cst₁, vec2[i], vec3[i])
        vec1[i] = muladd(cst₂, vec4[i], tmp)
    end
    for i in (nnz_m₃ + 1):nnz_m₁
        vec1[i] = muladd(cst₁, vec2[i], vec3[i])
    end

    return nothing
end

# ==============================================================================
# Compute L₁, L₂
# ==============================================================================
"""
    compute_L₁!(
        L₁, state, K_m₁xm₁, _2M_m₁xm₁, M_m₃xm₃, τFf₁, τFf₃, 
        τα, τα_q₅, _2q₁, τq₃, 
        vec_m₁, vec_m₃_1, vec_m₃_2)

Compute in-place

    L₁ = 2M_m₁xm₁⋅vⁿ⁻¹ - τα⋅K_m₁xm₁⋅dⁿ⁻¹ + τFf₁ + τα_q₅⋅[M_m₃xm₃⋅(2q₁⋅rⁿ⁻¹-τq₃⋅zⁿ⁻¹) + τFf₃; 0_(m₁-m₃)]

where `τα = τ⋅α(tₙ₋₁/₂)`, `τα_q₅ = τ⋅α(tₙ₋₁/₂)/q₅`, `τFf₁ = τ⋅F(f₁(tₙ₋₁/₂))`, and `τFf₃ = τ⋅F(f₃(tₙ₋₁/₂))`.
"""
function compute_L₁!(
        L₁::Vector{T},
        state::State{T},
        K_m₁xm₁::Symmetric{T, SparseMatrixCSC{T, I}},
        _2M_m₁xm₁::Symmetric{T, SparseMatrixCSC{T, I}},
        M_m₃xm₃::Symmetric{T, SparseMatrixCSC{T, I}},
        τFf₁::Vector{T},
        τFf₃::Vector{T},
        τα::T, τα_q₅::T, _2q₁::T, τq₃::T,
        vec_m₁::Vector{T}, vec_m₃_1::Vector{T}, vec_m₃_2::Vector{T}) where {T, I}
    vⁿ⁻¹ = state.v
    dⁿ⁻¹ = state.d
    rⁿ⁻¹ = state.r
    zⁿ⁻¹ = state.z

    m₁ = length(vⁿ⁻¹)
    m₃ = length(rⁿ⁻¹)

    # L₁ ← 2M_m₁xm₁ · vⁿ⁻¹
    mul!(L₁, _2M_m₁xm₁, vⁿ⁻¹)

    # vec_m₁ ← K_m₁xm₁⋅dⁿ⁻¹
    mul!(vec_m₁, K_m₁xm₁, dⁿ⁻¹)

    # vec_m₃_2 ← 2q₁⋅rⁿ⁻¹-τq₃⋅zⁿ⁻¹
    @. vec_m₃_2 = _2q₁ * rⁿ⁻¹ - τq₃ * zⁿ⁻¹

    # vec_m₃_1 ← M_m₃×m₃ · vec_m₃_2
    mul!(vec_m₃_1, M_m₃xm₃, vec_m₃_2)

    # Final assembly
    cst = -τα
    for i in 1:m₃
        tmp1 = muladd(cst, vec_m₁[i], τFf₁[i])
        tmp2 = vec_m₃_1[i] + τFf₃[i]
        L₁[i] += muladd(τα_q₅, tmp2, tmp1)
    end
    for i in (m₃ + 1):m₁
        L₁[i] += muladd(cst, vec_m₁[i], τFf₁[i])
    end

    return nothing
end

"""
    compute_L₂!(L₂, _2M_m₂xm₂, cⁿ⁻¹, τFf₂)

Compute in-place

    L₂ = 2M_m₂xm₂⋅cⁿ⁻¹ + τFf₂

where `τFf₂ = τ⋅F(f₂(tₙ₋₁/₂))`.
"""
function compute_L₂!(
        L₂::Vector{T},
        _2M_m₂xm₂::Symmetric{T, SparseMatrixCSC{T, I}},
        cⁿ⁻¹::Vector{T},
        τFf₂::Vector{T}) where {T, I}
    mul!(L₂, _2M_m₂xm₂, cⁿ⁻¹)
    @. L₂ += τFf₂
    return nothing
end

# ==============================================================================
# solve_nonlinear_system!
# ==============================================================================
"""
    solve_nonlinear_system!(
        cache, state, matrices, element_side_lengths,
        dof_map_m₁, dof_map_m₃, quad, input_data,
        τ, τ_2, τ²_2, τα;
        abstol, reltol, maxiter)

Solve the nonlinear system
```math
H(X) = Q X + [τα G^{m₁}(v̂ⁿ) + τ F^{m₁}(d̂ⁿ); τβ(\\mathbf{b}⋅ĉⁿ)K^{m_2\\times m_2}ĉⁿ] - [L₁;L₂] = 0
```
for ``X = [v̂ⁿ;ĉⁿ]`` via Newton's method, updating `cache.X` in-place. 
`cache.X` is warm-started from `[vⁿ⁻¹;cⁿ⁻¹]` and ``d̂ⁿ = (τ/2)v̂ⁿ + dⁿ⁻¹``.

Convergence is declared when either of the following holds at iteration ``k``:
1. Residual criterion:
   ``max|H(Xᵏ)| ≤ abstol``
2. Step criterion:
   ``max|Xᵏ⁺¹-Xᵏ| ≤ abstol + reltol ⋅ max|Xᵏ|``
"""
function solve_nonlinear_system!(
        cache::Scheme1Cache{T},
        state::State{T},
        matrices::SystemMatrices{T},
        element_side_lengths::NTuple{2, T},
        dof_map_m₁::DOFMap,
        dof_map_m₃::DOFMap,
        quad::QuadratureSetup,
        input_data::PDEInputData,
        τ::T, τ_2::T, τ²_2::T, τα::T;
        abstol::T = T(1e-15),
        reltol::T = T(1e-12),
        maxiter::Int = 6) where {T}
    Δx, Δy = element_side_lengths

    # Warm start: X ← [vⁿ⁻¹; cⁿ⁻¹]
    cache.v̂ⁿ .= state.v
    cache.ĉⁿ .= state.c
    @. cache.d̂ⁿ = muladd(τ_2, cache.v̂ⁿ, state.d)

    normH = zero(T)
    normΔX = zero(T)
    normX = zero(T)

    for _ in 1:maxiter
        bᵀĉⁿ = dot(matrices.b, cache.ĉⁿ)                      # b⋅ĉⁿ
        τβ_bᵀĉⁿ = τ * input_data.β(bᵀĉⁿ)                      # τβ(bᵀĉⁿ)
        τdβ_bᵀĉⁿ = τ * input_data.dβ(bᵀĉⁿ)                    # τβ'(bᵀĉⁿ)
        mul!(cache.K_m₂xm₂_ĉⁿ, matrices.K_m₂xm₂, cache.ĉⁿ)    # K_m₂xm₂⋅ĉⁿ

        assembly_nonlinearity_G!(
            cache.G, τα, input_data.g, cache.v̂ⁿ,
            dof_map_m₃, Δx, quad.xP, quad.ϕP, quad.W_ϕP)      # τα⋅G(v̂ⁿ)
        assembly_nonlinearity_F!(
            cache.F, τ, input_data.f, cache.d̂ⁿ, dof_map_m₁,
            element_side_lengths, quad.φP, quad.W_φP)         # τ⋅F(d̂ⁿ)

        compute_minusH!(
            cache.minusH, cache.v̂ⁿ, cache.ĉⁿ, τβ_bᵀĉⁿ, cache.K_m₂xm₂_ĉⁿ,
            cache.Q₁₁, matrices.τA_m₁xm₂, matrices.τA_m₂xm₁, matrices._2M_m₂xm₂,
            cache.G, cache.F, cache.L₁, cache.L₂,
            cache.vec_m₁, cache.vec_m₂)

        normH = maximum(abs, cache.minusH)
        normH ≤ abstol && return nothing                      # criterion 1

        JG = assembly_global_matrix_DG(
            τα, input_data.∂ₛg, cache.v̂ⁿ, dof_map_m₃,
            Δx, quad.xP, quad.ϕP, quad.W_ϕPϕP)                # τα⋅JG(v̂ⁿ)
        JF = assembly_global_matrix_DF(
            τ²_2, input_data.df, cache.d̂ⁿ, dof_map_m₁,
            element_side_lengths, quad.φP, quad.W_φPφP)       # (τ²/2)⋅JF(d̂ⁿ)

        compute_JH₁₁!(                                        # Q₁₁ + τα⋅JG(v̂ⁿ) + (τ²/2)⋅JF(d̂ⁿ)
            cache.JH₁₁, cache.Q₁₁, JG, JF)
        _muladd!(                                             # JH₂₂_sparse ← τβ(b⋅ĉⁿ)⋅K_m₂xm₂ + 2M_m₂xm₂
            cache.JH₂₂_sparse, τβ_bᵀĉⁿ, matrices.K_m₂xm₂, matrices._2M_m₂xm₂)

        scatter_symmetric!(cache.JH_sparse, cache.JH₁₁, cache.map11_direct, cache.map11_mirror)
        scatter_symmetric!(cache.JH_sparse, cache.JH₂₂_sparse, cache.map22_direct, cache.map22_mirror)

        solve_linear_system!(
            cache.linsolve, cache.minusH, τdβ_bᵀĉⁿ, cache.K_m₂xm₂_ĉⁿ,
            matrices.b, cache.vec_m₂, cache.vec_m₁_m₂)

        normΔX = maximum(abs, cache.linsolve.u)
        normX = maximum(abs, cache.X)
        cache.X .+= cache.linsolve.u
        @. cache.d̂ⁿ = muladd(τ_2, cache.v̂ⁿ, state.d)

        normΔX ≤ abstol + reltol * normX && return nothing     # criterion 2
    end

    @warn "solve_nonlinear_system! did not converge in $maxiter iterations " *
          "(‖H(Xᵏ)‖ = $(@sprintf("%.1e", normH)), " *
          "‖Xᵏ⁺¹-Xᵏ‖ = $(@sprintf("%.1e", normΔX)), " *
          "abstol+reltol*‖Xᵏ‖ = $(@sprintf("%.1e", abstol + reltol * normX)), " *
          "‖Xᵏ‖ = $(@sprintf("%.1e", normX)), " *
          "abstol = $(@sprintf("%.1e", abstol)), reltol = $(@sprintf("%.1e", reltol)))"
    return nothing
end

function solve_linear_system!(
        linsolve::LS.LinearCache,
        minusH::Vector{T},
        τdβ_bᵀĉⁿ::T,
        K_m₂xm₂_ĉⁿ::Vector{T},
        b::Vector{T},
        vec_m₂::Vector{T}, vec_m₁_m₂::Vector{T}) where {T}
    m₂ = length(K_m₂xm₂_ĉⁿ)
    m₁ = length(minusH) - m₂

    # Step 1: fresh factorization of JH_sparse, solve w₁ = JH_sparse \ (−H)
    # w₁ is copied into cache.vec_m₁_m₂ before step 2 overwrites linsolve.u
    linsolve.b .= minusH
    linsolve.isfresh = true
    LS.solve!(linsolve)
    vec_m₁_m₂ .= linsolve.u
    w₁ = vec_m₁_m₂

    # Step 2: reuse factorization, solve w₂ = JH_sparse \ u, with u = [0_m₁; b]
    for i in 1:m₁
        linsolve.b[i] = zero(T)
    end
    for i in 1:m₂
        linsolve.b[i + m₁] = b[i]
    end
    linsolve.isfresh = false
    LS.solve!(linsolve)
    w₂ = linsolve.u

    # Step 3: compute σ = (vᵀ w₁) / (1 + vᵀ w₂), with v = [0_m₁; τβ'(b·ĉⁿ)⋅K_m₂xm₂⋅ĉⁿ]
    v = vec_m₂
    @. v = τdβ_bᵀĉⁿ * K_m₂xm₂_ĉⁿ
    σ_num = dot(v, view(w₁, (m₁ + 1):(m₁ + m₂)))
    σ_den = one(T) + dot(v, view(w₂, (m₁ + 1):(m₁ + m₂)))
    abs(σ_den) < eps(T) && error("Sherman-Morrison breakdown: denominator ≈ 0")
    σ = σ_num / σ_den

    # Step 4: update solution
    @. linsolve.u = w₁ - σ * w₂

    return nothing
end

# ==============================================================================
# compute_minusH!
# ==============================================================================
"""
    compute_minusH!(
        minusH, v̂ⁿ, ĉⁿ, τβ_bᵀĉⁿ, K_m₂xm₂_ĉⁿ,
        Q₁₁, τA_m₁xm₂, τA_m₂xm₁, _2M_m₂xm₂,
        G, F, L₁, L₂, 
        vec_m₁, vec_m₂)

Compute in-place ``-H(X)`` where
```math
-H_{block1} = - Q₁₁⋅v̂ⁿ - τA_m₁xm₂⋅ĉⁿ - τα⋅G(v̂ⁿ) - τF(d̂ⁿ) + L₁\\
-H_{block2} = - τA_m₂xm₁⋅v̂ⁿ - 2M_m₂xm₂⋅ĉⁿ - τβ(b⋅ĉⁿ)⋅K_m₂xm₂⋅ĉⁿ + L₂
```
"""
function compute_minusH!(
        minusH, v̂ⁿ, ĉⁿ, τβ_bᵀĉⁿ, K_m₂xm₂_ĉⁿ,
        Q₁₁, τA_m₁xm₂, τA_m₂xm₁, _2M_m₂xm₂,
        G, F, L₁, L₂,
        vec_m₁, vec_m₂)
    m₁ = length(L₁)
    m₂ = length(L₂)
    m₃ = length(G)

    # Block 1
    mul!(view(minusH, 1:m₁), Q₁₁, v̂ⁿ)
    mul!(vec_m₁, τA_m₁xm₂, ĉⁿ)
    for i in 1:m₃
        minusH[i] = -minusH[i] - vec_m₁[i] - G[i] - F[i] + L₁[i]
    end
    for i in (m₃ + 1):m₁
        minusH[i] = -minusH[i] - vec_m₁[i] - F[i] + L₁[i]
    end

    # Block 2
    mul!(view(minusH, (m₁ + 1):(m₁ + m₂)), τA_m₂xm₁, v̂ⁿ)
    mul!(vec_m₂, _2M_m₂xm₂, ĉⁿ)
    for i in 1:m₂
        minusH[i + m₁] = -minusH[i + m₁] - vec_m₂[i] - τβ_bᵀĉⁿ * K_m₂xm₂_ĉⁿ[i] + L₂[i]
    end

    return nothing
end

# ==============================================================================
# Compute JH
# ==============================================================================
"""
    compute_JH₁₁!(JH₁₁, Q₁₁, JG, JF)

Compute in-place
```math
JH₁₁ = Q₁₁ + τα⋅JG(v̂ⁿ) + (τ²/2)⋅JF(d̂ⁿ).
```

Assumes `JH₁₁`, `Q₁₁`, and `JF` share the same sparsity pattern, and that `JG` occupies
the leading `nnz(JG)` entries of that pattern.
"""
function compute_JH₁₁!(
        JH₁₁::Symmetric{T, SparseMatrixCSC{T, I}},
        Q₁₁::Symmetric{T, SparseMatrixCSC{T, I}},
        JG::Symmetric{T, SparseMatrixCSC{T, I}},
        JF::Symmetric{T, SparseMatrixCSC{T, I}}) where {T, I}
    nnz_m₁ = nnz(JH₁₁.data)
    nnz_m₃ = nnz(JG.data)

    nzval_JH₁₁ = JH₁₁.data.nzval
    nzval_Q₁₁ = Q₁₁.data.nzval
    nzval_JF = JF.data.nzval
    nzval_JG = JG.data.nzval

    for i in 1:nnz_m₃
        nzval_JH₁₁[i] = nzval_Q₁₁[i] + nzval_JG[i] + nzval_JF[i]
    end
    for i in (nnz_m₃ + 1):nnz_m₁
        nzval_JH₁₁[i] = nzval_Q₁₁[i] + nzval_JF[i]
    end

    return nothing
end

# ==============================================================================
# solve_r̂ⁿ! — CrankNicolson
# ==============================================================================
"""
    compute_r̂ⁿ!(r̂ⁿ, v̂ⁿ, rⁿ⁻¹, zⁿ⁻¹, sol, cst1, cst2, cst3, cst4)

Compute in-place 

    r̂ⁿ = cst1⋅v̂ⁿ[1:m₃] + cst2⋅rⁿ⁻¹ + cst3⋅zⁿ⁻¹ + cst4⋅sol

where `cst1 = -τ⋅q₄/q₅`, `cst2 = 2⋅q₁/q₅`, `cst3 = -τ⋅q₃/q₅`, `cst4 = 1/q₅`, and `sol = inv(M_m₃xm₃)⋅τFf₃`
"""
function compute_r̂ⁿ!(
        r̂ⁿ, v̂ⁿ, rⁿ⁻¹, zⁿ⁻¹, sol, cst1, cst2, cst3, cst4)
    for i in eachindex(r̂ⁿ)
        r̂ⁿ[i] = cst1 * v̂ⁿ[i] + cst2 * rⁿ⁻¹[i] + cst3 * zⁿ⁻¹[i] + cst4 * sol[i]
    end
    return nothing
end