function solve_ode(
        cache::Scheme2Cache,
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

    # Step "1,0" (predictor) + step n=1 (corrector)
    perform_first_step!(
        cache, state, matrices, nel_per_dim, element_side_lengths,
        dof_map_m₁, dof_map_m₂, dof_map_m₃, quad, input_data, csts)
    apply!(
        callback, state, nel_per_dim, element_side_lengths,
        dof_map_m₁, dof_map_m₂, dof_map_m₃, quad, input_data)

    # Steps n ≥ 2
    for n in 2:(length(tspan) - 1)
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
        cache::Scheme2Cache{T},
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

    # d*ⁿ = (3*dⁿ⁻¹ - dⁿ⁻²)/2, c*ⁿ = (3*cⁿ⁻¹ - cⁿ⁻²)/2
    @. cache.dᵃⁿ .= T(0.5) * muladd(3, state.d, -cache.dⁿ⁻²)
    @. cache.cᵃⁿ .= T(0.5) * muladd(3, state.c, -cache.cⁿ⁻²)

    solve_v̂ⁿ(
        cache, state, matrices, nel_per_dim, element_side_lengths,
        dof_map_m₁, dof_map_m₃, quad, input_data,
        τ, τα, τα_q₅, _2q₁, τq₃, τ²α_2, τ²q₄α_q₅, Δx, τ_jac_1d, τ_jac_2d, t_half)
    solve_ĉⁿ(
        cache, state, matrices, nel_per_dim, element_side_lengths,
        dof_map_m₂, quad, input_data, τ, τ_jac_2d, t_half)

    LS.solve!(cache.linsolve_m₃)
    compute_r̂ⁿ!(cache.r̂ⁿ, cache.v̂ⁿ, state.r, state.z, cache.linsolve_m₃.u,
        _τq₄_q₅, _2q₁_q₅, _τq₃_q₅, _q₅)

    # dⁿ⁻² ← dⁿ⁻¹, cⁿ⁻² ← cⁿ⁻¹, update state solution
    cache.dⁿ⁻² .= state.d
    cache.cⁿ⁻² .= state.c
    update_state!(state, cache.v̂ⁿ, cache.linsolve_m₂.u, cache.r̂ⁿ, τ)

    return nothing
end

function perform_first_step!(
        cache::Scheme2Cache{T},
        state::State{T},
        matrices::SystemMatrices{T},
        nel_per_dim::NTuple{2, I},
        element_side_lengths::NTuple{2, T},
        dof_map_m₁::DOFMap,
        dof_map_m₂::DOFMap,
        dof_map_m₃::DOFMap,
        quad::QuadratureSetup,
        input_data::PDEInputData,
        csts::NTuple{13, T}
) where {T, I}
    τ, τ_2, τ²_2, _2q₁, τq₃, q₅, τ²q₄_q₅,
    _τq₄_q₅, _2q₁_q₅, _τq₃_q₅, _q₅, τ_jac_1d, τ_jac_2d = csts
    Δx = element_side_lengths[1]

    t_half = T(0.5) * τ
    α = input_data.α(t_half)

    τα = τ * α
    τα_q₅ = τα / q₅
    τ²α_2 = τ²_2 * α
    τ²q₄α_q₅ = τ²q₄_q₅ * α

    # --- Predictor: d*¹ = d⁰,  c*¹ = c⁰  ---
    cache.dᵃⁿ .= state.d
    cache.cᵃⁿ .= state.c
    solve_v̂ⁿ(
        cache, state, matrices, nel_per_dim, element_side_lengths,
        dof_map_m₁, dof_map_m₃, quad, input_data,
        τ, τα, τα_q₅, _2q₁, τq₃, τ²α_2, τ²q₄α_q₅, Δx, τ_jac_1d, τ_jac_2d, t_half)
    solve_ĉⁿ(
        cache, state, matrices, nel_per_dim, element_side_lengths,
        dof_map_m₂, quad, input_data, τ, τ_jac_2d, t_half)

    # --- Corrector: d*¹ = (τ/2)v̂^{"1,0"} + d⁰,  c*¹ = ĉ^{"1,0"} ---
    @. cache.dᵃⁿ = muladd(τ_2, cache.v̂ⁿ, state.d)
    cache.cᵃⁿ .= cache.linsolve_m₂.u
    solve_v̂ⁿ(
        cache, state, matrices, nel_per_dim, element_side_lengths,
        dof_map_m₁, dof_map_m₃, quad, input_data,
        τ, τα, τα_q₅, _2q₁, τq₃, τ²α_2, τ²q₄α_q₅, Δx, τ_jac_1d, τ_jac_2d, t_half)
    solve_ĉⁿ(
        cache, state, matrices, nel_per_dim, element_side_lengths,
        dof_map_m₂, quad, input_data, τ, τ_jac_2d, t_half)

    # --- Compute r̂ⁿ ---
    LS.solve!(cache.linsolve_m₃)
    compute_r̂ⁿ!(cache.r̂ⁿ, cache.v̂ⁿ, state.r, state.z, cache.linsolve_m₃.u,
        _τq₄_q₅, _2q₁_q₅, _τq₃_q₅, _q₅)

    # dⁿ⁻² ← dⁿ⁻¹, cⁿ⁻² ← cⁿ⁻¹, update state solution
    cache.dⁿ⁻² .= state.d
    cache.cⁿ⁻² .= state.c
    update_state!(state, cache.v̂ⁿ, cache.linsolve_m₂.u, cache.r̂ⁿ, τ)

    return nothing
end

"""
    solve_v̂ⁿ(
        cache, state, matrices, nel_per_dim, element_side_lengths,
        dof_map_m₁, dof_map_m₃, quad, input_data, 
        τ, τα, τα_q₅, _2q₁, τq₃, τ²α_2, τ²q₄α_q₅, Δx, τ_jac_1d, τ_jac_2d, t_half)

Sequentially:
1. Computes `Q₁(n)`
2. Computes `τ⋅F(f₁(tₙ₋₁/₂))`
3. Computes `τ⋅F(f₃(tₙ₋₁/₂))`
4. Computes `τ⋅F(dᵃⁿ)`
5. Computes `L₁(n,dᵃⁿ,cᵃⁿ)`
6. Solves the nonlinear system `Q₁⋅v̂ⁿ + τα⋅G(v̂ⁿ) - L₁ = 0`

**WARNING**: Assumes that `cache.dᵃⁿ` and `cache.cᵃⁿ` have already been populated.
"""
function solve_v̂ⁿ(
        cache::Scheme2Cache{T},
        state::State{T},
        matrices::SystemMatrices{T},
        nel_per_dim::NTuple{2, I},
        element_side_lengths::NTuple{2, T},
        dof_map_m₁::DOFMap,
        dof_map_m₃::DOFMap,
        quad::QuadratureSetup,
        input_data::PDEInputData,
        τ::T, τα::T, τα_q₅::T, _2q₁::T, τq₃::T, τ²α_2::T, τ²q₄α_q₅::T, Δx::T,
        τ_jac_1d::T, τ_jac_2d::T, t_half::T
) where {T, I}
    compute_Q₁!(
        cache.Q₁, matrices.K_m₁xm₁, matrices._2M_m₁xm₁, matrices.M_m₃xm₃, τ²α_2, τ²q₄α_q₅)
    assembly_rhs_2d!(
        cache.τFf₁, (x, y) -> input_data.f₁(x, y, t_half),
        τ_jac_2d, quad.W_φP, dof_map_m₁, nel_per_dim, element_side_lengths, quad.xP, quad.yP)
    assembly_rhs_1d!(
        cache.τFf₃, x -> input_data.f₃(x, t_half),
        τ_jac_1d, quad.W_ϕP, dof_map_m₃, Δx, quad.xP)
    assembly_nonlinearity_F!(
        cache.F, τ, input_data.f, cache.dᵃⁿ, dof_map_m₁, element_side_lengths, quad.φP, quad.W_φP)
    compute_L₁!(
        cache.L₁, state, cache.cᵃⁿ, matrices._2M_m₁xm₁,
        matrices.M_m₃xm₃, matrices.K_m₁xm₁, matrices.τA_m₁xm₂,
        cache.F, cache.τFf₁, cache.τFf₃,
        τα, τα_q₅, _2q₁, τq₃,
        cache.vec_m₁_1, cache.vec_m₁_2, cache.vec_m₃_1, cache.vec_m₃_2)
    solve_nonlinear_system!(
        cache, state, element_side_lengths, dof_map_m₃, quad, input_data, τα)
    return nothing
end

"""
    solve_ĉⁿ(
        cache, state, matrices, nel_per_dim, element_side_lengths,
        dof_map_m₂, quad, input_data, τ, τ_jac_2d, t_half)

Sequentially:
1. Computes `Q₂(cᵃⁿ)`
2. Computes `τ⋅F(f₂(tₙ₋₁/₂))`
3. Computes `L₂(n,v̂ⁿ)`
4. Solves the linear system `Q₂⋅ĉⁿ=L₂`

**WARNING**: Assumes that `cache.cᵃⁿ` and `cache.v̂ⁿ` have already been populated.
The solution `ĉⁿ` of the linear system is stored in `cache.linsolve_m₂.u`.
"""
function solve_ĉⁿ(
        cache::Scheme2Cache{T},
        state::State{T},
        matrices::SystemMatrices{T},
        nel_per_dim::NTuple{2, I},
        element_side_lengths::NTuple{2, T},
        dof_map_m₂::DOFMap,
        quad::QuadratureSetup,
        input_data::PDEInputData,
        τ::T, τ_jac_2d::T, t_half::T) where {T, I}
    bᵀcᵃⁿ = dot(matrices.b, cache.cᵃⁿ)
    τβ_bᵀcᵃⁿ = τ*input_data.β(bᵀcᵃⁿ)
    _muladd!(cache.Q₂_upper, τβ_bᵀcᵃⁿ, matrices.K_m₂xm₂, matrices._2M_m₂xm₂)
    scatter_symmetric!(cache.Q₂, cache.Q₂_upper, cache.map22_direct, cache.map22_mirror)

    assembly_rhs_2d!(
        cache.τFf₂, (x, y) -> input_data.f₂(x, y, t_half),
        τ_jac_2d, quad.W_φP, dof_map_m₂, nel_per_dim, element_side_lengths, quad.xP, quad.yP)

    compute_L₂!(
        cache.L₂, matrices._2M_m₂xm₂, matrices.τA_m₂xm₁, state.c,
        cache.v̂ⁿ, cache.τFf₂, cache.vec_m₂)

    solve_linear_system!(cache.linsolve_m₂)

    return nothing
end

# ==============================================================================
# compute_Q₁!
# ==============================================================================
"""
    compute_Q₁!(Q₁, K_m₁xm₁, _2M_m₁xm₁, M_m₃xm₃, τ²α_2, τ²q₄α_q₅)

Compute in-place

    Q₁ = τ²α_2*K_m₁xm₁ + _2M_m₁xm₁ + τ²q₄α_q₅*[M_m₃xm₃ 0_m₃x(m₁-m₃); 0_(m₁-m₃)xm₃ 0_m₃xm₃]

where `τ²α_2 = τ²⋅α(tₙ₋₁/₂)/2` and `τ²q₄α_q₅ = τ²⋅q₄⋅α(tₙ₋₁/₂)/q₅`
"""
function compute_Q₁!(
        Q₁::Symmetric{T, SparseMatrixCSC{T, I}},
        K_m₁xm₁::Symmetric{T, SparseMatrixCSC{T, I}},
        _2M_m₁xm₁::Symmetric{T, SparseMatrixCSC{T, I}},
        M_m₃xm₃::Symmetric{T, SparseMatrixCSC{T, I}},
        τ²α_2::T, τ²q₄α_q₅::T) where {T, I}
    vec1 = Q₁.data.nzval
    vec2 = K_m₁xm₁.data.nzval
    vec3 = _2M_m₁xm₁.data.nzval
    vec4 = M_m₃xm₃.data.nzval

    nnz_m₁ = length(vec1)
    nnz_m₃ = length(vec4)

    for i in 1:nnz_m₃
        tmp = muladd(τ²α_2, vec2[i], vec3[i])
        vec1[i] = muladd(τ²q₄α_q₅, vec4[i], tmp)
    end
    for i in (nnz_m₃ + 1):nnz_m₁
        vec1[i] = muladd(τ²α_2, vec2[i], vec3[i])
    end

    return nothing
end

# ==============================================================================
# Compute L₁, L₂
# ==============================================================================
"""
    compute_L₁!(
        L₁, state, cᵃⁿ, _2M_m₁xm₁, M_m₃xm₃, K_m₁xm₁, τA_m₁xm₂, τFdᵃⁿ, τFf₁, τFf₃, 
        τα, τα_q₅, _2q₁, τq₃,
        vec_m₁_1, vec_m₁_2, vec_m₃_1, vec_m₃_2)

Compute in-place

    L₁ = 2M_m₁xm₁⋅vⁿ⁻¹ - τα⋅K_m₁xm₁⋅dⁿ⁻¹ - τFdᵃⁿ - τA_m₁xm₂⋅cᵃⁿ  + τFf₁ + τα_q₅⋅[M_m₃xm₃⋅(2q₁⋅rⁿ⁻¹-τq₃⋅zⁿ⁻¹) + τFf₃; 0_(m₁-m₃)]

where `τα = τ⋅α(tₙ₋₁/₂)`, `τα_q₅ = τ⋅α(tₙ₋₁/₂)/q₅`, `τFdᵃⁿ = τF(dᵃⁿ)`, `τFf₁ = τ⋅F(f₁(tₙ₋₁/₂))`, and `τFf₃ = τ⋅F(f₃(tₙ₋₁/₂))`.
"""
function compute_L₁!(
        L₁::Vector{T},
        state::State{T},
        cᵃⁿ::Vector{T},
        _2M_m₁xm₁::Symmetric{T, SparseMatrixCSC{T, I}},
        M_m₃xm₃::Symmetric{T, SparseMatrixCSC{T, I}},
        K_m₁xm₁::Symmetric{T, SparseMatrixCSC{T, I}},
        τA_m₁xm₂::SparseMatrixCSC{T, I},
        τFdᵃⁿ::Vector{T},
        τFf₁::Vector{T},
        τFf₃::Vector{T},
        τα::T, τα_q₅::T, _2q₁::T, τq₃::T,
        vec_m₁_1::Vector{T}, vec_m₁_2::Vector{T}, vec_m₃_1::Vector{T}, vec_m₃_2::Vector{T}) where {
        T, I}
    vⁿ⁻¹ = state.v
    dⁿ⁻¹ = state.d
    rⁿ⁻¹ = state.r
    zⁿ⁻¹ = state.z

    m₁ = length(vⁿ⁻¹)
    m₃ = length(rⁿ⁻¹)

    # L₁ ← 2M_m₁xm₁ · vⁿ⁻¹
    mul!(L₁, _2M_m₁xm₁, vⁿ⁻¹)

    # vec_m₁_1 ← K_m₁xm₁⋅dⁿ⁻¹
    mul!(vec_m₁_1, K_m₁xm₁, dⁿ⁻¹)

    # vec_m₁_2 ← τA_m₁xm₂⋅cᵃⁿ
    mul!(vec_m₁_2, τA_m₁xm₂, cᵃⁿ)

    # vec_m₃_2 ← 2q₁⋅rⁿ⁻¹-τq₃⋅zⁿ⁻¹
    @. vec_m₃_2 = _2q₁ * rⁿ⁻¹ - τq₃ * zⁿ⁻¹

    # vec_m₃_1 ← M_m₃×m₃ · vec_m₃_2
    mul!(vec_m₃_1, M_m₃xm₃, vec_m₃_2)

    # Final assembly
    for i in 1:m₃
        tmp1 = muladd(τα, vec_m₁_1[i], τFdᵃⁿ[i])
        tmp2 = vec_m₃_1[i] + τFf₃[i]
        tmp3 = muladd(τα_q₅, tmp2, τFf₁[i])
        L₁[i] = L₁[i] - tmp1 - vec_m₁_2[i] + tmp3
    end
    for i in (m₃ + 1):m₁
        tmp1 = muladd(τα, vec_m₁_1[i], τFdᵃⁿ[i])
        L₁[i] = L₁[i] - tmp1 - vec_m₁_2[i] + τFf₁[i]
    end

    return nothing
end

"""
    compute_L₂!(L₂, _2M_m₂xm₂, τA_m₂xm₁, cⁿ⁻¹, v̂ⁿ, τFf₂, vec_m₂)

Compute in-place

    L₂ = 2M_m₂xm₂⋅cⁿ⁻¹ - τA_m₂xm₁⋅v̂ⁿ + τFf₂

where `τFf₂ = τ⋅F(f₂(tₙ₋₁/₂))`.
"""
function compute_L₂!(
        L₂::Vector{T},
        _2M_m₂xm₂::Symmetric{T, SparseMatrixCSC{T, I}},
        τA_m₂xm₁::SparseMatrixCSC{T, I},
        cⁿ⁻¹::Vector{T},
        v̂ⁿ::Vector{T},
        τFf₂::Vector{T},
        vec_m₂::Vector{T}) where {T, I}
    mul!(L₂, _2M_m₂xm₂, cⁿ⁻¹)
    mul!(vec_m₂, τA_m₂xm₁, v̂ⁿ)
    @. L₂ = L₂ + τFf₂ - vec_m₂
    return nothing
end

# ==============================================================================
# solve_nonlinear_system!
# ==============================================================================
"""
    solve_nonlinear_system!(
        cache, state, element_side_lengths,
        dof_map_m₃, quad, input_data,
        τα;
        abstol, reltol, maxiter)

Solve the nonlinear system
```math
H(X) = Q₁⋅X + τα⋅G(X) - L₁ = 0
```
for ``X = v̂ⁿ`` via Newton's method, updating `cache.v̂ⁿ` in-place. 
`cache.v̂ⁿ` is warm-started from `vⁿ⁻¹`.

Convergence is declared when either of the following holds at iteration ``k``:
1. Residual criterion:
   ``max|H(Xᵏ⁻¹)| ≤ abstol``
2. Step criterion:
   ``max|Xᵏ-Xᵏ⁻¹| ≤ abstol + reltol ⋅ max|Xᵏ⁻¹|``
"""
function solve_nonlinear_system!(
        cache::Scheme2Cache{T},
        state::State{T},
        element_side_lengths::NTuple{2, T},
        dof_map_m₃::DOFMap,
        quad::QuadratureSetup,
        input_data::PDEInputData,
        τα::T;
        abstol::T = T(1e-14),
        reltol::T = T(1e-9),
        maxiter::Int = 5) where {T}
    Δx, Δy = element_side_lengths

    # Warm start: X ← vⁿ⁻¹
    cache.v̂ⁿ .= state.v

    normH = zero(T)
    normΔX = zero(T)
    normX = zero(T)

    for _ in 1:maxiter
        assembly_nonlinearity_G!(
            cache.G, τα, input_data.g, cache.v̂ⁿ,
            dof_map_m₃, Δx, quad.xP, quad.ϕP, quad.W_ϕP)      # τα⋅G(v̂ⁿ)

        compute_minusH!(
            cache.minusH, cache.v̂ⁿ, cache.Q₁, cache.G, cache.L₁)

        normH = maximum(abs, cache.minusH)
        normH ≤ abstol && return nothing                      # criterion 1

        JG = assembly_global_matrix_DG(
            τα, input_data.∂ₛg, cache.v̂ⁿ, dof_map_m₃,
            Δx, quad.xP, quad.ϕP, quad.W_ϕPϕP)                # τα⋅JG(v̂ⁿ)

        compute_JH_upper!(cache.JH_upper, cache.Q₁, JG)
        scatter_symmetric!(cache.JH, cache.JH_upper, cache.map11_direct, cache.map11_mirror)

        solve_linear_system!(cache.linsolve_m₁)

        normΔX = maximum(abs, cache.linsolve_m₁.u)
        normX = maximum(abs, cache.v̂ⁿ)
        cache.v̂ⁿ .+= cache.linsolve_m₁.u

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

function solve_linear_system!(linsolve)
    linsolve.isfresh = true
    LS.solve!(linsolve)
    return nothing
end

# ==============================================================================
# compute_minusH!
# ==============================================================================
"""
    compute_minusH!(minusH, X, Q₁, G, L₁)

Compute in-place
```math
-H(X) = - Q₁⋅X - G + L₁.
```
where `G = τα(tₙ₋₁/₂)⋅G(X)`.
"""
function compute_minusH!(
        minusH::Vector{T},
        X::Vector{T},
        Q₁::Symmetric{T, SparseMatrixCSC{T, I}},
        G::Vector{T},
        L₁::Vector{T}) where {T, I}
    m₁ = length(L₁)
    m₃ = length(G)

    mul!(minusH, Q₁, X)
    for i in 1:m₃
        minusH[i] = -minusH[i] - G[i] + L₁[i]
    end
    for i in (m₃ + 1):m₁
        minusH[i] = -minusH[i] + L₁[i]
    end
    return nothing
end

# ==============================================================================
# Compute JH
# ==============================================================================
"""
    compute_JH_upper!(JH_upper, Q₁, JG)

Compute in-place

    JH_upper = Q₁ + JG.

where `JG = τα(tₙ₋₁/₂)⋅JG(X)`.
Assumes `JH_upper` and `Q₁` share the same sparsity pattern, and that `JG` occupies
the leading `nnz(JG)` entries of that pattern.
"""
function compute_JH_upper!(
        JH_upper::Symmetric{T, SparseMatrixCSC{T, I}},
        Q₁::Symmetric{T, SparseMatrixCSC{T, I}},
        JG::Symmetric{T, SparseMatrixCSC{T, I}}) where {T, I}
    nnz_m₁ = nnz(JH_upper.data)
    nnz_m₃ = nnz(JG.data)

    nzval_JH = JH_upper.data.nzval
    nzval_Q₁ = Q₁.data.nzval
    nzval_JG = JG.data.nzval

    for i in 1:nnz_m₃
        nzval_JH[i] = nzval_Q₁[i] + nzval_JG[i]
    end
    for i in (nnz_m₃ + 1):nnz_m₁
        nzval_JH[i] = nzval_Q₁[i]
    end
    return nothing
end