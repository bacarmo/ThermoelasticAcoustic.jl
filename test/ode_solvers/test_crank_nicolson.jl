@testitem "compute_Q!: CrankNicolson solver, Q matrix assembly" begin
    using ThermoelasticAcoustic: build_cache, CartesianMesh, compute_Q!, CrankNicolson,
                                 DOFMap, Lagrange, LeftRight, LeftRightBottomTop,
                                 LeftRightTop, specialize, SystemMatrices
    using LinearAlgebra: Symmetric
    using SparseArrays: SparseMatrixCSC

    # Setup
    pmin, pmax = (0.0, 0.0), (1.0, 1.0)
    Nx = (4, 3)
    τ = 1.0
    q₄ = 1.0
    q₅ = 1.0
    a = (1.0, 1.0)

    mesh1D = CartesianMesh((pmin[1],), (pmax[1],), (Nx[1],))
    mesh2D = CartesianMesh(pmin, pmax, Nx)

    fe1D = specialize(Lagrange{1}(), Val(1))
    fe2D = specialize(Lagrange{1}(), Val(2))

    dof_map_m₁ = DOFMap(mesh2D, fe2D, LeftRightTop())
    dof_map_m₂ = DOFMap(mesh2D, fe2D, LeftRightBottomTop())
    dof_map_m₃ = DOFMap(mesh1D, fe1D, LeftRight())

    m₁ = dof_map_m₁.m
    m₂ = dof_map_m₂.m
    m₃ = dof_map_m₃.m

    matrices = SystemMatrices(
        a, mesh1D, mesh2D, fe1D, fe2D, dof_map_m₁, dof_map_m₂, dof_map_m₃)
    cache = build_cache(CrankNicolson(), matrices, dof_map_m₁, dof_map_m₂, dof_map_m₃, τ)

    # Test α = 0: Q[1:m₁, 1:m₁] = 2M_m₁×m₁
    compute_Q!(cache, matrices, τ, 0.0, q₄, q₅)
    @test cache.M_m₁xm₁_vs2 ≈ matrices.M_m₁xm₁ * 2
    @test cache.Q[1:m₁, 1:m₁] ≈ Matrix(cache.M_m₁xm₁_vs2)

    # Test α = 1, q₄ = 0: Q[1:m₁, 1:m₁] = 2M_m₁×m₁ + (τ²/2)K_m₁×m₁
    compute_Q!(cache, matrices, τ, 1.0, 0.0, q₅)
    @test cache.Q[1:m₁, 1:m₁] ≈ Matrix(cache.M_m₁xm₁_vs2 + (τ^2 / 2) * matrices.K_m₁xm₁)

    # Test α = 1, q₄ = 1: Q[1:m₁, 1:m₁] = 2M_m₁×m₁ + (τ²/2)K_m₁×m₁ + (τ²q₄/q₅)[M_m₃×m₃  0; 0  0]
    compute_Q!(cache, matrices, τ, 1.0, q₄, q₅)
    Q_block_expected = Matrix(cache.M_m₁xm₁_vs2 + (τ^2 / 2) * matrices.K_m₁xm₁)
    Q_block_expected[1:m₃, 1:m₃] += Matrix(matrices.M_m₃xm₃ * (τ^2 * q₄ / q₅))
    @test cache.Q[1:m₁, 1:m₁] ≈ Q_block_expected

    # Test time-invariant blocks are preserved after compute_Q!
    @test cache.Q[1:m₁, (m₁ + 1):end] ≈ Matrix(τ * matrices.A_m₁xm₂)
    @test cache.Q[(m₁ + 1):end, 1:m₁] ≈ Matrix(τ * matrices.A_m₂xm₁)
    @test cache.Q[(m₁ + 1):end, (m₁ + 1):end] ≈ Matrix(cache.M_m₂xm₂_vs2)

    # Test matrix dimensions
    @test size(cache.Q) == (m₁ + m₂, m₁ + m₂)

    # Test type correctness
    @test cache.Q isa SparseMatrixCSC{Float64, Int64}

    # Test allocation behavior (FIXME: sparse block assignments currently allocate)
    @test (@allocated compute_Q!(cache, matrices, τ, 1.0, q₄, q₅)) > 0
end