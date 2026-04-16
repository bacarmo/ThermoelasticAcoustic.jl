@testitem "compute_Q₁₁!: CrankNicolson solver, Q₁₁ matrix assembly" begin
    using ThermoelasticAcoustic: build_cache, CartesianMesh, compute_Q₁₁!, CrankNicolson,
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

    # Test α = 0: Q₁₁ = 2M_m₁×m₁
    compute_Q₁₁!(cache, matrices, τ, 0.0, q₄, q₅)
    @test cache.M_m₁xm₁_vs2 ≈ matrices.M_m₁xm₁ * 2
    @test cache.Q₁₁ ≈ Matrix(cache.M_m₁xm₁_vs2)

    # Test α = 1, q₄ = 0: Q₁₁ = 2M_m₁×m₁ + (τ²/2)K_m₁×m₁
    compute_Q₁₁!(cache, matrices, τ, 1.0, 0.0, q₅)
    @test cache.Q₁₁ ≈ Matrix(cache.M_m₁xm₁_vs2 + (τ^2 / 2) * matrices.K_m₁xm₁)

    # Test α = 1, q₄ = 1: Q₁₁ = 2M_m₁×m₁ + (τ²/2)K_m₁×m₁ + (τ²q₄/q₅)[M_m₃×m₃  0; 0  0]
    compute_Q₁₁!(cache, matrices, τ, 1.0, q₄, q₅)
    Q₁₁_expected = Matrix(cache.M_m₁xm₁_vs2 + (τ^2 / 2) * matrices.K_m₁xm₁)
    Q₁₁_expected[1:m₃, 1:m₃] += Matrix(matrices.M_m₃xm₃ * (τ^2 * q₄ / q₅))
    @test cache.Q₁₁ ≈ Q₁₁_expected

    # Test matrix dimensions
    @test size(cache.Q₁₁) == (m₁, m₁)

    # Test type correctness
    @test cache.Q₁₁ isa Symmetric{Float64, SparseMatrixCSC{Float64, Int64}}

    # Test allocation behavior
    @test (@allocated compute_Q₁₁!(cache, matrices, τ, 1.0, q₄, q₅)) == 0
end

@testitem "compute_JH_sparse!: CrankNicolson solver, JH_sparse matrix assembly" begin
    using ThermoelasticAcoustic: assembly_global_matrix_DF, assembly_global_matrix_DG,
                                 build_cache, CartesianMesh, compute_JH_sparse!,
                                 compute_Q₁₁!,
                                 CrankNicolson, DOFMap, example1_manufactured, Lagrange,
                                 LeftRight, LeftRightBottomTop, LeftRightTop,
                                 QuadratureSetup, specialize, SystemMatrices
    using LinearAlgebra: dot, Symmetric, mul!
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

    quad = QuadratureSetup(fe1D, fe2D, mesh2D.Δx, pmin)
    input_data = example1_manufactured()

    @. cache.v̂ⁿ = 1.0
    @. cache.ĉⁿ = 1.0
    @. cache.d̂ⁿ = 1.0
    mul!(cache.K_m₂xm₂_vs_ĉⁿ, matrices.K_m₂xm₂, cache.ĉⁿ)

    # Test τ = α = 0: JH_sparse[1:m₁, 1:m₁] = 2M_m₁xm₁ (all τ-dependent terms vanish)
    compute_Q₁₁!(cache, matrices, 0.0, 0.0, q₄, q₅)
    compute_JH_sparse!(cache, matrices, dof_map_m₁, dof_map_m₂, dof_map_m₃,
        mesh1D, mesh2D, quad, 0.0, 0.0, input_data)
    @test cache.JH_sparse[1:m₁, 1:m₁] ≈ Matrix(cache.M_m₁xm₁_vs2)

    # Test α=1: JH_sparse[1:m₁, 1:m₁] = Q₁₁ + τα JG + (τ²/2) JF
    # Expected built independently via dense indexation
    compute_Q₁₁!(cache, matrices, τ, 1.0, q₄, q₅)

    JG = assembly_global_matrix_DG(
        τ, input_data.∂ₛg, @view(cache.v̂ⁿ[1:m₃]), mesh1D, dof_map_m₃, quad)
    JF = assembly_global_matrix_DF(
        τ^2 / 2, input_data.df, cache.d̂ⁿ, mesh2D, dof_map_m₁, quad)

    compute_JH_sparse!(cache, matrices, dof_map_m₁, dof_map_m₂, dof_map_m₃,
        mesh1D, mesh2D, quad, τ, τ, input_data)

    JH_block_expected = Matrix(cache.Q₁₁)
    JH_block_expected[1:m₃, 1:m₃] += Matrix(JG)
    JH_block_expected[1:m₁, 1:m₁] += Matrix(JF)

    @test cache.JH_sparse[1:m₁, 1:m₁] ≈ JH_block_expected

    # Test time-invariant off-diagonal blocks are preserved
    @test cache.JH_sparse[1:m₁, (m₁ + 1):end] ≈ Matrix(τ * matrices.A_m₁xm₂)
    @test cache.JH_sparse[(m₁ + 1):end, 1:m₁] ≈ Matrix(τ * matrices.A_m₂xm₁)

    # Test bottom-right m₂×m₂ block: 2M_m₂×m₂ + τβ K_m₂×m₂
    b_dot_ĉⁿ = dot(matrices.b, cache.ĉⁿ)
    @test b_dot_ĉⁿ ≈ prod(mesh2D.Δx) * m₂

    τβ = τ * input_data.β(b_dot_ĉⁿ)

    JH_β_expected = Matrix(cache.M_m₂xm₂_vs2) +
                    τβ * Matrix(matrices.K_m₂xm₂)
    @test cache.JH_sparse[(m₁ + 1):end, (m₁ + 1):end] ≈ JH_β_expected

    # Test matrix dimensions
    @test size(cache.JH_sparse) == (m₁ + m₂, m₁ + m₂)

    # Test type correctness
    @test cache.JH_sparse isa SparseMatrixCSC{Float64, Int64}

    # Test allocation behavior (FIXME: JG and JF allocate)
    @test (@allocated compute_JH_sparse!(
        cache, matrices, dof_map_m₁, dof_map_m₂, dof_map_m₃,
        mesh1D, mesh2D, quad, τ, τ, input_data)) > 0
end