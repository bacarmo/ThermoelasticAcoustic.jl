@testitem "compute_Q₁!: ModifiedCN solver, Q₁ matrix assembly" begin
    using ThermoelasticAcoustic: CartesianMesh, DOFMap, LeftRight, LeftRightBottomTop,
                                 LeftRightTop, Lagrange, ModifiedCN, SystemMatrices,
                                 build_cache, compute_Q₁!, specialize
    using SparseArrays: SparseMatrixCSC
    using LinearAlgebra: Symmetric

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

    m₃ = dof_map_m₃.m

    matrices = SystemMatrices(
        a, mesh1D, mesh2D, fe1D, fe2D, dof_map_m₁, dof_map_m₂, dof_map_m₃)
    cache = build_cache(ModifiedCN(), matrices, dof_map_m₁, dof_map_m₂, dof_map_m₃, τ)

    # Test α = 0: Q₁ = 2M_m₁×m₁
    compute_Q₁!(cache, matrices, τ, 0.0, q₄, q₅)
    @test cache.M_m₁xm₁_vs2 ≈ matrices.M_m₁xm₁ * 2
    @test cache.Q₁ ≈ cache.M_m₁xm₁_vs2

    # Test α = 1, q₄ = 0: Q₁ = 2M_m₁×m₁ + (τ²/2)K_m₁×m₁
    compute_Q₁!(cache, matrices, τ, 1.0, 0.0, q₅)
    @test cache.Q₁ ≈ cache.M_m₁xm₁_vs2 + (τ^2 / 2) * matrices.K_m₁xm₁

    # Test α = 1, q₄ = 1: Q₁ = 2M_m₁×m₁ + (τ²/2)K_m₁×m₁ + (τ²q₄/q₅)[M_m₃×m₃  0; 0  0]
    # Expected is built independently via dense indexation, bypassing .data access
    compute_Q₁!(cache, matrices, τ, 1.0, q₄, q₅)
    Q₁_expected = Matrix(cache.M_m₁xm₁_vs2 + (τ^2 / 2) * matrices.K_m₁xm₁)
    Q₁_expected[1:m₃, 1:m₃] += Matrix(matrices.M_m₃xm₃) * (τ^2 * q₄ / q₅)
    @test cache.Q₁ ≈ Q₁_expected

    # Test type correctness
    # cache.Q₁ is a Symmetric wrapper around a SparseMatrixCSC storing the upper triangle.
    # The production code uses .data.nzval to manipulate nonzeros directly,
    # bypassing the Symmetric wrapper and its index protection.
    @test cache.Q₁ isa Symmetric{Float64, SparseMatrixCSC{Float64, Int64}}
    @test cache.Q₁.data isa SparseMatrixCSC{Float64, Int64}
    @test cache.Q₁.data.nzval isa Vector{Float64}

    # Test allocation-free operation
    @test (@allocated compute_Q₁!(cache, matrices, τ, 1.0, q₄, q₅)) == 0
end

@testitem "compute_JH!: ModifiedCN solver, JH matrix assembly" begin
    using ThermoelasticAcoustic: assembly_global_matrix_DG, build_cache, CartesianMesh,
                                 compute_JH!, compute_Q₁!, DOFMap, example1_manufactured,
                                 Lagrange, LeftRight, LeftRightBottomTop, LeftRightTop,
                                 ModifiedCN, QuadratureSetup, specialize, SystemMatrices
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
    m₃ = dof_map_m₃.m

    matrices = SystemMatrices(
        a, mesh1D, mesh2D, fe1D, fe2D, dof_map_m₁, dof_map_m₂, dof_map_m₃)
    cache = build_cache(ModifiedCN(), matrices, dof_map_m₁, dof_map_m₂, dof_map_m₃, τ)

    quad = QuadratureSetup(fe1D, fe2D, mesh2D.Δx, pmin)
    input_data = example1_manufactured()

    # Test τα = 0: JH = Q₁ = 2M_m₁×m₁ (with α = 0)
    @. cache.v̂ⁿ = 1.0
    compute_Q₁!(cache, matrices, τ, 0.0, q₄, q₅)
    compute_JH!(cache, dof_map_m₃, mesh1D, quad, 0.0, input_data)
    @test cache.JH ≈ cache.Q₁
    @test cache.JH ≈ cache.M_m₁xm₁_vs2

    # Test τα = τ, α = 1: JH = Q₁ + τ[JG  0; 0  0]
    # Expected is built independently via dense indexation, bypassing .data access
    compute_Q₁!(cache, matrices, τ, 1.0, q₄, q₅)

    JG = assembly_global_matrix_DG(
        τ, input_data.∂ₛg, @view(cache.v̂ⁿ[1:m₃]), mesh1D, dof_map_m₃, quad)

    compute_JH!(cache, dof_map_m₃, mesh1D, quad, τ, input_data)

    JH_expected = Matrix(cache.M_m₁xm₁_vs2 + (τ^2 / 2) * matrices.K_m₁xm₁)
    JH_expected[1:m₃, 1:m₃] += Matrix(matrices.M_m₃xm₃) * (τ^2 * q₄ / q₅)
    JH_expected[1:m₃, 1:m₃] += Matrix(JG)

    @test cache.JH ≈ JH_expected

    # Test type correctness
    # cache.JH is a Symmetric wrapper around a SparseMatrixCSC storing the upper triangle.
    # The production code uses .data.nzval to manipulate nonzeros directly,
    # bypassing the Symmetric wrapper and its index protection.
    @test cache.JH isa Symmetric{Float64, SparseMatrixCSC{Float64, Int64}}
    @test cache.JH.data isa SparseMatrixCSC{Float64, Int64}
    @test cache.JH.data.nzval isa Vector{Float64}

    # Test upper triangular storage: lower triangle entries should be zero
    S = cache.JH.data
    @test all(S[i, j] == 0 for j in 1:(m₁ - 1) for i in (j + 1):m₁)
    @test size(cache.JH) == (m₁, m₁)

    # Test allocation behavior: assembly_global_matrix_DG allocates (FIXME)
    @test (@allocated compute_JH!(cache, dof_map_m₃, mesh1D, quad, τ, input_data)) > 0
end