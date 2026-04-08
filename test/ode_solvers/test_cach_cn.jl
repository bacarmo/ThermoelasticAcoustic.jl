# Key design facts verified here:
#   1. `linsolve.A` and `JH_sparse` share the same `nzval` array (no internal copy).
#   2. In-place mutation of `JH_sparse.nzval` is immediately visible to the solver.
#   3. `isfresh = true` triggers refactorization; `isfresh = false` suppresses it.
@testitem "cache_cn: KLU LinearCache.isfresh behavior" begin
    using ThermoelasticAcoustic: build_cache, CartesianMesh, CrankNicolson,
                                 DOFMap, Lagrange, LeftRight, LeftRightBottomTop,
                                 LeftRightTop, specialize, SystemMatrices
    import LinearSolve as LS

    # --- Setup ---
    pmin, pmax = (0.0, 0.0), (1.0, 1.0)
    mesh1D = CartesianMesh((pmin[1],), (pmax[1],), (4,))
    mesh2D = CartesianMesh(pmin, pmax, (4, 3))
    fe1D = specialize(Lagrange{1}(), Val(1))
    fe2D = specialize(Lagrange{1}(), Val(2))

    dof_map_m₁ = DOFMap(mesh2D, fe2D, LeftRightTop())
    dof_map_m₂ = DOFMap(mesh2D, fe2D, LeftRightBottomTop())
    dof_map_m₃ = DOFMap(mesh1D, fe1D, LeftRight())

    matrices = SystemMatrices(
        (1.0, 1.0), mesh1D, mesh2D, fe1D, fe2D,
        dof_map_m₁, dof_map_m₂, dof_map_m₃)
    cache = build_cache(
        CrankNicolson(), matrices, dof_map_m₁, dof_map_m₂, dof_map_m₃, 1.0)

    # Fact 1: the solver holds a reference to JH_sparse, not a copy.
    # In-place mutations to JH_sparse.nzval are immediately visible to the solver.
    @test cache.linsolve.A === cache.JH_sparse
    @test cache.linsolve.A.nzval === cache.JH_sparse.nzval

    # --- Baseline: solve A x = b with isfresh = true (first call always factorizes) ---
    cache.linsolve.b .= 1.0
    cache.linsolve.isfresh = true
    LS.solve!(cache.linsolve)
    sol_A = copy(cache.linsolve.u)
    sol_A_default = cache.JH_sparse \ cache.linsolve.b
    @test sol_A ≈ sol_A_default

    # --- Mutate A → 2A in-place; set isfresh = true to trigger refactorization ---
    # Expected: (2A) x = b  ⟹  x = sol_A / 2
    cache.JH_sparse.nzval .*= 2.0
    cache.linsolve.b .= 1.0
    cache.linsolve.isfresh = true
    LS.solve!(cache.linsolve)
    sol_2A = copy(cache.linsolve.u)
    @test sol_2A ≈ sol_A / 2

    # --- Change only b; set isfresh = false to reuse the factorization of 2A ---
    # Expected: (2A) x = 2b  ⟹  x = A⁻¹b = sol_A
    cache.linsolve.b .= 2.0
    cache.linsolve.isfresh = false
    LS.solve!(cache.linsolve)
    sol_2A_2b = copy(cache.linsolve.u)
    @test sol_2A_2b ≈ sol_A
end