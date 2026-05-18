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

@testitem "Synchronization problem - case study 1" begin
    using ThermoelasticAcoustic: CartesianMesh, DOFMap, Lagrange, LeftRight,
                                 LeftRightBottomTop, LeftRightTop, specialize,
                                 SystemMatrices
    using LinearAlgebra: UpperTriangular
    using SparseArrays: sparse, nnz, nzrange

    # --- Setup ---
    pmin, pmax = (0.0, 0.0), (1.0, 1.0)
    mesh1D = CartesianMesh((pmin[1],), (pmax[1],), (16,))
    mesh2D = CartesianMesh(pmin, pmax, (16, 16))
    fe1D = specialize(Lagrange{1}(), Val(1))
    fe2D = specialize(Lagrange{1}(), Val(2))

    dof_map_m₁ = DOFMap(mesh2D, fe2D, LeftRightTop())
    dof_map_m₂ = DOFMap(mesh2D, fe2D, LeftRightBottomTop())
    dof_map_m₃ = DOFMap(mesh1D, fe1D, LeftRight())

    matrices = SystemMatrices(
        (1.0, 1.0), mesh1D, mesh2D, fe1D, fe2D,
        dof_map_m₁, dof_map_m₂, dof_map_m₃)

    # ==========================================================================
    # Problem statement
    #
    # A : structurally symmetric sparse matrix — both triangles explicitly stored.
    # B : upper triangle of A, i.e. B = sparse(UpperTriangular(A)).
    # C : uninitialized sparse matrix with the same sparsity pattern as A.
    #
    # Goal: populate C.nzval from B.nzval such that C == A:
    #     C[i,j] = B[i,j]    for i ≤ j   (direct write)
    #     C[j,i] = B[i,j]    for i < j   (mirror write)
    #
    # Two strategies are compared:
    #   1. Naive     : locates the nzval positions of C[i,j] and C[j,i] by
    #                  linear search at every call — no preprocessing.
    #   2. Precomputed maps : same linear search performed once during
    #                  preprocessing; the hot loop uses the stored positions
    #                  directly — zero searches, zero allocations per call.
    # ==========================================================================
    A = sparse(matrices.M_m₁xm₁)
    B = sparse(UpperTriangular(A))
    C = similar(A)

    # ==========================================================================
    # Strategy 1 — naive linear search, no preprocessing
    #
    # For each entry (i,j) in B:
    #   - scans column j of C to locate C[i,j]   (direct position)
    #   - scans column i of C to locate C[j,i]   (mirror position)
    # ==========================================================================
    function sync_naive!(C, B)
        n = size(B, 2)
        for j in 1:n                          # Iterate over columns of B 
            for kb in nzrange(B, j)           # Iterate over nonzeros of B[:,j]
                i = B.rowval[kb]
                v = B.nzval[kb]
                for kc in nzrange(C, j)       # Find position of C[i,j]
                    if C.rowval[kc] == i
                        C.nzval[kc] = v
                        break
                    end
                end
                if i != j
                    for kcᵀ in nzrange(C, i)  # Find position of C[j,i]
                        if C.rowval[kcᵀ] == j
                            C.nzval[kcᵀ] = v
                            break
                        end
                    end
                end
            end
        end
        return nothing
    end

    # ==========================================================================
    # Strategy 2 — precomputed index maps
    #
    # Preprocessing: performs the same column scans as Strategy 1, but only once.
    # Stores the nzval positions in two index vectors:
    #   map_direct[kb] → position kc  in C.nzval of entry C[i,j]
    #   map_mirror[kb] → position kcᵀ in C.nzval of entry C[j,i]
    #
    # Hot loop (sync!): one sequential pass over B.nzval — zero allocations.
    # ==========================================================================
    function build_maps(C, B)
        nnz_B = nnz(B)
        map_direct = Vector{Int}(undef, nnz_B)
        map_mirror = Vector{Int}(undef, nnz_B)
        for j in 1:size(B, 2)
            for kb in nzrange(B, j)
                i = B.rowval[kb]
                for kc in nzrange(C, j)
                    if C.rowval[kc] == i
                        map_direct[kb] = kc
                        break
                    end
                end
                for kcᵀ in nzrange(C, i)
                    if C.rowval[kcᵀ] == j
                        map_mirror[kb] = kcᵀ
                        break
                    end
                end
            end
        end
        return map_direct, map_mirror
    end

    function sync!(C, B, map_direct, map_mirror)
        @inbounds for kb in eachindex(map_direct)
            v = B.nzval[kb]
            kc = map_direct[kb]
            kcᵀ = map_mirror[kb]
            C.nzval[kc] = v
            kc != kcᵀ && (C.nzval[kcᵀ] = v)   # skip diagonal (kc == kcᵀ)
        end
        return nothing
    end

    # ==========================================================================
    # Validation
    # ==========================================================================
    sync_naive!(C, B)
    @test C == A

    C = similar(A)
    map_direct, map_mirror = build_maps(C, B)
    alloc = @allocated sync!(C, B, map_direct, map_mirror)
    @test C == A
    @test alloc == 0
end

@testitem "Synchronization problem - case study 2" begin
    using ThermoelasticAcoustic: CartesianMesh, DOFMap, Lagrange, LeftRight,
                                 LeftRightBottomTop, LeftRightTop, specialize,
                                 SystemMatrices
    using LinearAlgebra: UpperTriangular
    using SparseArrays: sparse, nnz, nzrange

    function setup(fe)
        pmin, pmax = (0.0, 0.0), (1.0, 1.0)
        mesh1D = CartesianMesh((pmin[1],), (pmax[1],), (4,))
        mesh2D = CartesianMesh(pmin, pmax, (4, 4))
        fe1D = specialize(fe, Val(1))
        fe2D = specialize(fe, Val(2))

        dof_map_m₁ = DOFMap(mesh2D, fe2D, LeftRightTop())
        dof_map_m₂ = DOFMap(mesh2D, fe2D, LeftRightBottomTop())
        dof_map_m₃ = DOFMap(mesh1D, fe1D, LeftRight())

        matrices = SystemMatrices(
            (1.0, 1.0), mesh1D, mesh2D, fe1D, fe2D,
            dof_map_m₁, dof_map_m₂, dof_map_m₃)

        M̄₁₁ = sparse(matrices.M_m₁xm₁)
        M̄₂₂ = sparse(matrices.M_m₂xm₂)
        M₁₁ = sparse(UpperTriangular(M̄₁₁))
        M₂₂ = sparse(UpperTriangular(M̄₂₂))
        A₁₂ = matrices.A_m₁xm₂
        A₂₁ = matrices.A_m₂xm₁
        J = [similar(M̄₁₁) A₁₂; A₂₁ similar(M̄₂₂)]

        return J, M̄₁₁, M̄₂₂, M₁₁, M₂₂, A₁₂, A₂₁, dof_map_m₁.m
    end

    # ==========================================================================
    # Problem statement
    #
    # M̄₁₁, M̄₂₂ : structurally symmetric sparse matrix — both triangles stored.
    # M₁₁, M₂₂ : upper triangles, i.e. M₁₁ = sparse(UpperTriangular(M̄₁₁)), and analogously for M₂₂.
    # A₁₂, A₂₁  : general (non-symmetric) sparse matrices.
    #
    # The block Jacobian J has the structure:
    #
    #   J = [ similar(M̄₁₁)   A₁₂ ;
    #         A₂₁            similar(M̄₂₂) ]
    #
    # where the off-diagonal blocks are already populated and must not be modified.
    # Only the diagonal blocks of J need to be filled.
    #
    # Goal: populate J.nzval from M₁₁.nzval and M₂₂.nzval such that J == [M̄₁₁ A₁₂; A₂₁ M̄₂₂].
    #
    # The solution extends the precomputed-maps strategy from case study 1:
    # build index maps from M₁₁ → J and M₂₂ → J during preprocessing, then populate J.nzval.
    # ==========================================================================

    # Maps M₁₁.nzval positions → J.nzval positions for the (1,1) block.
    # Identical to case study 1: column indices in M₁₁ and J coincide.
    function build_maps_11(J, M₁₁)
        nnz_M₁₁ = nnz(M₁₁)
        map_direct = Vector{Int}(undef, nnz_M₁₁)
        map_mirror = Vector{Int}(undef, nnz_M₁₁)
        for j in 1:size(M₁₁, 2)        # Iterate over columns of M₁₁
            for kM in nzrange(M₁₁, j)  # Iterate over nonzeros of M₁₁[:,j]
                i = M₁₁.rowval[kM]
                for kJ in nzrange(J, j)# Iterate over nonzeros of J[:,j] and find position of J[i,j]
                    if J.rowval[kJ] == i
                        map_direct[kM] = kJ
                        break
                    end
                end
                for kᵀJ in nzrange(J, i)# Iterate over nonzeros of J[:,i] and find position of J[j,i]
                    if J.rowval[kᵀJ] == j
                        map_mirror[kM] = kᵀJ
                        break
                    end
                end
            end
        end
        return map_direct, map_mirror
    end

    # Maps M₂₂.nzval positions → J.nzval positions for the (2,2) block.
    #
    # Two index shifts are needed relative to case study 1:
    #   1. Column shift: column j of M₂₂ → column j+m₁ of J.
    #   2. Row offset within column: each column j+m₁ of J starts with
    #      length(nzrange(A₁₂, j)) entries from A₁₂ before the M̄₂₂ entries.
    #      start[j] skips those A₁₂ entries when searching for M̄₂₂[i,j] in J.
    function build_maps_22(J, M₂₂, A₁₂)
        m₁ = size(A₁₂, 1)

        # start[j]: index within nzrange(J, j+m₁) where the M̄₂₂ entries begin,
        # i.e. right after the A₁₂[:,j] entries.
        start = [length(nzrange(A₁₂, j)) + 1 for j in 1:size(A₁₂, 2)]

        nnz_M₂₂ = nnz(M₂₂)
        map_direct = Vector{Int}(undef, nnz_M₂₂)
        map_mirror = Vector{Int}(undef, nnz_M₂₂)
        for j in 1:size(M₂₂, 2)                         # Iterate over columns of M₂₂
            jJ = j + m₁
            for kM in nzrange(M₂₂, j)                   # Iterate over nonzeros of M₂₂[:,j]
                i = M₂₂.rowval[kM]
                iJ = i + m₁
                for kJ in nzrange(J, jJ)[start[j]:end]  # find J[i+m₁, j+m₁]
                    if J.rowval[kJ] == iJ
                        map_direct[kM] = kJ
                        break
                    end
                end
                for kᵀJ in nzrange(J, iJ)[start[i]:end] # find J[j+m₁, i+m₁]
                    if J.rowval[kᵀJ] == jJ
                        map_mirror[kM] = kᵀJ
                        break
                    end
                end
            end
        end
        return map_direct, map_mirror
    end

    function sync!(J, M, map_direct, map_mirror)
        @inbounds for kM in eachindex(map_direct)
            v = M.nzval[kM]
            kJ = map_direct[kM]
            kᵀJ = map_mirror[kM]
            J.nzval[kJ] = v
            kJ != kᵀJ && (J.nzval[kᵀJ] = v)   # skip diagonal (kJ == kᵀJ)
        end
        return nothing
    end

    # ==========================================================================
    # Validation
    # ==========================================================================
    @testset "$fe" for fe in (Lagrange{1}(), Lagrange{2}(), Lagrange{3}())
        J, M̄₁₁, M̄₂₂, M₁₁, M₂₂, A₁₂, A₂₁, m₁ = setup(fe)

        map_direct₁₁, map_mirror₁₁ = build_maps_11(J, M₁₁)
        map_direct₂₂, map_mirror₂₂ = build_maps_22(J, M₂₂, A₁₂)

        alloc1 = @allocated sync!(J, M₁₁, map_direct₁₁, map_mirror₁₁)
        alloc2 = @allocated sync!(J, M₂₂, map_direct₂₂, map_mirror₂₂)

        @test J[1:m₁, 1:m₁] == M̄₁₁
        @test J[(m₁ + 1):end, (m₁ + 1):end] == M̄₂₂
        @test J == [M̄₁₁ A₁₂; A₂₁ M̄₂₂]
        @test alloc1 == 0
        @test alloc2 == 0
    end
end