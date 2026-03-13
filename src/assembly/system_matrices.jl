# ==============================================================================
# SystemMatrices
# ==============================================================================

"""
    SystemMatrices{T, I}

Global FEM matrices for the coupled thermo-wave-acoustic system.
All matrices are assembled once before the time loop and remain constant
throughout the time integration.

# Type Parameters
- `T <: AbstractFloat`: floating-point precision (e.g. `Float64`).
- `I <: Integer`: integer type used by the sparse matrices (e.g. `Int64`).

# Fields
- `M_m₁xm₁`: mass matrix for the ``m_1``-space (wave field).
- `M_m₂xm₂`: mass matrix for the ``m_2``-space (heat field).
- `M_m₃xm₃`: mass matrix for the ``m_3``-space (acoustic field).
- `K_m₁xm₁`: stiffness matrix for the ``m_1``-space.
- `K_m₂xm₂`: stiffness matrix for the ``m_2``-space.
- `A_m₁xm₂`: coupling matrix with test functions in ``m_1``-space and trial functions
  in ``m_2``-space; arises from the term
  ``(\\varphi, (\\mathbf{a} \\cdot \\nabla)\\Theta^{\\ast n})``.
- `A_m₂xm₁`: coupling matrix with test functions in ``m_2``-space and trial functions
  in ``m_1``-space; arises from the term
  ``(\\psi, (\\mathbf{a} \\cdot \\nabla)\\hat{V}^n)``.
- `b`: global integral vector for the ``m_2``-space,
  ``b_i = \\int_{\\Omega} \\psi_i \\, d\\Omega``, ``i = 1, \\ldots, m_2``;
  see [`assembly_∫basis`](@ref).
"""
struct SystemMatrices{T <: AbstractFloat, I <: Integer}
    M_m₁xm₁::Symmetric{T, SparseMatrixCSC{T, I}}
    M_m₂xm₂::Symmetric{T, SparseMatrixCSC{T, I}}
    M_m₃xm₃::Symmetric{T, SparseMatrixCSC{T, I}}
    K_m₁xm₁::Symmetric{T, SparseMatrixCSC{T, I}}
    K_m₂xm₂::Symmetric{T, SparseMatrixCSC{T, I}}
    A_m₁xm₂::SparseMatrixCSC{T, I}
    A_m₂xm₁::SparseMatrixCSC{T, I}
    b::Vector{T}
end

"""
    SystemMatrices(a, mesh1D, mesh2D, fe1D, fe2D, dof_map_m₁, dof_map_m₂, dof_map_m₃)

Construct the global FEM matrices for the coupled thermo-wave-acoustic system.
See [`SystemMatrices`](@ref) for a description of the assembled fields.

# Arguments
- `a::NTuple{2, T}`: constant advection vector ``\\mathbf{a} = (a_1, a_2)``.
- `mesh1D::CartesianMesh{1}`: 1D Cartesian mesh (boundary ``\\Gamma_1``).
- `mesh2D::CartesianMesh{2}`: 2D Cartesian mesh (domain ``\\Omega``).
- `fe1D::DimensionalFEFamily`: 1D finite element family.
- `fe2D::DimensionalFEFamily`: 2D finite element family.
- `dof_map_m₁::DOFMap`: DOF map for the ``m_1``-space (left, right, and top boundaries prescribed).
- `dof_map_m₂::DOFMap`: DOF map for the ``m_2``-space (all boundaries prescribed).
- `dof_map_m₃::DOFMap`: DOF map for the ``m_3``-space (left and right boundaries prescribed).
"""
function SystemMatrices(
        a::NTuple{2, T},
        mesh1D::CartesianMesh{1},
        mesh2D::CartesianMesh{2},
        fe1D::F1,
        fe2D::F2,
        dof_map_m₁::DOFMap,
        dof_map_m₂::DOFMap,
        dof_map_m₃::DOFMap
) where {F1 <: DimensionalFEFamily, F2 <: DimensionalFEFamily, T <: AbstractFloat}
    # Local matrices — 2D
    Me_m₁xm₁ = Symmetric(assembly_local_matrix_ϕxϕ(mesh2D, fe2D))
    Ke_m₁xm₁ = Symmetric(assembly_local_matrix_∇ϕx∇ϕ(mesh2D, fe2D))
    Ae = assembly_local_matrix_ϕxc∇ϕ(a, mesh2D, fe2D)

    # Local matrices — 1D
    Me_m₃xm₃ = Symmetric(assembly_local_matrix_ϕxϕ(mesh1D, fe1D))

    # Local integral vector — m₂-space
    be = assembly_local_∫basis(fe2D, mesh2D)

    return SystemMatrices(
        # Mass matrices
        assembly_global_matrix(Me_m₁xm₁, dof_map_m₁),        # M_m₁xm₁
        assembly_global_matrix(Me_m₁xm₁, dof_map_m₂),        # M_m₂xm₂
        assembly_global_matrix(Me_m₃xm₃, dof_map_m₃),        # M_m₃xm₃
        # Stiffness matrices
        assembly_global_matrix(Ke_m₁xm₁, dof_map_m₁),        # K_m₁xm₁
        assembly_global_matrix(Ke_m₁xm₁, dof_map_m₂),        # K_m₂xm₂
        # Coupling matrices
        assembly_global_matrix(Ae, dof_map_m₁, dof_map_m₂),  # A_m₁xm₂
        assembly_global_matrix(Ae, dof_map_m₂, dof_map_m₁),  # A_m₂xm₁
        # Integral vector
        assembly_∫basis(be, dof_map_m₂)                      # b
    )
end