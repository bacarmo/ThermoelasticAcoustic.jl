"""
    SystemMatrices{T, I}

Global FEM matrices for the coupled thermo-wave-acoustic system.

## Fields
- `M_m₁xm₁`, `M_m₂xm₂`, `M_m₃xm₃`: Mass matrices
- `K_m₁xm₁`, `K_m₂xm₂`: Stiffness matrices
"""
struct SystemMatrices{T <: AbstractFloat, I <: Integer}
    M_m₁xm₁::Symmetric{T, SparseMatrixCSC{T, I}}
    M_m₂xm₂::Symmetric{T, SparseMatrixCSC{T, I}}
    M_m₃xm₃::Symmetric{T, SparseMatrixCSC{T, I}}
    K_m₁xm₁::Symmetric{T, SparseMatrixCSC{T, I}}
    K_m₂xm₂::Symmetric{T, SparseMatrixCSC{T, I}}
    A_m₁xm₂::SparseMatrixCSC{T, I}
    A_m₂xm₁::SparseMatrixCSC{T, I}
end

"""
    SystemMatrices(a, mesh1D, mesh2D, fe1D, fe2D, dof_map_m₁, dof_map_m₂, dof_map_m₃)

Assemble and return the global FEM matrices.

## Arguments
- `a::NTuple{2,T}`: Constant vector
- `mesh1D::CartesianMesh{1}`: 1D Cartesian mesh
- `mesh2D::CartesianMesh{2}`: 2D Cartesian mesh
- `fe1D::F1`: Specialized 1D finite element family
- `fe2D::F2`: Specialized 2D finite element family
- `dof_map_m₁::DOFMap`: DOF mapping with `LeftRightTop` boundary
- `dof_map_m₂::DOFMap`: DOF mapping with `LeftRightBottomTop` boundary
- `dof_map_m₃::DOFMap`: DOF mapping with `LeftRight` boundary
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
    Me_m₁xm₁ = Symmetric(assembly_local_matrix_ϕxϕ(mesh2D, fe2D))
    Me_m₃xm₃ = Symmetric(assembly_local_matrix_ϕxϕ(mesh1D, fe1D))
    Ke_m₁xm₁ = Symmetric(assembly_local_matrix_∇ϕx∇ϕ(mesh2D, fe2D))
    Ae = assemble_local_matrix_ϕxc∇ϕ(a, mesh2D, fe2D)

    return SystemMatrices(
        assembly_global_matrix(Me_m₁xm₁, dof_map_m₁),
        assembly_global_matrix(Me_m₁xm₁, dof_map_m₂),
        assembly_global_matrix(Me_m₃xm₃, dof_map_m₃),
        assembly_global_matrix(Ke_m₁xm₁, dof_map_m₁),
        assembly_global_matrix(Ke_m₁xm₁, dof_map_m₂),
        assembly_global_matrix(Ae, dof_map_m₁, dof_map_m₂),
        assembly_global_matrix(Ae, dof_map_m₂, dof_map_m₁)
    )
end