# ==============================================================================
# assembly_local_∫basis
# ==============================================================================

"""
    assembly_local_∫basis(fe, mesh)
 
Compute the local integral vector
```math
b^e_a = \\int_{\\Omega^e} \\varphi_a^e(x, y) \\, d\\Omega,
\\quad a = 1, \\ldots, Nb2,
```
where ``\\{\\varphi_a^e\\}`` are the local basis functions of `fe` on a single
element ``\\Omega^e`` of the Cartesian mesh `mesh`.
 
Since all elements are identical on a uniform Cartesian mesh, this local
vector is the same for every element and needs to be computed only once.
 
The full vector is evaluated as
```math
\\mathbf{b}^e = \\frac{\\Delta x \\, \\Delta y}{4}
\\sum_{i,j} W_i W_j \\, \\boldsymbol{\\varphi}(P_i, P_j),
```
where ``\\boldsymbol{\\varphi}(P_i, P_j) \\in \\mathbb{R}^{Nb2}`` contains all local
basis functions evaluated at ``(P_i, P_j)``. The quadrature uses
`Npg = polynomial_degree(fe) + 1` Gauss–Legendre points per direction.
 
Returns an `SVector{Nb2, T}` where `Nb2 = num_local_dof(fe)`.
 
# Arguments
- `fe::DimensionalFEFamily`: 2D finite element family defining the local basis.
- `mesh::CartesianMesh{2}`: 2D uniform Cartesian mesh.
"""
function assembly_local_∫basis(
        fe::DimensionalFEFamily,
        mesh::CartesianMesh{2, I}
) where {I <: Integer}
    Δx, Δy = mesh.Δx
    Nb2 = num_local_dof(fe)
    Npg = polynomial_degree(fe) + 1
    jacobian = Δx * Δy / 4

    P, W = legendre(Npg)

    result = zeros(Nb2)
    for j in 1:Npg, i in 1:Npg
        φ = basis_functions(fe, P[i], P[j])
        result .= muladd(W[i] * W[j], φ, result)
    end

    return SVector{Nb2}(result) * jacobian
end