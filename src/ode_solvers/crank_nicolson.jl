function ode_solve(
        ::CrankNicolson,
        state::FEMState{T},
        matrices::SystemMatrices,
        dof_map_m₁::DOFMap,
        dof_map_m₂::DOFMap,
        dof_map_m₃::DOFMap,
        mesh1D::CartesianMesh{1, I},
        mesh2D::CartesianMesh{2, I},
        quad::QuadratureSetup,
        tspan::StepRangeLen{T},
        input_data::PDEInputData,
        callback::AbstractCallback
) where {T <: AbstractFloat, I <: Integer}
    return callback
end