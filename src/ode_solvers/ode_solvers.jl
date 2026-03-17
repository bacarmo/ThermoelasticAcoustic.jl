# ========================================
# Type definitions
# ========================================
"""
    ODESolver

Abstract type for time integration schemes used in [`pde_solve`](@ref).

# Concrete subtypes
- [`CrankNicolson`](@ref)
- [`ModifiedCN`](@ref)
"""
abstract type ODESolver end

"""
    CrankNicolson <: ODESolver

Standard Crank-Nicolson scheme for the coupled thermo-wave-acoustic system.

See *Scheme 1, Strategy 2* in the method documentation for the full problem formulation.

# Example
```julia
pde_solve(Nx, fe, tspan, input_data, CrankNicolson(), callback)
```
"""
struct CrankNicolson <: ODESolver end

"""
    ModifiedCN <: ODESolver

Modified Crank-Nicolson scheme for the coupled thermo-wave-acoustic system.

See *Scheme 2* in the method documentation for the full problem formulation.

# Example
```julia
pde_solve(Nx, fe, tspan, input_data, ModifiedCN(), callback)
```
"""
struct ModifiedCN <: ODESolver end

# ========================================
# Generic interfaces
# ========================================
"""
    build_cache(solver, matrices, dof_map_m₁, dof_map_m₂, dof_map_m₃, τ)

Allocate and return the solver-specific cache for the time integration.
Dispatches on `solver` to construct the appropriate cache type.

This is the only allocation point for the time integration loop; all
subsequent calls to `ode_solve` are allocation-free.

# Arguments
- `solver::ODESolver`: time integration scheme (e.g. `ModifiedCN()`).
- `matrices::SystemMatrices{T,I}`: global FEM matrices.
- `dof_map_m₁::DOFMap`: DOF map for the `m₁`-space (wave field).
- `dof_map_m₂::DOFMap`: DOF map for the `m₂`-space (heat field).
- `dof_map_m₃::DOFMap`: DOF map for the `m₃`-space (acoustic field).
- `τ::T`: time-step size; used by solver-specific constructors that require it
  to pre-fill time-invariant matrix blocks (e.g. [`CrankNicolson`](@ref)).
 

# Extended help
To add support for a new solver `MySolver`, define:
```julia
build_cache(::MySolver, matrices, dof_map_m₁, dof_map_m₂, dof_map_m₃, τ) =
    MySolverCache(matrices, dof_map_m₁, dof_map_m₂, dof_map_m₃, τ)
```
"""
function build_cache end

"""
    ode_solve(solver, cache, state, matrices, dof_map_m₁, dof_map_m₂, dof_map_m₃,
              mesh1D, mesh2D, quad, tspan, input_data, callback)

Advance `state` over `tspan` using the time integration scheme `solver`.
Dispatches on `solver` to the appropriate implementation.

No allocation occurs; all workspace is provided via `cache`.

# Arguments
- `solver::ODESolver`: time integration scheme.
- `cache`: pre-allocated workspace; see [`build_cache`](@ref).
- `state::FEMState`: solution state; modified in-place at each time step.
- `matrices::SystemMatrices`: global FEM matrices.
- `dof_map_m₁::DOFMap`: DOF map for the `m₁`-space (wave field).
- `dof_map_m₂::DOFMap`: DOF map for the `m₂`-space (heat field).
- `dof_map_m₃::DOFMap`: DOF map for the `m₃`-space (acoustic field).
- `mesh1D::CartesianMesh{1}`: 1D Cartesian mesh.
- `mesh2D::CartesianMesh{2}`: 2D Cartesian mesh.
- `quad::QuadratureSetup`: precomputed quadrature data.
- `tspan::StepRangeLen{T}`: time grid defined as `t₀:τ:t_end`.
- `input_data::PDEInputData`: problem data (source terms, coefficients, boundary data).
- `callback::AbstractCallback`: invoked after each accepted time step.
"""
function ode_solve end

# ========================================
# Implementations
# ========================================
include("cache_cn.jl")
include("crank_nicolson.jl")

include("cache_mcn.jl")
include("modified_crank_nicolson.jl")