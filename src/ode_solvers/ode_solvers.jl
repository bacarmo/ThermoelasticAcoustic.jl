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
# ODE solvers
# ========================================
include("crank_nicolson.jl")
include("modified_crank_nicolson.jl")