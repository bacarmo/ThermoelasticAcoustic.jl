# ==============================================================================
# Convergence studies
# ==============================================================================
"""
    ConvergenceResults{T,I}

Stores the output of a mesh-refinement convergence study.

Fields
- `info`  : description of the run
- `Nx`    : number of elements per spatial direction at each refinement level
- `h`     : element diameter `√(Δx²+Δy²)` at each level
- `τ`     : time-step size at each level
- `errors`: L∞(L²) error norms; size `(n_levels, n_fields)`
- `rates` : convergence rates log(eᵢ₋₁/eᵢ)/log(δᵢ₋₁/δᵢ) where δ = h (spatial/coupled study) or δ = τ (temporal study); same size as `errors`
"""
struct ConvergenceResults{T <: AbstractFloat, I <: Integer}
    info::String
    Nx::Vector{I}
    h::Vector{T}
    τ::Vector{T}
    errors::Matrix{T}  # (n_levels, n_fields)
    rates::Matrix{T}   # (n_levels, n_fields); first row is always zero
end

"""
    run_cases(solver = Scheme1(); Nx_values, τ_fixed, Nx_fixed, τ_values, t_end = 1.0) -> Vector{ConvergenceResults}

Runs the spatial convergence study (`Nx_values` with `τ_fixed`) and the temporal one (`τ_values` with `Nx_fixed`) 
writing each [`ConvergenceResults`](@ref) to a file named `convergence_studies_<solver>_<YYYY-MM-DD_HH-MM-SS>.txt` as it becomes available.
"""
function run_cases(
        solver::AbstractODESolver = Scheme1();
        Nx_values = [2^i for i in 2:5], τ_fixed = 2.0^-13,
        Nx_fixed = 2^9, τ_values = [2.0^-i for i in 2:5],
        t_end = 1.0
)
    cases_space = (
        (fe = Lagrange{1}, id = example1_manufactured(1.76)),
        (fe = Lagrange{1}, id = example1_manufactured(2.4)),
        (fe = Lagrange{2}, id = example1_manufactured(2.58)),
        (fe = Lagrange{2}, id = example1_manufactured(3.4)),
        (fe = Lagrange{3}, id = example1_manufactured(3.51)),
        (fe = Lagrange{3}, id = example1_manufactured(4.4)))

    cases_time = (
        (fe = Lagrange{1}, id = example1_manufactured(2.4)),
    )

    filename = "convergence_studies_$(typeof(solver))_" *
               Dates.format(now(), "yyyy-mm-dd_HH-MM-SS") * ".txt"

    open(filename, "a") do io
        versioninfo(io)
    end

    num_cases_space = length(cases_space)
    num_cases_time = length(cases_time)
    results = Vector{ConvergenceResults}(undef, num_cases_space+num_cases_time)

    for (i, case) in enumerate(cases_space)
        # Warmup: one lightweight solve to trigger JIT compilation before timing
        Base.with_logger(Base.NullLogger()) do
            run_convergence_study(case.fe, [2, 4], 1 / 4, solver, case.id, t_end)
        end

        results[i] = run_convergence_study(
            case.fe, Nx_values, τ_fixed, solver, case.id, t_end)

        write_result(filename, results[i])
    end

    for (i, case) in enumerate(cases_time)
        # Warmup: one lightweight solve to trigger JIT compilation before timing
        Base.with_logger(Base.NullLogger()) do
            run_convergence_study(case.fe, 4, [1/2, 1/4], solver, case.id, t_end)
        end

        results[i + num_cases_space] = run_convergence_study(
            case.fe, Nx_fixed, τ_values, solver, case.id, t_end)

        write_result(filename, results[i + num_cases_space])
    end

    return results
end

# ==============================================================================
"""
    write_result(filename, result)

Appends `result` to `filename` in plain-text format.
"""
function write_result(filename::String, result::ConvergenceResults)
    open(filename, "a") do io
        println(io, "="^100)
        show(io, result)
        println(io)
        flush(io)
    end
    return nothing
end

# ==============================================================================
"""
    run_convergence_study(fe, Nx_values, τ_values, solver, input_data, t_end) -> ConvergenceResults

Runs one PDE solve per entry in `Nx_values`/`τ_values`, collects L∞(L²) errors, and returns a [`ConvergenceResults`](@ref).

# Example
```julia-repl
julia> using Thermoelastic

julia> result = run_convergence_study(
           Lagrange{1}, [2, 4, 8, 16], sqrt(2) ./ [2, 4, 8, 16],
           Scheme1(), example1_manufactured(2.4), 1.0)
```
"""
function run_convergence_study(
        ::Type{fe},
        Nx_values::Vector{I},
        τ_values::Vector{T},
        solver::AbstractODESolver,
        input_data::PDEInputData,
        t_end::T
) where {T <: AbstractFloat, I <: Integer, fe <: AbstractFEBasis}
    pmin, pmax = input_data.pmin, input_data.pmax

    n_levels = length(Nx_values)
    n_fields = 5

    errors = zeros(T, n_levels, n_fields)
    rates = zeros(T, n_levels, n_fields)
    h_values = zeros(T, n_levels)
    Δx_values = zeros(T, n_levels)
    times = zeros(n_levels)

    for i in 1:n_levels
        Nx = Nx_values[i]
        Δx, Δy = (pmax .- pmin) ./ (Nx, Nx)
        Δx_values[i] = Δx
        h_values[i] = sqrt(Δx^2 + Δy^2)

        tspan = zero(T):τ_values[i]:t_end
        callback = L2ErrorCallback(tspan)
        times[i] = @elapsed solve_pde(fe, (Nx, Nx), tspan, input_data, solver, callback)
        errors[i, 1] = maximum(callback.v_errors)
        errors[i, 2] = maximum(callback.d_errors)
        errors[i, 3] = maximum(callback.c_errors)
        errors[i, 4] = maximum(callback.r_errors)
        errors[i, 5] = maximum(callback.z_errors)
    end

    δ_2D = allequal(Nx_values) ? τ_values : h_values
    δ_1D = allequal(Nx_values) ? τ_values : Δx_values
    compute_rates!(view(rates, :, 1:3), view(errors, :, 1:3), δ_2D)
    compute_rates!(view(rates, :, 4:5), view(errors, :, 4:5), δ_1D)

    total_time = sum(times)
    time_str = total_time ≥ 1.0 ? @sprintf("%.2f s", total_time) :
               @sprintf("%.1f ms", total_time*1e3)

    info = string(
        "ConvergenceResults: t_end=", t_end, ", ", input_data.name, ", ", fe, ", ", solver,
        ", elapsed=", time_str)
    return ConvergenceResults(info, Nx_values, h_values, τ_values, errors, rates)
end

# Convenience overloads: fix Nx, vary τ; or fix τ, vary Nx
function run_convergence_study(
        ::Type{fe}, Nx::I, τ_values::Vector{T},
        solver::AbstractODESolver, input_data::PDEInputData,
        t_end::T
) where {T <: AbstractFloat, I <: Integer, fe <: AbstractFEBasis}
    run_convergence_study(
        fe, fill(Nx, length(τ_values)), τ_values, solver, input_data, t_end)
end

function run_convergence_study(
        ::Type{fe}, Nx_values::Vector{I}, τ::T,
        solver::AbstractODESolver, input_data::PDEInputData,
        t_end::T
) where {T <: AbstractFloat, I <: Integer, fe <: AbstractFEBasis}
    run_convergence_study(
        fe, Nx_values, fill(τ, length(Nx_values)), solver, input_data, t_end)
end

# ==============================================================================
"""
    compute_rates!(rates, errors, δ)

Fills `rates` in-place: `rates[i,j] = log(e[i-1,j]/e[i,j]) / log(δ[i-1]/δ[i])` for `i ≥ 2`.
Row 1 is left as zero (no previous level available).
"""
function compute_rates!(rates::AbstractArray, errors::AbstractArray, δ::Vector)
    n_levels, n_fields = size(errors)
    for j in 1:n_fields, i in 2:n_levels

        rates[i, j] = log(errors[i - 1, j] / errors[i, j]) / log(δ[i - 1] / δ[i])
    end
    return nothing
end

# ==============================================================================
"""
    Base.show(io, r::ConvergenceResults)

Prints a compact table of mesh sizes, errors, and convergence rates.
Invoked automatically by `print`, `display`, and the REPL.
"""
function Base.show(io::IO, r::ConvergenceResults)
    n_levels, n_fields = size(r.errors)
    println(io, r.info)
    println(io,
        "  Nx   log₂h   log₂τ    L∞L²_V    rate     L∞L²_U    rate     L∞L²_Θ    rate     L∞L²_R    rate     L∞L²_Z    rate")
    for i in 1:n_levels
        row = @sprintf("%4d  %6.2f  %6.2f", r.Nx[i], log2(r.h[i]), log2(r.τ[i]))
        for j in 1:n_fields
            row *= @sprintf("    %.2e % .3f", r.errors[i, j], r.rates[i, j])
        end
        println(io, row)
    end
end