# ==============================================================================
# run_convergence_studies.jl
#
# Runs coupled and spatial convergence studies for the ThermoelasticAcoustic
# package and saves the results to a timestamped plain-text file in the package root.
#
# HOW TO RUN
# ----------
# 0. (First use only) Download the package and install all dependencies.
#
#    Option A — Using Git:
#
#       git clone https://github.com/bacarmo/ThermoelasticAcoustic.jl.git
#
#    Option B — Without Git:
#       Go to https://github.com/bacarmo/ThermoelasticAcoustic.jl,
#       click "Code" > "Download ZIP", and extract the folder.
#
#    Then, open a terminal at the package root and run:
#
#       julia
#       ] activate .
#       ] instantiate
#       <backspace>  (to exit Pkg mode)
#
#    This step is only required once per machine.
#
# 1. Open a terminal at the package root:
#
#       ThermoelasticAcoustic/   <-- open the terminal here
#       ├── Project.toml
#       └── src/
#           └── postprocessing/
#               └── run_convergence_studies.jl
#
# 2. Start Julia and activate the package environment:
#
#       julia
#       ] activate .
#       <backspace>  (to exit Pkg mode)
#
# 3. Load and run the script:
#
#       include("src/postprocessing/run_convergence_studies.jl")
#
#    To re-run after modifying CASES or the package source, just call:
#
#       include("src/postprocessing/run_convergence_studies.jl")
#
# 4. Results are saved to a timestamped .txt file in the package root, e.g.:
#
#       convergence_results_YYYY-MM-DD__HH-MM-SS.txt
# ==============================================================================

using Revise
using Dates
using Printf
using ThermoelasticAcoustic

# ── Configuration ──────────────────────────────────────────────────────────────

# Cases to run
#! format: off
const CASES = (
    (input_data = example1_manufactured(1.76), solver = ModifiedCN(),    fe = Lagrange{1}(), Nx_exp_range = 3:6, τ_fixed = 2.0^(-15)),
    (input_data = example1_manufactured(2.4 ), solver = ModifiedCN(),    fe = Lagrange{1}(), Nx_exp_range = 3:6, τ_fixed = 2.0^(-15)),
    (input_data = example1_manufactured(1.76), solver = CrankNicolson(), fe = Lagrange{1}(), Nx_exp_range = 3:6, τ_fixed = 2.0^(-15)),
    (input_data = example1_manufactured(2.4 ), solver = CrankNicolson(), fe = Lagrange{1}(), Nx_exp_range = 3:6, τ_fixed = 2.0^(-15)),
)
#! format: on

# Studies to run - set any entry to false to skip it
const RUN_COUPLED = true
const RUN_SPATIAL = true

# Output file - timestamp prevents overwriting previous runs
const OUTPUT_FILE = "convergence_results_$(Dates.format(now(), "yyyy-mm-dd__HH-MM-SS")).txt"

# ── Helpers ────────────────────────────────────────────────────────────────────

# Shared loop: warm-up, run, print, flush, error handling
function run_study!(io, label, cases, warmup_fn, study_fn)
    println("# $label convergence study")
    println("# Run started: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))\n")
    flush(io)

    for case in cases
        try
            warmup_fn(case)
            elapsed = @elapsed results = study_fn(case)
            print_convergence_table(results)
            @printf("  Elapsed: %.2f s\n", elapsed)
        catch e
            println("\n[ERROR] $(case.input_data.name) · $(typeof(case.solver)) · " *
                    "$(typeof(case.fe)):\n  $e\n")
        end
        flush(io)
    end

    println("\n# Run finished: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
end

# ── Study definitions ──────────────────────────────────────────────────────────

# Coupled study — refines Nx and τ simultaneously
function warmup_coupled(case)
    #! format: off
    convergence_study_coupled(
        input_data   = case.input_data,
        solver       = case.solver,
        fe           = case.fe,
        Nx_exp_range = 2:2)
    #! format: on
end

function study_coupled(case)
    #! format: off
    convergence_study_coupled(
        input_data   = case.input_data,
        solver       = case.solver,
        fe           = case.fe,
        Nx_exp_range = case.Nx_exp_range)
    #! format: on
end

# Spatial study — refines Nx with τ fixed
function warmup_spatial(case)
    #! format: off
    convergence_study_spatial(
        input_data   = case.input_data,
        solver       = case.solver,
        fe           = case.fe,
        Nx_exp_range = 2:2,
        τ_fixed      = 2.0^(-3))
    #! format: on
end

function study_spatial(case)
    #! format: off
    convergence_study_spatial(
        input_data   = case.input_data,
        solver       = case.solver,
        fe           = case.fe,
        Nx_exp_range = case.Nx_exp_range,
        τ_fixed      = case.τ_fixed)
    #! format: on
end

# ── Main ───────────────────────────────────────────────────────────────────────

# open: creates the file and closes it on exit, even if an exception is raised
open(OUTPUT_FILE, "w") do io
    # redirect_stdout: redirects print/println/@printf to `io` instead of the terminal
    redirect_stdout(io) do
        RUN_COUPLED && run_study!(io, "Coupled", CASES, warmup_coupled, study_coupled)
        RUN_SPATIAL && run_study!(io, "Spatial", CASES, warmup_spatial, study_spatial)
    end
end