# ThermoelasticAcoustic 
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://bacarmo.github.io/ThermoelasticAcoustic.jl/stable/) 
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://bacarmo.github.io/ThermoelasticAcoustic.jl/dev/) 
[![Build Status](https://github.com/bacarmo/ThermoelasticAcoustic.jl/a-tions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/bacarmo/ThermoelasticAcoustic.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)

This repository provides all the code required to reproduce the numerical results presented in the article "Acoustic Boundary Conditions for a Nonlinear Coupled Thermoelastic System: Numerical Analysis".

## Installation
1. Install Julia v1.12.6
2. Get the repository, either by cloning it:
    ```bash
    git clone https://github.com/bacarmo/ThermoelasticAcoustic.jl.git
    ```
    or by downloading and extracting the [ZIP archive](https://github.com/bacarmo/ThermoelasticAcoustic.jl/archive/refs/heads/main.zip).

3. From the repository's root folder, instantiate the Julia environment
    ```bash
    julia --project=. -e 'using Pkg; Pkg.instantiate()'
    ```

## Reproducing the Results
Each Julia-REPL command below reproduces the numerical results reported in the corresponding table of the article.
### Spatial convergence order study (Table 2)
```julia-repl
julia> using ThermoelasticAcoustic

julia> cases1 = (
    (fe = Lagrange{1}, id = example1_manufactured(1.76), Nx = [2^i for i in 2:5], τ = 2.0^-13),
    (fe = Lagrange{1}, id = example1_manufactured(2.4 ), Nx = [2^i for i in 2:5], τ = 2.0^-13),
    (fe = Lagrange{2}, id = example1_manufactured(2.58), Nx = [2^i for i in 2:5], τ = 2.0^-13),
    (fe = Lagrange{2}, id = example1_manufactured(3.4 ), Nx = [2^i for i in 2:5], τ = 2.0^-13),
    (fe = Lagrange{3}, id = example1_manufactured(3.51), Nx = [2^i for i in 2:5], τ = 2.0^-13),
    (fe = Lagrange{3}, id = example1_manufactured(4.4 ), Nx = [2^i for i in 2:5], τ = 2.0^-13)
    );

julia> results1 = run_cases(Scheme1(), cases1)
```
### Temporal convergence order study (Table 3)
```julia-repl
julia> cases2 = (
    (fe = Lagrange{1}, id = example1_manufactured(2.4 ), Nx = 2^9, τ = [2.0^-i for i in 2:5]),
    );

julia> results2 = run_cases(Scheme1(), cases2)
```

## Declarations
### Acknowledgment 
This work was supported by CNPq - National Council for Scientific and Technological Development (Grant Number 151197/2025-3).

### AI use disclosure statement
In compliance with CNPq Portaria nº 2.664/2026 (Política de Integridade na Atividade Científica do CNPq), the following is disclosed.

During the implementation and documentation phase of the package, the author (B. A. Carmo) used Claude (Anthropic, free browser-based chat; exact model version not fixed) to: 
(i) review selected code sections to identify and fix bugs;
(ii) suggest more idiomatic or efficient Julia implementations;
(iii) refine docstrings and documentation text.
No IDE autocomplete or local AI agent was used.
All suggestions were reviewed by the author, who is fully responsible for the final code and results.

### Disclaimer
Everything is provided as is and without warranty. Use at your own risk.
