Documentation for [ThermoelasticAcoustic](https://github.com/bacarmo/ThermoelasticAcoustic.jl).

### Problem Definition
```@docs
PDEInputData
example1_manufactured
example1_zero_source
example2_manufactured
example2_zero_source
example3_manufactured
example3_zero_source
```

### Finite Element Family
```@docs
Lagrange
Hermite
```

### Time Integration
```@docs
CrankNicolson
ModifiedCN
```

### Solver
```@docs
pde_solve
```

### Callbacks
```@docs
L2ErrorCallback
```

### Postprocessing
```@docs
convergence_study_coupled
convergence_study_spatial
convergence_study_temporal
print_convergence_table
```

## Developer Reference
```@index
```
```@autodocs
Modules = [ThermoelasticAcoustic]
Filter = t -> !(t in [
    PDEInputData, example1_manufactured, example1_zero_source,
    example2_manufactured, example2_zero_source,
    example3_manufactured,
    example3_zero_source,
    Lagrange, Hermite,
    CrankNicolson, ModifiedCN,
    pde_solve,
    L2ErrorCallback,
    convergence_study_coupled, convergence_study_spatial,
    convergence_study_temporal, print_convergence_table
])
```