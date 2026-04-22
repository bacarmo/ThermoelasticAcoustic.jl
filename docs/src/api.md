## Index
```@index
```

## Public API

### Problem Definition
```@docs
PDEInputData
example0_manufactured
example0_zero_source
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
SolutionCallback
```

### Postprocessing
```@docs
convergence_study_coupled
convergence_study_spatial
convergence_study_temporal
print_convergence_table
```

## Internal API
```@autodocs
Modules = [ThermoelasticAcoustic]
Filter = t -> !(t in [
    PDEInputData, example0_manufactured, example0_zero_source,
    example1_manufactured, example1_zero_source,
    example2_manufactured, example2_zero_source,
    example3_manufactured,
    example3_zero_source,
    Lagrange, Hermite,
    CrankNicolson, ModifiedCN,
    pde_solve,
    L2ErrorCallback, SolutionCallback,
    convergence_study_coupled, convergence_study_spatial,
    convergence_study_temporal, print_convergence_table
])
```