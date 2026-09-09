# Scheme 1
```julia-repl
julia> using ThermoelasticAcoustic

julia> cases1a = (
    (fe = Lagrange{1}, id = example1_manufactured(1.76), Nx = [2^i for i in 2:5], τ = 2.0^-13),
    (fe = Lagrange{1}, id = example1_manufactured(2.4 ), Nx = [2^i for i in 2:5], τ = 2.0^-13),
    (fe = Lagrange{2}, id = example1_manufactured(2.58), Nx = [2^i for i in 2:5], τ = 2.0^-13),
    (fe = Lagrange{2}, id = example1_manufactured(3.4 ), Nx = [2^i for i in 2:5], τ = 2.0^-13),
    (fe = Lagrange{3}, id = example1_manufactured(3.51), Nx = [2^i for i in 2:5], τ = 2.0^-13),
    (fe = Lagrange{3}, id = example1_manufactured(4.4 ), Nx = [2^i for i in 2:5], τ = 2.0^-13)
    );

julia> cases1b = (
    (fe = Lagrange{1}, id = example1_manufactured(1.76), Nx = [2^i for i in 2:6], τ = 2.0^-15),
    (fe = Lagrange{1}, id = example1_manufactured(2.4 ), Nx = [2^i for i in 2:6], τ = 2.0^-15),
    (fe = Lagrange{2}, id = example1_manufactured(2.58), Nx = [2^i for i in 2:6], τ = 2.0^-15),
    (fe = Lagrange{2}, id = example1_manufactured(3.4 ), Nx = [2^i for i in 2:6], τ = 2.0^-15),
    (fe = Lagrange{3}, id = example1_manufactured(3.51), Nx = [2^i for i in 2:6], τ = 2.0^-15),
    (fe = Lagrange{3}, id = example1_manufactured(4.4 ), Nx = [2^i for i in 2:6], τ = 2.0^-15)
    );

julia> results1a = run_cases(Scheme1(), cases1a)
6-element Vector{ConvergenceResults}:
ConvergenceResults: t_end=1.0, example1_manufactured(1.76), Lagrange{1}, Scheme1(), sum(walltime)=172.84 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -13.00    1.07e-02  0.00    3.08e-02  0.00    9.00e-03  0.00    6.68e-02  0.00    4.64e-02  0.00    1.7634e+00
   8   -2.50  -13.00    2.75e-03  1.97    7.90e-03  1.96    2.30e-03  1.97    1.51e-02  2.14    1.07e-02  2.12    6.4987e+00
  16   -3.50  -13.00    7.08e-04  1.96    2.04e-03  1.95    5.87e-04  1.97    3.56e-03  2.09    2.56e-03  2.06    2.7660e+01
  32   -4.50  -13.00    1.89e-04  1.91    5.38e-04  1.92    1.50e-04  1.97    8.56e-04  2.06    6.28e-04  2.03    1.3692e+02

ConvergenceResults: t_end=1.0, example1_manufactured(2.4), Lagrange{1}, Scheme1(), sum(walltime)=171.50 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -13.00    1.92e-02  0.00    4.51e-02  0.00    1.86e-02  0.00    6.68e-02  0.00    4.95e-02  0.00    1.5818e+00
   8   -2.50  -13.00    4.95e-03  1.96    1.11e-02  2.02    4.62e-03  2.01    1.50e-02  2.15    1.15e-02  2.11    6.4173e+00
  16   -3.50  -13.00    1.25e-03  1.99    2.77e-03  2.01    1.15e-03  2.00    3.61e-03  2.06    2.79e-03  2.04    2.7670e+01
  32   -4.50  -13.00    3.13e-04  2.00    6.92e-04  2.00    2.88e-04  2.00    8.88e-04  2.02    6.92e-04  2.01    1.3583e+02

ConvergenceResults: t_end=1.0, example1_manufactured(2.58), Lagrange{2}, Scheme1(), sum(walltime)=630.04 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -13.00    5.92e-04  0.00    1.45e-03  0.00    6.13e-04  0.00    4.14e-03  0.00    3.75e-03  0.00    3.0705e+00
   8   -2.50  -13.00    8.37e-05  2.82    1.93e-04  2.91    8.14e-05  2.91    6.15e-04  2.75    5.23e-04  2.84    1.3250e+01
  16   -3.50  -13.00    1.06e-05  2.98    2.53e-05  2.93    1.07e-05  2.93    8.37e-05  2.88    6.84e-05  2.93    7.3971e+01
  32   -4.50  -13.00    1.41e-06  2.91    3.30e-06  2.94    1.39e-06  2.94    1.10e-05  2.93    8.77e-06  2.96    5.3975e+02

ConvergenceResults: t_end=1.0, example1_manufactured(3.4), Lagrange{2}, Scheme1(), sum(walltime)=634.54 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -13.00    1.49e-03  0.00    3.34e-03  0.00    1.63e-03  0.00    5.06e-03  0.00    3.70e-03  0.00    2.9789e+00
   8   -2.50  -13.00    1.82e-04  3.03    4.18e-04  3.00    2.04e-04  3.00    6.40e-04  2.98    4.84e-04  2.93    1.4945e+01
  16   -3.50  -13.00    2.32e-05  2.97    5.23e-05  3.00    2.55e-05  3.00    8.32e-05  2.94    6.16e-05  2.98    8.3547e+01
  32   -4.50  -13.00    2.90e-06  3.00    6.53e-06  3.00    3.19e-06  3.00    1.06e-05  2.97    7.74e-06  2.99    5.3306e+02

ConvergenceResults: t_end=1.0, example1_manufactured(3.51), Lagrange{3}, Scheme1(), sum(walltime)=1615.00 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -13.00    3.15e-05  0.00    7.15e-05  0.00    3.54e-05  0.00    4.05e-04  0.00    1.37e-04  0.00    6.9103e+00
   8   -2.50  -13.00    2.07e-06  3.93    4.89e-06  3.87    2.42e-06  3.87    4.05e-05  3.32    1.13e-05  3.60    3.4375e+01
  16   -3.50  -13.00    1.37e-07  3.92    3.29e-07  3.89    1.63e-07  3.90    2.54e-06  3.99    7.14e-07  3.98    1.9626e+02
  32   -4.50  -13.00    9.26e-09  3.89    2.19e-08  3.91    1.08e-08  3.91    1.49e-07  4.09    4.31e-08  4.05    1.3775e+03

ConvergenceResults: t_end=1.0, example1_manufactured(4.4), Lagrange{3}, Scheme1(), sum(walltime)=1818.06 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -13.00    1.04e-04  0.00    2.21e-04  0.00    1.19e-04  0.00    5.88e-04  0.00    1.47e-04  0.00    6.7762e+00
   8   -2.50  -13.00    6.59e-06  3.97    1.40e-05  3.98    7.52e-06  3.99    1.01e-04  2.55    2.41e-05  2.61    3.6991e+01
  16   -3.50  -13.00    3.79e-07  4.12    8.76e-07  4.00    4.70e-07  4.00    7.05e-06  3.84    1.69e-06  3.84    2.3734e+02
  32   -4.50  -13.00    2.39e-08  3.99    5.48e-08  4.00    2.94e-08  4.00    4.16e-07  4.08    1.00e-07  4.08    1.5370e+03

julia> results1b = run_cases(Scheme1(), cases1b)
6-element Vector{ConvergenceResults}:
ConvergenceResults: t_end=1.0, example1_manufactured(1.76), Lagrange{1}, Scheme1(), sum(walltime)=3172.05 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -15.00    1.07e-02  0.00    3.08e-02  0.00    9.00e-03  0.00    6.68e-02  0.00    4.64e-02  0.00    6.3320e+00
   8   -2.50  -15.00    2.75e-03  1.97    7.90e-03  1.96    2.30e-03  1.97    1.51e-02  2.14    1.07e-02  2.12    2.5672e+01
  16   -3.50  -15.00    7.08e-04  1.96    2.04e-03  1.95    5.87e-04  1.97    3.56e-03  2.09    2.56e-03  2.06    1.1013e+02
  32   -4.50  -15.00    1.89e-04  1.91    5.38e-04  1.92    1.50e-04  1.97    8.56e-04  2.06    6.28e-04  2.03    5.2864e+02
  64   -5.50  -15.00    6.23e-05  1.60    1.50e-04  1.85    3.87e-05  1.95    2.12e-04  2.01    1.57e-04  2.00    2.5013e+03

ConvergenceResults: t_end=1.0, example1_manufactured(2.4), Lagrange{1}, Scheme1(), sum(walltime)=3583.06 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -15.00    1.92e-02  0.00    4.51e-02  0.00    1.86e-02  0.00    6.68e-02  0.00    4.95e-02  0.00    6.3349e+00
   8   -2.50  -15.00    4.95e-03  1.96    1.11e-02  2.02    4.62e-03  2.01    1.50e-02  2.15    1.15e-02  2.11    2.5713e+01
  16   -3.50  -15.00    1.25e-03  1.99    2.77e-03  2.01    1.15e-03  2.00    3.61e-03  2.06    2.79e-03  2.04    1.1110e+02
  32   -4.50  -15.00    3.13e-04  2.00    6.92e-04  2.00    2.88e-04  2.00    8.88e-04  2.02    6.92e-04  2.01    5.4323e+02
  64   -5.50  -15.00    7.85e-05  2.00    1.73e-04  2.00    7.19e-05  2.00    2.21e-04  2.01    1.73e-04  2.00    2.8967e+03

ConvergenceResults: t_end=1.0, example1_manufactured(2.58), Lagrange{2}, Scheme1(), sum(walltime)=18840.08 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -15.00    5.92e-04  0.00    1.45e-03  0.00    6.13e-04  0.00    4.14e-03  0.00    3.75e-03  0.00    1.1022e+01
   8   -2.50  -15.00    8.37e-05  2.82    1.93e-04  2.91    8.14e-05  2.91    6.15e-04  2.75    5.23e-04  2.84    5.3380e+01
  16   -3.50  -15.00    1.06e-05  2.98    2.53e-05  2.93    1.07e-05  2.93    8.37e-05  2.88    6.84e-05  2.93    2.9471e+02
  32   -4.50  -15.00    1.41e-06  2.91    3.30e-06  2.94    1.39e-06  2.94    1.10e-05  2.93    8.77e-06  2.96    2.1062e+03
  64   -5.50  -15.00    1.96e-07  2.85    4.31e-07  2.94    1.79e-07  2.95    1.42e-06  2.95    1.12e-06  2.97    1.6375e+04

ConvergenceResults: t_end=1.0, example1_manufactured(3.4), Lagrange{2}, Scheme1(), sum(walltime)=19401.58 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -15.00    1.49e-03  0.00    3.34e-03  0.00    1.63e-03  0.00    5.06e-03  0.00    3.70e-03  0.00    1.1065e+01
   8   -2.50  -15.00    1.82e-04  3.03    4.18e-04  3.00    2.04e-04  3.00    6.40e-04  2.98    4.84e-04  2.93    5.3089e+01
  16   -3.50  -15.00    2.32e-05  2.97    5.23e-05  3.00    2.55e-05  3.00    8.32e-05  2.94    6.16e-05  2.98    2.9423e+02
  32   -4.50  -15.00    2.90e-06  3.00    6.53e-06  3.00    3.19e-06  3.00    1.06e-05  2.97    7.74e-06  2.99    2.1278e+03
  64   -5.50  -15.00    3.63e-07  3.00    8.17e-07  3.00    3.99e-07  3.00    1.34e-06  2.99    9.68e-07  3.00    1.6915e+04

ConvergenceResults: t_end=1.0, example1_manufactured(3.51), Lagrange{3}, Scheme1(), sum(walltime)=50124.23 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -15.00    3.15e-05  0.00    7.15e-05  0.00    3.54e-05  0.00    4.05e-04  0.00    1.37e-04  0.00    2.2086e+01
   8   -2.50  -15.00    2.07e-06  3.93    4.89e-06  3.87    2.42e-06  3.87    4.05e-05  3.32    1.13e-05  3.60    1.1759e+02
  16   -3.50  -15.00    1.37e-07  3.92    3.29e-07  3.89    1.63e-07  3.90    2.54e-06  3.99    7.14e-07  3.98    7.5159e+02
  32   -4.50  -15.00    9.21e-09  3.89    2.19e-08  3.91    1.08e-08  3.91    1.49e-07  4.09    4.30e-08  4.05    5.4924e+03
  64   -5.50  -15.00    9.39e-10  3.29    1.45e-09  3.92    4.21e-09  1.36    8.81e-09  4.08    2.62e-09  4.04    4.3741e+04

ConvergenceResults: t_end=1.0, example1_manufactured(4.4), Lagrange{3}, Scheme1(), sum(walltime)=50323.17 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -15.00    1.04e-04  0.00    2.21e-04  0.00    1.19e-04  0.00    5.88e-04  0.00    1.47e-04  0.00    2.3416e+01
   8   -2.50  -15.00    6.59e-06  3.97    1.40e-05  3.98    7.52e-06  3.99    1.01e-04  2.55    2.41e-05  2.61    1.1675e+02
  16   -3.50  -15.00    3.79e-07  4.12    8.76e-07  4.00    4.70e-07  4.00    7.05e-06  3.84    1.69e-06  3.84    7.5529e+02
  32   -4.50  -15.00    2.38e-08  4.00    5.48e-08  4.00    2.94e-08  4.00    4.16e-07  4.08    1.00e-07  4.08    5.5210e+03
  64   -5.50  -15.00    1.56e-09  3.93    3.43e-09  4.00    4.34e-09  2.76    2.43e-08  4.10    5.84e-09  4.10    4.3907e+04

julia> cases2a = (
    (fe = Lagrange{1}, id = example2_manufactured(1.76), Nx = [2^i for i in 2:5], τ = 2.0^-13),
    (fe = Lagrange{1}, id = example2_manufactured(2.4 ), Nx = [2^i for i in 2:5], τ = 2.0^-13),
    (fe = Lagrange{2}, id = example2_manufactured(2.58), Nx = [2^i for i in 2:5], τ = 2.0^-13),
    (fe = Lagrange{2}, id = example2_manufactured(3.4 ), Nx = [2^i for i in 2:5], τ = 2.0^-13),
    (fe = Lagrange{3}, id = example2_manufactured(3.51), Nx = [2^i for i in 2:5], τ = 2.0^-13),
    (fe = Lagrange{3}, id = example2_manufactured(4.4 ), Nx = [2^i for i in 2:5], τ = 2.0^-13)
    );

julia> cases2b = (
    (fe = Lagrange{1}, id = example2_manufactured(1.76), Nx = [2^i for i in 2:6], τ = 2.0^-15),
    (fe = Lagrange{1}, id = example2_manufactured(2.4 ), Nx = [2^i for i in 2:6], τ = 2.0^-15),
    (fe = Lagrange{2}, id = example2_manufactured(2.58), Nx = [2^i for i in 2:6], τ = 2.0^-15),
    (fe = Lagrange{2}, id = example2_manufactured(3.4 ), Nx = [2^i for i in 2:6], τ = 2.0^-15),
    (fe = Lagrange{3}, id = example2_manufactured(3.51), Nx = [2^i for i in 2:6], τ = 2.0^-15),
    (fe = Lagrange{3}, id = example2_manufactured(4.4 ), Nx = [2^i for i in 2:6], τ = 2.0^-15)
    );

julia> results2a = run_cases(Scheme1(), cases2a)
6-element Vector{ConvergenceResults}:
ConvergenceResults: t_end=1.0, example2_manufactured(1.76), Lagrange{1}, Scheme1(), sum(walltime)=171.25 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -13.00    1.11e-02  0.00    3.08e-02  0.00    9.01e-03  0.00    2.21e-02  0.00    2.60e-02  0.00    1.5947e+00
   8   -2.50  -13.00    2.80e-03  1.99    7.96e-03  1.95    2.31e-03  1.97    5.05e-03  2.13    6.06e-03  2.10    6.4138e+00
  16   -3.50  -13.00    7.08e-04  1.98    2.09e-03  1.93    5.90e-04  1.97    1.21e-03  2.07    1.48e-03  2.03    2.7685e+01
  32   -4.50  -13.00    2.11e-04  1.75    5.69e-04  1.87    1.52e-04  1.96    3.10e-04  1.96    3.78e-04  1.97    1.3556e+02

ConvergenceResults: t_end=1.0, example2_manufactured(2.4), Lagrange{1}, Scheme1(), sum(walltime)=172.63 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -13.00    1.90e-02  0.00    4.45e-02  0.00    1.87e-02  0.00    2.21e-02  0.00    2.80e-02  0.00    1.5773e+00
   8   -2.50  -13.00    4.86e-03  1.96    1.10e-02  2.02    4.63e-03  2.01    5.09e-03  2.12    6.58e-03  2.09    6.4313e+00
  16   -3.50  -13.00    1.23e-03  1.99    2.73e-03  2.01    1.16e-03  2.00    1.23e-03  2.04    1.62e-03  2.03    2.7702e+01
  32   -4.50  -13.00    3.08e-04  1.99    6.82e-04  2.00    2.89e-04  2.00    3.05e-04  2.02    4.02e-04  2.01    1.3692e+02

ConvergenceResults: t_end=1.0, example2_manufactured(2.58), Lagrange{2}, Scheme1(), sum(walltime)=620.59 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -13.00    5.88e-04  0.00    1.45e-03  0.00    6.13e-04  0.00    1.41e-03  0.00    2.44e-03  0.00    2.8525e+00
   8   -2.50  -13.00    8.33e-05  2.82    1.93e-04  2.91    8.14e-05  2.91    2.03e-04  2.80    3.31e-04  2.88    1.3270e+01
  16   -3.50  -13.00    1.06e-05  2.97    2.53e-05  2.93    1.07e-05  2.93    2.73e-05  2.89    4.26e-05  2.96    7.3285e+01
  32   -4.50  -13.00    1.48e-06  2.84    3.31e-06  2.94    1.39e-06  2.94    3.59e-06  2.93    5.41e-06  2.98    5.3118e+02

ConvergenceResults: t_end=1.0, example2_manufactured(3.4), Lagrange{2}, Scheme1(), sum(walltime)=632.29 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -13.00    1.47e-03  0.00    3.34e-03  0.00    1.63e-03  0.00    1.66e-03  0.00    2.35e-03  0.00    2.9788e+00
   8   -2.50  -13.00    1.81e-04  3.03    4.18e-04  3.00    2.04e-04  3.00    2.10e-04  2.98    3.12e-04  2.92    1.5037e+01
  16   -3.50  -13.00    2.31e-05  2.96    5.23e-05  3.00    2.55e-05  3.00    2.64e-05  2.99    3.96e-05  2.98    8.2812e+01
  32   -4.50  -13.00    2.89e-06  3.00    6.53e-06  3.00    3.19e-06  3.00    3.31e-06  3.00    4.96e-06  2.99    5.3146e+02

ConvergenceResults: t_end=1.0, example2_manufactured(3.51), Lagrange{3}, Scheme1(), sum(walltime)=1625.79 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -13.00    3.12e-05  0.00    7.15e-05  0.00    3.54e-05  0.00    4.64e-05  0.00    6.93e-05  0.00    6.6219e+00
   8   -2.50  -13.00    2.06e-06  3.92    4.90e-06  3.87    2.42e-06  3.87    3.12e-06  3.90    4.29e-06  4.02    3.5333e+01
  16   -3.50  -13.00    1.36e-07  3.91    3.29e-07  3.89    1.63e-07  3.90    2.05e-07  3.93    2.69e-07  4.00    1.9867e+02
  32   -4.50  -13.00    9.20e-09  3.89    2.19e-08  3.91    1.08e-08  3.91    1.34e-08  3.93    1.69e-08  3.99    1.3852e+03

ConvergenceResults: t_end=1.0, example2_manufactured(4.4), Lagrange{3}, Scheme1(), sum(walltime)=1831.76 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -13.00    1.02e-04  0.00    2.21e-04  0.00    1.19e-04  0.00    4.92e-05  0.00    5.59e-05  0.00    7.0161e+00
   8   -2.50  -13.00    6.55e-06  3.96    1.40e-05  3.98    7.52e-06  3.99    2.94e-06  4.06    3.38e-06  4.05    3.7481e+01
  16   -3.50  -13.00    3.81e-07  4.11    8.76e-07  4.00    4.70e-07  4.00    1.78e-07  4.05    2.09e-07  4.01    2.4142e+02
  32   -4.50  -13.00    2.39e-08  4.00    5.48e-08  4.00    2.94e-08  4.00    1.09e-08  4.03    1.30e-08  4.00    1.5458e+03

julia> results2b = run_cases(Scheme1(), cases2b)
6-element Vector{ConvergenceResults}:
ConvergenceResults: t_end=1.0, example2_manufactured(1.76), Lagrange{1}, Scheme1(), sum(walltime)=3177.63 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -15.00    1.11e-02  0.00    3.08e-02  0.00    9.01e-03  0.00    2.21e-02  0.00    2.60e-02  0.00    6.3038e+00
   8   -2.50  -15.00    2.80e-03  1.99    7.96e-03  1.95    2.31e-03  1.97    5.05e-03  2.13    6.06e-03  2.10    2.5639e+01
  16   -3.50  -15.00    7.08e-04  1.98    2.09e-03  1.93    5.90e-04  1.97    1.21e-03  2.07    1.48e-03  2.03    1.1006e+02
  32   -4.50  -15.00    2.11e-04  1.75    5.69e-04  1.87    1.52e-04  1.96    3.10e-04  1.96    3.78e-04  1.97    5.2887e+02
  64   -5.50  -15.00    1.17e-04  0.85    1.71e-04  1.73    3.97e-05  1.93    9.99e-05  1.63    1.09e-04  1.80    2.5068e+03

ConvergenceResults: t_end=1.0, example2_manufactured(2.4), Lagrange{1}, Scheme1(), sum(walltime)=3621.61 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -15.00    1.90e-02  0.00    4.45e-02  0.00    1.87e-02  0.00    2.21e-02  0.00    2.80e-02  0.00    6.3088e+00
   8   -2.50  -15.00    4.86e-03  1.96    1.10e-02  2.02    4.63e-03  2.01    5.09e-03  2.12    6.58e-03  2.09    2.5710e+01
  16   -3.50  -15.00    1.23e-03  1.99    2.73e-03  2.01    1.16e-03  2.00    1.23e-03  2.04    1.62e-03  2.03    1.1099e+02
  32   -4.50  -15.00    3.08e-04  1.99    6.82e-04  2.00    2.89e-04  2.00    3.05e-04  2.02    4.02e-04  2.01    5.4233e+02
  64   -5.50  -15.00    7.72e-05  2.00    1.70e-04  2.00    7.21e-05  2.00    7.59e-05  2.01    1.00e-04  2.00    2.9363e+03

ConvergenceResults: t_end=1.0, example2_manufactured(2.58), Lagrange{2}, Scheme1(), sum(walltime)=18770.62 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -15.00    5.88e-04  0.00    1.45e-03  0.00    6.13e-04  0.00    1.41e-03  0.00    2.44e-03  0.00    1.0974e+01
   8   -2.50  -15.00    8.33e-05  2.82    1.93e-04  2.91    8.14e-05  2.91    2.03e-04  2.80    3.31e-04  2.88    5.3446e+01
  16   -3.50  -15.00    1.06e-05  2.97    2.53e-05  2.93    1.07e-05  2.93    2.73e-05  2.89    4.26e-05  2.96    2.9397e+02
  32   -4.50  -15.00    1.48e-06  2.84    3.31e-06  2.94    1.39e-06  2.94    3.59e-06  2.93    5.41e-06  2.98    2.1122e+03
  64   -5.50  -15.00    2.75e-07  2.43    4.40e-07  2.91    1.79e-07  2.95    4.81e-07  2.90    6.92e-07  2.97    1.6300e+04

ConvergenceResults: t_end=1.0, example2_manufactured(3.4), Lagrange{2}, Scheme1(), sum(walltime)=19404.17 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -15.00    1.47e-03  0.00    3.34e-03  0.00    1.63e-03  0.00    1.66e-03  0.00    2.35e-03  0.00    1.1013e+01
   8   -2.50  -15.00    1.80e-04  3.03    4.18e-04  3.00    2.04e-04  3.00    2.10e-04  2.98    3.12e-04  2.92    5.3362e+01
  16   -3.50  -15.00    2.31e-05  2.96    5.23e-05  3.00    2.55e-05  3.00    2.64e-05  2.99    3.96e-05  2.98    2.9575e+02
  32   -4.50  -15.00    2.89e-06  3.00    6.53e-06  3.00    3.19e-06  3.00    3.31e-06  3.00    4.96e-06  2.99    2.1247e+03
  64   -5.50  -15.00    3.62e-07  3.00    8.17e-07  3.00    3.99e-07  3.00    4.13e-07  3.00    6.21e-07  3.00    1.6919e+04

ConvergenceResults: t_end=1.0, example2_manufactured(3.51), Lagrange{3}, Scheme1(), sum(walltime)=49570.72 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -15.00    3.12e-05  0.00    7.15e-05  0.00    3.54e-05  0.00    4.64e-05  0.00    6.93e-05  0.00    2.2230e+01
   8   -2.50  -15.00    2.06e-06  3.92    4.90e-06  3.87    2.42e-06  3.87    3.12e-06  3.90    4.29e-06  4.02    1.1737e+02
  16   -3.50  -15.00    1.36e-07  3.91    3.29e-07  3.89    1.63e-07  3.90    2.05e-07  3.93    2.69e-07  4.00    7.3353e+02
  32   -4.50  -15.00    9.19e-09  3.89    2.19e-08  3.91    1.08e-08  3.91    1.34e-08  3.93    1.69e-08  3.99    5.4702e+03
  64   -5.50  -15.00    9.87e-10  3.22    1.45e-09  3.92    4.21e-09  1.36    8.74e-10  3.94    1.07e-09  3.99    4.3227e+04

ConvergenceResults: t_end=1.0, example2_manufactured(4.4), Lagrange{3}, Scheme1(), sum(walltime)=50392.89 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -15.00    1.02e-04  0.00    2.21e-04  0.00    1.19e-04  0.00    4.92e-05  0.00    5.59e-05  0.00    2.4659e+01
   8   -2.50  -15.00    6.55e-06  3.96    1.40e-05  3.98    7.52e-06  3.99    2.94e-06  4.06    3.38e-06  4.05    1.2392e+02
  16   -3.50  -15.00    3.81e-07  4.11    8.76e-07  4.00    4.70e-07  4.00    1.78e-07  4.05    2.09e-07  4.01    7.7437e+02
  32   -4.50  -15.00    2.38e-08  4.00    5.48e-08  4.00    2.94e-08  4.00    1.09e-08  4.03    1.30e-08  4.00    5.5611e+03
  64   -5.50  -15.00    1.55e-09  3.94    3.43e-09  4.00    4.34e-09  2.76    6.74e-10  4.01    8.15e-10  4.00    4.3909e+04

julia> cases3 = (
    (fe = Lagrange{1}, id = example1_manufactured(2.4 ), Nx = 2^9, τ = [2.0^-i for i in 2:5]),
    (fe = Lagrange{1}, id = example2_manufactured(2.4 ), Nx = 2^9, τ = [2.0^-i for i in 2:5])
    );

julia> results3 = run_cases(Scheme1(), cases3)
2-element Vector{ConvergenceResults}:
ConvergenceResults: t_end=1.0, example1_manufactured(2.4), Lagrange{1}, Scheme1(), sum(walltime)=5157.03 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
 512   -8.50   -2.00    1.09e-02  0.00    3.47e-03  0.00    1.45e-03  0.00    7.12e-03  0.00    5.41e-03  0.00    3.9053e+02
 512   -8.50   -3.00    2.72e-03  2.00    8.95e-04  1.96    3.00e-04  2.27    1.79e-03  1.99    1.36e-03  1.99    7.3023e+02
 512   -8.50   -4.00    6.76e-04  2.01    2.27e-04  1.98    6.17e-05  2.28    4.50e-04  2.00    3.41e-04  2.00    1.4229e+03
 512   -8.50   -5.00    1.68e-04  2.01    5.80e-05  1.97    1.60e-05  1.95    1.13e-04  2.00    8.54e-05  2.00    2.6134e+03

ConvergenceResults: t_end=1.0, example2_manufactured(2.4), Lagrange{1}, Scheme1(), sum(walltime)=5168.19 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
 512   -8.50   -2.00    1.12e-02  0.00    3.90e-03  0.00    1.47e-03  0.00    1.44e-03  0.00    6.56e-04  0.00    3.8859e+02
 512   -8.50   -3.00    2.82e-03  2.00    9.80e-04  1.99    3.03e-04  2.28    4.24e-04  1.76    1.62e-04  2.02    7.2944e+02
 512   -8.50   -4.00    7.01e-04  2.01    2.47e-04  1.99    6.47e-05  2.23    1.12e-04  1.93    4.07e-05  1.99    1.4272e+03
 512   -8.50   -5.00    1.74e-04  2.01    6.26e-05  1.98    1.67e-05  1.96    2.87e-05  1.96    1.05e-05  1.96    2.6230e+03

julia> versioninfo()
Julia Version 1.12.7
Commit 6d172b025e4 (2026-08-15 08:05 UTC)
Build Info:
  Official https://julialang.org release
Platform Info:
  OS: Linux (x86_64-linux-gnu)
  CPU: 14 × Intel(R) Core(TM) Ultra 5 225H
  WORD_SIZE: 64
  LLVM: libLLVM-18.1.7 (ORCJIT, arrowlake)
  GC: Built with stock GC
Threads: 1 default, 1 interactive, 1 GC (on 14 virtual cores)
```