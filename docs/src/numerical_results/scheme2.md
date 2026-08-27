# Scheme 2
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

julia> results1a = run_cases(Scheme2(), cases1a)
6-element Vector{ConvergenceResults}:
 ConvergenceResults: t_end=1.0, example1_manufactured(1.76), Lagrange{1}, Scheme2(), sum(walltime)=129.08 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -13.00    1.07e-02  0.00    3.08e-02  0.00    9.00e-03  0.00    6.68e-02  0.00    4.64e-02  0.00    1.4815e+00
   8   -2.50  -13.00    2.75e-03  1.97    7.90e-03  1.96    2.30e-03  1.97    1.51e-02  2.14    1.07e-02  2.12    5.8984e+00
  16   -3.50  -13.00    7.08e-04  1.96    2.04e-03  1.95    5.87e-04  1.97    3.56e-03  2.09    2.56e-03  2.06    2.3915e+01
  32   -4.50  -13.00    1.89e-04  1.91    5.38e-04  1.92    1.50e-04  1.97    8.56e-04  2.06    6.28e-04  2.03    9.7788e+01

 ConvergenceResults: t_end=1.0, example1_manufactured(2.4), Lagrange{1}, Scheme2(), sum(walltime)=130.58 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -13.00    1.92e-02  0.00    4.51e-02  0.00    1.86e-02  0.00    6.68e-02  0.00    4.95e-02  0.00    1.4870e+00
   8   -2.50  -13.00    4.95e-03  1.96    1.11e-02  2.02    4.62e-03  2.01    1.50e-02  2.15    1.15e-02  2.11    5.9075e+00
  16   -3.50  -13.00    1.25e-03  1.99    2.77e-03  2.01    1.15e-03  2.00    3.61e-03  2.06    2.79e-03  2.04    2.3980e+01
  32   -4.50  -13.00    3.13e-04  2.00    6.92e-04  2.00    2.88e-04  2.00    8.88e-04  2.02    6.92e-04  2.01    9.9202e+01

 ConvergenceResults: t_end=1.0, example1_manufactured(2.58), Lagrange{2}, Scheme2(), sum(walltime)=197.15 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -13.00    5.91e-04  0.00    1.45e-03  0.00    6.13e-04  0.00    4.14e-03  0.00    3.75e-03  0.00    1.6047e+00
   8   -2.50  -13.00    8.34e-05  2.83    1.93e-04  2.91    8.14e-05  2.91    6.15e-04  2.75    5.23e-04  2.84    6.7761e+00
  16   -3.50  -13.00    1.06e-05  2.97    2.53e-05  2.93    1.07e-05  2.93    8.37e-05  2.88    6.84e-05  2.93    3.0855e+01
  32   -4.50  -13.00    1.46e-06  2.86    3.31e-06  2.94    1.39e-06  2.94    1.10e-05  2.93    8.77e-06  2.96    1.5792e+02

 ConvergenceResults: t_end=1.0, example1_manufactured(3.4), Lagrange{2}, Scheme2(), sum(walltime)=199.56 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -13.00    1.49e-03  0.00    3.34e-03  0.00    1.63e-03  0.00    5.06e-03  0.00    3.70e-03  0.00    1.6082e+00
   8   -2.50  -13.00    1.82e-04  3.03    4.18e-04  3.00    2.04e-04  3.00    6.40e-04  2.98    4.84e-04  2.93    6.7754e+00
  16   -3.50  -13.00    2.32e-05  2.97    5.23e-05  3.00    2.55e-05  3.00    8.32e-05  2.94    6.16e-05  2.98    3.1009e+01
  32   -4.50  -13.00    2.90e-06  3.00    6.53e-06  3.00    3.19e-06  3.00    1.06e-05  2.97    7.74e-06  2.99    1.6017e+02

 ConvergenceResults: t_end=1.0, example1_manufactured(3.51), Lagrange{3}, Scheme2(), sum(walltime)=399.58 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -13.00    3.05e-05  0.00    7.14e-05  0.00    3.53e-05  0.00    3.97e-04  0.00    1.35e-04  0.00    1.9968e+00
   8   -2.50  -13.00    2.06e-06  3.89    4.89e-06  3.87    2.42e-06  3.87    4.03e-05  3.30    1.13e-05  3.59    9.3074e+00
  16   -3.50  -13.00    1.37e-07  3.91    3.29e-07  3.89    1.63e-07  3.89    2.54e-06  3.99    7.13e-07  3.98    4.8891e+01
  32   -4.50  -13.00    1.01e-08  3.76    2.19e-08  3.91    1.08e-08  3.91    1.49e-07  4.09    4.30e-08  4.05    3.3939e+02

 ConvergenceResults: t_end=1.0, example1_manufactured(4.4), Lagrange{3}, Scheme2(), sum(walltime)=408.18 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -13.00    1.04e-04  0.00    2.21e-04  0.00    1.19e-04  0.00    5.48e-04  0.00    1.38e-04  0.00    2.0011e+00
   8   -2.50  -13.00    6.59e-06  3.97    1.40e-05  3.98    7.52e-06  3.99    1.00e-04  2.45    2.39e-05  2.52    9.2822e+00
  16   -3.50  -13.00    3.79e-07  4.12    8.76e-07  4.00    4.70e-07  4.00    7.04e-06  3.83    1.69e-06  3.83    4.9325e+01
  32   -4.50  -13.00    2.39e-08  3.99    5.48e-08  4.00    2.94e-08  4.00    4.16e-07  4.08    1.00e-07  4.08    3.4757e+02

julia> results1b = run_cases(Scheme2(), cases1b)
6-element Vector{ConvergenceResults}:
 ConvergenceResults: t_end=1.0, example1_manufactured(1.76), Lagrange{1}, Scheme2(), sum(walltime)=2148.41 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -15.00    1.07e-02  0.00    3.08e-02  0.00    9.00e-03  0.00    6.68e-02  0.00    4.64e-02  0.00    5.9226e+00
   8   -2.50  -15.00    2.75e-03  1.97    7.90e-03  1.96    2.30e-03  1.97    1.51e-02  2.14    1.07e-02  2.12    2.3505e+01
  16   -3.50  -15.00    7.08e-04  1.96    2.04e-03  1.95    5.87e-04  1.97    3.56e-03  2.09    2.56e-03  2.06    9.5078e+01
  32   -4.50  -15.00    1.89e-04  1.91    5.38e-04  1.92    1.50e-04  1.97    8.56e-04  2.06    6.28e-04  2.03    3.8875e+02
  64   -5.50  -15.00    6.23e-05  1.60    1.50e-04  1.85    3.87e-05  1.95    2.12e-04  2.01    1.57e-04  2.00    1.6352e+03

 ConvergenceResults: t_end=1.0, example1_manufactured(2.4), Lagrange{1}, Scheme2(), sum(walltime)=2149.68 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -15.00    1.92e-02  0.00    4.51e-02  0.00    1.86e-02  0.00    6.68e-02  0.00    4.95e-02  0.00    5.9112e+00
   8   -2.50  -15.00    4.95e-03  1.96    1.11e-02  2.02    4.62e-03  2.01    1.50e-02  2.15    1.15e-02  2.11    2.3520e+01
  16   -3.50  -15.00    1.25e-03  1.99    2.77e-03  2.01    1.15e-03  2.00    3.61e-03  2.06    2.79e-03  2.04    9.5044e+01
  32   -4.50  -15.00    3.13e-04  2.00    6.92e-04  2.00    2.88e-04  2.00    8.88e-04  2.02    6.92e-04  2.01    3.8871e+02
  64   -5.50  -15.00    7.85e-05  2.00    1.73e-04  2.00    7.19e-05  2.00    2.21e-04  2.01    1.73e-04  2.00    1.6365e+03

 ConvergenceResults: t_end=1.0, example1_manufactured(2.58), Lagrange{2}, Scheme2(), sum(walltime)=4051.93 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -15.00    5.91e-04  0.00    1.45e-03  0.00    6.13e-04  0.00    4.14e-03  0.00    3.75e-03  0.00    6.3051e+00
   8   -2.50  -15.00    8.34e-05  2.83    1.93e-04  2.91    8.14e-05  2.91    6.15e-04  2.75    5.23e-04  2.84    2.6206e+01
  16   -3.50  -15.00    1.06e-05  2.97    2.53e-05  2.93    1.07e-05  2.93    8.37e-05  2.88    6.84e-05  2.93    1.1436e+02
  32   -4.50  -15.00    1.46e-06  2.86    3.31e-06  2.94    1.39e-06  2.94    1.10e-05  2.93    8.77e-06  2.96    5.5922e+02
  64   -5.50  -15.00    2.30e-07  2.67    4.41e-07  2.91    1.79e-07  2.95    1.43e-06  2.94    1.12e-06  2.97    3.3458e+03

 ConvergenceResults: t_end=1.0, example1_manufactured(3.4), Lagrange{2}, Scheme2(), sum(walltime)=4083.87 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -15.00    1.49e-03  0.00    3.34e-03  0.00    1.63e-03  0.00    5.06e-03  0.00    3.70e-03  0.00    6.2977e+00
   8   -2.50  -15.00    1.82e-04  3.03    4.18e-04  3.00    2.04e-04  3.00    6.40e-04  2.98    4.84e-04  2.93    2.6152e+01
  16   -3.50  -15.00    2.32e-05  2.97    5.23e-05  3.00    2.55e-05  3.00    8.32e-05  2.94    6.16e-05  2.98    1.1452e+02
  32   -4.50  -15.00    2.90e-06  3.00    6.53e-06  3.00    3.19e-06  3.00    1.06e-05  2.97    7.74e-06  2.99    5.5872e+02
  64   -5.50  -15.00    3.63e-07  3.00    8.17e-07  3.00    3.99e-07  3.00    1.34e-06  2.99    9.68e-07  3.00    3.3782e+03

 ConvergenceResults: t_end=1.0, example1_manufactured(3.51), Lagrange{3}, Scheme2(), sum(walltime)=10221.69 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -15.00    3.05e-05  0.00    7.14e-05  0.00    3.53e-05  0.00    3.97e-04  0.00    1.35e-04  0.00    7.4039e+00
   8   -2.50  -15.00    2.06e-06  3.89    4.89e-06  3.87    2.42e-06  3.87    4.03e-05  3.30    1.13e-05  3.59    3.3045e+01
  16   -3.50  -15.00    1.37e-07  3.91    3.29e-07  3.89    1.63e-07  3.89    2.54e-06  3.99    7.13e-07  3.98    1.6705e+02
  32   -4.50  -15.00    9.98e-09  3.78    2.19e-08  3.91    1.08e-08  3.91    1.49e-07  4.09    4.30e-08  4.05    1.0759e+03
  64   -5.50  -15.00    9.23e-10  3.44    1.46e-09  3.91    7.16e-10  3.92    8.80e-09  4.08    2.62e-09  4.04    8.9383e+03

 ConvergenceResults: t_end=1.0, example1_manufactured(4.4), Lagrange{3}, Scheme2(), sum(walltime)=10171.95 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -15.00    1.04e-04  0.00    2.21e-04  0.00    1.19e-04  0.00    5.48e-04  0.00    1.38e-04  0.00    7.3927e+00
   8   -2.50  -15.00    6.59e-06  3.97    1.40e-05  3.98    7.52e-06  3.99    1.00e-04  2.45    2.39e-05  2.52    3.3061e+01
  16   -3.50  -15.00    3.79e-07  4.12    8.76e-07  4.00    4.70e-07  4.00    7.04e-06  3.83    1.69e-06  3.83    1.6601e+02
  32   -4.50  -15.00    2.38e-08  4.00    5.48e-08  4.00    2.94e-08  4.00    4.16e-07  4.08    9.99e-08  4.08    1.0723e+03
  64   -5.50  -15.00    1.53e-09  3.96    3.43e-09  4.00    1.84e-09  4.00    2.43e-08  4.10    5.84e-09  4.10    8.8931e+03

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

julia> results2a = run_cases(Scheme2(), cases2a)
6-element Vector{ConvergenceResults}:
 ConvergenceResults: t_end=1.0, example2_manufactured(1.76), Lagrange{1}, Scheme2(), sum(walltime)=128.18 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -13.00    1.11e-02  0.00    3.08e-02  0.00    9.01e-03  0.00    2.21e-02  0.00    2.60e-02  0.00    1.4703e+00
   8   -2.50  -13.00    2.80e-03  1.99    7.96e-03  1.95    2.31e-03  1.97    5.05e-03  2.13    6.06e-03  2.10    5.8825e+00
  16   -3.50  -13.00    7.08e-04  1.98    2.09e-03  1.93    5.90e-04  1.97    1.21e-03  2.07    1.48e-03  2.03    2.3718e+01
  32   -4.50  -13.00    2.11e-04  1.75    5.69e-04  1.87    1.52e-04  1.96    3.10e-04  1.96    3.78e-04  1.97    9.7112e+01

 ConvergenceResults: t_end=1.0, example2_manufactured(2.4), Lagrange{1}, Scheme2(), sum(walltime)=128.32 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -13.00    1.90e-02  0.00    4.45e-02  0.00    1.87e-02  0.00    2.21e-02  0.00    2.80e-02  0.00    1.4735e+00
   8   -2.50  -13.00    4.86e-03  1.96    1.10e-02  2.02    4.63e-03  2.01    5.09e-03  2.12    6.58e-03  2.09    5.8711e+00
  16   -3.50  -13.00    1.23e-03  1.99    2.73e-03  2.01    1.16e-03  2.00    1.23e-03  2.04    1.62e-03  2.03    2.3747e+01
  32   -4.50  -13.00    3.08e-04  1.99    6.82e-04  2.00    2.89e-04  2.00    3.05e-04  2.02    4.02e-04  2.01    9.7232e+01

 ConvergenceResults: t_end=1.0, example2_manufactured(2.58), Lagrange{2}, Scheme2(), sum(walltime)=176.52 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -13.00    5.87e-04  0.00    1.45e-03  0.00    6.13e-04  0.00    1.41e-03  0.00    2.44e-03  0.00    1.5677e+00
   8   -2.50  -13.00    8.31e-05  2.82    1.93e-04  2.91    8.14e-05  2.91    2.03e-04  2.79    3.31e-04  2.88    6.5154e+00
  16   -3.50  -13.00    1.08e-05  2.94    2.54e-05  2.93    1.07e-05  2.93    2.73e-05  2.89    4.27e-05  2.96    2.8719e+01
  32   -4.50  -13.00    1.65e-06  2.71    3.34e-06  2.92    1.39e-06  2.94    3.63e-06  2.91    5.44e-06  2.97    1.3972e+02

 ConvergenceResults: t_end=1.0, example2_manufactured(3.4), Lagrange{2}, Scheme2(), sum(walltime)=176.82 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -13.00    1.47e-03  0.00    3.34e-03  0.00    1.63e-03  0.00    1.66e-03  0.00    2.35e-03  0.00    1.5714e+00
   8   -2.50  -13.00    1.81e-04  3.03    4.18e-04  3.00    2.04e-04  3.00    2.10e-04  2.98    3.12e-04  2.92    6.5030e+00
  16   -3.50  -13.00    2.31e-05  2.96    5.23e-05  3.00    2.55e-05  3.00    2.64e-05  2.99    3.96e-05  2.98    2.8649e+01
  32   -4.50  -13.00    2.89e-06  3.00    6.53e-06  3.00    3.19e-06  3.00    3.31e-06  3.00    4.96e-06  2.99    1.4010e+02

 ConvergenceResults: t_end=1.0, example2_manufactured(3.51), Lagrange{3}, Scheme2(), sum(walltime)=321.20 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -13.00    3.02e-05  0.00    7.14e-05  0.00    3.53e-05  0.00    4.62e-05  0.00    6.94e-05  0.00    1.8553e+00
   8   -2.50  -13.00    2.05e-06  3.88    4.89e-06  3.87    2.42e-06  3.87    3.11e-06  3.89    4.29e-06  4.02    8.2168e+00
  16   -3.50  -13.00    1.38e-07  3.90    3.29e-07  3.89    1.63e-07  3.89    2.04e-07  3.93    2.69e-07  4.00    4.1478e+01
  32   -4.50  -13.00    1.06e-08  3.70    2.19e-08  3.91    1.08e-08  3.91    1.34e-08  3.93    1.70e-08  3.99    2.6965e+02

 ConvergenceResults: t_end=1.0, example2_manufactured(4.4), Lagrange{3}, Scheme2(), sum(walltime)=320.96 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -13.00    1.02e-04  0.00    2.21e-04  0.00    1.19e-04  0.00    4.93e-05  0.00    5.60e-05  0.00    1.8523e+00
   8   -2.50  -13.00    6.56e-06  3.96    1.40e-05  3.98    7.52e-06  3.99    2.94e-06  4.06    3.38e-06  4.05    8.2320e+00
  16   -3.50  -13.00    3.81e-07  4.11    8.76e-07  4.00    4.70e-07  4.00    1.78e-07  4.05    2.09e-07  4.01    4.1742e+01
  32   -4.50  -13.00    2.39e-08  4.00    5.48e-08  4.00    2.94e-08  4.00    1.11e-08  4.00    1.30e-08  4.00    2.6914e+02

julia> results2b = run_cases(Scheme2(), cases2b)
6-element Vector{ConvergenceResults}:
 ConvergenceResults: t_end=1.0, example2_manufactured(1.76), Lagrange{1}, Scheme2(), sum(walltime)=2152.26 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -15.00    1.11e-02  0.00    3.08e-02  0.00    9.01e-03  0.00    2.21e-02  0.00    2.60e-02  0.00    5.8978e+00
   8   -2.50  -15.00    2.80e-03  1.99    7.96e-03  1.95    2.31e-03  1.97    5.05e-03  2.13    6.06e-03  2.10    2.4699e+01
  16   -3.50  -15.00    7.08e-04  1.98    2.09e-03  1.93    5.90e-04  1.97    1.21e-03  2.07    1.48e-03  2.03    9.5093e+01
  32   -4.50  -15.00    2.11e-04  1.75    5.69e-04  1.87    1.52e-04  1.96    3.10e-04  1.96    3.78e-04  1.97    3.8875e+02
  64   -5.50  -15.00    1.17e-04  0.85    1.71e-04  1.73    3.97e-05  1.93    9.99e-05  1.63    1.09e-04  1.80    1.6378e+03

 ConvergenceResults: t_end=1.0, example2_manufactured(2.4), Lagrange{1}, Scheme2(), sum(walltime)=2150.39 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -15.00    1.90e-02  0.00    4.45e-02  0.00    1.87e-02  0.00    2.21e-02  0.00    2.80e-02  0.00    5.8840e+00
   8   -2.50  -15.00    4.86e-03  1.96    1.10e-02  2.02    4.63e-03  2.01    5.09e-03  2.12    6.58e-03  2.09    2.3483e+01
  16   -3.50  -15.00    1.23e-03  1.99    2.73e-03  2.01    1.16e-03  2.00    1.23e-03  2.04    1.62e-03  2.03    9.5037e+01
  32   -4.50  -15.00    3.08e-04  1.99    6.82e-04  2.00    2.89e-04  2.00    3.05e-04  2.02    4.02e-04  2.01    3.8965e+02
  64   -5.50  -15.00    7.72e-05  2.00    1.70e-04  2.00    7.21e-05  2.00    7.59e-05  2.01    1.00e-04  2.00    1.6363e+03

 ConvergenceResults: t_end=1.0, example2_manufactured(2.58), Lagrange{2}, Scheme2(), sum(walltime)=4049.20 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -15.00    5.87e-04  0.00    1.45e-03  0.00    6.13e-04  0.00    1.41e-03  0.00    2.44e-03  0.00    6.2856e+00
   8   -2.50  -15.00    8.31e-05  2.82    1.93e-04  2.91    8.14e-05  2.91    2.03e-04  2.79    3.31e-04  2.88    2.6009e+01
  16   -3.50  -15.00    1.08e-05  2.94    2.54e-05  2.93    1.07e-05  2.93    2.73e-05  2.89    4.27e-05  2.96    1.1455e+02
  32   -4.50  -15.00    1.65e-06  2.71    3.34e-06  2.92    1.39e-06  2.94    3.63e-06  2.91    5.44e-06  2.97    5.5784e+02
  64   -5.50  -15.00    4.59e-07  1.84    4.68e-07  2.84    1.80e-07  2.95    5.18e-07  2.81    7.17e-07  2.92    3.3445e+03

 ConvergenceResults: t_end=1.0, example2_manufactured(3.4), Lagrange{2}, Scheme2(), sum(walltime)=4055.04 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -15.00    1.47e-03  0.00    3.34e-03  0.00    1.63e-03  0.00    1.66e-03  0.00    2.35e-03  0.00    6.2897e+00
   8   -2.50  -15.00    1.81e-04  3.03    4.18e-04  3.00    2.04e-04  3.00    2.10e-04  2.98    3.12e-04  2.92    2.6000e+01
  16   -3.50  -15.00    2.31e-05  2.96    5.23e-05  3.00    2.55e-05  3.00    2.64e-05  2.99    3.96e-05  2.98    1.1475e+02
  32   -4.50  -15.00    2.89e-06  3.00    6.53e-06  3.00    3.19e-06  3.00    3.31e-06  3.00    4.96e-06  2.99    5.5974e+02
  64   -5.50  -15.00    3.62e-07  3.00    8.17e-07  3.00    3.99e-07  3.00    4.13e-07  3.00    6.21e-07  3.00    3.3483e+03

 ConvergenceResults: t_end=1.0, example2_manufactured(3.51), Lagrange{3}, Scheme2(), sum(walltime)=10209.52 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -15.00    3.02e-05  0.00    7.14e-05  0.00    3.53e-05  0.00    4.62e-05  0.00    6.94e-05  0.00    7.3918e+00
   8   -2.50  -15.00    2.05e-06  3.88    4.89e-06  3.87    2.42e-06  3.87    3.11e-06  3.89    4.29e-06  4.02    3.2823e+01
  16   -3.50  -15.00    1.37e-07  3.90    3.29e-07  3.89    1.63e-07  3.89    2.04e-07  3.93    2.69e-07  4.00    1.6644e+02
  32   -4.50  -15.00    1.04e-08  3.72    2.19e-08  3.91    1.08e-08  3.91    1.35e-08  3.92    1.71e-08  3.98    1.0761e+03
  64   -5.50  -15.00    1.13e-09  3.21    1.47e-09  3.90    7.16e-10  3.92    9.26e-10  3.87    1.14e-09  3.91    8.9268e+03

 ConvergenceResults: t_end=1.0, example2_manufactured(4.4), Lagrange{3}, Scheme2(), sum(walltime)=10157.28 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -15.00    1.02e-04  0.00    2.21e-04  0.00    1.19e-04  0.00    4.93e-05  0.00    5.60e-05  0.00    7.4097e+00
   8   -2.50  -15.00    6.56e-06  3.96    1.40e-05  3.98    7.52e-06  3.99    2.94e-06  4.06    3.38e-06  4.05    3.2970e+01
  16   -3.50  -15.00    3.81e-07  4.11    8.76e-07  4.00    4.70e-07  4.00    1.78e-07  4.05    2.09e-07  4.01    1.6571e+02
  32   -4.50  -15.00    2.38e-08  4.00    5.48e-08  4.00    2.94e-08  4.00    1.09e-08  4.03    1.30e-08  4.00    1.0715e+03
  64   -5.50  -15.00    1.53e-09  3.96    3.43e-09  4.00    1.84e-09  4.00    6.86e-10  3.99    8.15e-10  4.00    8.8797e+03

julia> cases3 = (
    (fe = Lagrange{1}, id = example1_manufactured(2.4 ), Nx = 2^9, τ = [2.0^-i for i in 2:5]),
    (fe = Lagrange{1}, id = example2_manufactured(2.4 ), Nx = 2^9, τ = [2.0^-i for i in 2:5])
    );

julia> results3 = run_cases(Scheme2(), cases3)
2-element Vector{ConvergenceResults}:
 ConvergenceResults: t_end=1.0, example1_manufactured(2.4), Lagrange{1}, Scheme2(), sum(walltime)=816.58 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
 512   -8.50   -2.00    1.15e-02  0.00    2.67e-03  0.00    1.45e-03  0.00    8.41e-03  0.00    5.85e-03  0.00    7.0936e+01
 512   -8.50   -3.00    2.21e-03  2.38    5.32e-04  2.33    3.01e-04  2.27    2.09e-03  2.01    1.50e-03  1.97    1.2001e+02
 512   -8.50   -4.00    4.45e-04  2.32    1.21e-04  2.13    6.95e-05  2.11    5.20e-04  2.00    3.77e-04  1.99    2.2188e+02
 512   -8.50   -5.00    1.02e-04  2.12    2.97e-05  2.03    1.77e-05  1.97    1.30e-04  2.00    9.48e-05  1.99    4.0375e+02

 ConvergenceResults: t_end=1.0, example2_manufactured(2.4), Lagrange{1}, Scheme2(), sum(walltime)=613.43 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
 512   -8.50   -2.00    1.17e-02  0.00    2.81e-03  0.00    1.48e-03  0.00    3.96e-03  0.00    1.55e-03  0.00    5.0392e+01
 512   -8.50   -3.00    2.26e-03  2.37    5.54e-04  2.35    3.03e-04  2.29    9.45e-04  2.07    4.34e-04  1.84    8.7771e+01
 512   -8.50   -4.00    4.58e-04  2.30    1.24e-04  2.16    6.83e-05  2.15    2.31e-04  2.04    1.13e-04  1.94    1.6353e+02
 512   -8.50   -5.00    1.05e-04  2.12    3.01e-05  2.04    1.76e-05  1.96    5.75e-05  2.00    2.88e-05  1.97    3.1173e+02

julia> versioninfo()
Julia Version 1.12.6
Commit 15346901f00 (2026-04-09 19:20 UTC)
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