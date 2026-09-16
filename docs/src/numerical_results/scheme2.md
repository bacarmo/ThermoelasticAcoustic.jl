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
 ConvergenceResults: t_end=1.0, example1_manufactured(1.76), Lagrange{1}, Scheme2(), sum(walltime)=182.46 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -13.00    1.07e-02  0.00    3.08e-02  0.00    8.99e-03  0.00    6.68e-02  0.00    4.64e-02  0.00    2.1035e+00
   8   -2.50  -13.00    2.74e-03  1.97    7.88e-03  1.97    2.30e-03  1.97    1.51e-02  2.14    1.07e-02  2.12    8.5534e+00
  16   -3.50  -13.00    7.03e-04  1.96    2.03e-03  1.96    5.86e-04  1.97    3.56e-03  2.09    2.56e-03  2.06    3.4105e+01
  32   -4.50  -13.00    1.84e-04  1.93    5.29e-04  1.94    1.49e-04  1.97    8.54e-04  2.06    6.27e-04  2.03    1.3770e+02

 ConvergenceResults: t_end=1.0, example1_manufactured(2.4), Lagrange{1}, Scheme2(), sum(walltime)=183.76 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -13.00    1.93e-02  0.00    4.51e-02  0.00    1.86e-02  0.00    6.68e-02  0.00    4.95e-02  0.00    2.0999e+00
   8   -2.50  -13.00    4.95e-03  1.96    1.11e-02  2.02    4.62e-03  2.01    1.50e-02  2.15    1.15e-02  2.11    8.4428e+00
  16   -3.50  -13.00    1.25e-03  1.99    2.77e-03  2.01    1.15e-03  2.00    3.61e-03  2.06    2.79e-03  2.04    3.4035e+01
  32   -4.50  -13.00    3.14e-04  2.00    6.92e-04  2.00    2.88e-04  2.00    8.88e-04  2.02    6.92e-04  2.01    1.3919e+02

 ConvergenceResults: t_end=1.0, example1_manufactured(2.58), Lagrange{2}, Scheme2(), sum(walltime)=316.60 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -13.00    5.92e-04  0.00    1.45e-03  0.00    6.13e-04  0.00    4.14e-03  0.00    3.75e-03  0.00    3.0380e+00
   8   -2.50  -13.00    8.39e-05  2.82    1.93e-04  2.91    8.14e-05  2.91    6.15e-04  2.75    5.23e-04  2.84    1.2490e+01
  16   -3.50  -13.00    1.06e-05  2.98    2.53e-05  2.93    1.07e-05  2.93    8.37e-05  2.88    6.84e-05  2.93    5.3456e+01
  32   -4.50  -13.00    1.40e-06  2.92    3.29e-06  2.94    1.39e-06  2.94    1.10e-05  2.93    8.76e-06  2.96    2.4761e+02

 ConvergenceResults: t_end=1.0, example1_manufactured(3.4), Lagrange{2}, Scheme2(), sum(walltime)=319.75 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -13.00    1.49e-03  0.00    3.34e-03  0.00    1.63e-03  0.00    5.06e-03  0.00    3.70e-03  0.00    3.0292e+00
   8   -2.50  -13.00    1.82e-04  3.03    4.18e-04  3.00    2.04e-04  3.00    6.40e-04  2.98    4.84e-04  2.93    1.2443e+01
  16   -3.50  -13.00    2.32e-05  2.97    5.23e-05  3.00    2.55e-05  3.00    8.32e-05  2.94    6.16e-05  2.98    5.3679e+01
  32   -4.50  -13.00    2.90e-06  3.00    6.53e-06  3.00    3.19e-06  3.00    1.06e-05  2.97    7.74e-06  2.99    2.5060e+02

 ConvergenceResults: t_end=1.0, example1_manufactured(3.51), Lagrange{3}, Scheme2(), sum(walltime)=606.71 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -13.00    3.15e-05  0.00    7.15e-05  0.00    3.54e-05  0.00    4.05e-04  0.00    1.37e-04  0.00    4.3618e+00
   8   -2.50  -13.00    2.07e-06  3.93    4.89e-06  3.87    2.42e-06  3.87    4.05e-05  3.32    1.13e-05  3.60    1.8778e+01
  16   -3.50  -13.00    1.37e-07  3.92    3.29e-07  3.89    1.63e-07  3.90    2.54e-06  3.99    7.14e-07  3.98    8.6350e+01
  32   -4.50  -13.00    9.32e-09  3.88    2.19e-08  3.91    1.08e-08  3.91    1.49e-07  4.09    4.31e-08  4.05    4.9722e+02

 ConvergenceResults: t_end=1.0, example1_manufactured(4.4), Lagrange{3}, Scheme2(), sum(walltime)=614.89 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -13.00    1.04e-04  0.00    2.21e-04  0.00    1.19e-04  0.00    5.88e-04  0.00    1.47e-04  0.00    4.3744e+00
   8   -2.50  -13.00    6.59e-06  3.97    1.40e-05  3.98    7.52e-06  3.99    1.01e-04  2.55    2.41e-05  2.61    1.8782e+01
  16   -3.50  -13.00    3.79e-07  4.12    8.76e-07  4.00    4.70e-07  4.00    7.05e-06  3.84    1.69e-06  3.84    8.7123e+01
  32   -4.50  -13.00    2.39e-08  3.99    5.48e-08  4.00    2.94e-08  4.00    4.16e-07  4.08    1.00e-07  4.08    5.0461e+02

julia> results1b = run_cases(Scheme2(), cases1b)
6-element Vector{ConvergenceResults}:
 ConvergenceResults: t_end=1.0, example1_manufactured(1.76), Lagrange{1}, Scheme2(), sum(walltime)=3002.64 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -15.00    1.07e-02  0.00    3.08e-02  0.00    8.99e-03  0.00    6.68e-02  0.00    4.64e-02  0.00    8.3970e+00
   8   -2.50  -15.00    2.74e-03  1.97    7.88e-03  1.97    2.30e-03  1.97    1.51e-02  2.14    1.07e-02  2.12    3.3745e+01
  16   -3.50  -15.00    7.03e-04  1.96    2.03e-03  1.96    5.86e-04  1.97    3.56e-03  2.09    2.56e-03  2.06    1.3550e+02
  32   -4.50  -15.00    1.84e-04  1.93    5.29e-04  1.94    1.49e-04  1.97    8.54e-04  2.06    6.27e-04  2.03    5.4868e+02
  64   -5.50  -15.00    5.44e-05  1.76    1.44e-04  1.88    3.83e-05  1.96    2.10e-04  2.03    1.56e-04  2.01    2.2763e+03

 ConvergenceResults: t_end=1.0, example1_manufactured(2.4), Lagrange{1}, Scheme2(), sum(walltime)=3000.95 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -15.00    1.93e-02  0.00    4.51e-02  0.00    1.86e-02  0.00    6.68e-02  0.00    4.95e-02  0.00    8.3761e+00
   8   -2.50  -15.00    4.95e-03  1.96    1.11e-02  2.02    4.62e-03  2.01    1.50e-02  2.15    1.15e-02  2.11    3.3625e+01
  16   -3.50  -15.00    1.25e-03  1.99    2.77e-03  2.01    1.15e-03  2.00    3.61e-03  2.06    2.79e-03  2.04    1.3539e+02
  32   -4.50  -15.00    3.14e-04  2.00    6.92e-04  2.00    2.88e-04  2.00    8.88e-04  2.02    6.92e-04  2.01    5.4870e+02
  64   -5.50  -15.00    7.85e-05  2.00    1.73e-04  2.00    7.19e-05  2.00    2.21e-04  2.01    1.73e-04  2.00    2.2749e+03

 ConvergenceResults: t_end=1.0, example1_manufactured(2.58), Lagrange{2}, Scheme2(), sum(walltime)=5980.93 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -15.00    5.92e-04  0.00    1.45e-03  0.00    6.13e-04  0.00    4.14e-03  0.00    3.75e-03  0.00    1.1981e+01
   8   -2.50  -15.00    8.39e-05  2.82    1.93e-04  2.91    8.14e-05  2.91    6.15e-04  2.75    5.23e-04  2.84    4.8771e+01
  16   -3.50  -15.00    1.06e-05  2.98    2.53e-05  2.93    1.07e-05  2.93    8.37e-05  2.88    6.84e-05  2.93    2.0481e+02
  32   -4.50  -15.00    1.40e-06  2.92    3.29e-06  2.94    1.39e-06  2.94    1.10e-05  2.93    8.76e-06  2.96    9.1784e+02
  64   -5.50  -15.00    1.87e-07  2.91    4.27e-07  2.95    1.79e-07  2.95    1.42e-06  2.95    1.11e-06  2.98    4.7975e+03

 ConvergenceResults: t_end=1.0, example1_manufactured(3.4), Lagrange{2}, Scheme2(), sum(walltime)=5992.80 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -15.00    1.49e-03  0.00    3.34e-03  0.00    1.63e-03  0.00    5.06e-03  0.00    3.70e-03  0.00    1.1983e+01
   8   -2.50  -15.00    1.82e-04  3.03    4.18e-04  3.00    2.04e-04  3.00    6.40e-04  2.98    4.84e-04  2.93    4.8709e+01
  16   -3.50  -15.00    2.32e-05  2.97    5.23e-05  3.00    2.55e-05  3.00    8.32e-05  2.94    6.16e-05  2.98    2.0538e+02
  32   -4.50  -15.00    2.90e-06  3.00    6.53e-06  3.00    3.19e-06  3.00    1.06e-05  2.97    7.74e-06  2.99    9.1840e+02
  64   -5.50  -15.00    3.63e-07  3.00    8.17e-07  3.00    3.99e-07  3.00    1.34e-06  2.99    9.68e-07  3.00    4.8083e+03

 ConvergenceResults: t_end=1.0, example1_manufactured(3.51), Lagrange{3}, Scheme2(), sum(walltime)=13355.81 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -15.00    3.15e-05  0.00    7.15e-05  0.00    3.54e-05  0.00    4.05e-04  0.00    1.37e-04  0.00    1.6832e+01
   8   -2.50  -15.00    2.07e-06  3.93    4.89e-06  3.87    2.42e-06  3.87    4.05e-05  3.32    1.13e-05  3.60    7.0612e+01
  16   -3.50  -15.00    1.37e-07  3.92    3.29e-07  3.89    1.63e-07  3.90    2.54e-06  3.99    7.14e-07  3.98    3.1577e+02
  32   -4.50  -15.00    9.24e-09  3.89    2.19e-08  3.91    1.08e-08  3.91    1.49e-07  4.09    4.30e-08  4.05    1.6932e+03
  64   -5.50  -15.00    6.13e-10  3.91    1.45e-09  3.92    7.16e-10  3.92    8.81e-09  4.08    2.62e-09  4.04    1.1259e+04

 ConvergenceResults: t_end=1.0, example1_manufactured(4.4), Lagrange{3}, Scheme2(), sum(walltime)=13397.63 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -15.00    1.04e-04  0.00    2.21e-04  0.00    1.19e-04  0.00    5.88e-04  0.00    1.47e-04  0.00    1.6823e+01
   8   -2.50  -15.00    6.59e-06  3.97    1.40e-05  3.98    7.52e-06  3.99    1.01e-04  2.55    2.41e-05  2.61    7.0843e+01
  16   -3.50  -15.00    3.79e-07  4.12    8.76e-07  4.00    4.70e-07  4.00    7.05e-06  3.84    1.69e-06  3.84    3.1719e+02
  32   -4.50  -15.00    2.38e-08  4.00    5.48e-08  4.00    2.94e-08  4.00    4.16e-07  4.08    1.00e-07  4.08    1.6923e+03
  64   -5.50  -15.00    1.53e-09  3.96    3.43e-09  4.00    1.84e-09  4.00    2.43e-08  4.10    5.84e-09  4.10    1.1300e+04

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
 ConvergenceResults: t_end=1.0, example2_manufactured(1.76), Lagrange{1}, Scheme2(), sum(walltime)=181.45 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -13.00    1.11e-02  0.00    3.08e-02  0.00    9.01e-03  0.00    2.21e-02  0.00    2.60e-02  0.00    2.0913e+00
   8   -2.50  -13.00    2.80e-03  1.98    7.92e-03  1.96    2.30e-03  1.97    5.04e-03  2.13    6.06e-03  2.10    8.3997e+00
  16   -3.50  -13.00    7.07e-04  1.99    2.06e-03  1.94    5.89e-04  1.97    1.20e-03  2.07    1.48e-03  2.03    3.3827e+01
  32   -4.50  -13.00    1.87e-04  1.92    5.53e-04  1.90    1.51e-04  1.97    3.01e-04  2.00    3.73e-04  1.99    1.3713e+02

 ConvergenceResults: t_end=1.0, example2_manufactured(2.4), Lagrange{1}, Scheme2(), sum(walltime)=181.48 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -13.00    1.90e-02  0.00    4.45e-02  0.00    1.87e-02  0.00    2.21e-02  0.00    2.80e-02  0.00    2.0900e+00
   8   -2.50  -13.00    4.86e-03  1.96    1.10e-02  2.02    4.63e-03  2.01    5.09e-03  2.12    6.58e-03  2.09    8.4067e+00
  16   -3.50  -13.00    1.23e-03  1.99    2.73e-03  2.01    1.16e-03  2.00    1.23e-03  2.04    1.62e-03  2.03    3.3808e+01
  32   -4.50  -13.00    3.08e-04  1.99    6.82e-04  2.00    2.89e-04  2.00    3.05e-04  2.02    4.02e-04  2.01    1.3717e+02

 ConvergenceResults: t_end=1.0, example2_manufactured(2.58), Lagrange{2}, Scheme2(), sum(walltime)=296.77 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -13.00    5.88e-04  0.00    1.45e-03  0.00    6.13e-04  0.00    1.41e-03  0.00    2.44e-03  0.00    2.9935e+00
   8   -2.50  -13.00    8.34e-05  2.82    1.93e-04  2.91    8.14e-05  2.91    2.03e-04  2.80    3.31e-04  2.88    1.2164e+01
  16   -3.50  -13.00    1.06e-05  2.98    2.53e-05  2.93    1.07e-05  2.93    2.73e-05  2.89    4.26e-05  2.96    5.1345e+01
  32   -4.50  -13.00    1.43e-06  2.89    3.30e-06  2.94    1.39e-06  2.94    3.57e-06  2.93    5.40e-06  2.98    2.3027e+02

 ConvergenceResults: t_end=1.0, example2_manufactured(3.4), Lagrange{2}, Scheme2(), sum(walltime)=297.02 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -13.00    1.47e-03  0.00    3.34e-03  0.00    1.63e-03  0.00    1.66e-03  0.00    2.35e-03  0.00    2.9899e+00
   8   -2.50  -13.00    1.80e-04  3.03    4.18e-04  3.00    2.04e-04  3.00    2.10e-04  2.98    3.12e-04  2.92    1.2151e+01
  16   -3.50  -13.00    2.31e-05  2.96    5.23e-05  3.00    2.55e-05  3.00    2.64e-05  2.99    3.96e-05  2.98    5.1245e+01
  32   -4.50  -13.00    2.89e-06  3.00    6.53e-06  3.00    3.19e-06  3.00    3.31e-06  3.00    4.96e-06  2.99    2.3063e+02

 ConvergenceResults: t_end=1.0, example2_manufactured(3.51), Lagrange{3}, Scheme2(), sum(walltime)=524.52 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -13.00    3.13e-05  0.00    7.15e-05  0.00    3.54e-05  0.00    4.64e-05  0.00    6.93e-05  0.00    4.1981e+00
   8   -2.50  -13.00    2.06e-06  3.92    4.90e-06  3.87    2.42e-06  3.87    3.12e-06  3.90    4.29e-06  4.02    1.7679e+01
  16   -3.50  -13.00    1.37e-07  3.91    3.29e-07  3.89    1.63e-07  3.90    2.05e-07  3.93    2.68e-07  4.00    7.9112e+01
  32   -4.50  -13.00    9.29e-09  3.88    2.19e-08  3.91    1.08e-08  3.91    1.35e-08  3.93    1.69e-08  3.99    4.2353e+02

 ConvergenceResults: t_end=1.0, example2_manufactured(4.4), Lagrange{3}, Scheme2(), sum(walltime)=524.43 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -13.00    1.02e-04  0.00    2.21e-04  0.00    1.19e-04  0.00    4.92e-05  0.00    5.59e-05  0.00    4.1889e+00
   8   -2.50  -13.00    6.55e-06  3.96    1.40e-05  3.98    7.52e-06  3.99    2.94e-06  4.06    3.38e-06  4.05    1.7639e+01
  16   -3.50  -13.00    3.81e-07  4.11    8.76e-07  4.00    4.70e-07  4.00    1.78e-07  4.05    2.09e-07  4.01    7.9846e+01
  32   -4.50  -13.00    2.39e-08  4.00    5.48e-08  4.00    2.94e-08  4.00    1.11e-08  4.00    1.30e-08  4.00    4.2275e+02

julia> results2b = run_cases(Scheme2(), cases2b)
6-element Vector{ConvergenceResults}:
 ConvergenceResults: t_end=1.0, example2_manufactured(1.76), Lagrange{1}, Scheme2(), sum(walltime)=2999.06 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -15.00    1.11e-02  0.00    3.08e-02  0.00    9.01e-03  0.00    2.21e-02  0.00    2.60e-02  0.00    8.3817e+00
   8   -2.50  -15.00    2.80e-03  1.98    7.92e-03  1.96    2.30e-03  1.97    5.04e-03  2.13    6.06e-03  2.10    3.3612e+01
  16   -3.50  -15.00    7.07e-04  1.99    2.06e-03  1.94    5.89e-04  1.97    1.20e-03  2.07    1.48e-03  2.03    1.3533e+02
  32   -4.50  -15.00    1.87e-04  1.92    5.53e-04  1.90    1.51e-04  1.97    3.01e-04  2.00    3.73e-04  1.99    5.4869e+02
  64   -5.50  -15.00    9.10e-05  1.04    1.60e-04  1.79    3.91e-05  1.95    8.89e-05  1.76    1.02e-04  1.87    2.2730e+03

 ConvergenceResults: t_end=1.0, example2_manufactured(2.4), Lagrange{1}, Scheme2(), sum(walltime)=2997.71 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -15.00    1.90e-02  0.00    4.45e-02  0.00    1.87e-02  0.00    2.21e-02  0.00    2.80e-02  0.00    8.3782e+00
   8   -2.50  -15.00    4.86e-03  1.96    1.10e-02  2.02    4.63e-03  2.01    5.09e-03  2.12    6.58e-03  2.09    3.3603e+01
  16   -3.50  -15.00    1.23e-03  1.99    2.73e-03  2.01    1.16e-03  2.00    1.23e-03  2.04    1.62e-03  2.03    1.3528e+02
  32   -4.50  -15.00    3.08e-04  1.99    6.82e-04  2.00    2.89e-04  2.00    3.05e-04  2.02    4.02e-04  2.01    5.4833e+02
  64   -5.50  -15.00    7.71e-05  2.00    1.70e-04  2.00    7.21e-05  2.00    7.58e-05  2.01    1.00e-04  2.00    2.2721e+03

 ConvergenceResults: t_end=1.0, example2_manufactured(2.58), Lagrange{2}, Scheme2(), sum(walltime)=5974.13 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -15.00    5.88e-04  0.00    1.45e-03  0.00    6.13e-04  0.00    1.41e-03  0.00    2.44e-03  0.00    1.1971e+01
   8   -2.50  -15.00    8.34e-05  2.82    1.93e-04  2.91    8.14e-05  2.91    2.03e-04  2.80    3.31e-04  2.88    4.8535e+01
  16   -3.50  -15.00    1.06e-05  2.98    2.53e-05  2.93    1.07e-05  2.93    2.73e-05  2.89    4.26e-05  2.96    2.0493e+02
  32   -4.50  -15.00    1.43e-06  2.89    3.30e-06  2.94    1.39e-06  2.94    3.57e-06  2.93    5.40e-06  2.98    9.2167e+02
  64   -5.50  -15.00    2.10e-07  2.77    4.31e-07  2.93    1.79e-07  2.95    4.69e-07  2.93    6.84e-07  2.98    4.7870e+03

 ConvergenceResults: t_end=1.0, example2_manufactured(3.4), Lagrange{2}, Scheme2(), sum(walltime)=5987.25 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -15.00    1.47e-03  0.00    3.34e-03  0.00    1.63e-03  0.00    1.66e-03  0.00    2.35e-03  0.00    1.1959e+01
   8   -2.50  -15.00    1.80e-04  3.03    4.18e-04  3.00    2.04e-04  3.00    2.10e-04  2.98    3.12e-04  2.92    4.8678e+01
  16   -3.50  -15.00    2.31e-05  2.96    5.23e-05  3.00    2.55e-05  3.00    2.64e-05  2.99    3.96e-05  2.98    2.0488e+02
  32   -4.50  -15.00    2.89e-06  3.00    6.53e-06  3.00    3.19e-06  3.00    3.31e-06  3.00    4.96e-06  2.99    9.2081e+02
  64   -5.50  -15.00    3.62e-07  3.00    8.17e-07  3.00    3.99e-07  3.00    4.13e-07  3.00    6.21e-07  3.00    4.8009e+03

 ConvergenceResults: t_end=1.0, example2_manufactured(3.51), Lagrange{3}, Scheme2(), sum(walltime)=13351.90 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -15.00    3.13e-05  0.00    7.15e-05  0.00    3.54e-05  0.00    4.64e-05  0.00    6.93e-05  0.00    1.6796e+01
   8   -2.50  -15.00    2.06e-06  3.92    4.90e-06  3.87    2.42e-06  3.87    3.12e-06  3.90    4.29e-06  4.02    7.0744e+01
  16   -3.50  -15.00    1.37e-07  3.91    3.29e-07  3.89    1.63e-07  3.90    2.05e-07  3.93    2.69e-07  4.00    3.1603e+02
  32   -4.50  -15.00    9.21e-09  3.89    2.19e-08  3.91    1.08e-08  3.91    1.34e-08  3.93    1.69e-08  3.99    1.6938e+03
  64   -5.50  -15.00    6.16e-10  3.90    1.45e-09  3.92    7.16e-10  3.92    8.75e-10  3.94    1.06e-09  3.99    1.1255e+04

 ConvergenceResults: t_end=1.0, example2_manufactured(4.4), Lagrange{3}, Scheme2(), sum(walltime)=13351.41 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
   4   -1.50  -15.00    1.02e-04  0.00    2.21e-04  0.00    1.19e-04  0.00    4.92e-05  0.00    5.59e-05  0.00    1.6814e+01
   8   -2.50  -15.00    6.55e-06  3.96    1.40e-05  3.98    7.52e-06  3.99    2.94e-06  4.06    3.38e-06  4.05    7.0688e+01
  16   -3.50  -15.00    3.81e-07  4.11    8.76e-07  4.00    4.70e-07  4.00    1.78e-07  4.05    2.09e-07  4.01    3.1710e+02
  32   -4.50  -15.00    2.38e-08  4.00    5.48e-08  4.00    2.94e-08  4.00    1.09e-08  4.03    1.30e-08  4.00    1.6914e+03
  64   -5.50  -15.00    1.53e-09  3.96    3.43e-09  4.00    1.84e-09  4.00    6.86e-10  3.99    8.15e-10  4.00    1.1255e+04

julia> cases3 = (
    (fe = Lagrange{1}, id = example1_manufactured(2.4 ), Nx = 2^9, τ = [2.0^-i for i in 2:5]),
    (fe = Lagrange{1}, id = example2_manufactured(2.4 ), Nx = 2^9, τ = [2.0^-i for i in 2:5])
    );

julia> results3 = run_cases(Scheme2(), cases3)
2-element Vector{ConvergenceResults}:
 ConvergenceResults: t_end=1.0, example1_manufactured(2.4), Lagrange{1}, Scheme2(), sum(walltime)=896.91 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
 512   -8.50   -2.00    1.15e-02  0.00    2.67e-03  0.00    1.45e-03  0.00    8.41e-03  0.00    5.85e-03  0.00    7.7026e+01
 512   -8.50   -3.00    2.21e-03  2.38    5.32e-04  2.33    3.01e-04  2.27    2.09e-03  2.01    1.50e-03  1.97    1.3066e+02
 512   -8.50   -4.00    4.45e-04  2.32    1.21e-04  2.13    6.95e-05  2.11    5.20e-04  2.00    3.77e-04  1.99    2.4558e+02
 512   -8.50   -5.00    1.02e-04  2.12    2.97e-05  2.03    1.77e-05  1.97    1.30e-04  2.00    9.48e-05  1.99    4.4364e+02

 ConvergenceResults: t_end=1.0, example2_manufactured(2.4), Lagrange{1}, Scheme2(), sum(walltime)=698.77 s
  Nx   log₂h   log₂τ    L∞L²_V    rate    L∞L²_U    rate    L∞L²_Θ    rate    L∞L²_R    rate    L∞L²_Z    rate    walltime[s]
 512   -8.50   -2.00    1.17e-02  0.00    2.81e-03  0.00    1.48e-03  0.00    3.96e-03  0.00    1.55e-03  0.00    5.6894e+01
 512   -8.50   -3.00    2.26e-03  2.37    5.54e-04  2.35    3.03e-04  2.29    9.45e-04  2.07    4.34e-04  1.84    9.9650e+01
 512   -8.50   -4.00    4.58e-04  2.30    1.24e-04  2.16    6.83e-05  2.15    2.31e-04  2.04    1.13e-04  1.94    1.8480e+02
 512   -8.50   -5.00    1.05e-04  2.12    3.01e-05  2.04    1.76e-05  1.96    5.74e-05  2.01    2.88e-05  1.97    3.5743e+02

julia> versioninfo()
Julia Version 1.13.0
Commit d1c37793dd2 (2026-09-09 19:00 UTC)
Build Info:
  Official https://julialang.org release
Platform Info:
  OS: Linux (x86_64-linux-gnu)
  CPU: 14 × Intel(R) Core(TM) Ultra 5 225H
  WORD_SIZE: 64
  LLVM: libLLVM-20.1.8 (ORCJIT, arrowlake)
  GC: Built with stock GC
Threads: 1 default, 1 interactive, 1 GC (on 14 virtual cores)
```