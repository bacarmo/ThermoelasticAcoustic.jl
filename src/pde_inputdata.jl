"""
    PDEInputData{Tα, Tβ, Tdβ, Tf, Tdf, Tg, T∂ₛg,
                 Tu₀, T∂ₓu₀, T∂ᵧu₀,
                 Tv₀, T∂ₓv₀, T∂ᵧv₀,
                 Tθ₀, T∂ₓθ₀, T∂ᵧθ₀,
                 Tz₀, Tr₀,
                 Tf₁, Tf₂, Tf₃,
                 Tu, Tv, Tθ, Tz, Tr}

Input data container for a coupled wave–thermal–acoustic PDE system posed on a
rectangular domain ``Ω ⊂ ℝ²``, with acoustic interaction prescribed on the lower
boundary ``Γ₁``.

# Mathematical Model

Wave equation in Ω:
```math
\\frac{∂²u}{∂t²}(x,y,t) - α(t)\\,Δu(x,y,t) + f(u(x,y,t)) + (\\mathbf{a} \\cdot \\nabla)θ(x,y,t) = f₁(x,y,t)
```

Thermal equation in Ω:
```math
\\frac{∂θ}{∂t}(x,y,t) - β\\!\\left(\\int_Ω θ(x,y,t)\\,\\mathrm{d}x\\,\\mathrm{d}y\\right) Δθ(x,y,t)
+ (\\mathbf{a} \\cdot \\nabla)\\frac{∂u}{∂t}(x,y,t) = f₂(x,y,t)
```

Acoustic boundary equation on Γ₁:
```math
q_1\\,\\frac{∂²z}{∂t²}(x,t) + q_2\\,\\frac{∂z}{∂t}(x,t) + q_3\\,z(x,t)
+ q_4\\,\\frac{∂u}{∂t}(x,y_{\\min},t) = f₃(x,t)
```

Nonlinear wave–acoustic coupling on Γ₁:
```math
-\\frac{∂u}{∂y}(x,y_{\\min},t) = \\frac{∂z}{∂t}(x,t) - g\\left(x,\\,\\frac{∂u}{∂t}(x,y_{\\min},t)\\right)
```

Boundary conditions:
```math
u = 0 \\text{ on } Γ_0, \\qquad θ = 0 \\text{ on } Γ
```

Initial conditions:
```math
u(x,y,0) = u_0(x,y), \\quad \\frac{∂u}{∂t}(x,y,0) = v_0(x,y), \\quad
θ(x,y,0) = θ_0(x,y), \\quad z(x,0) = z_0(x), \\quad \\frac{∂z}{∂t}(x,0) = r_0(x)
```

## Notes
- The rectangular domain is ``Ω = ]x_{\\min}, x_{\\max}[ × ]y_{\\min}, y_{\\max}[``.
- The lower boundary is ``Γ_1 = \\{(x, y_{\\min}) : x \\in ]x_{\\min}, x_{\\max}[\\}``.
- ``Γ`` denotes the full boundary of ``Ω`` and ``Γ_0 = Γ \\setminus Γ_1``.

# Fields

## Domain Geometry
- `pmin::NTuple{2,Float64}`: Bottom-left corner ``(x_{\\min}, y_{\\min})`` of ``Ω``.
- `pmax::NTuple{2,Float64}`: Top-right corner ``(x_{\\max}, y_{\\max})`` of ``Ω``.

## Physical Parameters
- `a::NTuple{2,Float64}`: Wave-thermal coupling strength
- `q₁::Float64`: Acoustic acceleration coefficient
- `q₂::Float64`: Acoustic velocity coefficient
- `q₃::Float64`: Acoustic displacement coefficient
- `q₄::Float64`: Wave-acoustic coupling strength

## Coefficient Functions
- `α::Tα`: Time-dependent wave diffusion coefficient α(t) > 0.
- `β::Tβ`, `dβ::Tdβ`: Nonlinear thermal diffusion coefficient β(s) and its derivative β′(s).
- `f::Tf`, `df::Tdf`: Nonlinear wave term f(s) and its derivative f′(s).
- `g::Tg`, `∂ₛg::T∂ₛg`: Nonlinear coupling function g(x, s) and its partial derivative ∂ₛg(x, s).

## Wave Initial Conditions (functions on Ω)
- `u₀::Tu₀`: Initial displacement u(x, y, 0).
- `∂ₓu₀::T∂ₓu₀`, `∂ᵧu₀::T∂ᵧu₀`: Partial derivatives of u₀.
- `v₀::Tv₀`: Initial velocity ∂ₜu(x, y, 0).
- `∂ₓv₀::T∂ₓv₀`, `∂ᵧv₀::T∂ᵧv₀`: Partial derivatives of v₀.

## Thermal Initial Conditions (functions on Ω)
- `θ₀::Tθ₀`: Initial temperature θ(x, y, 0).
- `∂ₓθ₀::T∂ₓθ₀`, `∂ᵧθ₀::T∂ᵧθ₀`: Partial derivatives of θ₀.

## Acoustic Initial Conditions (functions on Γ₁)
- `z₀::Tz₀`: Initial acoustic displacement z(x, 0).
- `r₀::Tr₀`: Initial acoustic velocity ∂ₜz(x, 0). Must satisfy ``r₀(x) = -\\frac{∂u₀}{∂y}(x, y_{\\min}) + g(x, v₀(x, y_{\\min}))``

## Source Terms
- `f₁::Tf₁`: Wave source f₁(x, y, t) on Ω.
- `f₂::Tf₂`: Thermal source f₂(x, y, t) on Ω.
- `f₃::Tf₃`: Acoustic source f₃(x, t) on Γ₁.

## Analytical Solutions
Optional fields for manufactured-solution convergence studies. Set to `nothing`
for physical simulations without known closed-form solutions.
- `u::Tu`, `v::Tv`: Analytical wave displacement u(x, y, t) and velocity v(x, y, t).
- `θ::Tθ`: Analytical thermal solution θ(x, y, t).
- `z::Tz`, `r::Tr`: Analytical acoustic displacement z(x, t) and velocity r(x, t).
"""
struct PDEInputData{
    Tα, Tβ, Tdβ, Tf, Tdf, Tg, T∂ₛg,
    Tu₀, T∂ₓu₀, T∂ᵧu₀,
    Tv₀, T∂ₓv₀, T∂ᵧv₀,
    Tθ₀, T∂ₓθ₀, T∂ᵧθ₀,
    Tz₀, Tr₀,
    Tf₁, Tf₂, Tf₃,
    Tu, Tv, Tθ, Tz, Tr}

    # Domain geometry
    pmin::NTuple{2, Float64}
    pmax::NTuple{2, Float64}

    # Physical parameters
    a::NTuple{2, Float64}
    q₁::Float64
    q₂::Float64
    q₃::Float64
    q₄::Float64

    # Coefficient functions
    α::Tα
    β::Tβ
    dβ::Tdβ
    f::Tf
    df::Tdf
    g::Tg
    ∂ₛg::T∂ₛg

    # Wave initial conditions
    u₀::Tu₀
    ∂ₓu₀::T∂ₓu₀
    ∂ᵧu₀::T∂ᵧu₀
    v₀::Tv₀
    ∂ₓv₀::T∂ₓv₀
    ∂ᵧv₀::T∂ᵧv₀

    # Thermal initial conditions
    θ₀::Tθ₀
    ∂ₓθ₀::T∂ₓθ₀
    ∂ᵧθ₀::T∂ᵧθ₀

    # Acoustic initial conditions
    z₀::Tz₀
    r₀::Tr₀

    # Source terms
    f₁::Tf₁
    f₂::Tf₂
    f₃::Tf₃

    # Analytical solutions (nothing for physical simulations)
    u::Tu
    v::Tv
    θ::Tθ
    z::Tz
    r::Tr
end

# ============================================================================
# Example 1
# ============================================================================
"""
    example1_manufactured(p::Float64=2.4) -> PDEInputData

Manufactured solution with oscillatory wave–acoustic coupling
``g(x, s) = (1 + e^{-x^2})(\\sin(s) + 2s)``:
```math
\\begin{alignat*}{2}
& u(x,y,t)   &&= (x^p - x)(y^p - 1)(4 + t^2), \\\\
& θ(x,y,t)   &&= \\sin(πx)\\sin(πy)\\,e^{-t}, \\\\
& z(x,t)     &&= \\sin(\\pi x) + (1+e^{-x^2})
                  \\left[\\frac{\\cos\\big(-2t(x^p-x)\\big)-1}{2(x^p-x)}
                  - 2t^2(x^p-x)\\right],
\\end{alignat*}
```
where the acoustic displacement ``z(x,t)`` is obtained by integrating
```math
\\frac{∂z}{∂t}(x,t) = -\\frac{∂u}{∂y}(x,y_{\\min},t) + g\\!\\left(x,\\,\\frac{∂u}{∂t}(x,y_{\\min},t)\\right).
```

# Arguments
- `p::Float64=2.4`: Smoothness parameter controlling solution regularity.

# Returns
`PDEInputData` with analytical solution for convergence study.
"""
function example1_manufactured(p::Float64 = 2.4)
    # Precompute exponent-related constants
    p_minus_1 = p - 1.0
    p_minus_2 = p - 2.0
    p_vs_p_minus_1 = p * p_minus_1

    # Physical parameters
    a = (0.0, 0.0)
    q₁ = q₂ = q₃ = q₄ = 1.0
    ymin = 0.0

    # Coefficient functions
    α = t -> 1.0 + exp(-t)

    β = s -> 1.0 + exp(-s * s)
    dβ = s -> -2.0 * s * exp(-s * s)

    f = function (s)
        s_abs = abs(s)
        s_abs3 = s_abs * s_abs * s_abs
        return s * s_abs3
    end
    df = function (s)
        s_abs = abs(s)
        s_abs3 = s_abs * s_abs * s_abs
        return 4.0 * s_abs3
    end

    g = (x, s) -> (1.0 + exp(-x * x)) * muladd(2.0, s, sin(s))
    ∂ₛg = (x, s) -> (1.0 + exp(-x * x)) * (2.0 + cos(s))

    # Analytical solutions
    u = function (x, y, t)
        xp = x^p
        yp = y^p
        return (xp - x) * (yp - 1.0) * (4.0 + t * t)
    end
    v = function (x, y, t)
        xp = x^p
        yp = y^p
        return (xp - x) * (yp - 1.0) * (2.0 * t)
    end

    θ = (x, y, t) -> sinpi(x) * sinpi(y) * exp(-t)

    z = function (x, t)
        xp = x^p
        xp_minus_x = xp - x
        two_t = 2.0 * t
        exp_term = 1.0 + exp(-x * x)
        return sinpi(x) +
               exp_term * ((cos(-xp_minus_x * two_t) - 1.0) /
                (2.0 * xp_minus_x) - xp_minus_x * (two_t * t))
    end
    r = function (x, t)
        xp = x^p
        xp_minus_x = xp - x
        exp_term = 1.0 + exp(-x * x)
        return exp_term * (sin(-xp_minus_x * (2.0 * t)) - xp_minus_x * (4.0 * t))
    end

    # Auxiliary functions for manufactured source terms
    ∂ₜₜu = function (x, y, t)
        xp = x^p
        yp = y^p
        return (xp - x) * (yp - 1.0) * 2.0
    end
    Δu = function (x, y, t)
        xp = x^p
        yp = y^p
        xp_minus_2 = x^p_minus_2
        yp_minus_2 = y^p_minus_2
        time_term = 4.0 + t * t
        return ((p_vs_p_minus_1 * xp_minus_2) * (yp - 1.0) +
                (xp - x) * (p_vs_p_minus_1 * yp_minus_2)) * time_term
    end
    vₓ = function (x, y, t)
        yp = y^p
        xp_minus_1 = x^p_minus_1
        return (p * xp_minus_1 - 1.0) * (yp - 1.0) * 2.0 * t
    end
    vᵧ = function (x, y, t)
        xp = x^p
        yp_minus_1 = y^p_minus_1
        return (xp - x) * (p * yp_minus_1) * 2.0 * t
    end
    θₓ = (x, y, t) -> π * cospi(x) * sinpi(y) * exp(-t)
    θᵧ = (x, y, t) -> π * sinpi(x) * cospi(y) * exp(-t)
    ∫θ = t -> 4.0 * exp(-t) / π^2
    ∂ₜₜz = function (x, t)
        xp = x^p
        xp_minus_x = xp - x
        exp_term = 1.0 + exp(-x * x)
        return -2.0 * exp_term * xp_minus_x * (2.0 + cos(2.0 * t * xp_minus_x))
    end

    # Manufactured source terms
    f₁ = (x, y, t) -> ∂ₜₜu(x, y, t) - α(t) * Δu(x, y, t) + f(u(x, y, t)) +
                      a[1] * θₓ(x, y, t) + a[2] * θᵧ(x, y, t)
    f₂ = (x, y, t) -> (-1.0 + 2.0 * pi * pi * β(∫θ(t))) * θ(x, y, t) +
                      a[1] * vₓ(x, y, t) + a[2] * vᵧ(x, y, t)
    f₃ = (x, t) -> q₁ * ∂ₜₜz(x, t) + q₂ * r(x, t) + q₃ * z(x, t) +
                   q₄ * v(x, ymin, t)

    # Initial conditions
    u₀ = function (x, y)
        xp = x^p
        yp = y^p
        return (xp - x) * (yp - 1.0) * 4.0
    end
    ∂ₓu₀ = function (x, y)
        yp = y^p
        xp_minus_1 = x^p_minus_1
        return (p * xp_minus_1 - 1.0) * (yp - 1.0) * 4.0
    end
    ∂ᵧu₀ = function (x, y)
        xp = x^p
        yp_minus_1 = y^p_minus_1
        return (xp - x) * (p * yp_minus_1) * 4.0
    end

    v₀ = (x, y) -> 0.0
    ∂ₓv₀ = (x, y) -> 0.0
    ∂ᵧv₀ = (x, y) -> 0.0

    θ₀ = (x, y) -> sinpi(x) * sinpi(y)
    ∂ₓθ₀ = (x, y) -> π * cospi(x) * sinpi(y)
    ∂ᵧθ₀ = (x, y) -> π * sinpi(x) * cospi(y)

    z₀ = x -> sinpi(x)
    r₀ = x -> 0.0

    return PDEInputData(
        (0.0, 0.0),
        (1.0, 1.0),
        a,
        q₁, q₂, q₃, q₄,
        α, β, dβ, f, df, g, ∂ₛg,
        u₀, ∂ₓu₀, ∂ᵧu₀,
        v₀, ∂ₓv₀, ∂ᵧv₀,
        θ₀, ∂ₓθ₀, ∂ᵧθ₀,
        z₀, r₀,
        f₁, f₂, f₃,
        u, v, θ, z, r
    )
end

"""
    example1_zero_source(p::Float64=2.4) -> PDEInputData

Same configuration as `example1_manufactured` but with f₁ = f₂ = f₃ = 0.
No analytical solution available.

# Arguments
- `p::Float64=2.4`: Smoothness parameter for initial conditions.

# Returns
`PDEInputData` with analytical solutions set to `nothing`.
"""
function example1_zero_source(p::Float64 = 2.4)
    # Precompute exponent-related constants
    p_minus_1 = p - 1.0

    # Physical parameters
    a = (0.0, 0.0)
    q₁ = q₂ = q₃ = q₄ = 1.0

    # Coefficient functions
    α = t -> 1.0 + exp(-t)

    β = s -> 1.0 + exp(-s * s)
    dβ = s -> -2.0 * s * exp(-s * s)

    f = function (s)
        s_abs = abs(s)
        s_abs3 = s_abs * s_abs * s_abs
        return s * s_abs3
    end
    df = function (s)
        s_abs = abs(s)
        s_abs3 = s_abs * s_abs * s_abs
        return 4.0 * s_abs3
    end

    g = (x, s) -> (1.0 + exp(-x * x)) * muladd(2.0, s, sin(s))
    ∂ₛg = (x, s) -> (1.0 + exp(-x * x)) * (2.0 + cos(s))

    # Zero source terms
    f₁ = (x, y, t) -> 0.0
    f₂ = (x, y, t) -> 0.0
    f₃ = (x, t) -> 0.0

    # Initial conditions
    u₀ = function (x, y)
        xp = x^p
        yp = y^p
        return (xp - x) * (yp - 1.0) * 4.0
    end
    ∂ₓu₀ = function (x, y)
        yp = y^p
        xp_minus_1 = x^p_minus_1
        return (p * xp_minus_1 - 1.0) * (yp - 1.0) * 4.0
    end
    ∂ᵧu₀ = function (x, y)
        xp = x^p
        yp_minus_1 = y^p_minus_1
        return (xp - x) * (p * yp_minus_1) * 4.0
    end

    v₀ = (x, y) -> 0.0
    ∂ₓv₀ = (x, y) -> 0.0
    ∂ᵧv₀ = (x, y) -> 0.0

    θ₀ = (x, y) -> sinpi(x) * sinpi(y)
    ∂ₓθ₀ = (x, y) -> π * cospi(x) * sinpi(y)
    ∂ᵧθ₀ = (x, y) -> π * sinpi(x) * cospi(y)

    z₀ = x -> sinpi(x)
    r₀ = x -> 0.0

    return PDEInputData(
        (0.0, 0.0),
        (1.0, 1.0),
        a,
        q₁, q₂, q₃, q₄,
        α, β, dβ, f, df, g, ∂ₛg,
        u₀, ∂ₓu₀, ∂ᵧu₀,
        v₀, ∂ₓv₀, ∂ᵧv₀,
        θ₀, ∂ₓθ₀, ∂ᵧθ₀,
        z₀, r₀,
        f₁, f₂, f₃,
        nothing, nothing, nothing, nothing, nothing
    )
end

# ============================================================================
# Example 2: Linear Coupling Test Case
# ============================================================================
"""
    example2_manufactured(p::Float64=2.4) -> PDEInputData

Manufactured solution with linear coupling ``g(x,s) = (1+e^{-x^2})s``:
```math
\\begin{alignat*}{2}
& u(x,y,t)   &&= (x^p - x)(y^p - 1)(4 + t^2), \\\\
& θ(x,y,t)   &&= \\sin(πx)\\sin(πy)\\,e^{-t}, \\\\
& z(x,t)     &&= \\sin(\\pi x) - (1+e^{-x^2})(x^p-x)t^2,
\\end{alignat*}
```
where the acoustic displacement ``z(x,t)`` is obtained by integrating
```math
\\frac{∂z}{∂t}(x,t) = -\\frac{∂u}{∂y}(x,y_{\\min},t) + g\\!\\left(x,\\,\\frac{∂u}{∂t}(x,y_{\\min},t)\\right).
```

# Arguments
- `p::Float64=2.4`: Smoothness parameter controlling solution regularity.

# Returns
`PDEInputData` with analytical solution for convergence study.
"""
function example2_manufactured(p::Float64 = 2.4)
    # Precompute exponent-related constants
    p_minus_1 = p - 1.0
    p_minus_2 = p - 2.0
    p_vs_p_minus_1 = p * p_minus_1

    # Physical parameters
    a = (0.0, 0.0)
    q₁ = q₂ = q₃ = q₄ = 1.0
    ymin = 0.0

    # Coefficient functions
    α = t -> 1.0 + exp(-t)

    β = s -> 1.0 + exp(-s * s)
    dβ = s -> -2.0 * s * exp(-s * s)

    f = function (s)
        s_abs = abs(s)
        s_abs3 = s_abs * s_abs * s_abs
        return s * s_abs3
    end
    df = function (s)
        s_abs = abs(s)
        s_abs3 = s_abs * s_abs * s_abs
        return 4.0 * s_abs3
    end

    g = (x, s) -> (1.0 + exp(-x * x)) * s
    ∂ₛg = (x, s) -> 1.0 + exp(-x * x)

    # Analytical solutions
    u = function (x, y, t)
        xp = x^p
        yp = y^p
        return (xp - x) * (yp - 1.0) * (4.0 + t * t)
    end
    v = function (x, y, t)
        xp = x^p
        yp = y^p
        return (xp - x) * (yp - 1.0) * (2.0 * t)
    end

    θ = (x, y, t) -> sinpi(x) * sinpi(y) * exp(-t)

    z = function (x, t)
        xp = x^p
        exp_term = 1.0 + exp(-x * x)
        return sinpi(x) - exp_term * (xp - x) * (t * t)
    end

    r = function (x, t)
        xp = x^p
        exp_term = 1.0 + exp(-x * x)
        return -exp_term * (xp - x) * (2.0 * t)
    end

    # Auxiliary functions for manufactured source terms
    ∂ₜₜu = function (x, y, t)
        xp = x^p
        yp = y^p
        return (xp - x) * (yp - 1.0) * 2.0
    end
    Δu = function (x, y, t)
        xp = x^p
        yp = y^p
        xp_minus_2 = x^p_minus_2
        yp_minus_2 = y^p_minus_2
        time_term = 4.0 + t * t
        return ((p_vs_p_minus_1 * xp_minus_2) * (yp - 1.0) +
                (xp - x) * (p_vs_p_minus_1 * yp_minus_2)) * time_term
    end
    vₓ = function (x, y, t)
        yp = y^p
        xp_minus_1 = x^p_minus_1
        return (p * xp_minus_1 - 1.0) * (yp - 1.0) * 2.0 * t
    end
    vᵧ = function (x, y, t)
        xp = x^p
        yp_minus_1 = y^p_minus_1
        return (xp - x) * (p * yp_minus_1) * 2.0 * t
    end
    θₓ = (x, y, t) -> π * cospi(x) * sinpi(y) * exp(-t)
    θᵧ = (x, y, t) -> π * sinpi(x) * cospi(y) * exp(-t)
    ∫θ = t -> 4.0 * exp(-t) / π^2
    ∂ₜₜz = @inline function (x, t)
        xp = x^p
        exp_term = 1.0 + exp(-x * x)
        return -2.0 * exp_term * (xp - x)
    end

    # Manufactured source terms
    f₁ = (x, y, t) -> ∂ₜₜu(x, y, t) - α(t) * Δu(x, y, t) + f(u(x, y, t)) +
                      a[1] * θₓ(x, y, t) + a[2] * θᵧ(x, y, t)
    f₂ = (x, y, t) -> (-1.0 + 2.0 * pi * pi * β(∫θ(t))) * θ(x, y, t) +
                      a[1] * vₓ(x, y, t) + a[2] * vᵧ(x, y, t)
    f₃ = (x, t) -> q₁ * ∂ₜₜz(x, t) + q₂ * r(x, t) + q₃ * z(x, t) +
                   q₄ * v(x, ymin, t)

    # Initial conditions
    u₀ = function (x, y)
        xp = x^p
        yp = y^p
        return (xp - x) * (yp - 1.0) * 4.0
    end
    ∂ₓu₀ = function (x, y)
        yp = y^p
        xp_minus_1 = x^p_minus_1
        return (p * xp_minus_1 - 1.0) * (yp - 1.0) * 4.0
    end
    ∂ᵧu₀ = function (x, y)
        xp = x^p
        yp_minus_1 = y^p_minus_1
        return (xp - x) * (p * yp_minus_1) * 4.0
    end

    v₀ = (x, y) -> 0.0
    ∂ₓv₀ = (x, y) -> 0.0
    ∂ᵧv₀ = (x, y) -> 0.0

    θ₀ = (x, y) -> sinpi(x) * sinpi(y)
    ∂ₓθ₀ = (x, y) -> π * cospi(x) * sinpi(y)
    ∂ᵧθ₀ = (x, y) -> π * sinpi(x) * cospi(y)

    z₀ = x -> sinpi(x)
    r₀ = x -> 0.0

    return PDEInputData(
        (0.0, 0.0),
        (1.0, 1.0),
        a,
        q₁, q₂, q₃, q₄,
        α, β, dβ, f, df, g, ∂ₛg,
        u₀, ∂ₓu₀, ∂ᵧu₀,
        v₀, ∂ₓv₀, ∂ᵧv₀,
        θ₀, ∂ₓθ₀, ∂ᵧθ₀,
        z₀, r₀,
        f₁, f₂, f₃,
        u, v, θ, z, r
    )
end

"""
    example2_zero_source(p::Float64=2.4) -> PDEInputData

Same configuration as `example2_manufactured` but with f₁ = f₂ = f₃ = 0.
No analytical solution available.

# Arguments
- `p::Float64=2.4`: Smoothness parameter for initial conditions.

# Returns
`PDEInputData` with analytical solutions set to `nothing`.
"""
function example2_zero_source(p::Float64 = 2.4)
    # Precompute exponent-related constants
    p_minus_1 = p - 1.0

    # Physical parameters
    a = (0.0, 0.0)
    q₁ = q₂ = q₃ = q₄ = 1.0

    # Coefficient functions
    α = t -> 1.0 + exp(-t)

    β = s -> 1.0 + exp(-s * s)
    dβ = s -> -2.0 * s * exp(-s * s)

    f = function (s)
        s_abs = abs(s)
        s_abs3 = s_abs * s_abs * s_abs
        return s * s_abs3
    end
    df = function (s)
        s_abs = abs(s)
        s_abs3 = s_abs * s_abs * s_abs
        return 4.0 * s_abs3
    end

    g = (x, s) -> (1.0 + exp(-x * x)) * s
    ∂ₛg = (x, s) -> 1.0 + exp(-x * x)

    # Manufactured source terms
    f₁ = (x, y, t) -> 0.0
    f₂ = (x, y, t) -> 0.0
    f₃ = (x, t) -> 0.0

    # Initial conditions
    u₀ = function (x, y)
        xp = x^p
        yp = y^p
        return (xp - x) * (yp - 1.0) * 4.0
    end
    ∂ₓu₀ = function (x, y)
        yp = y^p
        xp_minus_1 = x^p_minus_1
        return (p * xp_minus_1 - 1.0) * (yp - 1.0) * 4.0
    end
    ∂ᵧu₀ = function (x, y)
        xp = x^p
        yp_minus_1 = y^p_minus_1
        return (xp - x) * (p * yp_minus_1) * 4.0
    end

    v₀ = (x, y) -> 0.0
    ∂ₓv₀ = (x, y) -> 0.0
    ∂ᵧv₀ = (x, y) -> 0.0

    θ₀ = (x, y) -> sinpi(x) * sinpi(y)
    ∂ₓθ₀ = (x, y) -> π * cospi(x) * sinpi(y)
    ∂ᵧθ₀ = (x, y) -> π * sinpi(x) * cospi(y)

    z₀ = x -> sinpi(x)
    r₀ = x -> 0.0

    return PDEInputData(
        (0.0, 0.0),
        (1.0, 1.0),
        a,
        q₁, q₂, q₃, q₄,
        α, β, dβ, f, df, g, ∂ₛg,
        u₀, ∂ₓu₀, ∂ᵧu₀,
        v₀, ∂ₓv₀, ∂ᵧv₀,
        θ₀, ∂ₓθ₀, ∂ᵧθ₀,
        z₀, r₀,
        f₁, f₂, f₃,
        nothing, nothing, nothing, nothing, nothing
    )
end