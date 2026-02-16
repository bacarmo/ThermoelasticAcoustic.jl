## Strong Formulation
We seek functions ``u``, ``\theta``, and ``z`` satisfying the following system of equations:
```math
\begin{align*}
& u^{\prime\prime}(x,t) 
- \alpha(t)\Delta u(x,t) 
+ f\big(u(x,t)\big) 
= f_1(x,t),
\quad(x,t)\in\Omega\times(0,+\infty),
\\
& \theta^{\prime}(x,t) 
- \beta\Big(\int_\Omega\theta(t)dx\Big)\Delta\theta(x,t)
+ (\mathbf{a}\cdot\nabla)u^\prime(x,t)
= f_2(x,t),
\quad(x,t)\in\Omega\times(0,+\infty),
\\[5pt]
& q_1 z^{\prime\prime}(x,t)
+ q_2 z^\prime(x,t)
+ q_3 z(x,t)
+ q_4 u^\prime(x,t)
= f_3(x,t),
\quad(x,t)\in\Gamma_1\times(0,+\infty),
\\
& \frac{\partial u}{\partial\nu}(x,t)
= z^\prime(x,t)
- g\big(x,u^\prime(x,t)\big),
\quad(x,t)\in\Gamma_1\times(0,+\infty),
\\[5pt]
& u(x,t) = 0,\quad (x,t)\in\Gamma_0\times(0,+\infty),
\\[5pt]
& \theta(x,t) = 0,\quad (x,t)\in\Gamma\times(0,+\infty),
\end{align*}
```
with initial conditions
```math
\begin{align*}
& u(x,0) = u_0(x),\quad u^\prime(x,0)=v_0(x),\quad x\in\Omega,
\\[5pt]
& \theta_0(x,0)=\theta_0(x),\quad x\in\Omega,
\\
& z(x,0) = z_0(x),\quad
  z^\prime(x,0) = r_0(x) \equiv \frac{\partial u_0}{\partial\nu}(x) + g\big(x,v_0(x)\big),\quad x\in\Gamma_1,
\end{align*}
```
where ``\Omega`` is a bounded open subset of ``\mathbb{R}^n``, ``n\geq 2``, with smooth boundary ``\Gamma=\Gamma_0\cup\Gamma_1`` and disjoint ``\Gamma_0``, ``\Gamma_1``.

The existence, uniqueness, and asymptotic behavior of global solutions to this problem were investigated in [Braz e Silva et al. (2017)](https://link.springer.com/article/10.1007/s40314-015-0236-1). 
When the source terms $f_1$, $f_2$, and $f_3$ are identically zero and the nonlinear term ``f`` takes the form $f(s) = s|s|^\rho$, the authors establish that the total energy 
```math
\begin{align*}
E(t) 
=&  
  \frac{1}{2}\Big[
  \|u^\prime(t)\|^2
+ \|\theta(t)\|^2
+ \frac{1}{\rho+2}\|u(t)\|_{L^{\rho+2}(\Omega)}^{\rho+2}
+
  \alpha(t)\big(
  \|\nabla u(t)\|^2
+ \|\sqrt{q_1}z^\prime(t)\|_{\Gamma_1}^2
+ \|\sqrt{q_3}z(t)\|_{\Gamma_1}^2
\big)
\Big]
\end{align*}
```
associated with the strong global solution exhibits exponential decay.