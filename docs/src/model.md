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
& q_1(x) z^{\prime\prime}(x,t)
+ q_2(x) z^\prime(x,t)
+ q_3(x) z(x,t)
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
When the source terms $f_1$, $f_2$, and $f_3$ are identically zero and the nonlinear term ``f`` takes the form $f(s) = \lambda s|s|^\rho$, the authors establish that the total energy 
```math
\begin{align*}
E(t) 
=&  
  \frac{1}{2}\Big[
  \|u^\prime(t)\|^2
+ \|\theta(t)\|^2
+ \frac{\lambda}{\rho+2}\|u(t)\|_{L^{\rho+2}(\Omega)}^{\rho+2}
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

!!! details "Assumptions for existence and uniqueness"
    * ``\Omega \subset \mathbb{R}^n``, ``n\geq 2``, is a bounded, open, and connected domain, such that its boundary ``\Gamma`` is locally situated on one side.
    * ``\Gamma`` is a ``C^2``-class boundary, partitioned into two disjoint, connected components ``\Gamma_0`` and ``\Gamma_1``, each having positive measure.
    * ``\alpha \in L^\infty_{loc}(\mathbb{R}^+)``, ``\alpha'\in L^1(\mathbb{R}^+)``, and ``\alpha(t) \geq \alpha_0 > 0,\;\forall t \geq 0``.
    * The constant ``\mathbf{a}\in\mathbb{R}^n``, and ``(\mathbf{a}\cdot\nabla)`` is the operator ``\sum_{i=1}^n \mathbf{a}_i\frac{\partial}{\partial x_i}``.
    * ``f(s)= \lambda |s|^\rho s``, with ``\lambda > 0``. For ``n = 2``, ``\rho > 0``; for ``n \geq 3``, ``\frac{1}{n} \leq \rho \leq \frac{2}{n-2}``.
    * ``\beta^\prime\in L_{loc}^\infty(\mathbb{R})`` and ``\beta(s)\geq\beta_0>0``.
    * ``x\mapsto g(x,s)`` is continuous ``\forall s\in\mathbb{R}`` fixed; ``|g(x,s)-g(x,r)|\leq g_0|s-r|``, for all ``x\in\Gamma_1,\,s,r\in\mathbb{R}``; ``\;[g(x,s)-g(x,r)](s-r)\geq g_1(s-r)^2``; ``\;g(x,0)=0``;  ``\;g_0,g_1>0``.
    * ``q_1,q_2,q_3\in C^0(\Gamma_1)``, with ``q_1(x),q_3(x)>0``, and ``q_2(x)\geq 0``, for all ``x\in\Gamma_1``.
    * ``u_0\in H_{\Gamma_0}^1(\Omega)\cap H_\Delta(\Omega)``, ``\;v_0\in H_{\Gamma_0}^1(\Omega)``, ``\;\theta_0\in H_0^1(\Omega)``, ``\;z_0\in L^2(\Gamma_1)``.

    !!! note "Beware!"
        The source functions ``f_1``, ``f_2``, and ``f_3``, and the constant term ``q_4`` in the model are not present in the original formulation analyzed in [Braz e Silva et al. (2017)](https://link.springer.com/article/10.1007/s40314-015-0236-1).
        The inclusion of the sources terms serves exclusively to construct an exact solution for the model, which is required to verify the convergence order of the corresponding approximate solution.
        The constant ``q_4`` is introduced to facilitate validation: setting ``q_4=0`` decouples the acoustic solution ``z`` from the remaining equations, allowing independent verification.

!!! details "Supplementary assumptions for energy decay"
    * ``\Gamma_0`` and ``\Gamma_1`` have a special geometry: they are closed, connected and disjoints sets (``\overline{\Gamma_0}\cap\overline{\Gamma_1}=\empty``) with positive measure. Furthermore, there exists a function ``m(x) = x - x_0``, for some fixed ``x_0 \in \mathbb{R}^n``, such that
    ```math
    \Gamma_0 = \{ x \in \Gamma ; m(x) \cdot \nu(x) \leq 0 \}, \quad
    \Gamma_1 = \{ x \in \Gamma ; m(x) \cdot \nu(x) > 0 \}.
    ```
    * ``q_2(x)>0``, ``\;\forall x\in\Gamma_1``.
    * ``g`` has the special form:
    ```math
    \begin{aligned}
    & g(x,s) = [m(x)\cdot\nu(x)]\widetilde{g}(s),\quad\forall x\in\Gamma_1,
    \\[5pt]
    & |\widetilde{g}(s)| \leq \widetilde{g}_0|s|,\quad \forall s\in\mathbb{R},
    \\[5pt]
    & [\widetilde{g}(s)-\widetilde{g}(r)](s-r)\geq\widetilde{g}_1(s-r)^2,\quad\forall s,r\in\mathbb{R},\hbox{ for some }\widetilde{g}_0,\widetilde{g}_1>0.
    \end{aligned}
    ```
    * ``\alpha_0\leq\alpha(t)\leq\alpha_1\; \hbox{ a.e. in }\mathbb{R}^+``, and ``|\alpha^\prime(t)|\leq\alpha_0\epsilon,\;\forall t\in\mathbb{R}``, where 
    ```math
    \begin{aligned}
    & \epsilon
    = \min\Big\{
    \frac{\beta_0}{\kappa_2+\frac{1}{2\lambda_1}},
    \frac{\widetilde{g}_1}{\kappa_3},
    \frac{\overline{q}_2}{\kappa_4\widehat{q}_1}
    \Big\};
    \quad
    \overline{q}_i = \underset{x\in\overline{\Gamma}_1}{\min\,}q_i(x);
    \quad
    \widehat{q}_i = \underset{x\in\overline{\Gamma}_1}{\max\,}q_i(x);
    \\
    & \lambda_1
      \hbox{ is the first eigenvalue of the Laplace operator, } \lambda_1\|v\|^2\leq\|\nabla v\|^2,\,\forall v\in H_{\Gamma_0}^1(\Omega);
    \\
    & \kappa_2 
    = \frac{|\mathbf{a}|_\infty^2}{\alpha_0}
    \Big(
      16\widehat{m}^2 + \frac{2n}{\lambda_1}(n-\frac{1}{2})^2
    \Big);
    \;\; \widehat{m} = \underset{x\in\overline{\Omega}}{\max\,}|m(x)|_{\mathbb{R}^n};
    \;\; |\mathbf{a}|_\infty = \underset{1\leq i\leq n}{\max\,}|\mathbf{a}_i|;
    \\
    & \kappa_3
    = 2\widehat{m}^2\widetilde{g}_0^2 
    + \frac{1}{\alpha_0}
    + 2C_1(n-\frac{1}{2})^2\widetilde{g}_0^2
    + \alpha_0
    + \frac{3}{2}\alpha_1;
    \\
    & \kappa_4
    = \frac{2\widehat{m}^2 + 2C_1(n-\frac{1}{2})^2}{\zeta\overline{q}_1}
    + (\alpha_0 + \frac{3}{2}\alpha_1)\frac{\widehat{q}_2}{\overline{q}_1}
    + \alpha_0\epsilon;
    \;\;
    \zeta=\underset{x\in\overline{\Gamma}_1}{\min}\{m(x)\cdot\nu(x)\};
    \\
    &
    C_1\in\mathbb{R}^+\hbox{ such that }
    \|\sqrt{m\cdot\nu}\,v\|_{\Gamma_1}^2\leq\frac{C_1}{2}\|\nabla v\|^2,\;\forall v\in H_{\Gamma_0}^1(\Omega) \hbox{(mapping trace continuity)}.
    \end{aligned}
    ```
    !!! warning "Beware!"
        The definition of $\epsilon$ depends on $\kappa_4$, which in turn depends on $\epsilon$. How to solve this circular dependency?


!!! details "Notation: Functional spaces, inner products, and norms"
    We consider the function spaces
    ```math
        \begin{aligned}
        H^1(\Omega) &= \{ v \in L^2(\Omega) : \nabla v \in (L^2(\Omega))^n \},\\
        H_0^1(\Omega) &= \{ v \in H^1(\Omega) : v|_{\Gamma} = 0 \}, \\
        H_{\Gamma_0}^1(\Omega) &= \{ v \in H^1(\Omega) : v|_{\Gamma_0} = 0 \}, \\
        H_\Delta(\Omega) &= \{ v \in H^1(\Omega) : \Delta v \in L^2(\Omega) \}.
        \end{aligned}
    ```
    We denote the inner products and norms in $L^2(\Omega)$ and $L^2(\Gamma_1)$ as follows:
    ```math
    (\cdot, \cdot),\quad
    \|\cdot\|,
    \qquad\text{and}\qquad
    (\cdot, \cdot)_{\Gamma_1},\quad
    \|\cdot\|_{\Gamma_1}.
    ```

## Weak Formulation
We seek functions ``u(t)\in H_{\Gamma_0}^1(\Omega)``, ``\theta(t)\in H_0^1(\Omega)``, and ``z(t)\in L^2(\Gamma_1)`` such that
```math
\begin{align*}
& \big(\varphi,u^{\prime\prime}(t)\big)
+ \alpha(t)\Big[
    \big(\nabla\varphi,\nabla u(t)\big)
  - \big(\varphi,z^\prime(t))_{\Gamma_1}
  + \big(\varphi,g\big(u^\prime(t)\big)\big)_{\Gamma_1}\Big]
+ \big(\varphi,f\big(u(t)\big)\big) 
= \big(\varphi,f_1(t)\big),
\quad\forall\varphi\in H_{\Gamma_0}^1(\Omega),
\\
& \big(\psi,\theta^{\prime}(t)\big)
+ \beta\Big(\int_\Omega\theta(t)dx\Big)\big(\nabla\psi,\nabla\theta(t)\big)
+ \big(\psi,(\mathbf{a}\cdot\nabla)u^\prime(t)\big)
= \big(\psi,f_2(t)\big),
\quad\forall\psi\in H_0^1(\Omega),
\\[5pt]
& \big(\phi,q_1z^{\prime\prime}(t)
+ q_2z^\prime(t)
+ q_3z(t)
+ q_4u^\prime(t)\big)_{\Gamma_1}
= \big(\phi,f_3(t)\big)_{\Gamma_1},
\quad\forall\phi\in L^2(\Gamma_1),
\end{align*}
```
with 
``u(0)=u_0``, 
``u^\prime(0)=v_0``, 
``\theta(0)=\theta_0``,
``z(0)=z_0``, and 
``z^\prime(0) = r_0 \equiv \frac{\partial u_0}{\partial\nu} + g(v_0)``. 

By introducing the auxiliary variables ``v(t)=u^\prime(t)`` and ``r(t)=z^\prime(t)``, we obtain the equivalent first-order system: find functions ``u(t),v(t)\in H_{\Gamma_0}^1(\Omega)``, ``\theta(t)\in H_0^1(\Omega)``,  and ``z(t),r(t)\in L^2(\Gamma_1)`` such that
```math
\begin{align*}
& \big(\varphi,v^{\prime}(t)\big)
+ \alpha(t)\Big[
    \big(\nabla\varphi,\nabla u(t)\big)
  - \big(\varphi,r(t)\big)_{\Gamma_1}
  + \big(\varphi,g\big(v(t)\big)\big)_{\Gamma_1}\Big]
+ \big(\varphi,f\big(u(t)\big)\big) 
= \big(\varphi,f_1(t)\big),
\quad\forall\varphi\in H_{\Gamma_0}^1(\Omega),
\\
& \big(\psi,\theta^{\prime}(t)\big)
+ \beta\Big(\int_\Omega\theta(t)dx\Big)\big(\nabla\psi,\nabla\theta(t)\big)
+ \big(\psi,(\mathbf{a}\cdot\nabla)u^\prime(t)\big)
= \big(\psi,f_2(t)\big),
\quad\forall\psi\in H_0^1(\Omega),
\\[5pt]
& \big(\phi,q_1r^{\prime}(t)
+ q_2r(t)
+ q_3z(t)
+ q_4v(t)\big)_{\Gamma_1}
= \big(\phi,f_3(t)\big)_{\Gamma_1},
\quad\forall\phi\in L^2(\Gamma_1),
\\[5pt]
& u^\prime(t)=v(t),\quad z^\prime(t)=r(t),
\end{align*}
```
with initial conditions
``u(0)=u_0``, 
``v(0)=v_0``,
``\theta(0)=\theta_0``, 
``z(0)=z_0``, and 
``r(0) = r_0 \equiv \frac{\partial u_0}{\partial\nu} + g(v_0)``.