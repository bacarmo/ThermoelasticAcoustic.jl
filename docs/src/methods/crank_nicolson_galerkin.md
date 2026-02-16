## Fully discrete problem definition
The numerical scheme present is based on the Crank-Nicolson Galerkin method, which consists of finding ``U^n, V^n \in \mathcal{V}_{m_1}``, ``\Theta^n\in \mathcal{V}_{m_2}``, and ``Z^n, R^n \in \mathcal{V}_{m_3}`` such that
```math
\begin{align*}
& \big(\varphi,\bar{\partial}V^n\big)
+ \alpha(t_{n-\frac{1}{2}})\Big[
    \big(\nabla\varphi,\nabla\widehat{U}^n\big)
  - \big(\varphi,\widehat{R}^n)_{\Gamma_1}
  + \big(\varphi,g(\widehat{V}^n)\big)_{\Gamma_1}\Big]
+ \big(\varphi,f(\widehat{U}^n)\big) 
= \big(\varphi,f_1(t_{n-\frac{1}{2}}\big),
\quad\forall\varphi\in \mathcal{V}_{m_1},
\\
& \big(\psi,\bar{\partial}\Theta^n\big)
+ \beta\Big(\int_\Omega\hat{\Theta}^ndx\Big)\big(\nabla\psi,\nabla\widehat{\Theta}^n\big)
+ \big(\psi,(\mathbf{a}\cdot\nabla)\widehat{V}^n\big)
= \big(\psi,f_2(t_{n-\frac{1}{2}})\big),
\quad\forall\psi\in \mathcal{V}_{m_2},
\\[5pt]
& \big(\phi,q_1\bar{\partial}R^n
+ q_2\widehat{R}^n
+ q_3\widehat{Z}^n
+ q_4\widehat{V}^n\big)_{\Gamma_1}
= \big(\phi,f_3(t_{n-\frac{1}{2}})\big)_{\Gamma_1},
\quad\forall\phi\in \mathcal{V}_{m_3},
\\[5pt]
& 
\bar{\partial}U^n = \widehat{V}^n,\quad 
\bar{\partial}Z^n = \widehat{R}^n,
\end{align*}
```
with ``U^0, V^0 \in \mathcal{V}_{m_1}``, ``\Theta^0 \in \mathcal{V}_{m_2}``, and ``Z^0, R^0 \in \mathcal{V}_{m_3}`` given as approximations of the initial solutions ``u_0,\, v_0,\, \theta_0,\, z_0``, and ``r_0``.
We define these approximations by solving the following projections problems:
find  ``U^0,\, V^0 \in \mathcal{V}_{m_1}``, ``\Theta^0 \in \mathcal{V}_{m_2}``, and ``Z^0,\,R^0 \in \mathcal{V}_{m_3}`` such that
```math
\begin{align*}
& \big(\nabla\varphi,\nabla U^0\big) 
= \big(\nabla\varphi,\nabla u_0\big), \;\;
\forall\varphi\in \mathcal{V}_{m_1};
\\
& \big(\nabla\varphi,\nabla V^0\big) 
= \big(\nabla\varphi,\nabla v_0\big),\;\;
\forall\varphi\in \mathcal{V}_{m_1};
\\
& \big(\nabla\psi,\nabla\Theta^0\big) 
= \big(\nabla\psi,\nabla\theta_0\big),\;\;
\forall\psi\in \mathcal{V}_{m_2};
\\
& \big(\phi,Z^0\big)_{\Gamma_1} 
= \big(\phi,z_0\big)_{\Gamma_1},\;\;
\forall\phi\in \mathcal{V}_{m_3};
\\
& \big(\phi,R^0\big)_{\Gamma_1} 
= \big(\phi,r_0\big)_{\Gamma_1},\;\;
\forall\phi\in \mathcal{V}_{m_3}.
\end{align*}
```

!!! details "Notation"
    - ``\mathcal{V}_{m_1} \subset H_{\Gamma_0}^1(\Omega)``: Subspace of dimension ``m_1`` with basis ``\{\varphi_j\}_{j=1}^{m_1}``.
    - ``\mathcal{V}_{m_2} \subset H_{0}^1(\Omega)``: Subspace of dimension ``m_2`` with basis ``\{\psi_j\}_{j=1}^{m_2}``.
    - ``\mathcal{V}_{m_3} = \mathcal{V}_{m_1}|_{\Gamma_1}``: Subspace of dimension ``m_3`` with basis ``\{\phi_j\}_{j=1}^{m_3}``.
    - ``\displaystyle w^n:=w(t_n),\quad \bar{\partial}w^n := \frac{w^n - w^{n-1}}{\tau}\approx w^\prime(t_{n-\frac{1}{2}}),\quad \widehat{w}^n := \frac{w^n + w^{n-1}}{2}\approx w(t_{n-\frac{1}{2}}),``
    where ``\tau`` denotes the time step, ``t_n = n\tau`` the discrete times, ``t_{n-\frac{1}{2}}`` the midpoint of ``[t_{n-1},t_{n}]``, and ``w`` an arbitrary time-dependent function.