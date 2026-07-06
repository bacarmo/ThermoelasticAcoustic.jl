# Scheme 1
The numerical scheme present is based on the Crank-Nicolson Galerkin method, which consists of finding ``U^n, V^n \in \mathcal{V}_{m_1}``, ``\Theta^n\in \mathcal{V}_{m_2}``, and ``Z^n, R^n \in \mathcal{V}_{m_3}`` such that
```math
\begin{align*}
& \big(\varphi,\bar{\partial}V^n\big)
+ \alpha^{n-\frac{1}{2}}\Big[
    \big(\nabla\varphi,\nabla\widehat{U}^n\big)
  - \big(\varphi,\widehat{R}^n)_{\Gamma_1}
  + \big(\varphi,g(\widehat{V}^n)\big)_{\Gamma_1}\Big]
+ \big(\varphi,f(\widehat{U}^n)\big) 
+ \big(\varphi,(\mathbf{a}\cdot\nabla)\widehat\Theta^n\big) 
= \big(\varphi,f_1^{n-\frac{1}{2}}\big),
\quad\forall\varphi\in \mathcal{V}_{m_1},
\\
& \big(\psi,\bar{\partial}\Theta^n\big)
+ \beta\Big(\int_\Omega\widehat{\Theta}^ndx\Big)\big(\nabla\psi,\nabla\widehat{\Theta}^n\big)
+ \big(\psi,(\mathbf{a}\cdot\nabla)\widehat{V}^n\big)
= \big(\psi,f_2^{n-\frac{1}{2}}\big),
\quad\forall\psi\in \mathcal{V}_{m_2},
\\[5pt]
& \big(\phi,q_1\bar{\partial}R^n
+ q_2\widehat{R}^n
+ q_3\widehat{Z}^n
+ q_4\widehat{V}^n\big)_{\Gamma_1}
= \big(\phi,f_3^{n-\frac{1}{2}}\big)_{\Gamma_1},
\quad\forall\phi\in \mathcal{V}_{m_3},
\\[5pt]
& 
\bar{\partial}U^n = \widehat{V}^n,\quad 
\bar{\partial}Z^n = \widehat{R}^n,
\end{align*}
```
with ``U^0, V^0 \in \mathcal{V}_{m_1}``, ``\Theta^0 \in \mathcal{V}_{m_2}``, and ``Z^0, R^0 \in \mathcal{V}_{m_3}`` given as approximations of the initial solutions ``u_0``, ``\, v_0``, ``\,\theta_0``, ``\, z_0``, and ``r_0``.
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
    - ``\mathcal{V}_{m_3} = \mathcal{V}_{m_1}|_{\Gamma_1}``: Subspace of dimension ``m_3`` with basis ``\{\phi_j\}_{j=1}^{m_3}``, where the basis functions of ``\mathcal{V}_{m_1}``
      are indexed so that
    ```math
        \phi_j = \varphi_j|_{\Gamma_1}, \quad j = 1, \dots, m_3.
    ```
    - For a sufficiently regular time-dependent function ``w``, ``\tau`` the time step, and
      ``t_n = n\tau``:
      ```math
          w^n := w(t_n), \qquad
          \bar{\partial}w^n := \frac{w^n - w^{n-1}}{\tau} \approx w^\prime(t_{n-\frac{1}{2}}), \qquad
          \widehat{w}^n := \frac{w^n + w^{n-1}}{2} \approx w(t_{n-\frac{1}{2}}),
      ```
      where ``t_{n-\frac{1}{2}}`` denotes the midpoint of ``[t_{n-1}, t_n]``.

## Matrix Formulation
Representing the approximate solutions in terms of the basis functions,
```math
U^n      = \sum_{j=1}^{m_1} d_j^n\varphi_j,\;
V^n      = \sum_{j=1}^{m_1} v_j^n\varphi_j,\;
\Theta^n = \sum_{j=1}^{m_2} c_j^n\psi_j,\;
Z^n      = \sum_{j=1}^{m_3} z_j^n\phi_j,\;
R^n      = \sum_{j=1}^{m_3} r_j^n\phi_j,
```
and choosing test functions ``\varphi=\varphi_i`` for ``i = 1, \ldots, m_1``, ``\psi=\psi_i`` for ``i = 1, \ldots, m_2`` and ``\phi=\phi_i`` for ``i = 1, \ldots, m_3``, we obtain the system
```math
\begin{align*}
& M^{m_1\times m_1}\bar{\partial}v^n
+ \alpha^{n-\frac{1}{2}}\Big[  
    K^{m_1\times m_1}\widehat{d}^n
  - M^{m_1\times m_3}\widehat{r}^n
  + G^{m_1}(\widehat{v}^n)\Big]
+ F^{m_1}(\widehat{d}^n)
+ A^{m_1\times m_2}\widehat{c}^n
= \mathcal{F}^{m_1}(f_1^{n-\frac{1}{2}}),
\\
& M^{m_2\times m_2}\bar{\partial}c^n
+ \beta(\mathbf{b}\cdot\widehat{c}^n)K^{m_2\times m_2}\widehat{c}^n
+ A^{m_2\times m_1}\widehat{v}^n
= \mathcal{F}^{m_2}(f_2^{n-\frac{1}{2}}),
\\[3pt]
& M^{m_3\times m_3}\big[
  q_1\bar{\partial}r^n
+ q_2\widehat{r}^n
+ q_3\widehat{z}^n\big]
+ q_4M^{m_3\times m_1}\widehat{v}^n
= \mathcal{F}^{m_3}(f_3^{n-\frac{1}{2}}),
\\[5pt]
& \bar{\partial}d^n = \widehat{v}^n,\quad
\bar{\partial}z^n = \widehat{r}^n,
\end{align*}
```
with initial solution given by
```math
K^{m_1\times m_1} d^0 = \kappa^{m_1}(u_0),\quad
K^{m_1\times m_1} v^0 = \kappa^{m_1}(v_0),\quad
K^{m_2\times m_2} c^0 = \kappa^{m_2}(\theta_0),\quad
M^{m_3\times m_3} z^0 = \mathcal{F}^{m_3}(z_0),\quad
M^{m_3\times m_3} r^0 = \mathcal{F}^{m_3}(r_0).
```

!!! details "Matrix and vector definitions"
    ```math
    \begin{aligned}
    &
    M^{m_1\times m_1}_{i,j} = (\varphi_i,\varphi_j),\quad 
    M^{m_1\times m_3}_{i,j} = (\varphi_i,\phi_j)_{\Gamma_1},\quad 
    M^{m_3\times m_1}_{i,j} = (\phi_i,\varphi_j)_{\Gamma_1},\quad 
    M^{m_3\times m_3}_{i,j} = (\phi_i,\phi_j)_{\Gamma_1},
    \\[5pt]
    &
    M^{m_2\times m_2}_{i,j} = (\psi_i,\psi_j),\quad 
    K^{m_1\times m_1}_{i,j} = (\nabla\varphi_i,\nabla\varphi_j),\quad
    K^{m_2\times m_2}_{i,j} = (\nabla\psi_i,\nabla\psi_j),
    \\[5pt]
    &
    A^{m_1\times m_2}_{i,j} = (\varphi_i,(\mathbf{a}\cdot\nabla)\psi_j),\quad
    A^{m_2\times m_1}_{i,j} = (\psi_i,(\mathbf{a}\cdot\nabla)\varphi_j),
    \\[5pt]
    &
    \mathcal{F}_i^{m_1}(w) = \big(\varphi_i,w\big),\quad
    \mathcal{F}_i^{m_2}(w) = \big(\psi_i,w\big),\quad
    \mathcal{F}_i^{m_3}(w) = \big(\phi_i,w\big)_{\Gamma_1},
    \\[5pt]
    &
    \kappa^{m_1}(w) = \big(\nabla\varphi_i,\nabla w\big),\quad
    \kappa^{m_2}(w) = \big(\nabla\psi_i,\nabla w\big),
    \\[5pt]
    &
    G^{m_1}_i(\widehat{v}^n) 
    = \big(\varphi_i,g(\widehat{V}^n)\big)_{\Gamma_1}
    \equiv\int_{\Gamma_1}\varphi_i(x)g\Big(x,\sum_{\ell=1}^{m_1}\widehat{v}_\ell^n\varphi_\ell(x)\Big)d\Gamma,
    \\[5pt]
    &
    F_i^{m_1}(\widehat{d}^n) 
    = \big(\varphi_i,f(\widehat{U}^n)\big)
    \equiv\int_{\Omega}\varphi_i(x)f\Big(\sum_{\ell=1}^{m_1}\widehat{d}_\ell^n\varphi_\ell(x)\Big)dx,
    \\[5pt]
    &
    \beta\Big(\int_\Omega\widehat{\Theta}^ndx\Big)
    = \beta\Big(\int_\Omega\sum_{\ell=1}^{m_2}\widehat{c}_\ell^n\psi_\ell(x)dx\Big)
    = \beta\Big(\mathbf{b}\cdot \widehat{c}^n\Big),\quad\text{with}\;
    \mathbf{b}_i = \int_\Omega \psi_i(x)dx.
    \end{aligned}
    ```
    !!! note "Beware!"
        - Although ``q_1``, ``q_2``, and ``q_3`` are functions in the original model, we treat
          them as constants from this point on. This is a deliberate simplification adopted to reduce the complexity of the numerical implementation.
        - The index convention imposed on the basis functions of ``\mathcal{V}_{m_1}`` and ``\mathcal{V}_{m_3}``, namely
          ```math
                \phi_j = \varphi_j|_{\Gamma_1}, \quad j = 1, \dots, m_3,
          ```
          allows us to rewrite the matrices ``M^{m_1\times m_3}`` and ``M^{m_3\times m_1}`` as
          ```math
                M^{m_1\times m_3} = \begin{bmatrix} M^{m_3\times m_3} \\ 0^{(m_1-m_3)\times m_3} \end{bmatrix}
                \quad\text{and}\quad
                M^{m_3\times m_1} = \begin{bmatrix} M^{m_3\times m_3} & 0^{m_3\times(m_1-m_3)} \end{bmatrix}.
          ```
          Furthermore, the nonlinear term ``G^{m_1}(\widehat{v}^n) `` can be rewritten as
          ```math
          G^{m_1}(\widehat{v}^n) 
          = \begin{bmatrix} G^{m_3}(\widehat{v}_{1:m_3}^n)\\ 0^{(m_1-m_3)} \end{bmatrix},
          ```
          with
          ```math
          G^{m_3}_i(\widehat{v}_{1:m_3}^n) 
          =\int_{\Gamma_1}\phi_i(x)g\Big(x,\sum_{\ell=1}^{m_3}\widehat{v}_\ell^n\phi_\ell(x)\Big)d\Gamma.
          ```

## Solving the Nonlinear Algebraic System
In what follows, we show that the matrix formulation is equivalent to finding
``X \in \mathbb{R}^{m_1+m_2}`` such that ``H(X) = 0``, 
which we solve via Newton's method.
Starting from the initial guess 
``X_0 = [v^{n-1}; c^{n-1}]``, 
the method generates the sequence 
``X_{\kappa+1} = X_\kappa + S_\kappa``, 
where ``S_\kappa`` is obtained by solving the linear system
```math
JH(X_\kappa)\, S_\kappa = -H(X_\kappa),
```
with ``JH(X_\kappa)`` denoting the Jacobian of ``H`` evaluated at ``X_\kappa``.

### Formulation in Term of ``X=\begin{bmatrix}\hat{v}^n\\ \hat{c}^n\end{bmatrix}``
```math
Q(n)
\begin{bmatrix}
\hat{v}^n
\\[10pt] 
\hat{c}^n
\end{bmatrix}
+ 
\begin{bmatrix}
\tau\alpha^{n-\frac{1}{2}} G^{m_1}(\widehat{v}^n) + \tau F^{m_1}(\frac{\tau}{2}\widehat{v}^n + d^{n-1})
\\[10pt]
\tau\beta(\mathbf{b}\cdot\widehat{c}^n)K^{m_2\times m_2}\widehat{c}^n
\end{bmatrix}
-
\begin{bmatrix}
L_1(n)
\\[10pt]
L_2(n)
\end{bmatrix}
= 0.
```

Once ``\hat{v}^n`` and ``\hat{c}^n`` have been determined, we compute ``\hat{r}^n`` via
```math
\hat{r}^n
=
- \frac{\tau q_4}{q_5}\hat{v}_{1:m_3}^n
+ \frac{2q_1}{q_5} r^{n-1}
- \frac{\tau q_3}{q_5} z^{n-1}
+ \frac{\tau}{q_5} \Big(M^{m_3\times m_3}\Big)^{-1}
  \mathcal{F}^{m_3}(f_3^{n-\frac{1}{2}}).
```

The remaining quantities ``d^n``, ``v^n``, ``c^n``, ``z^n``, and ``r^n`` are then recovered by
```math
\begin{align*}
& v^n = 2\hat{v}^n - v^{n-1},
\quad c^n = 2\hat{c}^n - c^{n-1},
\quad r^n = 2\hat{r}^n - r^{n-1},
\quad d^n = \tau\hat{v}^n + d^{n-1},
\quad z^n = \tau\hat{r}^n + z^{n-1}.
\end{align*}
```

!!! details "Details"
    To rewrite the matrix formulation in terms of ``\widehat{v}^n``, ``\widehat{c}^n``, and ``\widehat{r}^n``, we make use of the following identities:
    ```math
    \bar\partial v^n = \frac{2}{\tau}\widehat{v}^n - \frac{2}{\tau}v^{n-1},
    \quad
    \widehat{d}^n = \frac{\tau}{2} \widehat{v}^n + d^{n-1},
    \quad
    \bar\partial r^n = \frac{2}{\tau}\widehat{r}^n - \frac{2}{\tau}z^{n-1},
    \quad
    \widehat{z}^n = \frac{\tau}{2} \widehat{r}^n + z^{n-1},
    \quad
    \bar\partial c^n = \frac{2}{\tau}\widehat{c}^n - \frac{2}{\tau}c^{n-1}.
    ```
    Applying these identities into the matrix formulation, we obtain:
    ```math
    \begin{align*}
    & M^{m_1\times m_1} \big(\frac{2}{\tau}\widehat{v}^n-\frac{2}{\tau}v^{n-1}\big)
    + \alpha^{n-\frac{1}{2}}\Big[  
        K^{m_1\times m_1} \big(\frac{\tau}{2}\widehat{v}^n+d^{n-1}\big)
      - M^{m_1\times m_3}\widehat{r}^n
      + G^{m_1}(\widehat{v}^n)\Big]
    + F^{m_1}(\frac{\tau}{2}\widehat{v}^n + d^{n-1})
    + A^{m_1\times m_2}\widehat{c}^n
    = \mathcal{F}^{m_1}(f_1^{n-\frac{1}{2}}),
    \\[10pt]
    & M^{m_2\times m_2}\big(\frac{2}{\tau}\widehat{c}^n-\frac{2}{\tau}c^{n-1}\big)
    + \beta(\mathbf{b}\cdot\widehat{c}^n)K^{m_2\times m_2}\widehat{c}^n
    + A^{m_2\times m_1}\widehat{v}^n
    = \mathcal{F}^{m_2}(f_2^{n-\frac{1}{2}}),
    \\[10pt]
    & M^{m_3\times m_3}\big[
      q_1\big(\frac{2}{\tau}\widehat{r}^n-\frac{2}{\tau}r^{n-1}\big)
    + q_2\widehat{r}^n
    + q_3\big(\frac{\tau}{2}\widehat{r}^n+z^{n-1}\big)\big]
    + q_4M^{m_3\times m_1}\widehat{v}^n
    = \mathcal{F}^{m_3}(f_3^{n-\frac{1}{2}}).
    \end{align*}
    ```
    Isolating ``\widehat{r}^n`` in the third equation:
    ```math
    (2q_1+\tau q_2+\frac{\tau^2}{2}q_3)M^{m_3\times m_3} \widehat{r}^n
    =
    - \tau q_4M^{m_3\times m_3} \widehat{v}_{1:m_3}^n
    + M^{m_3\times m_3}\big(2q_1r^{n-1} - \tau q_3z^{n-1}\big)
    + \tau \mathcal{F}^{m_3}(f_3^{n-\frac{1}{2}}).
    ```
    Denoting ``q_5 = 2q_1+\tau q_2+\frac{\tau^2}{2}q_3`` and using the result above in the first equation, we obtain:
    ```math
    \begin{align*}
    & M^{m_1\times m_1} \big(\frac{2}{\tau}\widehat{v}^n-\frac{2}{\tau}v^{n-1}\big)
    + \alpha^{n-\frac{1}{2}}\Big[  
        K^{m_1\times m_1} \big(\frac{\tau}{2}\widehat{v}^n+d^{n-1}\big)
      + G^{m_1}(\widehat{v}^n)\Big]
    + F^{m_1}(\frac{\tau}{2}\widehat{v}^n + d^{n-1})
    + A^{m_1\times m_2}\widehat{c}^n
    \\[10pt]
    &\qquad
    + \frac{\alpha^{n-\frac{1}{2}}}{q_5}
    \begin{bmatrix}
      \tau q_4M^{m_3\times m_3} \widehat{v}_{1:m_3}^n
      - M^{m_3\times m_3}\big( 2q_1r^{n-1} - \tau q_3z^{n-1} \big)
      - \tau \mathcal{F}^{m_3}(f_3^{n-\frac{1}{2}})
      \\[5pt]
      0^{(m_1-m_3)}
    \end{bmatrix}
    = \mathcal{F}^{m_1}(f_1^{n-\frac{1}{2}}),
    \\[10pt]
    & M^{m_2\times m_2}\big(\frac{2}{\tau}\widehat{c}^n-\frac{2}{\tau}c^{n-1}\big)
    + \beta(\mathbf{b}\cdot\widehat{c}^n)K^{m_2\times m_2}\widehat{c}^n
    + A^{m_2\times m_1}\widehat{v}^n
    = \mathcal{F}^{m_2}(f_2^{n-\frac{1}{2}}).
    \end{align*}
    ```
    Isolating ``\widehat{v}^n`` and ``\widehat{c}^n``:
    ```math
    \begin{align*}
    &
    \Big(
        2M^{m_1\times m_1}
        + \frac{\tau^2}{2}\alpha^{n-\frac{1}{2}} K^{m_1\times m_1}
        + \frac{\tau^2q_4}{q_5}\alpha^{n-\frac{1}{2}}
          \begin{bmatrix}
          M^{m_3\times m_3}       & 0^{m_3\times(m_1-m_3)}\\[5pt]
          0^{(m_1-m_3)\times m_3} & 0^{(m_1-m_3)\times(m_1-m_3)}
          \end{bmatrix}
    \Big) \widehat{v}^n
    + \tau A^{m_1\times m_2} \widehat{c}^n
    \\[10pt]
    &\qquad
    + \tau\alpha^{n-\frac{1}{2}} G^{m_1}(\widehat{v}^n) 
    + \tau F^{m_1}(\frac{\tau}{2}\widehat{v}^n + d^{n-1})
    - L_1(n) = 0,
    \\[10pt]
    &
    \tau A^{m_2\times m_1} \widehat{v}^n
    + 2M^{m_2\times m_2} \widehat{c}^n
    + \tau\beta(\mathbf{b}\cdot\widehat{c}^n)K^{m_2\times m_2}\widehat{c}^n
    - L_2(n) = 0.
    \end{align*}
    ```
    !!! note "Beware!"
        ```math
        \begin{align*}
        &
        \bar\partial v^n 
        = \frac{v^n-v^{n-1}}{\tau}
        = \frac{2\frac{v^n+v^{n-1}}{2}-v^{n-1}-v^{n-1}}{\tau}
        = \frac{2\widehat{v}^n-2v^{n-1}}{\tau}
        \quad\Rightarrow\quad
        \bar\partial v^n = \frac{2}{\tau}\widehat{v}^n - \frac{2}{\tau}v^{n-1}
        \\[10pt]
        &
        \bar\partial d^n = \widehat{v}^n
        \;\;\Rightarrow\;\;
        \frac{d^n-d^{n-1}}{\tau} = \widehat{v}^n
        \;\;\Rightarrow\;\;
        d^n = \tau \widehat{v}^n + d^{n-1}
        \;\;\Rightarrow\;\;
        \frac{d^n+d^{n-1}}{2} = \frac{\tau}{2} \widehat{v}^n + d^{n-1}
        \quad\Rightarrow\quad
        \widehat{d}^n = \frac{\tau}{2} \widehat{v}^n + d^{n-1}
        \end{align*}
        ```

!!! details "Matrix and vector definitions"
    ```math
    \begin{align*}
    Q =& 
    \begin{bmatrix}
    2M^{m_1\times m_1} 
    + \frac{\tau^2}{2}\alpha^{n-\frac{1}{2}}K^{m_1\times m_1}
    + \frac{\tau^2q_4}{q_5}\alpha^{n-\frac{1}{2}}
    \begin{bmatrix}
    M^{m_3\times m_3}       & 0^{m_3\times(m_1-m_3)}\\[5pt]
    0^{(m_1-m_3)\times m_3} & 0^{(m_1-m_3)\times(m_1-m_3)}
    \end{bmatrix}
    & 
    \tau A^{m_1\times m_2} 
    \\[30pt]
    \tau A^{m_2\times m_1} 
    &
    2M^{m_2\times m_2} 
    \end{bmatrix}
    \\[30pt]
    L_1 =&
    2M^{m_1\times m_1}v^{n-1}
    - \tau \alpha^{n-\frac{1}{2}}K^{m_1\times m_1} d^{n-1}
    + \tau \mathcal{F}^{m_1}(f_1^{n-\frac{1}{2}})
    \\
    &
    + \frac{\tau}{q_5}\alpha^{n-\frac{1}{2}}
    \begin{bmatrix}
      M^{m_3\times m_3}\big( 2q_1r^{n-1} - \tau q_3z^{n-1} \big)
      + \tau \mathcal{F}^{m_3}(f_3^{n-\frac{1}{2}})
      \\[5pt]
      0^{(m_1-m_3)}
    \end{bmatrix}
    \\[20pt]
    L_2 =& 
    2M^{m_2\times m_2} c^{n-1} 
    + \tau \mathcal{F}^{m_2}(f_2^{n-\frac{1}{2}})
    \end{align*}
    ```
    where ``q_5 = 2q_1+\tau q_2+\frac{\tau^2}{2}q_3``.

!!! details " Jacobian matrix calculation"
    Initially, note that:
    ```math
    H_i(X) = \sum_{\ell=1}^{m_1+m_2}Q_{i,\ell}X_\ell
    + 
    \begin{cases}\displaystyle
    \tau\alpha^{n-\frac{1}{2}}
    G_i^{m_1}\big(X_{1:m_1}\big) 
    + \tau F_i^{m_1}\big(\frac{\tau}{2}X_{1:m_1}+d^{n-1}\big)
    - [L_1]_i
    & \text{if } i\in\{1, \ldots, m_1\}
    \\[10pt]\displaystyle
    \tau\beta\big(\mathbf{b}\cdot X_{(m_1+1):(m_1+m_2)}\big)
    \sum_{\ell=1}^{m_2}K_{(i-m_1),\ell}^{m_2\times m_2}X_{m_1+\ell}
    -[L_2]_{(i-m_1)}
    & \text{if } i\in\{m_1+1, \ldots, m_1+m_2\}
    \end{cases}
    ```
    In this way,
    ```math
    \frac{\partial H_i}{\partial X_j}(X) 
    = Q_{i,j} 
    +
    \begin{cases}\displaystyle
    \tau\alpha^{n-\frac{1}{2}}
    \frac{\partial}{\partial X_j}G_i^{m_1}\big(X_{1:m_1}\big) 
    + \tau \frac{\partial}{\partial X_j}F_i^{m_1}\big(\frac{\tau}{2}X_{1:m_1}+d^{n-1}\big)
    & \text{if } i,j\in\{1, \ldots, m_1\}
    \\[10pt]
    0
    & \text{if } i\in\{1, \ldots, m_1\},\; j\in\{m_1+1, \ldots, m_1+m_2\}
    \\[10pt]
    0
    & \text{if } i\in\{m_1+1, \ldots, m_1+m_2\},\; j\in\{1, \ldots, m_1\}
    \\[10pt]\displaystyle
    \tau\beta^\prime\big(\mathbf{b}\cdot X_{(m_1+1):(m_1+m_2)}\big)\mathbf{b}_{j-m_1}
    \sum_{\ell=1}^{m_2}K_{(i-m_1),\ell}^{m_2\times m_2}X_{m_1+\ell}
    + \tau\beta\big(\mathbf{b}\cdot X_{(m_1+1):(m_1+m_2)}\big)
    K_{(i-m_1),(j-m_1)}^{m_2\times m_2}
    & \text{if } i,j\in\{m_1+1, \ldots, m_1+m_2\}
    \end{cases}
    ```
    where
    ```math
    G_i^{m_1}(v) 
    = \int_{\Gamma_1}\varphi_i(x)g\big(x,\sum_{\ell=1}^{m_1}v_\ell\varphi_\ell(x)\big)d\Gamma
    \;
    \Rightarrow
    \;
    \frac{\partial}{\partial X_j} G_i^{m_1}\big(X_{1:m_1}\big)
    = \int_{\Gamma_1}\varphi_i(x)\varphi_j(x)\frac{\partial g}{\partial s}\big(x,\sum_{\ell=1}^{m_1}X_\ell\varphi_\ell(x)\big)d\Gamma
    = JG\big(X_{1:m_1}\big)
    ```
    and
    ```math
    F_i^{m_1}(d) 
    = 
    \int_{\Omega}\varphi_i(x)f\Big(\sum_{\ell=1}^{m_1}d_\ell\varphi_\ell(x)\Big)dx
    \;
    \Rightarrow
    \;
    \frac{\partial}{\partial X_j}F_i^{m_1}\big(\frac{\tau}{2}X_{1:m_1}+d^{n-1}\big)
    = 
    \frac{\tau}{2}\int_{\Omega}\varphi_i(x)\varphi_j(x)f^\prime\Big(\sum_{\ell=1}^{m_1}
    \big[\frac{\tau}{2}X_\ell+d_\ell^{n-1}\big]
    \varphi_\ell(x)\Big)dx
    = \frac{\tau}{2} JF\big(\frac{\tau}{2}X_{1:m_1}+d^{n-1}\big).
    ```
    Consequently, the Jacobian matrix can be expressed as:
    ```math
    JH(X) 
    = Q
    + 
    \begin{bmatrix}
    \Big[
      \tau\alpha^{n-\frac{1}{2}} JG\big(X_{1:m_1}\big)
    + \frac{\tau^2}{2} JF\big(\frac{\tau}{2}X_{1:m_1}+d^{n-1}\big)
    \Big]^{m_1\times m_1}
    & 0^{m_1\times m_2}
    \\[10pt]
    0^{m_2\times m_1} &
    \tau\beta^\prime\big(\mathbf{b}\cdot X_{(m_1+1):(m_1+m_2)}\big)
    \Big[
    \big[K^{m_2\times m_2}X_{(m_1+1):(m_1+m_2)}\big]^{m_2\times 1}
    \times \big[\mathbf{b}^T\big]^{1\times m_2}
    \Big]^{m_2\times m_2}
    + \tau\beta\big(\mathbf{b}\cdot X_{(m_1+1):(m_1+m_2)}\big)K^{m_2\times m_2}
    \end{bmatrix}
    ```
    Given that ``\hat{v}^n=X_{1:m_1}`` and ``\hat{c}^n=X_{(m_1+1):(m_1+m_2)}``, we have
    ```math
    JH(X) 
    = Q
    + 
    \begin{bmatrix}
    \Big[
      \tau\alpha^{n-\frac{1}{2}} JG\big(\hat{v}^n\big)
    + \frac{\tau^2}{2} JF\big(\frac{\tau}{2}\hat{v}^n+d^{n-1}\big)
    \Big]^{m_1\times m_1}
    & 0^{m_1\times m_2}
    \\[10pt]
    0^{m_2\times m_1} &
    \tau\beta^\prime\big(\mathbf{b}\cdot \hat{c}^n\big)
    \Big[
    \big[K^{m_2\times m_2}\hat{c}^n\big]^{m_2\times 1}
    \times \big[\mathbf{b}^T\big]^{1\times m_2}
    \Big]^{m_2\times m_2}
    + \tau\beta\big(\mathbf{b}\cdot \hat{c}^n\big)K^{m_2\times m_2}
    \end{bmatrix}
    ```
    !!! note "Beware!"
        ```math
        JG^{m_1\times m_1}(v) =
        \begin{bmatrix}\displaystyle
        JG^{m_3\times m_3}(v_{1:m_3}) & 0^{m_3\times(m_1-m_3)}
        \\[10pt]
        0^{(m_1-m_3)\times m_3}       & 0^{(m_1-m_3)\times(m_1-m_3)}
        \end{bmatrix}
        ```