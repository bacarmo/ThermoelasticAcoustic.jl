# Scheme 2
The following fully discrete scheme combines the Crank-Nicolson Galerkin method with selective linearization of nonlinear terms and a decoupling strategy between the wave and heat equations. 
It consists of finding ``U^n, V^n \in \mathcal{V}_{m_1}``, ``\Theta^n\in \mathcal{V}_{m_2}``, and ``Z^n, R^n \in \mathcal{V}_{m_3}`` such that
```math
\begin{align*}
& \big(\varphi,\bar{\partial}V^n\big)
+ \alpha^{n-\frac{1}{2}}\Big[
    \big(\nabla\varphi,\nabla\widehat{U}^n\big)
  - \big(\varphi,\widehat{R}^n)_{\Gamma_1}
  + \big(\varphi,g(\widehat{V}^n)\big)_{\Gamma_1}\Big]
+ \big(\varphi,f(U^{\ast n})\big) 
+ \big(\varphi,(\mathbf{a}\cdot\nabla)\Theta^{\ast n}\big) 
= \big(\varphi,f_1^{n-\frac{1}{2}}\big),
\quad\forall\varphi\in \mathcal{V}_{m_1},
\\
& \big(\psi,\bar{\partial}\Theta^n\big)
+ \beta\Big(\int_\Omega\Theta^{\ast n}dx\Big)\big(\nabla\psi,\nabla\widehat{\Theta}^n\big)
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
for ``n = \text{“1,0''},\,1,\,2,\,\ldots``, with ``U^0, V^0 \in \mathcal{V}_{m_1}``, ``\Theta^0 \in \mathcal{V}_{m_2}``, and ``Z^0, R^0 \in \mathcal{V}_{m_3}`` given as approximations of the initial solutions ``u_0``, ``\, v_0``, ``\,\theta_0``, ``\, z_0``, and ``r_0``.
These approximations are defined in the same manner as in Scheme 1.

!!! details "Notation"
    In addition to the operators 
    ``\displaystyle\bar{\partial}w^n=\frac{w^n - w^{n-1}}{\tau}`` and 
    ``\displaystyle\widehat{w}^n = \frac{w^n + w^{n-1}}{2}``, consider
    ```math
    \bar{\partial}w^{\text{“1,0''}} = \frac{w^{\text{“1,0''}} - w^0}{\tau},\quad
    \widehat{w}^{\text{“1,0''}} = \frac{w^{\text{“1,0''}} + w^0}{2},
    \quad\text{and}\quad
    w^{*n} = 
    \begin{cases}\displaystyle
    w^0,                             & \text{if } n = \text{“1,0''},
    \\[10pt] \displaystyle
    \frac{w^{\text{“1,0''}}+w^0}{2}, & \text{if } n = 1,
    \\[10pt] \displaystyle
    \frac{3w^{n-1}-w^{n-2}}{2},      & \text{if } n \geq 2.
    \end{cases}
    ```

## Matrix Formulation
```math
\begin{align*}
& M^{m_1\times m_1}\bar{\partial}v^n
+ \alpha^{n-\frac{1}{2}}\Big[  
    K^{m_1\times m_1}\widehat{d}^n
  - M^{m_1\times m_3}\widehat{r}^n
  + G^{m_1}(\widehat{v}^n)\Big]
+ F^{m_1}(d^{\ast n})
+ A^{m_1\times m_2}c^{\ast n}
= \mathcal{F}^{m_1}(f_1^{n-\frac{1}{2}}),
\\
& M^{m_2\times m_2}\bar{\partial}c^n
+ \beta(\mathbf{b}\cdot c^{\ast n})K^{m_2\times m_2}\widehat{c}^n
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
\bar{\partial}z^n = \widehat{r}^n.
\end{align*}
```
## Solving the Algebraic Systems
The matrix formulation can be rewritten as:
```math
\begin{align*}
&
\big[Q_1(n)\big]^{m_1\times m_1} \widehat{v}^n
+ \tau\alpha^{n-\frac{1}{2}} G^{m_1}(\widehat{v}^n) 
- L_1(n,d^{\ast n},c^{\ast n}) = 0,
\\[10pt]
&
\big[Q_2(c^{\ast n})\big]^{m_2\times m_2}\widehat{c}^n
= L_2(n,\widehat{v}^n),
\\[10pt]
& \hat{r}^n
=
- \frac{\tau q_4}{q_5}\hat{v}_{1:m_3}^n
+ \frac{2q_1}{q_5} r^{n-1}
- \frac{\tau q_3}{q_5} z^{n-1}
+ \frac{\tau}{q_5} \Big(M^{m_3\times m_3}\Big)^{-1}
  \mathcal{F}^{m_3}(f_3^{n-\frac{1}{2}}).
\end{align*}
```
The three systems are solved sequentially. 
The first is a nonlinear system in ``\widehat{v}^n``, independent of ``\widehat{c}^n`` and ``\widehat{r}^n``, and is solved via Newton's method.
Once ``\widehat{v}^n`` is available, the second equation reduces to a linear system in ``\widehat{c}^n``, which is solved directly.
Finally, ``\widehat{r}^n`` is recovered explicitly from the third equation.

With ``\hat{v}^n``, ``\hat{c}^n``, and ``\hat{r}^n`` determined, the remaining unknowns are updated via
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
    Rewriting the matrix formulation in terms of ``\widehat{v}^n``, ``\widehat{c}^n``, and ``\widehat{r}^n``, we obtain:
    ```math
    \begin{align*}
    & M^{m_1\times m_1} \big(\frac{2}{\tau}\widehat{v}^n-\frac{2}{\tau}v^{n-1}\big)
    + \alpha^{n-\frac{1}{2}}\Big[  
        K^{m_1\times m_1} \big(\frac{\tau}{2}\widehat{v}^n+d^{n-1}\big)
      - M^{m_1\times m_3}\widehat{r}^n
      + G^{m_1}(\widehat{v}^n)\Big]
    + F^{m_1}(d^{\ast n})
    + A^{m_1\times m_2} c^{\ast n}
    = \mathcal{F}^{m_1}(f_1^{n-\frac{1}{2}}),
    \\[10pt]
    & M^{m_2\times m_2}\big(\frac{2}{\tau}\widehat{c}^n-\frac{2}{\tau}c^{n-1}\big)
    + \beta(\mathbf{b}\cdot c^{\ast n})K^{m_2\times m_2}\widehat{c}^n
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
    + F^{m_1}(d^{\ast n})
    + A^{m_1\times m_2}c^{\ast n}
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
    + \beta(\mathbf{b}\cdot c^{\ast n})K^{m_2\times m_2}\widehat{c}^n
    + A^{m_2\times m_1}\widehat{v}^n
    = \mathcal{F}^{m_2}(f_2^{n-\frac{1}{2}}),
    \\[10pt]
    & \hat{r}^n
    =
    - \frac{\tau q_4}{q_5}\hat{v}_{1:m_3}^n
    + \frac{2q_1}{q_5} r^{n-1}
    - \frac{\tau q_3}{q_5} z^{n-1}
    + \frac{\tau}{q_5} \Big(M^{m_3\times m_3}\Big)^{-1}
      \mathcal{F}^{m_3}(f_3^{n-\frac{1}{2}}).
    \end{align*}
    ```
    Isolating ``\widehat{v}^n`` and ``\widehat{c}^n``:
    ```math
    \begin{align*}
    &
    Q_1(n) \widehat{v}^n
    + \tau\alpha^{n-\frac{1}{2}} G^{m_1}(\widehat{v}^n) 
    - L_1(n,d^{\ast n},c^{\ast n}) = 0,
    \\[10pt]
    &
    Q_2(c^{\ast n})\widehat{c}^n
    = L_2(n,\widehat{v}^n),
    \\[10pt]
    & \hat{r}^n
    =
    - \frac{\tau q_4}{q_5}\hat{v}_{1:m_3}^n
    + \frac{2q_1}{q_5} r^{n-1}
    - \frac{\tau q_3}{q_5} z^{n-1}
    + \frac{\tau}{q_5} \Big(M^{m_3\times m_3}\Big)^{-1}
      \mathcal{F}^{m_3}(f_3^{n-\frac{1}{2}}).
    \end{align*}
    ```

!!! details "Matrix and vector definitions"
    ```math
    \begin{align*}
    Q_1 =& 
    2M^{m_1\times m_1} 
    + \frac{\tau^2}{2}\alpha^{n-\frac{1}{2}}K^{m_1\times m_1}
    + \frac{\tau^2q_4}{q_5}\alpha^{n-\frac{1}{2}}
    \begin{bmatrix}
    M^{m_3\times m_3}       & 0^{m_3\times(m_1-m_3)}\\[5pt]
    0^{(m_1-m_3)\times m_3} & 0^{(m_1-m_3)\times(m_1-m_3)}
    \end{bmatrix}
    \\[10pt]
    Q_2 =& 2M^{m_2\times m_2} + \tau\beta\big(\mathbf{b}\cdot c^{\ast n}\big)K^{m_2\times m_2}
    \\[10pt]
    L_1 =&
    - \tau F^{m_1}(d^{\ast n})
    - \tau A^{m_1\times m_2} c^{\ast n}
    + 2M^{m_1\times m_1}v^{n-1}
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
    - \tau A^{m_2\times m_1} \widehat{v}^n
    + 2M^{m_2\times m_2} c^{n-1} 
    + \tau \mathcal{F}^{m_2}(f_2^{n-\frac{1}{2}})
    \end{align*}
    ```

!!! details " Jacobian matrix calculation"
    Initially, note that:
    ```math
    H_i(X) 
    = \sum_{\ell=1}^{m_1}[Q_1]_{i,\ell}X_\ell
    + \tau\alpha^{n-\frac{1}{2}} G_i^{m_1}\big(X\big) 
    - \big[L_1\big]_i.
    ```
    In this way,
    ```math
    JH(X) 
    = Q_1 + \tau\alpha^{n-\frac{1}{2}}
      \begin{bmatrix}\displaystyle
      JG^{m_3\times m_3}(X_{1:m_3}) & 0^{m_3\times(m_1-m_3)}
      \\[10pt]
      0^{(m_1-m_3)\times m_3}       & 0^{(m_1-m_3)\times(m_1-m_3)}
      \end{bmatrix} ,
    ```
    where
    ```math
    \big[JG^{m_3\times m_3}(X_{1:m_3})\big]_{i,j}
    =
    \int_{\Gamma_1} \psi_i(x)\psi_j(x)
    \frac{\partial g}{\partial s}\Big(x,\sum_{\ell=1}^{m_3}X_\ell\psi_\ell(x)\Big)d\Gamma.
    ```