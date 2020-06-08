\appendix

# Steady-state Fock distribution for the resonator in RWA {#sec:qdlaser:pn-analytical-derivation}

In this Appendix, I derive the analytical expression in RWA for the steady-state 
Fock distribution $p_n$ of the harmonic oscillator in the single-atom laser setup (see Ch. -@sec:qdlaser), in the fully polarized case case ($P=1$), and at resonance ($\Delta\epsilon = \omega_0$). I further assume $\Gamma_L = \Gamma_R = \Gamma$. The derivation follows the seminal works of Scully and Lamb [@Scully1967;@Scully1997]. Starting from the master equation (-@eq:qdlaser:lindblad-equation), I label the matrix elements of $\rho$ as  $\rho^{nn'}_{ss'}$, with the upper index referring to the oscillator states
and the lower index referring to the quantum dot states
($s, s' \in \{0, \downarrow, \uparrow\}$).
The three equations for the diagonal dot part of the density matrix read

\begin{align}
    \label{eq:qdlaser:threeeq-1}
    \begin{split}
    \dot{\rho}^{n n'}_{\uparrow \uparrow} = &-i g \left(\sqrt{n+1} \rho^{n+1, n'}_{\downarrow \uparrow}- \sqrt{n'+ 1} \rho^{n, n'+1}_{\uparrow \downarrow}\right) \\
    &+ \frac{\kappa}{2}(1+n_B) \left( 2 \sqrt{(n+1)(n'+1)} \rho^{n+1, n'+1}_{\uparrow \uparrow} - (n+ n') \rho^{n n'}_{\uparrow \uparrow}\right) \\
    &+ \frac{\kappa}{2} n_B\left( 2\sqrt{n n'} \rho^{n-1, n'-1}_{\uparrow \uparrow} - (n+ n'+2) \rho^{n n'}_{\uparrow \uparrow} \right) \\
    & + \Gamma \rho^{n n'}_{00} ,
    \end{split} \\
    \label{eq:qdlaser:threeeq-2}
    \begin{split}
    \dot{\rho}^{n n'}_{\downarrow \downarrow} = &-i g \left(\sqrt{n} \rho^{n-1, n'}_{\uparrow \downarrow}- \sqrt{n'} \rho^{n, n'-1}_{\downarrow \uparrow}\right) \\
    &+ \frac{\kappa}{2}(1+n_B) \left( 2 \sqrt{(n+1)(n'+1)} \rho^{n+1, n'+1}_{\downarrow \downarrow} - (n+ n') \rho^{n n'}_{\downarrow \downarrow}\right) \\
    &+ \frac{\kappa}{2} n_B \left( 2\sqrt{n n'} \rho^{n-1, n'-1}_{\downarrow \downarrow} - (n+ n'+2) \rho^{n n'}_{\downarrow \downarrow} \right) \\
    &- \Gamma \rho^{n n'}_{\downarrow \downarrow}, 
    \end{split} \\
    \label{eq:qdlaser:threeeq-3}
    \begin{split}
    \dot{\rho}^{n n'}_{00} &= \frac{\kappa}{2}(1+n_B) \left( 2 \sqrt{(n+1)(n'+1)} \rho^{n+1, n'+1}_{00} - (n+ n') \rho^{n n'}_{00}\right) \\
    &+ \frac{\kappa}{2} n_B \left( 2\sqrt{n n'} \rho^{n-1, n'-1}_{00} - (n+ n'+2) \rho^{n n'}_{00} \right) \\
    &- \Gamma \rho^{n n'}_{00} + \Gamma \rho^{nn'}_{\downarrow \downarrow}.
    \end{split}
\end{align}

The equations couple to the off-diagonal elements of the dot ($\uparrow \downarrow$ and $\downarrow \uparrow$) and the $n+1, n'+1$ and $n-1, n'-1$ elements of the resonators, and thus cannot form a closed set of equations. 
However, one can find a closed form using the fact that $\{\Gamma; g\} \gg \kappa$: the electron dynamics is much faster than the decay and excitation processes of the resonator through the thermal bath. We therefore neglect the terms proportional to $\kappa$ in Eqs. (-@eq:qdlaser:threeeq-1)-(-@eq:qdlaser:threeeq-3). When calculating the equations of motion for 
$\rho^{n+1,n'}_{\downarrow \uparrow}$ and $\rho^{n,n'+1}_{\uparrow \downarrow}$ we obtain the system:

\begin{align}
 \label{eq:qdlaser:system_1}
 \dot{\rho}^{n n'}_{\uparrow \uparrow} = &-i g \left[\sqrt{n+1} \rho^{n+1, n'}_{\downarrow \uparrow}- \sqrt{n'+ 1} \rho^{n, n'+1}_{\uparrow \downarrow}\right] + \Gamma \rho^{n n'}_{00} , \\
  \label{eq:qdlaser:system_2}
   \dot{\rho}^{n, n'+1}_{\uparrow \downarrow} = &-i g \left[\sqrt{n+1} \rho^{n+1, n'+1}_{\downarrow \downarrow}- \sqrt{n'+1} \rho^{n, n'}_{\uparrow \uparrow}\right] - \frac{\Gamma}{2} \rho^{n, n'+1}_{\uparrow \downarrow}, \\
    \label{eq:qdlaser:system_3}
   \dot{\rho}^{n+1, n'}_{\downarrow \uparrow} = &-i g \left[\sqrt{n+1} \rho^{n, n'}_{\uparrow \uparrow}- \sqrt{n'+1} \rho^{n+1, n'+1}_{\downarrow \downarrow}\right] - \frac{\Gamma}{2} \rho^{n+1, n'}_{\downarrow \uparrow}, \\
 \label{eq:qdlaser:system_4}
 \dot{\rho}^{n+1, n'+1}_{\downarrow \downarrow} = &-i g \left[\sqrt{n+1} \rho^{n, n'+1}_{\uparrow \downarrow}- \sqrt{n'+1} \rho^{n+1, n'}_{\downarrow \uparrow}\right] - \Gamma \rho^{n+1, n'+1}_{\downarrow \downarrow}, \\
  \label{eq:qdlaser:system_5} 
 \dot{\rho}^{n n'}_{00} = &- \Gamma \rho^{n n'}_{00} + \Gamma \rho^{nn'}_{\downarrow \downarrow}.
\end{align}

Since we are interested in the stationary state, we put all time derivatives to zero. Now, from the last equation, we obtain  $\rho^{n n'}_{00} = \rho^{nn'}_{\downarrow \downarrow}$, and from the relation 
$\rho^{nn'} = \rho^{n n'}_{00} + \rho^{nn'}_{\downarrow \downarrow} + \rho^{nn'}_{\uparrow \uparrow}$ we obtain $\rho^{n n'}_{00} = \frac{1}{2} \rho^{nn'}  - \frac{1}{2} \rho^{nn'}_{\uparrow \uparrow},$ which can be inserted in Eq. (-@eq:qdlaser:system_1). We are left with a linear system of four unknowns 
which has the matrix form

$$
\mathbf{\underline{\underline{M}}}\mathbf{R} - \mathbf{A} = \mathbf{0}.
$$

We have defined:

$$
\mathbf{R} = \begin{pmatrix} \rho^{nn'}_{\uparrow \uparrow} \\  \dot{\rho}^{n, n'+1}_{\uparrow \downarrow} \\ \dot{\rho}^{n+1, n'}_{\downarrow \uparrow} \\ \dot{\rho}^{n+1, n'+1}_{\downarrow \downarrow} \end{pmatrix}, \quad \mathbf{A} = \dfrac{\Gamma}{2} \begin{pmatrix} \rho^{n n'} \\ 0 \\ 0 \\ 0 \end{pmatrix}, \quad
 \underline{\underline{\mathbf{M}}} = \begin{pmatrix}
\frac{\Gamma}{2} & - i g \sqrt{n' + 1} & i g \sqrt{n+1} & 0 \\
- i g \sqrt{n' + 1} & \frac{\Gamma}{2} & 0 & i g \sqrt{n+1} \\
i g \sqrt{n+1} & 0 & \frac{\Gamma}{2} & - i g \sqrt{n' + 1}  \\
0 & i g \sqrt{n+1} & - i g \sqrt{n' + 1}  & \Gamma 
\end{pmatrix}.
$$

The steady-state equation for the matrix elements of the oscillator reads:

$$
\begin{split}
0 &= - \left[ \frac{\mathscr{N}' _{n n'} \kappa (g/g_\mathrm{thr})^2}{1 + \mathscr{N}_{n n'} (A_s^2)^{-1} (g/g_\mathrm{thr})^2}\right] \rho^{n n'} + \left[ \frac{\sqrt{n n'} \kappa (g/g_\mathrm{thr})^2}{1 + \mathscr{N}_{n-1, n'-1}(A_s^2)^{-1} (g/g_\mathrm{thr})^2 }\right] \rho^{n-1, n'-1} \\
&+ \frac{\mathscr{C}_d}{2} \left[ 2 \sqrt{(n+1)(n'+1)} \rho^{n+1, n'+1} - (n + n')\rho^{n n'} \right] \\
&+ \frac{\mathscr{C}_e}{2} \left[2 \sqrt{n n'} \rho^{n-1, n'-1} - (n+1+n'+1) \rho^{n n'}\right],
\end{split}
$$

with the following quantities:

\begin{gather}
g_\mathrm{thr} = \frac{1}{2} \sqrt{\Gamma \kappa}\quad \text{(threshold coupling)}, \qquad A_s^2 = \frac{\Gamma}{3 \kappa}\quad \text{(saturation number)} \\
\mathscr{N}_{nn'} = \frac{1}{2}(n+1 + n'+1) + \frac{(n-n')^2}{36} \frac{1}{A_s^2}\left(\frac{g}{g_\mathrm{thr}} \right)^2, \\
\mathscr{N}'_{nn'} = \frac{1}{2}(n+1 + n'+1) + \frac{(n-n')^2}{12} \frac{1}{A_s^2}\left(\frac{g}{g_\mathrm{thr}} \right)^2, \\
\mathscr{C}_d = \kappa(1 + n_\mathrm{B}), \qquad \mathscr{C}_e = \kappa n_B.
\end{gather}

We are interested in the diagonal elements, and thus we put $n = n'$, obtaining the equation for the distribution $p_n$:

$$
\begin{split}
&- \left[ \frac{(n+1) \kappa \left(\frac{g}{g_\mathrm{thr}}\right)^2}{1 + (n+1) (A_s^2)^{-1}\left( \frac{g}{g_\mathrm{thr}}\right)^2} \right]p_n + \left[ \frac{n \kappa \left(\frac{g}{g_\mathrm{thr}}\right)^2}{1 + n\ (A_s^2)^{-1} \left(\frac{g}{g_\mathrm{thr}}\right)^2} \right]p_{n-1} \\
&-  \mathscr{C}_d n p_n +  \mathscr{C}_d (n+1)p_{n+1} - \mathscr{C}_e (n+1)p_n +  \mathscr{C}_e n p_{n-1}\\
&= 0.
\end{split}
$$
{#eq:qdlaser:equation_distribution}

To find the steady-state distribution which solves Eq. (-@eq:qdlaser:equation_distribution) we first rearrange the terms:

$$
\begin{split}
0 = &-\left\{\left[ \frac{(n+1) \kappa \left(\frac{g}{g_\mathrm{thr}}\right)^2}{1 + (n+1) (A_s^2)^{-1}\left( \frac{g}{g_\mathrm{thr}}\right)^2} \right]p_n -  \mathscr{C}_d (n+1)p_{n+1} + \mathscr{C}_e (n+1)p_n\right\} \\
&+  \left[ \frac{n \kappa \left(\frac{g}{g_\mathrm{thr}}\right)^2}{1 + n\ (A_s^2)^{-1} \left(\frac{g}{g_\mathrm{thr}}\right)^2} \right]p_{n-1} -  \mathscr{C}_d n p_n +  \mathscr{C}_e n p_{n-1}.
\end{split}
$$

Then, we notice that the equation is in the form  $0 = f_n - f_{n-1}$ and thus $f$ has to be zero (detailed balance condition).
We obtain the recursive equation:

$$
\left[ \frac{n \kappa\tilde{g}^2}{1 + n \frac{\tilde{g}^2}{A_s^2}} + \mathscr{C}_e n \right] p_{n-1} = \mathscr{C}_d p_n,
$$
{#eq:qdlaser:pn-solution-appendix}

with the shorthand $\tilde{g} = g/g_\mathrm{thr}$. Equation (-@eq:qdlaser:pn-solution-appendix) is equivalent to Eq. (-@eq:qdlaser:rwa-pn).

The solution can be written in terms of a finite product series 
or rising factorials:

$$
p_n = p_0 \prod_{k=1}^n \left[\frac{ \frac{\tilde{g}^2}{n_B+1}}{1 + k \frac{ \tilde{g}^2}{A_s^2}} + \frac{n_B}{n_B+1} \right] 
= p_0 \frac{\left(\dfrac{A_s^2}{n_B} + \frac{A_s^2}{\tilde{g}^2}+1\right)_n}{\left( \dfrac{A_s^2}{\tilde{g}^2} + 1 \right)_n} \left(\frac{n_B}{n_B+1} \right)^n,
$$
{#eq:qdlaser:steadystate_distribution}

equivalent to Eq. (-@eq:qdlaser:rwa-pn-solution).
The derivation outlined here assumes $P=1$ and equal tunneling rates $\Gamma_L = \Gamma_R$, but the general case of finite polarization and arbitrary rates simply renormalizes the saturation number $A_s^2$ and the threshold coupling $g_\mathrm{thr}$, yielding the expressions in Eq. (-@eq:qdlaser:rwa-pn-saturation-thr). The form of the solution for $p_n$ remains as in Eq. (-@eq:qdlaser:steadystate_distribution). 


# Operators $\mathcal{G}(\Delta_0, n)$, $\mathcal{F}(\Delta_0, n)$ and matrix element $\langle 1 | H_\text{eff} | 2 \rangle$ {#sec:jjcavity:matrix-elements}

I describe here the nonlinear operators $\mathcal{G}(\Delta_0, n)$ and
$\mathcal{F}(\Delta_0, n)$ appearing in the effective Hamiltonian (@eq:jjcavity:heff) of the
JJ-resonator system described in Ch. -@sec:jjcavity. They are defined by
the relations:
\begin{align}
\mathcal{G}(\Delta_0, n) n &= \sum_{p=1}^{\infty} \frac{4p}{p^2 - 1} [A_p,
A^\dagger_p], \\
a \mathcal{F}(\Delta_0, n) &= \sum_{p=0}^{\infty} \frac{(-1)^p}{2p + 1} [A_p,
A^\dagger_{p+1}],
\end{align}
where $A_p = (a^\dagger)^p K_p$, and $K_p$ is an hermitian operator defined as
$$
K_p = : \frac{J_p (2 \Delta_0 \sqrt{n})}{n^{p/2}} : = \sum_{m=0}^{\infty}
\frac{(-1)^m \Delta_0^{2m + p} (a^\dagger)^m a^m}{m! (m+p)!}.
$$
Effective Hamiltonian (@eq:jjcavity:heff) can be written in the Fock basis as
$$
H_\text{eff} = \sum_{q=0}^\infty (q \delta + \delta E_q) |q \rangle \langle q| +
i \sum_{q=0}^\infty \left[M_{q,q+1} |q \rangle \langle q +1 | e^{-2i\phi} - \text{H.c.} \right].
$$
The matrix elements are:
$$
\begin{aligned}
&\delta E_q = \frac{\tilde{E}_J^2}{4 \omega_J} \left \{ \sum_{p=1}^q
\frac{4p}{4p^2 - 1} \left[ \frac{\kappa_p^2(q - p)q!}{(q-p)!} \right] -
\sum_{p=1}^\infty \frac{4p}{4p^2 - 1} \left[ \frac{\kappa_p^2 (q) (q+p)!}{q!}\right]\right \},\\
&\begin{split}
    M_{q,q+1} = \frac{\tilde{E}_J^2}{4 \omega_J} & \left \{ \sum_{p=0}^q \frac{(-1)^p}{2p+1}
   \frac{\sqrt{q!(q+1)!}}{(q-p)!} \kappa_p (q - p)\kappa_{p+1}(q - p) \right. \\
              -&\left.\,\sum_{p=0}^\infty\frac{(-1)^p}{2p + 1}\frac{(q + p + 1)!}{\sqrt{q!(q+1)!}} \kappa_{p+1} (q)\kappa_{p}(q+1) \right \}.
\end{split}
\end{aligned}
$$

$\kappa_p (q)$ represents the $q$-th eigenvalue of $K_p$ and is given by
$$
\kappa_p (q) = q! \sum_{n=0}^q \frac{(-1)^n \Delta_0^{2n + p}}{n! (n+p)! (q - n)!}.
$$
In this representation, specific matrix elements of the effective Hamiltonian
can be written in closed form. We report here the element $M_{12}$ which is of
particular relevance as explained in @sec:jjcavity:g2:
$$
M_{12} = \frac{\tilde{E}_J^2}{4\sqrt{2}\omega_J} \left[ \Delta_0 e^{-\Delta_0^2}
\left(\frac{2}{3} \Delta_0^4 - \frac{10}{3} \Delta_0^2 + \frac{3}{2} \right) +
\sqrt{\pi} \text{erf} (\Delta_0) \left(\frac{2}{3} \Delta_0^6 - 3\Delta_0^4 + \frac{7}{2}\Delta_0^2-
\frac{3}{4} \right) \right],
$$
where $\text{erf}(x)$ is the error function. The first zero of $M_{12}$ occurs
at $\Delta_0 \approx 1.07$.