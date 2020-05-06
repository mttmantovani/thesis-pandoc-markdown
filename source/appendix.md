\appendix

# Operators $\mathcal{G}(\Delta_0, n)$, $\mathcal{F}(\Delta_0, n)$ and matrix element $\langle 1 | H_\text{eff} | 2 \rangle$ {#sec:jjcavity:matrix-elements}

We describe here the nonlinear operators $\mathcal{G}(\Delta_0, n)$ and
$\mathcal{F}(\Delta_0, n)$ appearing in the effective Hamiltonian (@eq:jjcavity:heff) of the
JJ-resonator system described in Chapter {@jjcavity}. They are defined by
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
