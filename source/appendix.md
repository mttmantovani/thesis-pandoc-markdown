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

# Full-counting statistics of charge transport {#sec:qdlaser:full-counting-statistics}


The full-counting statistics (FCS) of charge transport is a method for evaluating the full probability distribution of electron detected by a collector, in order to gain additional information on the transport, beyond the average current. Originally, the FCS method was devised within quantum optics, in the contest of photodetection [@Mandel1979;@Cook1981]. Later, it has been extendend to the study of quantum transport in mesoscopic conductors [@Levitov1993;@Levitov1996;@Blanter2000;@Bagrets2003;@Flindt2005]. Here, I briefly outline the theory of FCS and its application to the calculation of the current noise used in @sec:qdlaser:current.

## Basic concepts 

Assuming that the system dynamics is governed by a Markovian master equation as described in Ch. -@sec:theory, 

$$
\dot{\rho} = \mathcal{L} \rho,
$$
{#eq:qdlaser:fcs-master-equation}

we are interested in resolving the dynamics with respect to the detection or no-detection of electrons from, e.g., the right lead of a system.
We can write the total Liouvillian $\mathcal{L}$ as $\mathcal{L} = \mathcal{L}_0 + \mathcal{L}_1$, where $\mathcal{L}_0$ contains the detection-free evolution (the isolated dynamics of the system and the un-monitored detections) and $\mathcal{L}_1$ containes the detections. 
The propagator solving Eq. (-@eq:qdlaser:fcs-master-equation) has the general form

$$
\mathcal{P}(t) = e^{\mathcal{L}t}.
$$

To make explicit the appearance of the quantum jumps (electron being detected in the right lead) in the evolution, we consider the Laplace transform  $\mathcal{P}(z) = \int_0^\infty \od t \mathcal{P}(t) e^{-zt}$ of the propagator, obtaining

$$
\mathcal{P}(z) =\frac{1}{z \mathds{1} - \mathcal{L}_0 - \mathcal{L}_1}= \frac{1}{(z\mathds{1} - \mathcal{L}_0)\left(\mathds{1} - \dfrac{\mathcal{L}_1}{z\mathds{1} - \mathcal{L}_0}\right)}.
$$
{#eq:full_propagator_lapl}

Introducing the free propagator

$$
\mathcal{P}_0(z) = \frac{1}{z\mathds{1} - \mathcal{L}_0} = \int_0^\infty \od t e^{\mathcal{L}_0 t} e^{-z t},
$$

we can expand Eq. (-@eq:full_propagator_lapl) as

$$
\mathcal{P}(z) = \sum_{n=0}^\infty [\mathcal{P}_0 (z) \mathcal{L}_1]^n \mathcal{P}_0 (z) = \mathcal{P}_0 (z) + \mathcal{P}_0 (z) \mathcal{L}_1\mathcal{P}_0 (z)  + \mathcal{P}_0 (z) \mathcal{L}_1\mathcal{P}_0 (z) \mathcal{L}_1\mathcal{P}_0 (z) + \dots
$$

Using the convolution theorem, it can be shown that in the time domain we obtain

$$
\mathcal{P}(t) = e^{\mathcal{L}_0 t }+ \int_0^t \od t_1 e^{\mathcal{L}_0 (t - t_1)} \mathcal{L}_1 e^{\mathcal{L}_0 t_1} + \int_0^t \od t_2 \int_0^{t_2} \od t_1  e^{\mathcal{L}_0 (t - t_2)}\mathcal{L}_1 e^{\mathcal{L}_0 (t_2 -t_1)}\mathcal{L}_1 e^{\mathcal{L}_0 t_1} + \dots
$$
{#eq:full_propagator_time}

which can be interpreted as periods of free evolution interrupted by single jump events. Since $\rho (t) = \mathcal{P}(t)\rho_0$, where $\rho_0$ is the density matrix at $t=0$, we can write the probability for $n$ jump events to occur during time $t$:

$$
P_n (t) = \Tr  \left\{\int_0^t \od t_n  \dots \int_0^{t_2} \od t_1 e^{\mathcal{L}_0 (t - t_n) }\mathcal{L}_1 \dots \mathcal{L}_1 e^{\mathcal{L}_0 (t_2 - t_1)} \mathcal{L}_1 e^{\mathcal{L}_0 t_1 } \rho_0 \right\}.
$$
{#eq:p_n_t_probability}

In Laplace space, the above equation reads simply

$$
P_n (z) = \Tr \left\{ [\mathcal{P}_0 (z) \mathcal{L}_1]^n \mathcal{P}_0 (z) \rho_0\right\}.
$$

It is possible to calculate $P_n (t)$ even if we know only the full propagator, $\mathcal{P}(t)$. To this scope, we introduce the counting variable $\chi$, conjugate to $n$, satisfying the orthonormality relation

$$
\label{eq:orthonormality_relation}
\frac{1}{2\pi} \int_{-\pi}^{\pi} e^{in\chi} e^{-i m \chi} \od \chi = \delta_{nm}.
$$

Now, we make the replacement

$$
\mathcal{L}_1 \rightarrow \mathcal{L}_1 e^{i \chi}, \qquad \mathcal{L} \rightarrow \mathcal{L}(\chi) = \mathcal{L}_0 + \mathcal{L}_1 e^{i \chi},
$$

in the full propagator $\mathcal{P}(\chi, t) = e^{\mathcal{L}(\chi) t}$, now depending explicitly on the counting variable $\chi$.
Using Eq. (-@eq:orthonormality_relation) in Eq. (-@eq:p_n_t_probability) we obtain

$$
\label{eq:p_n_t_chi}
P_n (t) = \frac{1}{2\pi} \int_{-\pi}^{\pi} \od \chi \Tr \left\{\mathcal{P}(\chi, t) \rho_0\right\} e^{-i n \chi}, \qquad P_n (z) = \frac{1}{2\pi} \int_{-\pi}^{\pi} \od \chi \Tr \left\{\mathcal{P}(\chi, z) \rho_0\right\} e^{-i n \chi}.
$$

## Vector representation

For practical calculations, it is common to use a representation in which $\mathcal{L}$ is a matrix and $\rho$ is a column vector. 
To avoid confusion with the usual quantum mechanical operators acting in the Hilbert space of the usual state vectors (bras and kets) of the system, 
we call Liouville space the linear space of the usual operators. 
Operators acting in the Liouville space are called superoperators. 
In the following, superoperators are denoted by calligraphic symbols; operators acting in the Hilbert space will be distinguished from the usual bras and kets through double brackets, i.e., with the correspondence $O \leftrightarrow |o \rrangle$.
We can write a matrix representation for the superoperator $\mathcal{O}$ in a similar way to the Hilbert space case:

$$
\mathcal{O} = \sum_{kl} \sum_{mn} | k l \rrangle \llangle kl | \mathcal{O} | mn \rrangle =  \sum_{kl} \sum_{mn} O_{kl,mn} | kl \rrangle \llangle mn |.
$$

Here, $\{|m n \rrangle \equiv |m \rangle \langle n | \}$ is the system of all the projectors onto states forming an orthonormal basis of the Hilbert space,  which in turn is an orthonormal basis for the Liouville space, with respect to the scalar product $\llangle a | b \rrangle = \Tr (A^\dagger B)$.
The trace is performed on the system. 
We denote with $O_{kl,mn}$ the matrix elements of the superoperator, having four indices. 
The operator $O = \sum_{m, n} O_{mn} |m \rangle \langle n |$ represented 
by the matrix $O_{mn}$ corresponds to the vector $| o \rrangle =\sum_{m, n} O_{mn} |m  n \rrangle$, represented by a column vector with entries $\mathbf{o}= (O_{11}, O_{12}, O_{13}, \dots, O_{N1}, O_{N2}, \dots, O_{NN})^T$, where $N$ is the dimension of the Hilbert space. 
The exact ordering depends on the chosen ordering of the indices. 
We can now write the total Liovullian $\mathcal{L}$ with its eigendecomposition:

$$
	\mathcal{L} = \sum_{k=0}^{N^2-1} \lambda_k |\phi_k \rrangle \llangle \phi_k|.
$$

Since the Liouvillian is in general non-hermitian, the left eigenvector $|\phi_k \rrangle$ and the right eigenvector $\llangle \phi_k|$ are in general different. However, since we assume that the system reaches a stationary state, there is an eigenvalue $\lambda_0 = 0$, with right eigenvector $|\phi_0 \rrangle \equiv | 0 \rrangle = \rho_\mathrm{st}$. 
The stationary state is assumed unique and normalized. From the orthonormality condition, $\llangle \phi_j |  \phi_k \rrangle = \delta_{jk}$ we deduce that the left eigenvector corresponding to the eigenvalue $\lambda_0 = 0$ corresponds to the identity matrix in the Liouville space, or equivalently to performing the trace $\Tr (\rho_\mathrm{st}) = 1$ in the Hilbert space.
We denote with $\llangle \tilde{0} |$ the left eigenvalue.
Similarly, for the $\chi$-dependent Liouvillian we have the decomposition:

$$
\mathcal{L}(\chi) = \sum_{k=0}^{N^2-1} \lambda_k(\chi) |\phi_k (\chi)\rrangle \llangle \phi_k(\chi)|,
$$

with the assumption that $\lambda_0 (\chi) \rightarrow 0$ for $\chi \rightarrow 0$, and all the other eigenvalues have a larger negative real part near $\chi = 0$.

## Moments and cumulants

We can identify in Eq. (-@eq:p_n_t_chi) the moment-generating function (MGF):

$$
	\label{eq:mgf}
	M(\chi, t) = \Tr \left\{e^{\mathcal{L}(\chi) t} \rho_0 \right\},
$$

since it is defined as the Fourier transform of the probability distribution. 

From the MGF it is possible to compute all the moments by differentiation 
with respect to the counting field:

$$
	\langle n^k \rangle_t = \sum_{n=0}^{\infty} n^k P_n (t) = 
					(-i)^n \left.\frac{\partial M (\chi, t)}{\partial \chi^n} \right|_{\chi = 0},
$$

where $\langle n^k \rangle_t$ denotes the $k$-th moment of the probability 
distribution at time $t$. 
Note that the MGF can also be written as

$$
	M(\chi, t) = \sum_{n=0}^{\infty} e^{in \chi} P_n(t),
$$

which is the inverse Fourier transform of Eq. (-@eq:p_n_t_chi). 
The cumulant-generating function (CGF) is obtained by taking the logarithm of the MGF, 

$$
	C(\chi, t) = \ln M(\chi, t),
$$

and all the cumulants are recovered from 

$$
	\llangle n^k \rrangle_t = (-i)^n \left.\frac{\partial C (\chi, t)}{\partial \chi^n} \right|_{\chi = 0}.
$$

We recall that, in contrast to moments, higher cumulants are inert 
when a trivial transformation such as a translation is performed on the probability distribution.
It is possible to show that in the long-time limit, $t\rightarrow \infty$, 
the CGF is linear in time and depends on the lowest eigenvalue of the 
$\chi$-dependent Liouvillian $\mathcal{L}(\chi, t)$. 
Evaluating Eq. (-@eq:mgf) on the stationary state, and using the vector 
representation, we can write

$$
	M(\chi, t) = e^{C(\chi, t)} = 
			\sum_k \llangle \tilde{0} | \phi_k (\chi) \rrangle \llangle \phi_k (\chi) | 0 \rrangle e^{\lambda_k (\chi) t}.
$$

In the long-time limit, the only term that survives is that with $\lambda_0 (\chi)$ in the exponent, since all other terms yield exponentially damped contributions. 
Therefore, we have

$$
	e^{C(\chi, t)} \approx c(\chi) e^{\lambda_0 (\chi) t},
$$

with $c(\chi) = \llangle \tilde{0} | \phi_0 (\chi) \rrangle \llangle \phi_0 (\chi) | 0 \rrangle$. 
Taking the logarithm, we obtain $C(\chi, t) = \lambda_0(\chi)t + \ln c$, and therefore, to exponential accuracy

$$
	\label{eq:cgf_long_time}
	C(\chi, t) \approx \lambda_0(\chi)t.
$$

We can therefore find every cumulant by differentiation of the $\chi$-dependent eigenvalue $\lambda_0$.

## Expression for the average current and the current noise

Equation (-@eq:cgf_long_time) expresses the CGF for the distribution $P_n(t)$ 
in the long-time limit, which is the probability of $n$ jumps events occurring
in a time interval $t$. 
By identifying the jump processes with electrons entering 
the right lead, we can derive Eq. (-@eq:cgf_long_time) with respect to time and find the current cumulants in the stationary limit:

$$
\llangle I^m \rrangle = e^{m} \dt{} \lim_{t \rightarrow \infty} \llangle n^m \rrangle (t),
$$

with $e$ being the electron charge. 
To achieve this, we resolve the density matrix $\rho(t)$ and the master equation with respect to $n$. We introduce the superoperator $\mathcal{J}^{(e)}$ which is relative to the tunneling of electrons from the system to the right lead, and the superoperator $\mathcal{J}^{(h)}$, relative to the tunneling from the right lead to the system (viz., holes entering the collector). 
The $n$-resolved master equation has the form:

$$
	\dot{\rho}^{(n)}(t) = (\mathcal{L} -\mathcal{J}^{(e)} -\mathcal{J}^{(h)})\rho(t) + \mathcal{J}^{(e)}\rho^{(n-1)}(t) + \mathcal{J}^{(h)}\rho^{(n+1)}(t).
$$

It is possible to find the expression for the average current as

$$
I = 	\llangle I \rrangle = e \llangle \tilde{0} | \mathcal{J}^{(e)} - \mathcal{J}^{(h)} |0 \rrangle 
		=e \Tr [(  - \mathcal{J}^{(e)} - \mathcal{J}^{h}) \rho_\mathrm{st}],
$$

while for the second cumulant, i.e., the zero-frequency current noise, 
one finds [@Flindt2005]

$$
\label{eq:noise_fcs}
S(0) = 	\llangle I^2 \rrangle = e\Tr [(\mathcal{J}^{(e)} +\mathcal{J}^{(h)})\rho_\mathrm{st}] 
			- e^2  \Tr [(\mathcal{J}^{(e)} -\mathcal{J}^{(h)}) \mathcal{R} (\mathcal{J}^{(e)} -\mathcal{J}^{(h)}) \rho_\mathrm{st}].
$$

We have introduced the pseudoinverse of the Liouvillian $\mathcal{R} = \mathcal{Q} \mathcal{L}^{-1} \mathcal{Q}$, 
where $\mathcal{Q}$ is the projector out of the null-space of $\mathcal{L}$, 
which is spanned by $\rho_\mathrm{st}$. 
If $\mathcal{P} \equiv | 0 \rrangle \llangle \tilde{0} |$ is the projector onto 
the stationary state then $\mathcal{Q} = \mathds{1} - \mathcal{P}$. 
The pseudoinverse $\mathcal{R}$ is well defined, since the inversion is performed in 
the subspace spanned by $\mathcal{Q}$, where $\mathcal{L}$ is regular.
In the infinite-bias limit, transport from right to left is blocked, 
thus $\mathcal{J}^{(h)} = 0$ and $\mathcal{J}^{(e)} \pm \mathcal{J}^{(h)}  = \mathcal{J}^{(e)} \equiv \mathcal{J}$.
Then we have

$$
I = e \Tr [\mathcal{J} \rho_\mathrm{st}], \qquad 
S(0) = eI - 2e^2 \Tr[ \mathcal{J} \mathcal{R} \mathcal{J} \rho_\mathrm{st}],
$$


## Note on the numerical evaluation of the noise

The quantum dot + resonator system described in Ch. -@sec:qdlaser requires working with a very large Hilbert space, due to the necessity of including many Fock states (up to $n \approx 500$) to capture the lasing behavior with low quality factors. As a consequence, the actual calculation of $\mathcal{R}$ 
is unfeasible. 
However, it is not necessary to evaluate $\mathcal{R}$, but rather its action onto a vector, $\mathcal{R}|y \rrangle$ [@Flindt2004;@Flindt2010].
The trick to evaluate the second part of Eq. (-@eq:noise_fcs) is the following. 
We first write the linear system

$$
	\label{eq:linear_system}
	\mathcal{L} | x \rrangle = \mathcal{Q} | y \rrangle,
$$

where $|x \rrangle$ is an unknown vector and $| y \rrangle$ is known. 
The r.h.s. ensures that the solution lies in the range of $\mathcal{L}$. 
Since $\mathcal{L}$ is singular, the system has infinitely many solutions which we write in the form

$$
	\label{eq:linear_system_solution}
	|x \rrangle = |x_0 \rrangle + c |0 \rrangle,\quad c \in \mathds{C}.
$$

$|x_0 \rrangle$ is a particular solution which is found numerically. 
Now, we apply $\mathcal{R}$ to both sides of Eq. (-@eq:linear_system)
and obtain

$$
	\begin{split}
		\mathcal{R}\mathcal{L} | x \rrangle &= \mathcal{R}\mathcal{Q} | y \rrangle \\
		\mathcal{Q}\mathcal{L}^{-1}\mathcal{Q}\mathcal{L} | x \rrangle &= \mathcal{Q}\mathcal{L}^{-1}\mathcal{Q}^2 | y \rrangle \\
		\mathcal{Q}\mathcal{L}^{-1}\mathcal{L}\mathcal{Q} | x \rrangle &=  \mathcal{Q}\mathcal{L}^{-1}\mathcal{Q} | y \rrangle \\
		\mathcal{Q}| x \rrangle &= \mathcal{R}| y \rrangle.
	\end{split}
$$

We have used the definition of $\mathcal{R}$ and the identities 
$\mathcal{Q}^2 = \mathcal{Q}$ 
and 
$\mathcal{Q} \mathcal{L} = \mathcal{L}\mathcal{Q}$. 
Therefore, we find numerically the solution $|x_0 \rrangle$ to the system  (-@eq:linear_system) and then we apply the projector $\mathcal{Q}$ to Eq. (-@eq:linear_system_solution), 
yielding $\mathcal{Q}| x_0 \rrangle = \mathcal{R}| y \rrangle$, since $\mathcal{Q} |0 \rrangle = 0$. 
In our case, we need to choose $| y \rrangle = \mathcal{I} |0 \rrangle$.
The solution of the linear system Eq. (-@eq:linear_system_solution) revealed itself to be a formidable numerical problem due to the size of the matrices involved. We have used the iterative GMRES method [@Saad1986] to solve the  preconditioned system

$$
\mathcal{M}\mathcal{L} | x \rrangle = \mathcal{M}\mathcal{Q} | y \rrangle.
$$

The matrix $\mathcal{M}$ is an approximation of the pseudoinverse 
of $\mathcal{L}$ and it has been crucial for the convergence of the GMRES method. 
A wrong or absent preconditioner made the convergence impossible or extremely slow. 
The preconditioner was built in our case with a sparse incomplete 
LU decomposition (spILU) for the Liouvillian $\mathcal{L}$ [@Saad2003].


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