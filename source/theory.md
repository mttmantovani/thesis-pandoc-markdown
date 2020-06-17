# Theory of the quantum master equation {#sec:theory}

## Introduction

Any realistic investigation of a quantum system requires taking into account its
coupling to the rest of the world. Therefore, the derivation of a suitable
equation describing the dynamics of the density matrix of a system is
the primary ingredient both in the theoretical analysis and in the validation of
experimental data. In this Chapter, I will focus on a class of master
equations that can describe faithfully the solid-state devices presented in the
following Chapters of the Thesis. Typically, these systems follow a
Markovian dynamics, since they are weakly coupled to an environment with many
degrees of freedom, and system-environment correlations are lost quickly when
compared to the timescale on which the state of the system varies appreciably.

The most general form of a Markovian master equation can be rigorously
formulated in terms of the generator of a quantum dynamical semigroup, giving rise to the standard
Gorini-Kossakowski-Sudarshan-Lindblad (GKLS or, more commonly, Lindblad) equation [@Gorini1976;@Lindblad1976]. Here, I will instead focus on a derivation of the master equation which follows a well-developed microscopic approach and leads to the Bloch-Redfield equation [@Redfield1957;@Carmichael1993;@Breuer2002]. I will discuss critically the validity regime of the Bloch-Redfield equation, showing the condition under which it has a Lindblad form.
After presenting the master equation, I will briefly summarize a few commonly used methods (analytical and numerical) to solve it.

## Microscopic derivation of the quantum master equation

![Schematic picture of an open quantum system. The universe may be subdivided into a system weakly coupled to an environment through the interaction Hamiltonian $H_I$.](figures/theory-open-system.pdf){#fig:theory:open-system}

The general recipe to treat an open quantum system is to identify a *system*,
representing the portion of the universe that one wants to
closely investigate, which lives in an Hilbert space $\mathcal{H}_S$ and is described by a Hamiltonian $H_S$ and a density matrix $\rho_S$. The system is coupled to an *environment*, described by a state $\rho_E$ in an Hilbert space $\mathcal{H}_E$, and Hamiltonian $H_E$. The coupling, realized through an interaction Hamiltonian $H_I$, induces system-environment correlations, such that the *reduced* dynamics of $\rho_S$ can no longer be regarded as unitary (see @fig:theory:open-system).
The total Hamiltonian of the universe may be written as

$$
    H(t) = H_S \otimes \mathds{1}_E + \mathds{1}_S \otimes H_E + H_I(t),
$$
{#eq:theory:total-hamiltonian}

where $\mathds{1}_{\{E, S\}}$ are identity operators acting on the corresponding Hilbert space. The total Hamiltonian can be, in principle, time-dependent.

In order to make the problem tractable, many assumptions and approximations must be in order: One needs to isolate the relevant parts to be comprised in the model, discarding effects and interactions which can be safely neglected, and the system-environment coupling must be considered weak in
a sense that will be later clarified. In a typical solid-state or QED
device, the system is composed of a few-level atom, an artificial atom (a
quantum dot or a superconducting qubit), a mechanical or electromagnetic
resonator, as well as any coupled arrangement of the above. The environment
usually consists of metallic contacts (fermionic baths of quasicontinuous
states) or a bath of harmonic oscillators (substrate phonons, thermal
excitations, quasiparticles) in thermal equilibrium, to which the
system is naturally coupled, giving rise to
relaxation and decoherence.

### Liouville-von Neumann equation

The unitary evolution of the closed, isolated universe under Hamiltonian
\eqref{eq:theory:total-hamiltonian} is given by the Liouville-von Neumann
equation for $\rho$ [@Blum2012],

$$
    \dot{\rho}(t) = -i [H(t), \rho(t)] \equiv \mathcal{L}(t)\rho(t).
$$
{#eq:theory:liouville}

We have introduced the Liovullian superoperator $\mathcal{L}$, which will be of
central importance throughout this work. The general solution to Eq. (-@eq:theory:liouville) is given by

$$
    \dot{\rho}(t) = \mathcal{T} \exp \left[\int_0^t d \tau \mathcal{L}(\tau) \right] \rho(0),
$$
{#eq:theory:liouville-evolution}

where $\mathcal{T}$ is the time-ordering operator and $\rho(0)$ denotes the
state at the initial time. For a time-independent Liouvillian, Eq. (-@eq:theory:liouville-evolution) reduces to

$$
    \dot{\rho}(t) = \exp \left[\mathcal{L} t \right] \rho(0).
$$
{#eq:theory:liouville-evolution-time-independent}

Equation (-@eq:theory:liouville-evolution) is equivalent to

$$
    \dot{\rho}(t) = U(t,0) \rho(0) U^\dagger(t,0),
$$
{#eq:theory:liouville-evolution-unitary-operator}

with the time-evolution operator 

$$
    U(t,0) = \mathcal{T} \exp \left[-i \int_0^t d \tau H(\tau) \right],
$$
{#eq:theory:evolution-operator-solution}

which is unitary, $U(t,0)U^\dagger(t,0) = U^\dagger(t,0)U(t,0) = \mathds{1}$,
and satisfies

$$
   i \frac{\partial}{\partial t}U(t,0) = H(t) U(t,0).
$$

### Interaction picture {#sec:theory:interaction-picture}

Equations (-@eq:theory:liouville)-(-@eq:theory:liouville-evolution-unitary-operator)
are written in the Schrödinger picture, where the time evolution of the system
falls entirely on the state. The first step to obtain an equation for the
reduced system density matrix $\rho_S$ is move to the interaction picture with
the respect to the coupling Hamiltonian $H_I$. Let us write the total Hamiltonian as $H(t) = H_0 + H_I(t)$, where the evolution (i.e., the spectrum) of $H_0 = H_S + H_E$ is known and given by 

$$
    U_0(t,0) = \exp[-iH_0 t].
$$

The goal of the interaction picture is to split the global evolution as

$$
    U(t,0) = U_0(t,0) U_I(t,0),
$$
{#eq:theory:evolution-operator-interaction-picture}

where the subscript $(\cdot)_I$ stands for "interaction picture". For an operator $O$, the evolution of the expectation value on the full state $\rho$
can be written as 

$$
    \begin{split}
\langle O(t) \rangle &= \Tr [O \rho (t)] = \Tr [O  U(t,0) \rho(0) U^\dagger(t,0)]
= \Tr [O U_0(t,0) U_I(t,0) \rho(0) U_I^\dagger(t,0) U_0^\dagger(t,0)] \\
    &= \Tr [O_I(t) \rho_I(t)].
    \end{split}
$$

In the last line, we have used the cyclic property of the trace and defined

\begin{align}
    O_I(t) &= U_0^\dagger(t,0)O U_0(t,0),  \label{eq:theory:operator-interaction-picture}\\
    \rho_I(t) &= U_I^\dagger (t,0) \rho(0) U_I (t,0) = U_0^\dagger(t,0) \rho(t)
    U_0(t,0). \label{eq:theory:rho-interaction-picture}
\end{align}

Differentiating Eq. (-@eq:theory:rho-interaction-picture) with respect to time and using Eq. (-@eq:theory:evolution-operator-interaction-picture) leads to the Liouville-von Neumann equation in the interaction picture:

$$
    \dot{\rho}_I (t) = -i [H_I(t), \rho_I(t)].
$$
{#eq:theory:liouville-equation-interaction-picture}

Notice that $H_I(t)$ is the interaction Hamiltonian transformed to the
interaction picture according to Eq. (-@eq:theory:operator-interaction-picture). 

To start the formal derivation of the master equation, it will be necessary to perform a partial trace over the degrees of freedom of the environment, in order to get an equation for the reduced density matrix of the system, $\rho_S = \Tr_E [\rho]$. To this end, we will consider a specific form of the interaction Hamiltonian $H_I$, given by

$$
H_I = \sum_\alpha A_\alpha \otimes B_\alpha,
$$
{#eq:theory:interaction-hamiltonian-form}

where $A_\alpha$ are system operators and $B_\alpha$ are bath operators. In the interaction picture, Eq. (-@eq:theory:interaction-hamiltonian-form) transforms into

$$
H_I(t) = \sum_\alpha A_{\alpha,I}(t) \otimes B_{\alpha,I}(t),
$$

using Eq. (-@eq:theory:operator-interaction-picture). We will assume throughout for simplicity hermitian operators $A_\alpha = A_\alpha^\dagger$ and $B_\alpha = B_\alpha^\dagger$, although this is not a strict requirement, as long as $H_I$ is hermitian as a whole (see Appendix -@sec:theory:aux:1).

The standard derivation of the master equation stems from perturbation theory applied to Eq. (-@eq:theory:liouville-equation-interaction-picture), under the assumption of *weak coupling* between the system and the environment, followed by the subsequent application of the Born, Markov and secular approximations  [@Breuer2002]. Let us first integrate formally Eq. (-@eq:theory:liouville-equation-interaction-picture) and then insert the solution back on its right-hand side. After tracing over the environment, we obtain

$$
\dot{\rho}_S = -i \Tr_E \{ [H_I(t), \rho(0)]\} - \int_0^t \Tr_E \{ [H_I(t), [H_I(t'), \rho(t')]] \} dt'.
$$
{#eq:theory:liouville-equation-iteration}

Equation (-@eq:theory:liouville-equation-iteration) is exact, but depends on the full density matrix $\rho$ at all previous times. 

### Born approximation

The Born approximation relies on the fact that the environment is large, such that it is barely perturbed by the system, and the system-environment coupling is small: We assume that $H_I(t) \sim \mathcal{O}(\varepsilon)$, where $\varepsilon$ is a small perturbative dimensionless parameter. The density matrix at all times is then assumed to be of the form

$$
\rho(t) = \rho_S(t) \otimes \rho_E + \mathcal{O}(\varepsilon),
$$

with the factorized initial density matrix $\rho(0) = \rho_S(0) \otimes \rho_B$. Notice that the presence of the $\mathcal{O}(\varepsilon)$ term is crucial to induce correlations between the system and the environment, which allow the state of the system to evolve. After the Born approximation we obtain a perturbative expansion of Eq. (-@eq:theory:liouville-equation-iteration) which is valid to second order in $\varepsilon$:

$$
\dot{\rho}_S =  -i \Tr_E \{ [H_I(t), \rho_S(0) \otimes \rho_E]\} - \int_0^t \Tr_E \{ [H_I(t), [H_I(t'), \rho_S(t') \otimes \rho_E]] \} dt' + \mathcal{O}(\varepsilon^3).
$$
{#eq:theory:master-equation-born-approximation}

<!-- The truncation of the Lioville-von Neumann equation to second order requires a justification, which will be made clear later. -->

We now make the assumption

$$
\Tr \{ B_\alpha (t) \rho_E \} = 0.
$$
{#eq:theory:single-particle-exp-val}

This is not a restrictive condition, as it is always possible to modify the system Hamiltonian and the bath operators $B_\alpha$ to let the trace vanish (see Appendix -@sec:theory:aux:2). Using the interaction Hamiltonian decomposition (-@eq:theory:interaction-hamiltonian-form) and condition (-@eq:theory:single-particle-exp-val), Eq. (-@eq:theory:master-equation-born-approximation) becomes

$$
\dot{\rho}_S = -\sum_{\alpha \beta} \int_0^t dt' \left\{C_{\alpha\beta} (t,t') [A_\alpha (t), A_\beta (t') \rho_S(t')] + C_{\beta\alpha}(t',t) [\rho_S(t') A_\beta(t'), A_\alpha(t)] \right\},
$$
{#eq:theory:master-equation-born-approximation-corr-fun}

where we have defined the environmental correlation functions

$$
C_{\alpha\beta} (t, t') = \Tr \{ B_\alpha(t) B_\beta (t') \rho_E\}.
$$
{#eq:theory:environment-correlation-function}


### Markov approximation

Equation (-@eq:theory:master-equation-born-approximation-corr-fun) is closed (depends only on $\rho_S$), but it is non-Markovian, as $\rho_S$ must be known at all previous times. However, under the weak-coupling and large-reservoir approximations, it is possible to obtain a Markovian expression.

The next assumption in our treatment consists in assuming that the state of the environment is stationary, i.e., it is in thermal equilibrium,

$$
\rho_E = \frac{e^{-\beta H_E}}{\Tr \{e^{-\beta H_E}\}},
$$

where $\beta$ is the environmental inverse temperature. As a consequence, the environment Hamiltonian commutes with the bath density operator, $[H_E, \rho_E] = 0$. The correlation functions for a stationary bath satisfy the property

$$
C_{\alpha\beta}(t,t') = C_{\alpha\beta}(\tau \equiv t - t'),
$$

with

$$
C_{\alpha\beta}(\tau) = \Tr \{ e^{iH_E \tau }B_\alpha e^{-i H_E \tau }B_\beta \rho_E\}.
$$

Assuming hermitian bath operators $B_\alpha$, we will have additionally $C_{\alpha \beta}(\tau) = C^*_{\beta\alpha}(-\tau)$. The central idea behind the Markov approximation is that, when the environment is large and its spectrum becomes quasi-dense, the correlation functions $C_{\alpha\beta}(\tau)$ will by strongly peaked around $\tau = 0$ and will decay to zero *much faster* than the rate of variation of $\rho_S$. The consequence of this is twofold: first, we are allowed to replace $\rho_S(t')$ with $\rho_S(t)$ in Eq. (-@eq:theory:master-equation-born-approximation-corr-fun); second, we can push the integration limits to $t \rightarrow \infty$. In order to do this, we just make the change of variable $\tau = t - t'$ in Eq. (-@eq:theory:master-equation-born-approximation-corr-fun), obtaining

$$
\dot{\rho}_S = - \sum_{\alpha\beta} \int_0^\infty d\tau \{ C_{\alpha\beta} [A_\alpha(t), A_\beta(t - \tau) \rho_S (t)] + C_{\beta\alpha}(-\tau) [\rho_S(t)A_\beta(t - \tau),A_\alpha(t)]\}. 
$$
{#eq:theory:bloch-redfield}

After transforming back to the Schrödinger picture, we obtain

$$
\begin{split}
\dot{\rho}_S = &-i [H_S, \rho_S(t)] \\
&- \sum_{\alpha\beta} \int_0^\infty d\tau \{ C_{\alpha\beta} [A_\alpha,e^{-i H_S \tau} A_\beta e^{i H_S \tau} \rho_S (t)] + C_{\beta\alpha}(-\tau) [\rho_S(t) e^{-i H_S \tau} A_\beta e^{i H_S \tau},A_\alpha]\}.
\end{split}
$$
{#eq:theory:bloch-redfield-schroedinger}

Equation (-@eq:theory:bloch-redfield-schroedinger) is known as the Bloch-Redfield master equation [@Redfield1957;@Carmichael1993;@Breuer2002]; it is local in time and trace-preserving, but it is not guaranteed to preserve the positivity of the density matrix, because it cannot be generally put in a Lindblad form. To deal with this, the secular approximation must be employed.

### Secular approximation

To apply the secular approximation, it is crucial to start with the *interaction picture* master equation (-@eq:theory:bloch-redfield). Let us define the energy eigenbasis of the system Hamiltonian

$$
H_S |n \rangle = E_n |n \rangle,
$$

allowing the possibility of degeneration ($E_n = E_m$ for $n \neq m$). We introduce the *spectral decomposition* of the system operators $A_\alpha$ through the relations

\begin{align}
A_\alpha (\omega) &= \sum_{nm} \delta (\omega_{mn} - \omega) |n \rangle \langle n | A_\alpha |m \rangle \langle m |, \\
A_\alpha^\dagger (\omega) &= \sum_{nm} \delta (\omega_{nm} - \omega) |n \rangle \langle n | A_\alpha^\dagger |m \rangle \langle m |,
\end{align}

with the Bohr frequencies $\omega_{mn} = E_m - E_n$. The $A_\alpha(\omega)$ satisfy $\sum_\omega A_\alpha(\omega) = A_\alpha$, and in the interaction picture they become

$$
e^{i H_S t} A_\alpha(\omega) e^{-i H_S t} = e^{-i\omega t} A_{\alpha} (\omega).
$$
{#eq:theory:spectral-decomposition-int-pic}

Introducing Eq. (-@eq:theory:spectral-decomposition-int-pic) into Eq. (-@eq:theory:bloch-redfield) we arrive at

$$
\dot{\rho}_S = \sum_{\alpha \beta} \sum_{\omega \omega'} e^{i(\omega' - \omega)t} \Gamma_{\alpha \beta} (\omega) \left[A_\beta(\omega) \rho_S(t) A^\dagger_\alpha (\omega') - A^\dagger_\alpha(\omega') A_\beta (\omega) \rho_S(t) \right] + \mathrm{H.c.},
$$
{#eq:theory:bloch-redfield-after-spectral-decomposition}

with the one-sided Fourier transforms of the bath correlation functions

$$
\Gamma_{\alpha \beta} (\omega) = \int_0^\infty C_{\alpha \beta} (\tau) e^{i\omega \tau} d\tau.
$$

The secular approximation consists in neglecting the terms with $\omega'\neq \omega$ in Eq. (-@eq:theory:bloch-redfield-after-spectral-decomposition). For this to make sense, one must assume that the timescale on which $\rho_S$ varies appreciably is *much larger* than the typical values of $|\omega' - \omega|^{-1}$, such that the fast oscillating exponential terms in Eq. (-@eq:theory:bloch-redfield-after-spectral-decomposition) can be averaged out to zero if $\omega' \neq \omega$. Upon the secular approximation, we obtain

$$
\dot{\rho}_S = \sum_{\alpha \beta} \sum_{\omega} \Gamma_{\alpha \beta} (\omega) \left[A_\beta(\omega) \rho_S(t) A^\dagger_\alpha (\omega) - A^\dagger_\alpha(\omega) A_\beta (\omega) \rho_S(t) \right] + \mathrm{H.c.}
$$
{#eq:theory:bloch-redfield-after-secular}

Now, we split the function $\Gamma_{\alpha \beta} (\omega)$ into

$$
\begin{aligned}
 \Gamma_{\alpha \beta}(\omega) &=\frac{1}{2} \gamma_{\alpha \beta}(\omega)+\frac{1}{2} \sigma_{\alpha \beta}(\omega), \\ \Gamma_{\beta \alpha}^{*}(\omega) &=\frac{1}{2} \gamma_{\alpha \beta}(\omega)-\frac{1}{2} \sigma_{\alpha \beta}(\omega),
\end{aligned}
$$
{#eq:theory:damping-coefficients-1}

where

$$
\begin{array}{l}\gamma_{\alpha \beta}(\omega)=\Gamma_{\alpha \beta}(\omega)+\Gamma_{\beta \alpha}^{*}(\omega)=\int_{-\infty}^{+\infty} C_{\alpha \beta}(\tau) e^{+i \omega \tau} d \tau, \\ \sigma_{\alpha \beta}(\omega)=\Gamma_{\alpha \beta}(\omega)-\Gamma_{\beta \alpha}^{*}(\omega)=\int_{-\infty}^{+\infty} C_{\alpha \beta}(\tau) \operatorname{sgn}(\tau) \mathrm{e}^{+\mathrm{i} \omega \tau} \mathrm{d} \tau.\end{array}
$$
{#eq:theory:damping-coefficients-2}

With the help of Eq. (-@eq:theory:damping-coefficients-1), rearranging the terms in Eq. (-@eq:theory:bloch-redfield-after-secular) yields the Lindblad form

$$
\dot{\rho}_S = -i [H_{LS}, \rho_S(t)] + \sum_\omega \sum_{\alpha,\beta} \gamma_{\alpha\beta}(\omega) \left[ A_\beta(\omega) \rho_S(t) A^\dagger_\alpha(\omega) - \frac{1}{2} \{ A^\dagger_\alpha(\omega) A_\beta(\omega), \rho_S(t)\} \right].
$$

The notation $\{ \cdot, \cdot \}$ refers to the anticommutator. The operator

$$
H_{LS} = \sum_{\omega} \sum_{\alpha\beta} \frac{1}{2i} \sigma_{\alpha\beta}(\omega) A^\dagger_\alpha(\omega) A_\beta (\omega)
$$

is called Lamb-shift Hamiltonian, and its effect is a renormalization of the unperturbed energy levels of the system, arising from the system-environment coupling. In most practical applications (including the ones presented in this work) the Lamb shift is usually neglected. It can be shown that the matrix $\gamma_{\alpha\beta} (\omega)$ is hermitian and positive [@Breuer2002], hence it is diagonalizable with positive damping coefficients yielding a diagonal form of the master equation. In the special case where the energy dependence of the damping coefficients can be neglected, the Lindblad form is particularly simple and given by

$$
\dot{\rho}_S = -i[H_{LS}, \rho_S(t)] - \sum_{\alpha} \gamma_\alpha \mathcal{D}(A_\alpha)\rho_S(t),
$$
{#eq:theory:lindblad-equation-diagonal}

where I have introduced the *Lindblad dissipator* with operator $A$ acting on the state $\rho$,

$$
\mathcal{D}(A)\rho = A \rho A^\dagger - \frac{1}{2}\{ A^\dagger A, \rho \}.
$$


### Master equation in the energy representation

A sometimes useful representation of the master equation (without employing the spectral decomposition) is found by projecting Eq. (-@eq:theory:bloch-redfield) in the energy eigenbasis $|n\rangle$, which yields

$$
\dot{\rho}_S = \sum_{\alpha\beta} \sum_{nmpq} \Gamma_{\alpha\beta}(\omega_{mn})e^{i[\omega_{mn} - \omega_{qp}]t} (A_\beta)_{nm} (A_\alpha)_{pq}^* \left\{ L_{nm} \rho_S (t) L_{pq}^\dagger - L^\dagger_{pq} L_{nm} \rho_S(t) \right\} + \mathrm{H.c.}
$$

We have introduced the shorthands $(A_\alpha)_{nm} = \langle n | A_\alpha | m \rangle$ and $L_{nm} = |n \rangle \langle m|$ (the latter is the Lindblad jump operator in the energy representation). The secular approximation amounts here to neglecting all terms where $\omega_{nm} \neq \omega_{qp}$. The Lindblad form of the master equation in the Schrödinger picture reads

$$
\dot{\rho}_S = -i [H_S + H_{LS}, \rho_S(t)] + \sum_{nmpq} \gamma_{nm,pq} \left[ L_{nm} \rho_S(t) L_{pq}^\dagger - \frac{1}{2} \{L_{pq}^\dagger L_{nm}, \rho_S(t) \}\right].
$$

The Lamb-shift Hamiltonian is $H_{LS} = \sum_{nm} \sigma_{nm} L_{nm}$, with matrix elements

$$
\sigma_{n m}=-\frac{i}{2} \delta_{\omega_{mn}}\sum_{\alpha \beta} \sum_{p}  \sigma_{\alpha \beta}\left(\omega_{mp}\right) (A_\beta)_{pm} (A_\alpha)_{pn}^*,
$$

and the damping coefficients in the dissipative part read explicitly

$$
\gamma_{nm,pq} = \delta_{\omega_{mn},\omega_{qp}} \sum_{\alpha\beta} \gamma_{\alpha\beta}(\omega_{mn})  (A_\beta)_{nm} (A_\alpha)^*_{pq}.
$$

### Equivalence to rate equations {#sec:theory:rate-equations}

A further, relevant simplification of the master equation can be made if the spectrum of the system Hamiltonian $H_S$ is nondegenerate. In this case, the Kronecker delta functions appearing in the definitions of the damping coefficients simplify into $\delta_{\omega_{mn}} \rightarrow \delta_{mn}$. By looking at the equation for the populations of the density matrix, $\rho_{nn} = \langle n | \rho_S | n \rangle$, we obtain

$$
\dot{\rho}_{nn} = \sum_m \gamma_{nm,nm} \rho_{mm} - \left[ \sum_m \gamma_{mn,mn} \right] \rho_{nn},
$$
{#eq:theory:pauli-rate-equation}

i.e., populations only couple to the populations, while the coherences ($\rho_{nm}$ with $n \neq m$) decay to zero. The coefficients $\gamma_{nm,nm}$ correspond to the transition rates from state $m$ to $n$, and are equivalent to those obtained by Fermi's golden rule:

$$
\gamma_{nm,nm} = w_{n \leftarrow m} = \sum_{\alpha\beta} \gamma_{\alpha\beta} (\omega_{mn}) (A_\beta)_{nm} (A_\alpha)_{nm}^*.
$$

Equation (-@eq:theory:pauli-rate-equation) is known as Pauli master equation. It allows to compute the evolution of an open system by taking into account only the populations of the eigenstates (which are $N$ for a $N$-dimensional Hilbert space of the system) and not the full density matrix (composed of $N^2$ elements). Therefore, whenever the conditions to use it apply, it is of great numerical advantage for a large system.

<!-- ### Discussion on the assumptions made -->



<!-- We now consider the formal expression for $\rho_I(t + \delta t)$:

$$
    \rho_I (t+\delta t) = \rho_I(t) - i \int_t^{t+\delta t} d\tau [H_I(\tau), \rho_I(\tau)],
$$

which gives the density matrix after a short time $\delta t$, when
$\rho_I(t)$ is assumed to be known. We now iterate for two times and denote
the variation
$\delta \rho_I(t) = \rho_I(t+\delta t) - \rho_I(t)$. We obtain

$$
    \begin{split}
    \delta \rho_I (t) = &-i \int_t^{t+\delta t} d\tau [H_I(\tau), \rho_I(t)] -
    \int_t^{t+\delta t} d\tau \int_t^{\tau} d\tau' [H_I(\tau), [H_I(\tau'),
    \rho_I(t)]] \\
    &+i \int_t^{t+\delta t} d\tau \int_t^{\tau} d\tau' \int_t^{\tau'} d\tau''
    [H_I(\tau), [H_I(\tau'), [H_I(\tau''), \rho_I(\tau'')]]].
    \end{split}
$$
{#eq:theory:liouville-interaction-picture-third-order}

Equation (-@eq:theory:liouville-interaction-picture-third-order) is formally exact. In the time integrals, the time ordering follows $t + \delta t \geq \tau \geq \tau' \geq \tau'' \geq t$.  -->

## Solution of the master equation

I briefly present here a number of analytical and numerical tools than can be used to solve the master equation, some of which are employed in the following Chapters of the Thesis. For simplicity, I will restrict to the case of a master equation with time-independent Liouvillian, $\dot{\rho}(t)= \mathcal{L} \rho (t)$, which can be put in the diagonal form of Eq. (-@eq:theory:lindblad-equation-diagonal). I will drop the subscript $S$ and reference to the system density matrix as $\rho$.


### Analytical techniques

#### Equation of motion

The master equation for the density matrix can directly used to compute the dynamical evolution of the relevant expectation values of the system, by simply exploiting the relation $\langle O_i \rangle = \Tr [ O_i \rho]$, where $\{O_i\}$ is a family of system observables. Usually, this leads to a system of first-order linear differential equations involving the observables. If the system cannot be closed, further approximations may be required. It is worth stressing that for infinitely-dimensional Hilbert spaces (e.g., systems involving harmonic oscillators), this method is extremely useful. We write

$$
\begin{split}
\langle \dot{O}_i \rangle &= \Tr [O_i \dot {\rho}] = \Tr [O_i \mathcal{L}\rho] \\
    &= -i \langle [O_i, H + H_{LS}] \rangle + \sum_\alpha \gamma_\alpha \left[ \langle A_\alpha^\dagger O_i A_\alpha \rangle - \frac{1}{2} \langle \{O_i, A^\dagger_\alpha A_\alpha \} \rangle \right].
\end{split}
$$

For a finite-dimensional system, there exists a finite number of linearly independent observables $\{O_i \rangle\}$, hence one will end up with a linear system of the form

$$
\langle \dot{O}_i \rangle = \sum_j M_{ij} \langle O_j \rangle,
$$
{#eq:theory:exp-val-system}

where $M_{ij}$ are the entries of a matrix $\mathbf{M}$.


#### Adjoint master equation and quantum regression theorem

Similarly to the Heisenberg evolution of closed systems, it is possible to derive a master equation in which the state is kept constant while the system operators are left to evolve. This might be simply implied from Eq. (-@eq:theory:exp-val-system), yielding 

$$
\begin{split}
\dot{O}_i(t) &= -i[O_i(t), H+ H_{LS}] + \sum_\alpha \gamma_\alpha \left[A^\dagger_\alpha O_i(t) A_\alpha - \frac{1}{2} \{O_i(t), A^\dagger_\alpha A_\alpha \} \right] \equiv \mathcal{L}^\dagger O_i(t) \\
&= \sum_j M_{ij} O_j (t)
\end{split}
$$
{#eq:theory:adjoint-equation}

which is known as adjoint master equation with the hermitian conjugate of the Liouvillian, $\mathcal{L}^\dagger$ [@Breuer2002]. The time dependence in the operators implies we are working in the Heisenberg picture. Formally, this follows from the relation

$$
\langle O_i \rangle (t) =  \Tr_S [ O_i \rho(t) ] = \Tr_S [O_i e^{\mathcal{L} t} \rho] = \Tr_S [ (e^{\mathcal{L}^\dagger t} O_i) \rho],
$$

from which we define $O_i (t) = e^{\mathcal{L}^\dagger t} O_i$. Differentiating with respect to time yields Eq. (-@eq:theory:adjoint-equation).
The adjoint master equation becomes useful to calculate two-point correlation functions at different times, by knowing only single-point correlation functions (expectation values), a result known as quantum regression theorem [@Carmichael1999;@Gardiner2004]. From

$$
\frac{d}{d\tau} O_i (t + \tau) = \frac{d}{d\tau} e^{\mathcal{L}^\dagger (t + \tau)}O_i = \mathcal{L}^\dagger O_i(t + \tau) = \sum_j M_{ij}O_j (t + \tau),
$$

it follows

$$
\frac{d}{d\tau} \langle O_i (t + \tau) O_j (t) \rangle = \sum_j M_{ij} \langle O_i (t + \tau) O_j (t) \rangle,
$$

i.e., the same matrix coefficients stemming from the dynamical equations of single operators can be used. Two-time correlation functions are extremely useful when working with open systems, as they provide more statistical information than the expectation values. For example, the frequency spectrum of photons emitted from a cavity is given by the Fourier transform of the correlation function

$$
g^{(1)}(\tau) = \langle a^\dagger (t + \tau) a (t) \rangle,
$$

where $a$ is the bosonic mode of the cavity.

### Laplace transform

The Laplace transform method [@Schiff1999] is especially useful if one wants to compute the stationary state of a system, satisfying $\dot{\rho} = 0$. The Laplace transform is defined as

$$
\tilde{\rho}(z) = \int_0^\infty dt e^{-zt} \rho (t),
$$

from which one deducts

$$
\tilde{\rho}(z) = \frac{1}{z \mathds{1} - \mathcal{L}} \rho(0).
$$

with the identity matrix in Liouville space, $\mathds{1}$. The calculation of the inverse of the matrix $z \mathds{1} - \mathcal{L}$ is much simpler that exponentiating $\mathcal{L}t$, but one is then left with the problem of inverting the Laplace-transformed $\tilde{\rho}(z)$ back to real time, which is usually achieved through contour integration in complex space after identification of the poles of $(z\mathds{1} - \mathcal{L})^{-1}$.
The stationary state of the system is formally obtained by taking the limit 

$$
\rho_\mathrm{st} = \lim_{t\rightarrow \infty} \rho(t) = \lim_{z\rightarrow 0} z \tilde{\rho}(z).
$$

### Numerical techniques

#### Numerical integration

To propagate numerically a master equation $\dot{\rho} = \mathcal{L} \rho$ the general idea is time discretization into small time steps $\Delta t$. Explicit schemes rely on the recipe

$$
\frac{\rho(t + \Delta t) - \rho(t)}{\Delta t} = \mathcal{L} \rho(t),
$$

where $\mathcal{L}$ is applied to the "known" $\rho(t)$, while implicit methods are, e.g.,

$$
\frac{\rho(t + \Delta t) - \rho(t)}{\Delta t} = \mathcal{L} \frac{1}{2} [ \rho(t) + \rho(t + \Delta t)].
$$

Other solution schemes exist, such as semi-implicit methods and splitting strategies. As an example of an explicit method widely used for its stability, I briefly outline the fourth-order Runge-Kutta method [@Ascher1998]. For a choice of step size $\Delta t$, the density matrix at time $t + \Delta t$ ($\rho_{n+1}$) is calculated from the density matrix at time $t$ ($\rho_n$) according to

$$
\rho_{n+1} = \rho_n + \frac{1}{6} \Delta t (k_1 + 2k_2 + 2k_3 + k_4),
$$

with

$$
\begin{aligned}
k_1 &= \mathcal{L} \rho_n, \\
k_2 &= \mathcal{L} \left(\rho_n + \frac{1}{2} \Delta t k_1 \right), \\
k_3 &= \mathcal{L} \left(\rho_n + \frac{1}{2}\Delta t k_2\right), \\
k_4 &= \mathcal{L} \left(\rho_n + \Delta t k_3\right).
\end{aligned}
$$

The method applies the Liouvillian four times for each time step, and is indeed equivalent to the fourth-order expansion 

$$
\rho_{n+1} = \left[ \mathds{1} + \Delta t \mathcal{L} + \frac{1}{2!}(\Delta t)^2 \mathcal{L}^2 + \frac{1}{3!}(\Delta t)^3 \mathcal{L}^3 + \frac{1}{4!}(\Delta t)^4 \mathcal{L}^4 \right]\rho_n + \mathcal{O}\{(\Delta t)^5\}.
$$


#### Quantum trajectories (piecewise deterministic processes)

The basic idea underlying the quantum trajectories theory consists of rewriting
the master equation as an ensemble average over the individual stochastic
trajectories of an initial pure state of the system [@Dalibard1992;@Molmer1993;@Carmichael2008]. Specifically, an effective *deterministic* evolution is interrupted by stochastic jumps of the system, due to the interaction with the environment. This means that the master equation can be represented by a piecewise deterministic process (PDP) [@Breuer2002]. It is
possible to achieve this result in several ways, which correspond to different stochastic
*unravelings* of the master equation. The choice of the unraveling corresponds
to a set of probabilistic decisions over time, which defines the evolution of
the system (the "trajectory") and is inevitably lost in the ensemble-average
representation of the density matrix. The numerical advantage of the quantum trajectories method lies in the fact that it does not require to store the full density matrix (consisting of $N^2$ elements for a $N$-dimensional Hilbert space), but consists instead in realizing a large number $M \gg 1$ of trajectories starting from a pure state, which requires only $N$ complex observables to be stored. For larger systems, $M \ll N^2$ trajectories are in general sufficient to guarantee a good ensemble statistics.

A possible formulation of the quantum trajectories method can be devised in
terms of photon counting from a cavity using a detector.
This corresponds to the "quantum jump" approach to damping, which I will consider here. 
An equivalent derivation of a quantum trajectory follows from the description 
of the conditional evolution of a system under a continuous
measurement, conditioned on a stochastic record of measurements, and is generally referred
to as the "diffusive limit" of the master equation unraveling [@Walls2008;@Wiseman2010]. A common application of this limit is the
homodyne or heterodyne detection of the field leaving the cavity.

The idea behind the quantum jump approach is to consider the following nonlinear equation for a state $|\psi(t)\rangle$ of the system:

$$
| \dot{\psi} \rangle = -i \left(H - \frac{i}{2} \sum_\alpha \gamma_\alpha A^\dagger_\alpha A_\alpha \right)|\psi \rangle + \frac{1}{2} \left( \sum_\alpha \gamma_\alpha \langle \psi | A^\dagger_\alpha A_\alpha | \psi \rangle \right) | \psi \rangle.
$$

Its solution is given by [@Walls2008]

$$
| \psi(t) \rangle = \frac{e^{-i  H_\mathrm{eff} t} }{\sqrt{\langle \psi(0) | e^{i H_\mathrm{eff}^\dagger t} e^{-i H_\mathrm{eff} t} |\psi(0)\rangle}} |\psi(0) \rangle,
$$
{#eq:theory:quantum-trajectories-deterministic}

where I have introduced the nonhermitian Hamiltonian:

$$
H_\mathrm{eff} = H -\frac{i}{2} \sum_\alpha \gamma_\alpha A^\dagger_\alpha A_\alpha.
$$

Clearly, Eq. (-@eq:theory:quantum-trajectories-deterministic) gives a deterministic evolution. To reproduce the correct Lindblad dynamics, it is necessary to interrupt the deterministic dynamics with stochastic events (jumps). Defining the total probability of a jump to occur within the time interval $\delta t$,

$$
p_\mathrm{jump} = \delta t \sum_\alpha \gamma_\alpha,
$$

one needs to decide randomly which jump has occurred by updating the state according to

$$
|\tilde{\psi}(t + \delta t) \rangle  = A_\alpha |\psi(t) \rangle,
$$

to be normalized. The probability for the $\alpha$-jump to occur is given by

$$
p_\alpha = \frac{\gamma_\alpha \langle \psi(t) | A^\dagger_\alpha A_\alpha | \psi(t)  \rangle}{p_\mathrm{jump}/\delta t}.
$$

![Numerical scheme for the propagation of one time step in the quantum trajectories (PDP) method.](figures/theory-quantum-trajectories.pdf){#fig:theory:quantum-trajectories}

If the jump has not occurred, the state is propagated according to Eq. (-@eq:theory:quantum-trajectories-deterministic). The numerical scheme is summarized in @fig:theory:quantum-trajectories, for the simplified case of a single jump operator $a$.  At time $t$, $p_\mathrm{jump}$ is computed. Then, a random number $r$ is generated. If $p_\mathrm{jump} < r$, a jump is assumed to have taken place. The jump operator $a$ is applied to the state, and it is then normalized, giving $|\psi(t + \delta t)\rangle$. Conversely, if $p_\mathrm{jump} \geq 1$, no jump has occurred and the nonunitary effective evolution with $H_\mathrm{eff}$ is assumed. 

Formally, the procedure described above is described by a stochastic Schrödinger equation, given by [@Wiseman2010]

$$
|d \psi \rangle =\left( -i H_\mathrm{eff} + \frac{1}{2} \sum_\alpha \gamma_\alpha \langle \psi | A^\dagger_\alpha A_\alpha | \psi \rangle \right) | \psi \rangle dt + \sum_\alpha \left(  \frac{A_\alpha |\psi}{\sqrt{\langle \psi|A^\dagger_\alpha A_\alpha |\psi\rangle}} - |\psi\right) d N_\alpha.
$$
{#eq:theory:sse}

The $d N_\alpha$ are Poisson increments, satisfying

$$
d N_\alpha d N_\beta = \delta_{\alpha \beta} d N_\alpha, \quad \mathrm{E} [d N_\alpha] = \delta t \gamma_\alpha \langle \psi | A^\dagger_\alpha A_\alpha| \psi\rangle,
$$

where $\mathrm{E}[\cdot]$ corresponds to the classical expectation value. 
It is possible to show that the ensemble average $\mathrm{E} [|\psi \rangle \langle \psi| ]$ evolved with Eq. (-@eq:theory:sse), satisfies the Lindblad master equation (-@eq:theory:lindblad-equation-diagonal) [@Breuer2002;@Wiseman2010].




#### Steady-state solution

Of particular interest for this work is a method to look for the stationary state of the system $\rho_\mathrm{st}$, i.e., the state reached for $t \rightarrow \infty$. Instead of integrating numerically the master equation up to long times, one may focus on the steady equation

$$
\dot{\rho}  = \mathcal{L}\rho  \overset{!}{=} 0.
$$

In this way, we modify the problem into a linear algebra one, which consists in finding the null space (or kernel) of the Liouvillian operator $\mathcal{L}$, i.e., the right eigenvector associated to the zero eigenvalue of the matrix $\mathcal{L}$. Being a standard problem, many methods exist to solve it [@Saad2003]. The direct method, consisting in simply computing the inverse $\mathcal{L}^{-1}$, cannot be used as $\mathcal{L}$ is singular, and it would be a formidable numerical problem for a large system. An alternative is given by the Arnoldi iteration scheme to compute the zero-eigenvector of $\mathcal{L}$ by finding an orthonormal basis for the Krylov subspace from some initial guess state $\rho_0$. The Arnoldi scheme is the basis for the GMRES (generalized minimal residual) method, which is an efficient iterative procedure to compute $\rho_\mathrm{st}$ [@Saad1986].

Note that these methods are valid to solve in general the linear system $Ax = b$, and can thus be used also to find the vector of steady populations which satisfies the stationary Pauli rate equations presented in @sec:theory:rate-equations. In this case, one needs to find the vector $\rho_\mathrm{st}$ satisfying

$$
\sum_m w_{n \leftarrow m} \rho_{mm} - \sum_m w_{m \leftarrow n} \rho_{nn} = 0,
$$

i.e., the kernel of the reduced Liouvillian $\mathcal{W}$ with entries 

$$
W_{nm} = w_{n \leftarrow m} - \delta_{nm} \sum_p w_{p \leftarrow m}.
$$


<!--



## Quantum trajectories method

The basic idea underlying the quantum trajectories theory consists of rewriting
the master equation as an ensemble average over the individual stochastic
trajectories of an initial pure state of the system. It is
possible to achieve this result in several ways, which correspond to different stochastic
*unravelings* of the master equation. The choice of the unraveling corresponds
to a set of probabilistic decisions over time, which defines the evolution of
the system (the ``trajectory'') and is inevitably lost in the ensemble-average
representation of the density matrix. 

A possible formulation of the quantum trajectories method can be devised in
terms of photon counting from a cavity using a detector.
This corresponds to the "quantum jump" approach to damping. 
An equivalent derivation of a quantum trajectory follows from the description 
of the conditional evolution of a system under a continuous
measurement, conditioned on a stochastic record of measurements, and is generally referred
to as the "diffusive limit" of the master equation unraveling
\cite{Walls2008,Wiseman2010}. A common application of this limit is the
homodyne or heterodyne detection of the field leaving the cavity.


### Stochastic master equation

Consider for simplicity a cavity with decay rate $\kappa$, coupled to a zero-temperature
bath. The master equation for the intracavity field
is given by
\begin{equation}
    \label{eq:theory:qt-cavity-me}
    \dot{\rho} = -i[H,\rho] + \frac{\kappa}{2}(2 a \rho
    a^\dagger - \rho a^\dagger a - a^\dagger a \rho),
\end{equation}
Here, we recall that 
the state of the cavity $\rho$ is
obtained by tracing over the environment, $\rho = \Tr_E[\rho_\text{tot}]$.
The master equation gives, therefore, a local description of the system. The
trajectories method aims at disentangling system and reservoir, by keeping
count of the photons lost from the cavity and detected from a recording
apparatus, thus giving a nonlocal description of system and environment.
Integrating Eq. \eqref{eq:theory:qt-cavity-me} in a small time interval $\delta
t$ gives
\begin{equation}
    \label{eq:theory:qt-cavity-me-int}
    \rho(t + \delta t) = \{ \rho(t) - i[H, \rho(t)] \delta t - \frac{\kappa}{2}[a^\dagger a \rho(t) +
    \rho(t)a^\dagger a] \delta t \} + \kappa a \rho (t) a^\dagger \delta t.
\end{equation}
Suppose the cavity is at time $t$ in a Fock state with $n$ photons, represented by the density
matrix $\rho(t) = |n \rangle \langle n|$. Then, at time $t + \delta t$, the
state of the system is
\begin{equation}
    \label{eq:theory:qt-cavity-me-int-2}
    \rho(t + \delta t) = (1 - \kappa n \delta t) |n\rangle \langle n| + \kappa
    n \delta t |n - 1\rangle \langle n -1 |.
\end{equation}
It is clear that Eq. \eqref{eq:theory:qt-cavity-me-int-2} gives the
*unconditional* state of the system under two possible events: The first term
describes a situation where, after a time $\delta t$, no photons are lost from
the cavity. The probability for this is $1 - \kappa n \delta t$. Conversely,
the second term tells us that a photon is subtracted with probability $\kappa n
\delta t$, leaving the cavity in the Fock state $|n-1\rangle$. We will call
this event a
*quantum jump*, and it represents a sudden change in the knowledge of the state.
Equation \eqref{eq:theory:qt-cavity-me-int-2}
represents a statistical mixture. However, we want to gain information on the
state of the system *conditioned* upon a measurement of the photocounter, which
can, at each time interval $\delta t$, detect a jump or no jump. The latter
situation is
referred to as a null measurement. The state of the system conditioned upon a
null measurement in the time interval $[t, t+\delta t$)  is obtained by considering only the first term in Eq.
\eqref{eq:theory:qt-cavity-me-int-2}:
\begin{equation}
    \label{eq:theory:qt-cavity-nojump}
    \rho_{\text{no-jump}}(t + \delta t) = \frac{\rho(t) - i [H, \rho(t)] \delta
    t -  \frac{\kappa}{2}[a^\dagger a \rho(t) +
    \rho(t)a^\dagger a] \delta t }{\Tr \{ \rho(t) - i [H, \rho(t)] \delta
    t -  \frac{\kappa}{2}[a^\dagger a \rho(t) +
    \rho(t)a^\dagger a] \delta t \}},
\end{equation}
where the state has been normalized. Similarly, we could write the state
conditioned upon a jump.  We now wish to write a general, stochastic master equation for the 
density matrix, which reduces to Eq. \eqref{eq:theory:qt-cavity-me} when taking
a classical ensemble average. First, assume $\delta
t \ll \kappa^{-1}$, such that the jump probability $p_\text{jump}$ remains
small. The jump process is modeled by a classical Poisson point process
$dN(t)$ \cite{Karatzas1998}, with the properties
\begin{align}
    dN(t)^2 &= dN(t), \label{eq:theory:qt-poisson} \\
    \text{E} [ dN(t) ] &= \kappa \langle a^\dagger a \rangle \delta t.
    \label{eq:theory:qt-poisson-expect}
\end{align}
Equation \eqref{eq:theory:qt-poisson} refers to the fact that $dN(t)$ is either
0 or 1, since the jump probability is really small, and only no count or one
count is assumed possible within time $\delta t$. Equation \eqref{eq:theory:qt-poisson-expect} 
denotes the classical expectation value of the Poisson process, thus the jump
probability. Starting from a pure state $|\psi(t) \rangle$, if $dN(t) = 1$ a
jump as occurred and state changes to
\begin{equation}
    \label{eq:theory:qt-state-jump}
    |\psi_\text{jump}(t + \delta t) \rangle = \frac{a |\psi(t)
    \rangle}{\sqrt{\langle a^\dagger a \rangle}_t},
\end{equation}
where $\langle \cdot \rangle_t$ denotes the expectation value on $|\psi(t)
\rangle$. For no detection ($dN(t) = 0$), we have
\begin{align}
   |\psi_{\text{no-jump}}(t + \delta t) \rangle &= \frac{a |\psi(t)
    \rangle}{\sqrt{\langle a^\dagger a
    \rangle}_t} \nonumber \\
    &\approx \left\{ 1 - i H \delta t - \frac{\kappa}{2} (a^\dagger a + \langle
    a^\dagger a \rangle_t) \delta t \right\} |\psi(t) \rangle.
    \label{eq:theory:qt-state-nojump}
\end{align}

We have expanded the denominator to first order, as $\kappa \delta t \ll 1$.
Taking into account the stochastic term $dN(t)$ we arrive at the nonlinear stochastic
Schrödinger equation \cite{Wiseman2010}:
\begin{equation}
    \label{eq:qt-stochastic-se}
    d |\psi(t) \rangle = \left [ d N(t) \left ( \right) + \kappa \delta t
    \left( \right) \right] |\psi(t)\rangle.
\end{equation}
In writing Eq. \eqref{eq:qt-stochastic-se} we have neglected a term $dN(t)
\delta t$ which is $o(t)$. The stochastic equation has the simple
interpretation of a...


### Simulating quantum trajectories

![Numerical scheme for the evolution of one time step in the quantum trajectories method. \label{fig:theory:quantum-trajectories}](quantum-trajectories.pdf)

For a Hilbert space of dimension $N$, a density matrix $\rho$ representing a mixed
state has $N^2$ elements, and the total Liouvillian $\mathcal{L}$ describing its dynamics 
has $N^4$ elements. Although usually many of the elements of $\mathcal{L}$ are zero,
making computations tractable with sparse-matrices algorithms, the full
numerical solution of the quantum master equation can become extremely
demanding for large systems. On the other hand, a pure state requires only $N$
elements to be stored, and thus solving the stochastic Schrödinger equation is a much
less challenging task.

The idea behind the numerical simulation of a quantum trajectory follows a
Monte Carlo scheme and is summarized
in Fig. \ref{fig:theory:quantum-trajectories}. Suppose the system is initially 
in the state $|\psi(t)\rangle$, and we want to
compute the state after a short time interval $\delta t$, with $\delta t \ll
\kappa^{-1}$. Within this time $\delta t$, the evolution of the state is
determined by the nonhermitian effective Hamiltonian
\begin{equation}
    \label{eq:theory:qt-nonhermitian-hamiltonian}
    H_\text{eff} = H - i\frac{\kappa}{2}a^\dagger a,
\end{equation}
which gives the unnormalized state
\begin{equation}
    \label{eq:theory:qt-nonhermitian-evolved}
    |\tilde{\psi}(t + \delta t) \rangle = e^{- i H_\text{eff} \delta t}
    |\psi(t)\rangle.
\end{equation}
The increment can be calculated using standard numerical techniques. The
effective Hamiltonian $H_\text{eff}$ causes the state norm to decay over
time, as it keeps track only of the evolution without jumps.
It is easy to see that the norm of this state is $p_{\text{no-jump}} = 1 -
p_{\text{jump}}$ with 
\begin{equation}
    \label{eq:theory:qt-pjump}
    p_\text{jump} = \kappa \delta t \langle \psi(t) | a^\dagger a |\psi(t)
    \rangle.
\end{equation}
By keeping both $\gamma$ and $\delta t$ very small, we ensure that
$p_\text{jump} \ll 1$. We now compare this probability to a random number
$r$, uniformly distributed between 0 and 1. Most of the times, $p_\text{jump}
< r$, and we update the state of the system with the (normalized) state of Eq.
\eqref{eq:theory:qt-nonhermitian-evolved}, i.e.:
\begin{equation}
    \label{eq:theory:qt-nojump}
    |\psi(t + \delta t) \rangle = \frac{e^{- i H_\text{eff} \delta t}
    |\psi(t)\rangle}{\sqrt{1 - p_\text{jump}}}.
\end{equation}
If, however, $p_\text{jump} \geq r$, we assume that a jump has taken place.
The new, normalized state at time $t + \delta t$ is therefore
\begin{equation}
    \label{eq:theory:qt-jump}
    |\psi(t + \delta t) \rangle = \frac{ a
    |\psi(t)\rangle}{\sqrt{\frac{p_\text{jump}}{\kappa \delta t}}}.
\end{equation}
By keeping track of the times at which the jumps occur, we can reconstruct the
history of jump records $h[t]$. To have access to the ensemble average, we run the
simulation $N_\text{traj}$ times, always starting with the same initial
state, thus obtaining the mixed state 
\begin{equation}
    \label{eq:theory:qt-ensemble}
    \tilde{\rho}(t) = \frac{1}{N_\text{traj}} \sum_{k=1}^{N_\text{traj}}
    |\psi_k (t) \rangle.
\end{equation}
For $N_\text{traj} \gg 1$, one has $\tilde{\rho}(t) \approx \rho(t)$, where
$\rho(t)$ is the master equation solution, with an error which scales as
$1/\sqrt{N_\text{traj}}$. Although requiring $N_\text{traj} \gg 1$ to get 
good statistics, the computational scheme
scales only quadratically in the Hilbert space dimension $N$ and is, therefore,
of crucial advantage for large systems.

It is worth noticing here that the quantum
jumps are not really discontinuous in time: similarly to the master equation,
the evolution is assumed to be coarse-grained, such that the time step $\delta
t$ used in the computation remains much larger than the typical transition time
between the states interested. However, we require also that $\delta t$ is
small such that the jump probability for each time step is less than unity. The
separation of timescales follows again from the Markov approximation leading to
the master equation, thus enforcing the equivalence between the two approaches
when one averages over the trajectories.

It is easy to include several jump processes in the method, by using different
jump counters.
Finally, the scheme presented here can be generalized to the finite temperature
case, by performing a double average: the initial state $\psi(t)$ is each time
sampled from a thermal equilibrium distribution, and many trajectories
are simulated. The process is then repeated for many states, and the resulting mixed
density matrix is obtained by averaging with the Gibbs weights
of the initial states.




<!--
## Connection to thermodynamics
* Pioneering works and reviews
  \cite{Geusic1967,Alicki1979,Kosloff2013,Benenti2017};
* Other works \cite{Boukobza2006,Boukobza2006a,Boukobza2007};
-->

<!-- ## Applications of the master equation

### Damped harmonic oscillator
As a first application of the QME, we consider a single harmonic mode of a
resonator of frequency $\omega_0$, weakly coupled to a thermal bath. This is a standard way to treat a
damped electromagnetic or mechanical cavity with a finite quality factor $Q$.
The thermal bath is described by a bosonic reservoir of harmonic oscillators with a
quasicontinuous spectrum of frequencies $\omega_k$, linearly coupled to the cavity mode.
The total Hamiltonian is therefore written as follows:
\begin{equation}
    \label{eq:theory:damped-oscillator-hamiltonian}
H = \underbrace{\omega_0 b^\dagger b}_{H_{\text{osc}}} + \underbrace{\sum_k
\omega_k b^\dagger_k b_k}_{H_{\text{bath}}} + \underbrace{\sum_k g_k (b +
b^\dagger)(b_k + b^\dagger_k)}_{H_I}.
\end{equation}

### Relaxation and decoherence of a qubit
* Two-level system coupled to environment
	- Dephasing and relaxations -> useful for Chap. \ref{ch:jjcavity} and Chap. \ref{ch:spinboson}
	- Rabi oscillations: $T_1$ measurements;
	- Ramsey interferometry: $T_2$ and $T_2^*$ measurements.

### Transport through a double quantum dot
* Tunneling between dots and leads (discuss especially the form of the fermionic
  correlation functions \cite{Schoeller2000,Bevilacqua2016} and why Born-Markov
  is valid)
* Double dot coupled to leads; -->