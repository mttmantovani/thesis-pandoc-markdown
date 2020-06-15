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
Gorini-Kossakowski-Sudarshan-Lindblad (GKLS or, more commonly, Lindblad) equation [@Gorini1976;@Lindblad1976]. Here, I will instead focus on a derivation of the master equation which follows a well-developed microscopic approach and leads to the Bloch-Redfield equation [@Redfield1957;@Carmichael1993;@Breuer2002]. I will discuss critically the validity regime of the Bloch-Redfield equation, showing the condition under which it has a Lindblad form. Later, I will apply it to a class of pedagogical models that will be instrumental to the analysis of the systems presented in the following Chapters. 

In addition to the Bloch-Redfield equation, I will also discuss the quantum trajectories method [@Scully1997;@Carmichael2008;@Walls2008;@Wiseman2010], or Monte Carlo wavefunction method [@Dalibard1992;@Molmer1993]. Besides being of great numerical advantage when dealing with the dynamics of large systems (made up by several few-level systems or highly excited harmonic resonators), it also offers a purely stochastic representation of the evolution of a system in contact to an environment or to a measurement apparatus, which goes beyond the ensemble-average representation given by the density-matrix theory.

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

The truncation of the Lioville-von Neumann equation to second order requires a justification, which will be made clear later.

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

Now, we split the function $\Gamma_{\alpha \beta} (\omega)$ into

\begin{align}
 \Gamma_{\alpha \beta}(\omega) &=\frac{1}{2} \gamma_{\alpha \beta}(\omega)+\frac{1}{2} \sigma_{\alpha \beta}(\omega), \\ \Gamma_{\beta \alpha}^{*}(\omega) &=\frac{1}{2} \gamma_{\alpha \beta}(\omega)-\frac{1}{2} \sigma_{\alpha \beta}(\omega),
\end{align}

where

$$
\begin{array}{l}\gamma_{\alpha \beta}(\omega)=\Gamma_{\alpha \beta}(\omega)+\Gamma_{\beta \alpha}^{*}(\omega)=\int_{-\infty}^{+\infty} C_{\alpha \beta}(\tau) e^{+i \omega \tau} d \tau, \\ \sigma_{\alpha \beta}(\omega)=\Gamma_{\alpha \beta}(\omega)-\Gamma_{\beta \alpha}^{*}(\omega)=\int_{-\infty}^{+\infty} C_{\alpha \beta}(\tau) \operatorname{sgn}(\tau) \mathrm{e}^{+\mathrm{i} \omega \tau} \mathrm{d} \tau.\end{array}
$$

### Master equation in the energy eigenbasis

A sometimes useful representation of the master equation (without employing the spectral decomposition) is found by projecting Eq. (-@eq:theory:bloch-redfield) in the energy eigenbasis $|n\rangle$, which yields

$$
\dot{\rho}_S = \sum_{\alpha\beta} \sum_{nmpq} \Gamma_{\alpha\beta}(\omega_{mn})e^{i[\omega_{mn} - \omega_{qp}]t} (A_\beta)_{nm} (A_\alpha)_{pq}^* \left\{ L_{nm} \rho_S (t) L_{pq}^\dagger - L^\dagger_{pq} L_{nm} \rho_S(t) \right\} + \mathrm{H.c.}
$$

We have introduced the shorthands $(A_\alpha)_{nm} = \langle n | A_\alpha | m \rangle$ and $L_{nm} = |n \rangle \langle m|$ (the latter is the Lindblad jump operator in the energy representation). The secular approximation amounts here to neglecting all terms where $\omega_{nm} \neq \omega_{qp}$.

### Equivalence to rate equations



### Discussion on the assumptions made



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

### Analytical techniques

#### Equation of motion

#### Quantum regression theorem

### Numerical techniques

#### Steady-state solution

#### Numerical integration

#### Quantum trajectories (piecewise deterministic processes)

![Numerical scheme for the propagation of one time step in the quantum trajectories (PDP) method.](figures/theory-quantum-trajectories.pdf){#fig:theory:quantum-trajectories}

<!-- 

### Born approximation

Weak-coupling approximation: truncation to second order.

:   We now perform the Born (or weak-coupling) approximation
    \cite{Cohen-Tannoudji2004}. We assume that the system-environment coupling is
    *weak*, such that we can neglect the third-order term in
    Eq. \eqref{eq:theory:liouville-interaction-picture-third-order}. The condition
    under which this approximation is valid will be clarified later. We
    remain with
    \begin{equation}
        \label{eq:theory:liouville-interaction-picture-second-order}
        \delta \rho_I (t) = -i \int_t^{t+\delta t} d\tau [H_I(\tau), \rho_I(t)] -
        \int_t^{t+\delta t} d\tau \int_t^{\tau} d\tau' [H_I(\tau), [H_I(\tau'),
        \rho_I(t)]],
    \end{equation}
    which is valid to *second order* in the interaction Hamiltonian $H_I(t)$. 
    As we are interested in the dynamics of the *reduced* density matrix
    $\rho_{S,I}(t)$, we trace over the Hilbert space of the environment
    $\mathcal{H}_E$, obtaining
    \begin{equation}
        \label{eq:theory:system-liouville-interaction-picture-second-order}
        \delta \rho_{S,I} (t) = -i \int_t^{t+\delta t} d\tau \Tr_E [H_I(\tau), \rho_I(t)] -
        \int_t^{t+\delta t} d\tau \int_t^{\tau} d\tau' \Tr_E [H_I(\tau), [H_I(\tau'),
        \rho_I(t)]].
    \end{equation}
    Equation \eqref{eq:theory:system-liouville-interaction-picture-second-order}
    still contains the *full* density matrix at time $t$. 

Large-reservoir approximation: timescale separation.

:   Next, we make a key assumption on the state of the environment, which is
    obtained by tracing over the system degrees of freedom, $\rho_{E,I} (t) = \Tr_S
    \rho_I(t)$: The environment has usually a very large number of degrees of
    freedom, and is weakly coupled to the system. Hence, we may assume that the
    state $\rho_{E,I}$ is not appreciably changed by the system, i.e., it is
    constant in the interaction picture:
    \begin{equation}
        \rho_{E,I}(t)
        \approx \rho_{E,I}(0).
    \end{equation}
    Furthermore, we consider that the $\rho_{E}$ is stationary, i.e., it commutes
    with the environment Hamiltonian, $[\rho_{E}, H_E] = 0$. Equivalently, we may
    say that the environment is in thermodynamic equilibrium. The consequence
    of this large-reservoir approximation is that one can identify a separation
    of timescales in the evolution of $\rho_{S,I}(t)$. We will call $\tau_B$
    the typical timescale during which environmental correlations in the
    reservoir exist. When this time has elapsed, the reservoir state has lost
    dependence on its initial state. A second timescale, the relaxation time $\tau_R$,
    characterizes the typical duration over which the state of the system
    varies appreciably because of system-environment interaction, and is
    roughly specified by
    \begin{equation}
        \frac{\delta \rho_{S,I}(t)}{\delta t} \approx \frac{1}{\tau_R}
        \rho_{S,I}(t).
    \end{equation}




### Markov approximation and Bloch-Redfield equation
Equation ... is still nonlocal in time. We now employ the Markov approximation
to make it local, allowing us to write to proceed further in the calculations.

### Secular approximation and Lindblad form of the master equation


Main references:

* Recent literature \cite{Timm2008,Hussein2014a,Cattaneo2019};
* QuTiP \mcite{qutip, *Johansson2012, *Johansson2013}.


### Full-counting statistics for the evaluation of particle and heat currents

* Early work in quantum optics \cite{Mandel1979,Cook1981};
* Applications to quantum transport \cite{Levitov1993,Levitov1996,Blanter2000,Bagrets2003,Hussein2014a}

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
* Double dot coupled to leads; --> -->