# Theory of double Cooper-pair tunneling and light emission mediated by a resonator {#sec:jjcavity}

*The results presented in this Chapter are published in Ref.
\onlinecite{Morley2019}. I am responsible for the numerical
calculations that led to the results of  [@fig:jjcavity:photon-charge;@fig:jjcavity:comparisons]; moreover, 
I have independently implemented the code written by the first
author, W. T. Morley, to reproduce the results presented in [@fig:jjcavity:fano;@fig:jjcavity:g2].*

## Introduction

The purpose of this Chapter is to present an example of how
resonators coupled to mesoscopic conductors can influence dramatically the
charge transport, even changing the effective charge that tunnels through.
In standard mesoscopic conductors or in scanning tunneling microscopy (STM), the
charge tunneling often comes together with photon emission, and this mechanism
can be enhanced through the coupling to a resonator
[@Berndt1991;@Berndt1993;@Holst1994;@Atay2004;@Schneider2010;@Hofheinz2011]. 
Conversely, it has been shown how a resonator can mediate correlated tunneling
of *two* electrons across a junction, in order to generate a cavity photon at a
energy higher than the single-electron one. This fenomenon is known as overbias
emission [@Schneider2010;@Peters2017;@Xu2014;@Kaasbjerg2015;@Xu2016]. 
Here, I consider a solid-state device consisting of a Josephson
junction (JJ) coupled to a LC-resonator, and theoretically analyze the correlated
tunneling of *two* Cooper pairs across the JJ [@Morley2019]. This higher-order tunneling
process is enhanced and becomes resonant when the bias voltage applied to the
junction is such that two Cooper pairs are required to tunnel to cooperatively
emit a photon at the cavity frequency. I will present an effective Hamiltonian
formulation of the problem and discuss how the double-Cooper-pair tunneling (DCPT) rate
can exceed the single-Cooper-pair one. Further, I will reveal how this
unconventional transport channel leads to an unusual form of photon blockade
[@Imamoglu1997;@Biella2015;@Vaneph2018;@Chang2014] that makes the system an
ideal single-photon source [@Grimm2019].

## Josephson juction coupled to a microwave cavity

![(a) Model circuit. A damped  microwave resonator of frequency $\omega_0 = (LC)^{-1/2}$ 
is in series to a Josephson junction with
applied bias voltage $V$. Low-frequency voltage fluctuations in the circuit
may be modeled by an impedance $Z(\omega)$ in series with the junction.
(b) When the bias voltage is such that $\omega_J = \omega_0/2$, photon emission
in the cavity is mediated by tunneling of two Cooper pairs.](./figures/jjcavity.pdf){#fig:jjcavity:system}

The system considered in Ref. \onlinecite{Morley2019} consists of a LC-resonator
of frequency $\omega_0 = 1/\sqrt{LC}$ (realized with a superconducting microwave
cavity [@Hofheinz2011;@Chen2014;@Cassidy2017] or
a lumped-element oscillator [@Rolland2019]) in series with a voltage-biased Josephson
junction. It is schematically pictured in @fig:jjcavity:system(a). We
assume that the resonator capacitance is much larger than the capacitance of the Josephson
junction [@Armour2013;@Meister2015]. The system can then be described by the
following time-dependent Hamiltonian [@Gramich2013;@Vool2017;@Lorch2018]:

$$   
    H(t) = \omega_0 a^\dagger a - E_J \cos [\omega_J t - \phi + \Delta_0 (a + a^\dagger)],
$$  {#eq:jjcavity:hamiltonian}

where $a$ is the bosonic annihilation operator for the cavity mode, $\omega_J =
2eV$ is the Jospehson frequency ($V$ is the voltage applied to the junction and
$e$ is the electron's charge), $E_J$ is the Josephson energy and $\Delta_0 =
\sqrt{2 e^2} (L/C)^{1/4}$ is the zero-point flux fluctuations amplitude of the
resonator, given in units of the flux quantum. Further, we have included the
phase $\phi$ which is conjugate to the Cooper-pair number $N$ passing through
the junction. The operator

$$
e^{i p \phi} = \sum_{N=0}^{\infty} |N + p \rangle \langle N |, \qquad p \in \mathbb{N},
$$ {#eq:jjcavity:phase-operator}

describes transfer of $p$ Cooper pairs across the junction.
The value of $\Delta_0$ is determined by the impedance of the resonator, and can
range from $\Delta_0 \ll 1$ [@Hofheinz2011;@Chen2014;@Cassidy2017] to values
close to unity [@Rolland2019].

Losses in the resonator due to the resistance give rise to a
quality factor $Q = \omega_0 / \kappa$, with $\kappa$ the decay rate. The
presence of low-frequency impedances in the circuit causes voltage fluctuations
that dephase the charge on the junction, with rate $\gamma_\phi$ [@Wang2017]. These
decoherence effects can be described by a standard Lindblad master equation
[@Gramich2013]

$$
    \dot{\rho} = -i[ H(t), \rho] + \frac{\kappa}{2}\mathcal{D}[a]\rho +
    \frac{\gamma_\phi}{2} \mathcal{D}[N]\rho,
$$ {#eq:jjcavity:me}

where we have assumed the limit of low temperature. Typically, $\gamma_\phi \ll
\kappa$, hence the phase $\phi$ can be treated as constant neglecting the
dephasing term [@Armour2013;@Wang2017].