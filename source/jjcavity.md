# Theory of double Cooper-pair tunneling and light emission mediated by a resonator {#sec:jjcavity}

*The results presented in this Chapter are published in Ref.
\onlinecite{Morley2019}. I am responsible for the numerical
calculations that led to the results of  [@fig:jjcavity:photon-charge; @fig:jjcavity:comparisons]; moreover, 
I have independently implemented the code written by the first
author, W. T. Morley, to reproduce the results presented in [@fig:jjcavity:fano; @fig:jjcavity:g2].*

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
$$ {#eq:jjcavity:hamiltonian}

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

## Effective Hamiltonian description

![(a) Time-averaged cavity photon number (solid blue line) and Cooper-pair tunneling
rate (solid orange line) as a function of $\delta = \omega_0 - 2\omega_J$,
calculated using the full Hamiltonian (-@eq:jjcavity:hamiltonian) and the master
equation (-@eq:jjcavity:me). Parameters: $\Delta_0 = 0.15,\ Q = 1500,\
E_J=0.5\omega_0,\ \gamma_\phi=0$. Notice that, close to the resonance,
$\frac{\Gamma_{\text{CP}}}{\kappa \langle n \rangle} \approx 2$, i.e., two Cooper pairs
tunnel for each emitted photon. (b) Schematic diagram of the energy levels for
$n$ cavity photons and $N$ tunneling Cooper pairs, at the resonance $2\omega_J
\approx \omega_0$, illustrating a few effective second-order processes included
in the effective Hamiltonian (-@eq:jjcavity:heff).](./figures/jjcavity-photon-charge.pdf){#fig:jjcavity:photon-charge}

The study of photon emission and tunneling of individual Cooper pairs across a
voltage-biased JJ has been extensively investigated, both in theory [@Padurariu2012;@Leppakangas2013;@Armour2013;@Gramich2013;@Armour2015;@Trif2015;@Dambach2015;@Leppakangas2015;@Souquet2016;@Leppakangas2018;@Lorch2018;@Mendes2019;@Arndt2019] and
experiment [@Hofheinz2011;@Chen2014;@Cassidy2017;@Westig2017;@Jebari2018;@Grimm2019;@Rolland2019], and it is resonantly enhanced when the Josephson frequency
$\omega_J$ matches the resonator frequency $\omega_0$.
The goal of this Section is instead to derive an effective Hamiltonian for the system,
by considering the resonance $\omega_0 \approx 2\omega_J$. As Fig.
\ref{fig:jjcavity:system}(b) illustrates, for an applied bias voltage $V = \frac{\omega_0}{4e}$ to
the junction, two Cooper pairs are expected to contribute to the emission of a
single photon at frequency $\omega_0$. We then take Hamiltonian Eq.
(-@eq:jjcavity:hamiltonian) and move to a frame rotating at frequency
$2\omega_J$, obtaining

$$
\tilde{H} = \delta a^\dagger a - \frac{\tilde{E}_J}{2} \sum_{q=0}^{\infty} [O_q
e^{i(2q+1)\omega_J t} + \text{H.c.}],
$$ {#eq:jjcavity:hamiltonian-rot-frame}

with the detuning $\delta = \omega_0 - 2\omega_J$. We have defined the rescaled
Josephson energy $\tilde{E}_J = E_J e^{-\Delta_0^2/2}$ and the operator

$$
O_q =\ :\ i^q (a^\dagger)^{q+1} e^{-i\phi} \frac{J_q(2\Delta_0 \sqrt{n})}{n^{q/2}} + (-i)^{q+1}
(a^\dagger)^{(q+1)} e^{i\phi} \frac{J_{q+1}(2 \Delta_0 \sqrt{n})}{n^{(q+1)/2}}\ :.
$$ {#eq:jjcavity:operator-Oq}

Here, $J_q(z)$ denotes the Bessel function of order $q$, $n = a^\dagger a$ is
the photon number operator, and $: \cdots :$ implies normal ordering of the
operators. The next step in the derivation consists in performing a second-order
rotating-wave approximation to Hamiltonian
(-@eq:jjcavity:hamiltonian-rot-frame), in order to obtain an effective time-independent
operator [@James2007;@Gamel2010]. By averaging the operator over short time
scales of order
$\omega_J^{-1}$ we obtain

$$
H_\text{eff} = \delta a^\dagger a + \frac{\tilde{E}_J^2}{4 \omega_J}
\sum_{q=0}^\infty \frac{[O_q, O_q^\dagger]}{2q+1} = \left(\delta +
\frac{\tilde{E}_J^2\, \mathcal{G}}{4 \omega_J}\right) n - i \frac{\tilde{E}_J^2}{4\omega_J}
[\mathcal{F}a^\dagger e^{2i \phi} - \text{H.c.}].
$$ {#eq:jjcavity:heff}

We have introduced two highly nonlinear functions of $\Delta_0$ and of the
number operator, namely $\mathcal{F}(\Delta_0, n)$ and
$\mathcal{G}(\Delta_0, n)$. In Appendix -@sec:jjcavity:matrix-elements we provide their explicit expressions, and we
describe their matrix elements in the Fock state basis.

The effective Hamiltonian (-@eq:jjcavity:heff) provides an insight on the
relevant physical processes that we want to address. The first nonlinear term
describes a shift in the resonator frequency [@Meister2015] which can be
observed, e.g., by computing the average photon number in the cavity and the
Cooper-pair tunneling rate (see below, @sec:jjcavity:current), as is
reported in @fig:jjcavity:photon-charge(a). Similarly to the ac-Stark effect in
atomic levels [@James2007;@Gamel2010], this
frequency shift is due to the presence of a strong off-resonant drive (at
frequency $\omega_J$). However, the nonlinearity of $\mathcal{G}(\Delta_0, n)$
also up-converts the off-resonant drive into the resonant process of DCPT. Indeed, the factor $a^\dagger e^{2 i \phi}$ in the
second term of Hamiltonian (-@eq:jjcavity:heff) describes a coherent process in
which a photon is created in the cavity, and two Cooper pairs tunnel through the
junction \[see Eq. (-@eq:jjcavity:phase-operator)\]. To visualize this second-order
process, which is realized through intermediate virtual levels,
 a diagram of the energy levels involved is provided in
@fig:jjcavity:photon-charge(b).

When both $\Delta_0$ and the average photon number $\langle n \rangle$ are much
smaller than 1, an expansion of Eq. (-@eq:jjcavity:heff) to lowest order in
$\Delta_0$ yields the following effective Hamiltonian:

$$
H_{\text{eff}}^{(0)} = \delta'\, n + i \frac{\tilde{E}_J^2 \Delta_0^3}{3 \omega_J} [a e^{-2 i\phi} - a^\dagger e^{2i \phi}],
$$ {#eq:jjcavity:heff-zero-order}

with the linearized frequency shift

$$
\delta' = \delta + \frac{8 \tilde{E}_J^2 \Delta_0^4}{15 \omega_J}.
$$ {#eq:jjcavity:delta-prime}

Equations (-@eq:jjcavity:heff) and (-@eq:jjcavity:heff-zero-order), together with
the master equation (-@eq:jjcavity:me), can be used to compare the effective
Hamiltonian description to the numerical results obtained by time-averaging the
expectation values in the
long-time limit using the full time-dependent Hamiltonian $H(t)$ of Eq. (-@eq:jjcavity:hamiltonian).

## Results for the average photon number {#sec:jjcavity:photon-number}

Armed with the effective Hamiltonian description, we first analyze the photonic
properties on the system in the low-photon-number limit described by Eq. (-@eq:jjcavity:heff-zero-order). With the master equation (-@eq:jjcavity:me) we
write the equations of motion

\begin{align}
\frac{d}{dt}\langle a \rangle &= -\left(i \delta' +
\frac{\kappa}{2}\right)\langle a \rangle - X \langle e^{2 i \phi} \rangle, \\
\frac{d}{dt} \langle e^{2 i \phi} \rangle &= -2 \gamma_\phi\langle e^{2 i \phi},
\end{align}

with $X = \frac{\tilde{E}_J^2\Delta_0^3}{3 \omega_J}$. With the help of the
quantum regression theorem [@Carmichael1999] we obtain

\begin{align}
\langle a^\dagger(t) e^{2i\phi}(t + \tau) \rangle &= \langle a^\dagger (t) e^{2
i \phi} (t)\rangle e^{-2 \gamma_\phi \tau}, \\
    \langle a^\dagger(t) a (t + \tau)\rangle &= \langle a^\dagger(t) a (t) \rangle e^{-(\kappa/2 + i \delta')\tau} - X
    \frac{\langle a^\dagger (t) e^{2 i \phi}(t) \rangle}{\frac{\kappa}{2} - 2 \gamma_\phi + i \delta'}
    \left(e^{-2 \gamma_\phi \tau} - e^{-(\kappa/2 + i \delta')\tau} \right). \label{eq:jjcavity:regression}
\end{align}

For $t \rightarrow \infty$, we get the steady-state values

\begin{align}
    \langle n \rangle = \langle a^\dagger a \rangle &=
    \frac{X^2 \left(1 + \frac{4\gamma_\phi}{\kappa}
    \right)}{\frac{\kappa^2}{4}\left(1 +  \frac{4 \gamma_\phi}{\kappa}  \right)^2 +
    (\delta ')^2}, \label{eq:jjcavity:average-cavity-occupation} \\
    \langle a^\dagger e^{2 i \phi} \rangle &= \frac{- X}{\frac{\kappa}{2} + 2
    \gamma_\phi - i \delta'}. \label{eq:jjcavity:steady-aphi}
\end{align}

Plugging Eqs. (-@eq:jjcavity:average-cavity-occupation)-(-@eq:jjcavity:steady-aphi) into Eq. (-@eq:jjcavity:regression)  leads to the first-order coherence function for the field:

$$
\langle a^\dagger a(\tau) \rangle = \langle n \rangle e^{-(i \delta' +
\kappa/2)\tau} + \frac{X^2 \left( e^{-2 \gamma_\phi \tau} - e^{-(i \delta' + \kappa/2)\tau}\right)}{\left[i \delta' + \left(\frac{\kappa}{2} - 2
\gamma_\phi\right) \right] \left[-i \delta' + \left(\frac{\kappa}{2} + 2 \gamma_\phi \right) \right]},
$$

which in the low-dephasing limit ($\gamma_\phi \ll \kappa$) simplifies to

$$
\langle a^\dagger a(\tau) \rangle \approx \langle n \rangle e^{-2 \gamma_\phi \tau}.
$$ {#eq:jjcavity:first-order-coherence}

![Upper panel: time-averaged cavity occupation calculated with the full
Hamiltonian (-@eq:jjcavity:hamiltonian) (solid blue line), and with
$H_\text{eff}$ (dashed red line). The semiclassical analysis at low $E_J$ gives
the quadratic and quartic contributions (black dotted and green dash-dotted
lines, respectively). The crossover from single- to double-Cooper pair tunneling
occurs for $\tilde{E}_J \Delta_0^2 Q/\omega_0 = \sqrt{5/8}$ (vertical dotted
line). Lower panel: ratio between Cooper-pair tunneling rate and average photon
number, calculated using Eq. (-@eq:jjcavity:current-operator). The
effective-Hamiltonian prediction taken with Eq.
(-@eq:jjcavity:current-operator-eff)
gives exactly 2 (red dashed line). Parameters in both panels are: $\omega_J =
\omega_0/2,\ Q=1500,\ \Delta_0=0.15,\
\gamma_\phi=0$.](./figures/jjcavity-comparisons.pdf){#fig:jjcavity:comparisons}

Equation (-@eq:jjcavity:first-order-coherence) shows that the first-order
coherence function decays at a rate $2\gamma_\phi$, meaning that the resonator
linewidth is four times larger
than for the single-Cooper-pair resonance $\omega_J
\approx \omega_0$ [@Gramich2013]. Hence, signatures of DCPT can be easily found by looking at the resonator spectrum. The average
cavity occupation in Eq. (-@eq:jjcavity:average-cavity-occupation) shows a
lowest-order 
$E_J^4$-dependence and is instead
only weakly dependent on $\gamma_\phi$. In the following, we set $\gamma_\phi
= 0$ as voltage fluctuations can be negligible in experiments
[@Hofheinz2011;@Gramich2013;@Wang2017].

In @fig:jjcavity:comparisons (upper panel) we solve numerically [@Johansson2012;@Johansson2013] the Lindblad
equation Eq. (-@eq:jjcavity:me) with the full Hamiltonian (-@eq:jjcavity:hamiltonian)
and compare the average cavity occupation with the effective Hamiltonian
(-@eq:jjcavity:heff). It is apparent that the $E_J^4$-scaling is valid only at
intermediate values: At higher values of $E_J$, the lowest-order Hamiltonian
$H_\text{eff}^{(0)}$ is insufficient, and one must take into account the
higher-order contributions described by Eq. (-@eq:jjcavity:heff), while the nonlinear
components from the $\mathcal{F}(\Delta_0, n)$ and $\mathcal{G}(\Delta_0, n)$
functions lead to a saturation of the photon number. At the same time,
$H_\text{eff}$ is also not valid at very low $E_J$: It describes
adequately the resonant process at $2\omega_J \approx \omega_0$, but the
dominant process at low $E_J$ is the single-Cooper-pair tunneling oscillating at
the Josephson frequency $\omega_J$.

### Semiclassical interpretation

For $\gamma_\phi = 0$, the system can be mapped onto a nonlinearly driven
harmonic oscillator [@Armour2013;@Meister2015;@Armour2017]. At low-$E_J$, we can
write semiclassical equations of motion for the resonator to better understand
the competition between the $E_J^2$ and $E_J^4$ behaviors. In the limit
$\gamma_\phi \rightarrow 0$, setting $\phi = 0$ without loss of generality, and
assuming that the resonator is in a coherent state $|\alpha\rangle$, we arrive
at the equation of motion

$$
\frac{d \alpha}{d t} = - \left(i \omega_0 + \frac{\kappa}{2} \right)\alpha - i
\tilde{E}_J \Delta_0 \sin [\omega_J t + \Delta_0 (\alpha + \alpha^*)].
$$ {#eq:jjcavity:semiclassical-alpha}

Close to the resonance $2\omega_J \approx \omega_0$ and for very small $E_J$,
the resonator is off-resonantly driven at two frequencies $\pm \omega_J$, and in
the long-time limit $\alpha \approx \alpha_- e^{-i \omega_J t} + \alpha_+
e^{i\omega_J t}$, with the constants $\alpha_\pm$. When $E_J$ gets larger, however, nonlinearities come into play
and the oscillations will be up-converted into an effective drive at frequency
$2\omega_J$. We thus make the ansatz of a solution of the kind

$$
\alpha = \alpha_0 e^{-2 i \omega_J t} + \alpha_- e^{-i \omega_J t} + \alpha_+
e^{i\omega_J t},
$$
leading to:

\begin{align}
\alpha_0 &= -i \frac{\tilde{E}_J \Delta_0^2}{2} \frac{\alpha_- +
\alpha_+^*}{i(\omega_0 - 2 \omega_J) + \frac{\kappa}{2}},\\
    \alpha_- &= \frac{\tilde{E}_J \Delta_0}{2} \frac{1 - i \Delta_0
    \alpha_0}{i(\omega_0 - \omega_J) + \frac{\kappa}{2}},\\
    \alpha_+ &= -\frac{\tilde{E}_J \Delta_0}{2} \frac{1 + i \Delta_0
    \alpha_0^*}{i(\omega_0 + \omega_J) + \frac{\kappa}{2}}.
\end{align}

Setting $\omega_J \approx \omega_0/2$, and since $\kappa \ll \omega_J$, we
can write an expression for the average photon number to fourth order in
$E_J$:

$$
\overline{\langle n \rangle} = |\alpha_0|^2 + |\alpha_+|^2 + |\alpha_-|^2 \approx
\frac{\tilde{E}_J^2 \Delta_0^2} {2} \left[ \frac{\omega_0^2 +
\omega_J^2}{(\omega_0^2 - \omega_J^2)^2} \right] \left\{ 1 + \frac{\tilde{E}_J^2
\Delta_0^4}{2}
\frac{\omega_0^2}{[\omega_0^2 + \omega_J^2]\left[(\omega_0 - 2\omega_J)^2 + \frac{\kappa^2}{4} \right]}\right \}.
$$ {#eq:jjcavity:average-cavity-occupation-semiclassical}

The bar over $n$ implies a time average. It is clear that the $\pm\omega_J$
oscillations give a $\tilde{E}_J^2 \Delta_0^2$ contribution, while the
$2\omega_J$ oscillations give a $\tilde{E}_J^2 \Delta_0^6$ contribution. Both
are plotted as $\mathcal{O}(E_J^2)$ and $\mathcal{O}(E_J^4)$ in the upper panel
of @fig:jjcavity:comparisons. The two contributions have the same weight for
$\tilde{E}_J \Delta_0^2 Q = \sqrt{5/8}\omega_0$.


## Analysis of the charge transport {#sec:jjcavity:current}

### Average current and Cooper-pair tunneling rate

We now turn to a quantitative description of the Cooper-pair transport across
the junction. To begin with, we notice that the dissipative terms in the master
equation (-@eq:jjcavity:me) do not contribute to charge transfer, and is
therefore possible to define a Cooper-pair current operator simply as [@Gramich2013;@Armour2017]

$$
I_\text{CP}(t) = 2e \frac{d N}{d t} = 2 i e [H(t), N] = 2e E_J \sin [\omega_J t - \phi + \Delta_0 (a + a^\dagger)].
$$ {#eq:jjcavity:current-operator}

Averaging over a time $T \gg \omega_J^{-1}$ yields an average dc-current operator

$$
\overline{I}_\text{CP} = \frac{1}{T} \int_{t_0}^{t_0 + T} dt I_\text{CP}(t),
$$ {#eq:jjcavity:current-operator-average}

which leads to the definition of an effective, average Cooper-pair tunneling
rate

$$
\Gamma_\text{CP} = \frac{\overline{I}_\text{CP}}{2e}.
$$

The effective Hamiltonian (-@eq:jjcavity:heff) offers a direct, alternative
definition of the average current as

$$
\overline{I}_\text{CP} = 2ie [H_\text{eff}, N].
$$ {#eq:jjcavity:current-operator-eff}

Using the master equation (-@eq:jjcavity:me), Eq. (-@eq:jjcavity:current-operator-eff) predicts, in the steady state, the
remarkably simple relation 

$$
\frac{\Gamma_\text{CP}}{\kappa \langle n \rangle} = 2.
$$ {#eq:jjcavity:ratio-two}

Equation (-@eq:jjcavity:ratio-two) has a likewise simple interpretation: The
effective Hamiltonian describes a resonator oscillating at the frequency $2
\omega_J$, and photons are generated and destroy only when two Cooper pairs
tunnel. Accordingly, the lower panel of @fig:jjcavity:comparisons shows how Eq.
(-@eq:jjcavity:ratio-two) is in agreement with the full-Hamiltonian calculation
at sufficiently large $E_J$. For low $E_J$, oscillations at the Josephson
frequency cannot be neglected, and $\Gamma_\text{CP}$ drops below 2 as the 
single-Cooper-pair tunneling contribution grows.


### From incoherent to coherent double-Cooper-pair tunneling

![Current Fano factor as a function of $E_J$. The solid red line is obtained
with $H_\text{eff}$, while the orange dots are calculated with the full
Hamiltonian (-@eq:jjcavity:hamiltonian) and the current noise
Eq. (-@eq:jjcavity:current-noise). The dashed gray line marks the semiclassical
crossover at $\tilde{E}_J \Delta_0^2 Q = \sqrt{5/8}\omega_0$. Inset: Fano factor
of the resonator number, using the effective Hamiltonian $H_\text{eff}$.
Parameters are the same as in @fig:jjcavity:comparisons.](./figures/jjcavity-fano.pdf){#fig:jjcavity:fano}

Next, we characterize the transport statistics by defining the time-averaged current noise
[@Blanter2000]:

$$
S_\text{CP} = 2 \mathrm{Re}\ \int_0^{\infty} d\tau \int_{t_0}^{t_0 + T}
\frac{dt}{T} \left[\langle I_\text{CP}(t + \tau) I_\text{CP}(t) \rangle -
\langle I_\text{CP}(t + \tau)\rangle \langle I_\text{CP}(t) \rangle \right].
$$ {#eq:jjcavity:current-noise}

We also define the current Fano factor

$$
F_\text{CP} = \frac{S_\text{CP}}{2e \langle \overline{I}_\text{CP} \rangle},
$$ {#eq:jjcavity:current-fano-factor}

which relates the fluctuations in the current to the noise of a Poissonian
process of single-Cooper-pair tunneling [@Armour2017]. In @fig:jjcavity:fano, we
plot the Fano factor using both definitions of
Eqs. (-@eq:jjcavity:current-operator-average) and (-@eq:jjcavity:current-operator-eff)
in Eq. (-@eq:jjcavity:current-fano-factor). For small $E_J$, using the effective
Hamiltonian, $F_\text{CP} \approx 2$, meaning that *incoherent* DCPT takes place (i.e., an effective charge $4e$ is
transferred) [@Clerk2003]. At larger $E_J$, the value of $F_\text{CP}$ drops,
and since $\Gamma_\text{CP}/\kappa\langle n \rangle \approx 2$ in this regime
(see [@fig:jjcavity:comparisons]), this implies that DCPT is *coherent* [@Grabert2002]. In the inset of @fig:jjcavity:fano we
also consider the cavity Fano factor, $F_n = (\langle n^2 \rangle - \langle n
\rangle^2)/\langle n \rangle$, which drops below 1 witnessing a sub-Poissonian
photon emission statistics. This mechanism is in close analogy to the crossover
between incoherent and coherent *single*-Cooper-pair tunneling occurring when
$\omega_J \approx \omega_0$ [@Gramich2013;@Armour2017]. The orange dots in @fig:jjcavity:fano 
calculated using Eq. (-@eq:jjcavity:current-noise) show that $F_\text{CP}$ drops
below 2 at small $E_J$, again because of single-Cooper-pair tunneling
contributions.

By considering the validity regime of Eq. (-@eq:jjcavity:heff-zero-order) we can
pinpoint the parameter region dominated by incoherent DCPT, where $F_\text{CP} \approx 2$: the lower limit is determined by the
crossover from the $E_J^2$ and the $E_J^4$ behaviour, at $\tilde{E}_J \Delta_0^2
Q \approx \sqrt{5/8} \omega_0$, while the upper limit is fixed by the onset of
strong nonlinearities, when $\langle n \rangle \approx (4 \Delta_0^2)^{-1}$.
Therefore, incoherent behaviour dominates for 

$$
\sqrt{\frac{5}{8}} \ll \frac{\tilde{E}_J \Delta_0^2 Q}{\omega_0} \ll
\sqrt{\frac{3 Q}{8}},
$$

indicating that a high-$Q$ resonator is required to discern this regime from
other transport mechanism.

## Photon blockade and single-photon emission {#sec:jjcavity:g2}

![Second-order correlation function of the cavity field, calculated using Eqs.
(-@eq:jjcavity:heff) and (-@eq:jjcavity:g2), as a function of the zero-point amplitude $\Delta_0$. The
inset shows the behavior of the matrix element $M_{12} = |\langle 1 |
H_\text{eff} | 2 \rangle|$. The parameters used are: $\tilde{E}_J = 0.1
\omega_0,\ \omega_J = \omega_0/2,\ Q = 500,\ \gamma_\phi = 0$.](./figures/jjcavity-g2.pdf){#fig:jjcavity:g2}

For large resonator impedance, zero-point flux fluctuations become large and
$\Delta_0$ can approach unity [@Rolland2019]. In this regime, the light
emission from the resonator becomes strongly nonclassical, and its properties
can be investigated by looking at the second-order field correlation function
defined as [@Scully1997]

$$
g^{(2)} (0) = \frac{\langle a^\dagger a^\dagger a a \rangle}{\langle a^\dagger a \rangle^2},
$$ {#eq:jjcavity:g2}

where the average is performed on the state of the oscillator in the long-time
limit. When $g^{(2)}(0) \ll 1$, the photon-emission statistics is strongly
sub-Poissonian: photons are preferably emitted singularly and the resonator can
function as a single-photon source. This mechanism is in close relation to the
photon blockade of cQED experiments, which arises from the anharmonic spectrum of
the coupled light-matter Hamiltonian [@Imamoglu1997;@Birnbaum2005;@Chang2014;@Biella2015;@Vaneph2018;@Rolland2019].

In @fig:jjcavity:g2, we report
the behavior of $g^{(2)}(0)$ as a function of $\Delta_0$. The value of the
correlation function is close to zero for $\Delta_0 \approx 1.07$, indicating a
strong photon blockade leading to single-photon emission. It is possible to
interpret this result in terms of the effective Hamiltonian of Eq.
(-@eq:jjcavity:heff). Specifically, as shown in Appendix -@sec:jjcavity:matrix-elements,
the matrix element $\langle 1 |
H_\text{eff} | 2 \rangle$ can be calculated analytically in closed form. 
When it vanishes, the resonator remains trapped in the two-state basis spanned
by the Fock states $|0\rangle$ and $|1\rangle$ and further photon emission is
blocked. Consequently, the photon-photon
correlation function drops to zero. In the inset of @fig:jjcavity:g2,
we show that the matrix element $M_{12}$ has indeed a zero for $\Delta_0 \approx
1.07$, corresponding to the minimum of $g^{(2)}(0)$. The presence of vanishing matrix elements of
$H_\text{eff}$ in the Fock basis witnesses a destructive interference of many
processes, each contributing to the photon-blockade effect. It is worth noticing
here that the value $\Delta_0 \approx 1.07$ is significantly lower than the
matrix-element zero corresponding to the single-Cooper-pair tunneling resonance
($\omega_J \approx \omega_0$), occurring for $\Delta_0 = \sqrt{2}$
[@Gramich2013;@Souquet2016] and only slightly higher than the experimentally
achieved value [@Rolland2019].




## Conclusions {#sec:jjcavity:conclusions}

In this Chapter, I have analyzed the emergence of double-Cooper-pair tunneling
in a voltage-biased Josephson junction, when it is coupled to a microwave
resonator: When the Josephson frequency is half the resonator frequency, two
Cooper pairs are required to emit a photon at the cavity frequency. By
developing an effective Hamiltonian theory and performing numerical simulations,
I have shown a crossover from incoherent to coherent double-Cooper-pair
tunneling when increasing the Josephson energy $E_J$. Despite being higher-order
processes in both the Josephson energy and the resonator impedance, this
remarkable transport regime is readily accessible in experiment, as both
quantities can be tuned beyond the perturbative regime (both $E_J \gtrsim
\omega_0$ [@Chen2014] and $\Delta_0 \approx 1$ [@Rolland2019] have been
reached experimentally). 

When the
impedance resonator becomes large, the single-photon nonlinearity in the
spectrum leads to a
photon blockade mechanism which can be exploited to use this system as a
single-photon source. This is readily verified by computing second-order photon
correlation functions for the cavity. Beyond this practical advantage, the
results presented here will stimulate experimental and theoretical
investigation of hitherto unexplored higher-order tunneling processes in
JJ-resonator systems, and overbias photon emission mechanisms in superconducting STM
[@Ast2016].