# Dynamical multistability in a quantum-dot laser {#sec:qdlaser}


*The results presented in this Chapter are published in Ref. \onlinecite{Mantovani2019}.*

## Introduction

One of the most remarkable paradigms of light-matter interaction is the laser [@Haken1984;@Siegman1986;@Scully1997], or maser, when the emitted radiation wavelength falls in the microwave regime.
The lasing mechanism is generally based on a population inversion established in a gaining medium, which interacts with the electromagnetic field inside a cavity.
The gaining medium can be scaled down to a single emitter, realizing a single-atom laser [@Mu1992].
Their full quantum theory has been proposed in the pioneering works of Scully and Lamb [@Scully1966;@Scully1967;@Sargent1993].
Later on, single-atom lasing action has been first observed experimentally in the regime of strong atom-field coupling with a single caesium atom in a high-finesse optical cavity [@McKeever2003].
Unlike conventional lasers, one-atom amplifiers exhibit unique features such as thresholdless behaviour (due to strong coupling), self-quenching, and a sub-Poissonian photon statistics reflecting the quantum nature of both field and atoms [@Carmichael2003].
Within the contest of cavity QED, a related setup is the micromaser [@Meschede1985;@Filipowicz1986;@Lugiato1987;@Walther2006], where a stream of excited atoms is injected into a cavity at a low rate, such that at most one atom at a time resides in the cavity.
The micromaser can display a succession of sharp transitions in the cavity photon number and multistability, whereby two or
more stable amplitudes of oscillation coexist.

The study of single-emitter lasers and micromasers has been later extended to hybrid condensed matter systems, where a number of theoretical models have been proposed and realized experimentally on a variety of platforms, such as superconducting qubits [@Astafiev2007;@Andre2009a;@Ashhab2009], quantum dots embedded in optical photonic crystals [@Nomura2010], devices based on Josephson junctions [@Rodrigues2007;@Rodrigues2007a;@Chen2014;@Cassidy2017], and especially single or double quantum dots [@Jin2011;@Jin2013;@Kulkarni2014;@Liu2014;@Lambert2015;@Gullans2015;@Liu2015;@Liu2017;@Liu2017a;@Agarwalla2019;@Rastelli2019;@Tabatabaei2020].
Recent experiments on self-sustained oscillations in carbon nanotubes with ultra-high quality factors [@Urgell2019;@Wen2019;@Willick2020] are a key step towards the implementation of single-atom *phonon* lasers utilizing mechanical resonators [@Vahala2009;@Mahboob2013].
Furthermore, the engineerization of effective interactions between *spins* and oscillation (mechanical or electromagnetic) has attracted a lot of theoretical and experimental attention over the past few years [@Palyi2012;@Stadler2014;@Stadler2015;@Viennot2015;@Parafilo2016;@Benito2017;@Mi2017;@Mi2018;@Landig2018;@Landig2019;@Cubaynes2019]
(see Ch. -@sec:introduction).

In this Chapter, I present an implementation of a single-atom laser based on the interaction between the electron spin in a quantum dot and a resonant cavity [@Mantovani2019].
The quantum dot is in a spin-valve configuration: It is tunnel coupled to ferromagnetic conductors having noncollinear magnetizations [@Braun2004;@Sahoo2005;@Bordoloi2019].
As a result, the electron transport is highly spin-dependent and influenced by the spin-flip processes due to the coupling with the resonator.
Within this approach, a phonon laser based on the spin-valve mechanism has been proposed [@Khaetskii2013].
Here, I show that our setup represents not only a simple mapping of a widely-studied theoretical model---the one-atom laser---on a solid-state platform.
By contrast, it allows access to a novel and rich regime of multistable laser dynamics.
This phenomenon arises from the breakdown of the rotating-wave approximation (RWA), which is believed to hold when the effective coupling between light and matter $g$ remains small, when compared to the typical energy separation in the uncoupled subsystems (namely, the energy splitting $\Delta\epsilon$ of the two-level system and the cavity frequency $\omega_0$). 

After introducing the system and the model, I provide a semiclassical analysis to explain the multistable behavior, which agrees beautifully with the numerical simulations performed within a full quantum treatment of the resonator.
While the multistable regime of a mechanical oscillator or of a microwave resonator may be in principle detected through optical detection or quantum state tomography of a cavity, we show instead that it maps directly on the electron transport through the quantum dot:
The current displays a telegraph-like noise, associated with current jumps between different state of stability of the resonator. 

## Quantum-dot laser in the spin-valve setup


![Spin-valve-based quantum dot laser. (a) A microwave photon cavity or (b) a nanomechanical resonator mode are coupled to the spin of an electron in a quantum dot, in contact with ferromagnetic leads carrying spin-polarized current. (c) The system can be mapped onto a three-level single-atom laser, where electrons are pumped into the up-spin level, interact coherently emitting photons into the cavity, and decay into the empty-dot state by tunneling into the right lead.](figures/qdlaser-system-2.pdf){#fig:qdlaser:system}

The system I am considering is depicted in @fig:qdlaser:system. It consists of a quantum dot in the single-electron regime, with two nondegenerate spin levels ($|\uparrow\rangle$, $|\downarrow \rangle$), of energy difference $\Delta\epsilon$. 
The dot is embedded between two ferromagnetic contacts, which carry a spin-polarized current of opposite polarization.
The magnetization axis of the leads is chosen parallel to the quantization axis of the electronic spin in the dot.
The spin interacts with a single mode of a local resonator of frequency $\omega_0$, which can be either a microwave cavity or a nanomechanical resonator and is naturally damped through a thermal bath.
The Hamiltonian of the total system is given by

$$
    H_\text{tot} = H + H_\text{leads} + H_\text{tun} + H_\text{bath} + H_\text{osc-bath}.
$$
{#eq:qdlaser:total-hamiltonian}

Here, $H$ is the dot-resonator system Hamiltonian, $H_\text{leads}$ describes
the ferromagnetic contacts, $H_\text{tun}$ provides the electron tunneling
between dot and leads, $H_\text{bath}$ describes a generic
bosonic bath to which the resonator is damped through the coupling term $H_\text{osc-bath}$. The different contributions read:

\begin{align}
    H &= \epsilon_0 (n_\uparrow + n_\downarrow) + \frac{\Delta\epsilon}{2}
    (n_\uparrow - n_\downarrow) + U n_\uparrow n_\downarrow + \omega_0 b^\dagger b +
    g (b + b^\dagger) (d^\dagger_\uparrow d_\downarrow + \text{H.c.}),
    \label{eq:qdlaser:system-hamiltonian} \\
    H_\text{leads} &= \sum_{\nu=L,R} \sum_{k\sigma} (\varepsilon_{\nu k \sigma} -
    \mu_\nu) c^\dagger_{\nu k \sigma} c_{\nu k \sigma},
    \label{eq:qdlaser:leads-hamiltonian} \\
    H_\text{tun} &= \sum_{\nu=L,R} \sum_{k\sigma} V_{\nu \sigma} (c_{\nu k \sigma}
    d^\dagger_{\sigma} + \text{H.c.}),
    \label{eq:qdlaser:tunneling-hamiltonian} \\
    H_\text{osc-bath} &= b B + \text{H.c.}
    \label{eq:qdlaser:osc-bath-hamiltonian}
\end{align}

Let us first describe in detail the system Hamiltonian $H$: 
I have labeled with $\epsilon_0$ the average energy of the spin levels in the
dot, with $\Delta \epsilon$ their energy separation, and with $U>0$ the Coulomb
repulsion strength for the doubly-occupied state. 
The cavity mode (of frequency $\omega_0$) is coupled with strength $g$ to the quantum dot through the spin-oscillator interaction term. $d_\sigma$ is the fermionic operator annihilating an electron of spin $\sigma$ on the dot ($\sigma=\{ \uparrow, \downarrow \} = \{+, - \}$); 
$n_\sigma$ is the corresponding number operator. 
The cavity mode is described by the bosonic annihilation operator $b$.

The ferromagnetic leads are described in $H_\text{leads}$ by Fermi gas Hamiltonians, with $c_{\nu k \sigma}$ the annihilation operator for an excitation of momentum $k$ and energy $\epsilon_{\nu k \sigma}$ on the lead $\nu$ kept at chemical potential $\mu_\nu$. 
$V_{\nu \sigma}$ is the electron tunneling amplitude. 
Finally, $B$ describes a generic operator of the bosonic bath coupled linearly to the resonator.
The three reservoirs are assumed to be at separate thermal equilibria at
temperature $T$ (but can be at different chemical potentials).

## Quantum master equation

In this Section, I derive the Markovian master equation describing the dynamics of the system density matrix $\rho$, under the coherent evolution given by Eq. (-@eq:qdlaser:system-hamiltonian) and the dissipative mechanisms due to the electron tunneling and the resonator damping. 
The external environment is governed by the Hamiltonian $H_E = H_\text{leads} + H_\text{bath}$, and it interacts with the system through $H_\text{int} = H_\text{tun} + H_\text{osc-bath}$. 
Assuming that the weak-coupling and large reservoir limits, the Born-Markov approximation holds (see Ch. -@sec:theory). The Wangsness-Bloch-Redfield master equation for the system (in the Schrödinger picture) reads [@Redfield1957;@Timm2008]:

$$
    \dot{\rho}(t) = -i [ H, \rho(t)] - \int_0^\infty d\tau \text{Tr}_E \{
    [H_\text{int}, [H_\text{int}(-\tau), \rho (t) \rho_\text{leads}(0)
    \rho_\text{bath}(0) ] ] \} \equiv \mathcal{L} \rho(t).
$$
{#eq:qdlaser:bloch-redfield-equation}

Here, $\rho_\text{leads}(0)$ and $\rho_\text{bath}(0)$ represent the states of the reservoirs at the initial time ($t=0$); 
$\mathcal{L}$ is the total Liouvillian
superoperator. 
Its action on $\rho$ can be decomposed into the sum of a coherent part ($\mathcal{L}_H \rho = -i [H, \rho]$) and a dissipative part $\mathcal{L}_\text{leads}\rho + \mathcal{L}_\text{bath}\rho$.

In order to simplify Eq. (-@eq:qdlaser:bloch-redfield-equation), I make three important assumptions:

(i) A large bias voltage $V$ is applied symmetrically to the leads, such that $\mu_L = eV/2$, $\mu_R = -eV/2$, with $V \gg \{k_B T; \omega_0; \Delta\epsilon\}$. 
  The Fermi functions $f_\nu (x) = \{\exp [(x - \mu_\nu)/k_BT] + 1\}^{-1}$ are then $f_L \approx 1$ and $f_R \approx 0$, with no energy dependence. 
  This means that electron transport is allowed from the left to the right lead, while it is completely suppressed in the opposite direction.

(ii) The Coulomb repulsion strength $U$ is the largest energy scale in the system, such that the doubly-occupied state in the quantum dot becomes energetically inaccessible, and it cannot be thermally populated. 
  Thus, it can be projected out from the electronic Hilbert space. Formally, we have $U \gg eV \gg \{k_B T; \omega_0; \Delta\epsilon\}$. 
  Consequently, the coherent dynamics given by Hamiltonian (-@eq:qdlaser:system-hamiltonian) only involves the two spin states ($|\uparrow \rangle, |\downarrow \rangle$). 
  This allows us to map Eq. (-@eq:qdlaser:system-hamiltonian) into the Rabi Hamiltonian [@Rabi1936;@Rabi1937]:
  $$
  H = \frac{\Delta\epsilon}{2}\sigma_z + \omega_0 b^\dagger b + g (\sigma_+ +
      \sigma_-)(b + b^\dagger).
  $$
  {#eq:qdlaser:rabi-hamiltonian} 
  The Coulomb term in Eq. (-@eq:qdlaser:system-hamiltonian) has been removed, and we have mapped the dot operators to the Pauli algebra through $n_\uparrow-n_\downarrow \rightarrow \sigma_z$ and $d_\uparrow^\dagger d_\downarrow \rightarrow \sigma_+$, after projecting out the empty and the doubly-occupied states from the dot Hilbert space. 
  Notice that, however, the empty state $|0\rangle$ must be taken into account in the open dynamics: 
  In the single-electron regime, transport through the dot necessarily involves the third state $|0\rangle$.
  
(iii) The spin-resonator coupling strength $g$ is small compared to the typical bare energy scale of the two subsystems, i.e., 
  $g \ll \{\Delta\epsilon; \omega_0\}$.
  This allows us to safely treat Eq. (-@eq:qdlaser:bloch-redfield-equation) in the local basis $|n, \sigma\rangle$ of the uncoupled system ($n$ is the Fock number of the resonator), deriving a *local* master equation [@Cattaneo2019]. 

Under assumptions (i)-(iii), the time integrals in Eq. (-@eq:qdlaser:bloch-redfield-equation) can be performed and they lead to a quantum master equation in Lindblad form:

$$
    \dot{\rho} = -i[H, \rho] + \sum_{\sigma} [\Gamma_L^\sigma
    \mathcal{D}(F^\dagger_\sigma)\rho + \Gamma^\sigma_R \mathcal{D} (F_\sigma)\rho] + \kappa (1 + n_B)\mathcal{D}(b) +
    \kappa n_B \mathcal{D}(b^\dagger).
$$
{#eq:qdlaser:lindblad-equation}

Within the wide-band approximation, whereby the spectral densities of the
fermionic reservoirs are energy-independent, the electronic tunneling rates are given by
$\Gamma_\nu^\sigma = 2 \pi |V_{\nu \sigma}|^2 \rho_{\nu \sigma}$, where
$\rho_{\nu\sigma}$ is the spin-$\sigma$ density of states at the Fermi level of lead $\nu$. 
Since the leads are ferromagnetic, the spin dependence of the
tunneling rates can be written in terms of the polarizations $P_\nu$ of the
leads as $\Gamma_\nu^\sigma = \Gamma_\nu ( 1 + \sigma P_\nu)/2$, 
with $0 \leq P_\nu\leq 1$. 
Although for standard ferromagnetic contacts (Co and PdNi alloys) $P_\nu$ does not exceed 0.5 [@Sahoo2005;@Viennot2015;@Cubaynes2019], it can be equal to unity for half-metallic leads [@Ziese2002]. In the following, without losing generality in the results, I will assume a symmetric and opposite polarization of the leads, i.e., $P_L = -P_R = P$, allowing us to write the tunneling rates as

$$
\begin{aligned}
\Gamma_{L}^{\uparrow}&=\Gamma_{L}\left(\frac{1+P}{2}\right), \quad \Gamma_{L}^{\downarrow}=\Gamma_{L}\left(\frac{1-P}{2}\right), \\ \Gamma_{R}^{\uparrow}&=\Gamma_{R}\left(\frac{1-P}{2}\right), \quad \Gamma_{R}^{\downarrow}=\Gamma_{R}\left(\frac{1+P}{2}\right).
\end{aligned}
$$
{#eq:qdlaser:tunneling-rates}

In the electronic dissipators, we have replaced the fermionic operator $d_\sigma$ with $F_\sigma = (1 - n_{-\sigma})d_\sigma$, such that the population and the coherences involving the doubly-occupied state are forced to vanish, as required by the large Coulomb repulsion limit.

For the cavity, $\kappa = \omega_0 / Q$ is the decay rate ($Q$ is
the quality factor), and $n_B = [\exp(\omega_0 / k_B T) - 1]^{-1}$ is the average number of excitations in the thermal bath at the resonator frequency. 

The steady state $\rho_\text{st}$ of Eq. (-@eq:qdlaser:lindblad-equation), obtained by solving $\dot{\rho} = 0$, is the central object that I will employ to calculate the relevant stationary expectation values for the system, together with the cavity Fock distribution. 
The numerical steady solution of Eq. (-@eq:qdlaser:lindblad-equation) has been found using [qutip]{.smallcaps} [@Johansson2012].

## Single-atom laser within the rotating-wave approximation (RWA)

### RWA Hamiltonian and Jaynes-Cummings model

![Energy ladder of the uncoupled dot-resonator system (black lines). When $\Delta\epsilon \approx \omega_0$, the RWA Hamiltonian (-@eq:qdlaser:rwa-hamiltonian) describes the hybridization of the pairs $\{|\uparrow, n\rangle, |\downarrow,n+1\rangle\}$, giving rise to the Jaynes-Cummings doublet \[green lines, Eq. (-@eq:qdlaser:rwa-doublet)\]. The energy separation of each doublet grows proportionally to $g\sqrt{n+1}$. ](figures/qdlaser-rwa-diagram.pdf){#fig:qdlaser:rwa-diagram}

The theoretical treatment of the single-atom laser usually relies upon the so-called rotating-wave approximation (RWA). It consists in replacing the Rabi model Hamiltonian with the Jaynes-Cummings (JC) one [@Jaynes1963]:

$$
	H_\text{RWA} = \frac{\Delta \epsilon}{2} \sigma_z + \omega_0 b^\dagger b + g (\sigma_+ b + \sigma_- b^\dagger).
$$
{#eq:qdlaser:rwa-hamiltonian}

Equation (-@eq:qdlaser:rwa-hamiltonian) can be derived from the Rabi Hamiltonian (-@eq:qdlaser:rabi-hamiltonian) by moving to interaction picture and assuming a small detuning $\delta \ll \omega_0 + \Delta\epsilon$, with $\delta = |\omega_0 - \Delta\epsilon|$. Further, one must also ensure that $g\sqrt{n} \ll \{\omega_0; \Delta \epsilon\}$, where $n$ is the cavity Fock number. In @fig:qdlaser:rwa-diagram, I sketch the effective interaction described by the JC model, close to resonance, i.e., $\delta \approx 0$. The interaction (-@eq:qdlaser:rwa-hamiltonian) mixes *only* the states $|\uparrow, n\rangle$ and $|\downarrow,n+1\rangle$. Diagonalizing the interaction in the subspace spanned by these two states yields the JC doublet

$$
\begin{aligned}
    |n+\rangle &= \cos \left(\frac{\theta}{2}\right) |\uparrow, n\rangle + \sin \left( \frac{\theta}{2}\right) |\downarrow,n+1\rangle, \\
    |n-\rangle &= -\sin \left(\frac{\theta}{2}\right) |\uparrow, n\rangle + \cos \left(\frac{\theta}{2}\right)  |\downarrow,n+1\rangle, 
\end{aligned}
$$
{#eq:qdlaser:rwa-doublet}

with eigenenergies $E_{n\pm} = \omega_0 \left(n + \frac{1}{2} \right) \pm \frac{1}{2}\sqrt{4g^2(n+1) + \delta^2}$ and mixing angle $\theta = \arctan \left( \frac{2 g \sqrt{n+1}}{\delta}\right)$. Notice that the doublet splitting is proportional to $g\sqrt{n+1}$, i.e., the eigenvalue ladder of the JC Hamiltonian is anharmonic. If $g\sqrt{n+1}$ becomes comparable with the bare level separation of the uncoupled system, Eq. (-@eq:qdlaser:rwa-hamiltonian) is no longer valid, as doublets with different quantum number $n$ become closer in energy, and their interaction (which is not captured by the RWA Hamiltonian) cannot be neglected anymore. For atom-field cavity systems, $g/\omega_0 \sim 10^{-7}-10^{-6}$ [@Raimond2001], hence Eq. (-@eq:qdlaser:rwa-hamiltonian) remains valid also for photon numbers $n \gg 1$ characteristic of laser emission. However, as the technological development in solid-state devices pushes the ratio $g/\omega_0$ order of magnitudes higher, a more careful consideration on the validity of the RWA is necessary when designing systems with large photon number.


### Density matrix theory {#sec:qdlaser:rwa-density-matrix}

Using Eq. (-@eq:qdlaser:rwa-hamiltonian) in Eq. (-@eq:qdlaser:lindblad-equation), it is possible to obtain the approximate analytical expression for the steady-state Fock distribution of the cavity, $p_n$ (details of the calculations are provided in Appendix -@sec:qdlaser:pn-analytical-derivation). 
In particular, when close to the resonant condition $\Delta\epsilon = \omega_0$, the $p_n$ is obtained through the following recursive equation:

$$
	\left[\frac{n \kappa \left( \frac{g}{g_\text{thr}}\right)^2}{1 + \frac{n}{A_s^2}\left(\frac{g}{g_\text{thr}}\right)^2} + \kappa n_B \right] p_{n-1} = \kappa( 1+ n_B)p_n,
$$
{#eq:qdlaser:rwa-pn}

with the quantities

$$
A_s^2 = \frac{\Gamma_L \Gamma_R P}{(2\Gamma_L + \Gamma_R)\kappa}, \qquad g_\text{thr} = \sqrt{\frac{\Gamma_R \kappa}{4} \left[ \frac{2 \Gamma_L (1 + P^2)+ \Gamma_R(1 - P^2)}{4 \Gamma_L P}\right]}.
$$
{#eq:qdlaser:rwa-pn-saturation-thr}

I will refer to $A_s^2$ as the *saturation number* and to $g_\text{thr}$ as the *lasing threshold* coupling. 
The solution of Eq. (-@eq:qdlaser:rwa-pn) reads

$$
	p_n = p_0 \frac{\mathcal{N}_n}{\mathcal{D}_n} \left ( \frac{n_B}{n_B+1}\right)^n.
$$
{#eq:qdlaser:rwa-pn-solution}

Here, I have introduced the Pochhammer symbol (or rising factorial) $a_n = a (a+1) (a+2)\cdots (a + n - 1)$ and defined the quantities

$$
	\mathcal{N} = 1 + \frac{A_s^2}{n_B} + A_s^2 \frac{g_\text{thr}}{g}, \qquad \mathcal{D} = 1 + A_s^2 \frac{g_\text{thr}}{g}.
$$

The prefactor $p_0$ is the zero-Fock occupation number and is determined from the normalization condition $\sum_{n=0}^\infty p_n = 1$, leaving us with

$$
	p_0 = \frac{1}{_2 F_1 \left(1, \mathcal{N}; \mathcal{D}; \frac{n_B}{n_B+1}\right)},
$$

with the hypergeometric function $_2 F_1 (a, b; c; z)$. At zero temperature, one finds instead $p_0 = \left[_1 F_1 (1; A_s^2 g_\text{thr}^2 / g^2; n_B/(n_B+1)\right]^{-1}$, with the confluent hypergeometric function $_1 F_1(a; b; z)$. 
A more interesting quantity is represented by the average photon number, $\bar{n} = \sum_{n=0}^\infty n p_n$, which is found to be

$$
	\bar{n} = A_s^2 \left[1 -  \left(\frac{g_\text{thr}}{g} \right)^2\right] + n_B + n_B(n_B + 1) \left(\frac{g_\text{thr}}{g} \right)^2 p_0.
$$
{#eq:qdlaser:rwa-pn-nav}

The analytical expression Eq. (-@eq:qdlaser:rwa-pn-solution) allows also to calculate also the second moment of the distribution, $\overline{n^2}$, and from it the Fano factor [@Fano1947;@Scully1997], defined as 

$$
\mathcal{F} = \frac{\overline{n^2} - \bar{n}^2}{\bar{n}}.
$$

The Fano factor yields a measure of the dispersion of a probability distribution, compared to the variance of a Poissonian distribution (where it is equal to the mean). The Poissonian distribution has a Fano factor equal to 1, and is characteristic of a classical laser. A sub-Poissonian distribution ($\mathcal{F}<1$) can be instead signature of quantum behavior in the resonator, as the distribution tends to be more narrow and different from a classical field. The resulting expression for the Fano factor is

$$
\mathcal{F} = 1 + n_B + \frac{1}{\left(\frac{g}{g_\mathrm{thr}}\right)^2 + \frac{n_B + A_s^2}{n_B +1}} - \alpha p_0,
$$

where I have gathered in $\alpha$ an unimportant factor multiplying $p_0$. At zero temperature and far above threshold ($g \gg g_\mathrm{thr}$), the above expression becomes

$$
\mathcal{F} = 1 + \frac{1}{\left(\frac{g}{g_\mathrm{thr}}\right)^2 + A_s^2 -1}.
$$

The Fano factor is peaked around the threshold coupling, where quantum fluctuations increase, witnessing a dramatic change in the Fock distribution, analogously to a phase transition. Above threshold, the analytical expression predicts that $\mathcal{F}$ tends to the Poissonian value, $\mathcal{F} = 1$.

In @fig:qdlaser:rwa(a)-(b), I have collected the analytical results for both $\bar{n}$ and $\mathcal{F}$, together with their comparison to a numerical calculation using the full master equation and the RWA Hamiltonian and to the result from a semiclassical analysis (see below).


### Semiclassical theory {#sec:qdlaser:rwa-semiclassical}

Equation (-@eq:qdlaser:lindblad-equation) allows us to derive the equation of motion for the expectation value of a system operator $O$, defined as $\langle O \rangle = \text{Tr} (O \rho)$:

$$
\begin{aligned}
        \langle \dot{O} \rangle = &-i \langle [O, H] \rangle + \sum_\sigma \left[
         \Gamma_L^\sigma \left(\langle F_\sigma O F^\dagger_\sigma \rangle - \frac{1}{2}
        \langle \{ O, F_\sigma F_\sigma^\dagger \} \rangle \right) + \Gamma_R^\sigma \left( \langle F^\dagger_\sigma O F_\sigma \rangle - \frac{1}{2}
        \langle \{ O, F^\dagger_\sigma F_\sigma \} \rangle\right) \right] + \\
        &+ \kappa (1 + n_B) \left(\langle b^\dagger O b \rangle - \frac{1}{2}
        \langle \{ O, b^\dagger b \} \rangle \right) + \kappa n_B \left( \langle
        b O b^\dagger \rangle - \frac{1}{2}
        \langle \{ O, b b^\dagger \} \rangle\right). 
\end{aligned}
$$
{#eq:qdlaser:eq-motion}


Using the RWA Hamiltonian in Eq. (-@eq:qdlaser:eq-motion), one is able to derive a nonlinear set of dynamical equations which govern the system. 
For simplicity, let us set $P=1$ and $n_B = 0$. The equations read:

\begin{align}
	\langle \dot{n}_\uparrow \rangle &= - \Gamma_L \langle n_\uparrow - n_\downarrow \rangle - ig\langle b \sigma_+ - b^\dagger \sigma_- \rangle + \Gamma_L, \label{eq:qdlaser:rwa-semiclassical-system-1} \\
	\langle \dot{n}_\downarrow \rangle &= \Gamma_R \langle n_\downarrow \rangle + ig \langle b \sigma_+ - b^\dagger \sigma_-\rangle, \\
	\langle \dot{\sigma}_- \rangle &= \left(- i \Delta \epsilon - \frac{\Gamma_R}{2} \right) \langle \sigma_- \rangle + ig \langle (b + b^\dagger) (n_\uparrow - n_\downarrow) \rangle, \quad \text{and c.c.}, \\
	\langle \dot{b} \rangle &= \left( -i \omega_0 - \frac{\kappa}{2}\right) \langle b \rangle - ig\langle \sigma_- \rangle, \quad \text{and c.c.} \label{eq:qdlaser:rwa-semiclassical-system-4}
\end{align}

Equations (-@eq:qdlaser:rwa-semiclassical-system-1)-(-@eq:qdlaser:rwa-semiclassical-system-4) contain expectation values involving both spin and cavity operators: Writing equations of motions for these terms will involve higher order operators and lead to an infinite hierarchy of equations, which cannot be solved. 
To close the set of equations, the semiclassical approximation is needed. It consists in replacing the bosonic operator $b$ by the complex number $\alpha = Ae^{i\phi}$. $A$ and $\phi$ correspond to the semiclassical amplitude and phase of the resonator, respectively. 
The semiclassical approximation basically amounts to neglecting the quantum fluctuations in the resonator. 
Let us now switch to a rotating frame with the replacements $\langle \sigma_- \rangle \rightarrow \langle \sigma_- \rangle e^{-i \Delta\epsilon t}$ and $\alpha \rightarrow \alpha e^{-i\omega_0 t}$. Further, I set the spin quantities $S_x = \langle \sigma_+ + \sigma_-\rangle$, $S_y = -i\langle \sigma_+ - \sigma_- \rangle$, $S_z = \langle n_\uparrow - n_\downarrow \rangle$, and the total dot occupation $p_1 = \langle n_\uparrow + n_\downarrow \rangle$. Finally, we consider the special case where $\Gamma_L = \Gamma_R/2 = \Gamma$ in Eq. (-@eq:qdlaser:tunneling-rates): With this condition, the equation for $p_1$ decouples from the rest of the system and can thus be disregarded. However, it is worth remarking that this condition does not alter the physical results, but merely yields a simplification in the analytical calculations. At resonance $(\Delta \epsilon = \omega_0$) one obtains the set of equations

\begin{align}
\dot{S}_{x} &=-\Gamma S_{x}-2 g A \sin \phi S_{z}, \label{eq:qdlaser:rwa-semiclassical-system-alt-1} \\
\dot{S}_{y} &=-\Gamma S_{y}-2 g A \cos \phi S_{z}, \\ 
\dot{S}_{z} &= \Gamma-\Gamma S_{z}+2 g A\left(\sin \phi S_{x}+\cos \phi S_{y}\right), \\ 
\dot{A} &=-\frac{\kappa}{2} A+\frac{g}{2}\left(-\sin \phi S_{x}+\cos \phi S_{y}\right), \\ 
\dot{\phi}&=-\frac{g}{2 A}\left(\cos \phi S_{x}-\sin \phi S_{y}\right). \label{eq:qdlaser:rwa-semiclassical-system-alt-5}
\end{align}

The equations for $A$ and $\phi$ have been obtained by transforming the corresponding equations for $\alpha$ and $\alpha^*$, using $\alpha = A e^{i\phi}$. The system (-@eq:qdlaser:rwa-semiclassical-system-alt-1)-(-@eq:qdlaser:rwa-semiclassical-system-alt-1) is a nonlinear system describing lasing oscillations in the resonator, under the influence of the spin dynamics  (driven by the single-electron tunneling at rate $\Gamma$). It easy to study the stability of the resonator by setting the time derivatives to zero and looking for the stationary solution. The solution is independent on $\phi$, which can be set to zero. For finite polarization $P<1$, the dynamical equation for the resonator amplitude is

$$
\dot{A}=-\frac{A}{2}\left[\kappa-\frac{\frac{2 g^{2} P}{\Gamma}}{1+\left(\frac{2 g A}{\Gamma}\right)^{2}}\right]=-\frac{A}{2}\left[\kappa+\gamma_{\mathrm{RWA}}(A)\right].
$$
{#eq:qdlaser:rwa-resonator-amplitude}
In the latter equality, we have defined an effective, negative nonlinear damping $\gamma_\mathrm{RWA}(A)$: It encompasses the driving effect of the electron dynamics, opposing the natural damping of the resonator at decay rate $\kappa$. By setting $\dot{A} = 0$, we obtain a nonlinear algebraic equations with two steady solutions for the occupation number $\bar{n} = A^2$, i.e.:

$$
\bar{n} = 0 \quad \text{and} \quad
\bar{n}=A_{s}^{2}\left[1-\left(\frac{g_{\mathrm{thr}}}{g}\right)^{2}\right].
$$
{#eq:qdlaser:rwa-nav-semiclassical}

The saturation number and threshold coupling read

$$
A_s^2 = \frac{\Gamma P}{2 \kappa}, \quad g_\text{thr} = \sqrt{\frac{\Gamma \kappa}{2 P}},
$$

respectively. Notice that they are in full agreement with Eq. (-@eq:qdlaser:rwa-pn-saturation-thr), obtained with the full quantum approach. The stability analysis of Eq. (-@eq:qdlaser:rwa-resonator-amplitude) reveals a bifurcation point at $g= g_\mathrm{thr}$: for $g \leq g_\mathrm{thr}$, the only (stable) solution is $\bar{n} = 0$, corresponding to absence of lasing. For $g > g_\mathrm{thr}$, $\bar{n} = 0$ is an unstable solution, while the other solution with large $\bar{n}$ is the stable one. For a large coupling, $g \gg g_\mathrm{thr}$, the nonlinear damping $\gamma_\mathrm{RWA}$ essentially loses its dependence on $g$. This causes the average occupation number to saturate to the value $A_s^2$. Notice that the $A$-dependence of $\gamma_\mathrm{RWA}$ is crucial to obtain a finite steady value for $\bar{n}$. For completeness, I report here the more general expression of $\gamma_\mathrm{RWA}$ for arbitrary tunneling rates and $P=1$ (obtained by including the equation for $p_1$ in the stationary solution of the system), which is given by

$$
\gamma_{\mathrm{RWA}}(A)=-\frac{g^{2} \Gamma_{\mathrm{eff}}}{g^{2} A^{2}+\Gamma_{\mathrm{eff}} \Gamma_{R} / 4},
$$

where $\Gamma_{\mathrm{eff}}=\Gamma_{L} \Gamma_{R}/(2 \Gamma_{L}+\Gamma_{R})$. Accordingly, the expressions for saturation number and threshold coupling become

$$
g_{\mathrm{thr}}^{2}=\frac{\Gamma_{R}\kappa}{4}, \quad A_{s}=\sqrt{\frac{\Gamma_{\mathrm{eff}}}{\kappa} }
$$

The lasing threshold for fully polarized leads only depends on the right-tunneling rate $\Gamma_R$, as a consequence of the large Coulomb repulsion in the dot: Since only one electron can reside in the system at a given instant, a photon-emission event can only take place if the electron has tunneled out into the right contact. A large value of $\Gamma_R$ pushes the threshold for lasing to appear to higher values of $g$.

In @fig:qdlaser:rwa I summarize the analytical results obtained with the semiclassical approximation and with the density matrix approach of @sec:qdlaser:rwa-density-matrix, together with the numerical results. The analytical expression for $\bar{n}$ agrees well with the numerics. However, the Fano factor is slightly overestimated by the density matrix approach. Indeed, the numerics predicts a sub-Poissonian Fock distribution for the resonator, with $\mathcal{F}$ well below one. This can be observed more closely by inspecting the $p_n$ distributions in @fig:qdlaser:rwa(c)-(e) for the below threshold, slightly above threshold, and far above threshold couplings $g$. The sub-Poissonian statistics is a typical signature of the single-atom laser [@Carmichael2003;@McKeever2003;@Astafiev2007].

![Single-atom laser within the RWA. (a) Average phonon occupation and (b) Fano factor for the Fock distribution, on resonance ($\Delta\epsilon = \omega_0$), as a function of the spin-resonator coupling $g$. The solid red line is computed from the analytical expression for $p_n$, Eq. (-@eq:qdlaser:rwa-pn-nav), the filled blue circles represent the numerical calculation, and the green dashed line is the semiclassical result [Eq. (-@eq:qdlaser:rwa-nav-semiclassical)]. The vertical dotted grey lines correspond to the threshold coupling, $g_\mathrm{thr} = 0.005\omega_0$. (c)-(e) Probability distributions for the oscillator Fock number at three different values of $g$, calculated numerically (bars) and analytically (solid red curves). Parameters: $\Gamma = 0.1\omega_0,\ Q = \num{e3},\ T = 0,\ P = 1$.](figures/qdlaser-rwa.pdf){#fig:qdlaser:rwa}

## Beyond the rotating-wave approximation: Multistability

In the previous Section, I have analyzed the lasing behavior of the quantum-dot laser assuming the validity of the RWA Hamiltonian (-@eq:qdlaser:rwa-hamiltonian). However, the _a posteriori_ analysis of the results reveals a peculiarity of our system with respect to similar (atomic and solid-state) implementations of the single-atom laser: The spin-dependent transport yields a very large photon emission efficiency (for the case of full polarization and negligible spin relaxation processes, each electron passing through the dot emits one photon), such that the quantity $g\sqrt{\bar{n}}$ becomes comparable to the bare resonator frequency $\omega_0$, even for relatively low values of the ratio $g/\omega_0$. Consequently, without the technologically challenging need to enter the *ultrastrong-coupling* regime (characterized by $g/\omega_0 \sim 1$, see Ch. -@sec:introduction), the system is expected to show a clear and unique deviation from the RWA physics. 

**Write here a sentence to explain what we are going to do.**

### Semiclassical equations beyond RWA

The semiclassical approach utilized in @sec:qdlaser:rwa-semiclassical can be extended beyond the RWA by using the full Rabi Hamiltonian (-@eq:qdlaser:rabi-hamiltonian). First, Eqs. (-@eq:qdlaser:rwa-semiclassical-system-1)-(-@eq:qdlaser:rwa-semiclassical-system-4) are replaced by

\begin{align}
\left\langle\dot{ n}_{\uparrow}\right\rangle &=-\Gamma_{L}
\langle n_{\uparrow}\rangle-\Gamma_{L}\langle n_{\downarrow}\rangle-i g\langle( b+ b^{\dagger})(\sigma_{+}-\sigma_{-})\rangle+\Gamma_{L}, \\ \langle\dot{ n}_{\downarrow}\rangle &=-\Gamma_{R}\langle n_{\downarrow}\rangle+i g\langle( b+ b^{\dagger})(\sigma_{+}-\sigma_{-})\rangle, \\ \langle\dot{\sigma}_{-}\rangle &=\left(-i \Delta \epsilon-\frac{\Gamma_{R}}{2}\right)\langle\sigma_{-}\rangle+i g\langle( b+ b^{\dagger}) \sigma_{z}\rangle, \quad \text{and c.c.}, \\ \langle\dot{ b}\rangle &= \left(-i \omega_{0}-\frac{\kappa}{2}\right)\langle b\rangle-i g\langle \sigma_{+}+\sigma_{-}\rangle, \quad \text{and c.c.}
\end{align}

After the semiclassical approximation, setting again $\Gamma_L = \Gamma_R/2  = \Gamma$, $P=1$ and $\Delta\epsilon = \omega_0$, and moving to the rotating frame as above, the system can be written as

\begin{align} 
\dot{S}_{x} &=- \Gamma S_{x}-2 g A[\sin (2 \omega_{0} t-\phi)+\sin \phi] S_{z}, \label{eq:qdlaser:no-rwa-semiclassical-system-1} \\ 
\dot{S}_{y} &=-\Gamma S_{y}-2 g A[\cos (2 \omega_{0} t-\phi)+\cos \phi] S_{z},\label{eq:qdlaser:no-rwa-semiclassical-system-2} \\ 
\dot{S}_{z} &=-\Gamma S_{z}+2 g A \{[\sin (2 \omega_{0} t-\phi)+\sin \phi] S_{x} +[\cos (2 \omega_{0} t-\phi)+\cos \phi] S_{y}\}+\Gamma, \label{eq:qdlaser:no-rwa-semiclassical-system-3} \\ 
\dot{A} &= -\frac{\kappa}{2} A+\frac{g}{2}\left\{\left[\sin \left(2 \omega_{0} t-\phi\right)-\sin \phi\right] S_{x} + \left[\cos \left(2 \omega_{0} t-\phi\right)+\cos \phi\right] S_{y}\right\}, \label{eq:qdlaser:no-rwa-semiclassical-system-4} \\ 
\dot{\phi} &= -\frac{g}{2 A}\left\{\left[\cos \left(2 \omega_{0} t-\phi\right)+\cos \phi\right] S_{x} -\left[\sin \left(2 \omega_{0} t-\phi\right)+\sin \phi\right] S_{y}\right\}. \label{eq:qdlaser:no-rwa-semiclassical-system-5}
\end{align}

Contrary to the RWA case, this nonlinear set of equations does not have a stationary solution in the rotating frame: The resonator will approach a limit-cycle oscillating solution. A more thorough analysis is then necessary to gain further quantitative information on the long-time limit behavior of the resonator dynamics.
The basic idea consists in studying Eqs. (-@eq:qdlaser:no-rwa-semiclassical-system-1)-(-@eq:qdlaser:no-rwa-semiclassical-system-5) in Fourier space, deriving an effective expression for the nonlinear damping of the resonator, $\gamma_\mathrm{eff}$. To achieve this, it is crucial to exploit the fact that the dynamics of the resonator amplitude $A$ is *slow* compared to the electronic driving and to the oscillations themselves, a condition satisfied when $\kappa \ll \{\Gamma; g; \omega_0\}$. I also assume that, eventually, the phase of the oscillator evolves harmonically as $\phi = \omega_0 t$. This can be checked empirically by solving numerically the time-dependent equations (-@eq:qdlaser:no-rwa-semiclassical-system-1)-(-@eq:qdlaser:no-rwa-semiclassical-system-5) and looking at the long-time behavior of $A(t)$. It is then possible to average the charge dynamics \[determined by Eqs. (-@eq:qdlaser:no-rwa-semiclassical-system-1)-(-@eq:qdlaser:no-rwa-semiclassical-system-3)\] over a resonator period $\mathcal{T} = 2\pi/\omega_0$, during which $A$ can be considered constant. In this way, one can calculate the coarse-grained effect of the charge dynamics on the resonator amplitude. Averaging Eq. (-@eq:qdlaser:no-rwa-semiclassical-system-4) over $\mathcal{T}$, we obtain

$$
\dot{\bar{A}} = -\frac{\kappa}{2}\bar{A} + \frac{g}{\mathcal{T}} \int_0^\mathcal{T} dt' \cos (\omega_0 t') S_y(t', \bar{A}),
$$
{#eq:qdlaser:average-A-equation-1}

where $\bar{A} = 1/\mathcal{T}\int_0^\mathcal{T} dt' A(t')$. The cosine term in the integration acts as a Fourier filter, eliminating all Fourier components of $S_y$ except for the ones oscillating at $\pm\omega_0 t$. The Fourier transforms of the spin quantities are defined as

$$
\begin{aligned}
  S_k (t) &= \sum_{n=-\infty}^{+\infty} e^{-i\omega_0 n t}S_k^{(n)}, \\
  S_k^{(n)}&= \frac{1}{\mathcal{T}}\int_0^\mathcal{T} dt S_k(t) e^{i\omega_0 n t},
\end{aligned}
$$

with $k = x,y,z$. Because $S_k (t)$ is real, then $S_k^{(-n)} = S_k^{(n)*}$, and Eq. (-@eq:qdlaser:average-A-equation-1) becomes

$$
\dot{\bar{A}} = -\frac{\kappa}{2}\bar{A} + g\ \mathrm{Re} [S_y^{(1)}(\bar{A})] = -\frac{\bar{A}}{2} [ \kappa + \gamma_\mathrm{eff}(\bar{A})].
$$
The nonlinear damping is given by 

$$
\gamma_\mathrm{eff}(\bar{A}) = - \frac{2 g}{\bar{A}}\ \mathrm{Re} [S_y^{(1)}(\bar{A})].
$$

### Multistability analysis of the resonator

From the equation of the nonlinear damping to the multistability; plot of nonlinear damping and stability diagram, comparison with numerical results.

## Detecting lasing and multistability with transport measurements

In this Section, I describe a method to efficiently detect the lasing state in the resonator as well as the multistability.

### Full-counting statistics

### Current jumps and two-state model

## Experimental feasibility: Multistability in nonideal cases

The scope of this Section is to provide a realizability study for the system, by considering a number of processes which can hinder the onset of lasing and multistability in the cavity. Specifically, I will consider: (a) The effect of finite polarization ($P<1$) in the leads; (b) the effect of spin relaxation in the quantum dot; (c) for low-frequency mechanical oscillators, the effect of finite temperature and Duffing nonlinearity.

### Effect of finite temperature and finite polarization

![Effect of finite temperature and finite polarization. (a) Average occupation number in the oscillator at resonance as a function of the dot’s energy splitting $\Delta\epsilon$ and of the spin-resonator coupling strength $g$. (b) Stability diagram of the resonator. Italic numbers indicate the number of distinct peaks in the Fock distribution. Parameters: $Q = \num{e3}$, $\Gamma = 0.1\omega_0$, $P = 0.5$, and $T=10\omega_0$.](figures/qdlaser-finite-t-p.pdf){#fig:qdlaser:finite-t-p}

At finite polarization, a fraction of the total current passing through the quantum dot is elastic, i.e., no energy is exchanged with the resonator. This inevitably leads to a lower pumping efficiency, and it is natural to understand how sensitive lasing and multistability are to a decrease in ferromagnetic polarization. 
As a further detrimental effect, I consider here an environmental temperature which is large when compared to the resonator frequency, i.e., $T \gg \omega_0$. This condition is usually fulfilled for mechanical resonators such as carbon nanotubes (CNTs), which are a good candidate for an experimental realization of our system, see below. Indeed, their low mechanical frequency ($\omega_0/2\pi \approx \SI{100}{MHz}$) requires taking into account thermal fluctuations. Conversely, for microwave cavities in the GHz regime, $n_B = 0$ is usually a good approximation at cryogenic temperatures. 
In @fig:qdlaser:finite-t-p, I report the numerical calculation of the average occupation of the oscillator in the steady state, together with the stability diagram for a nonideal case ($T = 10\omega_0$ and $P=0.5$):  The qualitative picture is not destroyed, as lasing and multistability are retained. More specifically, the lasing threshold is pushed to a larger coupling strength, according to Eq. (-@eq:qdlaser:rwa-pn-saturation-thr), as well as the onset of bi- and multistability. The thermal noise smears out the transitions to the lasing state. 

### Effect of spin relaxation

The spin the quantum dot may be subject to spin relaxation processes due to interactions with bulk phonons in the substrate [@Khaetskii2000], or through spin-orbit coupling in CNTs [@Churchill2009;@Rice2013]. The overall effect of spin relaxation is to cause spin flip without photon emission into the cavity, as it does not come from coherent energy exchange. Its role is therefore similar to the effect of finite polarization. To model the spin relaxation, I assume a typical energy relaxation time $T_1$, but neglect a general inhomogeneous pure dephasing term of characteristic timescale $T_\phi$. This is justified as this term generally arises in CNTs from hyperfine coupling of the electronic spin to the nuclear spin of $^{13}$C atoms, whose natural abundance in carbon is less than 1% [@Churchill2009]. Taking the spin relaxation rate to be $\gamma_\mathrm{sr} = T_1^{-1}$, I simply add a dissipator $\mathcal{L}_\mathrm{sr} \rho = \gamma_\mathrm{sr} \mathcal{D}(\sigma_-)[\rho]$ to the master equation (-@eq:qdlaser:lindblad-equation) and search for the stationary state.
It is reasonable to assume that if $\gamma_\mathrm{sr}$ is much smaller than the generalized Rabi frequency $gA$ and of the
tunneling rates $\Gamma$, the resonator dynamics is expected to be unperturbed.
In @fig:qdlaser:spin-rel(a), we observe how lasing mechanics is noticeably suppressed for $\gamma_\mathrm{sr}/\omega_0 = \num{e-2}$. For the case of a CNTQD setup, the relaxation time in single-walled CNTs was reported to be $T_1 \approx \SI{100}{\micro\second}$ at $T = \SI{4}{\kelvin}$ corresponding to a relaxation rate of $\SI{10}{\kilo\hertz}$. At low temperature ($T \approx \SI{20}{\milli\kelvin}$), one can expect a substantial decrease of this value. In Figs. -@fig:qdlaser:spin-rel(b)-(c) I report the average cavity occupation and the stability diagram for $\gamma_\mathrm{sr} = \num{e-3}\omega_0$.

![Effect of spin relaxation in the quantum dot. (a) Average phonon occupation on resonance as a function of $g$, for increasing values of spin relaxation rate $\gamma_\mathrm{sr}$. (b) Average phonon occupation as a function of $\Delta\epsilon$ and $g$, for $\gamma_\mathrm{sr} = \num{e-3}\omega_0$. (c)Stability diagram of the resonator for $\gamma_\mathrm{sr} = \num{e-3}\omega_0$. The italic numbers indicate the number of distinct peaks in the Fock distribution. Parameters: $Q = \num{e3}$, $\Gamma = 0.1\omega_0$, $P=0.5$, and $T =10\omega_0$.](figures/qdlaser-spin-rel.pdf){#fig:qdlaser:spin-rel}

### Effect of Duffing nonlinearity

Nanomechanical resonators such as suspendend carbon nanotubes can be intrinsically weakly nonlinear, a property which is expected to play an important role when the amplitude of oscillations is amplified as described by our results. I include a  nonlinear term of Duffing type [@Strogatz2015] into Hamiltonian (-@eq:qdlaser:system-hamiltonian), which is modified into

$$
H = \frac{\Delta\epsilon}{2}\sigma_z + \omega_0 b^\dagger b + \frac{\tilde{\beta}}{4} (b+b^\dagger)^4 + g (\sigma_+ + \sigma_-) (b+b^\dagger).
$$
{#eq:qdlaser:duffing-hamiltonian}

I have introduced $\tilde{\beta} = \beta x_\mathrm{ZPM}^4$, where $\beta$ is the Duffing nonlinearity parameter and $x_\mathrm{ZPM} = \sqrt{\hbar/2m\omega_0}$ is the zero-point amplitude of the oscillator ($\hbar$ has been temporarily restored for clarity).
 To realistically estimate $\tilde{\beta}$, let us consider a carbon nanotube with mass approximately  given by $m = \pi L d/ 1315$, where $L$ and $d$ are length and diameter in meters, respectively [@Laurent2010;@Liu2008]. A typical mass of $m \approx \SI{e-21}{kg}$ gives zero-point fluctuations of order $x_\mathrm{ZPM} \approx \SI{10}{\pico\meter}$ for a frequency $\omega_0/2\pi = \SI{100}{MHz}$. Experimentally, the geometrical nonlinearity parameter for a nanotube is positive ($\beta>0$) and of order $\beta/m = \SI{e35}{\newton\per\kilogram\per\cubic\meter}$ [@Steele2009;@Meerwaldt2012]. It follows that $\tilde{\beta}/2\pi \approx \SI{1}{\kilo\hertz}$, i.e., $\tilde{\beta}/\omega_0 \approx \num{e-5}$.  We neglect the electrostatic nonlinearity arising from strong coupling effects between the leads and the nanotube and from single-electron tunneling, which is in general orders of magnitude smaller and is proportional to $\Gamma \ll \omega_0$. By solving the master equation with Hamiltonian (-@eq:qdlaser:duffing-hamiltonian) for the steady state, I report the average cavity occupation as a function of $g$ in @fig:qdlaser:duffing(a). Finally, Figs. -@fig:qdlaser:duffing(b)-(c) show that the combined presence of finite temperature, finite polarization, spin relaxation and Duffing nonlinearity can still preserve the main features of the system (lasing and bistability), despite the largely nonideal case.

![Effect of Duffing nonlinearity on the lasing, for the case of a nanomechanical resonator. (a) Average occupation for the oscillator as a function of $g$ for three different values of the Duffing nonlinearity parameter $\tilde{\beta}$, without spin relaxation ($\gamma_\mathrm{sr} = 0$). (b) Average occupation for the oscillator as a function $g$ and $\Delta\epsilon$ for $\gamma_\mathrm{sr} = \num{e-3}{\omega_0}$ and $\tilde{\beta} = \num{e-4}{\omega_0}$. (c) Stability diagram of the oscillator for $\gamma_\mathrm{sr} = \num{e-3}{\omega_0}$ and $\tilde{\beta} = \num{e-4}{\omega_0}$: the italic numbers indicate the number of distinct peaks in the Fock distribution. The rest of the parameters read: $Q = \num{e3}$, $\Gamma = 0.1\omega_0$, $P = 0.5$, and $T = 10\omega_0$.](figures/qdlaser-duffing.pdf){#fig:qdlaser:duffing}


### Implementations

The spin-resonator interaction lies at the heart of the physics described in this Chapter, and is the main ingredient to be sought after experimentally. Spin-valve-based carbon nanotube quantum dots (CNTQDs) constitute a promising route to implement the model. Indeed, spin-valve effect in CNTs has been demostrated with substantial spin polarization [@Sahoo2005;@Viennot2015], and CNTQDs have inherently small spin relaxation rate [@Rice2013;@Churchill2009] as well as huge quality factors [@Moser2014]. Very recently, suspended CNTs have witnessed a revival interest because of their ability to support self-sustained oscillations with large amplitudes and even bistability, when driven by electron tunneling [@Urgell2019;@Wen2019;@Willick2020;@Yang2020]. 
The spin-resonator interaction in suspended CNTQDs has been the object of theoretical investigations [@Palyi2012;@Stadler2014;@Stadler2015], where an interaction strength $g_\mathrm{exp} \approx \SI{1}{MHz}$ is predicted for a typical resonance frequency $\omega_0/2\pi = \SI{100}{MHz}$. For $Q=\num{e6}$, $P=0.5$, $\Gamma = 0.05\omega_0$, the threshold coupling can be as low as $g_\mathrm{thr} \approx \num{1.6e-4}\omega_0$, well below $g_\mathrm{exp}$.  


Another relevant platform that can be used to realize this single-atom laser is based on spin valves coupled to microwave cavity photons. In the past few years, a great technological effort has been put into reaching an effective strong coupling between spin and photons, using, e.g., CNTQDs and superconducting microwave cavities [@Viennot2015;@Cubaynes2019], silicon-based devices [@Mi2018;@Samkharadze2018;@Borjans2020], and quantum dots in GaAs/AlGaAs heterostructures [@Landig2018;@Landig2019;@Pan2020]. Hence, it is arguable that our device can stimulate further works in this direction.



## Conclusions

To conclude, I have analyzed here a model for a single-atom laser using a quantum-dot spin valve. This single-atom laser relies on a spin-polarized current pump from ferromagnetic leads into a quantum dot, where the spin is coupled to the motion of a resonator. 
I have shown how the spin-resonator coupling yields a very efficient photon emission in the resonator, such that the conventional rotating-wave approximation (RWA) breaks down, despite the relative weak bare coupling achievable between spin and cavity. 
As a consequence of RWA-breakdown, the resonator can develop a multistable behavior, which is characterized by the presence of multiple peaks in its Fock distribution, and can be detected through simple current measurements: 
As each stable amplitude carries a definite current value, time-monitoring of the current shows telegraph switching between different current plateaus, as long as each stable state is sufficiently long-lived. 
I have analyzed the stability of the system using numerical simulations, supporting the results with a semiclassical explanation;
similarly, the numerical experiments showing telegraph dynamics in the current have been confirmed by a simple two-state classical model, which assumes that the resonator state can be approximated by two well-distinct states with similar occupation probabilities. 
In view of an experimental realization, I have shown how the lasing and the multistability still survive when nonideal scenarios are considered, such as high environmental temperature, finite polarization of the leads, Duffing nonlinearities for the case of a mechanical resonators, and spin-relaxation processes in the quantum dots.
From the above analysis, it is arguable that this system represents a unique platform for investigating the coherent dynamics of a spin coupled to a resonator *beyond* the RWA, but without the need to access ultrastrong couplings. Further, our work raises a range of interesting theoretical and experimental questions about the extent to which multistable lasing can be controlled and exploited, e.g., in nonlinear amplifiers or force sensing devices.




