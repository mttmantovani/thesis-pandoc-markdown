# Single-photon pump by Cooper-pair splitting {#sec:cps}

*The results presented in this Chapter are published in Ref. \onlinecite{Mantovani2019a}.*

## Introduction 

One of the most remarkable features of hybrid quantum-dot-based systems is the
possibility of engineering correlations between distant subsystems and giving
rise to nonlocal effects, which are inherently quantum mechanical. A relevant
example is represented by the Cooper pairs in a superconductor, which are a
natural source of entangled electrons. Cooper pairs can be split using a
Cooper-pair splitter (CPS), a device widely studied both theoretically
[@Recher2001;@Chevallier2011;@Rech2012;@Scherubl2014;@Trocha2015;@Nigg2015;@Schroer2015;@Probst2016;@Dominguez2016;@Amitai2016;@Hussein2016;@Wrzesniewski2017;@Walldorf2018;@Bocian2018;@Wrzesniewski2020]
and
experimentally [@Hofstetter2009;@Herrmann2010;@Hofstetter2011;@Das2012;@Schindele2012;@Schindele2014;@Fulop2014;@Fulop2015;@Tan2015;@Borzenets2016;@Baba2018], while retaining their quantum correlations. It has been shown
recently that these correlations can also give rise to peculiar thermoelectric
effects [@Cao2015;@Sanchez2018;@Hussein2019;@Kirsanov2019]. On the other hand,
I have shown in the previous Chapters how mesoscopic QED devices have proven to be excellent
tools to couple electronic and bosonic degrees of freedom, enabling, in turn,
the possibility to control heat exchange [@Hofer2016;@Ronzani2018;@Dambach2019;@Senior2020;@Thomas2020a]. The purpose of this Chapter is to
conceive and analyze a device that bridges the existing gap between the study of
heat flows in quantum-dot-based [@Erdman2017;@Erdman2018;@Dutta2020] and circuit-QED
devices.

I consider a CPS consisting of a superconducting lead, which is put in the
proximity of two quantum dots. Each dot is further linearly coupled to a local
harmonic resonator, constituted by a microwave cavity or a nanomechanical
oscillator.  The system is sketched in @fig:cps:system. As I show below, it is
possible to obtain full control over heat and photon exchange between the two
resonators by merely tuning the energy level of the quantum dots to match
internal resonances of the coupled system. I will first present the model
Hamiltonian and the master equation to describe the behavior of the system
numerically. Next, I will develop an effective Hamiltonian formalism that
explains with simplicity the induced nonlocal coupling between the cavities due
to the superconductor, and the resulting heat exchange. Finally, I will quantify
and discuss the efficiency of the heat transfer by computing an appropriate
figure of merit.




## Double-quantum-dot Cooper-pair splitter coupled to resonators

### Cooper-pair splitter in quantum-dot setup

![Schematics of a Cooper-pair splitter with two quantum dots coupled to local resonators.](./figures/cps-system.pdf){#fig:cps:system}

In order to write the full Hamiltonian of the system, I first briefly derive the
effective Hamiltonian for a CPS in the limit of large superconducting gap $\Delta$ [@Rozhkov2000;@Eldridge2010]. A superconductor
is tunnel-coupled to two quantum dots with spin-degenerate levels, each in contact with a metallic lead
through a tunneling barrier (see @fig:cps:system). In principle, one can include
finite intradot ($U_\alpha$) and interdot ($U$) Coulomb repulsion strengths. The
Hamiltonian for the dots thus reads:

$$
H_{\mathrm{DQD}}=\sum_{\alpha=L, R}\left[\epsilon_{\alpha} \sum_{\sigma}
d_{\alpha \sigma}^{\dagger} d_{\alpha \sigma}+U_{\alpha} N_{\alpha \uparrow}
N_{\alpha \downarrow}\right]+U \sum_{\sigma \sigma^{\prime}} N_{L \sigma}
N_{R \sigma^{\prime}}.
$$ {#eq:cps:dots-hamiltonian}

The index $\alpha=L,R$ labels the left and right dot, respectively;
$d_{\alpha\sigma}$ is the fermionic annihilation operator for an electron in the
$\alpha$-dot with spin $\sigma$, energy $\epsilon_\alpha$ and corresponding number operator $N_{\alpha
\sigma}$. 
The three leads are indexed by $\chi = L, R, S$ (left, right, and
superconducting). We can write the Hamiltonian as follows:

$$
H_{\chi}=\sum_{k \sigma} (\varepsilon_{\chi k} - \mu_\chi) c_{\chi k \sigma}^{\dagger} c_{\chi k \sigma}-\delta_{\chi, S} \Delta \sum_{k}\left(c_{\chi-k_{\downarrow}} c_{\chi k \uparrow}+\mathrm{H.c.}\right).
$$ {#eq:cps:leads-hamiltonian}

The normal leads are described by Fermi liquids with fermionic annihilation
operators $c_{\chi k \sigma}$. They are kept at a fixed chemical potential $\mu_\chi$. The superconductor is described by the BCS term
with pairing potential $\Delta$ (assumed real and positive without loss of
generality), and its chemical potential is taken as the reference point, $\mu_S = 0$.
Next, we take into account the tunneling between the three leads and the quantum
dots:

\begin{align}
H_{\mathrm{tun},N} &=\sum_{k \sigma}\left[V_{N
L} c_{L k \sigma}^{\dagger} d_{L\sigma}+V_{N R} c_{R k \sigma}^{\dagger}
d_{R\sigma}+\mathrm{H.c.}\right], \label{eq:cps:normal-tunneling-hamiltonian}\\
H_{\mathrm{tun},S} &= \sum_{k \sigma}\left[V_{S L} c_{S k
\sigma}^{\dagger} d_{L \sigma}+V_{S R} c_{S k \sigma}^{\dagger} d_{R
\sigma}+\mathrm{H.c.}\right].  \label{eq:cps:sc-tunneling-hamiltonian}
\end{align}

The tunneling amplitudes $V_{N\alpha}$ and $V_{S\alpha}$ have been introduced.The CPS Hamiltonian is then given by $H_\mathrm{CPS}  = H_{\mathrm{DQD}} +
\sum_\chi H_\chi + H_{\mathrm{tun},N} + H_{\mathrm{tun},S}$.

### Large-gap and large-Coulomb-repulsion limit

Before introducing the master equation that governs the system dynamics, we consider the limits of large gap,
$\Delta \rightarrow \infty$, and large intradot Coulomb repulsion, $U_\alpha \rightarrow \infty$. 

The first approximation allows us to
integrate out the superconducting degrees of freedom in the Hamiltonian, as the
quasiparticles are inaccessible and no dissipation mechanism is introduced.
We then obtain for the CPS:

$$
H^\mathrm{eff}_\mathrm{CPS} = H_\mathrm{DQD} -\sum_{\alpha=L, R} \frac{\Gamma_{S
\alpha}}{2}\left(d_{\alpha\uparrow}^{\dagger} d_{\alpha
\downarrow}^{\dagger}+\mathrm{H.c}\right) - \frac{\sqrt{\Gamma_{S L} \Gamma_{S
R}}}{2}\left(d_{R \uparrow}^{\dagger} d_{L \downarrow}^{\dagger}-d_{R
\downarrow}^{\dagger} d_{L \uparrow}^{\dagger}+\mathrm{H.c.}\right).
$$ {#eq:cps:cps-large-gap-hamiltonian}

In Eq. (-@eq:cps:cps-large-gap-hamiltonian), we identify two contributions to the
DQD Hamiltonian: The first term describes local-Andreev reflection [@Rozhkov2000;@Gramich2015], whereby a
Cooper pair can be transferred from the superconductor to a quantum dot (and
vice-versa), with amplitude $\Gamma_{S\alpha} = 2\pi \rho_S |V_{S\alpha}|^2$
($\rho_S$ is the density of states at the Fermi level of the superconductor,
assumed energy-independent). The second term describes cross-Andreev reflection:
Cooper pairs can be split and recombined between the central superconductor
and both quantum dots, building up nonlocal correlations between the two. The
effective cross-Andreev coupling will be denoted by $\Gamma_S =
\sqrt{\Gamma_{SL} \Gamma_{SR}}$.

The large-$U_\alpha$ approximation is usually fulfilled in experiments since the
charging energy of the dots can reach several meV and is therefore larger than
all other energy scales in the system [@Herrmann2010]. Physically, this means that double
occupation in the dots is energetically inaccessible, and the local-Andreev
reflection contributions in Eq. (-@eq:cps:cps-large-gap-hamiltonian) can be
neglected. The CPS Hamiltonian becomes

$$
H^\mathrm{eff}_\mathrm{CPS} = H_\mathrm{DQD} - \frac{\Gamma_S}{2}\left(d_{R \uparrow}^{\dagger} d_{L \downarrow}^{\dagger}-d_{R
\downarrow}^{\dagger} d_{L \uparrow}^{\dagger}+\mathrm{H.c.}\right).
$$ {#eq:cps:cps-large-gap-large-coulomb-hamiltonian}

Notice that, for large $U_\alpha$, the relevant Hilbert spaces for the double
dot system is spanned by nine states: The empty state $|0\rangle$, four singly
occupied states $|\alpha \sigma \rangle = d^\dagger_{\alpha \sigma}|0\rangle$,
the singlet state $|S\rangle=\frac{1}{\sqrt{2}}\left(d_{R}^{\dagger}
d_{L\downarrow}^{\dagger}-d_{R \downarrow}^{\dagger}
d_{L\uparrow}^{\dagger}\right)|0\rangle$, the unpolarized triplet state $|T
0\rangle=\frac{1}{\sqrt{2}}\left(d_{R\uparrow}^{\dagger}
d_{L\downarrow}^{\dagger}+d_{R\downarrow}^{\dagger}
d_{L\uparrow}^{\dagger}\right)|0\rangle$ and two polarized triplet states $|T
\sigma\rangle=d_{R \sigma}^{\dagger} d_{L \sigma}^{\dagger}|0\rangle$. In the
subgap limit, the cross-Andreev reflection term hybridizes only the states
$|0\rangle$ and $|S\rangle$ of the
double dot: Only Cooper pairs in singlet states can be split and recombined with
the coherent amplitude $\Gamma_S$. The hybridization gives rise to the Andreev
bound states

$$
\begin{split}
|+\rangle &= \cos \frac{\theta}{2}|0\rangle + \sin \frac{\theta}{2} |S\rangle,\\
|- \rangle &=-\sin \frac{\theta}{2}|0\rangle +\cos\frac{\theta}{2} |S\rangle,
\end{split}
$$ {#eq:cps:andreev-states}

with the mixing angle $\theta = \arctan\left( \frac{\sqrt{2}\Gamma_S}{\epsilon_L +
\epsilon_R+U} \right)$. The energy splitting between the Andreev states is $\delta =
\sqrt{(\epsilon_L + \epsilon_R+U)^2 + 2 \Gamma_S^2}.$


### Capacitive coupling to resonators

Finally, we include a capacitive charge coupling
between each quantum dot and a local harmonic resonator of frequency
$\omega_\alpha$ and bosonic field operator $b_\alpha$
[@Stockklauser2015;@Deng2015;@Okazaki2016;@Mi2017;@Li2018]. Assuming a
negligible interdot Coulomb interaction $U$, the system
Hamiltonian that I will consider is then

$$
H = \sum_{\alpha \sigma} \epsilon_\alpha N_{\alpha \sigma} -
\frac{\Gamma_S}{2}(d^\dagger_{R\uparrow} d^\dagger_{L\downarrow} -
d^\dagger_{R\downarrow}d^\dagger_{L\uparrow} + \text{H.c.})
+ \sum_\alpha \omega_\alpha b^\dagger_\alpha b_\alpha + \sum_{\alpha\sigma}
g_\alpha (b_\alpha + b_\alpha^\dagger)N_{\alpha\sigma}.
$$ {#eq:cps:system-hamiltonian}

The strength of the linear interaction between dot and resonator is labeled by
$g_\alpha$.



## Master equation and stationary current

The electron tunneling into the normal leads and the dissipation processes in
the resonators will be treated to lowest order in perturbation theory, such that
the effective reduced dynamics for the coupled CPS-resonators system can be
described by rate equations involving only the populations of the eigenstates of
Hamiltonian (-@eq:cps:system-hamiltonian) 
[@Sauret2004;@Governale2008;@Eldridge2010;@Hussein2016]. From Eq.
(-@eq:cps:normal-tunneling-hamiltonian), we define the following dot-lead
tunneling rates: $\Gamma_L = 2\pi \rho_L |V_{NL}|^2$ and $\Gamma_R = 2\pi \rho_R
|V_{NR}|^2$, where $\rho_L$ and $\rho_R$ are the density of states in the leads
at the Fermi level, assumed constant. The resonators have quality factors
$Q_\alpha = \omega_\alpha / \kappa_\alpha$, where $\kappa_\alpha$ is the cavity
decay rate. To fulfill the sequential-tunneling and the low-damping regime, one
must ensure the conditions $\Gamma_\alpha \ll \{\Gamma_S; k_B T\}$ and
$\kappa_\alpha \ll \{\omega_\alpha; k_B T\}$, where $T$ is the environmental temperature of
the fermionic and bosonic reservoirs. Labeling with $|i\rangle$ the eigenstates
of
Hamiltonian (-@eq:cps:system-hamiltonian), the electronic and bosonic transition rates between two
eigenstates are given by Fermi's golden rule [@Benenti2017]:

\begin{align}
w^{\alpha,s}_{\mathrm{el},j\leftarrow i} &= \Gamma_\alpha f_\alpha^{(s)} (s
E_{ji}) \sum_\sigma |\langle j | d_{\alpha \sigma}^{(s)} | i \rangle |^2,
\label{eq:cps:fermionic-rates} \\
w^{\alpha,s}_{\mathrm{ph},j\leftarrow i} &= s\kappa_\alpha n_B (
E_{ji}) \sum_\sigma |\langle j | b_{\alpha}^{(s)} | i \rangle |^2.
\label{eq:cps:bosonic-rates}
\end{align}

We have introduced the generalized Fermi and Bose functions:

\begin{align}
f_\alpha^{(s)}(x) &= \frac{1}{\exp \left[\frac{s(x - \mu_\alpha)}{k_B T} \right] + 1},
\label{eq:cps:fermi-function} \\
n_B(x) &= \frac{1}{\exp \left[ \frac{x}{k_B T}\right] - 1}, \label{eq:cps:bose-function}
\end{align}

with $s = \pm$. $E_{ji} = E_j - E_i$ denotes the energy difference between two
eigenstates. Notice that it is reasonably implied that the spectrum of Hamiltonian
(-@eq:cps:system-hamiltonian) is nondegenerate. For the sake of compactness in Eqs.
(-@eq:cps:fermi-function)-(-@eq:cps:bose-function), I have provisionally used the
notation $d_{\alpha\sigma}^{(-)}$ [$d_{\alpha \sigma}^{(+)}$] for the fermionic
annihilation (creation) operators. A similar shorthand holds for the bosonic
operators. Let us define the total transition rates as

$$
w_{j\leftarrow i} = \sum_{\alpha s} ( w^{\alpha,s}_{\mathrm{el},j\leftarrow i} +
w^{\alpha,s}_{\mathrm{ph},j\leftarrow i} ).
$$ {#eq:cps:total-rates}

The populations $P_i$ of the system eigenstates will then obey the Pauli-type
master equation

$$
\dot{P}_i = \sum_j w_{i \leftarrow j} P_j - \sum_j w_{j \leftarrow i} P_i.
$$ {#eq:cps:master-equation}

Equation (-@eq:cps:master-equation) admits a stationary solution, represented by
the vector the stationay populations $P_i^\mathrm{st}$. It can be found by
setting $\dot{P}_i = 0$ and consequently solving the resulting linear system.
The stationary vector corresponds to the null space of the matrix $\textbf{W}$
with matrix elements $W_{ij} = w_{i\leftarrow j} - \delta_{ij} \sum_k w_{k \leftarrow j}$. It can be found numerically using standard
linear algebra methods [@Trefethen1997], and adding the normalization
condition for the populations, $\sum_i P_i = 1$. The matrix elements of
$\textbf{W}$ satisfy the sum rule $\sum_i W_{ij} = 0$, i.e., the
sum of the elements in each column of \textbf{W} vanishes.


The steady populations give us access to expression for particle and energy (or heat)
currents [@Benenti2017]. For electrons, they read

 \begin{align}
I_\alpha &= e \sum_{i, j, s} s w^{\alpha,s}_{\mathrm{el}, j \leftarrow i}P_i^\mathrm{st},
\label{eq:cps:electron-current} \\
\dot{Q}_\alpha &= \sum_{i,j,s} E_{ij}
w^{\alpha,s}_{\mathrm{el},j\leftarrow i} P_i^\mathrm{st} - \frac{\mu_\alpha}{e}I_\alpha.
 \label{eq:cps:electron-heat-current}
\end{align}

The corresponding expressions for photon particle and energy currents are
similarly

\begin{align}
I^\mathrm{ph}_\alpha &= \sum_{i, j, s} s w^{\alpha,s}_{\mathrm{ph}, j \leftarrow i}P_i^\mathrm{st},  \label{eq:cps:photon-current} \\
\dot{Q}^\mathrm{ph}_{\alpha} &= \sum_{i,j,s} E_{ij}
w^{\alpha,s}_{\mathrm{ph},j\leftarrow i} P_i^\mathrm{st}.
 \label{eq:cps:photon-heat-current}
\end{align}

For electrons, I have multiplied the particle current in Eq. (-@eq:cps:electron-current) 
by the electron's charge $e$ to get an electric
current.
Notice that, in contrast to the metallic leads, the bosonic reservoirs are at zero chemical potential, so the heat
current coincides with the energy current.

### Large bias limit and stationary current in the absence of resonators {#sec:cps:absence-resonators}

In the rest of the analysis, I will assume that the normal leads are negatively
biased with respect to the Fermi level of the superconductor, such that
$\mu_\alpha = - e V$, where $V>0$ is the applied bias voltage. Furthermore, the
bias voltage is assumed to be very large, such that $\{ U; \Delta \} \gg eV \gg
\{k_B T; \Gamma_S\}$. In this regime, electrons will flow unidirectionally from
the superconductor via the quantum dots into the leads. Transport in the
opposite direction will be blocked, as the rates $w^{\alpha,+}_{\mathrm{el},j
\leftarrow i}$ vanish. Furthermore, the temperature of the leads
becomes irrelevant. The electric current flowing in the lead is simply given by
$I_\alpha = e \Gamma_\alpha \sum_\sigma \langle N_{\alpha \sigma} \rangle$, with
the average performed on the steady state. It is worth noticing that the negative bias also reduces the number of available states in the double dot: The triplet states $|T0\rangle$ and $|T\sigma$ can be excluded from the computation, as they can only appear due to the blocked lead-to-dot tunneling.

In absence of resonators, the stationary solution to master equation (-@eq:cps:master-equation) can be
found analytically. Hence, we can find a simple expression for the stationary
current with largely negatively-biased leads, assuming $\Gamma_L = \Gamma_R = \Gamma$ without loss of generality:

$$
I_L = I_R = e\Gamma \frac{\Gamma_S^2}{(\epsilon_L + \epsilon_R)^2 + 2 \Gamma_S^2}.
$$ {#eq:cps:andreev-current-analytical}

The average current flows equally through both leads, and is pictured in
@fig:cps:andreev-current along the direction of the average dot level,
$\bar{\epsilon} = (\epsilon_L + \epsilon_R)/2$. It is
a lorentzian centered around $\bar{\epsilon} = 0$ with width proportional to $\Gamma_S$, and has a simple
physical meaning: When $\bar{\epsilon}$ lies closely to the Fermi level of the
superconductor, the hybridization between the singlet and empty states is
maximized. Cooper pairs injected in the separate dots can then tunnel into the
corresponding leads through single-electron events generating a current. The
transport window broadens with increasing $\Gamma_S$ as it leads to a stronger hybridization.

![Andreev stationary current through a normal lead, in absence of resonators, as a function
of the average dot level $\bar{\epsilon} = (\epsilon_L + \epsilon_R)/2$ for two
values of $\Gamma_S$.](./figures/cps-andreev-current.pdf){#fig:cps:andreev-current}

The transport mechanism can be better understood with the help of @fig:cps:andreev-cycle. When the leads are negatively biased, the double dot in the singlet state $|S\rangle$ decays into a singly-occupied state ($|\alpha\sigma\rangle$) through a single-electron tunneling event, giving up an electron to one of the leads, at the incoherent rate $\Gamma$. Similarly, the remaining electron will tunnel out, leaving the double dot in the empty state $|0\rangle$, see @fig:cps:andreev-cycle(a). However, the empty and the singlet state are coherently coupled with amplitude $\Gamma_S$, allowing the cycle to restart and giving rise to the steady current pictured in @fig:cps:andreev-current.

An interesting situation appears when the charge states are only *weakly* hybridized, i.e., whenever $|\epsilon| \gtrsim \Gamma_S$. In this case, for positive $\epsilon$, the Andreev bound state $|+\rangle$ is "mostly" the singlet state, while the state $|-\rangle$ is "mostly" empty, see @fig:cps:andreev-cycle(b). This leads to a strong asymmetry of the transition rates $|+\rangle \leftrightarrow |\alpha\sigma\rangle$ and $|-\rangle \leftrightarrow |\alpha\sigma\rangle$. The situation is reversed when $\epsilon <0$, as weak hybridization implies $|+\rangle \approx |0\rangle$ and $|-\rangle \approx |S\rangle$. This strong asymmetry lies at the heart of the photon transfer mechanism that I will describe below.

![(a) For negatively biased leads, subsequent incoherent tunneling events bring the singlet state into the empty state. The two states are coherently coupled due to the superconductors, leading to a stationary current. (b) The coherent coupling leads to the formation of the Andreev bound states $|\pm\rangle$, with energy separation $\delta$. For weak hybridization ($|\epsilon| \gtrsim \Gamma_S$) the transition rates $|+\rangle \leftrightarrow |\alpha\sigma\rangle$ and $|-\rangle \leftrightarrow |\alpha\sigma\rangle$ are strongly asymmetric. Here, the situation for $\epsilon > 0$ is pictured.](figures/cps-andreev-cycle.pdf){#fig:cps:andreev-cycle}

## Polaron transformation and effective Hamiltonian

Before showing the main features of the system through numerical calculations, I
perform a polaron transformation on Hamiltonian (-@eq:cps:system-hamiltonian), which facilitates the
unraveling of the physical processes I want to address. For a given operator
$O$, let us define the unitary transformation

$$
\overline{O} = e^{\xi} O e^{-\xi}, \quad \text{with} \quad \xi =
\sum_{\alpha\sigma}\Pi_{\alpha} N_{\alpha\sigma} \quad \text{and} \quad
\Pi_{\alpha} = \frac{g_\alpha}{\omega_\alpha}(b^\dagger_\alpha - b_\alpha).
$$ {#eq:cps:polaron-transformation}

Equation (-@eq:cps:polaron-transformation) is an example of polaron
transformation [@Lang1962], widely used also in the context of quantum transport
[@Brandes1999;@Brandes2005]. The application of Eq.
(-@eq:cps:polaron-transformation) to Eq. (-@eq:cps:system-hamiltonian) leads to
the polaron-transformed Hamiltonian

$$
\overline{H}=\sum_{\alpha \sigma} \bar{\epsilon}_{\alpha}\left( N_{\alpha \sigma} + \frac{|S\rangle \langle S|}{2}\right) -\frac{\Gamma_{S}}{\sqrt{2}}\left(|S\rangle\langle 0|X+| 0\rangle\langle S| X^{\dagger}\right)+\sum_{\alpha} \omega_{\alpha} b_{\alpha}^{\dagger} b_{\alpha},
$$ {#eq:cps:polaron-transformed-hamiltonian}

with $\bar{\epsilon}_\alpha = \epsilon_\alpha -
\frac{g_\alpha^2}{\omega_\alpha}$ and $X = \exp \left( \sum_\alpha \Pi_\alpha
\right)$. In this representation, the linear interaction term appearing in Eq.
(-@eq:cps:system-hamiltonian) now emerges as a transverse interaction between the
$|S\rangle$ and $|0\rangle$ states and the resonators, to _all_ orders in
$g_\alpha$, through the operator $X$. More intriguingly, this interaction is
purely nonlocal since it vanishes for $\Gamma_S = 0$.

At this point, I make the assumption of small coupling, $g_\alpha \ll \omega_\alpha$. The operators $X$ and $X^\dagger$ can be expanded up to second order in $g_\alpha$, such that the interaction term in Eq. (-@eq:cps:polaron-transformed-hamiltonian) becomes:

$$
\overline{H}_{\mathrm{int}}=-\frac{\Gamma_{S}}{\sqrt{2}}\left[i \sigma_{y} \Pi+\sigma_{x}\left(1+\frac{\Pi^{2}}{2}\right)\right]+\mathcal{O}\left(\Pi^{3}\right),
$$ {#eq:cps:polaron-transformed-interaction}

with the definitions: $i\Pi = \sum_\alpha i \Pi_\alpha$ (generalized total momentum), $\sigma_x = |0\rangle \langle S | + \text{H.c.}$, and $\sigma_y = -i |0\rangle \langle S | + \text{H.c.}$ By diagonalizing the interaction in the electronic subspace, we obtain the renormalized hybridized states 

\begin{align}
    |\bar{+}\rangle &=\cos \left(\frac{\bar{\theta}}{2}\right)|0\rangle+\sin \left(\frac{\bar{\theta}}{2}\right)|S\rangle, \\ 
    |\bar{-}\rangle&=-\sin \left(\frac{\bar{\theta}}{2}\right)|0\rangle+\cos \left(\frac{\bar{\theta}}{2}\right)|S\rangle,
\end{align}

with $\bar{\theta}=\arctan \left(\frac{\sqrt{2} \Gamma_{S}}{\bar{\epsilon}_{L}+\bar{\epsilon}_{R}}\right)$. The renormalized energy splitting then reads
$\bar{\delta} = \sqrt{(\bar{\epsilon}_L + \bar{\epsilon}_R)^2 + 2 \Gamma_S^2}$. Using this new basis, we define the Pauli algebra with matrices $\tau_+ = |\bar{+}\rangle\langle\bar{-}|$, $\tau_- = (\tau_+)^\dagger$, $\tau_x = \tau_+ + \tau_-$, $\tau_y =  -i(\tau_+ - \tau_-)$, $\tau_z = [\tau_+, \tau_-]$, through which Eq. (-@eq:cps:polaron-transformed-hamiltonian) up to second order becomes

$$
\overline{H} = \sum_{\alpha \sigma} \bar{\epsilon}_{\alpha} N_{\alpha \sigma}+\frac{\bar{\delta}}{2} \tau_{z}+\sum_{\alpha} \omega_{\alpha} b_{\alpha}^{\dagger} b_{\alpha} -\frac{\Gamma_{S}}{2 \sqrt{2}}\left[2 i \tau_{y} \Pi+\left(\sin \bar{\theta} \tau_{z}+\cos \bar{\theta} \tau_{x}\right) \Pi^{2}\right]+\mathcal{O}\left(\Pi^{3}\right).
$$

We now move to the interaction picture with respect to the free Hamiltonian $H_{0}=\sum_{\alpha \sigma} \bar{\epsilon}_{\alpha} N_{\alpha \sigma}+\frac{\delta}{2} \tau_{z}+ \sum_{\alpha} \omega_{\alpha} b_{\alpha}^{\dagger} b_{\alpha}.$  We obtain the time-dependent interaction Hamiltonian

$$
\begin{aligned} 
    \overline{H}_{\mathrm{int}}(t) =&-\sum_{\alpha} \frac{g_{\alpha} \Gamma_{S}}{\omega_{\alpha} \sqrt{2}}\left(e^{i \omega_{\alpha} t} b_{\alpha}^{\dagger}-e^{-i \omega_{\alpha} t} b_{\alpha}\right)\left(e^{i \bar{\delta} t} \tau_{+}-e^{-i \bar{\delta} t} \tau_{-}\right) \\ &-\frac{\Gamma_{S} g_{L} g_{R}}{\sqrt{2} \omega_{L} \omega_{R}}\left[e^{i \Omega t} b_{L}^{\dagger} b_{R}^{\dagger}+e^{-i \Omega t} b_{L} b_{R}-e^{i(\Delta \omega) t} b_{L}^{\dagger} b_{R}-e^{-i(\Delta \omega) t} b_{L} b_{R}^{\dagger}\right] \\ &\times \left[\sin (\bar{\theta}) \tau_{z}+\cos (\bar{\theta})\left(e^{i \bar{\delta} t} \tau_{+}+e^{-i \bar{\delta} t} \tau_{-}\right)\right] \\ &-\sum_{\alpha} \frac{\Gamma_{S} g_{\alpha}^{2}}{2 \sqrt{2} \omega_{\alpha}^{2}}\left[e^{2 i \omega_{\alpha} t}\left(b_{\alpha}^{\dagger}\right)^{2}+e^{-2 i \omega_{\alpha} t} b_{\alpha}^{2}-2 b_{\alpha}^{\dagger} b_{\alpha}-1\right] \\ & \times\left[\sin (\bar{\theta}) \tau_{z}+\cos (\bar{\theta})\left(e^{i \bar{\delta} t} \tau_{+}+e^{-i \bar{\delta} t} \tau_{-}\right)\right]+\mathcal{O}\left(g_{\alpha}^{3} / \omega_{\alpha}^{3}\right). 
\end{aligned}
$$ {#eq:cps:polaron-transformed-interaction-second-order}

I have introduced the sum, $\Omega = \omega_L + \omega_R$, and the difference, $\Delta\omega = \omega_L - \omega_R$, of the resonator frequencies. Equation (-@eq:cps:polaron-transformed-interaction-second-order)
is the central object that contains the resonant processes I will discuss below using a suitable rotating-wave approximation.

## Simultaneous ground-state cooling of nanoresonators {#sec:cps:local-cooling}

![Simultaneous cooling of resonators due to cross-Andreev reflection. (a) For two identical resonators, when $|\epsilon| \gtrsim \Gamma_S$ and $\epsilon > 0$, the weak charge hybridization leads to an asymmetry in the transition rates $|\bar{+}\rangle \leftrightarrow |\alpha\sigma \rangle$ and $|\bar{-}\rangle \leftrightarrow |\alpha\sigma \rangle$. Fast rates are depicted by blue solid arrows, slow transition with dashed light blue arrows. The coherent effective coupling \[Eq. (-@eq:cps:local-effective-hamiltonian), curved dashed arrows\] mixes states with different photon numbers close to the resonance $\bar{\delta} = \omega$. In the steady nonequilibrium state of the system, the cavities are effectively cooled down. (b) Current $I_\alpha$ as a function of the on-site dot energy $\epsilon$ at zero (red dashed line) and finite temperature (solid black line). (c) Average cavity occupation in one of the cavities, for $k_B T = 5\omega$. The horizontal dotted line corresponds to the thermal occupation. Inset: Photon occupation at the cooling resonance, $\epsilon = \epsilon_c$, as a function of $\Gamma_S$, for two different values of the electron tunneling rate $\Gamma$. The curves are rescaled to the thermal occupation. Parameters: $\Gamma  = 2 \times 10^{-4} \omega,\ g = 0.02\omega,\ Q = 10^5,\ \Gamma_S=0.2\omega$. ](figures/cps-local-cooling.pdf){#fig:cps:local-cooling}

The first main feature offered by our system is the possibility to obtain *simultaneous* ground-state cooling of the resonators, as well as simultaneous heating. To achieve this, we tune the dots' energy levels to the same value $\epsilon_L = \epsilon_R = \epsilon$ and assume two identical resonators, with $\omega_L = \omega_R = \omega$, $g_L = g_R = g$, and $Q_L = Q_R = Q$. Furthermore, we move close to the resonance condition $\bar{\delta} = \omega$. Notice that this is fulfilled by two values of $\epsilon$ of opposite sign, namely $\epsilon = \pm \sqrt{\omega^2 - 2 \Gamma_S^2}$. Close to the resonance, we can perform a rotating-wave approximation in Eq. (-@eq:cps:polaron-transformed-interaction-second-order), obtaining to first order in $g/\omega$ the simple, time-independent interaction

$$
\overline{H}_\text{loc} = \sum_\alpha \frac{g}{2} \sin \bar{\theta} (b_\alpha \tau_+ + b^\dagger_\alpha \tau_-).
$$ {#eq:cps:local-effective-hamiltonian}

Equation (-@eq:cps:local-effective-hamiltonian) describes hopping between the Andreev states $|\bar{+}\rangle$ and $|\bar{-}\rangle$ associated with one-photon loss and absorption in the cavities, through a Jaynes-Cummings interaction. Notice that the effective coupling is proportional to $\sin \bar{\theta} = \sqrt{2}\Gamma_S / \bar{\delta}$, and is therefore a direct consequence of nonlocal Andreev reflection. 

We can illustrate how this effective interaction leads to ground-state cooling of both cavities with the help of @fig:cps:local-cooling(a). The interaction (-@eq:cps:local-effective-hamiltonian) coherently mixes the states $|\bar{+}, n_{L}, n_{R}\rangle,\ |\bar{-}, n_{L}+1, n_{R}\rangle,$ and $|\bar{-}, n_{L}, n_{R}+1\rangle$, which are degenerate for $\overline{H}_\text{loc} = 0$. When $|\epsilon| \gtrsim \Gamma_S$, the hybridization of the states $|0\rangle$ and $|S\rangle$ is weak. What happens now depends on the sign of $\epsilon$: If $\epsilon < 0$, then $|\bar{+}\rangle \approx |0\rangle$ and $|\bar{-}\rangle \approx |S\rangle$. Conversely, for $\epsilon > 0$, $|\bar{+}\rangle \approx |S\rangle$ and $|\bar{-}\rangle \approx |0\rangle$. Let us consider the latter case. As discussed in @sec:cps:absence-resonators, the weak hybridization leads to the *fast* transitions $|\bar{+} \rangle \rightarrow |\alpha\sigma\rangle \rightarrow |\bar{-}\rangle$ and the *slow* transitions $|\bar{+}\rangle \leftarrow |\alpha\sigma \rangle \leftarrow |\bar{-} \rangle$, which conserve the photon number. When an electron reaches the state $|\bar{-}\rangle \approx |0\rangle$, the coherent cycle restarts due to the effective coupling, leading to further $|\bar{+}\rangle \rightarrow |\alpha\sigma \rangle \rightarrow |\bar{-} \rangle$ transitions. During each cycle, however, a boson is subtracted from both cavities, see @fig:cps:local-cooling(a). The process continues until the loss mechanisms in the resonators lead to a steady nonequilibrium state, which leaves the cavities cooled down. The opposite occurs for $\epsilon <0$: In this case, the cavities are heated as in each coherent cycle a photon is pumped into them.

To prove the qualitative description I have just presented, I show in @fig:cps:local-cooling(b)-(c) the electric current through one lead and the average photon occupation, $\bar{n}_\alpha = \langle b^\dagger_\alpha b_\alpha \rangle$, for one cavity as a function of $\epsilon$, finding the numerical steady solution of Eq. (-@eq:cps:master-equation) with the full Hamiltonian (-@eq:cps:system-hamiltonian). The current presents the characteristic broad resonance described in @sec:cps:absence-resonators. However, at zero temperature, the inelastic peak appears at $\epsilon = -\sqrt{\omega^2 - 2 \Gamma_S^2}$, which is associated to heating. At finite temperature, the second sideband peak emerges at positive $\epsilon$, associated with cavity cooling.

In the inset of @fig:cps:local-cooling(c) I further demonstrate how the cooling can be very effective, robustly bringing the resonators into their ground state, as witnessed by the value of $\bar{n}_\alpha \ll 1$. The inset shows an optimal cooling region over a wide range of values of $\Gamma_S$, which is the result of the interplay between the effective coupling $g/2\sin\bar{\theta}$, vanishing for $\Gamma_S \rightarrow 0$, and the hybridization between $|0\rangle$ and $|S\rangle$: A strong hybridization (achieved with large $\Gamma_S$) reduces the asymmetry between the transition rates $|\pm\rangle \leftrightarrow |\alpha\sigma\rangle$, deteriorating the cooling effect.

## Photon transfer between resonators

The local cooling mechanism presented in @sec:cps:local-cooling is due to one-photon resonances, and can be simply explained by analyzing the first-order term of Eq. (-@eq:cps:polaron-transformed-interaction-second-order). By keeping the second-order term, it is possible to describe the resonance $\bar{\delta} = \omega_L - \omega_R$. Let us assume, without loss of generality, $\omega_L > \omega_R$. For $\bar{\delta} = \omega_L - \omega_R$, the rotating-wave approximation yields the effective coupling Hamiltonian

$$
H_{\mathrm{RWA}}^{(-)}= \sum_{\alpha} \frac{\Gamma_{S} g_{\alpha}^{2}}{2 \sqrt{2} \omega_{\alpha}^{2}}\left(2 n_{\alpha}+1\right) \sin \bar{\theta} \tau_{z} +g_{\mathrm{NL}}\left(b_{L}^{\dagger} b_{R} \tau_{-}+\mathrm{H.c.}\right),
$$ {#eq:cps:nonlocal-effective-hamiltonian-difference}

with the effective nonlocal coupling strength

$$
g_{\mathrm{NL}}=\frac{\Gamma_{S} g_{L} g_{R}}{\sqrt{2} \omega_{L} \omega_{R}} \cos \bar{\theta}.
$$ {#eq:cps:nonlocal-effective-coupling}

Notice that, once again, Hamiltonian (-@eq:cps:nonlocal-effective-hamiltonian-difference) is purely nonlocal as it vanishes for $\Gamma_S = 0$. It is composed of two terms. Let us focus first on the second term. It describes coherent processes in which the superconductor mediates photon transfer between the resonators, by coupling the subspaces $\left|\bar{+}, n_{L}-1, n_{R}+1\right\rangle$ and $\left|\bar{-}, n_{L}, n_{R}\right\rangle$ \[see @fig:cps:nonlocal-transfer(a)\]. This term vanishes for identical resonators, as it would require $\bar{\delta} = 0$ and therefore $\Gamma_S = 0$. The mechanism which leads to a steady photon flow between the cavities is similar to the cooling mechanism described in @sec:cps:local-cooling: The unbalance in the transition rates due to the weak charge hybridization facilitates the processes whereby the left resonator loses a photon, which is transferred to the right cavity. In the steady state, a stationary heat flow is established between the oscillators \[see also the inset of @fig:cps:efficiency(b)\]. In @fig:cps:nonlocal-transfer(b), I also report the electronic current calculated numerically with the full Hamiltonian for the case of two resonators with different frequencies. A number of resonances appear, including the one-photon resonances discussed in @sec:cps:local-cooling, the two-photon resonance $\bar{\delta} = \omega_L - \omega_R$, as well as multiphoton resonances which are described in @sec:cps:higher-order-resonances and are a generalization of the process I just described.

The first term in Hamiltonian (-@eq:cps:nonlocal-effective-hamiltonian-difference) is proportional to $n_\alpha \tau_z$, and can be regarded as a dispersive shift of the cavity frequencies, depending on the Andreev bound state. This translates into a fine double-peak structure of the nonlocal resonance, see @fig:cps:efficiency(b), because the quantities calculated using Eq. (-@eq:cps:master-equation) are denstiy-matrix averages.

![(a) Sketch of the nonlocal photon transfer mechanism, occurring when $\bar{\delta} \approx \omega_L - \omega_R$. (b) Electronic current through one lead, as a function of the dots' level $\epsilon$, for two different values of $g = g_L = g_R$. A few resonances described by Eq. (-@eq:cps:higher-order-resonances) are indicated by arrows. Parameters: $\Gamma=10^{-4}\Gamma_S,\ \omega_L=5\Gamma_S,\ \omega_R=3\Gamma_S,\ Q_L=Q_R=10^5,\ k_B T = 5 \Gamma_S$.](figures/cps-nonlocal-transfer.pdf){#fig:cps:nonlocal-transfer}

## Efficiency

To understand the extent to which the system is able to cool the cavities
locally or transfer photons between them, we can quantify the efficiency of
these processes with the help of Eq. (-@eq:cps:photon-heat-current). The quantity
$\dot{Q}^\mathrm{ph}_\alpha$ measures the stationary heat current flowing from the bosonic
reservoir $\alpha$  to the corresponding resonator. It is negative (positive)
when the resonator is cooled (heated), and vanishes for an oscillator in thermal
equilibrium. 
For the case of local cooling, dividing this quantity by the resonator frequency $\omega_\alpha$
gives us an estimate of the energy quanta lost by the resonator on average per
unit time. We can then compare this number to the rate at which Cooper pairs are
injected into the system, which is given by $\frac{|I_S|}{2e}$. The Andreev
current through the superconductor $I_S = -I_L - I_R$ follows from current
conservation. The natural definition of the local cooling efficiency for the
$\alpha$-resonator is then

$$
\eta^{(\alpha)}_\mathrm{loc} =
\frac{|\dot{Q}^\mathrm{ph}_\alpha|}{\omega_\alpha} \frac{2e}{|I_S|}.
$$ {#eq:cps:local-efficiency}

In a similar fashion, when around the resonance $\bar{\delta} = \omega_L -
\omega_R$, the photon transfer efficiency is defined as

$$
\eta_\mathrm{NL} =
\frac{|\dot{Q}^\mathrm{ph}_L - \dot{Q}^\mathrm{ph}_R|}{\omega_L - \omega_R} \frac{2e}{|I_S|}.
$$ {#eq:cps:nonlocal-efficiency}

The results for Eqs. (-@eq:cps:local-efficiency) and
(-@eq:cps:nonlocal-efficiency) are shown in @fig:cps:efficiency, close to the
corresponding resonances. In both cases, the efficiency can reach 90\%: Almost one photon is absorbed from each cavity (local cooling) or transferred
from the left to the right cavity (nonlocal transfer) for each Cooper pair split
into the dots. 

The efficiency is essentially limited by two factors: (i) an elastic
contribution to the current (exemplified by the results of @sec:cps:absence-resonators) where electrons flow into the leads without
exchanging energy with the cavities; (ii) a finite fraction of the injected
electrons acting against the dominant process (cooling or photon transfer), as
illustrated by the dashed arrows in [@fig:cps:local-cooling(a);@fig:cps:nonlocal-transfer(a)]. Both processes become relevant when
$\Gamma_S$ is increased, and are a byproduct of the finite hybridization between the empty
and the singlet state. However, we remark that the hybridization is necessary to
achieve a nonzero efficiency.

![(a) Local cooling efficiency for the left cavity in the vicinity of the resonance $\bar{\delta} \approx \omega_L$. (b) Photon transfer efficiency, around $\bar{\delta} \approx \omega_L - \omega_R$. The inset shows the average cavity occupations, normalized to the thermal occupancy. The parameters are the same as in @fig:cps:nonlocal-transfer.](figures/cps-efficiency.pdf){#fig:cps:efficiency}

## Higher-order resonances {#sec:cps:higher-order-resonances}

The polaron-transformed Hamiltonian (-@eq:cps:polaron-transformed-hamiltonian) reveals that Eqs. (-@eq:cps:local-effective-hamiltonian) and (-@eq:cps:nonlocal-effective-hamiltonian-difference) are in fact two particular cases stemming from a more general resonance condition given by

$$
\bar{\delta} \approx\left|p \omega_{L} \pm q \omega_{R}\right|,
$$ {#eq:cps:higher-order-resonances}

where $p$ and $q$ are nonnegative integers. By expanding the operator $\Pi$ in Eq. (-@eq:cps:polaron-transformed-interaction) up to $n$-th order, one
obtains terms $(b_\alpha)^n$ and $(b_\alpha^\dagger)^n$, which, after moving to the interaction picture and performing a suitable rotating-wave approximation, will yield
$n$-photon local absorption or emission processes. The expansion
contains also terms of the form $(b^\dagger_\alpha)^p (b_{\bar{\alpha}})^q$ and $(b_{\alpha}^{\dagger})^{p}(b_{\bar{\alpha}}^{\dagger})^{q}$ together with their Hermitian conjugates, with the constraint $p+q = n$
($\bar{\alpha} = R$ if $\alpha = L$ and vice versa). The former term describes
the coherent transfer of $|p-q|$ photons between the cavities,
while the latter describes coherent emission and re-absorption
of $p$ and $q$ photons from the $\alpha$ and $\bar{\alpha}$ cavity, respectively. When either $p$ or $q$ is zero, the corresponding resonant process entails local cooling or heating of the cavities. An example for $p=0$ and $q=2$ is indicated in @fig:cps:nonlocal-transfer(b), at the resonance $\bar{\delta} = 2\omega_R$: This is a two-photon resonance at which two photons are absorbed (at $\epsilon <0$) or lost (at $\epsilon > 0$) by the $R$-cavity.


## Conclusions

In this Chapter, I have considered a double-quantum-dot-based CPS, with linear
capacitive coupling to local resonators. The main feature of this device is the
activation of a nonlocal coupling between the cavities, that allows to transfer
photon and control the heat exchange between them with single-photon precision.
This energy flow can also be channeled to cool or heat locally a single cavity.
The origin of this indirect coupling resides in the quantum correlations in the
Cooper pairs of the superconductor, which are maintained even when electrons are
split into the distant quantum dots, affecting the charge-vibration interaction.
As a result, this platform is a versatile tool to inspect quantum
thermodynamical effects involving electronic and bosonic degrees of freedom.
Because the effective interaction is at a single-photon level, one can also
devise few-photon or -phonon control and manipulation techniques
[@Okamoto2013;@Zhu2017;@Zhang2020] by implementing time-dependent protocols on the gate
voltages acting on the dots: By moving around the interesting resonances, one
can tune the strength of the nonlocal interactions dynamically. 

To conclude, I briefly point out the experimental practicability of the setup.
The capacitive coupling $g_\alpha/2\pi$ between quantum dots and microwave
resonators can be as high as 100 MHz with quality factors $Q \sim 10^4$ and
frequencies $\omega_\alpha/2\pi \approx 7$ GHz [@Bruhat2016;@Stockklauser2017].
For mechanical resonators, the typical coupling strength is only $\sim\!100$ kHz,
but for much lower resonance frequencies in the MHz regime and larger quality
factors up to $10^6$ [@Okazaki2016;@Moser2014]. The ratio
$g_\alpha/\omega_\alpha \sim 10^{-2}$ assumed in the analysis above is
then experimentally achievable. Finally, the magnitude of the cross-Andreev
reflection rate in a CPS is approximately $\Gamma_S \lesssim \sqrt{\Gamma_{SL}\Gamma_{SR}}$, when the distance between
the dots is much shorter than the coherence length in the superconducting
contact. The local Andreev reflection rates $\Gamma_{S\alpha}$ can be several tens of
microelectronvolts [@Gramich2015], which is orders of magnitude lower than the superconducting
gap $\Delta$ but comparable to the energy of microwave radiation, as assumed in
the calculations.