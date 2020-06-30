# Introduction {#sec:introduction}


## Coupling light and matter: Cavity quantum electrodynamics

One of the most succesful and ever-developing applications of the quantum theory
is the investigation of the interaction between *light* and *matter*.
Historically, the first experiments involved real atoms
interacting with the quantized electromagnetic field of optical and microwave
cavities, within the framework commonly
referred to as cavity quantum electrodynamics (QED)
[@Raimond2001;@Walther2006;@Haroche2006].
 
The paradigmatic model for the light-matter interaction is represented by the Rabi Hamiltonian [@Rabi1936;@Rabi1937], given by

$$
  H_\text{Rabi} = \frac{\Delta\epsilon}{2}\sigma_z + \omega_0 b^\dagger b + g
  (\sigma_+ + \sigma_-) (b + b^\dagger).
$$
{#eq:intro:rabi-hamiltonian}

It describes the coherent interaction between a two-level system (the *matter*, e.g., two levels of an atom or a spin-1/2) and a quantum harmonic oscillator (the *light*, e.g., a single mode of a microwave or optical cavity), see @fig:intro:cqed. The Rabi Hamiltonian follows from the interaction between the dipole moment of an atom and the electric field of the cavity, in the dipole approximation [@Scully1997]. In Eq. (-@eq:intro:rabi-hamiltonian), $\Delta \epsilon$ is the energy difference between the two levels and $\omega_0$ is the resonance frequency of the harmonic oscillator. The constant $g$ quantifies the strength of the interaction between light and matter.

![Sketch of a cavity QED setup. Two atomic levels of an atom, $|e\rangle$ and $|g\rangle$ with energy separation $\Delta \epsilon$ are coupled with strength $g$ to a single mode of a cavity with resonance frequency $\omega_0$. The atom and the cavity are subject to decay at rates $\gamma$ and $\kappa$, respectively, due to the interaction with the environment.](figures/intro-cqed.pdf){#fig:intro:cqed}

The two-level system is composed of an excited ($|e\rangle$) and a ground ($|g \rangle$) state, which can be represented by two column vectors $|e\rangle = (1, 0)^T$ and $|g\rangle = (0, 1)^T$. In the basis $\{|e\rangle, |g\rangle \}$, the qubit is described by the raising and lowering operators $\sigma_+ = |e \rangle \langle g|$ and $\sigma_- = |g \rangle \langle e |$, which generate the Pauli operators $\sigma_x = \sigma_+ + \sigma_-$, $\sigma_y = -i (\sigma_+ - \sigma_-)$ and $\sigma_z = [\sigma_+, \sigma_-]$, with matrix form

$$
\sigma_x = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}, \quad 
\sigma_y = \begin{pmatrix} 0 & -i \\ i &  0 \end{pmatrix}, \quad
\sigma_z = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}.
$$

The harmonic oscillator is described by the bosonic field annihilation and creation operators, $b$ and $b^\dagger$. 

The exact eigenenergies of Hamiltonian (-@eq:intro:rabi-hamiltonian) have been found only several decades after the first works of Rabi, and can only be written as the roots of a complex trascendental function [@Braak2011]. The primary way of simplifying the Rabi Hamiltonian is to assume a *small* value of the coupling strength, $g \ll \{\omega_0; \Delta\epsilon\}$ and consider the two systems to be close to resonance, i.e., $\Delta\epsilon \approx \omega_0$. This permits to apply the rotating-wave approximation (RWA) to obtain the Jaynes-Cummings (JC) Hamiltonian [@Jaynes1963]:

$$
H_\text{JC} = \frac{\Delta\epsilon}{2}\sigma_z + \omega_0 b^\dagger b + g
   (\sigma_+ b + b^\dagger \sigma_-),
$$
{#eq:intro:jc-hamiltonian}

which can be exactly solved, see Ch. -@sec:qdlaser.

The Rabi model has been extended to two-photon
interactions, where $b$ and $b^\dagger$ are replaced by $b^2$ and
$(b^\dagger)^2$ [@Travenec2012], and to $N$-level atoms [@Albert2012].
Similarly, extensions of the JC model with a single cavity mode coupled to many two-level systems (Dicke model [@Dicke1954] and Tavis-Cummings models [@Tavis1968;@Tavis1969]) have been devised and implemented.

When dealing with a cavity QED experiment, an important role is played by the surroundings of the light-matter system: Atoms and cavities are coupled not only to each other---which is desired---but also to the measurement apparati and to the thermal environment. This leads to inevitable relaxation and decoherence in the system, which in most cases can be captured by spontaneous emission rate $\gamma$ for the atom and decay rate $\kappa$ for the cavity (see Ch. -@sec:theory). Based on the intensity of the coupling strength $g$ compared to the decay rates and to the cavity and qubit frequencies, different regimes can be distinguished:

*Weak-coupling regime*, $g \lesssim \{\gamma; \kappa\}$: 

: In this case $g$ is not strong enough to observe coherence, because the energy is lost from the system before it can be exchanged between atom and cavity. Nevertheless, the coupling to the radiation field yields an enhancement of the spontaneous emission rate of the atom, a phenomenon known as Purcell effect [@Purcell1946].

*Strong-coupling regime*, $g \gg \{\gamma;\kappa\}$:

: When $g$ is strong enough to overcome the losses, the coherent exchange of energy between light and matter can be observed. This corresponds to the regime of validity of the JC model. The development of resonators with larger quality factors $Q = \omega_0/\kappa$, associated with lower decay rates, led to the observation of vacuum Rabi oscillations in superconducting microwave cavities coupled to ensembles of Rydberg atoms [@Kaluzny1983] and later to a single atom [@Meschede1985;@Brune1996;@Varcoe2000], vacuum Rabi splitting [@Thompson1992], collapse and revival of oscillations in the atomic population [@Rempe1987], and sub-Poissonian photon statistics [@Rempe1990]. The achievement of strong coupling is crucial to control and manipulate quantum systems and for quantum information processing.

*Ultrastrong-coupling regime*, $g \sim \{\omega_0; \Delta\epsilon \}$:

: When $g$ becomes comparable to the typical energy separation of atom and cavity, several nonperturbative effects come into play [@Kockum2019;@Forn-Diaz2019;@Boite2020]. In contrast to the JC physics, states with different number of excitations starts to be hybridized. Ultrastrong coupling has been first realized outside the cQED domain, in a superconducting flux qubit coupled to a coplanar-waveguide resonator [@Niemczyk2010].

*Deep strong-coupling regime*, $g > \{\omega_0; \Delta\epsilon \}$: 

: In this extreme regime, the light-matter coupling exceeds the subsystems bare energy [@Casanova2010]. Also in this case, superconducting qubits allowed observation of unprecedentedly high values of $g$ [@Forn-Diaz2017;@Yoshihara2017]. Currently, the highest ratio of $g/\omega_0 \approx 1.43$ has been measured in metamaterials coupled to cyclotron resonances [@Bayer2017].






## Light-matter coupling in hybrid systems: beyond the atomic QED architecture

<!-- * General overview of hybrid systems [@Xiang2013;@Treutlein2014;@Kurizki2015] -->
In the last decades, the technological advance in fabrication and control of solid-state devices has paved the way for the realization of the QED protocol in several different platforms. In 1992, Weisbuch and co-workers first observed vacuum Rabi splitting in GaAlAs/AlAs quantum wells coupled to an optical microcavity [@Weisbuch1992]. More than a decade later, strong coupling was observed between a single quantum dot and a photonic crystal [@Yoshie2004]. Gradually, the use of semiconducting artificial atoms in connection with microwave cavities and normal-metal or superconducting contacts became a well-established theoretical and experimental route to design mesoscopic QED systems [@Childress2004;@Cottet2015;@Viennot2016;@Cottet2017;@Burkard2020]. The works pertaining to this Thesis fall primarily into this domain. Simultaneously to quantum-dot-based devices, the increasing importance of superconducing qubits resulted in the parallel development of the field of circuit QED [@Blais2004;@Wallraff2004;@Paik2011;@Blais2020]. 
Fostered by outstanding experimental outcomes, these two main classes of *hybrid* systems stimulated theorists to apply this idea in more exotic contexts involving, e.g., topological superconductors [@Dmytruk2015] and topological waveguides [@Kim2020]. Finally, spin ensembles [@Imamoglu2009;@Kubo2010] and magnons [@Goryachev2014;@Zhang2015] coupled to microwave cavities are building blocks in *cavity magnomechanics* [@Zhang2016] and *cavity optomagnonics* [@Osada2016]. 

Each component of an hybrid system offers diverse intrinsic properties and features, which can then be combined in several ways in order to tailor the energy exchange mechanisms and let specific coherent interactions emerge. On the other hand, the same Hamiltonian can be implemented in radically different experimental contexts to explore broader regions of the parameter space. Finally, hybrid platforms allow to design, to some accuracy, the environmental surroundings of a system, uncovering novel phenomena arising from the complex open-system dynamics.

To better understand the working principles and the physics underlying a typical
solid-state hybrid circuit, I quickly
review below their main elements, and later discuss examples of how they can be coupled to realize light-matter coupling.

### Quantum dots

Quantum dots (QDs) are artificial structures with sizes
ranging from several nanometers to a few micrometers that can be realized on
many platforms [@Loss1998;@Kouwenhoven2001;@vanderWiel2002;@Hanson2007]. The goal is to trap a
small integer number of electrons in a potential well, by controlling an
electrostatic potential landscape using a set of gated electrodes. The confining
potential gives rise to discrete and tunable electronic energy levels, thus realizing an
artificial atom. By contacting the QDs with metallic electrodes made
out of normal metals, ferromagnets, or superconductors, it is possible to set up quantum transport experiments in a variety of configurations. Quantum dots can be fabricated in GaAs/AlGaAs heterostructures, in InSb, InAs, Si/SiGe semiconducting nanowires, carbon nanotubes [@Laird2015], Si [@Mi2018a] and Ge [@Hendrickx2020], and graphene [@Molitor2010;@Banszerus2018;@Kurzmann2019;@Kurzmann2019a;@Banszerus2020].
<!-- - All the applications (quantum computing, automated fine-tuning of dots using
  ML [@vanDiepen2018;@Mills2019a;@Lennon2019;@Zwolak2020;@Moon2020]); -->
Furthermore, quantum-dot devices have received great attention lately because of their potential applications as quantum computers [@Loss1998;@Chatterjee2020].

### Superconducting devices

#### Josephson junctions {.unnumbered}

The Josephson effect is arguably the most important consequence of
superconductivity [@Josephson1962;@Josephson1965;@Josephson1974].
It predicts an equilibrium Cooper-pair current flowing between two
superconductors when they are separated by a weak link, with no bias voltage
applied (dc-Josephson effect). The current is related to the phase difference
$\phi$ of the superconducting condensates, and is given by
\begin{equation}
    \label{eq:intro:josephson-dc-current}
    I  = I_c \sin \phi,
\end{equation}
where the critical current $I_c$ corresponds to the maximum current. A constant
bias voltage $V$ applied across the junction leads to the ac-Josephson effect:
The phase difference will evolve in time according to
\begin{equation}
    \label{eq:intro:josephson-ac-current}
    \dot{\phi} = -\frac{2eV}{\hbar}, 
\end{equation}
which leads to an ac-current of amplitude $I_c$ and Josephson frequency $\omega_J =
2eV/\hbar$.

From Eqs. (-@eq:intro:josephson-dc-current) and (-@eq:intro:josephson-ac-current), the electrical work done by a current source to change the phase difference $\Delta \phi$ is given by $I(t) V = \frac{\hbar}{2e} I(t) \dot{(\Delta \phi)}$. Integration in time yields the free energy

$$
F = \mathrm{const.} - \frac{\hbar I_c}{2e} \cos \phi,
$$
{#eq:intro:free-energy-jj}

which is minimum for $\Delta \phi = 0$. The Josephson effect can be formally derived using the Ginzburg-Landau (GL) theory, by considering two bulk superconductors separated by a short link of length $L \ll \xi$, where $\xi$ is the coherence length in the superconductors [@deGennes1999]. Solving the GL equation leads to the system free energy

$$
F = \frac{\hbar I_c}{2e} (1 - \cos \phi),
$$

In agreement with Eq. (-@eq:intro:free-energy-jj). The critical current $I_c$ is proportional to $A/L$, where $A$ is the cross-section surface area of the junction.

To study the dynamics of a Josephson junction in presence of a bias current and finite voltage, it is useful to resort to RCSJ model, in which one considers a circuit made up by an ideal Josephson junction shunted in parallel by a resistance $R$ and a capacitance $C$ [@Tinkham2004]. Given a bias current $I$, the phase difference across the junction obeys the differential equation

$$
\frac{d^2\phi}{d\tau^2} + \frac{1}{Q} \frac{d \phi}{d \tau} + \sin \Delta \phi = \frac{I}{I_c}.
$$
{#eq:intro:rcsj-equation}

I have introduced the dimensionless variable $\tau = \omega_p t$, with the plasma frequency and quality factor of the junction defined, respectively, as

$$
\omega_p = \sqrt{\frac{2 e I_c}{\hbar C}}, \qquad Q = \omega_p RC.
$$


In full analogy with a mechanical system, Eq. (-@eq:intro:rcsj-equation) is equivalent to the equation of motion of a particle of *mass* $(\hbar/2e)^2 C$, moving along the *coordinate* $\phi$, in a *potential*

$$
U(\phi) = -E_J \cos \phi - \frac{\hbar I}{2e}\phi,
$$
{#eq:intro:tilted-washboard-potential}

with a velocity-dependent *damping force* given by $(\hbar/2e)^2(1/R)d\phi/d\tau$. The quantity $E_J = \hbar I_c / 2e$ is known as *Josephson energy*. The potential (-@eq:intro:tilted-washboard-potential) is usually referred to as *tilted washboard potential*.

While the above analysis assumes that the phase difference, $\phi$, and the charge $Q = CV$ on the junction are classical variables, the variances $\Delta\phi$ and $\Delta N$ are related by the quantum uncertainty relation $\Delta \phi \Delta N \gtrsim 1$, where $N = Q/2e$ is the number of Cooper pairs transferred across the junction. This implies that $\phi$ and $N$ cannot be know with arbitrary precision, giving rise to quantum fluctuations both in the phase and charge. By promoting $\phi$ and $N$ to operators, one can write the Hamiltonian of an isolated Josephson junction at $T=0$, without bias current, as

$$
H = -E_J \cos \phi + \frac{Q^2}{2C} = -E_J \cos \phi - 4 E_C \frac{\partial^2}{\partial \phi^2},
$$
{#eq:intro:josephson-hamiltonian}

where I have made the canonical replacement $N \rightarrow i \partial/\partial\phi$, and defined the charging energy of the junction, $E_C = e^2/2C$. The ratio $4 E_C / E_J$ controls the interplay between charge and phase fluctuations. For a small junction, with $C \approx \SI{1}{\femto\farad}$, the junction behaves as a particle with very small *mass*. This implies $E_C \gg E_J$, and that fluctuations in phase $\phi$ are large. The ground-state wavefunction $\psi(\phi)$  of the junction approaches a constant value. Conversely, for large junctions, $E_C \ll E_J$, meaning that the wavefunction is peaked around the minimum of the potential $\phi = 0$, similarly to a heavy particle in a potential well. In this case, charge fluctuations dominate. Notice that, in this regime, the junction behaves as a nonlinear (anharmonic) oscillator, with the nonlinearity given by the cosine form of the potential. The intrinsic nonlinearity of a Josephson junction is key for the implementation of superconducting qubits, which I briefly review below.

<!-- When charge fluctuations are negligible (i.e., for small junctions), it is convenient to rewrite Hamiltonian (-@eq:intro:josephson-hamiltonian) in the basis of states with definite number of Cooper pairs, $|N \rangle$. -->

<!-- #### Andreev reflection {.unnumbered}


#### Andreev bound states {.unnumbered} -->


#### Superconducting qubits {.unnumbered}

Superconducting (SC) qubits are macroscopic circuit elements that behave quantum mechanically and exhibit observable coherence [@Wendin2017;@Gu2017;@Kjaergaard2020]. Typically, SC qubits are anharmonic oscillators, where the separation of energy levels is made nonuniform by introducing a nonlinearity through Josephson junctions. The anharmonicity is crucial to encode a qubit and perform gate operations involving only two states. A central parameter in the design of a SC qubit is the ratio between the charging energy $E_C$ and the Josephson energy $E_J$. According to this ratio and to the circuit topology, SC qubits are usually divided in three main classes: charge, flux (or persistent-current), and phase qubits.

The charge qubit, or Cooper-pair box, was the first superconducting circuit in which temporal coherence was first observed [@Nakamura1999]. It consists of a small superconducting island connected to a large superconducting reservoir via a Josephson junction, and it is characterized by $E_C / E_J \geq 1$, such that the charge on the island is a good quantum number. They can have large anharmonicity, but are generally affected by strong charge noise, which limits their coherence time. To circumvent this problem, a large shunt capacitor can be added to the circuit, realizing the *transmon* [@Koch2007], which is characterized by $E_J/ E_C \gtrsim 50$ and is insensitive to charge noise. Transmons are currently the most widely used qubits for high-fidelity quantum gates, and have been employed in the first experimental demonstration of quantum supremacy [@Arute2019].

In flux qubits [@Mooij1999;@vanderWal2000;@Friedman2000;@Chiorescu2003], a superconducting loop threaded by an external magnetic flux $\Phi_\mathrm{ext}$ is interrupted by one, three or more Josephson junctions. They are characterized by $E_J/E_C \sim 50$ and are generally affected by flux noise, since they need a large self-inductance. Flux qubits are the predominant elements in quantum annealing setups, such as the commercial system D-Wave [@Johnson2011].

Phase-qubit circuits consists of large Josephson junctions in the regime $E_J/E_C \sim \num{e6}$ [@Martinis1985;@Clarke1988;@Yu2002;@Martinis2009]. These circuits can be described by the tilted-washboard potential (-@eq:intro:tilted-washboard-potential), which is controlled by an applied bias current $I$ very close to the critical current $I_c$ of the junction. Phase qubits have small anharmonicity, but are insensitive to offset-charge noise because of the large $E_J/E_C$ ratio.

Many extensions and optimizations of the above qubits have been devised in order to mitigate detrimental effects. The first improvement to the original Cooper-pair box design led to the *quantronium* [@Cottet2002] with $E_J/E_C \sim 1$ (also called *charge-flux qubit*). To reduce the flux sensitivity in flux qubits, the loop size can be reduced, introducing however charge noise. For this reason, the small junction of a three-junction flux qubit can be shunted by a large capacitance [@You2007]. The shunt-capacitor design has been implemented also for the phase qubit [@Steffen2006]. Other variations of the three basis SC qubits include *four-junction flux qubit* [@Bertet2005],  *tunable-gap flux qubit* [@Paauw2009], *fluxonium* [@Manucharyan2009], *xmon*, *gmon* and *gatemon* qubits [@Barends2013;@Chen2014a;@Larsen2015].
Finally, a recent realization of a $0-\pi$ qubit holds the promise for the implementation of fault-tolerant quantum processors [@Gyenis2019].

### Electromagnetic cavities and mechanical resonators

The other basic ingredient of a mesoscopic QED device is the cavity, which is
described by a quantized bosonic field that can exchange energy with the
artificial atoms. Typically, a single resonant mode of given frequency is
considered for each cavity. Any real cavity has a finite quality factor $Q$, which is
connected to its energy decay rate $\kappa$ through $Q = \omega/\kappa$
($\omega$ is the resonance frequency of the mode). 

Optical cavities have first been used in QED experiments and consist of
conventional Fabry-Perot cavities with
large $Q$ up to \num{e8} [@Hood2001]. They provide ideal coupling to atoms and
spins, through dipole and magnetic interaction, respectively. Beside
the Fabry-Perot cavities, other resonators working in the optical regime have
been envisaged, such as microsphere [@Vernooy1998], microtoroidal [@Armani2003], and photonic band-gap cavities [@Lev2004].

In the microwave regime, cavities are usually fabricated on-chip in superconducting
circuits. In the coplanar waveguide (transmission line) geometry, the resonator is built out of two ground
conductors on the side of a central superconducting wire, with two capacitors on
the ends playing the role of the mirrors of a conventional Fabry-Perot cavity [@Wallraff2004;@Sillanpaa2007;@Hofheinz2008;@Hofheinz2009].
To achieve microwave frequencies, the size of a waveguide is in the millimeter
range, with quality factors of order $Q \approx \num{e3}-\num{e4}$. Another example of microwave resonator is offered by LC
resonators, composed by an inductor and a capacitor, with resonance frequency
$\omega = 1/\sqrt{LC}$. 

Nanomechanical resonators have attracted considerable interest in the last
decade, because their ability to respond to electrical and magnetic fields can
lead to ultrasensitive force detectors [@Rugar2004;@Chaste2012]. On the other
hand, they can be coupled to two- or few-level systems to investigate, prepare and
detect quantum states of mechanical vibration, or serve as a hybrid quantum
information processing platform. Mechanical resonators can vibrate at
frequencies ranging from tens of MHz to a few GHz [@Poot2012] with ultrahigh quality
factors [@Moser2014]. Because of the low frequency, the thermal energy $k_B T$ is
usually larger than the vibrational energy $\hbar \omega$. Hence, resonators need to be
actively cooled in order to exhibit quantum behavior. To this end, great
theoretical [@Wilson-Rae2004;@Blencowe2005;@Liu2013;@Ojanen2014;@Stadler2014;@Stadler2016] and experimental [@Naik2006;@OConnell2010;@Teufel2011;@Chan2011;@Peterson2016;@Clark2017] effort has been devoted to the development of
cooling setups. Diverse types of mechanical resonators include cantilevers, doubly-clamped beams, carbon nanotubes, nanodrums and membranes [@Poot2012].

### Engineering the interaction

The building blocks described above can be coupled together in order to realize hybrid systems that enhance the properties and the advantages of its constituents. Utilizing the charge degree of freedom, quantum dots can be dipole-coupled to microwave cavities [@Frey2012;@Basset2013] achieving stroung coupling of single electrons to photons [@Mi2017;@Li2018;@Wang2019;@Wang2020]. The coupling can also be realized using mechanical resonators such as carbon nanotubes [@Steele2009;@Okazaki2016]. 

Other efforts have been focused on coupling electronic *spin* to mechanical resonators [@Rugar2004;@Palyi2012;@Stadler2014], or to realize an effective spin-photon coupling based on  spin-charge and charge-photon coupling with quantum dots and microwave cavities [@Viennot2015;@Mi2018;@Samkharadze2018;@Landig2018;@Cubaynes2019;@Pan2020].

Beyond the typical *light-matter* paradigm, other hybrid implementations have been oriented to couple different oscillators: The optomechanical setup aims at achieving stroung coupling between mechanical and optical or microwave resonators exploiting the radiation pressure mechanism [@Teufel2011a;@Aspelmeyer2014], while recently, indirect coupling of macroscopic mechanical resonators hase been demonstrated [@Deng2016;@Luo2018;@Zhang2020].

The mesoscopic QED architecture with quantum dots is not simply an alternative tool to reproduce atomic or circuit QED experiments. 
First, quantum
transport is an important ingredient as it can drive the system out of
equilibrium. Hence, microwave photons can help to probe and detect these
nonequilibrium states [@Viennot2014;@Stockklauser2015;@Bruhat2016]; on the
other hand, microwave or mechanical cavities can themselves exploit
the out-of-equilibrium charge dynamics and become highly excited
[@Liu2014;@Liu2015;@Liu2017;@Liu2017a;@Urgell2019;@Wen2019;@Willick2020]. The present Thesis
will predominantly focus on this aspect.
Secondly, the interface between the conductors and the quantum dots can be made
highly transparent, therefore increasing the tunnel coupling strength. For
superconducting contacts, this can lead to coherent splitting and recombination
of Cooper pairs [@Hofstetter2009;@Herrmann2010]. In other setups, strong induced correlation can
enable detection of exotic states of matter, such as Kondo resonances [@Desjardins2017;@Borzenets2020] or Majorana bound
states [@Dartiailh2017;@Trif2019].

<!-- * Relevant examples that anticipate the results of the thesis. Especially: single atom laser to be studied again as example in the quantum master equation chapter. -->



## Structure of this Thesis and goals

The main goal of this Thesis is to discuss and investigate a
family of hybrid setups where quantum dots (or Josephson junctions) are coupled to bosonic cavities (nanomechanical
resonators or superconducting microwave cavities), and drive the latter into highly nonequilibrium states through electron tunneling events. I outline below the content of each Chapter:

* In Ch. -@sec:theory I introduce the theoretical and numerical framework
  upon which the main results are obtained, i.e., the Markovian quantum master
  equation to study the evolution of typical mesoscopic QED architectures.
* In Ch. -@sec:qdlaser I discuss a solid-state implementation of a
  single-atom maser, which is realized using a quantum dot in a
  spin-valve configuration between ferromagnetic contacts; I show how the
  systems exhibits unique features, such as the breakdown of the rotating-wave
  approximation, widely used in the light-matter interaction context, and a
  multistable lasing regime.
* In Ch. -@sec:cps I present a single-photon pump device based on a
  quantum-dot implementation of a Cooper-pair splitter, where an effective
  coupling between two distant
  harmonic resonators can be activated through Cooper-pair transport from a
  superconducting contact. This enables an efficient photon transfer mechanism, as well
  as cavity ground-state cooling.
* In Ch. -@sec:jjcavity I examine the theory of double-Cooper-pair
  tunneling in a voltage-biased Josephson junction coupled to a microwave
  resonator. I analyze the transition regime between coherent and incoherent
  double-Cooper-pair tunneling, and further show how the system can be used as a
  single-photon source.
* In Ch. -@sec:conclusions I draw the conclusions and give an outlook on
  the work presented.

