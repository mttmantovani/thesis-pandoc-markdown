# Introduction {#sec:introduction}


## Coupling light and matter: Cavity quantum electrodynamics

One of the most succesful and ever-developing applications of the quantum theory
is the investigation of the interaction between *light* and *matter*.
Historically, the first experiments involved real atoms
interacting with the quantized electromagnetic field of optical and microwave
cavities, within the framework commonly
referred to as cavity quantum electrodynamics (cQED)
[@Raimond2001;@Walther2006;@Haroche2006].
 
The paradigmatic model for the light-matter interaction is represented by the Rabi Hamiltonian [@Rabi1936;@Rabi1937], given by

$$
  H_\text{Rabi} = \frac{\Delta\epsilon}{2}\sigma_z + \omega_0 b^\dagger b + g
  (\sigma_+ + \sigma_-) (b + b^\dagger).
$$
{#eq:intro:rabi-hamiltonian}

It describes the coherent interaction between a two-level system (the *matter*, e.g., two levels of an atom or a spin-1/2) and a quantum harmonic oscillator (the *light*, e.g., a single mode of a microwave or optical cavity). The Rabi Hamiltonian follows from the interaction between the dipole moment of an atom and the electric field of the cavity, in the dipole approximation [@Scully1997]. In Eq. (-@eq:intro:rabi-hamiltonian), $\Delta \epsilon$ is the energy difference between the two levels and $\omega_0$ is the resonance frequency of the harmonic oscillator. The constant $g$ quantifies the strength of the interaction between light and matter.

![(a) Sketch of a cavity QED setup. A two-level atom with energy separation $\Delta \epsilon$ is coupled with strength $g$ to a cavity with strength resonance frequency $\omega_0$. The atom ](figures/intro-cqed.pdf){#fig:intro:cqed}

The two-level system is composed of an excited ($|e\rangle$) and a ground ($|g \rangle$) state, which can be represented by two column vectors $|e\rangle = (1, 0)^T$ and $|g\rangle = (0, 1)^T$. In the basis $\{|e\rangle, |g\rangle \}$, the qubit is described by the raising and lowering operators $\sigma_+ = |e \rangle \langle g|$ and $\sigma_- = |g \rangle \langle e |$, which generate the Pauli operators $\sigma_x = \sigma_+ + \sigma_-$, $\sigma_y = -i (\sigma_+ - \sigma_-)$ and $\sigma_z = [\sigma_+, \sigma_-]$, with matrix form

$$
\sigma_x = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}, \quad 
\sigma_y = \begin{pmatrix} 0 & -i \\ i &  0 \end{pmatrix}, \quad
\sigma_z = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}.
$$

The harmonic oscillator is described by the bosonic field annihilation and creation operators, $b$ and $b^\dagger$. 

The exact eigenenergies of Eq. (-@eq:intro:rabi-hamiltonian) have been found only several decades after the first works of Rabi, and can only be written as the roots of a complex trascendental function [@Braak2011]. The primary way of dealing with the Rabi Hamiltonian is to consider a *small* value of the coupling strength, $g \ll \{\omega_0; \Delta\epsilon\}$ and to tune the two systems close to resonance, i.e., $\Delta\epsilon \approx \omega_0$. This leads to the rotating-wave approximation (RWA) and to the Jaynes-Cummings Hamiltonian [@Jaynes1963]:

$$
H_\text{JC} = \frac{\Delta\epsilon}{2}\sigma_z + \omega_0 b^\dagger b + g
   (\sigma_+ b + b^\dagger \sigma_-),
$$
{#eq:intro:jc-hamiltonian}

which can be exactly solved, see Ch. -@sec:qdlaser.

Later, the Rabi model has been extended to two-photon
interactions, where $b$ and $b^\dagger$ are replaced by $b^2$ and
$(b^\dagger)^2$ [@Travenec2012], and to $N$-level atoms [@Albert2012].
 
 * Jaynes-Cummings Hamiltonian [@Jaynes1963], extensions (Dicke [@Dicke1954] and
  Tavis-Cummings [@Tavis1968;@Tavis1969] models)

  

  - Rydberg atoms; micromaser.

When dealing with a cavity QED experiment, an important role is played by the surroundings of the light-matter system: Atoms and cavities are coupled not only to each other--which is desired--but also to the measurement apparati and to the thermal environment. 

* Weak-coupling regime: $g \lesssim \gamma, \kappa$, Purcell effect
  (enhancement of the spontaneous emission rate of an atom due to the coupling
  to the radiation field) [@Purcell1946];
* Strong-coupling regime: $g \gg \gamma, \kappa$, vacuum Rabi oscillations
  [@Brune1996;@Varcoe2000],
  vacuum Rabi splitting [@Thompson1992], collapse and revival of
  oscillations [@Rempe1987], sub-Poissonian
  photon statistics [@Rempe1990];
* Ultrastrong-coupling regime: $g \sim \omega, \Delta\epsilon$
  [@Niemczyk2010;@Kockum2019;@Forn-Diaz2019].
* Deep strong-coupling regime: $g > \omega, \Delta\epsilon$ [@Casanova2010;@Bayer2017;@Yoshihara2017].


## Light-matter coupling in hybrid systems: beyond the QED architecture

In the last decades, the technological advance in fabrication and control of solid-state devices has paved the way for the realization of the QED protocol in several different platforms. 

* Circuit QED [@Blais2004;@Wallraff2004;@Paik2011]
* Mesoscopic QED [@Childress2004;@Cottet2015;@Viennot2016;@Cottet2017;@Burkard2020]
* General overview of hybrid systems [@Xiang2013;@Treutlein2014;@Kurizki2015]

To better understand the working principles and the physics underlying a typical
solid-state hybrid circuit, we
review below their main components, and later discuss how they can be coupled to
yield an *ad-hoc* designed Hamiltonian.

### Quantum dots

Quantum dots (QDs) are artificial structures with sizes
ranging from several nanometers to a few micrometers that can be realized on
many platforms [@Loss1998;@Kouwenhoven2001;@vanderWiel2002;@Hanson2007]. The goal is to trap a
small integer number of electrons in a potential well, by controlling an
electrostatic potential landscape using a set of gated electrodes. The confining
potential gives rise to discrete and tunable electronic energy levels, thus realizing an
artificial atom. By contacting the QDs with metallic electrodes made
out of normal metals, ferromagnets, or superconductors, it is possible to set up
quantum transport experiments in a variety of configurations.



- Some examples: heterostructures, ...
- GaAs/AlGaAs
- InSb, InAs, Si/SiGe nanowires
- Si **[CITE]** and Ge [@Hendrickx2020]
- CNTQDs [@Laird2015]
- Graphene
- All the applications (quantum computing, automated fine-tuning of dots using
  ML [@vanDiepen2018;@Mills2019a;@Lennon2019;@Zwolak2020;@Moon2020]);
- QD arrays [@Mills2019]

### Josephson junctions and superconducting qubits

The Josephson effect is arguably the most important consequence of
superconductivity [@Josephson1962;@Josephson1965;@Josephson1974].
It predicts an equilibrium Cooper-pair current flowing between two
superconductors when they are separated by a weak link, with no bias voltage
applied (dc-Josephson effect). The current is related to the phase difference
$\Delta\phi$ of the superconducting condensates, and is given by
\begin{equation}
    \label{eq:theory:josephson-dc-current}
    I  = I_c \sin \Delta\phi,
\end{equation}
where the critical current $I_c$ corresponds to the maximum current. A constant
bias voltage $V$ applied across the junction leads to the ac-Josephson effect:
The phase difference will evolve in time according to
\begin{equation}
    \label{eq:theory:josephson-ac-current}
    \dot{\Delta\phi} = -\frac{2eV}{\hbar}, 
\end{equation}
which leads to an ac-current of amplitude $I_c$ and Josephson frequency $\omega_J =
2eV/\hbar$.

- Josephson junction
- Andreev reflection
- SC qubits (overview)

	Superconducting (SC) qubits are macroscopic circuit elements that behave quantum mechanically and thus exhibit observable coherence. Typically, SC qubits are anharmonic oscillators, where the separation of energy levels is made nonuniform by introducing a nonlinearity through Josephson junctions. The anharmonicity is crucial to encode a qubit and perform gate operations involving only two states.
	- Reviews [@Wendin2017;@Gu2017];
	- Cooper-pair box [@Nakamura1999];
	- 3-junction flux qubit [@Mooij1999;@vanderWal2000];
	- Transmon [@Koch2007];
	- Fluxonium [@Manucharyan2009];
	- $0-\pi$ qubit [@Gyenis2019].


### Electromagnetic cavities and mechanical resonators

The other basic ingredient of a mesoscopic QED device is the cavity, which is
described by a quantized bosonic field that can exchange energy with the
artificial atoms. Typically, a single resonant mode of given frequency is
considered for each cavity. Any real cavity has a finite quality factor $Q$, which is
connected to its energy decay rate $\kappa$ through $Q = \omega/\kappa$
($\omega$ is the resonance frequency of the mode). 

Optical cavities have first been used in QED experiments and consist of
conventional Fabry-Perot cavities (see Fig. \ref{fig:intro:fp-cavity}) with
large $Q$ up to \num{e8} **[CITE]**. They provide ideal coupling to atoms and
spins, through dipole and magnetic interaction, respectively. Beside
the Fabry-Perot cavities, other resonators working in the optical regime have
been envisaged **[EXAMPLES]**.

In the microwave regime, cavities are usually fabricated on-chip in superconducting
circuits. In the coplanar waveguide (transmission line) geometry pictured in
Fig. \ref{fig:intro:cpw-cavity}, the resonator is built out of two ground
conductors on the side of a central superconducting wire, with two capacitors on
the ends playing the role of the mirrors of a conventional Fabry-Perot cavity.
To achieve microwave frequencies, the size of a waveguide is in the millimeter
range, with quality factors of order $Q \approx \num{e3}-\num{e4}$ **[CHECK,
CITATIONS]**. Another example of microwave resonator is offered by LC
resonators, composed by an inductor and a capacitor, with resonance frequency
$\omega = 1/\sqrt{LC}$. 

Nanomechanical resonators have attracted considerable interest in the last
decade, because their ability to respond to electrical and magnetic fields can
lead to ultrasensitive force detectors [@Rugar2004;@Chaste2012]. On the other
hand, they can be coupled to two- or few-level systems to investigate, prepare and
detect quantum states of mechanical vibration, or serve as a hybrid quantum
information processing platform. Mechanical resonators can vibrate at
frequencies ranging from tens of MHz to a few GHz **[CITE]** with ultrahigh quality
factors [@Moser2014]. Because of the low frequency, the thermal energy $k_B T$ is
usually larger than the vibrational energy $\hbar \omega$. Hence, resonators need to be
actively cooled in order to exhibit quantum behavior. To this end, great
theoretical [@Wilson-Rae2004;@Blencowe2005;@Liu2013;@Ojanen2014;@Stadler2014;@Stadler2016] and experimental [@Naik2006;@OConnell2010;@Teufel2011;@Chan2011;@Peterson2016;@Clark2017] effort has been devoted to the development of
cooling setups.


- Poot [@Poot2012]
- High Quality factors
- Types of resonators: cantilevers, doubly-clamped beams, carbon nanotubes,
  nanodrums and membranes.

### Engineering the interaction Hamiltonian



#### Charge coupling with microwave resonators
* Examples of how to couple the different subsystems;
	- Coupling dots (charge, dipole, Andreev) and SC-qubits to cavities [@Frey2012;@Basset2013;@Janvier2015;@Mi2017;@Li2018;@Wang2019;@Wang2020].
  - Coupling charge to mechanical resonators [@Steele2009];
  - Coupling spins to mechanical resonators [@Rugar2004;@Palyi2012;@Stadler2014];
  - Coupling spins to microwave cavities (Petta/Kontos)
    [@Viennot2015;@Mi2018;@Samkharadze2018;@Landig2018;@Cubaynes2019;@Pan2020];
  - Coupling mechanical resonators to microwave/optical resonators
    (optomechanics) [@Teufel2011a;@Aspelmeyer2014]
  - Coupling mechanical resonators to mechanical resonators (Guo)
    [@Deng2016;@Luo2018;@Zhang2020]

* **What are the advantages of hybrid systems (with QDs)?** 
* The mesoscopic QED architecture with quantum dots, however, is not simply an
alternative tool to reproduce atomic or circuit QED experiments. 
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
of Cooper pairs **[CITE CPS]**. In other setups, strong induced correlation can
enable detection of exotic states of matter, such as Kondo resonances [@Desjardins2017;@Borzenets2020] or Majorana bound
states [@Dartiailh2017;@Trif2019].
* Relevant examples that anticipate the results of the thesis. Especially: single atom laser to be studied again as example in the quantum master equation chapter.
  

**Other references**

Superconducting qubits:

* Review [@Kjaergaard2020]
  
Quantum dots and semiconducting qubits:

* New review [@Chatterjee2020]

Circuit QED:

* New review [@Blais2020]

Hybrid systems:

* Topology [@Kim2020];
  

## Structure of this Thesis and goals

The main goal of this Thesis is to discuss and investigate a
family of hybrid setups where quantum dots (or Josephson junctions) are coupled to bosonic cavities (nanomechanical
resonators or superconducting microwave cavities), and drive the latter into highly
nonequilibrium states through electron tunneling events. I outline below the content of each Chapter:

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

