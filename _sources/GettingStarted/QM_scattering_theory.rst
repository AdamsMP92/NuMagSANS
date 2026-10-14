.. _quantum-scattering-theory:

Quantum Scattering Theory
=========================

.. math::

   \boxed{
   \begin{array}{c}
      \text{stationary Schrödinger equation} \\
      \Downarrow \\
      \text{free incident solution} \\
      \Downarrow \\
      \text{outgoing Green function} \\
      \Downarrow \\
      \text{Lippmann--Schwinger integral equation} \\
      \Downarrow \\
      \text{far-field asymptotics} \\
      \Downarrow \\
      \text{scattering amplitude} \\
      \Downarrow \\
      \text{incident and scattered probability currents} \\
      \Downarrow \\
      \text{differential cross section} \\
      \Downarrow \\
      \text{first Born approximation} \\
      \Downarrow \\
      \text{Fourier transform of the interaction potential} \\
      \Downarrow \\
      \text{nuclear and magnetic scattering operator} \\
      \Downarrow \\
      \text{spin projection and SANS cross sections}
   \end{array}}

This page derives the elastic nuclear and magnetic SANS cross sections used by
NuMagSANS from a single quantum-mechanical scattering problem.  The purpose is
not to derive the strong neutron--nucleus interaction from first principles.
Instead, the coherent nuclear scattering-length density and the sample
magnetization are taken as the effective material fields, and the route from
these fields to the spin-resolved cross sections is made explicit.

Fermi's golden rule provides an equivalent first-order route from the
interaction matrix element to the transition rate.  It is therefore an
alternative to the scattering-state argument at the cross-section step, not
an additional approximation that must be inserted into the chain above.


Scope and assumptions
---------------------

The derivation describes the model implemented by NuMagSANS under the
following assumptions:

- the sample is static and the scattering is elastic, so that
  :math:`|\mathbf{k}_{f}|=|\mathbf{k}_{i}|=k`;
- nuclear scattering is coherent and spin independent;
- the neutron is nonrelativistic and the first Born approximation is valid;
- magnetic scattering is produced by a magnetostatic magnetization field;
- the nuclear and magnetic amplitudes are added coherently before their
  absolute square is taken;
- the incident and analyzed spin states are ideal eigenstates of one chosen
  polarization axis.

NuMagSANS does not include nuclear spin-incoherent scattering, inelastic
sample dynamics, multiple scattering, refraction at interfaces, instrumental
resolution, finite polarization efficiency, or spin leakage.  These effects
require extensions of the theory or experimental corrections.

The small-angle approximation :math:`q_x\simeq0` is introduced only when the
general result is specialized to the NuMagSANS detector geometry.  It is not
needed for the scattering-theory steps themselves.


Notation used below
-------------------

The same symbols occur in ordinary three-dimensional space and in the
two-dimensional neutron spin space.  The following distinctions are useful:

- bold symbols such as :math:`\mathbf r`, :math:`\mathbf B`, and
  :math:`\mathbf M` are vectors in real space;
- a hat on a bold vector, for example
  :math:`\widehat{\mathbf q}=\mathbf q/|\mathbf q|`, denotes a unit vector;
- a hat on a non-bold symbol, for example :math:`\widehat U`, denotes an
  operator acting on the two spin components;
- :math:`\mathbb I_2` and :math:`\mathbb I_3` are the identity matrices in
  two-dimensional spin space and three-dimensional real space, respectively;
- a tilde denotes a spatial Fourier transform, an asterisk denotes complex
  conjugation, and a dagger denotes a complex-conjugate transpose;
- the subscripts :math:`i` and :math:`f` mean ``initial`` and ``final``;
- :math:`d\sigma/d\Omega` is the differential cross section of one scattering
  system, whereas :math:`d\Sigma/d\Omega=(1/V)d\sigma/d\Omega` is the
  macroscopic cross section per sample volume and has units of inverse length.

All real-space nuclear and magnetization fields considered here are real.
Their Fourier transforms can nevertheless be complex because spatial phase
information is encoded in the Fourier amplitudes.


1. Stationary Schrödinger equation and effective neutron Hamiltonian
--------------------------------------------------------------------

Consider a neutron of mass :math:`m_n` in a static sample.  Its stationary
wavefunction is a two-component spinor.  In any chosen orthonormal spin basis
:math:`(\chi_1,\chi_2)`, it can be written explicitly as

.. math::

   \Psi(\mathbf r)
   =\psi_1(\mathbf r)\chi_1+\psi_2(\mathbf r)\chi_2
   =\begin{pmatrix}
      \psi_1(\mathbf r)\\
      \psi_2(\mathbf r)
   \end{pmatrix}
   \in\mathbb C^2.

The complex functions :math:`\psi_1` and :math:`\psi_2` are probability
amplitudes for the two basis states.  Their local probability density is
:math:`\Psi^{\dagger}\Psi=|\psi_1|^2+|\psi_2|^2`.  The effective low-energy
Hamiltonian used throughout this derivation is

.. math::

   \boxed{
   \widehat H
   =-\frac{\hbar^2}{2m_n}\nabla^2\mathbb I_2
   +U_{\mathrm N}(\mathbf r)\mathbb I_2
   -\widehat{\boldsymbol\mu}_n\cdot\mathbf B(\mathbf r).}

Because the sample is static, :math:`\widehat H` has no explicit time
dependence.  The time-dependent Schrödinger equation

.. math::

   i\hbar\frac{\partial}{\partial t}\boldsymbol\psi(\mathbf r,t)
   =\widehat H\boldsymbol\psi(\mathbf r,t)

therefore admits states with a separated time dependence,

.. math::

   \boldsymbol\psi(\mathbf r,t)
   =e^{-iEt/\hbar}\Psi(\mathbf r).

Substitution removes the common time-dependent factor and gives the
stationary Schrödinger equation used for elastic scattering,

.. math::

   \widehat H\Psi(\mathbf r)=E\Psi(\mathbf r).

The kinetic-energy operator acts in the same way on both components.  The
interaction :math:`\widehat U(\mathbf r)`, in contrast, is a local
:math:`2\times2` matrix and may mix them.  At the level required here it
contains a scalar nuclear part and a magnetic Zeeman part,

.. math::

   \widehat U(\mathbf r)
   =U_{\mathrm N}(\mathbf r)\mathbb I_2
   -\widehat{\boldsymbol\mu}_n\cdot\mathbf B(\mathbf r).

The symbol :math:`\widehat{\boldsymbol\mu}_n` is a *vector of three spin
operators*, not an ordinary spatial vector.  Introduce the neutron angular
momentum operator

.. math::

   \widehat{\mathbf S}_n
   =\frac{\hbar}{2}\boldsymbol\sigma^{(S)},

where :math:`\boldsymbol\sigma^{(S)}` is the vector of Pauli matrices in an
angular-momentum basis.  The physical magnetic-moment operator is

.. math::

   \boxed{
   \widehat{\boldsymbol\mu}_n
   =\gamma_n\widehat{\mathbf S}_n
   =\mu_n\boldsymbol\sigma^{(S)},}
   \qquad
   \mu_n=\frac{\gamma_n\hbar}{2}<0.

The negative sign expresses that the neutron magnetic moment is antiparallel
to its spin angular momentum.  The dot product in the interaction means

.. math::

   \widehat{\boldsymbol\mu}_n\cdot\mathbf B
   =\widehat\mu_{n,x}B_x
   +\widehat\mu_{n,y}B_y
   +\widehat\mu_{n,z}B_z.

Each :math:`\widehat\mu_{n,j}` is a :math:`2\times2` matrix, while each
:math:`B_j(\mathbf r)` is a real-valued classical field component.  A magnetic
field not aligned with the selected spin basis produces off-diagonal matrix
elements and can therefore flip the neutron channel state.

The coherent nuclear interaction is represented by the Fermi
pseudopotential,

.. math::

   U_{\mathrm N}(\mathbf r)
   =\frac{2\pi\hbar^2}{m_n}N(\mathbf r),

where :math:`N(\mathbf r)` is the coherent nuclear scattering-length density.
It has units of inverse area.  For discrete nuclei at positions
:math:`\mathbf R_j`, the idealized field is

.. math::

   N(\mathbf r)=\sum_j b_j\delta(\mathbf r-\mathbf R_j),

where :math:`b_j` is the coherent nuclear scattering length.  A continuum
description replaces this delta-function field by the corresponding local or
coarse-grained scattering-length density.

This is already an effective low-energy description: the short-ranged strong
interaction with the nuclei has been replaced by its coherent scattering
length.

The scalar nuclear potential multiplies :math:`\mathbb I_2` and cannot change
the neutron channel in this model.  The magnetic term couples the physical
moment operator to the magnetic induction :math:`\mathbf B` generated by the
sample.  Its connection to :math:`\mathbf M` is stated when the NuMagSANS
scattering operator is introduced and derived from Maxwell's equations in
:ref:`transverse-magnetization-derivation`.


2. From the Schrödinger equation to the Lippmann--Schwinger equation
---------------------------------------------------------------------

The purpose of the Lippmann--Schwinger equation is to rewrite the differential
equation as ``known incident wave + wave generated by the sample``.  This is
still an exact reformulation.  Its advantage is that the effect of the sample
appears as an integral correction that can later be approximated
systematically.


The free incoming solution
^^^^^^^^^^^^^^^^^^^^^^^^^^

First separate free propagation from the sample interaction,

.. math::

   \widehat H_0=-\frac{\hbar^2}{2m_n}\nabla^2\mathbb I_2,
   \qquad
   \widehat H=\widehat H_0+\widehat U.

For a free neutron, :math:`E=\hbar^2k^2/(2m_n)`.  The symbol :math:`k` is the
magnitude of the elastic incident and final wavevectors.  With
:math:`\widehat U=0`, the stationary equation becomes

.. math::

   (E-\widehat H_0)\Psi_{\mathrm{in}}(\mathbf r)=0.

This is the homogeneous free-particle equation.  The experimental preparation
selects one of its solutions as the incident neutron beam.  For an incident
wavevector :math:`\mathbf k_i` and a normalized two-component spinor
:math:`\chi_i`, choose

.. math::

   \boxed{
   \Psi_{\mathrm{in}}(\mathbf r)
   =e^{i\mathbf k_i\cdot\mathbf r}\chi_i,}
   \qquad
   |\mathbf k_i|=k.

Thus :math:`\Psi_{\mathrm{in}}` is not a second unknown wavefunction.  It is
the known free wave sent toward the sample.


The sample interaction as a source
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Restoring :math:`\widehat U` and moving the free part to the left gives

.. math::

   \boxed{
   (E-\widehat H_0)\Psi^{(+)}(\mathbf r)
   =\widehat U(\mathbf r)\Psi^{(+)}(\mathbf r).}

The right-hand side can be read as a spatially distributed source.  At every
point :math:`\mathbf r` in the sample, the interaction acts on the neutron
wave already present there and generates a contribution to the scattered
wave.  This source is self-consistent because it contains the still unknown
complete wave :math:`\Psi^{(+)}`.


The outgoing free Green function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A Green function describes how a point source propagates under a given
differential operator.  For the free neutron, define the outgoing Green
function by

.. math::

   (E-\widehat H_0)
   \left[
      G_0^{(+)}(\mathbf r-\mathbf r';E)\mathbb I_2
   \right]
   =\delta(\mathbf r-\mathbf r')\mathbb I_2.

The identity matrix appears because free propagation acts identically on the
two spin components.  Since
:math:`E=\hbar^2k^2/(2m_n)`, the required spatial kernel follows from the
outgoing Green function of the Helmholtz operator,

.. math::

   (\nabla^2+k^2)
   \left[-\frac{e^{ikR}}{4\pi R}\right]
   =\delta(\mathbf R),
   \qquad
   R=|\mathbf r-\mathbf r'|.

Consequently,

.. math::

   \boxed{
   G_0^{(+)}(\mathbf r-\mathbf r';E)
   =-\frac{m_n}{2\pi\hbar^2}
     \frac{e^{ik|\mathbf r-\mathbf r'|}}
          {|\mathbf r-\mathbf r'|}.}

The factor :math:`e^{ikR}/R` is a spherical wave travelling away from its
source.  The superscript :math:`(+)` records this outgoing boundary condition.
In the equivalent resolvent notation used in scattering theory, the same
choice is encoded by replacing :math:`E` with :math:`E+i0^+`.


Construction of the integral equation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Applying the Green function to the distributed source gives a particular
solution of the inhomogeneous equation,

.. math::

   \Psi_{\mathrm{sc}}(\mathbf r)
   =\int
     G_0^{(+)}(\mathbf r-\mathbf r';E)
     \widehat U(\mathbf r')
     \Psi^{(+)}(\mathbf r')\,d^3r'.

The complete solution is the sum of the chosen solution of the homogeneous
equation and this sample-generated particular solution,

.. math::

   \Psi^{(+)}(\mathbf r)
   =\Psi_{\mathrm{in}}(\mathbf r)+\Psi_{\mathrm{sc}}(\mathbf r).

Inserting the explicit Green function produces the outgoing
Lippmann--Schwinger equation in position space,

.. math::

   \boxed{
   \Psi^{(+)}(\mathbf r)
   =\Psi_{\mathrm{in}}(\mathbf r)
   -\frac{m_n}{2\pi\hbar^2}
   \int
   \frac{e^{ik|\mathbf r-\mathbf r'|}}
        {|\mathbf r-\mathbf r'|}
   \widehat U(\mathbf r')\Psi^{(+)}(\mathbf r')\,d^3r'.}

The equivalence with the differential equation can be checked directly.
Applying :math:`E-\widehat H_0` to both sides gives

.. math::

   \begin{aligned}
   (E-\widehat H_0)\Psi^{(+)}(\mathbf r)
   &=(E-\widehat H_0)\Psi_{\mathrm{in}}(\mathbf r) \\
   &\quad+\int
      (E-\widehat H_0)
      \left[G_0^{(+)}(\mathbf r-\mathbf r';E)\mathbb I_2\right]
      \widehat U(\mathbf r')\Psi^{(+)}(\mathbf r')\,d^3r' \\
   &=0+\int\delta(\mathbf r-\mathbf r')
      \widehat U(\mathbf r')\Psi^{(+)}(\mathbf r')\,d^3r' \\
   &=\widehat U(\mathbf r)\Psi^{(+)}(\mathbf r).
   \end{aligned}

Because the free kernel depends only on :math:`\mathbf r-\mathbf r'`, the
integral has the form of a convolution of the Green function with the source
:math:`\widehat U\Psi^{(+)}`.  The equation is nevertheless implicit: the
unknown complete wave occurs on both sides.  It is exactly equivalent to the
stationary Schrödinger equation together with the specified incident wave and
the outgoing boundary condition.  No Born approximation has yet been made.


3. From the scattering state to the differential cross section
----------------------------------------------------------------

Solving the Schrödinger or Lippmann--Schwinger equation gives a wavefunction,
whereas an experiment records a rate of neutrons arriving at a detector.  A
cross section is the quantity that connects these two descriptions.  This
connection requires more than the local probability density
:math:`\Psi^\dagger\Psi`: a density states how much probability is present
near a point, but not how rapidly or in which direction it moves through a
detector surface.  The required quantity is the probability current.


Why the probability current appears
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The probability current is fixed by probability conservation rather than
introduced as an additional scattering assumption.  Return temporarily to
the time-dependent spinor :math:`\boldsymbol\psi(\mathbf r,t)` and define its
probability density by

.. math::

   \rho(\mathbf r,t)
   =\boldsymbol\psi^\dagger(\mathbf r,t)
    \boldsymbol\psi(\mathbf r,t).

Use the time-dependent Schrödinger equation and its Hermitian conjugate to
differentiate :math:`\rho`.  For the local Hermitian interaction
:math:`\widehat U=\widehat U^\dagger` used here, the two potential terms
cancel.  The kinetic terms give

.. math::

   \begin{aligned}
   \frac{\partial\rho}{\partial t}
   &=\frac{1}{i\hbar}
     \left[
       \boldsymbol\psi^\dagger\widehat H\boldsymbol\psi
       -(\widehat H\boldsymbol\psi)^\dagger
        \boldsymbol\psi
     \right] \\
   &=\frac{i\hbar}{2m_n}
     \left[
       \boldsymbol\psi^\dagger\nabla^2\boldsymbol\psi
       -(\nabla^2\boldsymbol\psi^\dagger)
        \boldsymbol\psi
     \right].
   \end{aligned}

The vector-calculus identity

.. math::

   \boldsymbol\psi^\dagger\nabla^2\boldsymbol\psi
   -(\nabla^2\boldsymbol\psi^\dagger)\boldsymbol\psi
   =\nabla\cdot
    \left[
       \boldsymbol\psi^\dagger\nabla\boldsymbol\psi
       -(\nabla\boldsymbol\psi^\dagger)\boldsymbol\psi
    \right]

therefore turns the Schrödinger equation into the continuity equation

.. math::

   \boxed{
   \frac{\partial\rho}{\partial t}
   +\nabla\cdot\mathbf j=0,}

with

.. math::

   \boxed{
   \mathbf j[\boldsymbol\psi]
   =\frac{\hbar}{2m_n i}
    \left[
       \boldsymbol\psi^\dagger\nabla\boldsymbol\psi
       -(\nabla\boldsymbol\psi^\dagger)\boldsymbol\psi
    \right].}

Integrating the continuity equation over a volume and applying the divergence
theorem shows the meaning of :math:`\mathbf j`: its outward flux through the
boundary equals the rate at which probability leaves that volume.  Thus
:math:`\mathbf j\cdot\mathbf n\,dA` is the probability per unit time crossing
an oriented surface element :math:`dA` for one neutron state.

This probability current is not yet the experimental neutron flux.  The
latter is a *particle-number current*: the expected number of neutrons crossing
unit area per unit time.  For :math:`\mathcal N` statistically independent
realizations of the same one-neutron preparation, expectation values add, so
that

.. math::

   n_{\mathrm n}=\mathcal N\rho,
   \qquad
   \mathbf J_{\mathrm n}=\mathcal N\mathbf j.

Here :math:`n_{\mathrm n}` is the expected neutron number density and
:math:`\mathbf J_{\mathrm n}` is the neutron-number current density measured
in neutrons per unit area and unit time.  Its magnitude is what is commonly
called the neutron flux.  A detector count rate is the flux of
:math:`\mathbf J_{\mathrm n}` through its active area.  The ensemble may
represent distinct dilute wave packets present in the beam or repeated
one-neutron trials accumulated during a measurement; it does not imply
multiple occupation of one fermionic mode.

For example, a single box-normalized plane wave in volume
:math:`\mathcal V` has

.. math::

   \rho=\frac{1}{\mathcal V},
   \qquad
   \mathbf j=\frac{\mathbf v}{\mathcal V},
   \qquad
   \mathbf v=\frac{\hbar\mathbf k}{m_n}.

If the beam ensemble contains a mean number :math:`\mathcal N` of such
independently prepared neutrons distributed through that volume, its
ensemble-averaged quantities are

.. math::

   n_{\mathrm n}=\frac{\mathcal N}{\mathcal V},
   \qquad
   \mathbf J_{\mathrm n}=n_{\mathrm n}\mathbf v,

which is the familiar classical expression ``particle density times
velocity``.  Quantum mechanics supplies the corresponding probability flow
for each neutron; the beam ensemble converts it into the experimentally
counted particle flow.

The plane waves used below are not normalizable to one over all space.  They
are ideal stationary-beam states, which can equivalently be understood as the
large-volume or broad-wave-packet limit.  This causes no ambiguity in a cross
section.  The same ensemble factor :math:`\mathcal N`, or equivalently the
same overall beam normalization, multiplies both the incident and scattered
number currents.  Hence

.. math::

   \frac{J_{\mathrm n,sc}}{J_{\mathrm n,in}}
   =\frac{\mathcal N j_{\mathrm{sc}}}
          {\mathcal N j_{\mathrm{in}}}
   =\frac{j_{\mathrm{sc}}}{j_{\mathrm{in}}}.

This cancellation is why the cross section can be derived from
single-particle probability currents even though the laboratory measures a
many-neutron count rate.


What the differential cross section measures
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Consider a detector element far from the sample.  If it subtends a solid
angle :math:`d\Omega` at distance :math:`r`, its area perpendicular to the
radial direction is

.. math::

   dA=r^2d\Omega.

Let :math:`j_{\mathrm{in}}` be the magnitude of the incident probability
current density and let :math:`j_{\mathrm{sc},i\rightarrow f}^{(r)}` be the
outward radial current density scattered from the prepared spin state
:math:`\chi_i` into the analyzed spin state :math:`\chi_f`.  The differential
cross section is *defined* by

.. math::

   \boxed{
   \frac{d\sigma_{i\rightarrow f}}{d\Omega}
   =\frac{r^2j_{\mathrm{sc},i\rightarrow f}^{(r)}}
          {j_{\mathrm{in}}}.}

The numerator is the scattered probability per unit time and unit solid
angle; division by the incident probability per unit area and unit time
leaves an area.  The cross section is therefore an effective area that
quantifies how strongly one scattering system redirects the incident beam
into the selected angular and spin channel.  To evaluate this definition, the
incident and scattered parts of the Lippmann--Schwinger solution must now be
examined in the far field.


The outgoing wave in the far field
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The separation into an incident and a scattered wave was already obtained
when constructing the Lippmann--Schwinger equation,

.. math::

   \Psi^{(+)}=\Psi_{\mathrm{in}}+\Psi_{\mathrm{sc}}.

It is not an arbitrary decomposition at this stage.  The experimental
preparation fixes :math:`\Psi_{\mathrm{in}}`, and the outgoing Green function
fixes :math:`\Psi_{\mathrm{sc}}` as the wave generated by the sample.  Far
from a localized sample, :math:`r\gg|\mathbf r'|`, one has

.. math::

   |\mathbf r-\mathbf r'|
   =r-\widehat{\mathbf r}\cdot\mathbf r'
    +O(|\mathbf r'|^2/r),

where :math:`r=|\mathbf r|` is the detector distance,
:math:`\widehat{\mathbf r}=\mathbf r/r` is the observation direction, and
:math:`\mathbf r'` is a point inside the sample.  Consequently,

.. math::

   \frac{e^{ik|\mathbf r-\mathbf r'|}}
        {|\mathbf r-\mathbf r'|}
   \sim
   \frac{e^{ikr}}{r}e^{-i\mathbf k_f\cdot\mathbf r'},
   \qquad
   \mathbf k_f=k\widehat{\mathbf r}.

Inserting this approximation into the sample-generated part of the
Lippmann--Schwinger equation gives

.. math::

   \Psi_{\mathrm{sc}}(\mathbf r)
   \sim\frac{e^{ikr}}{r}
   \left[
      -\frac{m_n}{2\pi\hbar^2}
      \int e^{-i\mathbf k_f\cdot\mathbf r'}
      \widehat U(\mathbf r')
      \Psi^{(+)}(\mathbf r')\,d^3r'
   \right].

The spatial dependence :math:`e^{ikr}/r` is an outgoing spherical wave.  Its
spinor-valued coefficient contains the dependence on the observation
direction and on the incident state.  The scattering-amplitude operator is
defined as precisely this coefficient,

.. math::

   \boxed{
   \widehat f(\mathbf k_f,\mathbf k_i)\chi_i
   =-\frac{m_n}{2\pi\hbar^2}
   \int e^{-i\mathbf k_f\cdot\mathbf r'}
   \widehat U(\mathbf r')\Psi^{(+)}(\mathbf r')\,d^3r'.}

Linearity of the scattering problem in the incident spinor makes
:math:`\widehat f` a :math:`2\times2` operator in spin space.  The complete
wavefunction therefore has the asymptotic form

.. math::

   \Psi^{(+)}(\mathbf r)
   \sim e^{i\mathbf k_i\cdot\mathbf r}\chi_i
   +\frac{e^{ikr}}{r}
    \widehat f(\mathbf k_f,\mathbf k_i)\chi_i.

The amplitude has units of length.  At this point it is not yet a cross
section; it is the coefficient of the outgoing spherical wave.


From the outgoing amplitude to the cross section
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If the spin analyzer selects the normalized final spinor :math:`\chi_f`, the
corresponding scalar spatial amplitude in the scattered wave is

.. math::

   \psi_{\mathrm{sc},i\rightarrow f}(\mathbf r)
   =\chi_f^\dagger\Psi_{\mathrm{sc}}(\mathbf r)
   \sim\frac{e^{ik_f r}}{r}f_{fi},
   \qquad
   f_{fi}=\chi_f^\dagger
           \widehat f(\mathbf k_f,\mathbf k_i)\chi_i.

For the incident plane wave and normalized :math:`\chi_i`, the current
defined above is

.. math::

   \mathbf j_{\mathrm{in}}
   =\frac{\hbar\mathbf k_i}{m_n},
   \qquad
   j_{\mathrm{in}}=\frac{\hbar k_i}{m_n}.

For the outgoing radial wave,

.. math::

   \frac{\partial}{\partial r}
   \left(\frac{e^{ik_f r}}{r}\right)
   =\left(ik_f-\frac{1}{r}\right)
    \frac{e^{ik_f r}}{r}.

Substitution into the same current formula gives the outward radial current

.. math::

   j_{\mathrm{sc},i\rightarrow f}^{(r)}
   =\frac{\hbar k_f}{m_n}
    \frac{|f_{fi}|^2}{r^2}.

The :math:`1/r^2` decrease is exactly compensated by the
:math:`r^2d\Omega` area of a detector element with fixed solid angle.
Inserting the incident and scattered currents into the operational definition
of the cross section gives

.. math::

   \frac{d\sigma_{i\rightarrow f}}{d\Omega}
   =\frac{k_f}{k_i}|f_{fi}|^2
   =\frac{k_f}{k_i}
    \left|
       \chi_f^\dagger
       \widehat f(\mathbf k_f,\mathbf k_i)
       \chi_i
    \right|^2.

The present derivation is elastic, so :math:`k_f=k_i=k`.  It therefore reduces
to

.. math::

   \boxed{
   \frac{d\sigma_{i\rightarrow f}}{d\Omega}
   =\left|
      \chi_f^{\dagger}
      \widehat f(\mathbf k_f,\mathbf k_i)
      \chi_i
    \right|^2.}

Thus the absolute square is not postulated as the cross section.  It results
from the radial probability current carried by the outgoing spherical wave,
divided by the incident probability current.

Because the current is quadratic in the wavefunction, the current of
:math:`\Psi_{\mathrm{in}}+\Psi_{\mathrm{sc}}` also contains an interference
term; it is not simply the sum of the two currents.  The differential
scattering cross section is defined from the flux carried by the outgoing
spherical component into a detector direction outside the directly
transmitted beam.  There the asymptotic boundary condition separates that
component from the prepared incident plane wave.  In the exact forward
direction, the interference term cannot be discarded: it accounts for the
attenuation of the transmitted beam and leads to the optical theorem.  That
forward-beam relation is distinct from the angular scattered intensity
calculated here.


4. First Born approximation
----------------------------

The integral equation makes successive scattering events explicit.  Repeatedly
substituting its right-hand side into itself gives the Born series.  Suppressing
the spatial integrations for compactness,

.. math::

   \Psi^{(+)}
   =\Psi_{\mathrm{in}}
   +G_0^{(+)}\widehat U\Psi_{\mathrm{in}}
   +G_0^{(+)}\widehat U G_0^{(+)}\widehat U\Psi_{\mathrm{in}}
   +\cdots.

Each additional factor of :math:`\widehat U` represents another interaction
with the sample.  The first correction describes one scattering event; the
following terms contain repeated scattering and the associated distortion of
the wave inside the sample.

The first Born approximation replaces the exact wavefunction *inside the
scattering integral* by the incoming plane wave,

.. math::

   \Psi^{(+)}(\mathbf r)
   \simeq e^{i\mathbf k_i\cdot\mathbf r}\chi_i.

With the Fourier-space argument

.. math::

   \boldsymbol\kappa=\mathbf k_f-\mathbf k_i,

the amplitude becomes

.. math::

   \boxed{
   \widehat f^{(1)}(\boldsymbol\kappa)
   =-\frac{m_n}{2\pi\hbar^2}
   \int e^{-i\boldsymbol\kappa\cdot\mathbf r}
   \widehat U(\mathbf r)\,d^3r.}

This is the central result of the spatial scattering problem: in first order,
the scattering amplitude is proportional to the Fourier transform of the
interaction potential.

.. important:: Momentum-transfer sign

   Some neutron-scattering literature defines the momentum transferred to the
   sample as :math:`\mathbf q=\mathbf k_i-\mathbf k_f`, whereas the Born
   integral above naturally contains
   :math:`\boldsymbol\kappa=\mathbf k_f-\mathbf k_i`.  Thus
   :math:`\boldsymbol\kappa=-\mathbf q`.  Reversing this convention replaces
   every Fourier amplitude by its value at :math:`-\mathbf q`.  For real
   real-space fields this is its complex conjugate.  The convention must be
   changed consistently, especially for the chiral term, which is odd under
   :math:`\mathbf q\mapsto-\mathbf q`.

The first Born approximation is therefore a single-scattering approximation.
It is appropriate when the sample changes the neutron wave only weakly within
the scattering volume, so that the higher terms in the Born series are small.
This is often a good description of conventional transmission SANS because
the neutron interaction is weak and samples are commonly chosen to maintain
high transmission and suppress multiple scattering.  Small scattering angles
alone do not guarantee the approximation.  Thick or strongly scattering
samples, strong refraction, or appreciable multiple scattering can invalidate
it.  In general its validity is controlled by the complete sample and cannot
be reduced to a bound on one local value of the potential.


5. The NuMagSANS scattering operator
-------------------------------------

From this point onward, :math:`\mathbf q` denotes the Fourier-space argument
used by NuMagSANS.  If momentum transfer is defined with the opposite sign in
an experimental convention, the conversion described in the note above must
be applied.

NuMagSANS uses the unitary Fourier transforms

.. math::

   \widetilde N(\mathbf q)
   =\frac{1}{(2\pi)^{3/2}}
   \int N(\mathbf r)e^{-i\mathbf q\cdot\mathbf r}\,d^3r,

.. math::

   \widetilde{\mathbf M}(\mathbf q)
   =\frac{1}{(2\pi)^{3/2}}
   \int \mathbf M(\mathbf r)e^{-i\mathbf q\cdot\mathbf r}\,d^3r.

For :math:`\mathbf q\ne0`, the magnetic induction generated by the sample is
the transverse part of the magnetization.  With the sign convention used by
NuMagSANS, introduce the Halpern--Johnson vector

.. math::

   \boxed{
   \widetilde{\mathbf Q}(\mathbf q)
   =\widehat{\mathbf q}
    \left(\widehat{\mathbf q}\cdot\widetilde{\mathbf M}\right)
   -\widetilde{\mathbf M}
   =-\left(\mathbb I_3-
      \widehat{\mathbf q}\widehat{\mathbf q}^{\mathsf T}\right)
      \widetilde{\mathbf M},}

so that

.. math::

   \widetilde{\mathbf B}=-\mu_0\widetilde{\mathbf Q},
   \qquad
   \widehat{\mathbf q}\cdot\widetilde{\mathbf Q}=0.

Only this result is needed to continue the scattering calculation.  Its
derivation from the magnetostatic Maxwell equations, together with the role
of :math:`\mathbf H` and the magnetization gauge freedom, is given in
:ref:`transverse-magnetization-derivation`.

The magnetic scattering-length coefficient is

.. math::

   b_{\mathrm H}
   =\frac{m_n|\mu_n|\mu_0}{2\pi\hbar^2}
   \simeq2.91\times10^8\ \mathrm{A^{-1}m^{-1}},

when :math:`\mathbf M` is expressed in :math:`\mathrm{A/m}`.

Since :math:`N` has units :math:`\mathrm{m^{-2}}`, its transform
:math:`\widetilde N` has units of length.  Likewise,
:math:`\widetilde{\mathbf M}` has units :math:`\mathrm{A\,m^2}`, so
:math:`b_{\mathrm H}\widetilde{\mathbf Q}` also has units of length.  The
nuclear and magnetic terms can therefore be added at the amplitude level.

The physical moment was written above in an angular-momentum basis as
:math:`\widehat{\boldsymbol\mu}_n
=\mu_n\boldsymbol\sigma^{(S)}` with :math:`\mu_n<0`.  For the NuMagSANS
channel labels it is convenient to introduce

.. math::

   \boldsymbol\sigma=-\boldsymbol\sigma^{(S)}.

The same physical operator is then

.. math::

   \widehat{\boldsymbol\mu}_n=|\mu_n|\boldsymbol\sigma.

Thus a :math:`+` eigenstate of the channel matrix
:math:`\boldsymbol\sigma\cdot\widehat{\mathbf P}` has its magnetic moment
along :math:`\widehat{\mathbf P}` and its angular momentum in the opposite
direction.  This relabeling changes no physics, but it must be retained when
assigning :math:`+` and :math:`-` cross sections.  Using
:math:`\widetilde{\mathbf B}=-\mu_0\widetilde{\mathbf Q}` in the Zeeman
interaction gives

.. math::

   \widetilde{\widehat U}(\mathbf q)
   =\frac{2\pi\hbar^2}{m_n}
   \left[
      \widetilde N(\mathbf q)\mathbb I_2
      +b_{\mathrm H}\boldsymbol\sigma\cdot
       \widetilde{\mathbf Q}(\mathbf q)
   \right].

The first-Born amplitude therefore equals a common factor
:math:`-(2\pi)^{3/2}` times the following coherent scattering-length
operator:

.. math::

   \boxed{
   \widehat{\mathcal A}(\mathbf q)
   =\widetilde N(\mathbf q)\mathbb I_2
   +b_{\mathrm H}\,
    \boldsymbol\sigma\cdot\widetilde{\mathbf Q}(\mathbf q).}

Here :math:`\boldsymbol\sigma` denotes the Pauli matrices in the channel basis
used below.  This one operator contains the nuclear amplitude, the magnetic
amplitude, and their relative phase.

.. important:: Sign and channel convention

   The neutron gyromagnetic ratio is negative.  Its sign can instead be
   carried explicitly in the magnetic-moment operator, but the labels of the
   two channel states must then be reversed.  NuMagSANS fixes the convention
   operationally by the boxed equation above together with
   :math:`\widetilde{\mathbf Q}=\widehat{\mathbf q}
   (\widehat{\mathbf q}\cdot\widetilde{\mathbf M})
   -\widetilde{\mathbf M}` and by defining the first superscript of
   :math:`\Sigma^{s_i s_f}` as the incident state.  Once these definitions are
   chosen, the sign of the nuclear--magnetic interference and the assignment
   of the two chiral spin-flip channels are fixed; they are not individually
   arbitrary.

For a scattering volume :math:`V`, the matrix element between normalized
initial and final channel spinors is converted into the macroscopic
differential cross section according to

.. math::

   \boxed{
   \frac{d\Sigma^{s_i s_f}}{d\Omega}
   =\frac{8\pi^3}{V}
   \left|
   \chi_{s_f}^{\dagger}
   \widehat{\mathcal A}(\mathbf q)
   \chi_{s_i}
   \right|^2.}

The factor :math:`8\pi^3` is not an additional physical interaction.  It is
the square of :math:`(2\pi)^{3/2}`, which appears when the unnormalized Born
integral is expressed through the unitary Fourier transforms above.


6. Spin projection
------------------

Let :math:`\widehat{\mathbf P}` be the real unit vector defining the neutron
quantization axis.  Complete it to a right-handed orthonormal basis
:math:`(\mathbf e_1,\mathbf e_2,\widehat{\mathbf P})`.  In this basis the
three channel Pauli matrices are

.. math::

   \sigma_1=
   \begin{pmatrix}0&1\\1&0\end{pmatrix},
   \qquad
   \sigma_2=
   \begin{pmatrix}0&-i\\i&0\end{pmatrix},
   \qquad
   \sigma_{\parallel}=
   \begin{pmatrix}1&0\\0&-1\end{pmatrix}.

The two normalized channel spinors are explicitly

.. math::

   \chi_+=\begin{pmatrix}1\\0\end{pmatrix},
   \qquad
   \chi_-=\begin{pmatrix}0\\1\end{pmatrix},

and satisfy

.. math::

   (\boldsymbol\sigma\cdot\widehat{\mathbf P})\chi_{\pm}
   =\pm\chi_{\pm}.

The generally complex Halpern--Johnson amplitude can be resolved in the same
real-space basis,

.. math::

   \widetilde{\mathbf Q}
   =Q_1\mathbf e_1+Q_2\mathbf e_2+Q_{\parallel}\widehat{\mathbf P},
   \qquad
   Q_{\parallel}=\widehat{\mathbf P}\cdot\widetilde{\mathbf Q}.

Consequently,

.. math::

   \boldsymbol\sigma\cdot\widetilde{\mathbf Q}
   =\sigma_1Q_1+\sigma_2Q_2+\sigma_{\parallel}Q_{\parallel}
   =\begin{pmatrix}
      Q_{\parallel}&Q_1-iQ_2\\
      Q_1+iQ_2&-Q_{\parallel}
   \end{pmatrix}.

In the :math:`(+,-)` basis the scattering operator becomes

.. math::

   \widehat{\mathcal A}
   =
   \begin{pmatrix}
      \widetilde N+b_{\mathrm H}Q_{\parallel}
      & b_{\mathrm H}(Q_1-iQ_2) \\
      b_{\mathrm H}(Q_1+iQ_2)
      & \widetilde N-b_{\mathrm H}Q_{\parallel}
   \end{pmatrix}.

Rows in this matrix label the final state in the order :math:`(+,-)`, and
columns label the initial state in the same order.  The lower-left element,
for example, is therefore the amplitude for :math:`+\rightarrow-`.
The diagonal entries preserve the channel state and the off-diagonal entries
flip it.  We use the NuMagSANS convention that the first sign in
:math:`s_i s_f` labels the incident neutron and the second sign labels the
scattered neutron.  Thus, for example, :math:`+-` means
:math:`+\rightarrow-`.

It is useful to define the following five real building blocks, with
:math:`K=8\pi^3/V`:

.. math::

   \frac{d\Sigma_{\mathrm N}}{d\Omega}
   =K|\widetilde N|^2,

.. math::

   \frac{d\Sigma_{\mathrm M}}{d\Omega}
   =K b_{\mathrm H}^2
   \widetilde{\mathbf Q}\cdot\widetilde{\mathbf Q}^{\ast},

.. math::

   \frac{d\Sigma_{\mathrm P}}{d\Omega}
   =K b_{\mathrm H}^2
   \left|\widehat{\mathbf P}\cdot
   \widetilde{\mathbf Q}\right|^2,

.. math::

   \frac{d\Sigma_{\mathrm{NM}}}{d\Omega}
   =K b_{\mathrm H}\widehat{\mathbf P}\cdot
   \left(
      \widetilde N\widetilde{\mathbf Q}^{\ast}
      +\widetilde N^{\ast}\widetilde{\mathbf Q}
   \right),

.. math::

   \frac{d\Sigma_{\chi}}{d\Omega}
   =-iK b_{\mathrm H}^2\widehat{\mathbf P}\cdot
   \left(
      \widetilde{\mathbf Q}\times
      \widetilde{\mathbf Q}^{\ast}
   \right).

The cross product in the last expression is purely imaginary, so the chiral
cross section is real.  The Pauli identity

.. math::

   (\boldsymbol\sigma\cdot\mathbf a)
   (\boldsymbol\sigma\cdot\mathbf b)
   =(\mathbf a\cdot\mathbf b)\mathbb I_2
   +i\boldsymbol\sigma\cdot(\mathbf a\times\mathbf b)

then gives all four channels.


Non-spin-flip channels
^^^^^^^^^^^^^^^^^^^^^^

For the diagonal elements one simply expands

.. math::

   |\widetilde N\pm b_{\mathrm H}Q_{\parallel}|^2
   =|\widetilde N|^2+b_{\mathrm H}^2|Q_{\parallel}|^2
   \pm b_{\mathrm H}
   (\widetilde NQ_{\parallel}^{\ast}
   +\widetilde N^{\ast}Q_{\parallel}).

The three terms are the nuclear, longitudinal magnetic, and interference
contributions.  Hence,

.. math::

   \boxed{
   \frac{d\Sigma^{++}}{d\Omega}
   =\frac{d\Sigma_{\mathrm N}}{d\Omega}
   +\frac{d\Sigma_{\mathrm{NM}}}{d\Omega}
   +\frac{d\Sigma_{\mathrm P}}{d\Omega},}

.. math::

   \boxed{
   \frac{d\Sigma^{--}}{d\Omega}
   =\frac{d\Sigma_{\mathrm N}}{d\Omega}
   -\frac{d\Sigma_{\mathrm{NM}}}{d\Omega}
   +\frac{d\Sigma_{\mathrm P}}{d\Omega}.}

The nuclear--magnetic term changes sign because reversing the incident and
analyzed channel changes the sign of the longitudinal magnetic amplitude but
not the nuclear amplitude.


Spin-flip channels
^^^^^^^^^^^^^^^^^^

Only the part of :math:`\widetilde{\mathbf Q}` perpendicular to
:math:`\widehat{\mathbf P}` can connect the two channel states.  Define

.. math::

   \frac{d\Sigma_{\mathrm{sf}}}{d\Omega}
   =\frac{d\Sigma_{\mathrm M}}{d\Omega}
   -\frac{d\Sigma_{\mathrm P}}{d\Omega}.

The difference above equals
:math:`Kb_{\mathrm H}^2(|Q_1|^2+|Q_2|^2)`.  The imbalance between the two
off-diagonal matrix elements follows from

.. math::

   |Q_1\pm iQ_2|^2
   =|Q_1|^2+|Q_2|^2
   \mp i(Q_1Q_2^{\ast}-Q_2Q_1^{\ast}),

with

.. math::

   \widehat{\mathbf P}\cdot
   (\widetilde{\mathbf Q}\times\widetilde{\mathbf Q}^{\ast})
   =Q_1Q_2^{\ast}-Q_2Q_1^{\ast}.

The off-diagonal matrix elements therefore give

.. math::

   \boxed{
   \frac{d\Sigma^{+-}}{d\Omega}
   =\frac{d\Sigma_{\mathrm{sf}}}{d\Omega}
   +\frac{d\Sigma_{\chi}}{d\Omega},}

.. math::

   \boxed{
   \frac{d\Sigma^{-+}}{d\Omega}
   =\frac{d\Sigma_{\mathrm{sf}}}{d\Omega}
   -\frac{d\Sigma_{\chi}}{d\Omega}.}

Consequently, the polarization-independent spin-flip and chiral terms are the
half-sum and half-difference of the measurable spin-flip channels,

.. math::

   \frac{d\Sigma_{\mathrm{sf}}}{d\Omega}
   =\frac{1}{2}
   \left(
      \frac{d\Sigma^{+-}}{d\Omega}
      +\frac{d\Sigma^{-+}}{d\Omega}
   \right),

.. math::

   \frac{d\Sigma_{\chi}}{d\Omega}
   =\frac{1}{2}
   \left(
      \frac{d\Sigma^{+-}}{d\Omega}
      -\frac{d\Sigma^{-+}}{d\Omega}
   \right).


Unpolarized and half-polarized measurements
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Averaging over the two equally populated incident states and summing over the
unobserved final state removes the nuclear--magnetic and chiral terms.  The
result is the sum of the nuclear and unpolarized magnetic contributions,

.. math::

   \frac{1}{2}
   \sum_{s_i,s_f}
   \frac{d\Sigma^{s_i s_f}}{d\Omega}
   =\frac{d\Sigma_{\mathrm N}}{d\Omega}
   +\frac{d\Sigma_{\mathrm M}}{d\Omega}.

If the incident beam is polarized but the outgoing spin is not analyzed, the
SANSPOL cross sections are

.. math::

   \frac{d\Sigma^+}{d\Omega}
   =\frac{d\Sigma^{++}}{d\Omega}
   +\frac{d\Sigma^{+-}}{d\Omega},
   \qquad
   \frac{d\Sigma^-}{d\Omega}
   =\frac{d\Sigma^{--}}{d\Omega}
   +\frac{d\Sigma^{-+}}{d\Omega}.

These are precisely the combinations exported by the corresponding
NuMagSANS output options.


7. NuMagSANS SANS geometry
--------------------------

NuMagSANS fixes the incident beam along the :math:`x`-axis.  In the
small-angle approximation the detector lies in the :math:`yz`-plane and

.. math::

   \mathbf q
   =q\begin{pmatrix}0\\ \sin\theta\\ \cos\theta\end{pmatrix},
   \qquad
   \widehat{\mathbf q}
   =\begin{pmatrix}0\\ \sin\theta\\ \cos\theta\end{pmatrix}.

Inserting this direction into the transverse projector gives the components
actually evaluated in the kernels,

.. math::

   \widetilde{\mathbf Q}
   =\begin{pmatrix}
      -\widetilde M_x \\
      (\widetilde M_z\sin\theta
      -\widetilde M_y\cos\theta)\cos\theta \\
      (\widetilde M_y\cos\theta
      -\widetilde M_z\sin\theta)\sin\theta
   \end{pmatrix}.

The general spin-projection formulas above remain valid for any real unit
polarization vector :math:`\widehat{\mathbf P}`.  Choosing a particular guide
field direction simply specifies which component of
:math:`\widetilde{\mathbf Q}` is longitudinal in spin space.

The complete list of output names, formulas, symmetries, and configuration
examples is given under :ref:`output-options` in the
:doc:`ParameterOverview`.


8. Equivalent derivation from Fermi's golden rule
-------------------------------------------------

For completeness, take the initial and final plane waves to be normalized in
a box of volume :math:`\mathcal V`,

.. math::

   \Psi_i(\mathbf r)
   =\frac{e^{i\mathbf k_i\cdot\mathbf r}}{\sqrt{\mathcal V}}\chi_i,
   \qquad
   \Psi_f(\mathbf r)
   =\frac{e^{i\mathbf k_f\cdot\mathbf r}}{\sqrt{\mathcal V}}\chi_f.

The interaction matrix element between these two waves is

.. math::

   U_{fi}
   =\frac{1}{\mathcal V}
    \int e^{-i\boldsymbol\kappa\cdot\mathbf r}
    \chi_f^{\dagger}\widehat U(\mathbf r)\chi_i\,d^3r,
   \qquad
   \boldsymbol\kappa=\mathbf k_f-\mathbf k_i.

To first order, Fermi's golden rule gives

.. math::

   w_{i\rightarrow f}
   =\frac{2\pi}{\hbar}
   |U_{fi}|^2
   \delta(E_f-E_i).

Using the final-state density
:math:`\mathcal Vd^3k_f/(2\pi)^3`, integrating the energy delta function, and
dividing by the incident flux yields the cross section.  This is the same
flux concept derived from the continuity equation in Section 3.  The
box-normalized incident wave has density :math:`1/\mathcal V` and therefore
incident current
:math:`j_i=\hbar k_i/(m_n\mathcal V)`.  One then obtains

.. math::

   \frac{d\sigma_{i\rightarrow f}}{d\Omega}
   =\frac{k_f}{k_i}
   \left(\frac{m_n}{2\pi\hbar^2}\right)^2
   \left|
      \int e^{-i\boldsymbol\kappa\cdot\mathbf r}
      \chi_f^{\dagger}\widehat U(\mathbf r)\chi_i\,d^3r
   \right|^2.

For elastic scattering :math:`k_f/k_i=1`, and this is exactly the modulus
squared of the first-Born amplitude obtained from the Lippmann--Schwinger
equation.  The two routes therefore contain the same first-order physics:
one extracts an outgoing wave, while the other computes its transition rate.


.. _transverse-magnetization-derivation:

9. Why only transverse magnetization scatters
---------------------------------------------

This section derives the transverse projection used in the NuMagSANS
scattering operator above.  The projection follows from Maxwell's equations
rather than from an additional assumption about the neutron.  It can be
obtained without imposing a pointwise relation among :math:`\mathbf B`,
:math:`\mathbf H`, and :math:`\mathbf M`.  In magnetostatics, the fundamental
equations for the magnetic induction are

.. math::

   \nabla\cdot\mathbf B=0,
   \qquad
   \nabla\times\mathbf B=\mu_0\mathbf J_{\mathrm{tot}},

where :math:`\mathbf J_{\mathrm{tot}}` is the total electric-current density.
A magnetization :math:`\mathbf M`, understood as a magnetic-dipole-moment
density, represents the bound magnetization-current density

.. math::

   \mathbf J_{\mathrm M}=\nabla\times\mathbf M.

If :math:`\mathbf M` is extended by zero outside a finite sample, this equation
is understood distributionally: a discontinuity at the sample boundary also
contains the corresponding surface-current contribution.  The total current
can therefore be separated as

.. math::

   \mathbf J_{\mathrm{tot}}
   =\mathbf J_{\mathrm f}+\nabla\times\mathbf M.

NuMagSANS considers the induction generated by the sample magnetization and
assumes that no additional free electric current flows through the sample,
so :math:`\mathbf J_{\mathrm f}=0`.  The relevant real-space Maxwell equations
are then

.. math::

   \boxed{
   \nabla\cdot\mathbf B=0,
   \qquad
   \nabla\times\mathbf B
   =\mu_0\nabla\times\mathbf M.}

Taking the curl of the second equation and using the vector-Laplacian identity

.. math::

   \boxed{
   \nabla\times(\nabla\times\mathbf F)
   =\nabla(\nabla\cdot\mathbf F)-\nabla^2\mathbf F}

gives a purely differential relation between the induction and the
magnetization:

.. math::

   -\nabla^2\mathbf B
   =\mu_0\left[
      \nabla(\nabla\cdot\mathbf M)-\nabla^2\mathbf M
      \right].

Equivalently,

.. math::

   \boxed{
   \nabla^2\mathbf B
   =\mu_0\left[
      \nabla^2\mathbf M-\nabla(\nabla\cdot\mathbf M)
      \right]
    =-\mu_0\nabla\times(\nabla\times\mathbf M).}

Unlike a pointwise constitutive equation, this differential equation leaves
harmonic field components to be fixed by boundary conditions.  In particular,
it does not constrain the spatially constant mode.

For the Fourier convention
:math:`\widetilde f(\mathbf q)\propto\int f(\mathbf r)
e^{-i\mathbf q\cdot\mathbf r}\,d^3r`, a spatial derivative transforms as
:math:`\nabla\rightarrow i\mathbf q`, provided boundary terms vanish or
periodic boundary conditions are used.  Hence, for :math:`\mathbf q\ne0`,

.. math::

   \mathbf q\cdot\widetilde{\mathbf B}=0,
   \qquad
   \boxed{
   \mathbf q\times\widetilde{\mathbf B}
   =\mu_0\mathbf q\times\widetilde{\mathbf M}.}

To solve for :math:`\widetilde{\mathbf B}`, take the cross product of the
second equation with :math:`\mathbf q` from the left.  The vector
triple-product identity

.. math::

   \boxed{
   \mathbf a\times(\mathbf b\times\mathbf c)
   =\mathbf b(\mathbf a\cdot\mathbf c)
   -\mathbf c(\mathbf a\cdot\mathbf b)}

is often remembered as the ``BAC--CAB`` rule.  Applying it to both sides gives

.. math::

   \mathbf q\times
   (\mathbf q\times\widetilde{\mathbf B})
   =\mathbf q(\mathbf q\cdot\widetilde{\mathbf B})
    -q^2\widetilde{\mathbf B}
   =-q^2\widetilde{\mathbf B},

where :math:`\mathbf q\cdot\widetilde{\mathbf B}=0` was used.  The right-hand
side is

.. math::

   \mu_0\mathbf q\times
   (\mathbf q\times\widetilde{\mathbf M})
   =\mu_0\left[
      \mathbf q(\mathbf q\cdot\widetilde{\mathbf M})
      -q^2\widetilde{\mathbf M}
   \right].

Dividing by :math:`-q^2` yields

.. math::

   \widetilde{\mathbf B}(\mathbf q)
   =\mu_0\left(\mathbb I_3
   -\widehat{\mathbf q}\widehat{\mathbf q}^{\mathsf T}\right)
   \widetilde{\mathbf M}(\mathbf q).

The matrix
:math:`\boldsymbol\Pi_{\perp}=\mathbb I_3-\widehat{\mathbf q}
\widehat{\mathbf q}^{\mathsf T}` is a projector: it removes the component of
:math:`\widetilde{\mathbf M}` parallel to :math:`\mathbf q` and leaves the two
perpendicular components unchanged.  Indeed,

.. math::

   \boldsymbol\Pi_{\perp}^2=\boldsymbol\Pi_{\perp},
   \qquad
   \boldsymbol\Pi_{\perp}\widehat{\mathbf q}=0.

The complementary longitudinal projector is

.. math::

   \boldsymbol\Pi_{\parallel}
   =\widehat{\mathbf q}\widehat{\mathbf q}^{\mathsf T},
   \qquad
   \boldsymbol\Pi_{\perp}+\boldsymbol\Pi_{\parallel}=\mathbb I_3.

Writing

.. math::

   \widetilde{\mathbf M}_{\parallel}
   =\widehat{\mathbf q}
    (\widehat{\mathbf q}\cdot\widetilde{\mathbf M}),
   \qquad
   \widetilde{\mathbf M}_{\perp}
   =\widetilde{\mathbf M}-\widetilde{\mathbf M}_{\parallel},

the field equation simply states
:math:`\widetilde{\mathbf B}=\mu_0\widetilde{\mathbf M}_{\perp}`.
Equivalently, the triple-product identity gives the double-cross-product form

.. math::

   \widehat{\mathbf q}\times
   (\widehat{\mathbf q}\times\widetilde{\mathbf M})
   =\widehat{\mathbf q}
    (\widehat{\mathbf q}\cdot\widetilde{\mathbf M})
    -\widetilde{\mathbf M}.

NuMagSANS uses the Halpern--Johnson vector

.. math::

   \boxed{
   \widetilde{\mathbf Q}(\mathbf q)
   =\widehat{\mathbf q}
    \left(\widehat{\mathbf q}\cdot\widetilde{\mathbf M}\right)
   -\widetilde{\mathbf M}
   =-\boldsymbol\Pi_{\perp}\widetilde{\mathbf M}.}

Consequently,

.. math::

   \widetilde{\mathbf B}=-\mu_0\widetilde{\mathbf Q},
   \qquad
   \widehat{\mathbf q}\cdot\widetilde{\mathbf Q}=0.

The longitudinal component of :math:`\widetilde{\mathbf M}` therefore cannot
contribute to the magnetic scattering amplitude.


Role of :math:`\mathbf H` and the scattering null space
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The auxiliary field :math:`\mathbf H` is not needed to derive the transverse
magnetic scattering amplitude.  Its useful role here is instead to identify
the part of :math:`\mathbf M` that lies in the null space of the scattering
map.  Introduce :math:`\mathbf H` through its differential source equations,

.. math::

   \nabla\times\mathbf H=\mathbf J_{\mathrm f},
   \qquad
   \nabla\cdot\mathbf H=-\nabla\cdot\mathbf M.

Together with the differential equations for :math:`\mathbf B`, these imply

.. math::

   \boxed{
   \nabla^2\mathbf B
   =\mu_0\nabla^2(\mathbf H+\mathbf M).}

Indeed, the difference
:math:`\mathbf B-\mu_0(\mathbf H+\mathbf M)` has both zero divergence and zero
curl and is therefore harmonic.  The differential equations do not require
this harmonic contribution to vanish; doing so would be an additional choice
of boundary condition or of the homogeneous mode.

For :math:`\mathbf J_{\mathrm f}=0` and :math:`\mathbf q\ne0`, the equations
for :math:`\mathbf H` become

.. math::

   \mathbf q\times\widetilde{\mathbf H}=0,
   \qquad
   \mathbf q\cdot\widetilde{\mathbf H}
   =-\mathbf q\cdot\widetilde{\mathbf M}.

The first equation makes :math:`\widetilde{\mathbf H}` longitudinal; the
second fixes its magnitude.  Hence

.. math::

   \boxed{
   \widetilde{\mathbf H}
   =-\boldsymbol\Pi_{\parallel}\widetilde{\mathbf M},
   \qquad
   \frac{\widetilde{\mathbf B}}{\mu_0}
   =\boldsymbol\Pi_{\perp}\widetilde{\mathbf M}.}

This is the longitudinal--transverse split relevant to magnetic SANS.  The
transverse part is represented by the SANS-active induction, whereas the
longitudinal part satisfies
:math:`\boldsymbol\Pi_{\parallel}\widetilde{\mathbf M}
=-\widetilde{\mathbf H}` and lies in the null space of the magnetic scattering
operator.  The minus sign follows from
:math:`\nabla\cdot\mathbf H=-\nabla\cdot\mathbf M`.


Gauge freedom of the magnetization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The transverse projection also exposes a non-uniqueness of the magnetization
that is relevant to magnetic scattering.  Let :math:`\chi(\mathbf r)` be a
sufficiently regular scalar field, with boundary conditions compatible with
the Fourier transform.  The transformation

.. math::

   \boxed{
   \mathbf M(\mathbf r)
   \longrightarrow
   \mathbf M'(\mathbf r)
   =\mathbf M(\mathbf r)+\nabla\chi(\mathbf r)}

does not change the bound magnetization-current density, because

.. math::

   \nabla\times\mathbf M'
   =\nabla\times\mathbf M
    +\underbrace{\nabla\times\nabla\chi}_{=0}
   =\nabla\times\mathbf M.

At fixed boundary conditions it therefore leaves the magnetostatic Maxwell
equations above, and hence every finite-:math:`q` component of
:math:`\mathbf B`, unchanged.  The auxiliary longitudinal field transforms in
the opposite direction,

.. math::

   \mathbf H' = \mathbf H-\nabla\chi,

so that its differential source equations remain satisfied.  The invariance
of the induction is especially transparent in Fourier space:

.. math::

   \widetilde{\mathbf M}'
   =\widetilde{\mathbf M}
    +i\mathbf q\widetilde\chi,

.. math::

   \begin{aligned}
   \boldsymbol\Pi_{\perp}\widetilde{\mathbf M}'
   &=\boldsymbol\Pi_{\perp}\widetilde{\mathbf M}
     +i\widetilde\chi\,
      \underbrace{\boldsymbol\Pi_{\perp}\mathbf q}_{=0} \\
   &=\boldsymbol\Pi_{\perp}\widetilde{\mathbf M}.
   \end{aligned}

Consequently,

.. math::

   \widetilde{\mathbf B}[\mathbf M']
   =\widetilde{\mathbf B}[\mathbf M],
   \qquad
   \widetilde{\mathbf Q}[\mathbf M']
   =\widetilde{\mathbf Q}[\mathbf M].

This is the magnetization gauge freedom meant here: at every
:math:`\mathbf q\ne0`, magnetic SANS distinguishes only the equivalence class

.. math::

   \mathbf M\sim\mathbf M+\nabla\chi,

not the longitudinal part of :math:`\mathbf M` separately.  This statement
does not mean that the microscopic magnetization is physically arbitrary.
It means that the scattering map from :math:`\mathbf M` to
:math:`\widetilde{\mathbf Q}` is not injective: magnetization fields that
differ only by such a longitudinal contribution produce the same magnetic
SANS signal.

The relation can also be expressed through the Helmholtz decomposition

.. math::

   \mathbf M=\mathbf M_{\perp}+\nabla\phi,
   \qquad
   \nabla\cdot\mathbf M_{\perp}=0.

Taking the divergence gives the Poisson equation

.. math::

   \nabla^2\phi=\nabla\cdot\mathbf M,

whose solution depends on the chosen boundary conditions.  Subtracting
:math:`\nabla\phi` is a possible gauge choice; it selects the transverse
representative.  In Fourier space this representative is uniquely given by

.. math::

   \widetilde{\mathbf M}_{\perp}
   =\boldsymbol\Pi_{\perp}\widetilde{\mathbf M},
   \qquad
   \mathbf q\cdot\widetilde{\mathbf M}_{\perp}=0.

This Helmholtz decomposition is the real-space counterpart of the
double-cross-product projector derived above.

All projector and gauge statements in this section apply to
:math:`\mathbf q\ne0`.  A magnetization :math:`\mathbf M_0` that is constant
over all space belongs entirely to the :math:`\mathbf q=0` mode.  The
differential Maxwell equations do not then enforce
:math:`\mathbf H_0=-\mathbf M_0`; the homogeneous components of
:math:`\mathbf B` and :math:`\mathbf H` must be fixed by boundary conditions.
By contrast, a magnetization that is constant only inside a finite sample,
:math:`\mathbf M(\mathbf r)=\mathbf M_0\mathbb 1_V(\mathbf r)`, is not a pure
zero mode.  Its boundary and shape function generate finite-:math:`q`
components that can contribute to SANS.

A spatially uniform external guide field likewise belongs to the separate
:math:`\mathbf q=0` mode.  It sets the spin quantization axis but is not part of
the finite-:math:`q` SANS signal.  If free currents are present, their
Fourier-space contribution must be added and the magnetization-only projector
is no longer the complete field solution.


10. What changes beyond this model?
-----------------------------------

The compact NuMagSANS formulas are a consequence of the assumptions stated at
the beginning.  Relaxing them changes the problem as follows:

- **Inelastic scattering:** the cross section depends on energy transfer and
  time-dependent correlation functions, and contains the factor
  :math:`k_f/k_i`.
- **Spin-dependent nuclear scattering:** the nuclear term is no longer a
  scalar multiple of :math:`\mathbb I_2`; nuclear spin-incoherent and
  nuclear--magnetic terms must be added.
- **Multiple scattering or strong refraction:** the incident plane wave cannot
  replace the exact wave inside the sample; the first Born approximation is
  insufficient.
- **Statistical ensembles:** amplitudes are added coherently within one
  realization, whereas the physically appropriate ensemble average is
  normally taken over intensities or correlation functions afterward.
- **Instrumental effects:** finite polarization, analyzer efficiency,
  wavelength spread, detector resolution, and background act on the ideal
  cross sections after the sample-scattering calculation.

Successful execution of NuMagSANS verifies the deterministic calculation of
the selected model.  Scientific validation still requires suitable
analytical limits, symmetry relations, convergence studies, and comparison
with reference data or experiment.


References
----------

- O. Halpern and M. H. Johnson, `On the Magnetic Scattering of Neutrons
  <https://doi.org/10.1103/PhysRev.55.898>`_, *Physical Review* **55**, 898
  (1939).
- B. A. Lippmann and J. Schwinger, `Variational Principles for Scattering
  Processes. I <https://doi.org/10.1103/PhysRev.79.469>`_, *Physical Review*
  **79**, 469 (1950).
- R. M. Moon, T. Riste, and W. C. Koehler, `Polarization Analysis of
  Thermal-Neutron Scattering <https://doi.org/10.1103/PhysRev.181.920>`_,
  *Physical Review* **181**, 920 (1969).
- M. P. Adams and A. Michels, `NuMagSANS: a GPU-accelerated open-source
  software package for the generic computation of nuclear and magnetic
  small-angle neutron scattering observables of complex systems
  <https://doi.org/10.1107/S160057672600258X>`_, *Journal of Applied
  Crystallography* **59**, 943--959 (2026).
