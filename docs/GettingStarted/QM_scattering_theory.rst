.. _quantum-scattering-theory:

Quantum Scattering Theory
=========================

.. math::

   \boxed{
   \begin{array}{c}
      \text{stationary Schrödinger equation} \\
      \Downarrow \\
      \text{Lippmann--Schwinger equation} \\
      \Downarrow \\
      \text{far-field asymptotics} \\
      \Downarrow \\
      \text{scattering amplitude} \\
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

The corresponding stationary Schrödinger equation is

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
sample.  Its connection to :math:`\mathbf M` and its explicit matrix form are
derived below.


2. Lippmann--Schwinger equation
--------------------------------

Using the Hamiltonian defined above, separate the free propagation from the
sample interaction as

.. math::

   \widehat H_0=-\frac{\hbar^2}{2m_n}\nabla^2\mathbb I_2,
   \qquad
   \widehat H=\widehat H_0+\widehat U.

For a free neutron, :math:`E=\hbar^2k^2/(2m_n)`.  The symbol :math:`k` is the
magnitude of the elastic incident and final wavevectors.

The stationary Schrödinger equation can then be rewritten exactly as
the outgoing Lippmann--Schwinger equation

.. math::

   |\Psi^{(+)}\rangle
   =|\Phi\rangle
   +\widehat G_0^{(+)}(E)\widehat U|\Psi^{(+)}\rangle,
   \qquad
   \widehat G_0^{(+)}(E)
   =\left(E-\widehat H_0+i0^+\right)^{-1}.

Here :math:`|\Psi^{(+)}\rangle` denotes the complete quantum state in
abstract bra--ket notation; :math:`\Psi^{(+)}(\mathbf r)` below is the same
state represented as a two-component position-space wavefunction.
:math:`\widehat G_0^{(+)}(E)` is the free-particle Green operator, or
resolvent.

The infinitesimal positive imaginary term :math:`i0^+` specifies the outgoing
boundary condition: disturbances created by the sample propagate away from
it rather than toward it.

For an incoming plane wave with wavevector :math:`\mathbf k_i` and normalized
spinor :math:`\chi_i`,

.. math::

   \Phi(\mathbf r)=e^{i\mathbf k_i\cdot\mathbf r}\chi_i.

The outgoing free Green function in three dimensions gives

.. math::

   \Psi^{(+)}(\mathbf r)
   =\Phi(\mathbf r)
   -\frac{m_n}{2\pi\hbar^2}
   \int
   \frac{e^{ik|\mathbf r-\mathbf r'|}}
        {|\mathbf r-\mathbf r'|}
   \widehat U(\mathbf r')\Psi^{(+)}(\mathbf r')\,d^3r'.

This integral equation is still exact within the effective Hamiltonian above.
In particular, no Born approximation has yet been made.


3. Far field and scattering amplitude
-------------------------------------

Far from a localized sample, :math:`r\gg|\mathbf r'|`, one has

.. math::

   |\mathbf r-\mathbf r'|
   =r-\widehat{\mathbf r}\cdot\mathbf r'
    +O(|\mathbf r'|^2/r),

where :math:`r=|\mathbf r|` is the detector distance,
:math:`\widehat{\mathbf r}=\mathbf r/r` is the observation direction, and
:math:`\mathbf r'` is a point inside the sample.  The notation
:math:`O(|\mathbf r'|^2/r)` collects the next-order corrections in the ratio
of sample size to detector distance.  Therefore,

.. math::

   \frac{e^{ik|\mathbf r-\mathbf r'|}}
        {|\mathbf r-\mathbf r'|}
   \sim
   \frac{e^{ikr}}{r}e^{-i\mathbf k_f\cdot\mathbf r'},
   \qquad
   \mathbf k_f=k\widehat{\mathbf r}.

The asymptotic wavefunction is the sum of the incident plane wave and an
outgoing spherical wave,

.. math::

   \Psi^{(+)}(\mathbf r)
   \sim e^{i\mathbf k_i\cdot\mathbf r}\chi_i
   +\frac{e^{ikr}}{r}
   \widehat f(\mathbf k_f,\mathbf k_i)\chi_i.

The exact scattering amplitude is determined by

.. math::

   \widehat f(\mathbf k_f,\mathbf k_i)\chi_i
   =-\frac{m_n}{2\pi\hbar^2}
   \int e^{-i\mathbf k_f\cdot\mathbf r}
   \widehat U(\mathbf r)\Psi^{(+)}(\mathbf r)\,d^3r.

It is a matrix in spin space: the sample can scatter the spatial wave and can
also change the neutron spin state.  The scattering amplitude has units of
length.  Comparing the outgoing and incident probability currents gives, for
normalized initial and final spinors and elastic scattering,

.. math::

   \frac{d\sigma_{i\rightarrow f}}{d\Omega}
   =\left|
      \chi_f^{\dagger}
      \widehat f(\mathbf k_f,\mathbf k_i)
      \chi_i
    \right|^2.

Thus the spin-resolved cross section is obtained by projecting the amplitude
onto the selected final spinor before taking its absolute square.


4. First Born approximation
----------------------------

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

The first Born approximation is a single-scattering approximation.  Its
validity is controlled by the size of higher terms in the Born series for the
complete sample and cannot in general be reduced to a bound on one local
value of the potential.


5. Why only transverse magnetization scatters
---------------------------------------------

The magnetic projection follows from Maxwell's equations rather than from an
additional assumption about the neutron.  It can be derived without imposing
a pointwise relation among :math:`\mathbf B`, :math:`\mathbf H`, and
:math:`\mathbf M`.  In magnetostatics, the fundamental equations for the
magnetic induction are

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

The magnetic scattering-length coefficient is

.. math::

   b_{\mathrm H}
   =\frac{m_n|\mu_n|\mu_0}{2\pi\hbar^2}
   \simeq2.91\times10^8\ \mathrm{A^{-1}m^{-1}},

when :math:`\mathbf M` is expressed in :math:`\mathrm{A/m}`.


6. The NuMagSANS scattering operator
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


7. Spin projection
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


8. NuMagSANS SANS geometry
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


9. Equivalent derivation from Fermi's golden rule
-------------------------------------------------

For completeness, take plane waves normalized in a box of volume
:math:`\mathcal V`,

.. math::

   \langle\mathbf r|\mathbf k,s\rangle
   =\frac{e^{i\mathbf k\cdot\mathbf r}}{\sqrt{\mathcal V}}\chi_s.

To first order, Fermi's golden rule gives

.. math::

   w_{i\rightarrow f}
   =\frac{2\pi}{\hbar}
   |\langle f|\widehat U|i\rangle|^2
   \delta(E_f-E_i).

Using the final-state density
:math:`\mathcal Vd^3k_f/(2\pi)^3`, integrating the energy delta function, and
dividing by the incident flux
:math:`j_i=\hbar k_i/(m_n\mathcal V)` yields

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
