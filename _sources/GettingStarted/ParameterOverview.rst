Parameter Overview
==================

NuMagSANS simulations are configured using the Python function

.. code-block:: python

    NuMagSANS.write_config(...)

This function generates a configuration file that is passed to the
NuMagSANS backend executable. The parameters listed below correspond
directly to entries in the generated configuration file and control
the physical model and numerical settings of the simulation.


Example Script
--------------

The following script demonstrates a minimal NuMagSANS workflow.

.. code-block:: python

   from NuMagSANS import NuMagSANS
   from pathlib import Path

   BASE_DIR = Path(__file__).resolve().parent

   # Create NuMagSANS object
   sim = NuMagSANS()

   # Temporary configuration file
   config = BASE_DIR / "NuMagSANSInput_temp.conf"

   # Write configuration file
   sim.write_config(
       config,
       MagData_activate=1,
       MagDataPath=str(BASE_DIR / "RealSpaceData" / "MagData"),
       foldernameSANSData=str(BASE_DIR / "NuMagSANS_Output"),
       Fourier_Approach="atomistic",
       User_Selection=[1, 2, 3],
       Scattering_Volume_V=2.618e-24,
       enable_outputs=["SpinFlip_2D", "SpinFlip_1D"]
   )

   # Run NuMagSANS simulation
   sim.run(config)

   # Remove temporary configuration file
   sim.config_clear(config)


Data Paths and Data Selection
-----------------------------

These parameters define where the simulation reads input data and where
results are written.

NuMagSANS separates local object data from optional assembly metadata:

.. list-table::
   :header-rows: 1

   * - Data layer
     - Purpose
     - Typical file layout
   * - ``MagData``
     - Local magnetic moment data of one or more objects.
     - ``MagData/Object_1/m_1.csv``, ``MagData/Object_2/m_1.csv``, ...
   * - ``NucData``
     - Local nuclear scattering-length-density data of one or more objects.
     - ``NucData/Object_1/n_1.csv``, ``NucData/Object_2/n_1.csv``, ...
   * - ``StructData``
     - Object-center translations. One row per object with ``x y z``.
     - ``StructData.csv``
   * - ``RotData``
     - Object-wise local rotations. One row per object with ZYZ Euler angles ``alpha beta gamma`` in radians.
     - ``RotData.csv`` or ``RotData/RotData_1.csv``, ``RotData/RotData_2.csv``, ...

This separation allows several system types to be represented with the same
backend: fully materialized micromagnetic or atomistic datasets can be imported
directly through ``MagData`` and/or ``NucData``, while object-based assemblies
can additionally use ``StructData`` and ``RotData`` to translate and rotate
local object data without rewriting the object files.

``NucDataPath``  
    Default path ``RealSpaceData/NucData``

``NucData_activate``
    Default value 0

``MagDataPath`` 
    Default path ``RealSpaceData/MagData``

``MagData_activate``
    Default value 0

``StructDataFilename``
    Default path ``RealSpaceData/StructData.csv``

``StructData_activate``
    Default value 0

``StructDataPath``
    Default path ``RealSpaceData/StructData``

    Folder used by ``StructDataLoop``. In this mode, structure files are
    expected to follow the naming convention ``StructData_1.csv``,
    ``StructData_2.csv``, ... inside ``StructDataPath``. Each file contains one
    object-center translation table.

``RotDataFilename``
    Default path ``RealSpaceData/RotData.csv``

``RotData_activate``
    Default value 0

    If activated, ``RotDataFilename`` is read as an object-wise rotation table.
    The file must contain one row per object and three columns:
    ``alpha beta gamma``. The convention is
    :math:`R = R_z(\alpha) R_y(\beta) R_z(\gamma)`, with angles in radians.

``RotDataPath``
    Default path ``RealSpaceData/RotData``

    Folder used by ``RotDataLoop``. In this mode, rotation files are expected
    to follow the naming convention ``RotData_1.csv``, ``RotData_2.csv``, ...
    inside ``RotDataPath``. Each file contains one object-wise rotation table.

``foldernameSANSData``
    Default ``NuMagSANS_Output``

``FastLoad``
    Default value 0

    If set to ``1``, NuMagSANS skips repeated dimension checks for all files
    beyond the first data file inside each object folder. This can accelerate
    import for large datasets. It assumes that, within each object folder, all
    magnetic or nuclear data files have the same number of rows. Different
    objects may still have different numbers of rows.

``Exclude_Zero_Moments``
    Default value 0

    If set to ``1``, rows with zero magnetic moment are excluded from imported
    magnetic datasets.


Supported Data-Layer Combinations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The atomistic backend supports pure magnetic, pure nuclear, and combined
nuclear-magnetic calculations. Each of these can be evaluated without assembly
metadata, with ``StructData``, with ``RotData``, or with both ``StructData`` and
``RotData``.

The main combinations are:

.. list-table::
   :header-rows: 1

   * - Active layers
     - Interpretation
   * - ``MagData``
     - Pure magnetic scattering from local real-space data.
   * - ``NucData``
     - Pure nuclear scattering from local real-space data.
   * - ``MagData + NucData``
     - Combined nuclear-magnetic scattering from local real-space data.
   * - ``... + StructData``
     - Local objects are translated by object-center positions.
   * - ``... + RotData``
     - Local objects are rotated individually before scattering is computed.
   * - ``... + StructData + RotData``
     - Local objects are individually rotated and then placed at object-center positions.

This makes it possible to use NuMagSANS both for large fully materialized
simulation datasets and for object-based assembly models where local objects
are reused with different positions and orientations.


Fourier Approach
----------------

Defines how the Fourier transform from real space to reciprocal space
is performed.

``Fourier_Approach``
    Default ``atomistic``
    Options: ``atomistic``, ``micromagnetic``

Loop Control
------------

These parameters control batch simulations or repeated calculations.

``Loop_Modus``
    Description: Enable looping over selected data-file indices.
    Default ``0``

``Loop_From``
    Description: First index used in loop simulations 
    Default ``1``

``Loop_To``
    Description: Last index used in loop simulations
    Default ``20``

``User_Selection``
    Description: Explicit list of selected indices
    Default ``[1]``

    If ``Loop_Modus = 0``, only the listed file indices are evaluated. If
    ``Loop_Modus = 1``, the index range from ``Loop_From`` to ``Loop_To`` is
    used.

``StructDataLoop``
    Description: Enable an inner loop over several structure-data files.
    Default ``0``

    If ``StructDataLoop = 1``, NuMagSANS loads the selected
    ``StructData_#.csv`` files from ``StructDataPath`` while keeping the
    currently selected magnetic and nuclear local object data fixed. If
    ``RotDataLoop`` is also active, the structure-data loop is the outer loop
    and the rotation-data loop is the inner loop.

``StructDataLoop_From``
    Description: First structure-data file index used by the inner StructData
    loop.
    Default ``1``

``StructDataLoop_To``
    Description: Last structure-data file index used by the inner StructData
    loop.
    Default ``1``

``StructData_User_Selection``
    Description: Optional explicit list of selected structure-data file
    indices.

    If this option is present, it selects the exact ``StructData_#.csv`` files
    used by the inner StructData loop, for example ``{1, 3, 4}``. If it is not
    present, the contiguous range from ``StructDataLoop_From`` to
    ``StructDataLoop_To`` is used.

``RotDataLoop``
    Description: Enable an inner loop over several rotation-data files.
    Default ``0``

    If ``RotDataLoop = 1``, NuMagSANS loads the selected ``RotData_#.csv``
    files from ``RotDataPath`` while keeping the currently selected magnetic,
    nuclear, and structure data fixed. This is useful for evaluating the same
    local object data under multiple object-wise rotation configurations.

``RotDataLoop_From``
    Description: First rotation-data file index used by the inner RotData loop.
    Default ``1``

``RotDataLoop_To``
    Description: Last rotation-data file index used by the inner RotData loop.
    Default ``1``

``RotData_User_Selection``
    Description: Optional explicit list of selected rotation-data file indices.

    If this option is present, it selects the exact ``RotData_#.csv`` files
    used by the inner RotData loop, for example ``{1, 3, 4}``. If it is not
    present, the contiguous range from ``RotDataLoop_From`` to
    ``RotDataLoop_To`` is used.

    For a fixed data-file index and a selected rotation-data index, output is
    written to a nested folder such as ``SANS_1/RotData_3/``. If
    ``StructDataLoop`` is active as well, the output is nested as
    ``SANS_1/StructData_2/RotData_3/``.

Constant Parameters
-------------------

``XYZ_Unit_Factor``
    Conversion factors for coordinate units. Scaling factor applied to spatial coordinates.
    This parameter is equal to ``1``, if the positional input data is in units of :math:`\mathrm{nm}`.
    Default value ``1``

``Scattering_Volume_V``
    Effective scattering volume in units of :math:`\mathrm{m}^3`
    Default value ``2.618e-24``

``RotMat_alpha``
    Global sample rotation angle in degree.
    Default value ``0.0``

``RotMat_beta``
    Global sample rotation angle in degree. Rotates the sample in the
    :math:`x-y`-plane.
    Default value ``0.0``

    These two angles define a global rotation of the complete imported system.
    They are independent of ``RotData``, which defines object-wise local
    rotations.

``Polarization``
    Polarization vector (Px, Py, Pz). Defines the polarization vector of the incoming neutron beam.
    Default ``(0,0,1)``

``Number_Of_q_Points``
    Number of sampled q values in Fourier space.
    Default value ``1000``

``Number_Of_theta_Points``
    Number of angular sampling points in Fourier space.
    Default value ``1000``

``Number_Of_r_Points``
    Number of sampled q values in real space (correlation).
    Default value ``1000``

``Number_Of_alpha_Points``
    Number of angular sampling points in real space (correlation).
    Default value ``1000``

``q_max``
    Maximum scattering vector in Fourier space in units of :math:`\mathrm{nm}^{-1}`
    Default value ``3.0``

``r_max``
    Maximum real-space correlation distance in units of :math:`\mathrm{nm}`.
    Default value ``15.0 ``


Micromagnetic Parameters
------------------------

Parameters describing the magnetic properties of discretized simulation
cells.

``Cell_Nuclear_SLD``
    Description: Nuclear scattering length density of the cell
    Default ``8e14``

``Cell_Magnetization``
    Description: Magnetization inside the discretized cell
    Default ``486e3``

``Cuboid_Cell_Size``
    Description: Dimensions of the discretization cell (x,y,z)
    Default ``(2,2,2)``

Output Options
--------------

2D SANS cross sections
^^^^^^^^^^^^^^^^^^^^^^

In NuMagSANS, the incoming neutron beam is defined to propagate along the
:math:`x`-axis, such that the corresponding wave vector satisfies
:math:`\mathbf{k}_0 \parallel \mathbf{e}_x`. The two-dimensional detector
plane is oriented perpendicular to the incident beam direction. The
scattering vector :math:`\mathbf{q}` on the detector can therefore be
written as

.. math::

    \mathbf{q} =
    \begin{bmatrix}
        q_x \\
        q_y \\
        q_z
    \end{bmatrix}
    =
    \begin{bmatrix}
        0 \\
        q \sin\theta \\
        q \cos\theta
    \end{bmatrix}.

The magnetic SANS amplitudes are expressed through the Halpern-Johnson vector
:math:`\widetilde{\mathbf{Q}}`, i.e. the component of the Fourier-transformed
magnetization that contributes to magnetic neutron scattering:

.. math::

    \widetilde{\mathbf{Q}}
    =
    \hat{\mathbf{q}}
    \left(\hat{\mathbf{q}}\cdot\widetilde{\mathbf{M}}\right)
    -
    \widetilde{\mathbf{M}},

with

.. math::

    \hat{\mathbf{q}}
    =
    \begin{bmatrix}
        0 \\
        \sin\theta \\
        \cos\theta
    \end{bmatrix},
    \qquad
    \widetilde{\mathbf{M}}
    =
    \begin{bmatrix}
        \widetilde{M}_x \\
        \widetilde{M}_y \\
        \widetilde{M}_z
    \end{bmatrix}.

For the detector geometry used in NuMagSANS this gives

.. math::

    \widetilde{\mathbf{Q}}
    =
    \begin{bmatrix}
        -\widetilde{M}_x \\
        \left(\widetilde{M}_z\sin\theta-\widetilde{M}_y\cos\theta\right)\cos\theta \\
        \left(\widetilde{M}_y\cos\theta-\widetilde{M}_z\sin\theta\right)\sin\theta
    \end{bmatrix}.

The same expression is applied to the real and imaginary parts of
:math:`\widetilde{\mathbf{M}}` in the atomistic kernels.

In NuMagSANS the 2D SANS cross sections are exported in units of :math:`\mathrm{cm}^{-1}`.

``Nuclear_2D``
    Nuclear SANS cross section on the 2D detector.

    .. math::

        \frac{d\Sigma_{\mathrm{N}}}{d\Omega}(q,\theta)
        =
        \frac{8\pi^3}{V}\, |\widetilde{N}|^2

``Unpolarized_2D``
    Unpolarized magnetic SANS cross section on the 2D detector.

    .. math::

        \frac{d\Sigma_{\mathrm{M}}}{d\Omega}(q,\theta)
        =
        \frac{8\pi^3}{V}\, b_{\mathrm{H}}^2 |\widetilde{\mathbf{Q}}|^2

``Polarized_2D``
    Auxiliary polarized magnetic SANS cross section on the 2D detector.

    .. math::

        \frac{d\Sigma_{\mathrm{P}}}{d\Omega}(q,\theta)
        =
        \frac{8\pi^3}{V}\, b_{\mathrm{H}}^2
        \left|\hat{\mathbf{P}}\cdot\widetilde{\mathbf{Q}}\right|^2

``NuclearMagnetic_2D``
    Nuclear–magnetic SANS cross section on the 2D detector.

    .. math::

        \frac{d\Sigma_{\mathrm{NM}}}{d\Omega}(q,\theta)
        =
        \frac{8\pi^3}{V}\, b_{\mathrm{H}}
        \hat{\mathbf{P}}\cdot\left[\widetilde{N}\widetilde{\mathbf{Q}}^{\ast}+\widetilde{N}^{\ast}\widetilde{\mathbf{Q}}\right]

``Chiral_2D``
    Chiral SANS cross section on the 2D detector.

    .. math::

        \frac{d\Sigma_{\chi}}{d\Omega}(q,\theta)
        =
        -i\frac{8\pi^3}{V}\, b_{\mathrm{H}}^2
        \hat{\mathbf{P}}\cdot
        \left(\widetilde{\mathbf{Q}}\times\widetilde{\mathbf{Q}}^{\ast}\right)

``SpinFlip_2D``
    Spin-flip SANS cross section on the 2D detector.

    .. math::

        \frac{d\Sigma_{\mathrm{sf}}}{d\Omega}(q,\theta)
        =
        \frac{d\Sigma_{\mathrm{M}}}{d\Omega}(q,\theta)
        -
        \frac{d\Sigma_{\mathrm{P}}}{d\Omega}(q,\theta)

``PM_SpinFlip_2D``
    Spin-flip SANS cross section on the 2D detector.

    .. math::

        \frac{d\Sigma^{+-}}{d\Omega}(q,\theta)
        =
        \frac{d\Sigma_{\mathrm{sf}}}{d\Omega}(q,\theta)
        +
        \frac{d\Sigma_{\chi}}{d\Omega}(q,\theta)

``MP_SpinFlip_2D``
    Spin-flip SANS cross section on the 2D detector.

    .. math::

        \frac{d\Sigma^{+-}}{d\Omega}(q,\theta)
        =
        \frac{d\Sigma_{\mathrm{sf}}}{d\Omega}(q,\theta)
        -
        \frac{d\Sigma_{\chi}}{d\Omega}(q,\theta)


``PP_NonSpinFlip_2D``
    Non-spin-flip SANS cross section on the 2D detector.

    .. math::

        \frac{d\Sigma^{++}}{d\Omega}(q,\theta)
        =
        \frac{d\Sigma_{\mathrm{N}}}{d\Omega}(q,\theta)
        +
        \frac{d\Sigma_{\mathrm{NM}}}{d\Omega}(q,\theta)
        +
        \frac{d\Sigma_{\mathrm{P}}}{d\Omega}(q,\theta)

``MM_NonSpinFlip_2D``
    Non-spin-flip SANS cross section on the 2D detector.

    .. math::

        \frac{d\Sigma^{--}}{d\Omega}(q,\theta)
        =
        \frac{d\Sigma_{\mathrm{N}}}{d\Omega}(q,\theta)
        -
        \frac{d\Sigma_{\mathrm{NM}}}{d\Omega}(q,\theta)
        +
        \frac{d\Sigma_{\mathrm{P}}}{d\Omega}(q,\theta)

``P_SANSPOL_2D``
    Non-Spin-flip SANS cross section on the 2D detector.

    .. math::

        \frac{d\Sigma^{+}}{d\Omega}(q,\theta)
        =
        \frac{d\Sigma^{++}}{d\Omega}(q,\theta)
        +
        \frac{d\Sigma^{+-}}{d\Omega}(q,\theta)

``M_SANSPOL_2D``
    Non-Spin-flip SANS cross section on the 2D detector.

    .. math::

        \frac{d\Sigma^{-}}{d\Omega}(q,\theta)
        =
        \frac{d\Sigma^{--}}{d\Omega}(q,\theta)
        +
        \frac{d\Sigma^{-+}}{d\Omega}(q,\theta)

2D correlation functions
^^^^^^^^^^^^^^^^^^^^^^^^

NuMagSANS can compute two-dimensional real-space correlation functions from
the detector-plane SANS cross sections. The real-space vector in the detector
plane is parameterized as

.. math::

    \mathbf r =
    \begin{bmatrix}
        0 \\
        r\sin\alpha \\
        r\cos\alpha
    \end{bmatrix},

while the scattering vector in reciprocal space is written as

.. math::

    \mathbf q =
    \begin{bmatrix}
        0 \\
        q\sin\theta \\
        q\cos\theta
    \end{bmatrix}.

Using these parametrizations the scalar product becomes

.. math::

    \mathbf q\cdot\mathbf r = qr\cos(\theta-\alpha).

The two-dimensional correlation functions are obtained from the SANS cross
sections via the cosine transform

.. math::

    C(r,\alpha)
    =
    \int_{0}^{\infty}\int_{0}^{2\pi}
    \frac{d\Sigma}{d\Omega}(q,\theta)
    \cos\left(qr\cos(\theta-\alpha)\right)
    q\,dq\,d\theta .

The following two-dimensional correlation functions can be calculated.

``Nuclear_Corr_2D``
    2D correlation function of the nuclear SANS cross section.

    .. math::

        C_{\mathrm N}(r,\alpha)
        =
        \int_{0}^{\infty}\int_{0}^{2\pi}
        \frac{d\Sigma_{\mathrm N}}{d\Omega}(q,\theta)
        \cos\left(qr\cos(\theta-\alpha)\right)
        q\,dq\,d\theta


``Unpolarized_Corr_2D``
    2D correlation function of the unpolarized magnetic SANS cross section.

    .. math::

        C_{\mathrm M}(r,\alpha)
        =
        \int_{0}^{\infty}\int_{0}^{2\pi}
        \frac{d\Sigma_{\mathrm M}}{d\Omega}(q,\theta)
        \cos\left(qr\cos(\theta-\alpha)\right)
        q\,dq\,d\theta


``Polarized_Corr_2D``
    2D correlation function of the polarized magnetic SANS cross section.

    .. math::

        C_{\mathrm P}(r,\alpha)
        =
        \int_{0}^{\infty}\int_{0}^{2\pi}
        \frac{d\Sigma_{\mathrm P}}{d\Omega}(q,\theta)
        \cos\left(qr\cos(\theta-\alpha)\right)
        q\,dq\,d\theta


``NuclearMagnetic_Corr_2D``
    2D correlation function of the nuclear–magnetic interference SANS cross section.

    .. math::

        C_{\mathrm{NM}}(r,\alpha)
        =
        \int_{0}^{\infty}\int_{0}^{2\pi}
        \frac{d\Sigma_{\mathrm{NM}}}{d\Omega}(q,\theta)
        \cos\left(qr\cos(\theta-\alpha)\right)
        q\,dq\,d\theta


``SpinFlip_Corr_2D``
    2D correlation function of the spin-flip SANS cross section.

    .. math::

        C_{\mathrm{sf}}(r,\alpha)
        =
        \int_{0}^{\infty}\int_{0}^{2\pi}
        \frac{d\Sigma_{\mathrm{sf}}}{d\Omega}(q,\theta)
        \cos\left(qr\cos(\theta-\alpha)\right)
        q\,dq\,d\theta


``Chiral_Corr_2D``
    2D correlation function of the chiral SANS cross section.

    .. math::

        C_{\chi}(r,\alpha)
        =
        \int_{0}^{\infty}\int_{0}^{2\pi}
        \frac{d\Sigma_{\chi}}{d\Omega}(q,\theta)
        \cos\left(qr\cos(\theta-\alpha)\right)
        q\,dq\,d\theta


``PM_SpinFlip_Corr_2D``
    2D correlation function of the :math:`(+,-)` spin-flip channel.

    .. math::

        C^{+-}(r,\alpha)
        =
        \int_{0}^{\infty}\int_{0}^{2\pi}
        \frac{d\Sigma^{+-}}{d\Omega}(q,\theta)
        \cos\left(qr\cos(\theta-\alpha)\right)
        q\,dq\,d\theta


``MP_SpinFlip_Corr_2D``
    2D correlation function of the :math:`(-,+)` spin-flip channel.

    .. math::

        C^{-+}(r,\alpha)
        =
        \int_{0}^{\infty}\int_{0}^{2\pi}
        \frac{d\Sigma^{-+}}{d\Omega}(q,\theta)
        \cos\left(qr\cos(\theta-\alpha)\right)
        q\,dq\,d\theta


``PP_NonSpinFlip_Corr_2D``
    2D correlation function of the :math:`(+,+)` non-spin-flip channel.

    .. math::

        C^{++}(r,\alpha)
        =
        \int_{0}^{\infty}\int_{0}^{2\pi}
        \frac{d\Sigma^{++}}{d\Omega}(q,\theta)
        \cos\left(qr\cos(\theta-\alpha)\right)
        q\,dq\,d\theta


``MM_NonSpinFlip_Corr_2D``
    2D correlation function of the :math:`(-,-)` non-spin-flip channel.

    .. math::

        C^{--}(r,\alpha)
        =
        \int_{0}^{\infty}\int_{0}^{2\pi}
        \frac{d\Sigma^{--}}{d\Omega}(q,\theta)
        \cos\left(qr\cos(\theta-\alpha)\right)
        q\,dq\,d\theta


``P_SANSPOL_Corr_2D``
    2D correlation function of the SANSPOL :math:`(+ )` channel.

    .. math::

        C^{+}(r,\alpha)
        =
        \int_{0}^{\infty}\int_{0}^{2\pi}
        \frac{d\Sigma^{+}}{d\Omega}(q,\theta)
        \cos\left(qr\cos(\theta-\alpha)\right)
        q\,dq\,d\theta


``M_SANSPOL_Corr_2D``
    2D correlation function of the SANSPOL :math:`(- )` channel.

    .. math::

        C^{-}(r,\alpha)
        =
        \int_{0}^{\infty}\int_{0}^{2\pi}
        \frac{d\Sigma^{-}}{d\Omega}(q,\theta)
        \cos\left(qr\cos(\theta-\alpha)\right)
        q\,dq\,d\theta

1D Azimuthally averaged SANS cross sections
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In NuMagSANS the 1D SANS cross sections are exported in units of :math:`\mathrm{cm}^{-1}`.

``Nuclear_1D``
    Azimuthally averaged nuclear SANS cross section.

    .. math::

        I_{\mathrm{N}}(q)
        =
        \frac{1}{2\pi}
        \int_{0}^{2\pi}
        \frac{d\Sigma_{\mathrm{N}}}{d\Omega}(q,\theta)
        \, d\theta

``Unpolarized_1D``
    Azimuthally averaged unpolarized magnetic SANS cross section.

    .. math::

        I_{\mathrm{M}}(q)
        =
        \frac{1}{2\pi}
        \int_{0}^{2\pi}
        \frac{d\Sigma_{\mathrm{M}}}{d\Omega}(q,\theta)
        \, d\theta

``Polarized_1D``
    Azimuthally averaged polarized magnetic SANS cross section.

    .. math::

        I_{\mathrm{P}}(q)
        =
        \frac{1}{2\pi}
        \int_{0}^{2\pi}
        \frac{d\Sigma_{\mathrm{P}}}{d\Omega}(q,\theta)
        \, d\theta

``NuclearMagnetic_1D``
    Azimuthally averaged nuclear–magnetic interference cross section.

    .. math::

        I_{\mathrm{NM}}(q)
        =
        \frac{1}{2\pi}
        \int_{0}^{2\pi}
        \frac{d\Sigma_{\mathrm{NM}}}{d\Omega}(q,\theta)
        \, d\theta

``Chiral_1D``
    Azimuthally averaged chiral SANS cross section.

    .. math::

        I_{\chi}(q)
        =
        \frac{1}{2\pi}
        \int_{0}^{2\pi}
        \frac{d\Sigma_{\chi}}{d\Omega}(q,\theta)
        \, d\theta

``SpinFlip_1D``
    Azimuthally averaged spin-flip SANS cross section.

    .. math::

        I_{\mathrm{sf}}(q)
        =
        \frac{1}{2\pi}
        \int_{0}^{2\pi}
        \frac{d\Sigma_{\mathrm{sf}}}{d\Omega}(q,\theta)
        \, d\theta

``PM_SpinFlip_1D``
    Azimuthally averaged spin-flip SANS cross section for the
    :math:`(+,-)` polarization channel.

    .. math::

        I^{+-}(q)
        =
        \frac{1}{2\pi}
        \int_{0}^{2\pi}
        \frac{d\Sigma^{+-}}{d\Omega}(q,\theta)
        \, d\theta

``MP_SpinFlip_1D``
    Azimuthally averaged spin-flip SANS cross section for the
    :math:`(-,+)` polarization channel.

    .. math::

        I^{-+}(q)
        =
        \frac{1}{2\pi}
        \int_{0}^{2\pi}
        \frac{d\Sigma^{-+}}{d\Omega}(q,\theta)
        \, d\theta

``PP_NonSpinFlip_1D``
    Azimuthally averaged non-spin-flip SANS cross section for the
    :math:`(+,+)` polarization channel.

    .. math::

        I^{++}(q)
        =
        \frac{1}{2\pi}
        \int_{0}^{2\pi}
        \frac{d\Sigma^{++}}{d\Omega}(q,\theta)
        \, d\theta

``MM_NonSpinFlip_1D``
    Azimuthally averaged non-spin-flip SANS cross section for the
    :math:`(-,-)` polarization channel.

    .. math::

        I^{--}(q)
        =
        \frac{1}{2\pi}
        \int_{0}^{2\pi}
        \frac{d\Sigma^{--}}{d\Omega}(q,\theta)
        \, d\theta

``P_SANSPOL_1D``
    Azimuthally averaged SANSPOL cross section for positive neutron
    polarization.

    .. math::

        I^{+}(q)
        =
        \frac{1}{2\pi}
        \int_{0}^{2\pi}
        \frac{d\Sigma^{+}}{d\Omega}(q,\theta)
        \, d\theta

``M_SANSPOL_1D``
    Azimuthally averaged SANSPOL cross section for negative neutron
    polarization.

    .. math::

        I^{-}(q)
        =
        \frac{1}{2\pi}
        \int_{0}^{2\pi}
        \frac{d\Sigma^{-}}{d\Omega}(q,\theta)
        \, d\theta

1D Correlation functions
^^^^^^^^^^^^^^^^^^^^^^^^

In addition to the azimuthally averaged SANS cross sections :math:`I(q)`,
NuMagSANS can compute the corresponding real-space correlation functions
:math:`c(r)`. These functions are obtained from the scattering intensity
via a spherical Bessel transform

.. math::

    c(r)
    =
    \int_{0}^{\infty}
    I(q)\, j_0(qr)\, q^2 \, dq ,

where :math:`j_0(x)` denotes the spherical Bessel function of order zero,

.. math::

    j_0(x) = \frac{\sin x}{x}.

In NuMagSANS the 1D correlation functions are exported in units of :math:`\mathrm{nm}^{-4}`.
The following correlation functions can be calculated.

``Nuclear_Corr_1D``
    Correlation function of the nuclear SANS cross section.

    .. math::

        c_{\mathrm{N}}(r)
        =
        \int_{0}^{\infty}
        I_{\mathrm{N}}(q)\, j_0(qr)\, q^2 \, dq

``Unpolarized_Corr_1D``
    Correlation function of the unpolarized magnetic SANS cross section.

    .. math::

        c_{\mathrm{M}}(r)
        =
        \int_{0}^{\infty}
        I_{\mathrm{M}}(q)\, j_0(qr)\, q^2 \, dq

``Polarized_Corr_1D``
    Correlation function of the polarized magnetic SANS cross section.

    .. math::

        c_{\mathrm{P}}(r)
        =
        \int_{0}^{\infty}
        I_{\mathrm{P}}(q)\, j_0(qr)\, q^2 \, dq

``NuclearMagnetic_Corr_1D``
    Correlation function of the nuclear–magnetic interference SANS cross section.

    .. math::

        c_{\mathrm{NM}}(r)
        =
        \int_{0}^{\infty}
        I_{\mathrm{NM}}(q)\, j_0(qr)\, q^2 \, dq

``Chiral_Corr_1D``
    Correlation function of the chiral SANS cross section.

    .. math::

        c_{\chi}(r)
        =
        \int_{0}^{\infty}
        I_{\chi}(q)\, j_0(qr)\, q^2 \, dq

``SpinFlip_Corr_1D``
    Correlation function of the spin-flip SANS cross section.

    .. math::

        c_{\mathrm{sf}}(r)
        =
        \int_{0}^{\infty}
        I_{\mathrm{sf}}(q)\, j_0(qr)\, q^2 \, dq

``PM_SpinFlip_Corr_1D``
    Correlation function of the :math:`(+,-)` spin-flip channel.

    .. math::

        c^{+-}(r)
        =
        \int_{0}^{\infty}
        I^{+-}(q)\, j_0(qr)\, q^2 \, dq

``MP_SpinFlip_Corr_1D``
    Correlation function of the :math:`(-,+)` spin-flip channel.

    .. math::

        c^{-+}(r)
        =
        \int_{0}^{\infty}
        I^{-+}(q)\, j_0(qr)\, q^2 \, dq

``PP_NonSpinFlip_Corr_1D``
    Correlation function of the :math:`(+,+)` non-spin-flip channel.

    .. math::

        c^{++}(r)
        =
        \int_{0}^{\infty}
        I^{++}(q)\, j_0(qr)\, q^2 \, dq

``MM_NonSpinFlip_Corr_1D``
    Correlation function of the :math:`(-,-)` non-spin-flip channel.

    .. math::

        c^{--}(r)
        =
        \int_{0}^{\infty}
        I^{--}(q)\, j_0(qr)\, q^2 \, dq

``P_SANSPOL_Corr_1D``
    Correlation function of the SANSPOL :math:`(+ )` channel.

    .. math::

        c^{+}(r)
        =
        \int_{0}^{\infty}
        I^{+}(q)\, j_0(qr)\, q^2 \, dq

``M_SANSPOL_Corr_1D``
    Correlation function of the SANSPOL :math:`(- )` channel.

    .. math::

        c^{-}(r)
        =
        \int_{0}^{\infty}
        I^{-}(q)\, j_0(qr)\, q^2 \, dq

1D Pair-distance distribution functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In addition to the correlation functions :math:`c(r)`, NuMagSANS can compute
the corresponding pair-distance distribution functions :math:`p(r)`. These
functions are obtained from the correlation functions via

.. math::

    p(r) = r^2 c(r).

Using the definition of :math:`c(r)`, the pair-distance distribution function
can also be written as

.. math::

    p(r)
    =
    r^2
    \int_{0}^{\infty}
    I(q)\, j_0(qr)\, q^2 \, dq .

In NuMagSANS the 1D pair-distance distribution functions are exported in units of :math:`\mathrm{nm}^{-2}`.
The following pair-distance distribution functions can be calculated.

``Nuclear_PairDist_1D``
    Pair-distance distribution function of the nuclear SANS cross section.

    .. math::

        p_{\mathrm{N}}(r)
        =
        r^2
        \int_{0}^{\infty}
        I_{\mathrm{N}}(q)\, j_0(qr)\, q^2 \, dq

``Unpolarized_PairDist_1D``
    Pair-distance distribution function of the unpolarized magnetic SANS cross section.

    .. math::

        p_{\mathrm{M}}(r)
        =
        r^2
        \int_{0}^{\infty}
        I_{\mathrm{M}}(q)\, j_0(qr)\, q^2 \, dq

``Polarized_PairDist_1D``
    Pair-distance distribution function of the polarized magnetic SANS cross section.

    .. math::

        p_{\mathrm{P}}(r)
        =
        r^2
        \int_{0}^{\infty}
        I_{\mathrm{P}}(q)\, j_0(qr)\, q^2 \, dq

``NuclearMagnetic_PairDist_1D``
    Pair-distance distribution function of the nuclear–magnetic interference SANS cross section.

    .. math::

        p_{\mathrm{NM}}(r)
        =
        r^2
        \int_{0}^{\infty}
        I_{\mathrm{NM}}(q)\, j_0(qr)\, q^2 \, dq

``Chiral_PairDist_1D``
    Pair-distance distribution function of the chiral SANS cross section.

    .. math::

        p_{\chi}(r)
        =
        r^2
        \int_{0}^{\infty}
        I_{\chi}(q)\, j_0(qr)\, q^2 \, dq

``SpinFlip_PairDist_1D``
    Pair-distance distribution function of the spin-flip SANS cross section.

    .. math::

        p_{\mathrm{sf}}(r)
        =
        r^2
        \int_{0}^{\infty}
        I_{\mathrm{sf}}(q)\, j_0(qr)\, q^2 \, dq

``PM_SpinFlip_PairDist_1D``
    Pair-distance distribution function of the :math:`(+,-)` spin-flip channel.

    .. math::

        p^{+-}(r)
        =
        r^2
        \int_{0}^{\infty}
        I^{+-}(q)\, j_0(qr)\, q^2 \, dq

``MP_SpinFlip_PairDist_1D``
    Pair-distance distribution function of the :math:`(-,+)` spin-flip channel.

    .. math::

        p^{-+}(r)
        =
        r^2
        \int_{0}^{\infty}
        I^{-+}(q)\, j_0(qr)\, q^2 \, dq

``PP_NonSpinFlip_PairDist_1D``
    Pair-distance distribution function of the :math:`(+,+)` non-spin-flip channel.

    .. math::

        p^{++}(r)
        =
        r^2
        \int_{0}^{\infty}
        I^{++}(q)\, j_0(qr)\, q^2 \, dq

``MM_NonSpinFlip_PairDist_1D``
    Pair-distance distribution function of the :math:`(-,-)` non-spin-flip channel.

    .. math::

        p^{--}(r)
        =
        r^2
        \int_{0}^{\infty}
        I^{--}(q)\, j_0(qr)\, q^2 \, dq

``P_SANSPOL_PairDist_1D``
    Pair-distance distribution function of the SANSPOL :math:`(+ )` channel.

    .. math::

        p^{+}(r)
        =
        r^2
        \int_{0}^{\infty}
        I^{+}(q)\, j_0(qr)\, q^2 \, dq

``M_SANSPOL_PairDist_1D``
    Pair-distance distribution function of the SANSPOL :math:`(- )` channel.

    .. math::

        p^{-}(r)
        =
        r^2
        \int_{0}^{\infty}
        I^{-}(q)\, j_0(qr)\, q^2 \, dq


Angular Spectrum
----------------

Settings for angular spectrum calculations.

``k_max``
    Maximum angular wave number.
    Default value ``10``

``Angular_Spec``
    Enable angular spectrum calculation. (binary selection)
    Default value ``0``
