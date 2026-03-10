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


Data Paths
----------

These parameters define where the simulation reads input data and where
results are written.

+---------------------+------------------------------------------------------------+------------------------------+
| Parameter           | Description                                                | Default                      |
+=====================+============================================================+==============================+
| ``NucDataPath``     | Path to nuclear real-space scattering data                 | ``RealSpaceData/NucData``    |
+---------------------+------------------------------------------------------------+------------------------------+
| ``MagDataPath``     | Path to magnetic real-space moment data                    | ``RealSpaceData/MagData``    |
+---------------------+------------------------------------------------------------+------------------------------+
| ``StructDataFilename`` | Structural data file describing atomic positions       | ``RealSpaceData/StructData.csv`` |
+---------------------+------------------------------------------------------------+------------------------------+
| ``foldernameSANSData`` | Output directory where scattering data is written      | ``NuMagSANS_Output``         |
+---------------------+------------------------------------------------------------+------------------------------+


Data Selection
--------------

These flags control which datasets are used by the simulation.

+-------------------------+-----------------------------------------------------------+---------+
| Parameter               | Description                                               | Default |
+=========================+===========================================================+=========+
| ``NucData_activate``    | Activate nuclear scattering data                          | 0       |
+-------------------------+-----------------------------------------------------------+---------+
| ``MagData_activate``    | Activate magnetic moment data                             | 0       |
+-------------------------+-----------------------------------------------------------+---------+
| ``StructData_activate`` | Activate structural data input                            | 0       |
+-------------------------+-----------------------------------------------------------+---------+
| ``Exclude_Zero_Moments``| Ignore atoms or cells with zero magnetic moment           | 0       |
+-------------------------+-----------------------------------------------------------+---------+


Fourier Approach
----------------

Defines how the Fourier transform from real space to reciprocal space
is performed.

+-------------------+----------------------------------------------------+-----------+
| Parameter         | Description                                        | Default   |
+===================+====================================================+===========+
| ``Fourier_Approach`` | Method used for Fourier transformation          | ``atomistic`` |
+-------------------+----------------------------------------------------+-----------+


Loop Control
------------

These parameters control batch simulations or repeated calculations.

+-------------------+---------------------------------------------------+---------+
| Parameter         | Description                                       | Default |
+===================+===================================================+=========+
| ``Loop_Modus``    | Enable looping over simulation indices            | 0       |
+-------------------+---------------------------------------------------+---------+
| ``Loop_From``     | First index used in loop simulations              | 1       |
+-------------------+---------------------------------------------------+---------+
| ``Loop_To``       | Last index used in loop simulations               | 20      |
+-------------------+---------------------------------------------------+---------+
| ``User_Selection``| Explicit list of selected indices                 | [1]     |
+-------------------+---------------------------------------------------+---------+


Units
-----

Conversion factors for coordinate units.

+-------------------+---------------------------------------------------+---------+
| Parameter         | Description                                       | Default |
+===================+===================================================+=========+
| ``XYZ_Unit_Factor`` | Scaling factor applied to spatial coordinates    | 1       |
+-------------------+---------------------------------------------------+---------+


Micromagnetic Parameters
------------------------

Parameters describing the magnetic properties of discretized simulation
cells.

+----------------------+----------------------------------------------------+-----------+
| Parameter            | Description                                        | Default   |
+======================+====================================================+===========+
| ``Cell_Nuclear_SLD`` | Nuclear scattering length density of the cell      | 8e14      |
+----------------------+----------------------------------------------------+-----------+
| ``Cell_Magnetization`` | Magnetization inside the discretized cell        | 486e3     |
+----------------------+----------------------------------------------------+-----------+
| ``Cuboid_Cell_Size`` | Dimensions of the discretization cell (x,y,z)      | (2,2,2)   |
+----------------------+----------------------------------------------------+-----------+


Scattering Volume
-----------------

Defines the physical scattering volume used for normalization.

+----------------------+----------------------------------------------------+-----------+
| Parameter            | Description                                        | Default   |
+======================+====================================================+===========+
| ``Scattering_Volume_V`` | Effective scattering volume                     | 2.618e-24 |
+----------------------+----------------------------------------------------+-----------+


Rotation of the Sample
----------------------

These parameters define the orientation of the sample relative to the
scattering geometry.

+-------------------+---------------------------------------------------+---------+
| Parameter         | Description                                       | Default |
+===================+===================================================+=========+
| ``RotMat_alpha``  | Rotation angle around x-axis                      | 0.0     |
+-------------------+---------------------------------------------------+---------+
| ``RotMat_beta``   | Rotation angle around y-axis                      | 0.0     |
+-------------------+---------------------------------------------------+---------+


Neutron Polarization
--------------------

Defines the polarization vector of the incoming neutron beam.

+-------------------+---------------------------------------------------+-----------+
| Parameter         | Description                                       | Default   |
+===================+===================================================+===========+
| ``Polarization``  | Polarization vector (Px, Py, Pz)                  | (0,0,1)   |
+-------------------+---------------------------------------------------+-----------+


Reciprocal Space Sampling
-------------------------

Controls the resolution of the scattering calculation.

+-----------------------------+--------------------------------------------+---------+
| Parameter                   | Description                                | Default |
+=============================+============================================+=========+
| ``Number_Of_q_Points``      | Number of sampled q values                 | 1000    |
+-----------------------------+--------------------------------------------+---------+
| ``Number_Of_theta_Points``  | Number of angular sampling points          | 1000    |
+-----------------------------+--------------------------------------------+---------+
| ``Number_Of_r_Points``      | Sampling points in real-space correlations | 1000    |
+-----------------------------+--------------------------------------------+---------+
| ``Number_Of_alpha_Points``  | Angular resolution for correlation analysis| 1000    |
+-----------------------------+--------------------------------------------+---------+
| ``q_max``                   | Maximum scattering vector                  | 3.0     |
+-----------------------------+--------------------------------------------+---------+
| ``r_max``                   | Maximum real-space distance                | 15.0    |
+-----------------------------+--------------------------------------------+---------+


Angular Spectrum
----------------

Settings for angular spectrum calculations.

+---------------+----------------------------------------------------+---------+
| Parameter     | Description                                        | Default |
+===============+====================================================+=========+
| ``k_max``     | Maximum angular wave number                        | 10      |
+---------------+----------------------------------------------------+---------+
| ``Angular_Spec`` | Enable angular spectrum calculation              | 0       |
+---------------+----------------------------------------------------+---------+


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

Azimuthally averaged SANS cross sections
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

Pair-distance distribution functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
