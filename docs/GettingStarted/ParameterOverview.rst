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

   # Create NuMagSANS facade object
   facade = NuMagSANS()

   # Temporary configuration file
   config = BASE_DIR / "NuMagSANSInput_temp.conf"

   # Write configuration file
   facade.write_config(
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
   facade.run(config)

   # Remove temporary configuration file
   facade.config_clear(config)


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


Output Selection
----------------

Simulation outputs are controlled through the ``enable_outputs`` argument.

Example:

.. code-block:: python

    sim.write_config(
        "simulation.conf",
        enable_outputs=[
            "Unpolarized_2D",
            "SpinFlip_2D"
        ]
    )

Each output corresponds to a specific scattering quantity calculated by
the backend. Internally these outputs are written to the output folder
defined by ``foldernameSANSData``.

Typical examples include:

+---------------------+---------------------------------------------+
| Output              | Description                                 |
+=====================+=============================================+
| ``Nuclear_2D``      | Nuclear scattering intensity (2D detector)  |
+---------------------+---------------------------------------------+
| ``Unpolarized_2D``  | Total unpolarized scattering cross section  |
+---------------------+---------------------------------------------+
| ``SpinFlip_2D``     | Spin-flip magnetic scattering               |
+---------------------+---------------------------------------------+
| ``Nuclear_1D``      | Radially averaged nuclear scattering        |
+---------------------+---------------------------------------------+
| ``SpinFlip_1D``     | Radially averaged spin-flip scattering      |
+---------------------+---------------------------------------------+


Generated Configuration File
----------------------------

The Python interface generates a configuration file containing the
parameters listed above.

Example snippet:

.. code-block:: text

    Cell_Magnetization = 486000;
    q_max = 3.0;
    Number_Of_q_Points = 1000;

This file is then passed to the NuMagSANS backend executable during
simulation execution.
