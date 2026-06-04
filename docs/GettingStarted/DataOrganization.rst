Data Organization
=================

The input data for NuMagSANS consist of discretized real-space datasets
that can be organized in several ways. At the highest level, NuMagSANS
expects a ``RealSpaceData`` directory.

Inside this directory, the input data are separated into four data layers:
``MagData``, ``NucData``, ``StructData``, and ``RotData``. ``MagData`` and
``NucData`` contain local real-space object data, while ``StructData`` and
``RotData`` provide optional assembly metadata for object translations and
object-wise rotations.

The possible combinations of these data layers are summarized in the table
below and are described in more detail in the following scenarios.

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
   * - ``MagData + StructData``
     - Magnetic local objects translated by object-center positions.
   * - ``MagData + RotData``
     - Magnetic local objects rotated individually.
   * - ``MagData + StructData + RotData``
     - Magnetic local objects rotated and translated as an assembly.
   * - ``NucData + StructData``
     - Nuclear local objects translated by object-center positions.
   * - ``NucData + RotData``
     - Nuclear local objects rotated individually.
   * - ``NucData + StructData + RotData``
     - Nuclear local objects rotated and translated as an assembly.
   * - ``MagData + NucData + StructData + RotData``
     - Combined nuclear-magnetic assembly with object-wise translations and rotations.


Simulation scenarios and their data organization
------------------------------------------------

NuMagSANS supports several simulation scenarios, and the corresponding data
handling can be difficult to understand all at once. Therefore, the input-data
organization is introduced scenario by scenario, starting with the simplest
setup.


Scenario 1: Single object with a single magnetic state
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The simplest NuMagSANS scenario is a single magnetic object with a single
magnetic state. This could be, for example, a spherical magnetic nanoparticle
with a uniform magnetization state.

For this scenario, the input dataset is organized as follows. At the highest
level, the ``RealSpaceData`` directory contains a ``MagData`` directory. Inside
``MagData``, the local object data are stored in the ``Object_1`` directory,
which contains a single ``m_1.csv`` file.

It is important that the object directory name ends with an underscore ``_``
followed by an enumeration starting from ``1``. The prefix before the
underscore can be chosen by the user. For example, ``Object_1``,
``Particle_1``, ``Angle_1``, or ``Orientation_1`` are possible names as long as
the numbering convention is preserved.

The same convention applies to the magnetic data file ``m_1.csv``. The prefix
``m`` identifies magnetic data, while the index ``1`` is the data-file index
selected through ``User_Selection`` or ``Loop_Modus``.

.. code-block:: text

   RealSpaceData/
     MagData/
       Object_1/
         m_1.csv

The ``m_1.csv`` file is a header-free, whitespace-separated file with six
columns. The first three columns contain the Cartesian position information
``x``, ``y``, ``z``. The last three columns contain the magnetic moment or
magnetization vector components ``mx``, ``my``, ``mz``.

.. code-block:: text

   x y z mx my mz

For atomistic systems, the coordinates are commonly given in nanometers and
the atomic magnetic moments in units of the Bohr magneton. The relevant minimal
configuration is:

.. code-block:: conf

   MagData_activate = 1;
   NucData_activate = 0;
   StructData_activate = 0;
   RotData_activate = 0;
   MagDataPath = RealSpaceData/MagData;
   User_Selection = {1};

In this scenario, no ``NucData``, ``StructData``, or ``RotData`` are required.
The object position, shape, and magnetic state are fully contained in the local
coordinates and magnetic vector components stored in ``m_1.csv``.





Consistency Rules
-----------------

The following consistency rules should be satisfied:

- Active ``MagData`` and ``NucData`` should use compatible object folders when
  combined in one calculation.
- If ``StructData`` is active, the number of rows in ``StructData.csv`` must
  match the number of local objects.
- If ``RotData`` is active without ``RotDataLoop``, the number of rows in
  ``RotData.csv`` must match the number of local objects.
- If ``RotDataLoop`` is active, every selected ``RotData_#.csv`` file must
  contain the same number of rows.
- If both ``StructData`` and ``RotData`` are active, both metadata layers must
  contain the same number of object entries.
- Files are whitespace-separated and do not use CSV headers.


Future Extension: Replication Import
------------------------------------

A planned extension is a replication-import mode. In this mode, one local
object would be read once and replicated internally in RAM before GPU upload.
The kernels would then receive the same in-memory structure as for a fully
materialized dataset, while avoiding thousands of duplicated object files on
disk.

This mode is not part of the current input format, but it follows naturally
from the same separation between local object data and assembly metadata.
