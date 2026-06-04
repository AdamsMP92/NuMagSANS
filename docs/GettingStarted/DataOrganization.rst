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

The scenarios below are ordered from simple to complex and build on each
other. Information introduced in an earlier scenario is not repeated in full
later. For example, the format of ``m_1.csv`` is explained in Scenario 1;
later scenarios only describe the additional object folders, data files,
metadata layers, or loop settings.


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


Scenario 2: Single object with several magnetic states
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The next scenario keeps the same single local object, but stores several
magnetic states for this object. This is useful, for example, when comparing
several analytically generated magnetization states or several states exported
from another simulation.

Only the number of magnetic data files changes. The object folder remains
``Object_1``, but it now contains ``m_1.csv``, ``m_2.csv``, ``m_3.csv``, and
so on.

.. code-block:: text

   RealSpaceData/
     MagData/
       Object_1/
         m_1.csv
         m_2.csv
         m_3.csv

The file format is the same as in Scenario 1. The selected file index is
controlled by ``User_Selection`` or by ``Loop_Modus``.

.. code-block:: conf

   MagData_activate = 1;
   MagDataPath = RealSpaceData/MagData;
   Loop_Modus = 0;
   User_Selection = {1, 2, 3};

Alternatively, a consecutive range can be selected through the loop settings:

.. code-block:: conf

   Loop_Modus = 1;
   Loop_From = 1;
   Loop_To = 3;

In this scenario, NuMagSANS repeatedly imports the same local object geometry
with different magnetic vector data.


Scenario 3: Several explicitly stored magnetic objects
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Several local objects can be stored explicitly by adding more object folders
inside ``MagData``. Each active object folder provides the selected magnetic
data file.

.. code-block:: text

   RealSpaceData/
     MagData/
       Object_1/
         m_1.csv
       Object_2/
         m_1.csv
       Object_3/
         m_1.csv

The object order is defined by the enumeration of the object folders. For a
selected data-file index ``k``, every active object folder should contain the
corresponding ``m_k.csv`` file.

Without ``StructData`` or ``RotData``, the coordinates in the local object
files are interpreted directly as the real-space coordinates used for the
calculation. This means that translations or relative object positions must
already be contained in the ``x``, ``y``, ``z`` columns of the object files.


Scenario 4: Magnetic and nuclear local object data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Nuclear scattering data are added through the ``NucData`` layer. The directory
structure follows the same object-folder and data-file indexing convention as
``MagData``.

.. code-block:: text

   RealSpaceData/
     MagData/
       Object_1/
         m_1.csv
     NucData/
       Object_1/
         n_1.csv

The nuclear file is header-free and whitespace-separated. It contains the
local coordinates and the nuclear scattering length density or nuclear
contrast value:

.. code-block:: text

   x y z nuc

For a combined nuclear-magnetic calculation, both data layers are activated:

.. code-block:: conf

   MagData_activate = 1;
   NucData_activate = 1;
   MagDataPath = RealSpaceData/MagData;
   NucDataPath = RealSpaceData/NucData;
   User_Selection = {1};

Pure nuclear calculations use the same ``NucData`` organization, but deactivate
``MagData``.


Scenario 5: Object translations with StructData
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``StructData`` adds an assembly layer on top of the local object data. Instead
of baking the object-center positions into every local object file, the local
objects can be stored in their own coordinate systems and translated by a
global structure file.

.. code-block:: text

   RealSpaceData/
     MagData/
       Object_1/
         m_1.csv
       Object_2/
         m_1.csv
     StructData.csv

The ``StructData.csv`` file is header-free and whitespace-separated. Each row
contains the center translation of one local object:

.. code-block:: text

   x y z

The row order follows the object-folder order. For the example above, row 1
belongs to ``Object_1`` and row 2 belongs to ``Object_2``.

.. code-block:: conf

   MagData_activate = 1;
   StructData_activate = 1;
   MagDataPath = RealSpaceData/MagData;
   StructDataFilename = RealSpaceData/StructData.csv;

The same metadata layer can be used with ``NucData`` or with combined
``MagData + NucData`` calculations.


Scenario 6: Object-wise rotations with RotData
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``RotData`` adds object-wise rotations. This is useful when the same local
object data should be evaluated in different orientations without rewriting
the local coordinates or vector fields on disk.

.. code-block:: text

   RealSpaceData/
     MagData/
       Object_1/
         m_1.csv
       Object_2/
         m_1.csv
     RotData.csv

The ``RotData.csv`` file is header-free and whitespace-separated. Each row
contains the three rotation angles of one local object:

.. code-block:: text

   alpha beta gamma

The angles are interpreted by the rotation kernels using the internal
object-wise rotation convention. The row order follows the object-folder order.

.. code-block:: conf

   MagData_activate = 1;
   RotData_activate = 1;
   MagDataPath = RealSpaceData/MagData;
   RotDataFilename = RealSpaceData/RotData.csv;

The same rotation layer can be used with ``NucData`` or with combined
``MagData + NucData`` calculations.


Scenario 7: Object translations and rotations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``StructData`` and ``RotData`` can be combined. In this case, the local object
data define the object-internal coordinates and fields, ``RotData`` defines
the object-wise rotations, and ``StructData`` defines the object-center
translations.

.. code-block:: text

   RealSpaceData/
     MagData/
       Object_1/
         m_1.csv
       Object_2/
         m_1.csv
     StructData.csv
     RotData.csv

Both metadata files must contain one row per local object, and the row order
must be consistent between them. This scenario is the standard organization
for an explicitly materialized assembly of local objects.

.. code-block:: conf

   MagData_activate = 1;
   StructData_activate = 1;
   RotData_activate = 1;
   MagDataPath = RealSpaceData/MagData;
   StructDataFilename = RealSpaceData/StructData.csv;
   RotDataFilename = RealSpaceData/RotData.csv;

For combined nuclear-magnetic assemblies, the corresponding ``NucData`` layer
is added using the same object-folder and file-index convention introduced
above.


Scenario 8: Sweeping several RotData files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For orientation sweeps, ``RotData`` can also be organized as a folder
containing several rotation datasets. The local object data and optional
``StructData`` are imported once, while NuMagSANS loops over selected
``RotData`` files.

.. code-block:: text

   RealSpaceData/
     MagData/
       Object_1/
         m_1.csv
       Object_2/
         m_1.csv
     StructData.csv
     RotData/
       RotData_1.csv
       RotData_2.csv
       RotData_3.csv

Each ``RotData_#.csv`` file follows the same row convention as the single-file
``RotData.csv`` case. The active rotation files can be selected explicitly:

.. code-block:: conf

   RotData_activate = 1;
   RotDataLoop = 1;
   RotDataPath = RealSpaceData/RotData;
   RotData_User_Selection = {1, 2, 3};

or through a consecutive loop range:

.. code-block:: conf

   RotDataLoop = 1;
   RotDataLoop_From = 1;
   RotDataLoop_To = 3;

This mode is intended for repeated calculations with fixed local object data
and changing object-wise orientations.





Consistency Rules
-----------------

The following consistency rules should be satisfied:

- Active ``MagData`` and ``NucData`` should use compatible object folders when
  combined in one calculation.
- Object folders and data files should end with an underscore followed by an
  enumeration starting from ``1``.
- For a selected data-file index ``k``, each active object folder should
  contain the corresponding ``m_k.csv`` or ``n_k.csv`` file.
- If ``StructData`` is active, the number of rows in ``StructData.csv`` must
  match the number of local objects.
- If ``RotData`` is active without ``RotDataLoop``, the number of rows in
  ``RotData.csv`` must match the number of local objects.
- If ``RotDataLoop`` is active, every selected ``RotData_#.csv`` file must
  contain one row per local object.
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
