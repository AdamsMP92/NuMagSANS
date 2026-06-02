Real-Space Input Data Organization
==================================

NuMagSANS reads real-space input data from a ``RealSpaceData`` directory.
The input organization separates local object data from optional assembly
metadata. This makes it possible to represent fully materialized datasets,
object-based assemblies, and rotation sweeps using the same backend.

The basic directory structure is:

.. code-block:: text

   RealSpaceData/
     MagData/
     NucData/
     StructData.csv
     RotData.csv
     RotData/
       RotData_1.csv
       RotData_2.csv
       ...

Only the data layers that are activated in the configuration file are required.
The examples below mostly use ``MagData`` for compactness. The same directory
logic applies to ``NucData`` and to combined ``MagData + NucData`` calculations.


Local Object Data
-----------------

Local object data are stored in object folders. Each object folder contains one
or more data files. The file index is selected by ``User_Selection`` or by
``Loop_Modus``.

Magnetic data:

.. code-block:: text

   RealSpaceData/
     MagData/
       Object_1/
         m_1.csv
         m_2.csv
         ...
       Object_2/
         m_1.csv
         m_2.csv
         ...

Each magnetic file contains rows with:

.. code-block:: text

   x y z mx my mz

Nuclear data:

.. code-block:: text

   RealSpaceData/
     NucData/
       Object_1/
         n_1.csv
         n_2.csv
         ...
       Object_2/
         n_1.csv
         n_2.csv
         ...

Each nuclear file contains rows with:

.. code-block:: text

   x y z nuc

The object folders define the local object ordering. For example,
``Object_1`` corresponds to the first object, ``Object_2`` to the second
object, and so on.


Mode 1: Fully Materialized Local Data
-------------------------------------

In the simplest mode, NuMagSANS reads only local magnetic and/or nuclear data.
No additional object-center translations or object-wise rotations are applied.

Typical layout:

.. code-block:: text

   RealSpaceData/
     MagData/
       Object_1/
         m_1.csv

or, for combined magnetic and nuclear data:

.. code-block:: text

   RealSpaceData/
     MagData/
       Object_1/
         m_1.csv
     NucData/
       Object_1/
         n_1.csv

Typical configuration:

.. code-block:: conf

   MagData_activate = 1;
   NucData_activate = 0;
   StructData_activate = 0;
   RotData_activate = 0;
   User_Selection = {1};

This mode is suitable for a single complete object or for a fully materialized
system that is stored directly as local real-space coordinates.


Mode 2: Materialized Multi-Object Data
--------------------------------------

Multiple local objects can be provided explicitly as separate object folders.
The kernels treat them as separate objects even without ``StructData`` or
``RotData``.

Typical layout:

.. code-block:: text

   RealSpaceData/
     MagData/
       Object_1/
         m_1.csv
       Object_2/
         m_1.csv
       Object_3/
         m_1.csv

Typical configuration:

.. code-block:: conf

   MagData_activate = 1;
   StructData_activate = 0;
   RotData_activate = 0;
   User_Selection = {1};

In this mode, positions and orientations are already contained in the local
coordinates stored in the individual object files.


Mode 3: StructData Assembly
---------------------------

``StructData`` stores object-center translations. It is used when local object
data should be placed at object-center positions during the scattering
calculation.

Typical layout:

.. code-block:: text

   RealSpaceData/
     MagData/
       Object_1/
         m_1.csv
       Object_2/
         m_1.csv
     StructData.csv

``StructData.csv`` contains one row per object:

.. code-block:: text

   x y z

Typical configuration:

.. code-block:: conf

   MagData_activate = 1;
   StructData_activate = 1;
   RotData_activate = 0;
   StructDataFilename = RealSpaceData/StructData.csv;
   User_Selection = {1};

The number of rows in ``StructData.csv`` must match the number of local
objects.


Mode 4: RotData Assembly
------------------------

``RotData`` stores object-wise local rotations. It is used when local object
data should be rotated individually during the scattering calculation.

Typical layout:

.. code-block:: text

   RealSpaceData/
     MagData/
       Object_1/
         m_1.csv
       Object_2/
         m_1.csv
     RotData.csv

``RotData.csv`` contains one row per object:

.. code-block:: text

   alpha beta gamma

The angles are interpreted as ZYZ Euler angles in radians:

.. math::

   R = R_z(\alpha) R_y(\beta) R_z(\gamma).

Typical configuration:

.. code-block:: conf

   MagData_activate = 1;
   StructData_activate = 0;
   RotData_activate = 1;
   RotDataFilename = RealSpaceData/RotData.csv;
   User_Selection = {1};

The number of rows in ``RotData.csv`` must match the number of local objects.


Mode 5: StructData and RotData Assembly
---------------------------------------

``StructData`` and ``RotData`` can be combined. In this case, local objects are
rotated individually and placed at object-center positions.

Typical layout:

.. code-block:: text

   RealSpaceData/
     MagData/
       Object_1/
         m_1.csv
       Object_2/
         m_1.csv
     StructData.csv
     RotData.csv

Typical configuration:

.. code-block:: conf

   MagData_activate = 1;
   StructData_activate = 1;
   RotData_activate = 1;
   StructDataFilename = RealSpaceData/StructData.csv;
   RotDataFilename = RealSpaceData/RotData.csv;
   User_Selection = {1};

The number of local objects, the number of ``StructData`` rows, and the number
of ``RotData`` rows must be identical.


Mode 6: Data-File Selection and Outer Data Loops
------------------------------------------------

The file index inside each object folder is selected globally. For example,
``m_1.csv`` and ``n_1.csv`` correspond to data-file index ``1``.

Explicit selection:

.. code-block:: conf

   Loop_Modus = 0;
   User_Selection = {1, 3, 5};

Loop mode:

.. code-block:: conf

   Loop_Modus = 1;
   Loop_From = 1;
   Loop_To = 20;

In loop mode, NuMagSANS evaluates the selected data-file indices one after
another. Output folders are generated as:

.. code-block:: text

   NuMagSANS_Output/
     SANS_1/
     SANS_2/
     ...


Mode 7: RotData Folder Loop
---------------------------

The RotData folder loop evaluates several rotation-data files for the same
selected magnetic, nuclear, and structure data.

Typical layout:

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

Typical configuration using a contiguous range:

.. code-block:: conf

   RotData_activate = 1;
   RotDataLoop = 1;
   RotDataPath = RealSpaceData/RotData;
   RotDataLoop_From = 1;
   RotDataLoop_To = 3;

Typical configuration using an explicit selection:

.. code-block:: conf

   RotData_activate = 1;
   RotDataLoop = 1;
   RotDataPath = RealSpaceData/RotData;
   RotData_User_Selection = {1, 3};

Each selected ``RotData_#.csv`` file must contain the same number of rows, and
this number must match the number of local objects. Output is written into
nested folders:

.. code-block:: text

   NuMagSANS_Output/
     SANS_1/
       RotData_1/
       RotData_3/


Supported Data-Layer Combinations
---------------------------------

The atomistic backend supports magnetic data, nuclear data, or both. Each of
these can be evaluated with no assembly metadata, with ``StructData``, with
``RotData``, or with both.

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
