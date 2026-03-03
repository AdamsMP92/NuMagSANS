Installation
============

NuMagSANS requires a working CUDA toolchain and a C++17-compatible compiler.

This section provides a quick system check followed by build and Python installation steps.


Test Your System
----------------

Before building NuMagSANS, verify that the required tools are available.

1) Check for a CUDA-capable GPU:

.. code-block:: bash

   nvidia-smi

You should see your GPU listed.


2) Check CUDA toolkit (nvcc):

.. code-block:: bash

   nvcc --version

This should print the CUDA compiler version.


3) Check CMake:

.. code-block:: bash

   cmake --version

CMake version 3.18 or newer is required.


4) Check your C++ compiler:

.. code-block:: bash

   g++ --version

Make sure it supports C++17 (GCC >= 9 recommended).

If all commands run without errors, your system is ready.


Build NuMagSANS
---------------

Clone the repository:

.. code-block:: bash

   git clone https://github.com/AdamsMP92/NuMagSANS.git
   cd NuMagSANS

Configure:

.. code-block:: bash

   cmake -S . -B build -DCMAKE_CUDA_ARCHITECTURES=native

Compile:

.. code-block:: bash

   cmake --build build -j

The executable will be created at:

::

   build/NuMagSANS


Python Interface Installation
-----------------------------

NuMagSANS provides a Python interface for configuration and execution.
It is strongly recommended to install it inside a virtual environment.


Create a Virtual Environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   python3 -m venv .venv

Activate it:

.. code-block:: bash

   source .venv/bin/activate

Upgrade pip (recommended):

.. code-block:: bash

   pip install --upgrade pip


Install NuMagSANS (Editable Mode)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   pip install -e .

This makes the package available as:

.. code-block:: python

   from NuMagSANS import NuMagSANS

Verify installation:

.. code-block:: bash

   python -c "from NuMagSANS import NuMagSANS; print('Installation successful')"


Running an Example
------------------

After building and installing:

.. code-block:: bash

   cd example
   python NuMagSANSrun.py


Common Issues
-------------

"No space left on device"
^^^^^^^^^^^^^^^^^^^^^^^^^

If compilation fails with an error like:

::

   No space left on device

your system temporary directory (/tmp) may be full.
This can happen on shared HPC systems.

You can fix this by setting a private temporary directory:

.. code-block:: bash

   mkdir -p $HOME/tmp
   export TMPDIR=$HOME/tmp

Then re-run CMake.


HPC Cluster Example
-------------------

On many clusters:

.. code-block:: bash

   module load gcc
   module load python
   module load cuda
   module load cmake

   cmake -S . -B build -DCMAKE_CUDA_ARCHITECTURES=native
   cmake --build build -j

If using the Python interface on HPC:

.. code-block:: bash

   cd NuMagSANS
   python3 -m venv .venv
   source .venv/bin/activate
   pip install -e .


Important Note on Paths
-----------------------

.. note::

   NuMagSANS expects absolute paths in the configuration file.

   The Python interface automatically resolves paths relative to the example script.
   This ensures deterministic behaviour independent of the current working directory.

