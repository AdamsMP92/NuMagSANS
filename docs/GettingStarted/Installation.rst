Installation
============

NuMagSANS requires a working CUDA toolchain and a C++17-compatible compiler.

This section provides a quick system check followed by build and Python installation steps.

Quick Build NuMagSANS
---------------------

The following command lines 
(1) download NuMagSANS from GitHub
(2) compile the C++/cuda source code using cmake
(3) create a virtual environment
(4) install the NuMagSANS python interface
(5) run the NuMagSANS example script

.. code-block:: bash

   git clone https://github.com/AdamsMP92/NuMagSANS.git
   cd NuMagSANS
   cmake -S . -B build -DCMAKE_CUDA_ARCHITECTURES=native
   cmake --build build -j
   python3 -m venv .venv
   source .venv/bin/activate
   pip install --upgrade pip
   pip install -e .
   python -c "from NuMagSANS import NuMagSANS; print('Installation successful')"
   cd example
   python NuMagSANSrun.py
   
Test Your System
----------------

Before building NuMagSANS, verify that the required tools are available on your system.

The following checks ensure that CUDA, the compiler toolchain, and the build system are correctly installed.

1) Check for a CUDA-capable GPU
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   nvidia-smi

This command should display your GPU and the installed driver version.

Typical output may look like:

::

   +-----------------------------------------------------------------------------+
   | NVIDIA-SMI 535.xx       Driver Version: 535.xx       CUDA Version: 12.x     |
   | GPU  Name                    Persistence-M | Bus-Id        Disp.A | Volatile |
   | 0    NVIDIA RTX 3080                         Off  | 00000000:01:00.0  Off     |
   +-----------------------------------------------------------------------------+

If this command fails:

- Ensure that an **NVIDIA GPU** is installed.
- Install the **NVIDIA driver** for your system.
- On HPC systems, you may need to load a module:

.. code-block:: bash

   module load cuda


2) Check the CUDA Toolkit
^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   nvcc --version

Typical output:

::

   nvcc: NVIDIA (R) Cuda compiler driver
   release 12.3, V12.3.x

NuMagSANS generally works with **CUDA 11 or newer**.

If ``nvcc`` is not found:

- Ensure the CUDA toolkit is installed.
- Add CUDA to your PATH, for example:

.. code-block:: bash

   export PATH=/usr/local/cuda/bin:$PATH


3) Check CMake
^^^^^^^^^^^^^^

.. code-block:: bash

   cmake --version

Expected output:

::

   cmake version 3.22.1

NuMagSANS requires **CMake 3.18 or newer**.

If CMake is missing:

Install it via your system package manager.

Examples:

.. code-block:: bash

   sudo apt install cmake

or on HPC systems:

.. code-block:: bash

   module load cmake


4) Check your C++ compiler
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   g++ --version

Example output:

::

   g++ (GCC) 11.4.0

Your compiler must support **C++17**.

Recommended versions:

- GCC ≥ 9
- Clang ≥ 10

If your compiler is too old:

- Install a newer GCC version
- or load a module on HPC systems:

.. code-block:: bash

   module load gcc


If all commands run successfully, your system is ready to build NuMagSANS.

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
