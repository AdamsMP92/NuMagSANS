Installation
============

NuMagSANS requires a working CUDA toolchain and a C++17-compatible compiler.

This section provides a quick system check followed by the build steps.


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

   cmake -S . -B build

Compile:

.. code-block:: bash

   cmake --build build -j

The executable will be created at:

::

   build/NuMagSANS


Common Issues
-------------

If compilation fails with an error like:

  "No space left on device"

your system temporary directory (/tmp) may be full.
This can happen on shared HPC systems.

You can fix this by setting a private temporary directory:

   mkdir -p $HOME/tmp
   export TMPDIR=$HOME/tmp

Then re-run CMake.


HPC Cluster Example
-------------------

On many clusters:

.. code-block:: bash

   module load gcc/13.2.0
   module load python/3.11.7
   module load cuda/12.6.2
   module load cmake/3.27.9

   cmake -S . -B build -DCMAKE_CUDA_ARCHITECTURES=80
   cmake --build build -j
