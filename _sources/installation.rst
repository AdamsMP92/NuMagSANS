Installation
============

NuMagSANS is a GPU-accelerated computational engine and requires
a working CUDA toolchain.

Prerequisites
-------------

The following software must be available on your system:

- **CMake** >= 3.18
- **CUDA Toolkit** (nvcc in PATH)
- A C++17 compatible compiler (e.g. GCC >= 9)
- A CUDA-capable NVIDIA GPU

To verify CUDA availability:

.. code-block:: bash

   nvcc --version

To verify CMake:

.. code-block:: bash

   cmake --version


Build from Source
-----------------

Clone the repository:

.. code-block:: bash

   git clone https://github.com/<your-username>/NuMagSANS
   cd NuMagSANS

Configure the project:

.. code-block:: bash

   cmake -S . -B build

Compile:

.. code-block:: bash

   cmake --build build -j

The executable will be located at:

::

   build/NuMagSANS


Specifying GPU Architecture (Optional)
--------------------------------------

If you want to manually specify the CUDA architecture:

.. code-block:: bash

   cmake -S . -B build -DCMAKE_CUDA_ARCHITECTURES=80

Replace ``80`` with the architecture matching your GPU.


HPC Cluster Example
-------------------

On many clusters, CUDA must be loaded as a module:

.. code-block:: bash

   module load cuda
   module load gcc

   cmake -S . -B build -DCMAKE_CUDA_ARCHITECTURES=80
   cmake --build build -j


Debug Build
-----------

To build in debug mode:

.. code-block:: bash

   cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
   cmake --build build


Installation (Optional)
-----------------------

To install the executable system-wide:

.. code-block:: bash

   cmake --install build

This will place the executable in the configured installation prefix.
