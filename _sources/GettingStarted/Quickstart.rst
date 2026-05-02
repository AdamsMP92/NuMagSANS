Quickstart Installation
=======================

NuMagSANS can be executed through the python interface, which allows
users to generate configuration files and run simulations programmatically.

The repository contains a minimal example together with a small dataset.

Example Installation in HPC Environments
----------------------------------------

Many HPC clusters provide software through environment modules.
Module names and partitions vary between systems.

The following examples illustrate typical setups on two HPC systems.
Adjust the module names according to your cluster.

Example: IRIS HPC (University of Luxembourg)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

https://hpc-docs.uni.lu

Login to a interactive gpu node:

.. code-block:: bash

   salloc -p interactive --gpus 1 --cpus-per-gpu 1 --mem 8G --time=60

Load required modules:

.. code-block:: bash

   module purge
   module load lib/GDRCopy/2.4-GCCcore-13.2.0
   module load lang/Python/3.11.5-GCCcore-13.2.0
   module load system/CUDA/12.6.0
   module load devel/CMake/3.27.6-GCCcore-13.2.0

Check if the required modules are activated:

.. code-block:: bash

   nvidia-smi
   nvcc --version
   cmake --version
   g++ --version
   python --version

4) Build NuMagSANS using the following commands:

1. download NuMagSANS from GitHub
2. compile the C++/CUDA source code using CMake
3. create a virtual environment
4. install the NuMagSANS python interface
5. run the NuMagSANS example script

.. code-block:: bash

   git clone https://github.com/AdamsMP92/NuMagSANS.git
   cd NuMagSANS
   cmake -S . -B build -DCMAKE_CUDA_ARCHITECTURES=native
   cmake --build build -j
   python3 -m venv .venv
   source .venv/bin/activate
   pip install --upgrade pip
   pip install -e .
   python3 -c "from NuMagSANS import NuMagSANS; print('Installation successful')"
   cd example
   python3 NuMagSANSrun.py

Example: MPSD HPC (Hamburg)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

https://computational-science.mpsd.mpg.de/about.html

We assume you are logged in at a mpsd cluster login-node.

1) Login to a interactive gpu node:

.. code-block:: bash

   salloc -p gpu-interactive --gpus 1 --cpus-per-gpu 1 --mem 8G --time=60

2) Load required modules:

.. code-block:: bash

   module purge
   mpsd-modules 25c native

   module load gcc/13.2.0
   module load python/3.11.7
   module load cuda/12.6.2
   module load cmake/3.27.9

3) Check if the required modules are activated:

.. code-block:: bash

   nvidia-smi
   nvcc --version
   cmake --version
   g++ --version
   python --version

4) Build NuMagSANS using the following commands:

1. download NuMagSANS from GitHub
2. compile the C++/CUDA source code using CMake
3. create a virtual environment
4. install the NuMagSANS python interface
5. run the NuMagSANS example script

.. code-block:: bash

   git clone https://github.com/AdamsMP92/NuMagSANS.git
   cd NuMagSANS
   cmake -S . -B build -DCMAKE_CUDA_ARCHITECTURES=native
   cmake --build build -j
   python3 -m venv .venv
   source .venv/bin/activate
   pip install --upgrade pip
   pip install -e .
   python3 -c "from NuMagSANS import NuMagSANS; print('Installation successful')"
   cd example
   python3 NuMagSANSrun.py



