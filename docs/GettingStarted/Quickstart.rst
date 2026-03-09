Quickstart
==========

NuMagSANS can be executed through the python interface, which allows
users to generate configuration files and run simulations programmatically.

The repository contains a minimal example together with a small dataset.

Quick System Check
------------------

The following command check your system. 

Linux
^^^^^

.. code-block:: bash

   nvidia-smi
   nvcc --version
   cmake --version
   g++ --version
   python --version

Windows PowerShell
^^^^^^^^^^^^^^^^^^

.. code-block:: powershell

   nvidia-smi
   nvcc --version
   cmake --version
   cl
   python --version

Quick Build NuMagSANS
---------------------

The following commands

1. download NuMagSANS from GitHub
2. compile the C++/CUDA source code using CMake
3. create a virtual environment
4. install the NuMagSANS python interface
5. run the NuMagSANS example script

Linux
^^^^^

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

Windows PowerShell
^^^^^^^^^^^^^^^^^^

.. code-block:: powershell
   
   git clone https://github.com/AdamsMP92/NuMagSANS.git
   cd NuMagSANS
   cmake -S . -B build -DCMAKE_CUDA_ARCHITECTURES=native
   cmake --build build -j
   python -m venv .venv
   .venv\Scripts\Activate.ps1
   pip install --upgrade pip
   pip install -e .
   python -c "from NuMagSANS import NuMagSANS; print('Installation successful')"
   cd example
   python NuMagSANSrun.py

Example Script
--------------

The following provides a quick overview on the NuMagSANS example that is executed in the last step of the Quick build. 
It demonstrates a minimal NuMagSANS workflow.

This example performs the following steps:

1. Create a ``NuMagSANS`` interface object.
2. Generate a temporary configuration file.
3. Specify simulation parameters and input data.
4. Execute the NuMagSANS simulation.
5. Remove the temporary configuration file.

Input data are provided in the repository:

::

   example/RealSpaceData/MagData

Simulation results are written to:

::

   example/NuMagSANS_Output

.. code-block:: python

   from NuMagSANS import NuMagSANS
   from pathlib import Path

   BASE_DIR = Path(__file__).resolve().parent

   # Create NuMagSANS interface object
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


