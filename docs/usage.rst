Usage
=====

NuMagSANS can be executed through the Python interface, which provides a
convenient way to generate configuration files and run simulations.

The repository contains an example dataset and script demonstrating a
minimal workflow.

---

## Quick Example

Navigate to the example directory:

```bash
cd example
```

Run the example script:

```bash
python NuMagSANSrun.py
```

The simulation will generate output data in the directory:

```
example/NuMagSANS_Output
```

---

## Example Script

The following Python script demonstrates a minimal NuMagSANS workflow.

```python
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
```

---

## Description of the Workflow

The example script performs the following steps:

1. **Create a NuMagSANS interface object**
2. **Generate a temporary configuration file**
3. **Specify simulation parameters and input data**
4. **Run the NuMagSANS executable**
5. **Remove the temporary configuration file**

The required magnetization data are provided in the repository:

```
example/RealSpaceData/MagData
```

Simulation results will be written to:

```
example/NuMagSANS_Output
```

---
