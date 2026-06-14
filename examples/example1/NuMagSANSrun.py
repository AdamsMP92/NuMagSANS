from pathlib import Path

from GenerateRealSpaceData import generate_real_space_data

from NuMagSANS import NuMagSANS

BASE_DIR = Path(__file__).resolve().parent

# Generate the example input data locally so the repository does not need to
# carry the materialized RealSpaceData files.
generate_real_space_data(BASE_DIR / "RealSpaceData")

# NuMagSANS object
sim = NuMagSANS()

# filename of the NuMagSANS-Config file
config = BASE_DIR / "NuMagSANSInput_temp.conf"

# Config schreiben
sim.write_config(
    config,
    MagData_activate=1,
    MagDataPath=str(BASE_DIR / "RealSpaceData" / "MagData"),
    foldernameSANSData=str(BASE_DIR / "NuMagSANS_Output"),
    Fourier_Approach="atomistic",
    User_Selection=[1, 2, 3],
    Scattering_Volume_V=2.618e-24,
    enable_outputs=["SpinFlip_2D", "SpinFlip_1D"],
)

# Run NuMagSANS
sim.run(config)

# Delete temporary config file
sim.config_clear(config)
