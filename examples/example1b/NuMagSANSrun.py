from pathlib import Path

from GenerateRealSpaceDataSystemDesigner import generate_real_space_data_system_designer

from NuMagSANS import NuMagSANS

BASE_DIR = Path(__file__).resolve().parent

# Generate the example input data through the SystemDesigner workflow. This
# example uses NuMagSANS replication import instead of materializing all objects
# as separate MagData folders.
summary = generate_real_space_data_system_designer(BASE_DIR)

# NuMagSANS object
sim = NuMagSANS()

# filename of the NuMagSANS config file
config = BASE_DIR / "NuMagSANSInput_temp.conf"

config_hints = summary["config_hints"]

sim.write_config(
    config,
    MagData_activate=1,
    MagDataPath=config_hints["MagDataPath"],
    MagData_ReplicationImport=config_hints["MagData_ReplicationImport"],
    MagData_NumberOfReplications=config_hints["MagData_NumberOfReplications"],
    StructData_activate=config_hints["StructData_activate"],
    StructDataFilename=config_hints["StructDataFilename"],
    RotData_activate=config_hints["RotData_activate"],
    RotDataPath=config_hints["RotDataPath"],
    RotDataLoop=config_hints["RotDataLoop"],
    RotDataLoop_From=config_hints["RotDataLoop_From"],
    RotDataLoop_To=config_hints["RotDataLoop_To"],
    RotData_User_Selection=config_hints["RotData_User_Selection"],
    foldernameSANSData=str(BASE_DIR / "NuMagSANS_Output"),
    Fourier_Approach="atomistic",
    Loop_Modus=config_hints["Loop_Modus"],
    Loop_From=config_hints["Loop_From"],
    Loop_To=config_hints["Loop_To"],
    User_Selection=config_hints["User_Selection"],
    Scattering_Volume_V=2.618e-24,
    enable_outputs=["SpinFlip_2D", "SpinFlip_1D"],
)

# Run NuMagSANS
sim.run(config)

# Delete temporary config file
sim.config_clear(config)
