from NuMagSANS import NuMagSANS

# NuMagSANSFacade object
facade = NuMagSANS()

# filename of the NuMagSANS-Config file
config = "NuMagSANSInput_temp.conf"

# Config schreiben
facade.write_config(
    config, 
    MagData_activate=1, 
    MagDataPath="RealSpaceData/MagData/",
    Fourier_Approach="atomistic",
    User_Selection=[1, 2, 3],
    Scattering_Volume_V=2.618e-24,
    enable_outputs=["SpinFlip_2D", "SpinFlip_1D"]
)

# Run NuMagSANS
facade.run(config)

# Delete temporary config file
facade.config_clear(config)
