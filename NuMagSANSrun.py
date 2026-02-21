from NuMagSANSFacade import NuMagSANSFacade

# NuMagSANSFacade object
facade = NuMagSANSFacade()

# Config schreiben
facade.write_config(
    "NuMagSANSInput.conf", 
    MagData_activate=1, 
    MagDataPath="Realspace/MagData/",
    Fourier_Approach="atomistic",
    User_Selection=[1, 2, 3],
    Scattering_Volume_V=2.618e-24,
    enable_outputs=["SpinFlip_2D"]
)

# Run starten
facade.run("NuMagSANSInput.conf")

