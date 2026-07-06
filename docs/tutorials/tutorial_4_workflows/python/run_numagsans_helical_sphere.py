"""Configure and run NuMagSANS on the generated helical sphere dataset."""

import math
from pathlib import Path

from NuMagSANS import NuMagSANS


def main():
    """Write a NuMagSANS config for HelicalSphere_R20 and start the simulation."""
    repo_root = Path(__file__).resolve().parents[4]
    dataset_dir = repo_root / "HelicalSphere_R20"
    real_space_dir = dataset_dir / "RealSpaceData"

    radius_nm = 20.0
    n_replications = 1
    scattering_volume = n_replications * (4.0 / 3.0) * math.pi * (radius_nm * 1e-9) ** 3

    if not real_space_dir.exists():
        raise FileNotFoundError(
            f"Could not find {real_space_dir}. Run helical_single_sphere.py first."
        )

    sim = NuMagSANS()
    config = dataset_dir / "NuMagSANSInput_helical_sphere.conf"

    sim.write_config(
        config,
        MagData_activate=1,
        MagDataPath=str(real_space_dir / "MagData"),
        MagData_ReplicationImport=1,
        MagData_NumberOfReplications=n_replications,
        RotData_activate=1,
        RotDataPath=str(real_space_dir / "RotData"),
        RotDataLoop=1,
        RotDataLoop_From=1,
        RotDataLoop_To=1,
        RotData_User_Selection=[1],
        User_Selection=[1],
        foldernameSANSData=str(dataset_dir / "NuMagSANS_Output"),
        SANSData_Output_Format="csv",
        Fourier_Approach="atomistic",
        Scattering_Volume_V=scattering_volume,
        enable_outputs=["SpinFlip_2D", "SpinFlip_1D"],
    )

    sim.run(config)
    sim.config_clear(config)


if __name__ == "__main__":
    main()
