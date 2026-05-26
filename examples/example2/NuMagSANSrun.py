from pathlib import Path
import shutil

import numpy as np

from NuMagSANS import NuMagSANS
from UniformSphere import (
    random_two_sphere_case,
    scaled_spin_flip_1d,
    system_scattering_volume_m3,
    write_uniform_sphere_case,
)


BASE_DIR = Path(__file__).resolve().parent
CONFIG = BASE_DIR / "NuMagSANSInput_temp.conf"
OUTPUT_DIR = BASE_DIR / "NuMagSANS_Output"

N_ITERATIONS = 5
SEED = 12345
N_Q = 160
N_THETA = 361
Q_MAX = 1.2
POLARIZATION = (0.0, 0.0, 1.0)


def read_spin_flip_1d(filename: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read q and I_sf from a NuMagSANS SANS1D.csv output file."""

    data = np.genfromtxt(filename, delimiter=",", names=True)
    return np.asarray(data["q"], dtype=float), np.asarray(data["I_sf"], dtype=float)


def relative_mse(reference: np.ndarray, prediction: np.ndarray) -> float:
    """Return the relative MSE without fitting or rescaling the prediction."""

    reference = np.asarray(reference, dtype=float)
    prediction = np.asarray(prediction, dtype=float)
    mse = float(np.mean((reference - prediction) ** 2))
    norm = float(np.mean(reference**2))
    return mse if norm == 0.0 else mse / norm


def run_single_case(sim: NuMagSANS, iteration: int, rng: np.random.Generator) -> float:
    """Generate one random two-sphere case, run NuMagSANS, and compare SpinFlip_1D."""

    spheres = random_two_sphere_case(rng)
    scattering_volume = system_scattering_volume_m3(spheres)
    write_uniform_sphere_case(BASE_DIR, spheres, dataset_index=1, clean=True)

    if OUTPUT_DIR.exists():
        shutil.rmtree(OUTPUT_DIR)

    sim.write_config(
        CONFIG,
        MagData_activate=1,
        StructData_activate=1,
        RotData_activate=1,
        FastLoad=1,
        MagDataPath=str(BASE_DIR / "RealSpaceData" / "MagData"),
        StructDataFilename=str(BASE_DIR / "RealSpaceData" / "StructData.csv"),
        RotDataFilename=str(BASE_DIR / "RealSpaceData" / "RotData.csv"),
        foldernameSANSData=str(OUTPUT_DIR),
        Fourier_Approach="atomistic",
        Loop_Modus=0,
        User_Selection=[1],
        Number_Of_q_Points=N_Q,
        Number_Of_theta_Points=N_THETA,
        q_max=Q_MAX,
        Scattering_Volume_V=scattering_volume,
        Polarization=POLARIZATION,
        enable_outputs=["SpinFlip_1D"],
    )

    try:
        sim.run(CONFIG)
    finally:
        sim.config_clear(CONFIG)

    q, numagsans_spin_flip = read_spin_flip_1d(OUTPUT_DIR / "SANS_1" / "SANS1D.csv")
    analytic_spin_flip = scaled_spin_flip_1d(
        q_values=q,
        n_theta=N_THETA,
        spheres=spheres,
        scattering_volume=scattering_volume,
        polarization=POLARIZATION,
    )
    mse = relative_mse(numagsans_spin_flip, analytic_spin_flip)

    print(f"iteration {iteration}: relative MSE = {mse:.6e}")
    return mse


def main() -> None:
    """Run five random non-overlapping two-sphere checks."""

    sim = NuMagSANS()
    rng = np.random.default_rng(SEED)
    mse_values = [run_single_case(sim, index, rng) for index in range(1, N_ITERATIONS + 1)]

    print("MSE values:")
    for mse in mse_values:
        print(f"{mse:.6e}")


if __name__ == "__main__":
    main()
