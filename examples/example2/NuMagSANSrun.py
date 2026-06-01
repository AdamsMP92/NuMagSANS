from pathlib import Path
import shutil
from dataclasses import replace

import numpy as np

from NuMagSANS import NuMagSANS
from UniformSphere import (
    random_two_sphere_case,
    random_zyz_angles,
    scaled_spin_flip_1d,
    system_scattering_volume_m3,
    write_uniform_sphere_case,
)


BASE_DIR = Path(__file__).resolve().parent
CONFIG = BASE_DIR / "NuMagSANSInput_temp.conf"
OUTPUT_DIR = BASE_DIR / "NuMagSANS_Output"
REAL_SPACE_DIR = BASE_DIR / "RealSpaceData"

N_ITERATIONS = 5
SEED = 12345
N_Q = 160
N_THETA = 361
Q_MAX = 1.2
POLARIZATION = (0.0, 0.0, 1.0)
N_ROTDATA = 4


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


def format_vector(vector: tuple[float, float, float]) -> str:
    """Return a compact vector string for terminal tables."""

    return "(" + ", ".join(f"{value: .3f}" for value in vector) + ")"


def write_rotdata_loop_cases(rotation_cases: list[list], real_space_dir: Path) -> None:
    """Write RotData_#.csv files for an inner NuMagSANS rotation-data loop."""

    rotdata_dir = real_space_dir / "RotData"
    rotdata_dir.mkdir(parents=True, exist_ok=True)

    for index, spheres in enumerate(rotation_cases, start=1):
        angles = np.asarray([sphere.angles for sphere in spheres], dtype=float)
        np.savetxt(rotdata_dir / f"RotData_{index}.csv", angles, fmt="%.10e", delimiter=" ")


def sample_rotation_cases(spheres: list, rng: np.random.Generator, n_cases: int) -> list[list]:
    """Return copies of the spheres with independently sampled object rotations."""

    rotation_cases = []
    for _ in range(n_cases):
        rotation_cases.append(
            [replace(sphere, angles=random_zyz_angles(rng)) for sphere in spheres]
        )
    return rotation_cases


def run_single_case(sim: NuMagSANS, iteration: int, rng: np.random.Generator) -> dict:
    """Generate one random two-sphere case, run NuMagSANS, and compare SpinFlip_1D."""

    spheres = random_two_sphere_case(rng)
    scattering_volume = system_scattering_volume_m3(spheres)
    write_uniform_sphere_case(BASE_DIR, spheres, dataset_index=1, clean=True)
    rotation_cases = sample_rotation_cases(spheres, rng, N_ROTDATA)
    write_rotdata_loop_cases(rotation_cases, REAL_SPACE_DIR)

    if OUTPUT_DIR.exists():
        shutil.rmtree(OUTPUT_DIR)

    sim.write_config(
        CONFIG,
        MagData_activate=1,
        StructData_activate=1,
        RotData_activate=1,
        RotDataLoop=1,
        RotDataLoop_From=1,
        RotDataLoop_To=N_ROTDATA,
        FastLoad=1,
        MagDataPath=str(BASE_DIR / "RealSpaceData" / "MagData"),
        StructDataFilename=str(BASE_DIR / "RealSpaceData" / "StructData.csv"),
        RotDataFilename=str(BASE_DIR / "RealSpaceData" / "RotData.csv"),
        RotDataPath=str(BASE_DIR / "RealSpaceData" / "RotData"),
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

    rotdata_results = []
    for rotdata_index, rotated_spheres in enumerate(rotation_cases, start=1):
        q, numagsans_spin_flip = read_spin_flip_1d(
            OUTPUT_DIR / "SANS_1" / f"RotData_{rotdata_index}" / "SANS1D.csv"
        )
        analytic_spin_flip = scaled_spin_flip_1d(
            q_values=q,
            n_theta=N_THETA,
            spheres=rotated_spheres,
            scattering_volume=scattering_volume,
            polarization=POLARIZATION,
        )
        mse = relative_mse(numagsans_spin_flip, analytic_spin_flip)
        rotdata_results.append(
            {
                "rotdata_index": rotdata_index,
                "mse": mse,
                "spheres": rotated_spheres,
            }
        )
        print(f"iteration {iteration}, RotData {rotdata_index}: relative MSE = {mse:.6e}")

    return {
        "iteration": iteration,
        "scattering_volume": scattering_volume,
        "rotdata_results": rotdata_results,
    }


def print_summary_table(results: list[dict]) -> None:
    """Print MSE values together with the generated two-sphere parameters."""

    print("\nSummary:")
    print(
        f"{'it':>2} {'rot':>3} {'MSE':>12} {'particle':>8} {'R':>8} {'a':>8} "
        f"{'position':>30} {'magnetization':>30} {'angles':>30}"
    )
    for result in results:
        for rotdata_result in result["rotdata_results"]:
            for particle_index, sphere in enumerate(rotdata_result["spheres"], start=1):
                mse = f"{rotdata_result['mse']:.6e}" if particle_index == 1 else ""
                iteration = str(result["iteration"]) if particle_index == 1 else ""
                rotdata_index = str(rotdata_result["rotdata_index"]) if particle_index == 1 else ""
                print(
                    f"{iteration:>2} {rotdata_index:>3} {mse:>12} {particle_index:>8} "
                    f"{sphere.radius:8.3f} {sphere.spacing:8.3f} "
                    f"{format_vector(sphere.center):>30} "
                    f"{format_vector(sphere.magnetization):>30} "
                    f"{format_vector(sphere.angles):>30}"
                )

    mse_values = np.asarray(
        [
            rotdata_result["mse"]
            for result in results
            for rotdata_result in result["rotdata_results"]
        ],
        dtype=float,
    )
    print(f"\nmean MSE: {np.mean(mse_values):.6e}")


def cleanup_generated_data() -> None:
    """Remove generated test input and output data."""

    for path in (REAL_SPACE_DIR, OUTPUT_DIR):
        if path.exists():
            shutil.rmtree(path)


def main() -> None:
    """Run five random non-overlapping two-sphere checks."""

    sim = NuMagSANS()
    rng = np.random.default_rng(SEED)
    try:
        results = [run_single_case(sim, index, rng) for index in range(1, N_ITERATIONS + 1)]
        print_summary_table(results)
    finally:
        cleanup_generated_data()


if __name__ == "__main__":
    main()
