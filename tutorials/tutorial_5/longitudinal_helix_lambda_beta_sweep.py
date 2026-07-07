"""Generate a longitudinal-helix lambda and beta-distribution sweep."""

import math
import sys
from pathlib import Path

import numpy as np

BASE_DIR = Path(__file__).resolve().parent
REPO_ROOT = BASE_DIR.parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from NuMagSANS.SystemDesigner.Workflows import write_spherical_replication_vectorfield_sweep  # noqa: E402

WRITE_DIAGNOSTIC_PNG = True


def main():
    """Create one sphere with 16 longitudinal-helix fields and 16 RotData files."""
    output_dir = BASE_DIR / "LongitudinalHelixSweep_R20"

    radius_nm = 20.0
    lattice_constant_nm = 1.0
    pitch_nm = 20.0
    chirality = 1.0
    phase = 0.0
    n_replications = 1024

    lambda_values = np.linspace(0.0, 10.0, 16)
    beta2_values = np.linspace(0.0, math.pi, 16)
    beta1_values = np.zeros_like(beta2_values)

    field_parameter_cases = [
        {
            "field_type": "longitudinal_helix",
            "k": 2.0 * math.pi / pitch_nm,
            "chirality": chirality,
            "phase": phase,
            "lambda": lambda_value,
            "normalize": True,
        }
        for lambda_value in lambda_values
    ]

    summary = write_spherical_replication_vectorfield_sweep(
        output_dir=output_dir,
        R=radius_nm,
        a=lattice_constant_nm,
        atomtype="Fe",
        n_replications=n_replications,
        field_parameter_cases=field_parameter_cases,
        n_rotdata=len(beta2_values),
        beta_min=beta1_values,
        beta_max=beta2_values,
        rotation_seed=1,
    )

    print("Generated longitudinal-helix field variants:")
    for index, lambda_value in enumerate(lambda_values, start=1):
        print(f"  m_{index}.csv: lambda = {lambda_value:g}")

    print("Generated beta-distribution RotData files:")
    for index, beta2_value in enumerate(beta2_values, start=1):
        print(f"  RotData_{index}.csv: beta_1 = 0 deg, beta_2 = {math.degrees(beta2_value):.1f} deg")

    if WRITE_DIAGNOSTIC_PNG:
        from NuMagSANS.SystemDesigner.AssemblyAnalyzer import plot_magnetization_file

        image_dir = output_dir / "images"
        image_dir.mkdir(parents=True, exist_ok=True)

        magdata_dir = Path(summary["real_space_dir"]) / "MagData" / "Object_1"
        for index, lambda_value in enumerate(lambda_values, start=1):
            mag_file = magdata_dir / f"m_{index}.csv"
            image_file = image_dir / f"longitudinal_helix_lambda_{lambda_value:g}.png"
            plot_magnetization_file(
                mag_file,
                vector_scale=1.6,
                center_vectors=True,
                seed=1,
                show_points=False,
                color_by="mx",
                cmap="coolwarm",
                clim=(-1.0, 1.0),
                show=False,
                output_png=image_file,
                arrow_tip_length=0.63,
                arrow_tip_radius=0.15,
                arrow_shaft_radius=0.06,
            )
            print(f"Wrote diagnostic image: {image_file}")


if __name__ == "__main__":
    main()
