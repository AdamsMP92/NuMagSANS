"""Generate the RealSpaceData input set for NuMagSANS example 1.

The example uses five identical local spherical magnetic objects and places
them through ``StructData.csv``. Each local object folder contains six data
files so the existing ``User_Selection`` indices can select equivalent test
datasets without storing them manually in the repository.

All CSV files are written without headers and with whitespace separators,
matching the NuMagSANS real-space input format.
"""

from __future__ import annotations

from pathlib import Path
import argparse
import shutil

import numpy as np


BASE_DIR = Path(__file__).resolve().parent
REAL_SPACE_DIR = BASE_DIR / "RealSpaceData"

DEFAULT_OBJECT_CENTERS = np.array(
    [
        (0.0, 0.0, 0.0),
        (10.0, 0.0, 0.0),
        (0.0, 10.0, 0.0),
        (0.0, 0.0, 10.0),
        (20.0, 0.0, 0.0),
    ],
    dtype=float,
)


def sample_spherical_grid(radius: float, spacing: float, points_per_axis: int) -> np.ndarray:
    """Return cubic-grid points clipped by a sphere centered at the origin."""

    if points_per_axis < 1:
        raise ValueError("points_per_axis must be at least 1.")

    if points_per_axis % 2 == 0:
        raise ValueError("points_per_axis must be odd so the grid contains the origin.")

    half_width = points_per_axis // 2
    axis = (np.arange(points_per_axis, dtype=float) - half_width) * spacing
    x_grid, y_grid, z_grid = np.meshgrid(axis, axis, axis, indexing="ij")
    mask = x_grid**2 + y_grid**2 + z_grid**2 <= radius**2
    positions = np.column_stack((x_grid[mask], y_grid[mask], z_grid[mask]))
    order = np.lexsort((positions[:, 1], positions[:, 0], positions[:, 2]))
    return positions[order]


def write_struct_data(real_space_dir: Path, object_centers: np.ndarray) -> None:
    """Write the global object-center positions."""

    np.savetxt(real_space_dir / "StructData.csv", object_centers, fmt="%.10g", delimiter=" ")


def write_magnetic_data(
    real_space_dir: Path,
    local_positions: np.ndarray,
    magnetization: tuple[float, float, float],
    n_objects: int,
    n_datasets: int,
) -> None:
    """Write local magnetic object files as x y z mx my mz."""

    mag_dir = real_space_dir / "MagData"
    moment_rows = np.repeat(np.asarray(magnetization, dtype=float)[None, :], len(local_positions), axis=0)
    magnetic_data = np.column_stack((local_positions, moment_rows))

    for object_index in range(1, n_objects + 1):
        object_dir = mag_dir / f"Angle_{object_index}"
        object_dir.mkdir(parents=True, exist_ok=True)

        for dataset_index in range(1, n_datasets + 1):
            np.savetxt(
                object_dir / f"m_{dataset_index}.csv",
                magnetic_data,
                fmt="%.10g",
                delimiter=" ",
            )


def write_nuclear_data(
    real_space_dir: Path,
    local_positions: np.ndarray,
    nuclear_sld: float,
    n_objects: int,
    n_datasets: int,
) -> None:
    """Write local nuclear object files as x y z rho."""

    nuc_dir = real_space_dir / "NucData"
    sld_rows = np.full((len(local_positions), 1), nuclear_sld, dtype=float)
    nuclear_data = np.column_stack((local_positions, sld_rows))

    for object_index in range(1, n_objects + 1):
        object_dir = nuc_dir / f"Angle_{object_index}"
        object_dir.mkdir(parents=True, exist_ok=True)

        for dataset_index in range(1, n_datasets + 1):
            np.savetxt(
                object_dir / f"n_{dataset_index}.csv",
                nuclear_data,
                fmt="%.10g",
                delimiter=" ",
            )


def generate_real_space_data(
    output_dir: Path = REAL_SPACE_DIR,
    radius: float = 5.0,
    spacing: float = 0.3554,
    points_per_axis: int = 29,
    magnetization: tuple[float, float, float] = (0.0, 0.0, 1.7),
    nuclear_sld: float = 1.0,
    n_datasets: int = 6,
    clean: bool = True,
) -> Path:
    """Generate the complete RealSpaceData folder for example 1."""

    output_dir = Path(output_dir)

    if clean and output_dir.exists():
        shutil.rmtree(output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)
    local_positions = sample_spherical_grid(radius, spacing, points_per_axis)

    write_struct_data(output_dir, DEFAULT_OBJECT_CENTERS)
    write_magnetic_data(output_dir, local_positions, magnetization, len(DEFAULT_OBJECT_CENTERS), n_datasets)
    write_nuclear_data(output_dir, local_positions, nuclear_sld, len(DEFAULT_OBJECT_CENTERS), n_datasets)

    return output_dir


def parse_args() -> argparse.Namespace:
    """Parse command line arguments for manual data generation."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=REAL_SPACE_DIR)
    parser.add_argument("--radius", type=float, default=5.0)
    parser.add_argument("--spacing", type=float, default=0.3554)
    parser.add_argument("--points-per-axis", type=int, default=29)
    parser.add_argument("--mx", type=float, default=0.0)
    parser.add_argument("--my", type=float, default=0.0)
    parser.add_argument("--mz", type=float, default=1.7)
    parser.add_argument("--nuclear-sld", type=float, default=1.0)
    parser.add_argument("--n-datasets", type=int, default=6)
    parser.add_argument("--no-clean", action="store_true")
    return parser.parse_args()


def main() -> None:
    """Generate the example data set from command line parameters."""

    args = parse_args()
    output_dir = generate_real_space_data(
        output_dir=args.output_dir,
        radius=args.radius,
        spacing=args.spacing,
        points_per_axis=args.points_per_axis,
        magnetization=(args.mx, args.my, args.mz),
        nuclear_sld=args.nuclear_sld,
        n_datasets=args.n_datasets,
        clean=not args.no_clean,
    )
    print(f"Generated RealSpaceData in: {output_dir}")


if __name__ == "__main__":
    main()
