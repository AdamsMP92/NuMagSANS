"""Workflows for replicated spherical objects with vector-field sweeps."""

import csv
import math
from pathlib import Path

import numpy as np

from ..AssemblyBaseDesigner.AssemblyBaseTemplates import write_monodisperse_spherical_nanoparticle_base
from ..AssemblyBaseMagnetizer.SphericalVectorFieldLib import unit_field
from ..AssemblyBaseOrchestrator import sobol_zyz_rotations, write_rotdata_loop


def _parse_value(value):
    """Convert metadata strings to int or float when possible."""
    try:
        integer = int(value)
    except ValueError:
        integer = None
    else:
        if str(integer) == value:
            return integer

    try:
        return float(value)
    except ValueError:
        return value


def _read_object_metadata(object_dir):
    """Read one ``meta.csv`` file from a materialized local object."""
    meta_path = Path(object_dir) / "meta.csv"
    if not meta_path.exists():
        raise FileNotFoundError(f"Missing metadata file: {meta_path}")

    metadata = {}
    with open(meta_path, newline="") as file:
        reader = csv.DictReader(file)
        for row in reader:
            metadata[row["key"]] = _parse_value(row["value"])

    return metadata


def _read_positions(object_dir):
    """Read local object coordinates from ``pos.csv``."""
    pos_path = Path(object_dir) / "pos.csv"
    if not pos_path.exists():
        raise FileNotFoundError(f"Missing position file: {pos_path}")

    positions = np.loadtxt(pos_path, delimiter=",", skiprows=1)
    if positions.ndim == 1:
        positions = positions.reshape(1, -1)

    if positions.shape[1] != 3:
        raise ValueError(f"Expected three coordinate columns in {pos_path}.")

    return positions[:, 0], positions[:, 1], positions[:, 2]


def _diameter_from_metadata(metadata):
    """Infer spherical field diameter from object metadata."""
    if "template_param.D" in metadata:
        return float(metadata["template_param.D"])

    if "template_param.R" in metadata:
        return 2.0 * float(metadata["template_param.R"])

    raise KeyError("Could not infer diameter. Expected template_param.R or template_param.D in object metadata.")


def _normalize_field_parameter_cases(field_parameter_cases):
    """Return a validated list of vector-field parameter dictionaries."""
    if field_parameter_cases is None:
        field_parameter_cases = [
            {
                "field_type": "vortex",
                "profile_type": "linear",
                "xi_type": "cylindrical_xi",
                "kappa": 1.0,
                "turns": 1.0,
            }
        ]

    if isinstance(field_parameter_cases, dict):
        field_parameter_cases = [field_parameter_cases]

    normalized = [dict(case) for case in field_parameter_cases]
    if not normalized:
        raise ValueError("field_parameter_cases must contain at least one parameter dictionary.")

    for index, field_params in enumerate(normalized, start=1):
        if "field_type" not in field_params:
            raise KeyError(f"field_parameter_cases[{index}] is missing 'field_type'.")

    return normalized


def _write_single_object_field_variants(
    object_dir,
    output_dir,
    field_parameter_cases,
    output_file_prefix="m",
):
    """Write ``m_1.csv``, ``m_2.csv``, ... for one local object."""
    metadata = _read_object_metadata(object_dir)
    x, y, z = _read_positions(object_dir)
    diameter = _diameter_from_metadata(metadata)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    output_files = []
    for index, field_params in enumerate(field_parameter_cases, start=1):
        mx, my, mz = unit_field(x, y, z, diameter, field_params)
        output_data = np.column_stack((x, y, z, mx, my, mz))

        output_path = output_dir / f"{output_file_prefix}_{index}.csv"
        np.savetxt(output_path, output_data, fmt="%.12g", delimiter=" ")
        output_files.append(str(output_path))

    return {
        "magdata_object_dir": str(output_dir),
        "n_field_variants": len(output_files),
        "output_files": output_files,
        "field_parameter_cases": field_parameter_cases,
    }


def _rotation_seed(base_seed, index):
    """Return a deterministic per-case seed or ``None``."""
    if base_seed is None:
        return None
    return int(base_seed) + int(index)


def write_spherical_replication_vectorfield_sweep(
    output_dir,
    R,
    a,
    atomtype,
    n_replications,
    field_parameter_cases=None,
    n_rotdata=1,
    beta_min=0.0,
    beta_max=math.pi,
    rotation_seed=None,
    scramble=True,
):
    """Create a single-object spherical sweep for NuMagSANS replication import.

    This workflow materializes one local spherical crystal object, writes one
    or more local magnetic vector-field variants as ``m_1.csv``, ``m_2.csv``,
    ..., and generates one or more angular Sobol distributions as
    ``RotData_1.csv``, ``RotData_2.csv``, ...

    The intended NuMagSANS usage is a two-dimensional sweep over
    ``User_Selection`` and ``RotDataLoop``:

    ``m_i`` selects a vector-field parameter case, while ``RotData_j`` selects
    an angular distribution realization for the replicated object ensemble.

    Parameters
    ----------
    output_dir : str or pathlib.Path
        Output directory containing the generated ``RealSpaceData`` tree.
    R : float
        Sphere radius used by the crystal template.
    a : float
        Simple-cubic lattice constant.
    atomtype : str
        Per-atom label written by the crystal template.
    n_replications : int
        Number of replicated objects used by NuMagSANS replication import.
    field_parameter_cases : dict or iterable of dict, optional
        Parameter dictionaries passed directly to
        ``SphericalVectorFieldLib.unit_field``.
    n_rotdata : int, optional
        Number of angular distribution files to write.
    beta_min, beta_max : float, optional
        Polar-angle range in radians. Sampling is isotropic in ``cos(beta)``.
    rotation_seed : int, optional
        Base seed for scrambled Sobol rotation distributions.
    scramble : bool, optional
        If ``True`` (default), use scrambled Sobol sequences.

    Returns
    -------
    dict
        Summary with generated file locations and NuMagSANS configuration hints.
    """
    n_replications = int(n_replications)
    if n_replications <= 0:
        raise ValueError("n_replications must be a positive integer.")

    n_rotdata = int(n_rotdata)
    if n_rotdata <= 0:
        raise ValueError("n_rotdata must be a positive integer.")

    output_dir = Path(output_dir)
    real_space_dir = output_dir / "RealSpaceData"
    field_parameter_cases = _normalize_field_parameter_cases(field_parameter_cases)

    object_summary = write_monodisperse_spherical_nanoparticle_base(
        R=R,
        a=a,
        atomtype=atomtype,
        n_objects=1,
        output_dir=output_dir,
        name="single replicated spherical nanoparticle object",
    )

    object_dir = real_space_dir / "Local_Objects" / "Object_1"
    magdata_object_dir = real_space_dir / "MagData" / "Object_1"
    field_summary = _write_single_object_field_variants(
        object_dir=object_dir,
        output_dir=magdata_object_dir,
        field_parameter_cases=field_parameter_cases,
    )

    rotation_cases = [
        sobol_zyz_rotations(
            n_samples=n_replications,
            beta_min=beta_min,
            beta_max=beta_max,
            isotropic=True,
            scramble=scramble,
            seed=_rotation_seed(rotation_seed, index),
        )
        for index in range(n_rotdata)
    ]
    rotdata_summary = write_rotdata_loop(
        rotation_cases,
        output_dir=real_space_dir / "RotData",
    )

    config_hints = {
        "MagData_activate": 1,
        "MagDataPath": str(real_space_dir / "MagData"),
        "MagData_ReplicationImport": 1,
        "MagData_NumberOfReplications": n_replications,
        "RotData_activate": 1,
        "RotDataPath": str(real_space_dir / "RotData"),
        "RotDataLoop": 1,
        "RotDataLoop_From": 1,
        "RotDataLoop_To": n_rotdata,
        "RotData_User_Selection": list(range(1, n_rotdata + 1)),
        "Loop_Modus": 1,
        "Loop_From": 1,
        "Loop_To": len(field_parameter_cases),
        "User_Selection": list(range(1, len(field_parameter_cases) + 1)),
    }

    return {
        "output_dir": str(output_dir),
        "real_space_dir": str(real_space_dir),
        "object_summary": object_summary,
        "field_summary": field_summary,
        "rotdata_summary": rotdata_summary,
        "config_hints": config_hints,
    }
