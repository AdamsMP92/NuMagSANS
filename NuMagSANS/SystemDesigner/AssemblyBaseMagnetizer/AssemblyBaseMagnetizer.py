"""Materialize local magnetic vector-field data for assembly-base objects."""

import csv
from pathlib import Path

import numpy as np

from .MagnetizationBase import sample_field_params
from .VectorFieldModels.SphericalVectorFieldLib import unit_field

SPHERICAL_VECTOR_FIELD_TEMPLATE = "spherical_unit_field"


def _object_sort_key(path):
    """Sort Object_1, Object_2, ... by their numeric suffix."""
    try:
        return int(path.name.split("_")[-1])
    except ValueError:
        return path.name


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


def _resolve_objects_dir(assembly_dir):
    """Return the Local_Objects directory for a materialized assembly base."""
    assembly_dir = Path(assembly_dir)

    if assembly_dir.name == "Local_Objects":
        return assembly_dir

    objects_dir = assembly_dir / "RealSpaceData" / "Local_Objects"
    if objects_dir.exists():
        return objects_dir

    raise FileNotFoundError(f"Could not find RealSpaceData/Local_Objects in {assembly_dir}.")


def _read_object_metadata(object_dir):
    """Read one Object_i/meta.csv file."""
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
    """Read local object coordinates from pos.csv."""
    pos_path = Path(object_dir) / "pos.csv"
    if not pos_path.exists():
        raise FileNotFoundError(f"Missing position file: {pos_path}")

    positions = np.loadtxt(pos_path, delimiter=",", skiprows=1)
    if positions.ndim == 1:
        positions = positions.reshape(1, -1)

    if positions.shape[1] != 3:
        raise ValueError(f"Expected three coordinate columns in {pos_path}.")

    return positions[:, 0], positions[:, 1], positions[:, 2]


def _field_params_for_object(sampled_params, index):
    """Collect one object's sampled field parameters."""
    params = {}
    for name, values in sampled_params.items():
        array = np.asarray(values)
        if array.ndim == 0:
            value = array.item()
        else:
            value = array[index]
            if hasattr(value, "item"):
                value = value.item()
        params[name] = value
    return params


def _diameter_from_metadata(metadata):
    """Infer spherical field diameter from object metadata."""
    if "template_param.D" in metadata:
        return float(metadata["template_param.D"])

    if "template_param.R" in metadata:
        return 2.0 * float(metadata["template_param.R"])

    raise KeyError("Could not infer diameter. Expected template_param.R or template_param.D in object metadata.")


def write_spherical_magnetization_data(
    assembly_dir,
    magnetization_base_data,
    field_template_name=SPHERICAL_VECTOR_FIELD_TEMPLATE,
    output_file="m_1.csv",
    seed=None,
):
    """Write local spherical magnetic vector-field data for each object.

    The function reads materialized local objects from
    ``assembly_dir/RealSpaceData/Local_Objects`` and writes magnetic data to
    ``assembly_dir/RealSpaceData/LocalMagData/Object_i/output_file``. The output
    files contain no header and use spaces as separators.

    Parameters
    ----------
    assembly_dir : str or pathlib.Path
        Existing assembly directory containing ``RealSpaceData/Local_Objects``.
    magnetization_base_data : dict
        Magnetization base created by ``magnetization_base`` or a high-level
        magnetization template.
    field_template_name : str, optional
        Field template to sample. Currently only ``"spherical_unit_field"`` is
        materialized by this writer.
    output_file : str, optional
        File name written inside each ``MagData/Object_i`` directory.
    seed : int, optional
        Random seed for sampled field parameters.

    Returns
    -------
    dict
        Summary containing output path, object count, and sampled field
        parameters.
    """
    if field_template_name != SPHERICAL_VECTOR_FIELD_TEMPLATE:
        raise ValueError(
            f"write_spherical_magnetization_data currently supports only {SPHERICAL_VECTOR_FIELD_TEMPLATE!r}."
        )

    assembly_dir = Path(assembly_dir)
    objects_dir = _resolve_objects_dir(assembly_dir)
    object_dirs = sorted(objects_dir.glob("Object_*"), key=_object_sort_key)
    if not object_dirs:
        raise ValueError(f"No Object_* directories found in {objects_dir}.")

    mag_dir = assembly_dir / "RealSpaceData" / "LocalMagData"
    mag_dir.mkdir(parents=True, exist_ok=True)

    sampled_params = sample_field_params(
        magnetization_base_data,
        field_template_name=field_template_name,
        size=len(object_dirs),
        seed=seed,
    )

    for index, object_dir in enumerate(object_dirs):
        object_name = object_dir.name
        metadata = _read_object_metadata(object_dir)
        x, y, z = _read_positions(object_dir)
        diameter = _diameter_from_metadata(metadata)
        field_params = _field_params_for_object(sampled_params, index)

        mx, my, mz = unit_field(x, y, z, diameter, field_params)
        output_data = np.column_stack((x, y, z, mx, my, mz))

        object_mag_dir = mag_dir / object_name
        object_mag_dir.mkdir(parents=True, exist_ok=True)
        np.savetxt(
            object_mag_dir / output_file,
            output_data,
            fmt="%.12g",
            delimiter=" ",
        )

    return {
        "mag_dir": str(mag_dir),
        "field_template_name": field_template_name,
        "output_file": output_file,
        "n_objects": len(object_dirs),
        "sampled_params": sampled_params,
    }
