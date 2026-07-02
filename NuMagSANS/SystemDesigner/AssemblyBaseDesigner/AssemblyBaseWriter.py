import csv
from pathlib import Path

import numpy as np

from ..CrystalDesigner.CrystalTemplates import get_crystal_template
from .AssemblyBase import sample_template_params


def _as_object_value(sampled_values, index):
    """Return one sampled parameter value for one object."""
    array = np.asarray(sampled_values)

    if array.ndim == 0:
        return array.item()

    return array[index].item() if hasattr(array[index], "item") else array[index]


def _object_params(sampled_params, index):
    """Collect all template parameters for one object index."""
    return {name: _as_object_value(values, index) for name, values in sampled_params.items()}


def _validate_crystal_positions(crystal):
    """Check that a generated crystal contains aligned x/y/z coordinate arrays."""
    for key in ("x", "y", "z"):
        if key not in crystal:
            raise KeyError(f"Generated crystal is missing coordinate key: {key}")

    lengths = {key: len(crystal[key]) for key in ("x", "y", "z")}
    if len(set(lengths.values())) != 1:
        raise ValueError(f"Generated crystal coordinate arrays must match: {lengths}")


def _write_position_data(path, crystal):
    """Write local object coordinates to pos.csv."""
    positions = np.column_stack((crystal["x"], crystal["y"], crystal["z"]))
    np.savetxt(
        path,
        positions,
        delimiter=",",
        header="x,y,z",
        comments="",
    )


def _write_atomtype_data(path, crystal):
    """Write per-atom type labels if the generated crystal provides them."""
    if "atomtype" not in crystal:
        return

    atomtypes = np.asarray(crystal["atomtype"], dtype=str)
    with open(path, "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["atom_index", "atomtype"])
        for index, atomtype in enumerate(atomtypes):
            writer.writerow([index, atomtype])


def _template_distribution_specs(assembly_base_data, template_name):
    """Return distribution specifications for one template in an assembly base."""
    return assembly_base_data["crystal_templates"][template_name]["parameters"]


def _write_metadata(
    path,
    object_id,
    template_name,
    template_params,
    distribution_specs,
    crystal,
):
    """Write object-level metadata to meta.csv."""
    metadata_rows = [
        ("object_id", object_id),
        ("template_name", template_name),
        ("n_atoms", len(crystal["x"])),
    ]
    metadata_rows.extend((f"template_param.{name}", value) for name, value in template_params.items())

    for param_name, spec in distribution_specs.items():
        for spec_key, spec_value in spec.items():
            metadata_rows.append((f"template_param_distribution.{param_name}.{spec_key}", spec_value))

    if "atomtype" in crystal:
        atomtypes = np.unique(crystal["atomtype"])
        metadata_rows.append(("unique_atomtypes", ";".join(map(str, atomtypes))))

    with open(path, "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["key", "value"])
        writer.writerows(metadata_rows)


def write_assembly_base_objects(
    assembly_base_data,
    template_name,
    n_objects,
    output_dir,
    seed=None,
):
    """Generate and store local crystal objects from an assembly base.

    This writer materializes only the local crystal objects. It does not assign
    object centers, translations, or rotations. A later AssemblyBaseOrganizer
    can operate on the materialized object data and write ``StructData.csv``
    with the object-center arrangement.

    The directory layout is:

        output_dir/
            RealSpaceData/
                Local_Objects/
                    Object_1/
                        pos.csv
                        meta.csv
                        atomtype.csv
                    Object_2/
                        ...

    Parameters
    ----------
    assembly_base_data : dict
        Dictionary returned by ``assembly_base``.
    template_name : str
        Name of the crystal template to materialize.
    n_objects : int
        Number of local object geometries to generate.
    output_dir : str or pathlib.Path
        Base directory where ``RealSpaceData/Local_Objects`` is created.
    seed : int, optional
        Random seed used for sampling template parameters.

    Returns
    -------
    dict
        Summary with output path, template name, object count, and sampled
        template parameters.
    """
    n_objects = int(n_objects)
    if n_objects < 1:
        raise ValueError("n_objects must be at least 1.")

    output_dir = Path(output_dir)
    objects_dir = output_dir / "RealSpaceData" / "Local_Objects"
    objects_dir.mkdir(parents=True, exist_ok=True)

    template = get_crystal_template(template_name)
    distribution_specs = _template_distribution_specs(assembly_base_data, template_name)
    sampled_params = sample_template_params(
        assembly_base_data,
        template_name=template_name,
        size=n_objects,
        seed=seed,
    )

    for index in range(n_objects):
        object_id = index + 1
        template_params = _object_params(sampled_params, index)
        crystal = template(**template_params)
        _validate_crystal_positions(crystal)

        object_dir = objects_dir / f"Object_{object_id}"
        object_dir.mkdir(parents=True, exist_ok=True)

        _write_position_data(object_dir / "pos.csv", crystal)
        _write_atomtype_data(object_dir / "atomtype.csv", crystal)
        _write_metadata(
            object_dir / "meta.csv",
            object_id,
            template_name,
            template_params,
            distribution_specs,
            crystal,
        )

    return {
        "objects_dir": str(objects_dir),
        "template_name": template_name,
        "n_objects": n_objects,
        "sampled_params": sampled_params,
    }
