import csv
import math
from pathlib import Path


def _resolve_objects_dir(path):
    """Resolve either a project output directory or a local object directory."""
    path = Path(path)

    if path.name in ("Local_Objects", "Objects"):
        return path

    local_objects_dir = path / "RealSpaceData" / "Local_Objects"
    if local_objects_dir.exists():
        return local_objects_dir

    legacy_objects_dir = path / "RealSpaceData" / "Objects"
    if legacy_objects_dir.exists():
        return legacy_objects_dir

    raise FileNotFoundError(
        "Could not find a Local_Objects directory. Pass either the "
        "Local_Objects path or a directory containing RealSpaceData/Local_Objects."
    )


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


def read_object_metadata(object_dir):
    """Read one Object_i/meta.csv file into a metadata dictionary."""
    meta_path = Path(object_dir) / "meta.csv"
    if not meta_path.exists():
        raise FileNotFoundError(f"Missing metadata file: {meta_path}")

    metadata = {}
    with open(meta_path, newline="") as file:
        reader = csv.DictReader(file)
        for row in reader:
            metadata[row["key"]] = _parse_value(row["value"])

    return metadata


def analyze_assembly_base_dataset(path):
    """Extract realized template-parameter distributions from object metadata.

    Parameters
    ----------
    path : str or pathlib.Path
        Either the dataset root containing ``RealSpaceData/Local_Objects`` or
        the ``Local_Objects`` directory itself. Legacy ``RealSpaceData/Objects``
        datasets are still accepted.

    Returns
    -------
    dict
        Analysis dictionary containing raw object records, template parameter
        arrays, and the subset of parameters that are numeric.
    """
    objects_dir = _resolve_objects_dir(path)
    object_dirs = sorted(objects_dir.glob("Object_*"), key=_object_sort_key)

    if not object_dirs:
        raise ValueError(f"No Object_* directories found in {objects_dir}")

    records = []
    for object_dir in object_dirs:
        metadata = read_object_metadata(object_dir)
        metadata["object_dir"] = str(object_dir)
        records.append(metadata)

    parameter_names = sorted(
        {
            key.replace("template_param.", "", 1)
            for record in records
            for key in record
            if key.startswith("template_param.")
        }
    )

    parameters = {}
    numeric_parameters = {}
    parameter_distributions = {}
    for parameter_name in parameter_names:
        key = f"template_param.{parameter_name}"
        values = [record.get(key) for record in records]
        parameters[parameter_name] = values

        try:
            numeric_values = [float(value) for value in values]
        except (TypeError, ValueError):
            continue

        if all(math.isfinite(value) for value in numeric_values):
            numeric_parameters[parameter_name] = numeric_values

        distribution_prefix = f"template_param_distribution.{parameter_name}."
        distribution_spec = {}
        for key, value in records[0].items():
            if key.startswith(distribution_prefix):
                spec_key = key.replace(distribution_prefix, "", 1)
                distribution_spec[spec_key] = value

        if distribution_spec:
            parameter_distributions[parameter_name] = distribution_spec

    return {
        "objects_dir": str(objects_dir),
        "n_objects": len(records),
        "records": records,
        "parameters": parameters,
        "numeric_parameters": numeric_parameters,
        "parameter_distributions": parameter_distributions,
    }


def write_parameter_table(analysis, output_path):
    """Write extracted template parameters to a compact CSV table."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    parameter_names = list(analysis["parameters"].keys())
    fieldnames = ["object_id", "template_name", "n_atoms", *parameter_names]

    with open(output_path, "w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()

        for record in analysis["records"]:
            row = {
                "object_id": record.get("object_id"),
                "template_name": record.get("template_name"),
                "n_atoms": record.get("n_atoms"),
            }
            for parameter_name in parameter_names:
                row[parameter_name] = record.get(f"template_param.{parameter_name}")
            writer.writerow(row)

    return output_path
