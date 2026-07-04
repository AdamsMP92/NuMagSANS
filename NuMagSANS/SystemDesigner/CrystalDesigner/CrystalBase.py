import numpy as np

CRYSTAL_BASE_KEYS = ("x", "y", "z", "atomtype")


def _as_1d_array(value, name):
    """Convert input values to a one-dimensional NumPy array."""
    array = np.asarray(value)

    if array.ndim == 0:
        array = array.reshape(1)

    if array.ndim != 1:
        raise ValueError(f"{name} must be a one-dimensional value or array.")

    return array


def crystal_base(x, y, z, atomtype):
    """Create a clean crystal base dictionary.

    A crystal base is the smallest structural unit used by the builder. It
    stores one or more atom positions together with a free-form atom type label.
    The atom type is intentionally not validated chemically, so values such as
    ``"Fe"``, ``"Vacancy"``, ``"A"``, or any other label are allowed.

    Parameters
    ----------
    x, y, z : scalar or array-like
        Cartesian coordinate values. All coordinate arrays must have the same
        length after conversion to one-dimensional NumPy arrays.
    atomtype : scalar or array-like
        Atom type label or labels. A single value is broadcast to all positions.
        If multiple values are provided, their length must match the coordinate
        length.

    Returns
    -------
    dict
        Dictionary with exactly the keys ``"x"``, ``"y"``, ``"z"``, and
        ``"atomtype"``. Each value is a one-dimensional NumPy array.
    """
    x_array = _as_1d_array(x, "x")
    y_array = _as_1d_array(y, "y")
    z_array = _as_1d_array(z, "z")

    coordinate_lengths = {
        "x": len(x_array),
        "y": len(y_array),
        "z": len(z_array),
    }
    if len(set(coordinate_lengths.values())) != 1:
        raise ValueError(f"Coordinate arrays must have equal length: {coordinate_lengths}")

    number_of_atoms = len(x_array)
    atomtype_array = _as_1d_array(atomtype, "atomtype")

    if len(atomtype_array) == 1:
        atomtype_array = np.repeat(atomtype_array, number_of_atoms)
    elif len(atomtype_array) != number_of_atoms:
        raise ValueError(
            "atomtype must be a single value or have the same length as the "
            f"coordinate arrays. Got {len(atomtype_array)} atom types for "
            f"{number_of_atoms} positions."
        )

    return {
        "x": x_array,
        "y": y_array,
        "z": z_array,
        "atomtype": atomtype_array,
    }
