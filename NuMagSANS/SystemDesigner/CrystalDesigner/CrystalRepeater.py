import numpy as np

CRYSTAL_BASE_KEYS = ("x", "y", "z", "atomtype")


def _copy_value(value):
    """Copy array-like values when possible, otherwise return the value itself."""
    return value.copy() if hasattr(value, "copy") else value


def _copy_crystal_base(crystal_base):
    """Create a shallow copy of the crystal base with copied array values."""
    return {key: _copy_value(value) for key, value in crystal_base.items()}


def _as_1d_array(value, name):
    """Convert a value to a one-dimensional NumPy array."""
    array = np.asarray(value)

    if array.ndim == 0:
        array = array.reshape(1)

    if array.ndim != 1:
        raise ValueError(f"{name} must be a one-dimensional value or array.")

    return array


def _as_vector3(value, name):
    """Convert a value to a three-dimensional vector."""
    vector = np.asarray(value, dtype=float)

    if vector.shape != (3,):
        raise ValueError(f"{name} must be a three-dimensional vector.")

    return vector


def _as_positive_integer(value, name):
    """Validate repetition counts such as N1, N2, and N3."""
    integer = int(value)

    if integer != value:
        raise ValueError(f"{name} must be an integer.")

    if integer < 1:
        raise ValueError(f"{name} must be at least 1.")

    return integer


def _validate_crystal_base(crystal_base):
    """Check that the crystal base has x/y/z/atomtype arrays of equal length."""
    missing_keys = [key for key in CRYSTAL_BASE_KEYS if key not in crystal_base]
    if missing_keys:
        raise KeyError(f"Crystal base is missing keys: {missing_keys}")

    lengths = {key: len(_as_1d_array(crystal_base[key], key)) for key in CRYSTAL_BASE_KEYS}
    if len(set(lengths.values())) != 1:
        raise ValueError(f"Crystal base arrays must have equal length: {lengths}")

    return lengths["x"]


def crystal_repeater(crystal_base, a1, a2, a3, N1, N2, N3, u1, u2, u3):
    """Expand a crystal base by repeated lattice translations.

    Each atom in the base is translated by all linear combinations

        n1 * a1 * u1 + n2 * a2 * u2 + n3 * a3 * u3

    with ``n1 = 0 ... N1 - 1``, ``n2 = 0 ... N2 - 1``, and
    ``n3 = 0 ... N3 - 1``. The resulting crystal keeps atom positions,
    atom-type labels, a copy of the original base, and the lattice metadata
    used to generate the structure.

    Parameters
    ----------
    crystal_base : dict
        Base dictionary with the keys ``"x"``, ``"y"``, ``"z"``, and
        ``"atomtype"``. All arrays must be one-dimensional and have equal
        length.
    a1, a2, a3 : float
        Lattice constants multiplying the basis vectors ``u1``, ``u2``, and
        ``u3``.
    N1, N2, N3 : int
        Number of repetitions along the three lattice directions. Repetitions
        start at zero, so ``N1=1`` means only the original base position along
        direction 1 is included.
    u1, u2, u3 : array-like
        Three-dimensional lattice direction vectors.

    Returns
    -------
    dict
        Expanded crystal dictionary. The arrays ``"x"``, ``"y"``, ``"z"``,
        and ``"atomtype"`` are aligned by index. Additional mapping arrays
        ``"base_atom_index"``, ``"n1"``, ``"n2"``, and ``"n3"`` record where
        each generated atom came from.
    """
    number_of_base_atoms = _validate_crystal_base(crystal_base)

    base_x = _as_1d_array(crystal_base["x"], "x").astype(float)
    base_y = _as_1d_array(crystal_base["y"], "y").astype(float)
    base_z = _as_1d_array(crystal_base["z"], "z").astype(float)
    base_atomtype = _as_1d_array(crystal_base["atomtype"], "atomtype")

    a1 = float(a1)
    a2 = float(a2)
    a3 = float(a3)

    N1 = _as_positive_integer(N1, "N1")
    N2 = _as_positive_integer(N2, "N2")
    N3 = _as_positive_integer(N3, "N3")

    u1 = _as_vector3(u1, "u1")
    u2 = _as_vector3(u2, "u2")
    u3 = _as_vector3(u3, "u3")

    positions = []
    atomtypes = []
    base_atom_indices = []
    n1_indices = []
    n2_indices = []
    n3_indices = []

    base_positions = np.column_stack((base_x, base_y, base_z))

    for n1 in range(N1):
        for n2 in range(N2):
            for n3 in range(N3):
                translation = n1 * a1 * u1 + n2 * a2 * u2 + n3 * a3 * u3
                translated_positions = base_positions + translation

                positions.append(translated_positions)
                atomtypes.append(base_atomtype)
                base_atom_indices.append(np.arange(number_of_base_atoms))
                n1_indices.append(np.full(number_of_base_atoms, n1, dtype=int))
                n2_indices.append(np.full(number_of_base_atoms, n2, dtype=int))
                n3_indices.append(np.full(number_of_base_atoms, n3, dtype=int))

    positions = np.vstack(positions)

    return {
        "x": positions[:, 0],
        "y": positions[:, 1],
        "z": positions[:, 2],
        "atomtype": np.concatenate(atomtypes),
        "base_atom_index": np.concatenate(base_atom_indices),
        "n1": np.concatenate(n1_indices),
        "n2": np.concatenate(n2_indices),
        "n3": np.concatenate(n3_indices),
        "crystal_base": _copy_crystal_base(crystal_base),
        "a1": a1,
        "a2": a2,
        "a3": a3,
        "N1": N1,
        "N2": N2,
        "N3": N3,
        "u1": u1,
        "u2": u2,
        "u3": u3,
    }


def CrystalRepeater(crystal_base, a1, a2, a3, N1, N2, N3, u1, u2, u3):
    """Expand a crystal base using the original prototype function name."""
    return crystal_repeater(crystal_base, a1, a2, a3, N1, N2, N3, u1, u2, u3)
