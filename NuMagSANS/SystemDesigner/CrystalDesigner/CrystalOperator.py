import numpy as np

# Coordinate fields expected in every crystal dictionary.
#
# The current data model is intentionally lightweight:
#     crystal = {
#         "x": np.ndarray,
#         "y": np.ndarray,
#         "z": np.ndarray,
#         "atomtype": np.ndarray,  # optional per-atom labels
#     }
#
# Any additional array-like field with the same length as the coordinates is
# treated as per-atom data and is transformed or filtered together with the
# coordinate arrays.
COORDINATE_KEYS = ("x", "y", "z")
METADATA_KEYS = ("crystal_base", "a1", "a2", "a3", "N1", "N2", "N3", "u1", "u2", "u3")

# Functions that are available inside implicit operator expressions.
#
# Keeping this list explicit makes expression evaluation more predictable than
# exposing all Python builtins. Example expression:
#     "x**2 + y**2 + z**2 - R**2"
ALLOWED_FUNCTIONS = {
    "abs": np.abs,
    "sqrt": np.sqrt,
    "sin": np.sin,
    "cos": np.cos,
    "tan": np.tan,
    "exp": np.exp,
    "log": np.log,
    "min": np.minimum,
    "max": np.maximum,
}


def _copy_crystal(crystal):
    """Return a shallow copy of a crystal, copying array-like values when possible."""
    return {key: value.copy() if hasattr(value, "copy") else value for key, value in crystal.items()}


def _validate_coordinates(crystal):
    """Check that the crystal has x/y/z arrays with matching lengths."""
    missing_keys = [key for key in COORDINATE_KEYS if key not in crystal]
    if missing_keys:
        raise KeyError(f"Crystal is missing coordinate keys: {missing_keys}")

    lengths = {key: len(crystal[key]) for key in COORDINATE_KEYS}
    if len(set(lengths.values())) != 1:
        raise ValueError(f"Coordinate arrays must have equal length: {lengths}")


def _apply_mask(crystal, mask):
    """Filter all per-atom fields with a boolean mask.

    Fields whose length matches the coordinate length are interpreted as
    per-atom arrays. Scalar values or global metadata fields are copied through
    unchanged. Known lattice metadata fields are never filtered, even if their
    length happens to match the number of atoms.
    """
    return {
        key: value[mask]
        if key not in METADATA_KEYS and hasattr(value, "__len__") and len(value) == len(mask)
        else value
        for key, value in crystal.items()
    }


def crystal_centering(crystal, inplace=False):
    """Center crystal coordinates around their geometric mean.

    Parameters
    ----------
    crystal : dict
        Crystal data containing at least the keys ``"x"``, ``"y"``, and ``"z"``.
        Each coordinate field must be array-like and have the same length.
    inplace : bool, optional
        If ``False`` (default), return a copied crystal and leave the input
        unchanged. If ``True``, modify and return the input dictionary.

    Returns
    -------
    dict
        Crystal data with centered coordinate arrays.
    """
    _validate_coordinates(crystal)

    centered = crystal if inplace else _copy_crystal(crystal)
    for key in COORDINATE_KEYS:
        centered[key] = np.asarray(centered[key]) - np.mean(centered[key])

    return centered


def rotation_matrix_zyz(alpha, beta, gamma):
    """Return the active Z-Y-Z rotation matrix.

    The convention is

        R = Rz(alpha) @ Ry(beta) @ Rz(gamma)

    and the returned matrix acts on column vectors ``[x, y, z]``.
    Angles are given in radians.
    """
    ca = np.cos(alpha)
    sa = np.sin(alpha)
    cb = np.cos(beta)
    sb = np.sin(beta)
    cg = np.cos(gamma)
    sg = np.sin(gamma)

    rz_alpha = np.array(
        [
            [ca, -sa, 0.0],
            [sa, ca, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    ry_beta = np.array(
        [
            [cb, 0.0, sb],
            [0.0, 1.0, 0.0],
            [-sb, 0.0, cb],
        ]
    )
    rz_gamma = np.array(
        [
            [cg, -sg, 0.0],
            [sg, cg, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )

    return rz_alpha @ ry_beta @ rz_gamma


def crystal_rotator_zyz(x, y, z, alpha, beta, gamma):
    """Rotate coordinate arrays with a Z-Y-Z Euler rotation.

    Parameters
    ----------
    x, y, z : array-like
        Coordinate arrays with matching shape.
    alpha, beta, gamma : float
        Rotation angles in radians. The convention is
        ``Rz(alpha) @ Ry(beta) @ Rz(gamma)``.

    Returns
    -------
    tuple of numpy.ndarray
        Rotated ``x``, ``y``, and ``z`` arrays.
    """
    x = np.asarray(x)
    y = np.asarray(y)
    z = np.asarray(z)

    if not (x.shape == y.shape == z.shape):
        raise ValueError("x, y, and z must have matching shapes.")

    positions = np.vstack((x.ravel(), y.ravel(), z.ravel()))
    rotated_positions = rotation_matrix_zyz(alpha, beta, gamma) @ positions

    return (
        rotated_positions[0].reshape(x.shape),
        rotated_positions[1].reshape(y.shape),
        rotated_positions[2].reshape(z.shape),
    )


def crystal_rotate_zyz(crystal, alpha, beta, gamma, inplace=False):
    """Rotate a crystal dictionary with a Z-Y-Z Euler rotation."""
    _validate_coordinates(crystal)

    target = crystal if inplace else _copy_crystal(crystal)
    target["x"], target["y"], target["z"] = crystal_rotator_zyz(
        target["x"],
        target["y"],
        target["z"],
        alpha,
        beta,
        gamma,
    )

    return target


def implicit_crystal_operator(crystal, operator, operator_params=None, inplace=False):
    """Filter a crystal with an implicit inequality.

    The operator string is evaluated once for all positions. Positions are kept
    where the expression value is less than or equal to zero:

        operator(x, y, z, params) <= 0

    Example:
        implicit_crystal_operator(crystal, "x**2 + y**2 + z**2 - R**2", {"R": 5})

    This example keeps all atoms inside or on a sphere with radius ``R``.

    Parameters
    ----------
    crystal : dict
        Crystal data containing at least ``"x"``, ``"y"``, and ``"z"`` arrays.
        Additional per-atom arrays, such as ``"atomtype"``, are filtered with
        the same mask.
    operator : str
        Python expression using ``x``, ``y``, ``z``, values from
        ``operator_params``, ``np``, and functions listed in
        ``ALLOWED_FUNCTIONS``. Use Python power syntax ``**`` instead of ``^``.
    operator_params : dict, optional
        Named scalar or array parameters used by the operator expression.
    inplace : bool, optional
        If ``False`` (default), return a filtered copy. If ``True``, modify and
        return the input dictionary.

    Returns
    -------
    dict
        Filtered crystal data.
    """
    _validate_coordinates(crystal)

    operator_params = operator_params or {}
    target = crystal if inplace else _copy_crystal(crystal)

    namespace = {
        "__builtins__": {},
        "np": np,
        "x": np.asarray(target["x"]),
        "y": np.asarray(target["y"]),
        "z": np.asarray(target["z"]),
        **ALLOWED_FUNCTIONS,
        **operator_params,
    }

    # Evaluate the expression in a restricted namespace. This is suitable for
    # trusted local formulas; avoid passing untrusted user input as operator text.
    values = eval(operator, namespace)
    mask = np.asarray(values) <= 0

    if mask.shape != np.asarray(target["x"]).shape:
        raise ValueError("Implicit operator must return one value per position.")

    filtered = _apply_mask(target, mask)

    if inplace:
        target.clear()
        target.update(filtered)
        return target

    return filtered


# Backwards-compatible aliases for the initial API names.
def CrystalCentering(crystal):
    """Center coordinates in-place using the original prototype function name."""
    return crystal_centering(crystal, inplace=True)


def ImplicitCrystalOperator(crystal, operator, operator_params=None):
    """Apply an implicit crystal filter in-place using the original function name."""
    return implicit_crystal_operator(crystal, operator, operator_params, inplace=True)
