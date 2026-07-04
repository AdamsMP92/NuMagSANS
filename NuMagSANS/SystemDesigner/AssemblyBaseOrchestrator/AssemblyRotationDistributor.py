"""Rotation-angle distributors for assembly-level orientation data.

The functions in this module are intentionally independent of the concrete
real-space data layout. They only generate Z-Y-Z Euler angles that can later be
written to one or more NuMagSANS ``RotData`` files.
"""

import math
from pathlib import Path

import numpy as np
from scipy.stats import qmc


def _validate_n_samples(n_samples):
    """Return ``n_samples`` as a positive integer."""
    n_samples = int(n_samples)
    if n_samples <= 0:
        raise ValueError("n_samples must be a positive integer.")
    return n_samples


def _validate_beta_range(beta_min, beta_max):
    """Validate and return a beta interval in radians."""
    beta_min = float(beta_min)
    beta_max = float(beta_max)

    if not 0.0 <= beta_min <= math.pi:
        raise ValueError("beta_min must be in the interval [0, pi].")
    if not 0.0 <= beta_max <= math.pi:
        raise ValueError("beta_max must be in the interval [0, pi].")
    if beta_max < beta_min:
        raise ValueError("beta_max must be greater than or equal to beta_min.")

    return beta_min, beta_max


def _next_power_of_two(value):
    """Return the smallest power of two greater than or equal to value."""
    return 1 << (int(value) - 1).bit_length()


def _unit_cube_to_zyz(unit_samples, beta_min, beta_max, isotropic=True):
    """Map unit-cube samples to Z-Y-Z Euler angles.

    Parameters
    ----------
    unit_samples : numpy.ndarray
        Array with shape ``(n_samples, 3)`` and values in ``[0, 1)``.
    beta_min, beta_max : float
        Lower and upper beta angle in radians.
    isotropic : bool, optional
        If ``True``, sample beta uniformly in ``cos(beta)``. If ``False``,
        sample beta directly and uniformly in angle space.
    """
    unit_samples = np.asarray(unit_samples, dtype=float)
    if unit_samples.ndim != 2 or unit_samples.shape[1] != 3:
        raise ValueError("unit_samples must have shape (n_samples, 3).")

    beta_min, beta_max = _validate_beta_range(beta_min, beta_max)

    alpha = 2.0 * math.pi * unit_samples[:, 0]
    gamma = 2.0 * math.pi * unit_samples[:, 2]

    if isotropic:
        cos_beta_min = math.cos(beta_min)
        cos_beta_max = math.cos(beta_max)
        cos_beta = cos_beta_min + unit_samples[:, 1] * (cos_beta_max - cos_beta_min)
        beta = np.arccos(np.clip(cos_beta, -1.0, 1.0))
    else:
        beta = beta_min + unit_samples[:, 1] * (beta_max - beta_min)

    return {
        "alpha": alpha,
        "beta": beta,
        "gamma": gamma,
    }


def sobol_zyz_rotations(
    n_samples,
    beta_min=0.0,
    beta_max=math.pi,
    isotropic=True,
    scramble=True,
    seed=None,
):
    """Generate Z-Y-Z Euler-angle samples from a Sobol sequence.

    The returned angles follow the NuMagSANS rotation convention

    ``R = Rz(alpha) @ Ry(beta) @ Rz(gamma)``.

    Parameters
    ----------
    n_samples : int
        Number of angle triples to generate.
    beta_min, beta_max : float, optional
        Beta range in radians. The default covers the full interval
        ``[0, pi]``.
    isotropic : bool, optional
        If ``True`` (default), beta is sampled uniformly in ``cos(beta)``.
        This is the appropriate choice for isotropic orientation axes. If
        ``False``, beta is sampled uniformly as an angle.
    scramble : bool, optional
        If ``True`` (default), use a scrambled Sobol sequence.
    seed : int, optional
        Random seed used by the scrambled Sobol generator.

    Returns
    -------
    dict
        Dictionary with ``"alpha"``, ``"beta"``, and ``"gamma"`` arrays.
    """
    n_samples = _validate_n_samples(n_samples)
    beta_min, beta_max = _validate_beta_range(beta_min, beta_max)

    sampler = qmc.Sobol(d=3, scramble=scramble, seed=seed)

    # Sobol balance properties are best for powers of two. Generate the next
    # balanced block and truncate to support arbitrary requested sample counts.
    n_balanced = _next_power_of_two(n_samples)
    unit_samples = sampler.random_base2(m=int(math.log2(n_balanced)))[:n_samples]

    return _unit_cube_to_zyz(
        unit_samples,
        beta_min=beta_min,
        beta_max=beta_max,
        isotropic=isotropic,
    )


def rotations_to_array(rotations):
    """Convert a rotation dictionary to an ``(n_samples, 3)`` array."""
    try:
        alpha = np.asarray(rotations["alpha"], dtype=float)
        beta = np.asarray(rotations["beta"], dtype=float)
        gamma = np.asarray(rotations["gamma"], dtype=float)
    except KeyError as error:
        raise KeyError("rotations must contain alpha, beta, and gamma arrays.") from error

    if not (alpha.shape == beta.shape == gamma.shape):
        raise ValueError("alpha, beta, and gamma arrays must have matching shapes.")

    return np.column_stack((alpha, beta, gamma))


def write_rotdata_loop(rotations_cases, output_dir, filename_prefix="RotData"):
    """Write a sequence of rotation-angle arrays as ``RotData_#.csv`` files.

    Parameters
    ----------
    rotations_cases : iterable
        Iterable of rotation dictionaries or ``(n_samples, 3)`` arrays.
    output_dir : str or pathlib.Path
        Directory where ``RotData_1.csv``, ``RotData_2.csv``, ... are written.
    filename_prefix : str, optional
        File prefix used before the numeric suffix.

    Returns
    -------
    dict
        Summary with output directory, file paths, and number of written cases.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    output_files = []
    for index, rotations in enumerate(rotations_cases, start=1):
        if isinstance(rotations, dict):
            rotations = rotations_to_array(rotations)
        else:
            rotations = np.asarray(rotations, dtype=float)

        if rotations.ndim != 2 or rotations.shape[1] != 3:
            raise ValueError("Each rotation case must have shape (n_samples, 3).")

        output_path = output_dir / f"{filename_prefix}_{index}.csv"
        np.savetxt(output_path, rotations, fmt="%.12g", delimiter=" ")
        output_files.append(str(output_path))

    return {
        "rotdata_dir": str(output_dir),
        "n_rotdata": len(output_files),
        "output_files": output_files,
    }
