"""Analytical magnetic vector fields with spherical or radial length scaling.

This module is intentionally only a local vector-field library. It does not
read object folders and it does not write MagData files. Higher-level
AssemblyBaseMagnetizer code can apply these fields to materialized
Local_Objects and store the resulting magnetic data.
"""

import warnings

import numpy as np


# ------------------------------------------------------------
#  Coordinate helpers
# ------------------------------------------------------------
def _cylindrical(x, y, z):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return rho, phi


def _spherical(x, y, z):
    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arctan2(np.sqrt(x**2 + y**2), z)
    phi = np.arctan2(y, x)
    return r, theta, phi


# ------------------------------------------------------------
#  Parameter checking
# ------------------------------------------------------------
def _check_param_range(name, value, min_val, max_val, fatal=True):
    if not (min_val <= value <= max_val):
        msg = f"Parameter '{name}' = {value} is outside the allowed range [{min_val}, {max_val}]"
        if fatal:
            raise ValueError(msg)
        else:
            warnings.warn(msg)
            return False
    return True


# ============================================================
#  Alpha-profile generator (radial transition functions)
# ============================================================
def alpha_profile(xi, params):
    """Compute tilt angle α(ξ) for various radial transition functions."""
    typ = params.get("profile_type", "linear").lower()
    kappa = params.get("kappa", 1.0)

    if typ == "linear":
        _check_param_range("kappa", kappa, 0.0, 1.0)
        return np.arccos(1 - kappa * xi)

    elif typ == "parabolic":
        _check_param_range("kappa", kappa, 0.0, 1.0)
        return np.arccos(1 - kappa * xi**2)

    elif typ == "trigonometric":
        _check_param_range("kappa", kappa, 0.0, 1.0)
        return 0.5 * np.pi * kappa * xi

    elif typ == "hyperbolic":
        _check_param_range("kappa", kappa, -10.0, 10.0, fatal=False)
        return np.arccos(np.sqrt(1 - np.tanh(kappa * xi) ** 2))

    elif typ == "argument_hyperbolic":
        _check_param_range("kappa", kappa, 0.0, 10.0, fatal=False)
        return 0.5 * np.pi * np.tanh(kappa * xi)

    elif typ == "lorentzian":
        n = params.get("exponent", 2.0)
        _check_param_range("kappa", kappa, 0.1, 5.0)
        _check_param_range("exponent", n, 1.0, 10.0)
        return np.arccos((1 + kappa * xi**2) ** (-n))

    elif typ == "gaussian":
        _check_param_range("kappa", kappa, 0.1, 2.0)
        return np.arccos(np.exp(-(xi**2) / kappa**2))

    elif typ == "arctan":
        _check_param_range("kappa", kappa, 0.0, 10.0, fatal=False)
        return 2.0 * np.arctan(kappa * xi)

    else:
        raise ValueError(f"Unknown profile_type: '{typ}'")


# ============================================================
#  Unified vector-field generator
# ============================================================
def unit_field(x, y, z, D, params):
    """
    Generate unit vector field for various core types.

    Parameters
    ----------
    x, y, z : ndarray
        Cartesian coordinates.
    D : float
        Particle diameter. [in nanometer]
    params : dict
        {
            "field_type": "vortex" | "hedgehog" | "artichoke" |
                          "skyrmion" | "poloidal_vortex" |
                          "transversal_helix",
            "profile_type": ... (see alpha_profile),
            "xi_type": "cylindrical_xi" | "spherical_xi",
            "kappa": float,
            "N": int (for skyrmion),
            "m": int, "gamma": float (for skyrmion base),
            "k": float, "chirality": float, "phase": float
                 (for transversal_helix),
            "turns": float
        }

    Returns
    -------
    mx, my, mz : ndarray
        Normalized components.
    """

    R = D / 2.0
    field_type = params.get("field_type", "vortex").lower()

    if params.get("coord_transform", True):
        # coordinate systems
        rho, phi_cyl = _cylindrical(x, y, z)
        r, theta, phi = _spherical(x, y, z)

    if params.get("coord_normalize", True):
        xi_type = params.get("xi_type", "cylindrical_xi").lower()
        # choose normalized radius ξ
        if xi_type == "spherical_xi":
            xi = np.clip(r / R, 0.0, 1.0)
        elif xi_type == "cylindrical_xi":
            xi = np.clip(rho / R, 0.0, 1.0)
        else:
            raise ValueError(f"Unknown xi_type '{xi_type}' (use 'spherical_xi' or 'cylindrical_xi')")

    # --- basis-dependent definition ---
    if field_type == "linearized_vortex":
        m_0 = params.get("m_0", 1.0)
        m_1 = params.get("m_1", 1.0)
        mx = -m_1 * y
        my = m_1 * x
        mz = np.full_like(x, m_0)

    elif field_type in ("transversal_helix", "transverse_helix"):
        k = params.get("k", 1.0)
        chirality = params.get("chirality", params.get("chi", 1.0))
        phase = params.get("phase", 0.0)
        argument = k * z + phase
        mx = np.cos(argument)
        my = chirality * np.sin(argument)
        mz = np.zeros_like(x)

    elif field_type == "vortex":
        # compute tilt α(ξ)
        alpha = alpha_profile(xi, params)
        alpha *= params.get("turns", 1.0)
        m_phi = np.sin(alpha)
        m_z = np.cos(alpha)
        mx = -m_phi * np.sin(phi)
        my = m_phi * np.cos(phi)
        mz = m_z

    elif field_type == "hedgehog":
        # compute tilt α(ξ)
        alpha = alpha_profile(xi, params)
        alpha *= params.get("turns", 1.0)
        m_r = np.sin(alpha) * np.cos(theta)
        m_z = np.cos(alpha)
        mx = m_r * np.sin(theta) * np.cos(phi)
        my = m_r * np.sin(theta) * np.sin(phi)
        mz = m_r * np.cos(theta) + m_z

    elif field_type == "artichoke":
        # compute tilt α(ξ)
        alpha = alpha_profile(xi, params)
        alpha *= params.get("turns", 1.0)
        m_theta = -np.sin(alpha)
        m_z = np.cos(alpha)
        mx = m_theta * np.cos(phi) * np.cos(theta)
        my = m_theta * np.sin(phi) * np.cos(theta)
        mz = -m_theta * np.sin(theta) + m_z

    elif field_type == "poloidal_vortex":
        # compute tilt α(ξ)
        alpha = alpha_profile(xi, params)
        alpha *= params.get("turns", 1.0)
        m_theta = np.sin(alpha)
        m_z = np.cos(alpha)
        mx = m_theta * np.cos(phi) * np.cos(theta)
        my = m_theta * np.sin(phi) * np.cos(theta)
        mz = m_theta * np.sin(theta) + m_z

    elif field_type == "skyrmion":
        N = params.get("N", 1.0)
        m = params.get("m", 1.0)
        gamma = params.get("gamma", 0.0)
        mx = np.sin(np.pi * N * xi) * np.cos(m * phi + gamma)
        my = np.sin(np.pi * N * xi) * np.sin(m * phi + gamma)
        mz = np.cos(np.pi * N * xi)

    else:
        raise ValueError(f"Unknown field_type '{field_type}'")

    if params.get("additional_normalize", True):
        # normalization
        norm = np.sqrt(mx**2 + my**2 + mz**2) + 1e-12
        mx /= norm
        my /= norm
        mz /= norm

    return mx, my, mz
