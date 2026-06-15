"""Utilities for two uniform spherical nanoparticles.

This module is intentionally only a small test helper. It can generate
NuMagSANS-compatible real-space input data for a set of uniform spheres and it
can compute the corresponding analytical Fourier magnetization.

No data are generated on import. Call ``write_uniform_sphere_case`` explicitly
when a concrete test setup should be materialized.
"""

from __future__ import annotations

import shutil
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.special import jv


@dataclass(frozen=True)
class UniformSphere:
    """Definition of one uniform spherical particle in local coordinates."""

    radius: float
    spacing: float
    magnetization: tuple[float, float, float]
    center: tuple[float, float, float]
    angles: tuple[float, float, float] = (0.0, 0.0, 0.0)


def rotation_z(angle: float) -> np.ndarray:
    """Return an active z-axis rotation matrix."""

    c = np.cos(angle)
    s = np.sin(angle)
    return np.array(
        [
            [c, -s, 0.0],
            [s, c, 0.0],
            [0.0, 0.0, 1.0],
        ],
        dtype=float,
    )


def rotation_y(angle: float) -> np.ndarray:
    """Return an active y-axis rotation matrix."""

    c = np.cos(angle)
    s = np.sin(angle)
    return np.array(
        [
            [c, 0.0, s],
            [0.0, 1.0, 0.0],
            [-s, 0.0, c],
        ],
        dtype=float,
    )


def rotation_zyz(alpha: float, beta: float, gamma: float) -> np.ndarray:
    """Return Rz(alpha) Ry(beta) Rz(gamma)."""

    return rotation_z(alpha) @ rotation_y(beta) @ rotation_z(gamma)


def sample_sphere(radius: float, spacing: float) -> np.ndarray:
    """Sample local Cartesian positions on a cubic grid clipped by a sphere."""

    axis = np.arange(-radius, radius + 0.5 * spacing, spacing, dtype=float)
    x_grid, y_grid, z_grid = np.meshgrid(axis, axis, axis, indexing="ij")
    mask = x_grid**2 + y_grid**2 + z_grid**2 <= radius**2
    return np.column_stack((x_grid[mask], y_grid[mask], z_grid[mask]))


def sphere_form_factor(q: np.ndarray, radius: float) -> np.ndarray:
    """Return integral exp(-i q*r) dV over a sphere of radius R."""

    q = np.asarray(q, dtype=float)
    u = q * radius
    volume = 4.0 * np.pi * radius**3 / 3.0
    out = np.empty_like(u, dtype=float)

    small = np.abs(u) < 1.0e-8
    u_small = u[small]
    out[small] = volume * (1.0 - u_small**2 / 10.0 + u_small**4 / 280.0)

    u_large = u[~small]
    out[~small] = volume * 3.0 * (np.sin(u_large) - u_large * np.cos(u_large)) / u_large**3
    return out


def sphere_volume_nm3(radius: float) -> float:
    """Return the volume of one sphere in nm^3."""

    return 4.0 * np.pi * radius**3 / 3.0


def system_scattering_volume_m3(spheres: list[UniformSphere]) -> float:
    """Return the total particle volume in m^3."""

    volume_nm3 = sum(sphere_volume_nm3(sphere.radius) for sphere in spheres)
    return volume_nm3 * 1.0e-27


def sphere_number_density_nm3(sphere: UniformSphere) -> float:
    """Return the point density represented by the sampled cubic grid."""

    return 1.0 / sphere.spacing**3


def analytical_fourier_magnetization(
    qx: np.ndarray,
    qy: np.ndarray,
    qz: np.ndarray,
    spheres: list[UniformSphere],
    global_angles: tuple[float, float, float] = (0.0, 0.0, 0.0),
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Analytical Fourier magnetization of uniform translated/rotated spheres."""

    qx = np.asarray(qx, dtype=float)
    qy = np.asarray(qy, dtype=float)
    qz = np.asarray(qz, dtype=float)
    q = np.sqrt(qx**2 + qy**2 + qz**2)

    global_rotation = rotation_zyz(*global_angles)
    mx_q = np.zeros_like(q, dtype=complex)
    my_q = np.zeros_like(q, dtype=complex)
    mz_q = np.zeros_like(q, dtype=complex)

    for sphere in spheres:
        object_rotation = rotation_zyz(*sphere.angles)
        center = np.asarray(sphere.center, dtype=float)
        moment = np.asarray(sphere.magnetization, dtype=float)

        center_global = global_rotation @ center
        moment_global = global_rotation @ object_rotation @ moment

        phase = np.exp(-1j * (qx * center_global[0] + qy * center_global[1] + qz * center_global[2]))
        shape = sphere_form_factor(q, sphere.radius) * sphere_number_density_nm3(sphere)

        mx_q += moment_global[0] * shape * phase
        my_q += moment_global[1] * shape * phase
        mz_q += moment_global[2] * shape * phase

    return mx_q, my_q, mz_q


def spin_flip_2d_from_fourier(
    q: np.ndarray,
    theta: np.ndarray,
    spheres: list[UniformSphere],
    polarization: tuple[float, float, float] = (0.0, 0.0, 1.0),
    global_angles: tuple[float, float, float] = (0.0, 0.0, 0.0),
) -> np.ndarray:
    """Return the unscaled magnetic spin-flip intensity on NuMagSANS q/theta points."""

    q = np.asarray(q, dtype=float)
    theta = np.asarray(theta, dtype=float)
    qx = np.zeros_like(q)
    qy = q * np.sin(theta)
    qz = q * np.cos(theta)

    mx_q, my_q, mz_q = analytical_fourier_magnetization(
        qx=qx,
        qy=qy,
        qz=qz,
        spheres=spheres,
        global_angles=global_angles,
    )

    sin_theta = np.sin(theta)
    cos_theta = np.cos(theta)
    qx_proj = -mx_q
    qy_proj = (mz_q * sin_theta - my_q * cos_theta) * cos_theta
    qz_proj = (my_q * cos_theta - mz_q * sin_theta) * sin_theta

    p = np.asarray(polarization, dtype=float)
    p = p / np.linalg.norm(p)

    unpolarized = np.abs(qx_proj) ** 2 + np.abs(qy_proj) ** 2 + np.abs(qz_proj) ** 2
    polarized = np.abs(p[0] * qx_proj + p[1] * qy_proj + p[2] * qz_proj) ** 2
    return np.real(unpolarized - polarized)


def spin_flip_1d(
    q_values: np.ndarray,
    n_theta: int,
    spheres: list[UniformSphere],
    polarization: tuple[float, float, float] = (0.0, 0.0, 1.0),
    global_angles: tuple[float, float, float] = (0.0, 0.0, 0.0),
) -> np.ndarray:
    """Azimuthally average the analytical spin-flip intensity like NuMagSANS."""

    q_values = np.asarray(q_values, dtype=float)
    dtheta = 2.0 * np.pi / (float(n_theta) - 1.0)
    theta = np.arange(n_theta, dtype=float) * dtheta

    out = np.zeros_like(q_values, dtype=float)
    for index, q_value in enumerate(q_values):
        q = np.full(n_theta, q_value, dtype=float)
        spin_flip_2d = spin_flip_2d_from_fourier(
            q=q,
            theta=theta,
            spheres=spheres,
            polarization=polarization,
            global_angles=global_angles,
        )
        out[index] = np.sum(spin_flip_2d[:-1] + spin_flip_2d[1:]) / (4.0 * np.pi) * dtheta

    return out


def _poly_mul(a: dict[int, complex], b: dict[int, complex]) -> dict[int, complex]:
    """Multiply two Fourier polynomials in exp(i n theta)."""

    out = {}
    for n_a, c_a in a.items():
        for n_b, c_b in b.items():
            n = n_a + n_b
            out[n] = out.get(n, 0.0j) + c_a * c_b
    return out


def _poly_pow(a: dict[int, complex], exponent: int) -> dict[int, complex]:
    """Raise a Fourier polynomial to a non-negative integer power."""

    out = {0: 1.0 + 0.0j}
    for _ in range(exponent):
        out = _poly_mul(out, a)
    return out


def _angular_phase_average(
    coeffs: dict[int, complex],
    q_value: float,
    dy: float,
    dz: float,
) -> complex:
    """Return <p(theta) exp(-i q (dy sin theta + dz cos theta))>_theta."""

    rho = q_value * np.hypot(dy, dz)
    if rho == 0.0:
        return coeffs.get(0, 0.0j)

    phi = np.arctan2(dy, dz)
    out = 0.0j
    for n, coeff in coeffs.items():
        out += coeff * np.exp(1j * n * phi) * (-1j) ** n * jv(n, rho)
    return out


def spin_flip_1d_exact(
    q_values: np.ndarray,
    spheres: list[UniformSphere],
    polarization: tuple[float, float, float] = (0.0, 0.0, 1.0),
    global_angles: tuple[float, float, float] = (0.0, 0.0, 0.0),
) -> np.ndarray:
    """Closed-form azimuthal SpinFlip_1D for uniform spheres and P || z."""

    p = np.asarray(polarization, dtype=float)
    p = p / np.linalg.norm(p)
    if not np.allclose(p, (0.0, 0.0, 1.0)):
        raise NotImplementedError("The closed-form 1D expression is implemented for P || z.")

    sin_poly = {-1: 0.5j, 1: -0.5j}
    cos_poly = {-1: 0.5, 1: 0.5}
    one = {0: 1.0 + 0.0j}
    cos4 = _poly_pow(cos_poly, 4)
    sin2_cos2 = _poly_mul(_poly_pow(sin_poly, 2), _poly_pow(cos_poly, 2))
    sin_cos3 = _poly_mul(sin_poly, _poly_pow(cos_poly, 3))

    global_rotation = rotation_zyz(*global_angles)
    centers = []
    moments = []

    for sphere in spheres:
        object_rotation = rotation_zyz(*sphere.angles)
        centers.append(global_rotation @ np.asarray(sphere.center, dtype=float))
        moments.append(global_rotation @ object_rotation @ np.asarray(sphere.magnetization, dtype=float))

    q_values = np.asarray(q_values, dtype=float)
    out = np.zeros_like(q_values, dtype=float)

    for q_index, q_value in enumerate(q_values):
        shapes = [
            sphere_form_factor(np.asarray([q_value]), sphere.radius)[0] * sphere_number_density_nm3(sphere)
            for sphere in spheres
        ]
        value = 0.0j

        for i, center_i in enumerate(centers):
            for j, center_j in enumerate(centers):
                moment_i = moments[i]
                moment_j = moments[j]
                dy = center_i[1] - center_j[1]
                dz = center_i[2] - center_j[2]

                value += (
                    shapes[i]
                    * shapes[j]
                    * (
                        moment_i[0] * moment_j[0] * _angular_phase_average(one, q_value, dy, dz)
                        + moment_i[1] * moment_j[1] * _angular_phase_average(cos4, q_value, dy, dz)
                        + moment_i[2] * moment_j[2] * _angular_phase_average(sin2_cos2, q_value, dy, dz)
                        - (moment_i[1] * moment_j[2] + moment_i[2] * moment_j[1])
                        * _angular_phase_average(sin_cos3, q_value, dy, dz)
                    )
                )

        out[q_index] = float(np.real(value))

    return out


def magnetic_sans_prefactor(scattering_volume: float) -> float:
    """Return the atomistic NuMagSANS magnetic SANS prefactor."""

    mu_b = 9.2740100783e-24
    b_h = 2.91e8
    return (mu_b * b_h) ** 2 / scattering_volume * 1.0e-2


def scaled_spin_flip_1d(
    q_values: np.ndarray,
    n_theta: int,
    spheres: list[UniformSphere],
    scattering_volume: float,
    polarization: tuple[float, float, float] = (0.0, 0.0, 1.0),
    global_angles: tuple[float, float, float] = (0.0, 0.0, 0.0),
) -> np.ndarray:
    """Return analytically scaled SpinFlip_1D in the NuMagSANS convention."""

    return magnetic_sans_prefactor(scattering_volume) * spin_flip_1d_exact(
        q_values=q_values,
        spheres=spheres,
        polarization=polarization,
        global_angles=global_angles,
    )


def random_unit_vector(rng: np.random.Generator) -> np.ndarray:
    """Sample a random unit vector with an isotropic direction."""

    vector = rng.normal(size=3)
    return vector / np.linalg.norm(vector)


def random_zyz_angles(rng: np.random.Generator) -> tuple[float, float, float]:
    """Sample ZYZ Euler angles with an isotropic middle angle."""

    alpha = rng.uniform(0.0, 2.0 * np.pi)
    beta = np.arccos(rng.uniform(-1.0, 1.0))
    gamma = rng.uniform(0.0, 2.0 * np.pi)
    return (float(alpha), float(beta), float(gamma))


def random_two_sphere_case(rng: np.random.Generator) -> list[UniformSphere]:
    """Create two random, non-overlapping uniform spheres."""

    r1 = rng.uniform(5.0, 9.0)
    r2 = rng.uniform(5.0, 9.0)
    direction = random_unit_vector(rng)
    distance = (r1 + r2) * rng.uniform(1.25, 1.8)
    center1 = -0.5 * distance * direction
    center2 = 0.5 * distance * direction

    m1 = rng.uniform(0.7, 1.2) * random_unit_vector(rng)
    m2 = rng.uniform(0.7, 1.2) * random_unit_vector(rng)

    return [
        UniformSphere(
            radius=float(r1),
            spacing=float(rng.uniform(0.9, 1.3)),
            magnetization=tuple(m1),
            center=tuple(center1),
            angles=random_zyz_angles(rng),
        ),
        UniformSphere(
            radius=float(r2),
            spacing=float(rng.uniform(0.9, 1.3)),
            magnetization=tuple(m2),
            center=tuple(center2),
            angles=random_zyz_angles(rng),
        ),
    ]


def write_uniform_sphere_case(
    output_dir: Path,
    spheres: list[UniformSphere],
    dataset_index: int = 1,
    clean: bool = True,
) -> Path:
    """Write RealSpaceData with MagData, StructData.csv, and RotData.csv."""

    output_dir = Path(output_dir)
    real_space_dir = output_dir / "RealSpaceData"

    if clean and real_space_dir.exists():
        shutil.rmtree(real_space_dir)

    mag_dir = real_space_dir / "MagData"
    mag_dir.mkdir(parents=True, exist_ok=True)

    struct_rows = []
    rot_rows = []

    for index, sphere in enumerate(spheres, start=1):
        object_dir = mag_dir / f"Object_{index}"
        object_dir.mkdir(parents=True, exist_ok=True)

        positions = sample_sphere(sphere.radius, sphere.spacing)
        moments = np.repeat(np.asarray(sphere.magnetization, dtype=float)[None, :], len(positions), axis=0)
        data = np.column_stack((positions, moments))

        np.savetxt(
            object_dir / f"m_{dataset_index}.csv",
            data,
            fmt="%.10e",
            delimiter=" ",
        )

        struct_rows.append(sphere.center)
        rot_rows.append(sphere.angles)

    np.savetxt(real_space_dir / "StructData.csv", np.asarray(struct_rows), fmt="%.10e", delimiter=" ")
    np.savetxt(real_space_dir / "RotData.csv", np.asarray(rot_rows), fmt="%.10e", delimiter=" ")

    return real_space_dir


def clear_uniform_sphere_case(output_dir: Path) -> None:
    """Delete the generated RealSpaceData directory for this example."""

    real_space_dir = Path(output_dir) / "RealSpaceData"
    if real_space_dir.exists():
        shutil.rmtree(real_space_dir)


def example_two_spheres() -> list[UniformSphere]:
    """Return a simple two-sphere setup for manual tests."""

    return [
        UniformSphere(
            radius=10.0,
            spacing=1.0,
            magnetization=(0.0, 0.0, 1.0),
            center=(-14.0, 0.0, 0.0),
            angles=(0.0, 0.0, 0.0),
        ),
        UniformSphere(
            radius=7.0,
            spacing=1.5,
            magnetization=(0.0, 0.0, 0.8),
            center=(14.0, 3.0, 0.0),
            angles=(0.4, 0.8, 0.2),
        ),
    ]
