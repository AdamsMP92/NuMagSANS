import numpy as np
import pytest

from NuMagSANS.SystemDesigner.AssemblyBaseMagnetizer import evaluate_vector_field
from NuMagSANS.SystemDesigner.AssemblyBaseMagnetizer.VectorFieldModels.OperatorKernelVectorFieldLib import (
    GaussianKernel,
    SechKernel,
    build_operator_kernel_model,
)
from NuMagSANS.SystemDesigner.AssemblyBaseMagnetizer.VectorFieldModels.SphericalVectorFieldLib import (
    unit_field as spherical_unit_field,
)


def _sample_coordinates():
    x = np.array([-1.0, 0.0, 1.0, 0.5])
    y = np.array([0.0, 1.0, 0.0, -0.5])
    z = np.array([0.0, 0.0, 1.0, 0.25])
    return x, y, z


def _assert_unit_vectors(mx, my, mz):
    field = np.column_stack((mx, my, mz))
    assert field.shape[1] == 3
    assert np.all(np.isfinite(field))
    assert np.allclose(np.linalg.norm(field, axis=1), 1.0)


def test_spherical_vortex_unit_field_is_finite_and_normalized():
    x, y, z = _sample_coordinates()
    params = {
        "field_type": "vortex",
        "profile_type": "linear",
        "xi_type": "cylindrical_xi",
        "kappa": 1.0,
    }

    mx, my, mz = spherical_unit_field(x, y, z, D=2.0, params=params)

    _assert_unit_vectors(mx, my, mz)


def test_longitudinal_helix_is_normalized_after_adding_lambda_component():
    x, y, z = _sample_coordinates()
    params = {
        "field_type": "longitudinal_helix",
        "k": 0.7,
        "chirality": 1.0,
        "lambda": 0.5,
    }

    mx, my, mz = spherical_unit_field(x, y, z, D=2.0, params=params)

    _assert_unit_vectors(mx, my, mz)
    assert np.allclose(mz, 0.5 / np.sqrt(1.25))


def test_longitudinal_helix_normalization_can_be_disabled():
    x, y, z = _sample_coordinates()
    params = {
        "field_type": "longitudinal_helix",
        "k": 0.7,
        "chirality": 1.0,
        "lambda_z": 0.5,
        "normalize": False,
    }

    mx, my, mz = spherical_unit_field(x, y, z, D=2.0, params=params)
    norm = np.sqrt(mx**2 + my**2 + mz**2)

    assert np.allclose(mz, 0.5)
    assert np.allclose(norm, np.sqrt(1.25))


def test_spherical_legacy_additional_normalize_key_is_still_accepted():
    x, y, z = _sample_coordinates()
    params = {
        "field_type": "longitudinal_helix",
        "k": 0.7,
        "chirality": 1.0,
        "lambda_z": 0.5,
        "additional_normalize": False,
    }

    mx, my, mz = spherical_unit_field(x, y, z, D=2.0, params=params)
    norm = np.sqrt(mx**2 + my**2 + mz**2)

    assert np.allclose(norm, np.sqrt(1.25))


def test_operator_kernel_gaussian_field_is_finite_and_normalized():
    x, y, z = _sample_coordinates()
    params = {
        "library": "operator_kernel",
        "kernel_type": "gaussian",
        "sigma": [0.4, 0.4, 0.4],
        "component_specs": ["-a * dy", "a * dx", "c"],
        "operator_parameters": {"a": 0.8, "c": 0.55},
        "normalize": True,
    }

    mx, my, mz = evaluate_vector_field(x, y, z, D=2.0, params=params)

    _assert_unit_vectors(mx, my, mz)


def test_operator_kernel_sech_field_accepts_zero_component():
    x, y, z = _sample_coordinates()
    params = {
        "library": "operator_kernel",
        "kernel_type": "sech",
        "a": [2.0, 2.0, 2.4],
        "component_specs": ["0", "a * dx", "c"],
        "operator_parameters": {"a": 0.8, "c": 0.55},
        "normalize": False,
    }

    mx, my, mz = evaluate_vector_field(x, y, z, D=2.0, params=params)

    assert np.allclose(mx, 0.0)
    assert np.all(np.isfinite(my))
    assert np.all(np.isfinite(mz))


def test_vector_field_registry_defaults_to_spherical_model():
    x, y, z = _sample_coordinates()
    params = {
        "field_type": "vortex",
        "profile_type": "linear",
        "xi_type": "cylindrical_xi",
        "kappa": 1.0,
    }

    registry_field = evaluate_vector_field(x, y, z, D=2.0, params=params)
    spherical_field = spherical_unit_field(x, y, z, D=2.0, params=params)

    assert all(np.allclose(a, b) for a, b in zip(registry_field, spherical_field))


def test_vector_field_registry_rejects_unknown_library():
    x, y, z = _sample_coordinates()

    with pytest.raises(ValueError, match="Unknown vector-field library"):
        evaluate_vector_field(x, y, z, D=2.0, params={"library": "unknown"})


def test_operator_kernel_model_real_and_fourier_shapes():
    kernel = GaussianKernel(np.diag([0.4, 0.5, 0.6]) ** 2)
    model = build_operator_kernel_model(
        kernel=kernel,
        component_specs=["-a * dy", "a * dx", "c"],
        operator_parameters={"a": 0.8, "c": 0.55},
    )
    points = np.array([[0.0, 0.0, 0.0], [0.2, -0.1, 0.3]])

    real_field = model.real(points)
    fourier_field = model.fourier(points)

    assert real_field.shape == (2, 3)
    assert fourier_field.shape == (2, 3)
    assert np.all(np.isfinite(real_field))
    assert np.all(np.isfinite(fourier_field))


def test_kernel_parameter_validation():
    with pytest.raises(ValueError, match="Sigma must have shape"):
        GaussianKernel(np.array([1.0, 2.0]))

    with pytest.raises(ValueError, match="a must contain three"):
        SechKernel([1.0, 2.0])
