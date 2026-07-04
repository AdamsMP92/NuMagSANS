import math

import numpy as np
import pytest

from NuMagSANS.SystemDesigner.AssemblyBaseOrchestrator import (
    rotations_to_array,
    sobol_zyz_rotations,
    write_rotdata_loop,
)


def test_sobol_zyz_rotations_shape_ranges_and_reproducibility():
    rotations_a = sobol_zyz_rotations(n_samples=7, beta_min=0.2, beta_max=1.1, seed=42)
    rotations_b = sobol_zyz_rotations(n_samples=7, beta_min=0.2, beta_max=1.1, seed=42)
    array = rotations_to_array(rotations_a)

    assert array.shape == (7, 3)
    assert np.allclose(array, rotations_to_array(rotations_b))
    assert np.all((0.0 <= array[:, 0]) & (array[:, 0] <= 2.0 * math.pi))
    assert np.all((0.2 <= array[:, 1]) & (array[:, 1] <= 1.1))
    assert np.all((0.0 <= array[:, 2]) & (array[:, 2] <= 2.0 * math.pi))


def test_sobol_zyz_rotations_rejects_invalid_input():
    with pytest.raises(ValueError, match="n_samples"):
        sobol_zyz_rotations(n_samples=0)

    with pytest.raises(ValueError, match="beta_min"):
        sobol_zyz_rotations(n_samples=4, beta_min=-0.1)

    with pytest.raises(ValueError, match="beta_max"):
        sobol_zyz_rotations(n_samples=4, beta_min=1.0, beta_max=0.5)


def test_rotations_to_array_validation():
    with pytest.raises(KeyError, match="alpha, beta, and gamma"):
        rotations_to_array({"alpha": [0.0], "beta": [0.0]})

    with pytest.raises(ValueError, match="matching shapes"):
        rotations_to_array({"alpha": [0.0, 1.0], "beta": [0.0], "gamma": [0.0]})


def test_write_rotdata_loop(tmp_path):
    rotations = [
        sobol_zyz_rotations(n_samples=3, seed=1),
        np.zeros((3, 3)),
    ]

    summary = write_rotdata_loop(rotations, tmp_path)

    assert summary["n_rotdata"] == 2
    assert len(summary["output_files"]) == 2
    for filename in summary["output_files"]:
        data = np.loadtxt(filename)
        assert data.shape == (3, 3)


def test_write_rotdata_loop_rejects_wrong_shape(tmp_path):
    with pytest.raises(ValueError, match="shape"):
        write_rotdata_loop([np.zeros((3, 2))], tmp_path)
