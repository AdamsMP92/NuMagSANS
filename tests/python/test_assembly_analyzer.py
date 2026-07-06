import numpy as np
import pytest

from NuMagSANS.SystemDesigner.AssemblyAnalyzer import (
    downsample_vectors,
    evaluate_color_field,
    load_magnetization_file,
)
from NuMagSANS.SystemDesigner.AssemblyAnalyzer.MagnetizationPlot import (
    _cut_mask,
    _save_plotter_png,
    _scalar_clim,
)


def test_load_magnetization_file(tmp_path):
    filename = tmp_path / "m_1.csv"
    np.savetxt(
        filename,
        np.array(
            [
                [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                [1.0, 0.0, 0.0, 0.0, 1.0, 0.0],
            ]
        ),
    )

    loaded = load_magnetization_file(filename)

    assert loaded["positions"].shape == (2, 3)
    assert loaded["vectors"].shape == (2, 3)
    assert loaded["data"].shape == (2, 6)


def test_load_magnetization_file_rejects_wrong_column_count(tmp_path):
    filename = tmp_path / "bad.csv"
    np.savetxt(filename, np.array([[0.0, 1.0, 2.0]]))

    with pytest.raises(ValueError, match="Expected six columns"):
        load_magnetization_file(filename)


def test_downsample_vectors_stride_and_random_selection():
    positions = np.arange(30, dtype=float).reshape(10, 3)
    vectors = positions + 100.0

    stride_positions, stride_vectors = downsample_vectors(positions, vectors, step=2)
    assert np.array_equal(stride_positions, positions[::2])
    assert np.array_equal(stride_vectors, vectors[::2])

    sample_positions, sample_vectors = downsample_vectors(positions, vectors, max_vectors=4, seed=1)
    assert sample_positions.shape == (4, 3)
    assert sample_vectors.shape == (4, 3)


def test_downsample_vectors_validation():
    with pytest.raises(ValueError, match="positions"):
        downsample_vectors(np.zeros((3, 2)), np.zeros((3, 3)))

    with pytest.raises(ValueError, match="same number"):
        downsample_vectors(np.zeros((2, 3)), np.zeros((3, 3)))

    with pytest.raises(ValueError, match="step"):
        downsample_vectors(np.zeros((2, 3)), np.zeros((2, 3)), step=0)


def test_evaluate_color_field_builtin_and_callable_components():
    positions = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    vectors = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]])

    assert np.allclose(evaluate_color_field("mz", positions, vectors)[1], [3.0, 6.0, 9.0])
    assert np.allclose(evaluate_color_field("m_rho", positions, vectors)[1], [1.0, 5.0, 0.0])
    assert np.allclose(evaluate_color_field("m_phi", positions, vectors)[1], [2.0, -4.0, 0.0])
    assert np.allclose(
        evaluate_color_field(lambda pos, vec: pos[:, 2] + vec[:, 0], positions, vectors)[1],
        [1.0, 4.0, 8.0],
    )


def test_evaluate_color_field_validation():
    positions = np.zeros((3, 3))
    vectors = np.zeros((3, 3))

    with pytest.raises(ValueError, match="Unknown color field"):
        evaluate_color_field("does_not_exist", positions, vectors)

    with pytest.raises(ValueError, match="one-dimensional"):
        evaluate_color_field(lambda pos, vec: np.zeros((3, 2)), positions, vectors)


def test_cut_mask_and_scalar_clim_helpers():
    positions = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]])

    assert np.array_equal(_cut_mask(positions, {"x": 0.5}, "upper"), [True, False, True, True])
    assert np.array_equal(_cut_mask(positions, {"z": 1.0}, "lower"), [False, False, False, True])
    assert _scalar_clim(np.array([1.0, 1.0])) == (0.95, 1.05)

    with pytest.raises(ValueError, match="cut_mode"):
        _cut_mask(positions, {"x": 0.0}, "bad")


def test_save_plotter_png_uses_screenshot_image(tmp_path):
    class FakePlotter:
        def screenshot(self, return_img=False):
            assert return_img is True
            return np.zeros((4, 5, 3), dtype=np.uint8)

    output = _save_plotter_png(FakePlotter(), tmp_path / "nested" / "plot.png")

    assert output == tmp_path / "nested" / "plot.png"
    assert output.exists()
