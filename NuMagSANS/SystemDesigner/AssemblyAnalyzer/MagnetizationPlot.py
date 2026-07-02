"""Local plotting helpers for materialized magnetization files.

These functions are intended for interactive inspection on a local workstation.
They are not used by the SystemDesigner workflows and should not be part of
standard HPC production runs.
"""

from pathlib import Path

import numpy as np


def load_magnetization_file(filename):
    """Load a NuMagSANS magnetization file with columns x y z mx my mz."""
    filename = Path(filename)
    if not filename.exists():
        raise FileNotFoundError(f"Magnetization file does not exist: {filename}")

    data = np.loadtxt(filename)
    if data.ndim == 1:
        data = data.reshape(1, -1)

    if data.shape[1] != 6:
        raise ValueError(f"Expected six columns x y z mx my mz in {filename}.")

    return {
        "positions": data[:, :3],
        "vectors": data[:, 3:6],
        "data": data,
    }


def downsample_vectors(positions, vectors, step=1, max_vectors=None, seed=None):
    """Downsample vector-field data by stride or random selection.

    Parameters
    ----------
    positions, vectors : array-like
        Arrays with shape ``(n_points, 3)``.
    step : int, optional
        Keep every ``step``-th vector before optional random sub-selection.
    max_vectors : int, optional
        If provided, randomly select at most this many vectors after striding.
    seed : int, optional
        Random seed used when ``max_vectors`` triggers random sub-selection.
    """
    positions = np.asarray(positions, dtype=float)
    vectors = np.asarray(vectors, dtype=float)

    if positions.ndim != 2 or positions.shape[1] != 3:
        raise ValueError("positions must have shape (n_points, 3).")
    if vectors.ndim != 2 or vectors.shape[1] != 3:
        raise ValueError("vectors must have shape (n_points, 3).")
    if len(positions) != len(vectors):
        raise ValueError("positions and vectors must contain the same number of rows.")

    step = int(step)
    if step <= 0:
        raise ValueError("step must be a positive integer.")

    positions = positions[::step]
    vectors = vectors[::step]

    if max_vectors is not None:
        max_vectors = int(max_vectors)
        if max_vectors <= 0:
            raise ValueError("max_vectors must be a positive integer.")

        if len(positions) > max_vectors:
            rng = np.random.default_rng(seed)
            indices = np.sort(rng.choice(len(positions), size=max_vectors, replace=False))
            positions = positions[indices]
            vectors = vectors[indices]

    return positions, vectors


def _cylindrical_projection(positions, vectors, component):
    """Return cylindrical magnetization components around the z-axis."""
    rho = np.hypot(positions[:, 0], positions[:, 1])
    values = np.zeros(len(positions), dtype=float)
    mask = rho > 0.0

    if component == "m_rho":
        values[mask] = (vectors[mask, 0] * positions[mask, 0] + vectors[mask, 1] * positions[mask, 1]) / rho[mask]
    elif component == "m_phi":
        values[mask] = (-vectors[mask, 0] * positions[mask, 1] + vectors[mask, 1] * positions[mask, 0]) / rho[mask]
    else:
        raise ValueError(f"Unknown cylindrical component: {component}")

    return values


def _evaluate_builtin_color_field(color_by, positions, vectors):
    """Evaluate a named scalar field for coloring magnetization glyphs."""
    color_by = color_by.lower()

    if color_by == "x":
        return positions[:, 0]
    if color_by == "y":
        return positions[:, 1]
    if color_by == "z":
        return positions[:, 2]
    if color_by == "rho":
        return np.hypot(positions[:, 0], positions[:, 1])
    if color_by == "r":
        return np.linalg.norm(positions, axis=1)
    if color_by == "mx":
        return vectors[:, 0]
    if color_by == "my":
        return vectors[:, 1]
    if color_by == "mz":
        return vectors[:, 2]
    if color_by in {"magnitude", "m_abs", "m"}:
        return np.linalg.norm(vectors, axis=1)
    if color_by in {"m_rho", "m_phi"}:
        return _cylindrical_projection(positions, vectors, color_by)

    raise ValueError(
        f"Unknown color field '{color_by}'. Supported names are x, y, z, rho, r, "
        "mx, my, mz, magnitude, m_abs, m, m_rho, and m_phi."
    )


def evaluate_color_field(color_by, positions, vectors):
    """Evaluate a scalar color field from a name or callable.

    Parameters
    ----------
    color_by : str or callable or None
        Built-in scalar field name or a callable with signature
        ``color_by(positions, vectors)``. If ``None``, no scalar field is used.
    positions, vectors : array-like
        Arrays with shape ``(n_points, 3)``.

    Returns
    -------
    tuple
        ``(field_name, values)``. Both entries are ``None`` if ``color_by`` is
        ``None``.
    """
    if color_by is None:
        return None, None

    positions = np.asarray(positions, dtype=float)
    vectors = np.asarray(vectors, dtype=float)

    if isinstance(color_by, str):
        field_name = color_by
        values = _evaluate_builtin_color_field(color_by, positions, vectors)
    elif callable(color_by):
        field_name = getattr(color_by, "__name__", "color_field")
        values = color_by(positions, vectors)
    else:
        raise TypeError("color_by must be None, a built-in field name, or a callable.")

    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or len(values) != len(positions):
        raise ValueError("The color field must return a one-dimensional array with one value per vector.")

    return field_name, values


def plot_magnetization_file(
    filename,
    vector_scale=1.0,
    center_vectors=False,
    arrow_tip_length=None,
    arrow_tip_radius=None,
    arrow_tip_resolution=None,
    arrow_shaft_radius=None,
    arrow_shaft_resolution=None,
    downsample_step=1,
    max_vectors=None,
    seed=None,
    show_points=True,
    point_size=4.0,
    vector_color="tomato",
    point_color="black",
    color_by=None,
    cmap="viridis",
    show_scalar_bar=True,
    scalar_bar_title=None,
    background="white",
    window_size=(1100, 850),
    show=True,
    return_plotter=False,
):
    """Plot one magnetization file as a PyVista vector-field glyph plot.

    Parameters
    ----------
    filename : str or pathlib.Path
        Magnetization file with whitespace-separated columns ``x y z mx my mz``.
    vector_scale : float, optional
        Uniform glyph length scale.
    center_vectors : bool, optional
        If ``True``, center each arrow glyph around its ``x y z`` coordinate.
        If ``False``, use the default PyVista arrow anchoring.
    arrow_tip_length, arrow_tip_radius : float, optional
        Arrow head length and radius. If any arrow-shape parameter is set, a
        custom PyVista arrow geometry is used.
    arrow_tip_resolution : int, optional
        Angular resolution of the arrow head.
    arrow_shaft_radius : float, optional
        Radius of the arrow shaft.
    arrow_shaft_resolution : int, optional
        Angular resolution of the arrow shaft.
    downsample_step : int, optional
        Plot every ``downsample_step``-th vector.
    max_vectors : int, optional
        Randomly keep at most this many vectors after stride downsampling.
    seed : int, optional
        Seed for random vector sub-selection.
    show_points : bool, optional
        If ``True``, also show the underlying point cloud.
    point_size : float, optional
        Point size used when ``show_points`` is active.
    vector_color, point_color : str, optional
        PyVista-compatible colors.
    color_by : str or callable, optional
        Scalar field used to color the vectors. Built-in names include ``mz``,
        ``m_phi``, ``m_rho``, ``magnitude``, ``rho``, and coordinate/vector
        components. A custom callable must use ``color_by(positions, vectors)``
        and return one scalar per vector.
    cmap : str, optional
        PyVista/Matplotlib color map used when ``color_by`` is active.
    show_scalar_bar : bool, optional
        If ``True``, show a scalar bar for the vector coloring.
    scalar_bar_title : str, optional
        Custom scalar bar title. Defaults to the selected color field name.
    background : str, optional
        Plot background color.
    window_size : tuple, optional
        PyVista render window size.
    show : bool, optional
        If ``True``, open the interactive PyVista window.
    return_plotter : bool, optional
        If ``True``, return the PyVista plotter after adding all actors.

    Returns
    -------
    pyvista.Plotter or None
        Plotter if ``return_plotter`` is true, otherwise ``None``.
    """
    try:
        import pyvista as pv
    except ImportError as error:
        raise ImportError(
            "plot_magnetization_file requires pyvista. Install it locally with "
            "`pip install pyvista` or install NuMagSANS with a visualization extra."
        ) from error

    loaded = load_magnetization_file(filename)
    positions, vectors = downsample_vectors(
        loaded["positions"],
        loaded["vectors"],
        step=downsample_step,
        max_vectors=max_vectors,
        seed=seed,
    )

    mesh = pv.PolyData(positions)
    mesh["magnetization"] = vectors
    mesh["mx"] = vectors[:, 0]
    mesh["my"] = vectors[:, 1]
    mesh["mz"] = vectors[:, 2]
    mesh["magnitude"] = np.linalg.norm(vectors, axis=1)

    color_field_name, color_values = evaluate_color_field(color_by, positions, vectors)
    if color_field_name is not None:
        mesh[color_field_name] = color_values

    arrow_kwargs = {
        "tip_length": arrow_tip_length,
        "tip_radius": arrow_tip_radius,
        "tip_resolution": arrow_tip_resolution,
        "shaft_radius": arrow_shaft_radius,
        "shaft_resolution": arrow_shaft_resolution,
    }
    arrow_kwargs = {key: value for key, value in arrow_kwargs.items() if value is not None}
    use_custom_arrow = center_vectors or bool(arrow_kwargs)

    if use_custom_arrow:
        arrow = pv.Arrow(
            start=(-0.5, 0.0, 0.0) if center_vectors else (0.0, 0.0, 0.0),
            direction=(1.0, 0.0, 0.0),
            scale=1.0,
            **arrow_kwargs,
        )
        glyphs = mesh.glyph(
            orient="magnetization",
            scale=False,
            factor=float(vector_scale),
            geom=arrow,
        )
    else:
        glyphs = mesh.glyph(orient="magnetization", scale=False, factor=float(vector_scale))

    plotter = pv.Plotter(window_size=window_size)
    plotter.set_background(background)

    if show_points:
        plotter.add_mesh(mesh, color=point_color, point_size=point_size, render_points_as_spheres=True)

    if color_field_name is None:
        plotter.add_mesh(glyphs, color=vector_color)
    else:
        plotter.add_mesh(
            glyphs,
            scalars=color_field_name,
            cmap=cmap,
            show_scalar_bar=show_scalar_bar,
            scalar_bar_args={"title": scalar_bar_title or color_field_name},
        )

    plotter.add_axes()
    plotter.show_bounds(grid="front", location="outer")

    if show:
        plotter.show()

    if return_plotter:
        return plotter

    return None
