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


def _mesh_from_vectors(pv, positions, vectors, color_field_name=None, color_values=None):
    """Create a PyVista point mesh with magnetization vectors and optional scalars."""
    mesh = pv.PolyData(positions)
    mesh["magnetization"] = vectors
    mesh["mx"] = vectors[:, 0]
    mesh["my"] = vectors[:, 1]
    mesh["mz"] = vectors[:, 2]
    mesh["magnitude"] = np.linalg.norm(vectors, axis=1)

    if color_field_name is not None:
        mesh[color_field_name] = color_values

    return mesh


def _glyph_from_mesh(mesh, pv, vector_scale, center_vectors=False, arrow_kwargs=None):
    """Create vector glyphs from a point mesh."""
    if mesh.n_points == 0:
        return pv.PolyData()

    arrow_kwargs = arrow_kwargs or {}
    use_custom_arrow = center_vectors or bool(arrow_kwargs)

    if use_custom_arrow:
        arrow = pv.Arrow(
            start=(-0.5, 0.0, 0.0) if center_vectors else (0.0, 0.0, 0.0),
            direction=(1.0, 0.0, 0.0),
            scale=1.0,
            **arrow_kwargs,
        )
        return mesh.glyph(
            orient="magnetization",
            scale=False,
            factor=float(vector_scale),
            geom=arrow,
        )

    return mesh.glyph(orient="magnetization", scale=False, factor=float(vector_scale))


def _axis_index(axis):
    """Return the Cartesian column index for an axis name."""
    axis = str(axis).lower()
    if axis == "x":
        return 0
    if axis == "y":
        return 1
    if axis == "z":
        return 2
    raise ValueError("cut axes must be selected from 'x', 'y', and 'z'.")


def _cut_mask(positions, cut_values, cut_mode):
    """Return a boolean mask for Cartesian cut-slider values."""
    mask = np.ones(len(positions), dtype=bool)
    for axis, value in cut_values.items():
        coordinate = positions[:, _axis_index(axis)]
        if cut_mode == "upper":
            mask &= coordinate <= value
        elif cut_mode == "lower":
            mask &= coordinate >= value
        else:
            raise ValueError("cut_mode must be either 'upper' or 'lower'.")
    return mask


def _scalar_clim(values):
    """Return a non-degenerate scalar color range."""
    value_min = float(np.min(values))
    value_max = float(np.max(values))
    if value_min == value_max:
        pad = 1.0 if value_min == 0.0 else abs(value_min) * 0.05
        return value_min - pad, value_max + pad
    return value_min, value_max


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
    enable_cut_sliders=False,
    cut_axes=("x", "y", "z"),
    cut_mode="upper",
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
    enable_cut_sliders : bool, optional
        If ``True``, add interactive Cartesian cut sliders. In the default
        ``cut_mode="upper"``, each slider shows points with coordinate values
        below the current threshold. The initial slider positions show the full
        vector field.
    cut_axes : tuple of str, optional
        Cartesian axes for which sliders are created. Allowed entries are
        ``"x"``, ``"y"``, and ``"z"``.
    cut_mode : {"upper", "lower"}, optional
        In ``"upper"`` mode, show points with coordinate ``<=`` threshold.
        In ``"lower"`` mode, show points with coordinate ``>=`` threshold.
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

    color_field_name, color_values = evaluate_color_field(color_by, positions, vectors)

    arrow_kwargs = {
        "tip_length": arrow_tip_length,
        "tip_radius": arrow_tip_radius,
        "tip_resolution": arrow_tip_resolution,
        "shaft_radius": arrow_shaft_radius,
        "shaft_resolution": arrow_shaft_resolution,
    }
    arrow_kwargs = {key: value for key, value in arrow_kwargs.items() if value is not None}

    plotter = pv.Plotter(window_size=window_size)
    plotter.set_background(background)

    active_cut_axes = tuple(dict.fromkeys(str(axis).lower() for axis in cut_axes))
    for axis in active_cut_axes:
        _axis_index(axis)

    if cut_mode not in {"upper", "lower"}:
        raise ValueError("cut_mode must be either 'upper' or 'lower'.")

    cut_values = {}
    for axis in active_cut_axes:
        axis_values = positions[:, _axis_index(axis)]
        cut_values[axis] = float(np.max(axis_values) if cut_mode == "upper" else np.min(axis_values))

    actors = {
        "points": None,
        "glyphs": None,
    }

    def add_vector_actors(mask):
        active_positions = positions[mask]
        active_vectors = vectors[mask]
        active_color_values = color_values[mask] if color_values is not None else None
        mesh = _mesh_from_vectors(
            pv,
            active_positions,
            active_vectors,
            color_field_name=color_field_name,
            color_values=active_color_values,
        )
        glyphs = _glyph_from_mesh(
            mesh,
            pv,
            vector_scale=vector_scale,
            center_vectors=center_vectors,
            arrow_kwargs=arrow_kwargs,
        )

        if show_points:
            actors["points"] = plotter.add_mesh(
                mesh,
                color=point_color,
                point_size=point_size,
                render_points_as_spheres=True,
            )

        if color_field_name is None or glyphs.n_points == 0:
            actors["glyphs"] = plotter.add_mesh(glyphs, color=vector_color)
        else:
            actors["glyphs"] = plotter.add_mesh(
                glyphs,
                scalars=color_field_name,
                cmap=cmap,
                clim=_scalar_clim(color_values),
                show_scalar_bar=show_scalar_bar,
                scalar_bar_args={"title": scalar_bar_title or color_field_name},
            )

    def update_cut(axis, value):
        cut_values[axis] = float(value)
        for actor_name in ("points", "glyphs"):
            if actors[actor_name] is not None:
                plotter.remove_actor(actors[actor_name])
                actors[actor_name] = None

        add_vector_actors(_cut_mask(positions, cut_values, cut_mode))
        plotter.render()

    add_vector_actors(_cut_mask(positions, cut_values, cut_mode))

    if enable_cut_sliders:
        slider_positions = {
            "x": ((0.05, 0.92), (0.35, 0.92)),
            "y": ((0.38, 0.92), (0.68, 0.92)),
            "z": ((0.71, 0.92), (0.95, 0.92)),
        }
        for axis in active_cut_axes:
            axis_values = positions[:, _axis_index(axis)]
            value_range = (float(np.min(axis_values)), float(np.max(axis_values)))
            pointa, pointb = slider_positions.get(axis, ((0.05, 0.92), (0.35, 0.92)))
            plotter.add_slider_widget(
                lambda value, selected_axis=axis: update_cut(selected_axis, value),
                rng=value_range,
                value=cut_values[axis],
                title=f"{axis} cut",
                pointa=pointa,
                pointb=pointb,
            )

    plotter.add_axes()
    plotter.show_bounds(grid="front", location="outer")

    if show:
        plotter.show()

    if return_plotter:
        return plotter

    return None
