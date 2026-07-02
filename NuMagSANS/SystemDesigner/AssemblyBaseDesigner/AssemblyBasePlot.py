from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def _select_parameters(analysis, parameters):
    """Return the parameter names that should be plotted."""
    available = list(analysis["parameters"].keys())
    if parameters is None:
        return available

    missing = [parameter for parameter in parameters if parameter not in available]
    if missing:
        raise KeyError(f"Unknown parameter names: {missing}")

    return list(parameters)


def _plot_numeric_histogram(ax, values, parameter_name, bins):
    """Plot a histogram for one numeric parameter."""
    values = np.asarray(values, dtype=float)
    ax.hist(
        values,
        bins=bins,
        density=True,
        color="#4a76b8",
        edgecolor="#202020",
        alpha=0.65,
        label="realized",
    )
    ax.axvline(np.mean(values), color="#b83232", linewidth=1.4, label="mean")
    ax.set_xlabel(parameter_name)
    ax.set_ylabel("density")
    ax.legend(loc="best")


def _normal_pdf(x, mean, sigma):
    """Evaluate a normal probability density."""
    if sigma <= 0:
        return np.zeros_like(x)

    return 1.0 / (sigma * np.sqrt(2.0 * np.pi)) * np.exp(-0.5 * ((x - mean) / sigma) ** 2)


def _lognormal_pdf(x, mean, sigma, mean_type):
    """Evaluate a log-normal probability density."""
    if sigma <= 0:
        return np.zeros_like(x)

    log_mean = mean
    if mean_type == "arithmetic":
        log_mean = np.log(mean) - 0.5 * sigma**2

    pdf = np.zeros_like(x)
    positive = x > 0
    pdf[positive] = (
        1.0
        / (x[positive] * sigma * np.sqrt(2.0 * np.pi))
        * np.exp(-0.5 * ((np.log(x[positive]) - log_mean) / sigma) ** 2)
    )
    return pdf


def _uniform_pdf(x, low, high):
    """Evaluate a uniform probability density."""
    pdf = np.zeros_like(x)
    if high <= low:
        return pdf

    pdf[(x >= low) & (x <= high)] = 1.0 / (high - low)
    return pdf


def _plot_expected_density(ax, values, distribution_spec):
    """Overlay the expected density implied by a distribution specification."""
    if not distribution_spec:
        return

    distribution = distribution_spec.get("distribution")
    values = np.asarray(values, dtype=float)
    value_min = np.min(values)
    value_max = np.max(values)
    value_span = value_max - value_min
    if value_span == 0:
        value_span = max(abs(value_min), 1.0)

    x_min = value_min - 0.15 * value_span
    x_max = value_max + 0.15 * value_span

    if distribution == "normal":
        mean = float(distribution_spec["mean"])
        sigma = float(distribution_spec["sigma"])
        x_min = min(x_min, mean - 4.0 * sigma)
        x_max = max(x_max, mean + 4.0 * sigma)
        x = np.linspace(x_min, x_max, 400)
        y = _normal_pdf(x, mean, sigma)
    elif distribution == "lognormal":
        mean = float(distribution_spec["mean"])
        sigma = float(distribution_spec["sigma"])
        mean_type = distribution_spec.get("mean_type", "arithmetic")
        x_min = max(0.0, x_min)
        x = np.linspace(x_min, x_max, 400)
        y = _lognormal_pdf(x, mean, sigma, mean_type)
    elif distribution == "uniform":
        low = float(distribution_spec["low"])
        high = float(distribution_spec["high"])
        x_min = min(x_min, low)
        x_max = max(x_max, high)
        x = np.linspace(x_min, x_max, 400)
        y = _uniform_pdf(x, low, high)
    elif distribution == "constant":
        ax.axvline(
            float(distribution_spec["value"]),
            color="#202020",
            linewidth=1.4,
            linestyle="--",
            label="expected",
        )
        ax.legend(loc="best")
        return
    else:
        return

    ax.plot(x, y, color="#202020", linewidth=1.8, linestyle="--", label="expected")
    ax.legend(loc="best")


def _plot_categorical_counts(ax, values, parameter_name):
    """Plot category counts for one non-numeric parameter."""
    values = np.asarray(values, dtype=str)
    labels, counts = np.unique(values, return_counts=True)
    ax.bar(labels, counts, color="#4a8f68", edgecolor="#202020", alpha=0.85)
    ax.set_xlabel(parameter_name)
    ax.set_ylabel("count")


def plot_parameter_distributions(
    analysis,
    output_path=None,
    parameters=None,
    bins=20,
    figsize_per_panel=(4.0, 3.2),
):
    """Plot realized template-parameter distributions from analyzer output.

    Numeric parameters are shown as histograms. Non-numeric parameters, such as
    atom type labels, are shown as category-count bar plots.

    Parameters
    ----------
    analysis : dict
        Dictionary returned by ``analyze_assembly_base_dataset``.
    output_path : str or pathlib.Path, optional
        If provided, save the figure to this path.
    parameters : iterable of str, optional
        Parameter names to plot. By default all extracted parameters are shown.
    bins : int, optional
        Number of histogram bins for numeric parameters.
    figsize_per_panel : tuple, optional
        Width and height per subplot in inches.

    Returns
    -------
    matplotlib.figure.Figure
        The created Matplotlib figure.
    """
    parameter_names = _select_parameters(analysis, parameters)
    if not parameter_names:
        raise ValueError("No template parameters available for plotting.")

    n_panels = len(parameter_names)
    fig_width = figsize_per_panel[0] * n_panels
    fig_height = figsize_per_panel[1]
    fig, axes = plt.subplots(1, n_panels, figsize=(fig_width, fig_height))

    if n_panels == 1:
        axes = [axes]

    numeric_parameters = analysis["numeric_parameters"]
    parameter_distributions = analysis.get("parameter_distributions", {})
    for ax, parameter_name in zip(axes, parameter_names):
        if parameter_name in numeric_parameters:
            _plot_numeric_histogram(
                ax,
                numeric_parameters[parameter_name],
                parameter_name,
                bins,
            )
            _plot_expected_density(
                ax,
                numeric_parameters[parameter_name],
                parameter_distributions.get(parameter_name),
            )
        else:
            _plot_categorical_counts(
                ax,
                analysis["parameters"][parameter_name],
                parameter_name,
            )

        ax.set_title(parameter_name)

    fig.suptitle(f"Assembly base parameter distributions, n={analysis['n_objects']}")
    fig.tight_layout()

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=300)

    return fig
