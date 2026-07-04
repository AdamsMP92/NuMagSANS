import matplotlib.pyplot as plt
import numpy as np


def set_equal_aspect_3d(ax, x, y, z):
    """Set equal visual scaling for all three axes in a 3D Matplotlib plot.

    Matplotlib's 3D axes do not automatically use the same scale for x, y, and
    z. This helper expands the smaller axis ranges to match the largest one and
    then requests an equal box aspect ratio.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        A 3D Matplotlib axis.
    x, y, z : array-like
        Coordinate arrays used to determine the axis limits.
    """
    x = np.asarray(x)
    y = np.asarray(y)
    z = np.asarray(z)

    x_center = 0.5 * (np.max(x) + np.min(x))
    y_center = 0.5 * (np.max(y) + np.min(y))
    z_center = 0.5 * (np.max(z) + np.min(z))

    max_range = max(
        np.max(x) - np.min(x),
        np.max(y) - np.min(y),
        np.max(z) - np.min(z),
    )
    half_range = 0.5 * max_range

    ax.set_xlim(x_center - half_range, x_center + half_range)
    ax.set_ylim(y_center - half_range, y_center + half_range)
    ax.set_zlim(z_center - half_range, z_center + half_range)
    ax.set_box_aspect((1, 1, 1))


def _atomtype_label(crystal):
    """Create a compact legend label from the atom types in a crystal."""
    if "atomtype" not in crystal:
        return "atoms"

    atomtypes = np.unique(crystal["atomtype"])
    if len(atomtypes) == 1:
        return str(atomtypes[0])

    return "mixed atom types"


def scatter_plot_crystal(crystal, title):
    """Plot crystal positions as a 3D scatter plot with equal axis scaling."""
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection="3d")

    ax.scatter(
        crystal["x"],
        crystal["y"],
        crystal["z"],
        s=36,
        c="#a83232",
        edgecolors="#2b2b2b",
        linewidths=0.35,
        alpha=0.9,
        label=_atomtype_label(crystal),
    )

    set_equal_aspect_3d(ax, crystal["x"], crystal["y"], crystal["z"])

    ax.set_title(title)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.legend(loc="upper right")

    plt.tight_layout()
    plt.show()
