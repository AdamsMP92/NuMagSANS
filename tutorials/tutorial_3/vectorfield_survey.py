import sys
from pathlib import Path

import numpy as np

BASE_DIR = Path(__file__).resolve().parent
REPO_ROOT = BASE_DIR.parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from NuMagSANS.SystemDesigner.AssemblyAnalyzer import plot_magnetization_file  # noqa: E402
from NuMagSANS.SystemDesigner.Workflows import write_spherical_replication_vectorfield_sweep  # noqa: E402

OUTPUT_DIR = BASE_DIR / "VectorFieldSurvey"
IMAGE_DIR = OUTPUT_DIR / "images"

# Set this to an integer, for example 3, to inspect only m_3.csv.
# Set it to "all" to walk through all generated vector-field variants.
PLOT_SELECTION = "all"

# Set to True for local interactive PyVista windows. The default writes PNG
# diagnostics, which also works for headless tutorial/smoke-test runs.
INTERACTIVE = False

FIELD_PARAMETER_CASES = [
    {
        "label": "hyperbolic vortex",
        "field_type": "vortex",
        "profile_type": "hyperbolic",
        "xi_type": "cylindrical_xi",
        "kappa": 6.0,
        "turns": 1.0,
    },
    {
        "label": "hyperbolic hedgehog",
        "field_type": "hedgehog",
        "profile_type": "hyperbolic",
        "xi_type": "spherical_xi",
        "kappa": 6.0,
        "turns": 1.0,
    },
    {
        "label": "hyperbolic artichoke",
        "field_type": "artichoke",
        "profile_type": "hyperbolic",
        "xi_type": "spherical_xi",
        "kappa": 6.0,
        "turns": 1.0,
    },
    {
        "label": "hyperbolic poloidal vortex",
        "field_type": "poloidal_vortex",
        "profile_type": "hyperbolic",
        "xi_type": "spherical_xi",
        "kappa": 6.0,
        "turns": 1.0,
    },
    {
        "label": "skyrmion",
        "field_type": "skyrmion",
        "xi_type": "cylindrical_xi",
        "N": 1,
        "m": 1,
        "gamma": 0.0,
    },
    {
        "label": "transversal helix",
        "field_type": "transversal_helix",
        "k": 0.6,
        "chirality": 1.0,
        "phase": 0.0,
    },
    {
        "label": "longitudinal helix",
        "field_type": "longitudinal_helix",
        "k": 0.6,
        "chirality": 1.0,
        "lambda": 0.5,
        "phase": 0.0,
    },
    {
        "label": "linearized vortex",
        "field_type": "linearized_vortex",
        "m_0": 1.0,
        "m_1": 0.05,
    },
    {
        "label": "operator-kernel Gaussian vortex",
        "library": "operator_kernel",
        "kernel_type": "gaussian",
        "sigma": [0.35, 0.35, 0.35],
        "component_specs": ["-a * dy", "a * dx", "c"],
        "operator_parameters": {"a": 0.8, "c": 0.55},
        "normalize": True,
    },
    {
        "label": "operator-kernel Gaussian skyrmion-like",
        "library": "operator_kernel",
        "kernel_type": "gaussian",
        "sigma": [0.42, 0.42, 0.32],
        "component_specs": ["-a * dy", "a * dx", "c + d * (dx^2 + dy^2)"],
        "operator_parameters": {"a": 0.7, "c": -1.6, "d": -0.5},
        "normalize": True,
    },
    {
        "label": "operator-kernel Sech vortex",
        "library": "operator_kernel",
        "kernel_type": "sech",
        "a": [2.0, 2.0, 2.4],
        "component_specs": ["-a * dy", "a * dx", "c"],
        "operator_parameters": {"a": 0.8, "c": 0.55},
        "normalize": True,
    },
]


def _plot_indices(selection, n_cases):
    if selection == "all":
        return range(1, n_cases + 1)

    index = int(selection)
    if index < 1 or index > n_cases:
        raise ValueError(f"PLOT_SELECTION must be in [1, {n_cases}] or 'all'.")

    return [index]


write_spherical_replication_vectorfield_sweep(
    output_dir=OUTPUT_DIR,
    R=20.0,
    a=1.0,
    atomtype="Fe",
    n_replications=256,
    field_parameter_cases=FIELD_PARAMETER_CASES,
    n_rotdata=1,
    beta_min=0.0,
    beta_max=np.pi,
    rotation_seed=1,
)

print("Generated vector-field variants:")
for index, field_params in enumerate(FIELD_PARAMETER_CASES, start=1):
    print(f"  m_{index}.csv: {field_params['label']}")

magdata_dir = OUTPUT_DIR / "RealSpaceData" / "MagData" / "Object_1"

IMAGE_DIR.mkdir(parents=True, exist_ok=True)

for index in _plot_indices(PLOT_SELECTION, len(FIELD_PARAMETER_CASES)):
    field_params = FIELD_PARAMETER_CASES[index - 1]
    print(f"Plotting m_{index}.csv: {field_params['label']}")

    if field_params.get("field_type") in ("transversal_helix", "longitudinal_helix"):
        cb = "mx"
    else:
        cb = "mz"

    mag_file = magdata_dir / f"m_{index}.csv"

    if INTERACTIVE:
        plot_magnetization_file(
            mag_file,
            vector_scale=1.5,
            center_vectors=True,
            arrow_tip_length=0.63,
            arrow_tip_radius=1.5 * 0.17,
            arrow_tip_resolution=24,
            arrow_shaft_radius=1.5 * 0.075,
            arrow_shaft_resolution=16,
            downsample_step=1,
            show_points=True,
            color_by=cb,
            scalar_bar_title=cb,
            enable_cut_sliders=True,
        )
    else:
        output_png = IMAGE_DIR / f"m_{index:02d}_{field_params['label'].replace(' ', '_')}.png"
        plot_magnetization_file(
            mag_file,
            vector_scale=1.5,
            center_vectors=True,
            arrow_tip_length=0.63,
            arrow_tip_radius=1.5 * 0.17,
            arrow_tip_resolution=24,
            arrow_shaft_radius=1.5 * 0.075,
            arrow_shaft_resolution=16,
            downsample_step=1,
            show_points=True,
            color_by=cb,
            scalar_bar_title=cb,
            show=False,
            output_png=output_png,
        )
