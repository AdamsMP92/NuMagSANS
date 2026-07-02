from pathlib import Path

import numpy as np

from NuMagSANS.SystemDesigner.AssemblyAnalyzer import plot_magnetization_file
from NuMagSANS.SystemDesigner.Workflows import write_spherical_replication_vectorfield_sweep

OUTPUT_DIR = Path("Playground_VectorFieldSurvey")

# Set this to an integer, for example 3, to inspect only m_3.csv.
# Set it to "all" to walk through all generated vector-field variants.
PLOT_SELECTION = "all"

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
        "label": "linearized vortex",
        "field_type": "linearized_vortex",
        "m_0": 1.0,
        "m_1": 0.05,
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

for index in _plot_indices(PLOT_SELECTION, len(FIELD_PARAMETER_CASES)):
    field_params = FIELD_PARAMETER_CASES[index - 1]
    print(f"Plotting m_{index}.csv: {field_params['label']}")

    if index == 6: 
        cb = "mx"
    else:
        cb = "mz"

    plot_magnetization_file(
        magdata_dir / f"m_{index}.csv",
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
    )
