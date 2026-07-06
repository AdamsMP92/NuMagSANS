"""Generate and inspect one helical spherical particle with SystemDesigner."""

import math
from pathlib import Path

from NuMagSANS.SystemDesigner.AssemblyAnalyzer import plot_magnetization_file
from NuMagSANS.SystemDesigner.Workflows import write_spherical_replication_vectorfield_sweep


def main():
    """Create a local helical sphere and open a diagnostic vector-field plot."""
    output_dir = Path("HelicalSphere_R20")

    radius_nm = 20.0
    lattice_constant_nm = 1.0
    pitch_nm = 20.0
    chirality = 1.0
    phase = 0.0

    summary = write_spherical_replication_vectorfield_sweep(
        output_dir=output_dir,
        R=radius_nm,
        a=lattice_constant_nm,
        atomtype="Fe",
        n_replications=1,
        field_parameter_cases={
            "field_type": "transversal_helix",
            "k": 2.0 * math.pi / pitch_nm,
            "chirality": chirality,
            "phase": phase,
        },
        n_rotdata=1,
        beta_min=0.0,
        beta_max=0.0,
        rotation_seed=1,
    )

    mag_file = Path(summary["real_space_dir"]) / "MagData" / "Object_1" / "m_1.csv"

    plot_magnetization_file(
        mag_file,
        vector_scale=1.6,
        center_vectors=True,
        seed=1,
        show_points=False,
        color_by="mx",
        cmap="coolwarm",
        enable_cut_sliders=True,
        arrow_tip_length=0.63,
        arrow_tip_radius=0.15,
        arrow_shaft_radius=0.06,
        #output_png="helical_sphere.png",
    )


if __name__ == "__main__":
    main()
