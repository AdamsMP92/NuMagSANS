"""Workflow helpers for spherical nanoparticle systems.

Workflows orchestrate the lower SystemDesigner layers into one callable unit.
They keep frequently used system-building recipes compact while still exposing
the physically relevant parameters at the function boundary.
"""

from ..AssemblyBaseDesigner.AssemblyBaseAnalyzer import (
    analyze_assembly_base_dataset,
    write_parameter_table,
)
from ..AssemblyBaseDesigner.AssemblyBaseTemplates import (
    write_monodisperse_spherical_nanoparticle_base,
)
from ..AssemblyBaseMagnetizer.AssemblyBaseMagnetizer import (
    write_spherical_magnetization_data,
)
from ..AssemblyBaseMagnetizer.MagnetizationBaseTemplates import (
    spherical_vortex_magnetization_base,
)


def write_monodisperse_spherical_vortex_system(
    output_dir,
    R=10.0,
    a=1.0,
    atomtype="Fe",
    n_objects=500,
    kappa=1.0,
    profile_type="linear",
    xi_type="cylindrical_xi",
    turns=1.0,
    bins=20,
    write_parameter_csv=True,
    plot_parameter_histograms=True,
    output_file="m_1.csv",
):
    """Generate a monodisperse spherical-particle system with local vortex data.

    The workflow creates local spherical crystal objects, analyzes the realized
    object metadata, optionally writes/plots parameter distributions, and then
    materializes local magnetization data on the generated objects.

    Parameters
    ----------
    output_dir : str or pathlib.Path
        Output directory containing the generated ``RealSpaceData`` tree.
    R : float, optional
        Sphere radius used by the crystal template.
    a : float, optional
        Simple-cubic lattice constant.
    atomtype : str, optional
        Per-atom label written by the crystal template.
    n_objects : int, optional
        Number of local objects to materialize.
    kappa, profile_type, xi_type, turns : optional
        Parameters passed to the spherical-vortex magnetization template.
    bins : int, optional
        Histogram bin count for parameter-distribution plots.
    write_parameter_csv : bool, optional
        If ``True``, write ``RealSpaceData/parameter_table.csv``.
    plot_parameter_histograms : bool, optional
        If ``True``, write ``RealSpaceData/parameter_distributions.png``.
    output_file : str, optional
        Magnetization filename written inside each ``Object_i`` directory.

    Returns
    -------
    dict
        Summary containing object-generation, analysis, optional output paths,
        and magnetization-generation results.
    """
    output_dir = str(output_dir)

    object_summary = write_monodisperse_spherical_nanoparticle_base(
        R=R,
        a=a,
        atomtype=atomtype,
        n_objects=n_objects,
        output_dir=output_dir,
        name="monodisperse spherical nanoparticle system",
    )

    analysis = analyze_assembly_base_dataset(output_dir)

    parameter_table_path = None
    if write_parameter_csv:
        parameter_table_path = write_parameter_table(
            analysis,
            f"{output_dir}/RealSpaceData/parameter_table.csv",
        )

    plot_path = None
    if plot_parameter_histograms:
        from ..AssemblyBaseDesigner.AssemblyBasePlot import plot_parameter_distributions

        plot_path = f"{output_dir}/RealSpaceData/parameter_distributions.png"
        plot_parameter_distributions(
            analysis,
            output_path=plot_path,
            bins=bins,
        )

    mag_summary = write_spherical_magnetization_data(
        assembly_dir=output_dir,
        magnetization_base_data=spherical_vortex_magnetization_base(
            kappa=kappa,
            profile_type=profile_type,
            xi_type=xi_type,
            turns=turns,
        ),
        output_file=output_file,
    )

    return {
        "output_dir": output_dir,
        "object_summary": object_summary,
        "analysis": analysis,
        "parameter_table_path": parameter_table_path,
        "parameter_plot_path": plot_path,
        "magnetization_summary": mag_summary,
    }


if __name__ == "__main__":
    write_monodisperse_spherical_vortex_system("AssemblyExample1")
