from .AssemblyBaseAnalyzer import (
    analyze_assembly_base_dataset,
    write_parameter_table,
)
from .AssemblyBasePlot import (
    plot_parameter_distributions,
)
from .AssemblyBaseTemplates import (
    write_gaussian_spherical_nanoparticle_base,
)

# Example 2:
# Generate a dilute object set made from spherical Fe nanoparticle templates.
#
# This example only creates the local object geometries:
#
#     RealSpaceData/Local_Objects/Object_1/pos.csv
#     RealSpaceData/Local_Objects/Object_1/meta.csv
#     ...
#
# It does not yet define object-center positions or rotations. That later
# organization step belongs to AssemblyBaseOrganizer and can write StructData.csv.
# After generation, the realized template-parameter distributions are extracted
# from the Object_i/meta.csv files and plotted.


OUTPUT_DIR = "AssemblyExample2"

# Generate a dilute ensemble of spherical nanoparticles with a Gaussian
# distribution of radii.
summary = write_gaussian_spherical_nanoparticle_base(
    R_mean=10.0,
    R_std=3.0,
    a=1.0,
    atomtype="Fe",
    n_objects=500,
    output_dir=OUTPUT_DIR,
    seed=123,
    name="dilute spherical Fe nanoparticles",
)

# Analyze the generated output.
analysis = analyze_assembly_base_dataset(OUTPUT_DIR)

parameter_table_path = write_parameter_table(
    analysis,
    f"{OUTPUT_DIR}/RealSpaceData/parameter_table.csv",
)

# Plot histograms of analyzer output.
plot_path = f"{OUTPUT_DIR}/RealSpaceData/parameter_distributions.png"
plot_parameter_distributions(
    analysis,
    output_path=plot_path,
    bins=20,
)
