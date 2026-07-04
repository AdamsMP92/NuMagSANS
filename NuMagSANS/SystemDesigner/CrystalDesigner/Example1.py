from .CrystalPlot import scatter_plot_crystal
from .CrystalTemplates import sc_sphere_crystal

crystal = sc_sphere_crystal(a=0.8, R=10.0, atomtype="Fe")

scatter_plot_crystal(crystal, "Simple Cubic Fe Lattice with Spherical Cut")
