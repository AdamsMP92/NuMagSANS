from .CrystalPlot import scatter_plot_crystal
from .CrystalTemplates import sc_cylinder_crystal

crystal = sc_cylinder_crystal(a=0.8, R=10.0, h=8.0, atomtype="Fe")

scatter_plot_crystal(crystal, "Simple Cubic Fe Lattice with Cylindrical Cut")
