# -- Project information -----------------------------------------------------
project = "NuMagSANS"
author = "Michael Philipp Adams"
copyright = "2025, Michael Philipp Adams"
release = "beta1"

# -- General configuration ---------------------------------------------------
extensions = [
    "myst_parser",
    "breathe",
    "exhale",
    "sphinx.ext.autodoc",
]

# Markdown und reStructuredText
source_suffix = [".rst", ".md"]
myst_enable_extensions = ["colon_fence"]

# Hauptseite (für GitHub Pages wichtig)
master_doc = "index"

# Doxygen XML → Breathe
breathe_projects = {"NuMagSANS": "./doxygen/xml"}
breathe_default_project = "NuMagSANS"

# -- Options for HTML output -------------------------------------------------
html_theme = "furo"

# Falls du später Custom CSS oder Logos einfügen möchtest:
# html_static_path = ["_static"]

html_theme_options = {
    "sidebar_hide_name": False,
}

# -- Autodoc: Python module path ---------------------------------------------
import os
import sys
sys.path.insert(0, os.path.abspath("../src"))
