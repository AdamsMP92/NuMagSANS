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

# Damit Markdown-Dateien (index.md) von Sphinx richtig verarbeitet werden
source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}

# Hauptseite explizit setzen, damit GitHub Pages nicht ins Leere zeigt
master_doc = 'index'

# Pfad zu den Doxygen-XML-Dateien (relativ zu docs/)
breathe_projects = {"NuMagSANS": "./doxygen/xml"}
breathe_default_project = "NuMagSANS"

# Quelle kann .md oder .rst sein
source_suffix = [".rst", ".md"]
myst_enable_extensions = ["colon_fence"]

# -- Options for HTML output -------------------------------------------------
html_theme = "furo"

html_theme_options = {
    "light_logo": "logo_light.png",
    "dark_logo": "logo_dark.png",
    "sidebar_hide_name": False,
    "light_css_variables": {
        "color-brand-primary": "#1e90ff",
        "color-brand-content": "#0040a0",
        "color-background-secondary": "#f5f5f8",
    },
    "dark_css_variables": {
        "color-brand-primary": "#84caff",
        "color-brand-content": "#c5e1ff",
        "color-background-secondary": "#0e1117",
    },
}

# Optional: Falls du später Python oder CUDA-Code referenzierst
import os, sys
sys.path.insert(0, os.path.abspath('../src'))
