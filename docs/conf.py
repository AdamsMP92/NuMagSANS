# -- Project information -----------------------------------------------------
project = "NuMagSANS"
author = "Michael Philipp Adams"
copyright = "2025, Michael Philipp Adams"
release = "beta1"

# -- General configuration ---------------------------------------------------
extensions = [
    "breathe",
    "myst_parser",
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
html_static_path = ["_static"]

html_theme_options = {
    "sidebar_hide_name": False,
    "light_logo": "logo.png",  # optional, falls du ein Logo hast
    "dark_logo": "logo.png",   # dito
}

# Optional: Falls du später Python oder CUDA-Code referenzierst
import os, sys
sys.path.insert(0, os.path.abspath('../src'))
