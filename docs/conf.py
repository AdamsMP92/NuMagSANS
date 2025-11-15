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
    "sphinx.ext.mathjax",   # ← wichtig für Rendern der Formeln
]

# Markdown + reStructuredText
source_suffix = [".rst", ".md"]

# MyST-Erweiterungen
myst_enable_extensions = [
    "colon_fence",   # deine bisherige Erweiterung
    "dollarmath",    # ← NEU: aktiviert $$…$$, $…$, \[…\], \(...\)
]

# Hauptdokument (für GitHub Pages wichtig)
master_doc = "index"

# ---------------------------------------------------------------------------
# Breathe + Exhale: Doxygen XML Integration
# ---------------------------------------------------------------------------

# Verweis auf den Ordner, in dem Doxygen das XML erzeugt
breathe_projects = {"NuMagSANS": "./doxygen/xml"}
breathe_default_project = "NuMagSANS"

# Exhale benötigt zusätzliche Imports
import os
import textwrap

# Exhale-Konfiguration
exhale_args = {
    "containmentFolder": "./api",
    "rootFileName": "library_root.rst",
    "rootFileTitle": "C++ API Reference",
    "doxygenStripFromPath": "..",
    "createTreeView": True,
    "exhaleExecutesDoxygen": False,
    "verboseBuild": False,
}

# -- Options for HTML output -------------------------------------------------
html_theme = "furo"

# -- Autodoc: Python modules from src/ ---------------------------------------
import sys
sys.path.insert(0, os.path.abspath("../src"))

