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

# Markdown + reStructuredText
source_suffix = [".rst", ".md"]
myst_enable_extensions = ["colon_fence"]

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

# Exhale-Konfiguration (vollständig, getestet)
exhale_args = {
    # Der Ordner in docs/, in dem Exhale seine generierten .rst-Dateien ablegt
    "containmentFolder": "./api",

    # Root-Datei für die API-Referenz
    "rootFileName": "library_root.rst",
    "rootFileTitle": "C++ API Reference",

    # Entfernt den oberen Teil des Pfades in Doxygen-Daten
    "doxygenStripFromPath": "..",

    # GUI-Tree-View in der Seitenleiste
    "createTreeView": True,

    # Doxygen wird *nicht* von Exhale selbst ausgeführt
    "exhaleExecutesDoxygen": False,

    # Optional: weniger Log-Ausgabe
    "verboseBuild": False,
}

# -- Options for HTML output -------------------------------------------------
html_theme = "furo"

# -- Autodoc: Python modules from src/ ---------------------------------------
import sys
sys.path.insert(0, os.path.abspath("../src"))

