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
    "sphinx.ext.mathjax",
]

# Markdown + reStructuredText
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

# MyST-Erweiterungen
myst_enable_extensions = [
    "colon_fence",
    "dollarmath",
]

# Hauptdokument
master_doc = "index"

# ---------------------------------------------------------------------------
# Breathe + Exhale
# ---------------------------------------------------------------------------
import os
import textwrap

breathe_projects = {"NuMagSANS": "./doxygen/xml"}
breathe_default_project = "NuMagSANS"

exhale_args = {
    "containmentFolder": "./api",
    "rootFileName": "library_root.rst",
    "rootFileTitle": "C++ API Reference",
    "doxygenStripFromPath": "..",
    "createTreeView": True,
    "exhaleExecutesDoxygen": False,
    "verboseBuild": False,
}

# -- Theme settings -----------------------------------------------------------
html_theme = "pydata_sphinx_theme"

html_theme_options = {
    "navbar_start": ["navbar-logo"],
    "navbar_center": ["navbar-nav"],
    "navbar_end": ["theme-switcher", "navbar-icon-links"],

    # Navigation bar links (top of page)
    "navbar_links": [
        {"name": "Home", "url": "index.html"},
        {"name": "Installation", "url": "installation.html"},
        {"name": "Usage", "url": "usage.html"},
        {"name": "API Reference", "url": "api_reference.html"},
    ],

    # Aktiviert linke Sidebar automatisch
    "show_nav_level": 2,
    "navigation_depth": 4,
}

# Optional: dein Logo (falls vorhanden)
#html_logo = "_static/logo.png"

# -- Autodoc --------------------------------------------------------------
import sys
sys.path.insert(0, os.path.abspath("../src"))


