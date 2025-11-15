# -- Project information -----------------------------------------------------
project = "NuMagSANS"
author = "Michael Philipp Adams"
copyright = "2025"

# -- General configuration ---------------------------------------------------
extensions = [
    "myst_parser",
    "sphinx_copybutton",
    "sphinx_design",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "breathe",
    "exhale",
    "sphinx.ext.mathjax",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

myst_enable_extensions = [
    "colon_fence",
    "dollarmath",
]

# -- Breathe / Exhale configuration ------------------------------------------
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

# -- HTML output -------------------------------------------------------------
html_theme = "pydata_sphinx_theme"
html_static_path = ["_static"]
html_css_files = ["css/custom.css"]   # ← WICHTIG!

#html_logo = "_static/logo.png"
#html_favicon = "_static/favicon.png"

html_show_sourcelink = False
html_sourcelink_suffix = ""

html_theme_options = {
    "external_links": [
        {"name": "NuMagSANS GitHub", "url": "https://github.com/AdamsMP92/NuMagSANS"},
    ],
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/AdamsMP92",
            "icon": "fab fa-github-square",
        },
    ],
    "header_links_before_dropdown": 6,
    "logo": {
        "text": "Documentation",
    },
    "secondary_sidebar_items": ["page-toc"],
}

html_sidebars = {
    "**": ["sidebar-nav-bs", "sidebar-ethical-ads"],
}

# -- Python path -------------------------------------------------------------
import os
import sys
sys.path.insert(0, os.path.abspath("../src"))
