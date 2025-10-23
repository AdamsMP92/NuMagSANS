# -- Project information -----------------------------------------------------
project = "NuMagSANS"
copyright = "©2025, Michael Philipp Adams"
author = "Michael Philipp Adams"
release = "beta1"

# -- General configuration ---------------------------------------------------
extensions = [
    "breathe",
    "myst_parser",       # erlaubt Markdown in der Dokumentation
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.todo",
]

breathe_projects = {"NuMagSANS": "../docs/doxygen/xml"}
breathe_default_project = "NuMagSANS"

templates_path = ["_templates"]
exclude_patterns = []
language = "en"

# -- Options for HTML output -------------------------------------------------
html_theme = "furo"   

html_static_path = ["_static"]

# -- Furo theme options ------------------------------------------------------
html_theme_options = {
    "sidebar_hide_name": False,
    "light_logo": "logo_light.png",
    "dark_logo": "logo_dark.png",
    "light_css_variables": {
        "color-brand-primary": "#4B8BBE",
        "color-brand-content": "#306998",
    },
    "dark_css_variables": {
        "color-brand-primary": "#FFE873",
        "color-brand-content": "#FFD43B",
    },
}

# -- Optional custom CSS (falls du Feinheiten ändern willst) -----------------
html_css_files = ["custom.css"]
