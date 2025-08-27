# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'NuMagSANS'
copyright = '2025, Michael Philipp Adams'
author = 'Michael Philipp Adams'
release = 'beta1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["breathe", "myst_parser"]
breathe_projects = {"NuMagSANS": "../docs/doxygen/xml"}
breathe_default_project = "NuMagSANS"
html_theme = "furo"

templates_path = ['_templates']
exclude_patterns = []

language = '[en]'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']
