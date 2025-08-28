# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'ExaGOOP'
copyright = '2025, Sreejith N.A., Hariswaran Sitaraman, Nick Deak, Marc Day'
author = 'Sreejith N.A., Hariswaran Sitaraman, Nick Deak, Marc Day'
release = '1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['breathe']
breathe_projects = {"ExaGOOP","_build/xml"}
breathe_default_project = "ExaGOOP"
html_theme = "furo"

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']
