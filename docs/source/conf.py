import os
import sys
sys.path.insert(0, os.path.abspath('../../analysis'))  # Point this to your code directory if needed
sys.path.insert(0, os.path.abspath('../../run'))      # Point this to your code directory if needed
sys.path.insert(0, os.path.abspath('../../include'))  # Point this to your code directory if needed
sys.path.insert(0, os.path.abspath('../../src'))      # Point this to your code directory if needed
sys.path.insert(0, os.path.abspath('../'))       # Points to the docs directory
sys.path.insert(0, os.path.abspath('.'))  # Points to the source directory where your .rst files are


# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'G4FastSim'
copyright = '2024, Auke Colijn'
author = 'Auke Colijn'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

html_extra_path = ['../build/html']
import subprocess
subprocess.call('cd ../.. ; doxygen docs/Doxyfile', shell=True)


extensions = [
    'sphinx.ext.autodoc',   # Auto-generates documentation from docstrings
    'sphinx.ext.napoleon',  # Supports Google-style and NumPy-style docstrings
    'sphinx.ext.viewcode',  # Links to the source code in the documentation
    'breathe'              # Integrates Doxygen with Sphinx
]
breathe_projects = {
    "G4FastSim": "../doxygen_output/xml"
}
breathe_default_project = "G4FastSim"

templates_path = ['_templates']
exclude_patterns = []

autodoc_mock_imports = ['pandas', 'shutil', 'json']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
