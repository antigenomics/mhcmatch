"""Sphinx configuration for mhcmatch."""
import datetime
import os
import sys

# Import the package from src/ without installing it (seqtree is mocked below).
sys.path.insert(0, os.path.abspath("../src"))

project = "mhcmatch"
author = "ISALGO laboratory"
copyright = f"{datetime.date.today().year}, {author}"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.githubpages",
]

# Heavy / native dependencies mocked at doc-build time (seqtree ships a C++ core; the rest are
# optional). Docs build with only sphinx + the theme installed.
autodoc_mock_imports = ["seqtree", "numpy", "logomaker", "matplotlib", "pandas"]
autosummary_generate = False
autodoc_member_order = "bysource"
autodoc_typehints = "description"
napoleon_google_docstring = True

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "pydata_sphinx_theme"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_title = "mhcmatch"
html_theme_options = {
    "github_url": "https://github.com/antigenomics/mhcmatch",
    "show_prev_next": False,
}
