"""Sphinx configuration for mhcmatch."""
import datetime
import os
import sys

# Import the package from src/ without installing it (seqtree is mocked below).
sys.path.insert(0, os.path.abspath("../src"))

project = "mhcmatch"
author = "ISALGO laboratory"
copyright = f"{datetime.date.today().year}, {author}"

# **The release, read from the package, never typed into a page.** Docs used to name the version in
# prose -- ">= 1.7.3" in three places, three different numbers, none of them the shipped pin -- and
# every one of them went stale at a release nobody remembered to grep. `|release|` substitutes this,
# so a page can say which version it documents and cannot be wrong about it. Read from the installed
# metadata, falling back to pyproject for a docs build in a bare checkout.
def _version() -> str:
    try:
        from importlib.metadata import version as _v
        return _v("mhcmatch")
    except Exception:
        import re
        src = open(os.path.join(os.path.dirname(__file__), "..", "pyproject.toml")).read()
        m = re.search(r'^version = "([^"]+)"', src, re.M)
        return m.group(1) if m else "0.0.0"


release = _version()
version = ".".join(release.split(".")[:2])
rst_prolog = f".. |release| replace:: {release}\n"

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
# Render __init__ docstrings alongside the class docstring: AnchorModel documents all 19 parameters in
# __init__, and without "both" Sphinx publishes only the short class docstring and drops them.
autoclass_content = "both"
autodoc_typehints = "description"
napoleon_google_docstring = True

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "pydata_sphinx_theme"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
# The version is in the brand, on every page: which release a page documents is the first thing a
# reader needs and the last thing anyone remembers to write down.
html_title = f"mhcmatch {release}"
html_theme_options = {
    "github_url": "https://github.com/antigenomics/mhcmatch",
    "show_prev_next": False,
    # **The sidebar was a flat wall of thirteen full page titles**, in source order, with no grouping
    # and nothing marking where the reader is -- so it was scenery, not navigation. The toctree in
    # `index.rst` is now four captioned groups with short link titles; these three options are what
    # make the theme render that structure rather than flatten it.
    "show_nav_level": 2,        # open each caption's pages, do not collapse to the caption alone
    "navigation_depth": 3,      # let a page's own sections show under it
    "collapse_navigation": False,
    "header_links_before_dropdown": 4,
    "show_toc_level": 2,        # the right-hand on-page TOC: subsections too, not just top level
}


# **The model tables are generated here, on every build, from the artifacts that ship.** They are
# the one thing these docs used to refuse to print: six pages carried their own copy of the EPIC
# coefficients and all six went stale together the first time the model was refitted, because
# nothing read them. Generating removes the class of error rather than the page -- `docs/models.rst`
# includes what this writes, and `docs/_generated/` is gitignored so there is no committed copy to
# drift. The README's summary block is the one committed copy, pinned by `tests/test_modeldoc.py`.
def setup(app):
    from mhcmatch import _modeldoc

    _modeldoc.write(os.path.abspath(".."))
    return {"parallel_read_safe": True}
