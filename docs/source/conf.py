"""Sphinx configuration for the itacart documentation.

The build reads the **working tree**, never an installed copy — see the
``sys.path`` bootstrap below. On a development branch the version pip has
installed is almost never the version under test, and documentation generated
from the wrong source is worse than no documentation at all.

The API reference is not written by hand. ``sphinx.ext.autosummary`` walks the
package recursively and generates one page per module and per public name out
of the docstrings themselves, so a module delivered by a later phase appears in
the documentation without editing anything here.
"""

from __future__ import annotations

import sys
from pathlib import Path

# --------------------------------------------------------------------------
# Import the package from the working tree
# --------------------------------------------------------------------------

SOURCE_DIR = Path(__file__).resolve().parent
ROOT = SOURCE_DIR.parents[1]

sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(SOURCE_DIR / "_ext"))

import itacart  # noqa: E402  (must follow the sys.path bootstrap)

_loaded_from = Path(itacart.__file__).resolve().parent
if _loaded_from != (ROOT / "src" / "itacart").resolve():
    raise RuntimeError(
        f"itacart was imported from {_loaded_from}, not from the working tree "
        f"at {ROOT / 'src' / 'itacart'}. An installed copy is shadowing it; "
        "uninstall it or build in a clean environment."
    )

# --------------------------------------------------------------------------
# Project metadata, read from the package rather than restated
# --------------------------------------------------------------------------

project = "ITACaRT"
copyright = "2025-2026, Instituto Tecnologico de Aeronautica"
author = "Israel Nunes da Silva"
release = itacart.__version__
version = release

# --------------------------------------------------------------------------
# Extensions
# --------------------------------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx.ext.mathjax",
    "sphinx.ext.coverage",
    "sphinx_autodoc_typehints",
    "sphinx_copybutton",
    "myst_nb",  # supersedes myst_parser: same MyST, plus execution
    "phase_figures",  # local, in _ext/
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
# "myst-nb" is the parser name myst_nb registers for Markdown; the older
# "markdown" belonged to myst_parser and no longer exists here. Naming it
# explicitly keeps the mapping visible instead of relying on the extension
# to register the suffix behind our back.
source_suffix = {".rst": "restructuredtext", ".md": "myst-nb"}

# --------------------------------------------------------------------------
# autosummary and autodoc — what makes the reference self-generating
# --------------------------------------------------------------------------

autosummary_generate = True

# The package facade re-exports 109 names from the submodules. Documenting them
# at package level as well would duplicate every entry and produce a page of
# cross-references to itself, so the facade shows only what it defines and the
# submodule pages carry the rest.
autosummary_imported_members = False

# ``itacart.constants`` exports both the REFINEMENT_ALPHABET table and the
# refinement_alphabet() accessor. autosummary derives a stub filename from the
# object name, so on a case-insensitive filesystem -- NTFS, APFS -- the two
# stubs are the same file: the second write overwrites the first and the build
# then warns that a stub is missing. Linux never sees it, and neither does the
# Ubuntu runner in CI, so a strict build fails only on a contributor's own
# machine. Mapping one of the two to a distinct filename removes the collision
# on every filesystem.
autosummary_filename_map = {
    "itacart.constants.refinement_alphabet": "itacart.constants.refinement_alphabet_fn",
}

# Deliberately no "members" here. A module page is an index: its own docstring
# plus the summary tables the template builds. The detail lives on one page per
# name, generated into api/generated/. Turning members on as well documents
# every object twice — once inline and once on its own page — and Sphinx says
# so, loudly, as "duplicate object description".
autodoc_default_options = {
    "show-inheritance": True,
    "member-order": "bysource",
}
autodoc_typehints = "description"
autodoc_preserve_defaults = True
# Only the optional extras are mocked. shapely is a hard runtime dependency,
# so mocking it would render every shapely type in a signature as a mock.
autodoc_mock_imports = ["geopandas", "joblib"]

# Google style throughout, as the repository already had it. The one file that
# used NumPy style, itacart_core/compositional_index.py, is an origin document
# and not part of the package. Flip numpy_docstring back on if a module ported
# from itacart_core ever lands with "Parameters / ----------" headings.
napoleon_google_docstring = True
napoleon_numpy_docstring = False
napoleon_use_admonition_for_notes = True
napoleon_use_rtype = False

always_document_param_types = True
typehints_defaults = "comma"

coverage_show_missing_items = True

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "shapely": ("https://shapely.readthedocs.io/en/stable", None),
    "geopandas": ("https://geopandas.org/en/stable", None),
}

# --------------------------------------------------------------------------
# MyST
# --------------------------------------------------------------------------

myst_enable_extensions = ["colon_fence", "deflist", "dollarmath"]
myst_heading_anchors = 3

# The worked examples are MyST text, not .ipynb: the repository ignores
# notebooks at any depth and that is deliberate. They are executed at build
# time instead, so a documented example that stops working stops the build
# rather than quietly lying to the reader. Nothing is cached between runs —
# the freshly executed output is the only output the reader ever sees.
nb_execution_mode = "auto"
nb_execution_raise_on_error = True
nb_execution_timeout = 120

# --------------------------------------------------------------------------
# HTML output
# --------------------------------------------------------------------------

html_theme = "sphinx_rtd_theme"
html_title = f"itacart {release}"
# Charts and images are maximised on purpose: the figures of the paper live
# in docs/_static and are referenced from the concept pages.
#
# Two entries, on purpose. "../_static" is the repository's own directory and
# is kept so that any page already pointing at _static/f1/... keeps resolving.
# "_static" is source-local, for assets that belong to the documentation
# rather than to a phase. Drop the first if no page references it directly:
# the phase galleries do not need it, they copy what they use.
html_static_path = ["_static", "../_static"]
html_copy_source = False
html_show_sphinx = False

# The documentation is published to GitHub Pages at a fixed address, so the
# canonical URL is knowable at build time and Sphinx can emit a correct
# <link rel="canonical"> on every page. The theme carries an older option of
# its own, ``canonical_url``, for the same purpose; it is deliberately left
# unset so that there is one source of truth rather than two.
html_baseurl = "https://itacart.github.io/itacart/"

html_theme_options = {
    # Concepts & Guides nests three levels deep, and the theme collapses every
    # branch but the current one by default, which would hide most of the tree
    # from the sidebar. Three is the depth the content needs. A fourth level
    # adds nothing a reader navigates by — it is the per-name pages generated
    # under the API reference — and the sidebar is repeated verbatim on every
    # page, so the depth is paid once per page rather than once per site.
    "collapse_navigation": False,
    "navigation_depth": 3,
    # A ``:hidden:`` toctree is hidden from the page body, not from this
    # theme's sidebar, which is built from the global toctree and shows
    # hidden entries unless told otherwise. This is what actually keeps
    # the phase galleries out of the navigation while leaving them built
    # and linkable from the pages that use their figures.
    "includehidden": False,
    "sticky_navigation": True,
    # The one piece of visual identity carried over from the previous theme.
    # This theme has no dark variant, so the darker of the two brand colours
    # is the one that survives.
    "style_nav_header_background": "#2f6f9f",
}

# No "Edit on GitHub" link for now. The repository that serves this site is
# not yet the upstream where these pages live, so any edit link would point at
# a branch that does not contain the page being read. To restore it, set
# ``display_github`` to True and add ``github_user``, ``github_repo``,
# ``github_version`` and ``conf_py_path`` ("/docs/source/").
html_context = {
    "display_github": False,
}

# --------------------------------------------------------------------------
# Phase figure galleries
# --------------------------------------------------------------------------
#
# The verification notebooks write their figures to docs/_static/fN/, which is
# outside this source directory: the notebooks themselves are not versioned but
# the figures they produce are. The phase_figures extension copies
# each directory into _generated/ at build time and writes one gallery page per
# phase, so a new phase's figures need no edit here either.

phase_figures_root = ROOT / "docs" / "_static"
phase_figures_titles = {
    "f0": "F0 — Bootstrap",
    "f1": "F1 — Geodesy",
    "f2": "F2 — Compositional index",
    "f3": "F3 — Resolutions and cells",
    "f4": "F4 — Boundary behaviour",
    "f5": "F5 — Hierarchy",
    "f6": "F6 — Topology",
    "f7": "F7 — Vector geometry",
    "f8": "F8 — Serialization",
    "f9": "F9 — Interoperability and conformance",
    # Not a phase: the eight figures published with the paper, kept at the
    # resolution they were authored at. The key starts with "f" because the
    # extension only scans directories that do, and sorts after "f9".
    "fpaper": "Figures of the paper",
}

# The paper's figures were authored by hand for the publication, not
# produced by a verification notebook, and the gallery must not claim
# otherwise. Every other gallery keeps the default.
phase_figures_provenance = {
    "fpaper": "as published with the paper",
}
