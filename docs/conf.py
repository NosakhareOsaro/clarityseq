"""
Sphinx configuration for GenomeForge documentation.

Builds autodoc from Python docstrings (Google style via napoleon),
Mermaid diagrams via sphinxcontrib-mermaid, and furo theme.
"""

project = "GenomeForge"
copyright = "2026, GenomeForge Contributors"
author = "GenomeForge Contributors"
release = "0.1.0"

extensions = [
    "sphinx.ext.autodoc",           # Auto-document from Python docstrings
    "sphinx.ext.napoleon",          # Google-style docstrings
    "sphinx.ext.viewcode",          # Source code links
    "sphinx.ext.intersphinx",       # Cross-project links
    "sphinx.ext.autosummary",       # Module summary tables
    "myst_parser",                  # Markdown support
    "sphinxcontrib.mermaid",        # Mermaid diagrams
]

# Napoleon config for Google-style docstrings
napoleon_google_docstring = True
napoleon_numpy_docstring = False
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_use_admonition_for_examples = True

# Furo theme (modern, accessible)
html_theme = "furo"
html_theme_options = {
    "sidebar_hide_name": False,
    "light_css_variables": {
        "color-brand-primary": "#2563eb",
        "color-brand-content": "#1d4ed8",
    },
}

# Autodoc settings
autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
    "member-order": "bysource",
}

# Source suffixes
source_suffix = {".rst": "restructuredtext", ".md": "markdown"}

# Intersphinx mapping to external docs
intersphinx_mapping = {
    "python": ("https://docs.python.org/3.12/", None),
    "pydantic": ("https://docs.pydantic.dev/latest/", None),
}
