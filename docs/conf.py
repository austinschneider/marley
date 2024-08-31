import sys, os, subprocess

from sphinx.highlighting import lexers
from pygments.lexers.web import PhpLexer

project = u'MARLEY'
copyright = u'2016-2024 Steven Gardiner'
master_doc = 'index'
templates_path = [ '_templates' ]
extensions = [ 'sphinxcontrib.bibtex', 'sphinxcontrib.newsfeed',
  'sphinx.ext.todo' ]
source_suffix = '.rst'
version = '1.2.0'
exclude_patterns = ['_build']

highlight_language = 'none'

# TODO notes

#todo_include_todos = True
todo_include_todos = False

# -- HTML theme settings ------------------------------------------------

html_favicon = 'mar.png'
html_show_sourcelink = False
html_sidebars = {
    '**': ['logo-text.html',
           'globaltoc.html',
           'localtoc.html',
           'searchbox.html']
}

import guzzle_sphinx_theme

extensions.append("guzzle_sphinx_theme")
html_theme_path = guzzle_sphinx_theme.html_theme_path()
html_theme = 'guzzle_sphinx_theme'

# Guzzle theme options (see theme.conf for more information)
html_theme_options = {
    "base_url": "www.marleygen.org",
}

#html_add_permalinks = None
html_static_path = [ '_static' ]
html_extra_path = [ '../CITATION.bib', './CNAME' ]

# Bibliography files used by the sphinxcontrib.bibtex extension
bibtex_bibfiles = [ 'marley_pubs.bib', 'external_pubs.bib' ]
