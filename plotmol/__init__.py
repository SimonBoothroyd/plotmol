"""
plotmol

Interactive plotting of data annotated with molecule structures.
"""

from . import _version
from ._plotmol import InputSizeError, default_tooltip_template, scatter

__version__ = _version.get_versions()["version"]
__all__ = ["InputSizeError", "default_tooltip_template", "scatter", "__version__"]
