"""
plotmol
Interactive plotting of data annotated with molecule structures.
"""

from ._version import get_versions
from .plotmol import scatter

__all__ = [scatter]

__version__ = get_versions()["version"]
del get_versions
