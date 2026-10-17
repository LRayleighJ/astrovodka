"""Re-export MARTINI's IllustrisTNG snapshot I/O helpers.

The implementation lives in :mod:`martini.sources._illustris_tools`; astrovodka
imports it rather than carrying a copy. This module exists so that the TNG
data-loading code in astrovodka — and the tutorial, which imports these names
directly — keeps one stable import path.
"""

from martini.sources._illustris_tools import (
    gcPath,
    getNumPart,
    getSnapOffsets,
    loadHeader,
    loadSingle,
    loadSubset,
    offsetPath,
    partTypeNum,
    snapPath,
)

__all__ = [
    "gcPath",
    "getNumPart",
    "getSnapOffsets",
    "loadHeader",
    "loadSingle",
    "loadSubset",
    "offsetPath",
    "partTypeNum",
    "snapPath",
]
