"""IRanges/GRanges port to Python (experiment e0b-no-r). See README.md."""

from ._sew import NO_HIT
from .iranges import IRanges
from .overlaps import Hits
from .rle import Rle

__all__: list[str] = ["IRanges", "Hits", "Rle", "NO_HIT"]
