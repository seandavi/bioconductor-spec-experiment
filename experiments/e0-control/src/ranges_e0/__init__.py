"""ranges_e0: a plain Python port of Bioconductor IRanges/GRanges (experiment E0)."""

from .hits import Hits
from .iranges import NO_MATCH, IRanges, concat
from .rle import Rle

__all__ = ["Hits", "IRanges", "NO_MATCH", "Rle", "concat"]
