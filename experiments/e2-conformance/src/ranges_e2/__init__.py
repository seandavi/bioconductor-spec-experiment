"""IRanges/GRanges port (experiment E2: conformance suite first)."""

from .granges import NA_INDEX, GRanges
from .hits import Hits
from .iranges import IRanges
from .rle import Factor, Rle
from .seqinfo import Seqinfo

__version__ = "0.1.0"

__all__ = ["GRanges", "IRanges", "Seqinfo", "Hits", "Rle", "Factor", "NA_INDEX"]
