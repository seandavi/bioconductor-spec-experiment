"""Shared type aliases."""

from __future__ import annotations

from typing import Any

import numpy as np
import numpy.typing as npt

IntArray = npt.NDArray[np.int64]
BoolArray = npt.NDArray[np.bool_]
# Anything accepted where a per-range integer argument is wanted (CONV.6).
IntLike = int | np.integer[Any] | npt.ArrayLike
# Integer result where NA positions are masked (CONV.4).
MaskedInts = np.ma.MaskedArray[Any, np.dtype[np.int64]]
