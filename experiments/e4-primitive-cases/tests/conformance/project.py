"""Project Python results onto the golden's primitive fields (mirror of project.R).

Dispatch is by class name so the projector works at every build stage.
Index vectors are 0-based in Python; projection adds 1 (ADR-0001).
"""

from __future__ import annotations

from typing import Any

import numpy as np
import polars as pl


def _vec(x: Any, plus_one: bool = False) -> list[Any]:
    """numpy / masked array -> list; masked -> None; optional +1 on non-null."""
    if isinstance(x, np.ma.MaskedArray):
        mask = np.ma.getmaskarray(x).tolist()
        vals = x.data.tolist()
        return [None if m else (v + 1 if plus_one else v) for v, m in zip(vals, mask, strict=True)]
    vals = np.asarray(x).tolist() if not isinstance(x, list) else x
    return [None if v is None else (v + 1 if plus_one else v) for v in vals]


def _mcols(df: pl.DataFrame) -> dict[str, list[Any]]:
    out: dict[str, list[Any]] = {}
    for name in df.columns:
        col = df[name]
        # `revmap` holds positions; every other column is data
        if col.dtype == pl.List:
            out[name] = [None if r is None else [v + 1 if name == "revmap" else v for v in r] for r in col.to_list()]
        else:
            out[name] = col.to_list()
    return out


def _seqinfo(si: Any) -> dict[str, Any]:
    return {
        "seqnames": list(si.seqnames),
        "seqlengths": list(si.seqlengths.values()),
        "isCircular": list(si.is_circular.values()),
        "genome": list(si.genome.values()),
    }


def project(x: Any, is_index: bool = False, kind_hint: str | None = None) -> dict[str, Any]:
    kind = type(x).__name__
    if kind == "GRanges":
        return {
            "kind": "GRanges",
            "seqnames": list(x.seqnames), "start": _vec(x.start), "end": _vec(x.end),
            "width": _vec(x.width), "strand": list(x.strand),
            "names": None if x.names is None else list(x.names),
            "mcols": _mcols(x.mcols), "seqinfo": _seqinfo(x.seqinfo),
        }
    if kind == "IRanges":
        return {
            "kind": "IRanges",
            "start": _vec(x.start), "end": _vec(x.end), "width": _vec(x.width),
            "names": None if x.names is None else list(x.names),
            "mcols": _mcols(x.mcols),
        }
    if kind == "Hits":
        return {
            "kind": "Hits", "query": _vec(x.query, True), "subject": _vec(x.subject, True),
            "query_length": x.query_length, "subject_length": x.subject_length,
        }
    if kind == "Rle":
        return {"kind": "Rle", "values": _vec(x.values), "lengths": _vec(x.lengths)}
    if kind == "Seqinfo":
        return {"kind": "Seqinfo"} | _seqinfo(x)
    # an empty dict is ambiguous (RleList vs named vector): the op says which
    is_rle_list = isinstance(x, dict) and all(type(v).__name__ == "Rle" for v in x.values()) and (x or kind_hint == "coverage")
    if is_rle_list:
        return {
            "kind": "RleList", "names": list(x),
            "elements": [{"values": _vec(r.values), "lengths": _vec(r.lengths)} for r in x.values()],
        }
    if isinstance(x, dict):  # named atomic vector (seqlengths, isCircular, genome)
        return {"kind": "vector", "values": _vec(list(x.values())), "names": list(x)}
    if isinstance(x, (int, np.integer)) and not isinstance(x, bool):
        return {"kind": "vector", "values": [int(x) + (1 if is_index else 0)], "names": None}
    if isinstance(x, (np.ndarray, list)):
        return {"kind": "vector", "values": _vec(x, is_index), "names": None}
    raise TypeError(f"cannot project {type(x).__name__}")
