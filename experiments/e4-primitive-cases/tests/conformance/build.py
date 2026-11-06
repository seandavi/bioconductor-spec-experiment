"""Build Python objects from primitive input specs (mirror of runners/r/build.R)."""

from __future__ import annotations

from typing import Any

import numpy as np

from ranges_e4 import IRanges


def _vec(x: Any, kind: str | None = None) -> Any:
    if x is None:
        return None
    if not isinstance(x, list):
        x = [x]
    if kind == "int":
        return [None if v is None else int(v) for v in x]
    if kind == "num":
        return [None if v is None else float(v) for v in x]
    return x


def _ranges_args(s: dict[str, Any]) -> dict[str, Any]:
    return {k: _vec(s.get(k), "int") for k in ("start", "end", "width")} | {
        "names": _vec(s.get("names")),
        "mcols": s.get("mcols"),
    }


def build_input(spec: dict[str, Any], inputs: dict[str, Any]) -> Any:
    (kind, s), = spec.items()
    match kind:
        case "IRanges":
            return IRanges(**_ranges_args(s))
        case "Seqinfo":
            return _seqinfo(s)
        case "GRanges":
            from ranges_e4 import GRanges

            si = s.get("seqinfo")
            if isinstance(si, str):
                si = inputs[si]
            elif si is not None:
                si = _seqinfo(si)
            return GRanges(
                seqnames=_vec(s["seqnames"]),
                strand=_vec(s.get("strand")),
                seqinfo=si,
                **_ranges_args(s),
            )
        case "int":
            return np.array(_vec(s, "int"), dtype=object if None in _vec(s, "int") else np.int64)
        case "num":
            return np.array(_vec(s, "num"), dtype=float)
        case "chr":
            return np.array(_vec(s), dtype=object)
        case "lgl":
            return np.array(_vec(s), dtype=bool)
    raise ValueError(f"unknown input kind: {kind}")


def _seqinfo(s: dict[str, Any]) -> Any:
    from ranges_e4 import Seqinfo

    return Seqinfo(
        seqnames=_vec(s["seqnames"]),
        seqlengths=_vec(s.get("seqlengths"), "int"),
        is_circular=_vec(s.get("isCircular")),
        genome=_vec(s.get("genome")),
    )


def build_inputs(specs: dict[str, Any]) -> dict[str, Any]:
    inputs: dict[str, Any] = {}
    for name, spec in specs.items():
        inputs[name] = build_input(spec, inputs)
    return inputs
