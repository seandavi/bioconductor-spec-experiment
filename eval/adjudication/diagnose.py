"""Which representation-level normalization makes a mismatch disappear (the "decisive difference").

Each normalizer is applied to both the R and the port projection; the comparator's own `diffs`
then decides equality. The first single normalizer (in LADDER order) that yields no diffs names
the cluster; otherwise the first pair; otherwise the first differing path/kind is used.
"""
from __future__ import annotations

import copy
import itertools
import re
import sys
from typing import Any, Callable

from load import EVAL

sys.path.insert(0, str(EVAL / "src"))
from bioc_eval.compare import diffs  # noqa: E402

Proj = Any


def _walk(o: Proj, f: Callable[[dict], None]) -> Proj:
    if isinstance(o, dict):
        f(o)
        for v in o.values():
            _walk(v, f)
    elif isinstance(o, list):
        for v in o:
            _walk(v, f)
    return o


def nan_as_na(o: Proj) -> Proj:
    def f(d: dict) -> None:
        if "values" in d and isinstance(d["values"], list):
            d["values"] = [None if v == "NaN" else v for v in d["values"]]
    return _walk(o, f)


def num_type(o: Proj) -> Proj:
    """int and dbl storage types treated alike (values still compared numerically)."""
    def f(d: dict) -> None:
        if d.get("type") in ("int", "dbl"):
            d["type"] = "num"
    return _walk(o, f)


def lgl_as_int(o: Proj) -> Proj:
    """logical vs integer/double storage (TRUE == 1)."""
    def f(d: dict) -> None:
        if d.get("type") in ("int", "dbl", "lgl"):
            if d["type"] == "lgl" and isinstance(d.get("values"), list):
                d["values"] = [None if v is None else int(v) for v in d["values"]]
            d["type"] = "num"
    return _walk(o, f)


def drop_names(o: Proj) -> Proj:
    def f(d: dict) -> None:
        for k in ("names", "rownames"):
            if k in d and d.get("kind") != "RleList":
                d[k] = None
    return _walk(o, f)


def _seqinfo_sorted(si: dict) -> dict:
    idx = sorted(range(len(si["seqnames"])), key=lambda i: si["seqnames"][i])
    return {**si, **{k: [si[k][i] for i in idx] for k in ("seqnames", "seqlengths", "isCircular", "genome")}}


def seqlevel_order(o: Proj) -> Proj:
    """Seqinfo compared as a set of rows (seqlevel order ignored); RleList elements too."""
    if not isinstance(o, dict):
        return o
    if o.get("kind") == "Seqinfo":
        return _seqinfo_sorted(o)
    if o.get("kind") == "GRanges":
        o["seqinfo"] = _seqinfo_sorted(o["seqinfo"])
    if o.get("kind") == "RleList" and o.get("names"):
        idx = sorted(range(len(o["names"])), key=lambda i: o["names"][i])
        o["names"] = [o["names"][i] for i in idx]
        o["elements"] = [o["elements"][i] for i in idx]
    return o


def drop_seqinfo(o: Proj) -> Proj:
    if isinstance(o, dict) and o.get("kind") == "GRanges":
        o["seqinfo"] = None
    return o


def seqinfo_unknowns(o: Proj) -> Proj:
    """Seqinfo rows kept, but seqlengths/isCircular/genome content ignored."""
    def f(d: dict) -> None:
        if d.get("kind") == "Seqinfo":
            for k in ("seqlengths", "isCircular", "genome"):
                d[k] = None
    return _walk(o, f)


def drop_mcols(o: Proj) -> Proj:
    if isinstance(o, dict):
        if "mcols" in o:
            o["mcols"] = []
        if o.get("kind") == "DataFrame":
            o["columns"] = []
    return o


def mcols_types(o: Proj) -> Proj:
    """mcols column storage types (int/dbl/lgl/chr/Factor) ignored; values compared as strings."""
    def col(c: dict) -> dict:
        if c.get("kind") == "Factor":
            return {"v": [None if k is None else c["levels"][k - 1] for k in c["codes"]]}
        if c.get("kind") == "vector":
            return {"v": [None if v is None else str(float(v)) if isinstance(v, (int, float)) and not isinstance(v, bool)
                          else str(v) for v in c["values"]]}
        return c
    if isinstance(o, dict):
        if "mcols" in o:
            o["mcols"] = [[n, col(c)] for n, c in o["mcols"]]
        if o.get("kind") == "DataFrame":
            o["columns"] = [[n, col(c)] for n, c in o["columns"]]
    return o


def _rows(o: dict, fields: tuple[str, ...]) -> list[tuple]:
    return list(zip(*(map(repr, o[k]) if isinstance(o[k], list) else [repr(o[k])] * len(o["start"]) for k in fields)))


def row_order(o: Proj) -> Proj:
    """Element order ignored (ranges, vectors, Hits pairs)."""
    if not isinstance(o, dict):
        return o
    k = o.get("kind")
    if k in ("IRanges", "GRanges"):
        fields = tuple(f for f in ("seqnames", "start", "width", "strand", "names") if f in o and isinstance(o[f], list))
        cols = [(n, c) for n, c in o["mcols"] if c.get("kind") == "vector"]
        rows = [tuple(o[f][i] for f in fields) + tuple(repr(c["values"][i]) for _, c in cols)
                for i in range(len(o["start"]))]
        order = sorted(range(len(rows)), key=lambda i: repr(rows[i]))
        for f in fields:
            o[f] = [o[f][i] for i in order]
        o["mcols"] = [[n, ({**c, "values": [c["values"][i] for i in order]} if c.get("kind") == "vector" else c)]
                      for n, c in o["mcols"]]
    elif k == "vector":
        pairs = sorted(zip(map(repr, o["values"]), o["values"], o["names"] or [None] * len(o["values"])))
        o["values"] = [p[1] for p in pairs]
        o["names"] = None if o["names"] is None else [p[2] for p in pairs]
    elif k == "Hits":
        pairs = sorted(zip(o["from"], o["to"]))
        o["from"], o["to"] = [a for a, _ in pairs], [b for _, b in pairs]
    elif k == "IntegerList":
        o["elements"] = [sorted(e, key=repr) for e in o["elements"]]
    return o


def mcol_order(o: Proj) -> Proj:
    """mcols column order ignored; revmap/IntegerList element order ignored."""
    if isinstance(o, dict) and "mcols" in o:
        o["mcols"] = sorted(([n, (row_order(copy.deepcopy(c)) if c.get("kind") == "IntegerList" else c)]
                             for n, c in o["mcols"]), key=lambda nc: nc[0])
    return o


def strand_star(o: Proj) -> Proj:
    if isinstance(o, dict) and o.get("kind") == "GRanges":
        o["strand"] = ["*"] * len(o["strand"])
    return o


LADDER: list[tuple[str, Callable[[Proj], Proj]]] = [
    ("nan_as_na", nan_as_na),
    ("num_type", num_type),
    ("lgl_vs_num", lgl_as_int),
    ("names", drop_names),
    ("seqlevel_order", seqlevel_order),
    ("seqinfo_fields", seqinfo_unknowns),
    ("seqinfo", drop_seqinfo),
    ("mcol_order", mcol_order),
    ("mcol_types", mcols_types),
    ("mcols", drop_mcols),
    ("row_order", row_order),
    ("strand", strand_star),
]


def _equal_under(r: Proj, p: Proj, fs: list[Callable[[Proj], Proj]]) -> bool:
    r, p = copy.deepcopy(r), copy.deepcopy(p)
    for f in fs:
        r, p = f(r), f(p)
    return not diffs(r, p)


def decisive(r: Proj, p: Proj) -> str:
    """Name of the normalization (or pair) under which R and port agree; else `diff:<kind>@<path>`."""
    for name, f in LADDER:
        if _equal_under(r, p, [f]):
            return name
    for (n1, f1), (n2, f2) in itertools.combinations(LADDER, 2):
        if _equal_under(r, p, [f1, f2]):
            return f"{n1}+{n2}"
    d = diffs(r, p)[0]
    path = re.sub(r"\[\d+\]", "[]", d["path"])
    return f"diff:{d['kind']}@{path}"
