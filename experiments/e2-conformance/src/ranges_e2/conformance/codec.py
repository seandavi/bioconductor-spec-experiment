"""JSON tree <-> ranges_e2 objects, per conformance/FORMAT.md.

`Vec` carries R-typed atomic vectors (kind, values with None = NA, optional names) so adapters
can see exactly what R would see. Encoding produces plain JSON-compatible trees (dict/list/...)
with the same key order as the R codec.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import numpy as np
import polars as pl

from ..granges import STRAND_LEVELS, GRanges
from ..hits import Hits
from ..iranges import IRanges
from ..rle import Factor, Rle
from ..seqinfo import Seqinfo

Tree = Any
TAGS = {"$int": "int", "$dbl": "dbl", "$chr": "chr", "$lgl": "lgl"}
KIND_TAG = {v: k for k, v in TAGS.items()}


@dataclass
class Vec:
    """R atomic vector: kind in int|dbl|chr|lgl, values with None for NA, optional names."""

    kind: str
    values: list[Any]
    names: list[str | None] | None = None

    def __len__(self) -> int:
        return len(self.values)

    def has_na(self) -> bool:
        return any(v is None for v in self.values)

    def ints(self) -> np.ndarray:
        """int64 array (doubles truncate); NA is an error."""
        if self.kind not in ("int", "dbl", "lgl"):
            raise TypeError(f"expected numeric, got {self.kind}")
        if self.has_na():
            raise ValueError("NA not allowed")
        return np.array([int(v) for v in self.values], dtype=np.int64)

    def numeric(self) -> np.ndarray:
        """int64 for int vectors, float64 for double vectors."""
        if self.has_na():
            raise ValueError("NA not allowed")
        if self.kind == "dbl":
            return np.array([float(v) for v in self.values], dtype=np.float64)
        return self.ints()

    def strings(self) -> list[str | None]:
        if self.kind != "chr":
            raise TypeError(f"expected character, got {self.kind}")
        return list(self.values)

    def bools(self) -> np.ndarray:
        if self.kind != "lgl":
            raise TypeError(f"expected logical, got {self.kind}")
        if self.has_na():
            raise ValueError("NA not allowed")
        return np.array(self.values, dtype=bool)


@dataclass
class DFrameOut:
    """A metadata frame with explicit row count (polars cannot hold 0 columns with rows)."""

    frame: pl.DataFrame | None
    nrow: int


# ---- decode ----------------------------------------------------------------


def _num(v: Any) -> Any:
    if isinstance(v, str):
        return {"NaN": math.nan, "Inf": math.inf, "-Inf": -math.inf}[v]
    return v


def _shorthand(tree: Tree) -> Vec:
    items = tree if isinstance(tree, list) else [tree]
    if not items:
        raise ValueError("ambiguous empty array: use a tag")
    present = [v for v in items if v is not None]
    if not present:
        return Vec("lgl", list(items))
    if any(isinstance(v, (list, dict)) for v in present):
        raise ValueError("bad shorthand value")
    if all(isinstance(v, str) for v in present):
        return Vec("chr", list(items))
    if all(isinstance(v, bool) for v in present):
        return Vec("lgl", list(items))
    if all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in present):
        if any(isinstance(v, float) for v in present):
            return Vec("dbl", [None if v is None else float(v) for v in items])
        return Vec("int", list(items))
    raise ValueError("mixed kinds in shorthand array")


def _tagged(tree: dict[str, Any]) -> Vec:
    tag = next(k for k in tree if k in TAGS)
    kind = TAGS[tag]
    items = tree[tag]
    if kind == "dbl":
        values: list[Any] = [None if v is None else float(_num(v)) for v in items]
    elif kind == "int":
        values = [None if v is None else int(v) for v in items]
    else:
        values = list(items)
    names = tree.get("names")
    return Vec(kind, values, None if names is None else list(names))


_POLARS_DTYPES = {"int": pl.Int32, "dbl": pl.Float64, "chr": pl.Utf8, "lgl": pl.Boolean}


def vec_to_series(name: str, vec: Vec) -> pl.Series:
    return pl.Series(name, vec.values, dtype=_POLARS_DTYPES[vec.kind], strict=False)


def decode(tree: Tree, fixtures: dict[str, Any] | None = None) -> Any:
    fx = fixtures or {}
    if tree is None:
        return None
    if not isinstance(tree, dict):
        return _shorthand(tree)
    if "$ref" in tree:
        if tree["$ref"] not in fx:
            raise KeyError(f"unknown fixture {tree['$ref']}")
        return fx[tree["$ref"]]
    if any(k in TAGS for k in tree):
        return _tagged(tree)
    kind = tree.get("$type")
    if kind is None:
        raise ValueError("object without $type")

    def get(key: str) -> Any:
        return decode(tree[key], fx) if tree.get(key) is not None else None

    if kind == "IRanges":
        return _decode_iranges(tree, fx)
    if kind == "Seqinfo":
        return Seqinfo(
            get("seqnames").strings(),
            None if get("seqlengths") is None else get("seqlengths").values,
            None if get("isCircular") is None else get("isCircular").values,
            None if get("genome") is None else get("genome").values,
        )
    if kind == "GRanges":
        sn = get("seqnames").strings()
        rng = _decode_iranges(tree, fx)
        si = get("seqinfo")
        strand = get("strand")
        return GRanges(
            sn,
            rng,
            None if strand is None else strand.strings(),
            seqinfo=si if si is not None else Seqinfo(list(dict.fromkeys(sn))),
        )
    if kind == "DFrame":
        return _decode_dframe(tree, fx)
    if kind == "Hits":
        return Hits(
            np.array(get("from").values, dtype=np.int64) - 1,
            np.array(get("to").values, dtype=np.int64) - 1,
            int(tree["nLnode"]),
            int(tree["nRnode"]),
        )
    if kind == "Factor":
        lv = get("levels").strings()
        codes = np.array([-1 if c is None else c - 1 for c in get("codes").values], dtype=np.int64)
        return Factor(tuple(lv), codes)
    if kind == "Rle":
        values = get("values")
        lengths = np.array(get("lengths").values, dtype=np.int64)
        if isinstance(values, Vec):
            dtype = {"int": np.int32, "dbl": np.float64, "lgl": np.bool_}.get(values.kind)
            if dtype is None:
                raise TypeError("Rle of character is not supported")
            values = np.array(values.values, dtype=dtype)
        return Rle(values, lengths)
    if kind == "RleList":
        return {k: decode(v, fx) for k, v in tree["elements"].items()}
    raise ValueError(f"unknown $type {kind}")


def _decode_iranges(tree: dict[str, Any], fx: dict[str, Any]) -> IRanges:
    def arr(key: str) -> Any:
        v = decode(tree[key], fx) if tree.get(key) is not None else None
        return None if v is None else v.numeric()

    names = decode(tree["names"], fx) if tree.get("names") is not None else None
    mcols = decode(tree["mcols"], fx) if tree.get("mcols") is not None else None
    frame = mcols.frame if isinstance(mcols, DFrameOut) else mcols
    return IRanges(arr("start"), arr("end"), arr("width"), None if names is None else names.strings(), frame)


def _decode_dframe(tree: dict[str, Any], fx: dict[str, Any]) -> DFrameOut:
    cols = []
    for name, col in tree["columns"].items():
        if isinstance(col, dict) and col.get("$type") == "IntegerList":
            lists = [[int(i) - 1 for i in el] for el in col["elements"]]
            cols.append(pl.Series(name, lists, dtype=pl.List(pl.Int32)))
        else:
            v = decode(col, fx)
            if not isinstance(v, Vec):
                raise TypeError(f"unsupported DFrame column {name}")
            cols.append(vec_to_series(name, v))
    nrow = tree.get("nrow", len(cols[0]) if cols else 0)
    return DFrameOut(pl.DataFrame(cols) if cols else None, int(nrow))


# ---- encode ----------------------------------------------------------------


def _enc_item(kind: str, v: Any) -> Any:
    if v is None:
        return None
    if kind == "dbl":
        f = float(v)
        if math.isnan(f):
            return "NaN"
        if math.isinf(f):
            return "Inf" if f > 0 else "-Inf"
        return f
    if kind == "int":
        return int(v)
    if kind == "lgl":
        return bool(v)
    return v


def enc_vec(kind: str, values: Any, names: list[str | None] | None = None) -> Tree:
    out: dict[str, Any] = {KIND_TAG[kind]: [_enc_item(kind, v) for v in values]}
    if names is not None:
        out["names"] = list(names)
    return out


def _enc_names(names: list[str | None] | None) -> Tree:
    return None if names is None else enc_vec("chr", names)


def _enc_dframe(frame: pl.DataFrame | None, nrow: int) -> Tree:
    cols: dict[str, Any] = {}
    if frame is not None:
        for name in frame.columns:
            s = frame[name]
            cols[name] = _enc_column(name, s)
    return {"$type": "DFrame", "nrow": nrow, "columns": cols}


def _enc_column(name: str, s: pl.Series) -> Tree:
    dt = s.dtype
    if dt.is_integer():
        return enc_vec("int", s.to_list())
    if dt.is_float():
        return enc_vec("dbl", s.to_list())
    if dt == pl.Boolean:
        return enc_vec("lgl", s.to_list())
    if dt in (pl.Utf8, pl.String):
        return enc_vec("chr", s.to_list())
    if isinstance(dt, pl.List) and name == "revmap":
        return {"$type": "IntegerList", "elements": [[i + 1 for i in el] for el in s.to_list()], "names": None}
    raise TypeError(f"cannot encode column {name} of type {dt}")


def _enc_seqinfo(si: Seqinfo) -> Tree:
    return {
        "$type": "Seqinfo",
        "seqnames": enc_vec("chr", si.seqnames),
        "seqlengths": enc_vec("int", si.seqlengths),
        "isCircular": enc_vec("lgl", si.is_circular),
        "genome": enc_vec("chr", si.genome),
    }


def encode_factor(f: Factor) -> Tree:
    return {
        "$type": "Factor",
        "levels": enc_vec("chr", f.levels),
        "codes": enc_vec("int", [None if c < 0 else int(c) + 1 for c in f.codes.tolist()]),
    }


def encode(x: Any) -> Tree:
    if x is None:
        return None
    if isinstance(x, GRanges):
        return {
            "$type": "GRanges",
            "seqnames": enc_vec("chr", x.seqname_strings()),
            "start": enc_vec("int", x.start.tolist()),
            "width": enc_vec("int", x.width.tolist()),
            "strand": enc_vec("chr", [STRAND_LEVELS[c] for c in x._strand.tolist()]),
            "names": _enc_names(x.names),
            "mcols": _enc_dframe(x._ranges._mcols, len(x)),
            "seqinfo": _enc_seqinfo(x.seqinfo),
        }
    if isinstance(x, IRanges):
        return {
            "$type": "IRanges",
            "start": enc_vec("int", x.start.tolist()),
            "width": enc_vec("int", x.width.tolist()),
            "names": _enc_names(x.names),
            "mcols": _enc_dframe(x._mcols, len(x)),
        }
    if isinstance(x, Seqinfo):
        return _enc_seqinfo(x)
    if isinstance(x, Hits):
        return {
            "$type": "Hits",
            "from": enc_vec("int", (x.query_hits + 1).tolist()),
            "to": enc_vec("int", (x.subject_hits + 1).tolist()),
            "nLnode": x.query_length,
            "nRnode": x.subject_length,
        }
    if isinstance(x, Factor):
        return encode_factor(x)
    if isinstance(x, Rle):
        if isinstance(x.values, Factor):
            values = encode_factor(x.values)
        else:
            values = encode(x.values)
        return {"$type": "Rle", "values": values, "lengths": enc_vec("int", x.lengths.tolist())}
    if isinstance(x, dict) and all(isinstance(v, Rle) for v in x.values()):
        return {"$type": "RleList", "elements": {k: encode(v) for k, v in x.items()}}
    if isinstance(x, DFrameOut):
        return _enc_dframe(x.frame, x.nrow)
    if isinstance(x, Vec):
        return enc_vec(x.kind, x.values, x.names)
    if isinstance(x, np.ndarray):
        kind = {"i": "int", "u": "int", "f": "dbl", "b": "lgl"}.get(x.dtype.kind)
        if kind is None:
            raise TypeError(f"cannot encode array of dtype {x.dtype}")
        return enc_vec(kind, x.tolist())
    if isinstance(x, list):
        return enc_vec("chr", x)
    if isinstance(x, bool):
        return enc_vec("lgl", [x])
    if isinstance(x, (int, np.integer)):
        return enc_vec("int", [int(x)])
    if isinstance(x, float):
        return enc_vec("dbl", [x])
    if isinstance(x, str):
        return enc_vec("chr", [x])
    raise TypeError(f"cannot encode {type(x).__name__}")


# ---- comparison -------------------------------------------------------------


def _group_hits(tree: Tree) -> Tree:
    """hits-query-grouped: within each run of equal `from`, sort the `to` values."""
    if isinstance(tree, dict):
        if tree.get("$type") == "Hits":
            frm = tree["from"]["$int"]
            to = list(tree["to"]["$int"])
            i = 0
            while i < len(frm):
                j = i
                while j < len(frm) and frm[j] == frm[i]:
                    j += 1
                to[i:j] = sorted(to[i:j])
                i = j
            return {**tree, "to": {**tree["to"], "$int": to}}
        return {k: _group_hits(v) for k, v in tree.items()}
    if isinstance(tree, list):
        return [_group_hits(v) for v in tree]
    return tree


def tree_equal(a: Tree, b: Tree) -> bool:
    if a is None or b is None:
        return a is None and b is None
    if isinstance(a, (dict, list)) or isinstance(b, (dict, list)):
        if not (isinstance(a, (dict, list)) and isinstance(b, (dict, list))):
            return False
        if len(a) == 0 and len(b) == 0:
            return True
        if isinstance(a, dict) != isinstance(b, dict):
            return False
        if isinstance(a, dict):
            return list(a.keys()) == list(b.keys()) and all(tree_equal(a[k], b[k]) for k in a)
        return len(a) == len(b) and all(tree_equal(x, y) for x, y in zip(a, b))
    if isinstance(a, bool) or isinstance(b, bool):
        return isinstance(a, bool) and isinstance(b, bool) and a == b
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return abs(a - b) <= 1e-9 * max(1.0, abs(b))
    return bool(a == b)


def compare(actual: Tree, expected: Tree, mode: str = "exact") -> bool:
    if mode == "hits-query-grouped":
        actual, expected = _group_hits(actual), _group_hits(expected)
    elif mode != "exact":
        raise ValueError(f"unknown compare mode {mode}")
    return tree_equal(actual, expected)
