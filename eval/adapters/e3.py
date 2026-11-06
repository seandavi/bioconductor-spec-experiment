"""Adapter for port e3 (ranges_e3, "spec and ADRs first").

Port docs used: experiments/e3-spec-adr/README.md, spec/00-conventions.md (CONV),
10-iranges-core.md, 80-order.md, 90-seqinfo.md, 95-granges.md, and the docstrings
and public names under src/ranges_e3/.

Port index base: 0-based element indices (Hits, select results, order, nearest /
precede / follow, revmap, subscripts; CONV.3); coordinates are 1-based closed.

Translations (details in e3.md):
- R `a.b` argument names -> `a_b`; `maxgap`/`minoverlap` -> `max_gap`/`min_overlap`;
  `ties.method` -> `ties_method`; `fromLast` -> `from_last`; `subsetByOverlaps(ranges=)` -> `subject`.
- R generics -> methods on the first operand (`union(x, y)` -> `x.union(y)`); `c` -> `concat`.
- Index-valued results (`index_valued` in ops.json) are shifted +1 here, at the op
  boundary, and returned as RVector so nested calls see R-valued integers:
  `order`, `select != "all"` results, `Hits` (shifted in `project`), `revmap` list columns
  (an IntegerList column named `revmap` is 0-based inside the port; shifted in both directions).
- Positive integer subscripts are shifted -1. Negative / zero / NA integers, NA names, and
  logical masks that are not full length (R recycles) are `Unsupported`: the port documents
  Python negative indexing and refuses short masks.
- mcols are polars frames: atomic columns -> typed Series (NA -> null), Factor -> pl.Enum,
  IntegerList -> List(Int64).
- `seqlengths` / `isCircular` / `genome` results carry the seqlevels as names (R returns named vectors).
- Arguments the port's method does not accept (or required port arguments absent from the
  R call) are reported as unsupported by `_invoke`; nothing is filled in.
"""

from __future__ import annotations

import inspect
from typing import Any

import numpy as np
import polars as pl

from adapter_base import (
    MISSING,
    AdapterBase,
    RVector,
    Unsupported,
    UnknownArgument,
    p_dataframe,
    p_factor,
    p_granges,
    p_hits,
    p_integerlist,
    p_iranges,
    p_rle,
    p_rlelist,
    p_seqinfo,
    p_vector,
    p_vector_from,
)
from ranges_e3 import GRanges, Hits, IRanges, Rle, Seqinfo, concat

_POLARS_DTYPE = {"int": pl.Int64, "dbl": pl.Float64, "chr": pl.String, "lgl": pl.Boolean}
_REVMAP = "revmap"  # mcols column holding 0-based index lists in the port


# ---------------------------------------------------------------- helpers

def _val(v: Any) -> Any:
    """RVector -> scalar (length 1) or masked ndarray; natives pass through.

    A length-1 NA stays a masked array: `None` is R NULL / the port's "absent" sentinel.
    """
    if isinstance(v, RVector):
        return v.scalar() if len(v) == 1 and v.values[0] is not None else v.np_masked()
    return v


def _list(v: Any) -> Any:
    """RVector -> list (None = NA); natives pass through."""
    return v.list() if isinstance(v, RVector) else v


def _args(**kw: Any) -> dict:
    """Present arguments only, converted with `_val`."""
    return {k: _val(v) for k, v in kw.items() if v is not MISSING}


def _rename(kw: dict, names: dict[str, str]) -> dict:
    return {names.get(k, k): v for k, v in kw.items()}


def _rng(x: Any, op: str) -> Any:
    if not isinstance(x, (IRanges, GRanges)):
        raise Unsupported(f"{op}: port has no such operation on {type(x).__name__}")
    return x


def _invoke_fn(fn: Any, name: str, kw: dict) -> Any:
    """Call a port function with keyword arguments the port itself declares."""
    params = inspect.signature(fn).parameters
    for k in kw:
        if k not in params:
            raise UnknownArgument(k)
    for pname, p in params.items():
        if p.default is p.empty and pname not in kw:
            raise Unsupported(f"{name}: required port argument '{pname}' absent from the R call")
    return fn(**kw)


def _invoke(obj: Any, method: str, /, **kw: Any) -> Any:
    fn = getattr(obj, method, None)
    if fn is None or not callable(fn):
        raise Unsupported(f"{type(obj).__name__}.{method} does not exist")
    return _invoke_fn(fn, method, kw)


def _to_column(name: str, native: Any) -> pl.Series:
    if isinstance(native, RVector):
        return pl.Series(name, native.list(), dtype=_POLARS_DTYPE[native.type])
    if isinstance(native, pl.Series):
        if isinstance(native.dtype, pl.List) and name == _REVMAP:
            return pl.Series(
                name, [[None if v is None else v - 1 for v in e] for e in native.to_list()], dtype=native.dtype
            )
        return native.alias(name)
    if isinstance(native, np.ndarray):  # public accessor result used as metadata; keep dtype, mask -> null
        data = pl.Series(name, np.ma.getdata(native))
        mask = np.ma.getmaskarray(native)
        return data.zip_with(pl.Series(~mask), pl.Series(name, [None] * len(data), dtype=data.dtype))
    raise Unsupported(f"mcols column of {type(native).__name__}")


def _frame(columns: list[tuple[str, Any]]) -> pl.DataFrame | None:
    if not columns:
        return None
    names = [n for n, _ in columns]
    if len(set(names)) != len(names):
        raise Unsupported("mcols: duplicate column names (polars frame)")
    return pl.DataFrame([_to_column(n, v) for n, v in columns])


def _column_projection(name: str, s: pl.Series) -> dict:
    dt = s.dtype
    if isinstance(dt, pl.Enum):
        return p_factor(dt.categories.to_list(), s.to_physical().to_list(), base=0)
    if isinstance(dt, pl.List):
        return p_integerlist(s.to_list(), shift=1 if name == _REVMAP else 0)
    if dt.is_integer():
        return p_vector("int", s.to_list())
    if dt.is_float():
        return p_vector("dbl", s.to_list())
    if dt == pl.Boolean:
        return p_vector("lgl", s.to_list())
    if dt == pl.String:
        return p_vector("chr", s.to_list())
    raise ValueError(f"cannot project polars dtype {dt}")


def _mcols_projection(df: pl.DataFrame) -> list:
    return [(c, _column_projection(c, df[c])) for c in df.columns]


def _index_result(res: Any) -> Any:
    """0-based port index result -> R 1-based. Hits stay native (shifted in project)."""
    if isinstance(res, Hits):
        return res
    mask = np.ma.getmaskarray(res)
    return RVector("int", [None if m else int(v) + 1 for v, m in zip(np.ma.getdata(res).tolist(), mask.tolist())])


def _subscript(i: Any, n: int) -> Any:
    """R subscript -> port subscript, for the forms whose port semantics are documented (CORE.9)."""
    if i is MISSING or i is None:
        raise Unsupported("subscript:missing")
    if isinstance(i, RVector):
        t, vals = i.type, i.list()
    elif isinstance(i, np.ndarray):
        t = {"i": "int", "u": "int", "f": "dbl", "b": "lgl", "U": "chr", "O": "chr"}[i.dtype.kind]
        vals = [None if m else v for v, m in zip(np.ma.getdata(i).tolist(), np.ma.getmaskarray(i).tolist())]
    else:
        raise Unsupported(f"subscript:{type(i).__name__}")
    if t == "dbl":
        raise Unsupported("subscript:double (port documents integer sequences)")
    if t == "int":
        if any(v is None or v <= 0 for v in vals):
            raise Unsupported("subscript:non-positive or NA integer")
        return np.array(vals, dtype=np.int64) - 1
    if t == "lgl":
        if any(v is None for v in vals) or len(vals) != n:
            raise Unsupported("subscript:logical with NA or not full length (R recycles)")
        return np.array(vals, dtype=np.bool_)
    if any(v is None for v in vals):
        raise Unsupported("subscript:NA name")
    return list(vals)


# ---------------------------------------------------------------- adapter

class Adapter(AdapterBase):
    port = "e3"
    port_packages = ("ranges_e3",)
    unsupported_ops = {
        "seqinfo<-": "no public seqinfo setter",
        "seqlevels<-": "no public seqlevels setter",
        "seqlengths<-": "no public seqlengths setter",
        "isCircular<-": "no public is_circular setter",
        "genome<-": "no public genome setter",
    }

    # ---- inputs
    def build_IRanges(self, p: dict) -> IRanges:
        return IRanges(
            start=np.array(p["start"], dtype=np.int64),
            width=np.array(p["width"], dtype=np.int64),
            names=p["names"],
            mcols=self._input_mcols(p),
        )

    def build_GRanges(self, p: dict) -> GRanges:
        return GRanges(
            p["seqnames"],
            IRanges(start=np.array(p["start"], dtype=np.int64), width=np.array(p["width"], dtype=np.int64)),
            strand=p["strand"],
            seqinfo=self.build_Seqinfo(p["seqinfo"]),
            names=p["names"],
            mcols=self._input_mcols(p),
        )

    def build_Seqinfo(self, p: dict) -> Seqinfo:
        return Seqinfo(p["seqnames"], p["seqlengths"], p["isCircular"], p["genome"])

    def build_DataFrame(self, p: dict) -> pl.DataFrame:
        if p["rownames"] is not None:
            raise Unsupported("input:DataFrame rownames")
        if not p["columns"] and p["nrow"]:
            raise Unsupported("input:DataFrame with rows and no columns (polars frame is 0 x 0)")
        return _frame([(n, self.build_input(c)) for n, c in p["columns"]]) if p["columns"] else pl.DataFrame()

    def build_Factor(self, p: dict) -> pl.Series:
        lv = p["levels"]
        return pl.Series("", [None if c is None else lv[c - 1] for c in p["codes"]], dtype=pl.Enum(lv))

    def build_IntegerList(self, p: dict) -> pl.Series:
        if p["names"] is not None:
            raise Unsupported("input:IntegerList names")
        return pl.Series("", p["elements"], dtype=pl.List(pl.Int64))

    def _input_mcols(self, p: dict) -> pl.DataFrame | None:
        return _frame([(n, self.build_input(c)) for n, c in p["mcols"]])

    # ---- construct
    def op_IRanges(self, *, start=MISSING, end=MISSING, width=MISSING, names=MISSING, _mcols=None):
        kw = _args(start=start, end=end, width=width)
        if names is not MISSING:
            kw["names"] = _list(names)
        if _mcols:
            kw["mcols"] = _frame(list(_mcols))
        return _invoke_fn(IRanges, "IRanges", kw)

    def op_GRanges(self, *, seqnames=MISSING, ranges=MISSING, strand=MISSING, seqinfo=MISSING, _mcols=None):
        kw = {k: v for k, v in (("seqnames", seqnames), ("ranges", ranges), ("seqinfo", seqinfo)) if v is not MISSING}
        if "seqnames" in kw:
            kw["seqnames"] = _list(seqnames)
        if strand is not MISSING:
            kw["strand"] = _list(strand)
        if _mcols:
            kw["mcols"] = _frame(list(_mcols))
        return _invoke_fn(GRanges, "GRanges", kw)

    def op_Seqinfo(self, *, seqnames=MISSING, seqlengths=MISSING, isCircular=MISSING, genome=MISSING):
        kw = {"seqnames": seqnames, "seqlengths": seqlengths, "is_circular": isCircular, "genome": genome}
        return _invoke_fn(Seqinfo, "Seqinfo", {k: _list(v) for k, v in kw.items() if v is not MISSING})

    # ---- accessors
    def op_start(self, *, x):
        return _rng(x, "start").start

    def op_end(self, *, x):
        return _rng(x, "end").end

    def op_width(self, *, x):
        return _rng(x, "width").width

    def op_names(self, *, x):
        return _rng(x, "names").names

    def op_length(self, *, x):
        if not isinstance(x, (IRanges, GRanges, Seqinfo)):
            raise Unsupported(f"length: {type(x).__name__}")
        return RVector("int", [len(x)])

    def op_seqnames(self, *, x):
        if not isinstance(x, GRanges):
            raise Unsupported(f"seqnames: {type(x).__name__}")
        return x.seqnames

    def op_strand(self, *, x):
        if not isinstance(x, GRanges):
            raise Unsupported(f"strand: {type(x).__name__}")
        return x.strand

    def op_ranges(self, *, x):
        if not isinstance(x, GRanges):
            raise Unsupported(f"ranges: {type(x).__name__}")
        return x.ranges

    def op_mcols(self, *, x):
        return _rng(x, "mcols").mcols

    # ---- metadata
    def op_set_names(self, *, x, value=MISSING):
        if value is MISSING:
            raise Unsupported("names<-: value absent")
        return _invoke(_rng(x, "names<-"), "with_names", names=_list(value))

    def op_set_mcols(self, *, x, value=MISSING):
        if value is MISSING:
            raise Unsupported("mcols<-: value absent")
        if value is not None and not isinstance(value, pl.DataFrame):
            raise Unsupported(f"mcols<-: value {type(value).__name__}")
        return _invoke(_rng(x, "mcols<-"), "with_mcols", mcols=value)

    # ---- seqinfo
    def _seq_attr(self, x: Any, op: str, gr_attr: str, si_attr: str) -> list:
        if isinstance(x, GRanges):
            return getattr(x, gr_attr)
        if isinstance(x, Seqinfo):
            return getattr(x, si_attr)
        raise Unsupported(f"{op}: {type(x).__name__}")

    def _levels(self, x: Any, op: str) -> list:
        return self._seq_attr(x, op, "seqlevels", "seqnames")

    def op_seqinfo(self, *, x):
        if not isinstance(x, GRanges):
            raise Unsupported(f"seqinfo: {type(x).__name__}")
        return x.seqinfo

    def op_seqlevels(self, *, x):
        return RVector("chr", self._levels(x, "seqlevels"))

    def op_seqlengths(self, *, x):
        return RVector("int", self._seq_attr(x, "seqlengths", "seqlengths", "seqlengths"), self._levels(x, "seqlengths"))

    def op_isCircular(self, *, x):
        return RVector("lgl", self._seq_attr(x, "isCircular", "is_circular", "is_circular"), self._levels(x, "isCircular"))

    def op_genome(self, *, x):
        return RVector("chr", self._seq_attr(x, "genome", "genome", "genome"), self._levels(x, "genome"))

    def op_merge(self, *, x, y=MISSING):
        if not isinstance(x, Seqinfo):
            raise Unsupported(f"merge: {type(x).__name__}")
        if y is MISSING:
            raise Unsupported("merge: y absent")
        return x.merge(y)

    def op_trim(self, *, x, **kw):
        return _invoke(_rng(x, "trim"), "trim", **_args(**kw))

    # ---- subset / concat
    def op_subset(self, *, x, i=MISSING):
        x = _rng(x, "[")
        return x[_subscript(i, len(x))]

    def op_c(self, *dots, x, **kw):
        for k in kw:
            raise UnknownArgument(k)
        objs = (x, *dots)
        for o in objs:
            _rng(o, "c")
        return concat(*objs)

    # ---- intra-range
    def op_shift(self, *, x, **kw):
        return _invoke(_rng(x, "shift"), "shift", **_args(**kw))

    def op_narrow(self, *, x, **kw):
        return _invoke(_rng(x, "narrow"), "narrow", **_args(**kw))

    def op_resize(self, *, x, **kw):
        fix = kw.get("fix")
        if isinstance(fix, RVector) and len(fix) > 1:
            raise Unsupported("argument:fix as a per-range vector (port takes one string)")
        return _invoke(_rng(x, "resize"), "resize", **_args(**kw))

    def op_flank(self, *, x, **kw):
        return _invoke(_rng(x, "flank"), "flank", **_args(**kw))

    def op_promoters(self, *, x, **kw):
        return _invoke(_rng(x, "promoters"), "promoters", **_args(**kw))

    # ---- inter-range
    def op_range(self, *, x, **kw):
        return _invoke(_rng(x, "range"), "range", **_args(**kw))

    def op_reduce(self, *, x, **kw):
        return _invoke(_rng(x, "reduce"), "reduce", **_args(**kw))

    def op_gaps(self, *, x, **kw):
        return _invoke(_rng(x, "gaps"), "gaps", **_args(**kw))

    def op_disjoin(self, *, x, **kw):
        return _invoke(_rng(x, "disjoin"), "disjoin", **_args(**kw))

    # ---- set operations
    def op_union(self, *, x, **kw):
        return _invoke(_rng(x, "union"), "union", **_args(**kw))

    def op_intersect(self, *, x, **kw):
        return _invoke(_rng(x, "intersect"), "intersect", **_args(**kw))

    def op_setdiff(self, *, x, **kw):
        return _invoke(_rng(x, "setdiff"), "setdiff", **_args(**kw))

    def op_punion(self, *, x, **kw):
        return _invoke(_rng(x, "punion"), "punion", **_args(**kw))

    def op_pintersect(self, *, x, **kw):
        return _invoke(_rng(x, "pintersect"), "pintersect", **_args(**kw))

    def op_psetdiff(self, *, x, **kw):
        return _invoke(_rng(x, "psetdiff"), "psetdiff", **_args(**kw))

    # ---- overlaps
    _OVERLAP_NAMES = {"maxgap": "max_gap", "minoverlap": "min_overlap"}

    def op_findOverlaps(self, *, query, **kw):
        res = _invoke(_rng(query, "findOverlaps"), "find_overlaps", **_args(**_rename(kw, self._OVERLAP_NAMES)))
        return _index_result(res)

    def op_countOverlaps(self, *, query, **kw):
        return _invoke(_rng(query, "countOverlaps"), "count_overlaps", **_args(**_rename(kw, self._OVERLAP_NAMES)))

    def op_overlapsAny(self, *, query, **kw):
        return _invoke(_rng(query, "overlapsAny"), "overlaps_any", **_args(**_rename(kw, self._OVERLAP_NAMES)))

    def op_subsetByOverlaps(self, *, x, **kw):
        names = {**self._OVERLAP_NAMES, "ranges": "subject"}
        return _invoke(_rng(x, "subsetByOverlaps"), "subset_by_overlaps", **_args(**_rename(kw, names)))

    # ---- nearest
    def op_nearest(self, *, x, **kw):
        return _index_result(_invoke(_rng(x, "nearest"), "nearest", **_args(**kw)))

    def op_precede(self, *, x, **kw):
        return _index_result(_invoke(_rng(x, "precede"), "precede", **_args(**kw)))

    def op_follow(self, *, x, **kw):
        return _index_result(_invoke(_rng(x, "follow"), "follow", **_args(**kw)))

    def op_distance(self, *, x, **kw):
        return _invoke(_rng(x, "distance"), "distance", **_args(**kw))

    # ---- coverage
    def op_coverage(self, *, x, **kw):
        return _invoke(_rng(x, "coverage"), "coverage", **_args(**kw))

    # ---- ordering
    def op_order(self, *dots, **kw):
        if len(dots) != 1:
            raise Unsupported("order: port orders one object (no secondary keys)")
        res = _invoke(_rng(dots[0], "order"), "order", **_args(**kw))
        return RVector("int", [int(v) + 1 for v in res.tolist()])

    def op_sort(self, *, x, **kw):
        return _invoke(_rng(x, "sort"), "sort", **_args(**kw))

    def op_rank(self, *, x, **kw):
        return _invoke(_rng(x, "rank"), "rank", **_args(**kw))

    def op_duplicated(self, *, x, **kw):
        return _invoke(_rng(x, "duplicated"), "duplicated", **_args(**_rename(kw, {"fromLast": "from_last"})))

    def op_unique(self, *, x, **kw):
        return _invoke(_rng(x, "unique"), "unique", **_args(**_rename(kw, {"fromLast": "from_last"})))

    # ---- projection
    def project(self, value: Any) -> dict:
        if isinstance(value, GRanges):
            return p_granges(
                value.seqnames, value.start, value.width, value.strand, value.names,
                _mcols_projection(value.mcols), self.project(value.seqinfo),
            )
        if isinstance(value, IRanges):
            return p_iranges(
                value.start, value.width, value.names, _mcols_projection(value.mcols)
            )
        if isinstance(value, Seqinfo):
            return p_seqinfo(value.seqnames, value.seqlengths, value.is_circular, value.genome)
        if isinstance(value, Hits):
            return p_hits(value.query_hits + 1, value.subject_hits + 1, value.query_length, value.subject_length)
        if isinstance(value, Rle):
            return p_rle(value.values, value.lengths)
        if isinstance(value, dict):  # GRanges.coverage: seqlevel -> Rle
            return p_rlelist(list(value), [(r.values, r.lengths) for r in value.values()])
        if isinstance(value, pl.DataFrame):
            return p_dataframe(value.height, _mcols_projection(value))
        if isinstance(value, pl.Series):
            return _column_projection("", value)
        if isinstance(value, np.ndarray):
            return p_vector_from(value)
        if isinstance(value, list):  # names
            return p_vector("chr", value)
        raise TypeError(f"cannot project {type(value).__name__}")
