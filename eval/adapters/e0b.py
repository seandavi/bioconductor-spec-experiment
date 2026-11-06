"""Adapter for port e0b (ranges_e0b). Thin: translation and projection only.

Port docs used: experiments/e0b-no-r/README.md, docstrings and public names in src/ranges_e0b.
Port index base: 0-based for positional indices (subscripts, Hits, select= results, order, revmap);
NO_HIT (-1) is the port's NA for index results. Coordinates are 1-based like R.
Translations (see e0b.md for the unsupported list):
  - R `a.b` argument names -> port `a_b` (done by the base); `y`/`subject`/`ranges`/`value` -> the
    port's `other` / `subject` / `with_*` parameter names.
  - R ops are generics, the port's are methods: x.reduce(...), x.find_overlaps(...), x.concat(...).
  - 1-length RVector -> python scalar, longer -> ndarray (numeric/logical) or list (character).
  - positive integer subscripts shift to 0-based; Hits / select / order results shift back to 1-based;
    NO_HIT -> NA. Other subscript forms are Unsupported (the port's negative/zero meaning differs).
  - mcols are polars frames: int/dbl/chr/lgl -> Int64/Float64/String/Boolean, Factor -> Enum,
    IntegerList -> List(Int64).
  - reduce(with_revmap=True) returns (ranges, list[ndarray]); projected as a `revmap` mcols column (+1).
  - GRanges.seqlengths/is_circular/genome dicts and Seqinfo tuples -> named vectors (names = seqlevels).
  - GRanges.coverage dict[str, Rle] -> RleList; IRanges.coverage Rle -> Rle.
  - distance: NO_HIT (-1, documented as R NA in GRanges.distance) -> NA.
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
    UnknownArgument,
    Unsupported,
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
from ranges_e0b import NO_HIT, GRanges, Hits, IRanges, Rle, Seqinfo

_PL_DTYPE = {"int": pl.Int64, "dbl": pl.Float64, "chr": pl.String, "lgl": pl.Boolean}


# ---------------------------------------------------------------- argument translation

def _arg(v: Any, name: str = "argument") -> Any:
    """RVector -> scalar (length 1; NA -> None) / ndarray (numeric, logical) / list (character)."""
    if not isinstance(v, RVector):
        return v
    if len(v) == 1:
        return v.scalar()
    if v.type == "chr":
        return v.list()
    try:
        return v.np()
    except ValueError:
        raise Unsupported(f"NA in vector argument {name}") from None


def _list(v: Any) -> list | None:
    """Character vector -> list (NULL -> None)."""
    return None if v is None else v.list()


def _scalar_or_list(v: Any) -> Any:
    """Seqinfo fields: the port broadcasts a scalar; longer vectors are plain lists."""
    if v is None:
        return None
    return v.scalar() if len(v) == 1 else v.list()


def _mapping_or_arg(v: Any, name: str) -> Any:
    """A named R vector is the port's `Mapping[str, int]` (per seqlevel); otherwise as `_arg`."""
    if isinstance(v, RVector) and v.names is not None:
        return dict(zip(v.names, v.list()))
    return _arg(v, name)


def _series(v: RVector) -> pl.Series:
    if v.names is not None:
        raise Unsupported("named vector as a metadata column")
    return pl.Series(v.list(), dtype=_PL_DTYPE[v.type])


def _mcols_mapping(pairs: Any) -> dict[str, pl.Series] | None:
    """`...` extra named arguments -> mcols mapping (the port builds the frame)."""
    if not pairs:
        return None
    out = {}
    for name, value in pairs:
        if not isinstance(value, RVector):
            raise Unsupported("metadata column that is not an atomic vector")
        out[name] = _series(value)
    return out


def _invoke(obj: Any, method: str, /, **kw: Any) -> Any:
    """Call obj.method(**kw) with port parameter names; an R argument the port lacks is unsupported."""
    fn = getattr(obj, method, None)
    if fn is None:
        raise Unsupported(f"{type(obj).__name__}.{method}: no such public method")
    params = inspect.signature(fn).parameters
    for k in kw:
        if k not in params:
            raise UnknownArgument(k)
    required = [p.name for p in params.values() if p.default is p.empty and p.name not in kw
                and p.kind in (p.POSITIONAL_OR_KEYWORD, p.KEYWORD_ONLY)]
    if required:
        raise Unsupported(f"{method}: port requires {', '.join(required)}")
    return fn(**{k: _arg(v, k) for k, v in kw.items()})


def _pass(obj: Any, method: str, kw: dict, /, **renamed: Any) -> Any:
    """_invoke with R arguments renamed to the port's parameter names (MISSING = absent)."""
    return _invoke(obj, method, **kw, **{k: v for k, v in renamed.items() if v is not MISSING})


def _prop(obj: Any, name: str) -> Any:
    try:
        return getattr(obj, name)
    except AttributeError:
        raise Unsupported(f"{type(obj).__name__}.{name}: no such public attribute") from None


def _one_based(values: Any) -> RVector:
    """0-based port indices (NO_HIT = NA) -> R integer vector."""
    return RVector("int", [None if v == NO_HIT else int(v) + 1 for v in np.asarray(values).tolist()])


def _hits_or_indices(res: Any) -> Any:
    return res if isinstance(res, Hits) else _one_based(res)


def _subscript(i: Any, n: int) -> Any:
    """R subscript -> what the port documents: 0-based int array, bool mask of length n, or names."""
    if not isinstance(i, RVector):
        raise Unsupported("subscript:missing or NULL")
    if i.type == "int":
        vals = i.list()
        if any(v is None or v <= 0 for v in vals):
            raise Unsupported("subscript:non-positive or NA integer")
        return np.array([v - 1 for v in vals], dtype=np.int64)
    if i.type == "lgl":
        vals = i.list()
        if any(v is None for v in vals):
            raise Unsupported("subscript:NA logical")
        if len(vals) != n:
            raise Unsupported("subscript:logical recycling")
        return np.array(vals, dtype=bool)
    if i.type == "chr":
        if any(v is None for v in i.list()):
            raise Unsupported("subscript:NA name")
        return i.list()
    raise Unsupported("subscript:double")


# ---------------------------------------------------------------- projection helpers

def _column(s: pl.Series) -> dict:
    dt = s.dtype
    if isinstance(dt, pl.Enum):
        return p_factor(dt.categories.to_list(), s.to_physical().to_list(), base=0)
    if isinstance(dt, pl.List) and dt.inner is not None and dt.inner.is_integer():
        return p_integerlist(s.to_list())
    for kind, test in (("int", dt.is_integer()), ("dbl", dt.is_float()), ("chr", dt == pl.String),
                       ("lgl", dt == pl.Boolean)):
        if test:
            return p_vector(kind, s.to_list())
    raise ValueError(f"cannot project polars dtype {dt}")


def _mcols_of(df: pl.DataFrame) -> list:
    return [(c, _column(df[c])) for c in df.columns]


def _project_seqinfo(si: Seqinfo) -> dict:
    return p_seqinfo(si.seqnames, si.seqlengths, si.is_circular, si.genome)


def _project_ranges(value: Any, revmap: Any = None) -> dict:
    mcols = _mcols_of(value.mcols)
    if revmap is not None:
        mcols.append(("revmap", p_integerlist(revmap, shift=1)))  # 0-based -> 1-based
    if isinstance(value, GRanges):
        return p_granges(value.seqnames, value.start, value.width, value.strand, value.names, mcols,
                         _project_seqinfo(value.seqinfo))
    return p_iranges(value.start, value.width, value.names, mcols)


# ---------------------------------------------------------------- adapter

class Adapter(AdapterBase):
    port = "e0b"
    port_packages = ("ranges_e0b",)
    unsupported_ops = {
        "seqlevels<-": "no public seqlevels setter (only GRanges.with_seqinfo)",
        "seqlengths<-": "no public seqlengths setter (only GRanges.with_seqinfo)",
        "isCircular<-": "no public isCircular setter (only GRanges.with_seqinfo)",
        "genome<-": "no public genome setter (only GRanges.with_seqinfo)",
    }

    # ---- inputs ------------------------------------------------------------
    @staticmethod
    def _column_series(col: dict) -> pl.Series:
        kind = col["kind"]
        if kind == "vector":
            return _series(RVector(col["type"], col["values"], col["names"]))
        if col.get("names") is not None:
            raise Unsupported(f"named {kind} column")
        if kind == "Factor":
            levels = col["levels"]
            return pl.Series([None if c is None else levels[c - 1] for c in col["codes"]], dtype=pl.Enum(levels))
        if kind == "IntegerList":
            return pl.Series(col["elements"], dtype=pl.List(pl.Int64))
        raise Unsupported(f"metadata column:{kind}")

    def _mcols_input(self, pairs: list) -> dict[str, pl.Series] | None:
        return {name: self._column_series(col) for name, col in pairs} if pairs else None

    def _iranges_input(self, p: dict) -> IRanges:
        return IRanges(start=np.array(p["start"], dtype=np.int64), width=np.array(p["width"], dtype=np.int64),
                       names=p["names"], mcols=self._mcols_input(p["mcols"]))

    def build_IRanges(self, p: dict) -> IRanges:
        return self._iranges_input(p)

    def build_Seqinfo(self, p: dict) -> Seqinfo:
        return Seqinfo(p["seqnames"], p["seqlengths"], p["isCircular"], p["genome"])

    def build_GRanges(self, p: dict) -> GRanges:
        ranges = IRanges(start=np.array(p["start"], dtype=np.int64), width=np.array(p["width"], dtype=np.int64),
                         names=p["names"])
        return GRanges(p["seqnames"], ranges, p["strand"], self._mcols_input(p["mcols"]),
                       self.build_Seqinfo(p["seqinfo"]))

    def build_DataFrame(self, p: dict) -> pl.DataFrame:
        if p["rownames"] is not None:
            raise Unsupported("DataFrame rownames")
        if not p["columns"] and p["nrow"]:
            raise Unsupported("DataFrame with rows but no columns")
        return pl.DataFrame([self._column_series(col).alias(name) for name, col in p["columns"]])

    def build_Factor(self, p: dict) -> pl.Series:
        return self._column_series(p)

    # ---- constructors --------------------------------------------------------
    def op_IRanges(self, *, start=MISSING, end=MISSING, width=MISSING, names=MISSING, _mcols=None):
        kw = {k: v for k, v in (("start", start), ("end", end), ("width", width)) if v is not MISSING}
        if names is not MISSING:
            kw["names"] = _list(names)
        return IRanges(**{k: _arg(v, k) for k, v in kw.items()}, mcols=_mcols_mapping(_mcols))

    def op_GRanges(self, *, seqnames=MISSING, ranges=MISSING, strand=MISSING, seqinfo=MISSING,
                   seqlengths=MISSING, _mcols=None):
        if seqlengths is not MISSING:
            raise UnknownArgument("seqlengths")
        if seqnames is MISSING or ranges is MISSING:
            raise Unsupported("GRanges: port requires seqnames and ranges")
        if isinstance(seqnames, pl.Series):
            raise Unsupported("argument:seqnames as a factor (the port has no factor seqnames)")
        kw = {} if strand is MISSING else {"strand": _arg(strand, "strand")}
        if seqinfo is not MISSING:
            kw["seqinfo"] = seqinfo
        return GRanges(_arg(seqnames, "seqnames"), ranges, mcols=_mcols_mapping(_mcols), **kw)

    def op_Seqinfo(self, *, seqnames=MISSING, seqlengths=MISSING, isCircular=MISSING, genome=MISSING):
        if seqnames is MISSING:
            raise Unsupported("Seqinfo: port requires seqnames")
        kw = {k: _scalar_or_list(v) for k, v in (("seqlengths", seqlengths), ("is_circular", isCircular),
                                                  ("genome", genome)) if v is not MISSING}
        return Seqinfo(_list(seqnames), **kw)

    # ---- accessors ---------------------------------------------------------------
    def op_start(self, *, x):
        return _prop(x, "start")

    def op_end(self, *, x):
        return _prop(x, "end")

    def op_width(self, *, x):
        return _prop(x, "width")

    def op_names(self, *, x):
        names = _prop(x, "names")
        return None if names is None else RVector("chr", names)

    def op_length(self, *, x):
        return RVector("int", [len(x)])

    def op_seqnames(self, *, x):
        return RVector("chr", _prop(x, "seqnames"))

    def op_strand(self, *, x):
        return RVector("chr", _prop(x, "strand"))

    def op_ranges(self, *, x):
        return _prop(x, "ranges")

    def op_mcols(self, *, x):
        return _prop(x, "mcols")

    # ---- metadata -----------------------------------------------------------------
    def op_set_names(self, *, x, value=MISSING):
        return _pass(x, "with_names", {}, names=value if value is MISSING else _list(value))

    def op_set_mcols(self, *, x, value=MISSING):
        return _pass(x, "with_mcols", {}, mcols=value)

    # ---- seqinfo --------------------------------------------------------------------
    def op_seqinfo(self, *, x):
        return _prop(x, "seqinfo")

    def op_seqlevels(self, *, x):
        if isinstance(x, Seqinfo):
            return RVector("chr", list(x.seqnames))
        return RVector("chr", _prop(x, "seqlevels"))

    @staticmethod
    def _per_seqlevel(x: Any, gr_attr: str, si_attr: str, type_: str) -> RVector:
        """Named R vector over seqlevels from GRanges' dict or Seqinfo's tuple."""
        if isinstance(x, Seqinfo):
            return RVector(type_, list(getattr(x, si_attr)), list(x.seqnames))
        d = _prop(x, gr_attr)
        return RVector(type_, list(d.values()), list(d.keys()))

    def op_seqlengths(self, *, x):
        return self._per_seqlevel(x, "seqlengths", "seqlengths", "int")

    def op_isCircular(self, *, x):
        return self._per_seqlevel(x, "is_circular", "is_circular", "lgl")

    def op_genome(self, *, x):
        return self._per_seqlevel(x, "genome", "genome", "chr")

    def op_set_seqinfo(self, *, x, value=MISSING, **kw):
        return _pass(x, "with_seqinfo", kw, seqinfo=value)

    def op_merge(self, *, x, y=MISSING):
        return _pass(x, "merge", {}, other=y)

    def op_trim(self, *, x, **kw):
        return _invoke(x, "trim", **kw)

    # ---- subsetting and combining -------------------------------------------------------
    def op_subset(self, *, x, i=MISSING):
        if not isinstance(x, (IRanges, GRanges)):
            raise Unsupported(f"subset of {type(x).__name__}")
        return x[_subscript(i, len(x))]

    def op_c(self, *dots, x, **kw):
        for k in kw:
            raise UnknownArgument(k)
        return _prop(x, "concat")(*dots)

    # ---- intra-range ---------------------------------------------------------------------
    def op_shift(self, *, x, **kw):
        return _invoke(x, "shift", **kw)

    def op_narrow(self, *, x, **kw):
        return _invoke(x, "narrow", **kw)

    def op_resize(self, *, x, **kw):
        return _invoke(x, "resize", **kw)

    def op_flank(self, *, x, **kw):
        return _invoke(x, "flank", **kw)

    def op_promoters(self, *, x, **kw):
        return _invoke(x, "promoters", **kw)

    # ---- inter-range -------------------------------------------------------------------------
    def op_range(self, *, x, **kw):
        return _invoke(x, "range", **kw)

    def op_reduce(self, *, x, **kw):
        return _invoke(x, "reduce", **kw)

    def op_gaps(self, *, x, start=MISSING, end=MISSING, **kw):
        bounds = {k: _mapping_or_arg(v, k) for k, v in (("start", start), ("end", end)) if v is not MISSING}
        return _invoke(x, "gaps", **kw, **bounds)

    def op_disjoin(self, *, x, **kw):
        return _invoke(x, "disjoin", **kw)

    # ---- set operations ------------------------------------------------------------------------
    def op_union(self, *, x, y=MISSING, **kw):
        return _pass(x, "union", kw, other=y)

    def op_intersect(self, *, x, y=MISSING, **kw):
        return _pass(x, "intersect", kw, other=y)

    def op_setdiff(self, *, x, y=MISSING, **kw):
        return _pass(x, "setdiff", kw, other=y)

    def op_punion(self, *, x, y=MISSING, **kw):
        return _pass(x, "punion", kw, other=y)

    def op_pintersect(self, *, x, y=MISSING, **kw):
        return _pass(x, "pintersect", kw, other=y)

    def op_psetdiff(self, *, x, y=MISSING, **kw):
        return _pass(x, "psetdiff", kw, other=y)

    # ---- overlaps ----------------------------------------------------------------------------------
    def op_findOverlaps(self, *, query, subject=MISSING, **kw):
        return _hits_or_indices(_pass(query, "find_overlaps", kw, subject=subject))

    def op_countOverlaps(self, *, query, subject=MISSING, **kw):
        return _pass(query, "count_overlaps", kw, subject=subject)

    def op_overlapsAny(self, *, query, subject=MISSING, **kw):
        return _pass(query, "overlaps_any", kw, subject=subject)

    def op_subsetByOverlaps(self, *, x, ranges=MISSING, **kw):
        return _pass(x, "subset_by_overlaps", kw, subject=ranges)

    # ---- nearest ---------------------------------------------------------------------------------------
    def op_nearest(self, *, x, subject=MISSING, **kw):
        return _hits_or_indices(_pass(x, "nearest", kw, subject=subject))

    def op_precede(self, *, x, subject=MISSING, **kw):
        return _hits_or_indices(_pass(x, "precede", kw, subject=subject))

    def op_follow(self, *, x, subject=MISSING, **kw):
        return _hits_or_indices(_pass(x, "follow", kw, subject=subject))

    def op_distance(self, *, x, y=MISSING, **kw):
        res = _pass(x, "distance", kw, other=y)
        return RVector("int", [None if v == NO_HIT else v for v in res.tolist()])  # NO_HIT is R NA here

    # ---- coverage ------------------------------------------------------------------------------------------
    def op_coverage(self, *, x, width=MISSING, **kw):
        bounds = {} if width is MISSING else {"width": _mapping_or_arg(width, "width")}
        return _invoke(x, "coverage", **kw, **bounds)

    # ---- ordering ------------------------------------------------------------------------------------------
    def op_order(self, *dots, **kw):
        if len(dots) != 1:
            raise Unsupported("order: the port orders one object")
        return _one_based(_invoke(dots[0], "order", **kw))

    def op_sort(self, *, x, **kw):
        return _invoke(x, "sort", **kw)

    def op_rank(self, *, x, **kw):
        return _invoke(x, "rank", **kw)

    def op_duplicated(self, *, x, **kw):
        return _invoke(x, "duplicated", **kw)

    def op_unique(self, *, x, **kw):
        return _invoke(x, "unique", **kw)

    # ---- results ------------------------------------------------------------------------------------------------
    def project(self, value: Any) -> dict:
        if isinstance(value, (IRanges, GRanges)):
            return _project_ranges(value)
        if isinstance(value, tuple) and len(value) == 2 and isinstance(value[0], (IRanges, GRanges)):
            return _project_ranges(value[0], revmap=value[1])  # reduce(with_revmap=True)
        if isinstance(value, Seqinfo):
            return _project_seqinfo(value)
        if isinstance(value, Hits):
            return p_hits(value.query_hits + 1, value.subject_hits + 1, value.query_length, value.subject_length)
        if isinstance(value, Rle):
            return p_rle(value.values, value.lengths)
        if isinstance(value, dict) and all(isinstance(r, Rle) for r in value.values()):
            return p_rlelist(list(value), [(r.values, r.lengths) for r in value.values()])
        if isinstance(value, pl.DataFrame):
            return p_dataframe(value.height, _mcols_of(value))
        if isinstance(value, pl.Series):
            return _column(value)
        if isinstance(value, np.ndarray):
            return p_vector_from(value)
        raise TypeError(f"cannot project {type(value).__name__}")
