"""Eval adapter for port e4 (ranges_e4).

Port docs used: experiments/e4-primitive-cases README.md ("Package design", "Known limits"),
docstrings and public names in src/ranges_e4 (iranges, granges, seqinfo, hits, rle).
Index base: 0-based (README "Index base"). Translated here: Hits (both columns), integer
results of select != "all" (masked arrays), `order`, `[` integer subscripts, and the
`revmap` mcols column of range/reduce/disjoin results.

Translations (representation only; see e4.md):
  * R `a.b` argument names -> port `a_b`; exceptions: fromLast -> from_last,
    subsetByOverlaps `ranges` -> `subject`, merge `y` -> `other`.
  * RVector -> ndarray (masked array when it holds NA); length-1 scalars -> Python scalars
    where the port takes a scalar; chr `weight` (mcols column name) -> str.
  * Port-side mcols are polars frames: vector -> Int64/Float64/String/Boolean series,
    Factor -> pl.Enum series (codes <-> levels), IntegerList -> List(Int64) series.
  * Named R vectors the port returns as dicts (seqlengths, isCircular, genome) -> RVector.
  * Rle / dict[str, Rle] -> Rle / RleList projections.
  * The port exposes no setter for seqinfo; `names<-`/`mcols<-` use with_names/with_mcols.
"""

from __future__ import annotations

import inspect
from typing import Any

import numpy as np
import polars as pl
from ranges_e4 import GRanges, Hits, IRanges, Rle, Seqinfo

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

_RANGES = (IRanges, GRanges)
_PORT_OBJECTS = (IRanges, GRanges, Seqinfo)
_SERIES_DTYPES = {"int": pl.Int64, "dbl": pl.Float64, "chr": pl.String, "lgl": pl.Boolean}


# ---------------------------------------------------------------- argument conversion

def _np(v: Any) -> Any:
    """RVector -> ndarray (masked where NA); anything else passes through."""
    if not isinstance(v, RVector):
        return v
    return v.np_masked() if any(x is None for x in v.values) else v.np()


def _scalar(name: str, v: Any) -> Any:
    """RVector of length 1 -> Python scalar (None = NA); R NULL stays None."""
    if v is None:
        return None
    if isinstance(v, RVector):
        if len(v) != 1:
            raise Unsupported(f"argument:{name}: port takes a scalar, case passes length {len(v)}")
        return v.scalar()
    return v


def _strs(v: Any) -> Any:
    """chr RVector -> list (None = NA); anything else passes through."""
    return v.list() if isinstance(v, RVector) else v


def _ints(values: list) -> np.ndarray:
    """Projection int list -> int64 array, or a masked array when it holds NA."""
    if any(x is None for x in values):
        return np.ma.masked_array([0 if x is None else x for x in values], mask=[x is None for x in values], dtype=np.int64)
    return np.array(values, dtype=np.int64)


def _series(name: str, col: Any) -> pl.Series:
    """One mcols/DataFrame column (projection dict or RVector) -> polars series."""
    if isinstance(col, pl.Series):
        return col.alias(name)
    if isinstance(col, RVector):
        return pl.Series(name, col.list(), dtype=_SERIES_DTYPES[col.type])
    if col["kind"] == "vector":
        return _series(name, RVector(col["type"], col["values"]))
    if col["kind"] == "Factor":
        levels = col["levels"]
        vals = [None if c is None else levels[c - 1] for c in col["codes"]]
        return pl.Series(name, vals, dtype=pl.Enum(levels))
    if col["kind"] == "IntegerList":
        return pl.Series(name, col["elements"], dtype=pl.List(pl.Int64))
    raise Unsupported(f"mcols column of kind {col['kind']}")


def _frame(pairs: Any) -> pl.DataFrame | None:
    """[(name, column)] -> polars frame; no columns -> None (port: no mcols)."""
    series = [_series(str(n), c) for n, c in pairs]
    return pl.DataFrame(series) if series else None


def _mcols_arg(v: Any) -> Any:
    if v is None or isinstance(v, pl.DataFrame):
        return v
    raise Unsupported(f"mcols value of type {type(v).__name__}")


# ---------------------------------------------------------------- result projection helpers

def _column_projection(s: pl.Series) -> dict:
    dt = s.dtype
    if isinstance(dt, pl.Enum):
        return p_factor(dt.categories.to_list(), s.to_physical().to_list(), base=0)
    if isinstance(dt, pl.List):
        return p_integerlist(s.to_list())
    if dt.is_integer():
        return p_vector("int", s.to_list())
    if dt.is_float():
        return p_vector("dbl", s.to_list())
    if dt == pl.Boolean:
        return p_vector("lgl", s.to_list())
    if dt == pl.String:
        return p_vector("chr", s.to_list())
    raise ValueError(f"cannot project mcols column of dtype {dt}")


def _mcols_projection(df: pl.DataFrame) -> list:
    return [(c, _column_projection(df[c])) for c in df.columns]


def _one_based(res: Any) -> Any:
    """0-based index results (select != "all", order) -> 1-based; Hits are shifted in project()."""
    return res if isinstance(res, Hits) else res + 1


class Adapter(AdapterBase):
    port = "e4"
    port_packages = ("ranges_e4",)
    unsupported_ops = {
        "seqinfo<-": "port has no public seqinfo setter (objects are immutable)",
        "seqlevels<-": "port has no public seqlevels setter",
        "seqlengths<-": "port has no public seqlengths setter",
        "isCircular<-": "port has no public isCircular setter",
        "genome<-": "port has no public genome setter",
    }

    # ------------------------------------------------------------ generic call into the port

    @staticmethod
    def _need(op: str, x: Any, classes: tuple = _RANGES) -> None:
        if not isinstance(x, classes):
            raise Unsupported(f"{op}: input is {type(x).__name__}, not a port object")

    def _invoke(self, op: str, x: Any, method: str, kw: dict, args: tuple = (), classes: tuple = _RANGES) -> Any:
        """x.method(*args, **kw), dropping MISSING.

        An argument the port method lacks is UnknownArgument; a port parameter without a
        default that the case did not supply is Unsupported (the port has no default).
        """
        self._need(op, x, classes)
        if not hasattr(x, method):
            raise Unsupported(f"{op}: {type(x).__name__} has no {method}")
        fn = getattr(x, method)
        kw = {k: v for k, v in kw.items() if v is not MISSING}
        params = inspect.signature(fn).parameters
        for k in kw:
            if k not in params:
                raise UnknownArgument(k)
        if args and not any(p.kind is p.VAR_POSITIONAL for p in params.values()):
            raise Unsupported(f"{op}: extra positional objects")
        for name, p in params.items():
            if p.default is p.empty and p.kind in (p.POSITIONAL_OR_KEYWORD, p.KEYWORD_ONLY) and name not in kw:
                raise Unsupported(f"argument:{name}: no default in port")
        return fn(*args, **kw)

    def _property(self, op: str, x: Any, attr: str, classes: tuple = _RANGES) -> Any:
        self._need(op, x, classes)
        if not hasattr(x, attr):
            raise Unsupported(f"{op}: {type(x).__name__} has no {attr}")
        return getattr(x, attr)

    @staticmethod
    def _revmap_one_based(res: Any) -> Any:
        """The port stores `revmap` 0-based; shift it so the native result is R-indexed."""
        if "revmap" not in res.mcols.columns:
            return res
        col = res.mcols["revmap"]
        shifted = pl.Series("revmap", [[v + 1 for v in m] for m in col.to_list()], dtype=pl.List(pl.Int64))
        return res.with_mcols(res.mcols.with_columns(shifted))

    # ------------------------------------------------------------ inputs

    def build_IRanges(self, p: dict) -> Any:
        return IRanges(start=_ints(p["start"]), width=_ints(p["width"]), names=p["names"], mcols=_frame(p["mcols"]))

    def build_Seqinfo(self, p: dict) -> Any:
        return Seqinfo(p["seqnames"], p["seqlengths"], p["isCircular"], p["genome"])

    def build_GRanges(self, p: dict) -> Any:
        ranges = IRanges(start=_ints(p["start"]), width=_ints(p["width"]), names=p["names"])
        return GRanges(
            p["seqnames"], ranges, p["strand"], self.build_Seqinfo(p["seqinfo"]), mcols=_frame(p["mcols"])
        )

    def build_DataFrame(self, p: dict) -> Any:
        if p["rownames"] is not None:
            raise Unsupported("input:DataFrame rownames: polars frames carry no row names")
        df = pl.DataFrame([_series(n, c) for n, c in p["columns"]])
        if not p["columns"] and p["nrow"]:
            # a zero-column frame of the requested height: add a column, then drop it
            df = pl.DataFrame({"_": range(p["nrow"])}).drop("_")
        return df

    def build_Factor(self, p: dict) -> Any:
        return _series("", {"kind": "Factor", "levels": p["levels"], "codes": p["codes"]})

    # ------------------------------------------------------------ constructors

    def op_IRanges(self, *, start=MISSING, end=MISSING, width=MISSING, names=MISSING, _mcols=()):
        kw = {"start": _np(start), "end": _np(end), "width": _np(width), "names": _strs(names)}
        kw = {k: v for k, v in kw.items() if v is not MISSING}
        return IRanges(**kw, mcols=_frame(_mcols))

    def op_GRanges(self, *, seqnames=MISSING, ranges=MISSING, strand=MISSING, seqinfo=MISSING, _mcols=()):
        kw = {"seqnames": _strs(seqnames), "ranges": ranges, "strand": _strs(strand), "seqinfo": seqinfo}
        kw = {k: v for k, v in kw.items() if v is not MISSING}
        if "seqnames" not in kw:
            raise Unsupported("argument:seqnames: no default in port")
        return GRanges(**kw, mcols=_frame(_mcols))

    def op_Seqinfo(self, *, seqnames=MISSING, seqlengths=MISSING, isCircular=MISSING, genome=MISSING):
        if seqnames is MISSING:
            raise Unsupported("argument:seqnames: no default in port")
        kw = {"seqlengths": _strs(seqlengths), "is_circular": _strs(isCircular), "genome": _strs(genome)}
        return Seqinfo(_strs(seqnames), **{k: v for k, v in kw.items() if v is not MISSING})

    # ------------------------------------------------------------ accessors

    def op_start(self, *, x):
        return self._property("start", x, "start")

    def op_end(self, *, x):
        return self._property("end", x, "end")

    def op_width(self, *, x):
        return self._property("width", x, "width")

    def op_names(self, *, x):
        names = self._property("names", x, "names")
        return None if names is None else RVector("chr", names.tolist())

    def op_length(self, *, x):
        self._need("length", x, _PORT_OBJECTS + (Hits, Rle))
        return RVector("int", [len(x)])

    def op_seqnames(self, *, x):
        return RVector("chr", list(self._property("seqnames", x, "seqnames", _PORT_OBJECTS)))

    def op_strand(self, *, x):
        return RVector("chr", list(self._property("strand", x, "strand")))

    def op_ranges(self, *, x):
        return self._property("ranges", x, "ranges")

    def op_mcols(self, *, x):
        return self._property("mcols", x, "mcols")

    # ------------------------------------------------------------ metadata setters

    def op_set_names(self, *, x, value=MISSING):
        if value is MISSING:
            raise Unsupported("argument:value: missing")
        return self._invoke("names<-", x, "with_names", {"names": _strs(value)})

    def op_set_mcols(self, *, x, value=MISSING):
        if value is MISSING:
            raise Unsupported("argument:value: missing")
        return self._invoke("mcols<-", x, "with_mcols", {"mcols": _mcols_arg(value)})

    # ------------------------------------------------------------ seqinfo family

    def op_seqinfo(self, *, x):
        return self._property("seqinfo", x, "seqinfo")

    def op_seqlevels(self, *, x):
        return RVector("chr", list(self._property("seqlevels", x, "seqlevels", _PORT_OBJECTS)))

    def op_seqlengths(self, *, x):
        d = self._property("seqlengths", x, "seqlengths", _PORT_OBJECTS)
        return RVector("int", list(d.values()), list(d.keys()))

    def op_isCircular(self, *, x):
        d = self._property("isCircular", x, "is_circular", _PORT_OBJECTS)
        return RVector("lgl", list(d.values()), list(d.keys()))

    def op_genome(self, *, x):
        d = self._property("genome", x, "genome", _PORT_OBJECTS)
        return RVector("chr", list(d.values()), list(d.keys()))

    def op_merge(self, *, x, y=MISSING):
        return self._invoke("merge", x, "merge", {"other": y}, classes=(Seqinfo,))

    def op_trim(self, *, x):
        return self._invoke("trim", x, "trim", {})

    # ------------------------------------------------------------ subset / concat

    def op_subset(self, *, x, i=MISSING):
        self._need("[", x)
        if i is MISSING:
            raise Unsupported("subscript:missing")
        if not isinstance(i, RVector):
            raise Unsupported("subscript:NULL" if i is None else f"subscript:{type(i).__name__}")
        if any(v is None for v in i.values):
            raise Unsupported("subscript:NA")
        if i.type == "int":
            if any(v <= 0 for v in i.values):
                raise Unsupported("subscript:non-positive integer")
            idx = i.np() - 1  # port is 0-based
        elif i.type == "lgl":
            if len(i) != len(x):
                raise Unsupported("subscript:logical recycling")
            idx = i.np()
        elif i.type == "chr":
            idx = i.np()
        else:
            raise Unsupported(f"subscript:{i.type}")
        return x[idx]

    def op_c(self, *dots, x):
        self._need("c", x)
        return type(x).concat(x, *dots)

    # ------------------------------------------------------------ intra-range

    def op_shift(self, *, x, shift=MISSING):
        return self._invoke("shift", x, "shift", {"shift": _np(shift)})

    def op_narrow(self, *, x, start=MISSING, end=MISSING, width=MISSING):
        return self._invoke("narrow", x, "narrow", {"start": _np(start), "end": _np(end), "width": _np(width)})

    def op_resize(self, *, x, width=MISSING, fix=MISSING, use_names=MISSING, ignore_strand=MISSING):
        kw = {"width": _np(width), "fix": _np(fix), "use_names": _scalar("use.names", use_names),
              "ignore_strand": _scalar("ignore.strand", ignore_strand)}
        return self._invoke("resize", x, "resize", kw)

    def op_flank(self, *, x, width=MISSING, start=MISSING, both=MISSING, use_names=MISSING, ignore_strand=MISSING):
        kw = {"width": _np(width), "start": _np(start), "both": _scalar("both", both),
              "use_names": _scalar("use.names", use_names), "ignore_strand": _scalar("ignore.strand", ignore_strand)}
        return self._invoke("flank", x, "flank", kw)

    def op_promoters(self, *, x, upstream=MISSING, downstream=MISSING, use_names=MISSING):
        kw = {"upstream": _np(upstream), "downstream": _np(downstream), "use_names": _scalar("use.names", use_names)}
        return self._invoke("promoters", x, "promoters", kw)

    # ------------------------------------------------------------ inter-range

    def op_range(self, *dots, x, ignore_strand=MISSING):
        res = self._invoke("range", x, "range", {"ignore_strand": _scalar("ignore.strand", ignore_strand)}, args=dots)
        return self._revmap_one_based(res)

    def op_reduce(self, *, x, drop_empty_ranges=MISSING, min_gapwidth=MISSING, with_revmap=MISSING, ignore_strand=MISSING):
        kw = {"drop_empty_ranges": _scalar("drop.empty.ranges", drop_empty_ranges),
              "min_gapwidth": _scalar("min.gapwidth", min_gapwidth),
              "with_revmap": _scalar("with.revmap", with_revmap),
              "ignore_strand": _scalar("ignore.strand", ignore_strand)}
        return self._revmap_one_based(self._invoke("reduce", x, "reduce", kw))

    def op_gaps(self, *, x, start=MISSING, end=MISSING, ignore_strand=MISSING):
        kw = {"start": _scalar("start", start), "end": _scalar("end", end),
              "ignore_strand": _scalar("ignore.strand", ignore_strand)}
        return self._invoke("gaps", x, "gaps", kw)

    def op_disjoin(self, *, x, with_revmap=MISSING, ignore_strand=MISSING):
        kw = {"with_revmap": _scalar("with.revmap", with_revmap), "ignore_strand": _scalar("ignore.strand", ignore_strand)}
        return self._revmap_one_based(self._invoke("disjoin", x, "disjoin", kw))

    # ------------------------------------------------------------ set operations

    def _setop(self, op: str, x: Any, y: Any, **kw: Any) -> Any:
        return self._invoke(op, x, op, {"y": y, **{k: _scalar(k, v) for k, v in kw.items() if v is not MISSING}})

    def op_union(self, *, x, y=MISSING, ignore_strand=MISSING):
        return self._setop("union", x, y, ignore_strand=ignore_strand)

    def op_intersect(self, *, x, y=MISSING, ignore_strand=MISSING):
        return self._setop("intersect", x, y, ignore_strand=ignore_strand)

    def op_setdiff(self, *, x, y=MISSING, ignore_strand=MISSING):
        return self._setop("setdiff", x, y, ignore_strand=ignore_strand)

    def op_punion(self, *, x, y=MISSING, fill_gap=MISSING, ignore_strand=MISSING):
        return self._setop("punion", x, y, fill_gap=fill_gap, ignore_strand=ignore_strand)

    def op_pintersect(self, *, x, y=MISSING, resolve_empty=MISSING, ignore_strand=MISSING):
        return self._setop("pintersect", x, y, resolve_empty=resolve_empty, ignore_strand=ignore_strand)

    def op_psetdiff(self, *, x, y=MISSING, ignore_strand=MISSING):
        return self._setop("psetdiff", x, y, ignore_strand=ignore_strand)

    # ------------------------------------------------------------ overlaps

    def _overlap_kw(self, subject: Any, maxgap: Any, minoverlap: Any, type: Any, **more: Any) -> dict:
        kw = {"subject": subject, "maxgap": _scalar("maxgap", maxgap), "minoverlap": _scalar("minoverlap", minoverlap),
              "type": _scalar("type", type)}
        kw.update({k: _scalar(k.replace("_", "."), v) for k, v in more.items()})
        return kw

    def op_findOverlaps(self, *, query, subject=MISSING, maxgap=MISSING, minoverlap=MISSING, type=MISSING,
                        select=MISSING, ignore_strand=MISSING):
        kw = self._overlap_kw(subject, maxgap, minoverlap, type, select=select, ignore_strand=ignore_strand)
        return _one_based(self._invoke("findOverlaps", query, "find_overlaps", kw))

    def op_countOverlaps(self, *, query, subject=MISSING, maxgap=MISSING, minoverlap=MISSING, type=MISSING,
                         ignore_strand=MISSING):
        kw = self._overlap_kw(subject, maxgap, minoverlap, type, ignore_strand=ignore_strand)
        return self._invoke("countOverlaps", query, "count_overlaps", kw)

    def op_overlapsAny(self, *, query, subject=MISSING, maxgap=MISSING, minoverlap=MISSING, type=MISSING,
                       ignore_strand=MISSING):
        kw = self._overlap_kw(subject, maxgap, minoverlap, type, ignore_strand=ignore_strand)
        return self._invoke("overlapsAny", query, "overlaps_any", kw)

    def op_subsetByOverlaps(self, *, x, ranges=MISSING, maxgap=MISSING, minoverlap=MISSING, type=MISSING,
                            invert=MISSING, ignore_strand=MISSING):
        kw = self._overlap_kw(ranges, maxgap, minoverlap, type, invert=invert, ignore_strand=ignore_strand)
        return self._invoke("subsetByOverlaps", x, "subset_by_overlaps", kw)

    # ------------------------------------------------------------ nearest

    def _nearest(self, op: str, method: str, x: Any, subject: Any, select: Any, ignore_strand: Any) -> Any:
        kw = {"subject": subject, "select": _scalar("select", select), "ignore_strand": _scalar("ignore.strand", ignore_strand)}
        return _one_based(self._invoke(op, x, method, kw))

    def op_nearest(self, *, x, subject=MISSING, select=MISSING, ignore_strand=MISSING):
        return self._nearest("nearest", "nearest", x, subject, select, ignore_strand)

    def op_precede(self, *, x, subject=MISSING, select=MISSING, ignore_strand=MISSING):
        return self._nearest("precede", "precede", x, subject, select, ignore_strand)

    def op_follow(self, *, x, subject=MISSING, select=MISSING, ignore_strand=MISSING):
        return self._nearest("follow", "follow", x, subject, select, ignore_strand)

    def op_distance(self, *, x, y=MISSING, ignore_strand=MISSING):
        return self._invoke("distance", x, "distance", {"y": y, "ignore_strand": _scalar("ignore.strand", ignore_strand)})

    # ------------------------------------------------------------ coverage

    def op_coverage(self, *, x, shift=MISSING, width=MISSING, weight=MISSING):
        w = _scalar("weight", weight) if isinstance(weight, RVector) and weight.type == "chr" else _np(weight)
        return self._invoke("coverage", x, "coverage", {"shift": _np(shift), "width": _np(width), "weight": w})

    # ------------------------------------------------------------ ordering

    def op_order(self, *dots, decreasing=MISSING):
        if len(dots) != 1:
            raise Unsupported(f"order: port orders one object, case passes {len(dots)}")
        res = self._invoke("order", dots[0], "order", {"decreasing": _scalar("decreasing", decreasing)})
        return _one_based(res)

    def op_sort(self, *, x, decreasing=MISSING):
        return self._invoke("sort", x, "sort", {"decreasing": _scalar("decreasing", decreasing)})

    def op_rank(self, *, x, ties_method=MISSING):
        return self._invoke("rank", x, "rank", {"ties_method": _scalar("ties.method", ties_method)})

    def op_duplicated(self, *, x, fromLast=MISSING):
        return self._invoke("duplicated", x, "duplicated", {"from_last": _scalar("fromLast", fromLast)})

    def op_unique(self, *, x):
        return self._invoke("unique", x, "unique", {})

    # ------------------------------------------------------------ projection

    def _project_seqinfo(self, x: Seqinfo) -> dict:
        return p_seqinfo(x.seqnames, list(x.seqlengths.values()), list(x.is_circular.values()), list(x.genome.values()))

    def project(self, value: Any) -> dict:
        if isinstance(value, GRanges):
            return p_granges(
                value.seqnames, value.start, value.width, value.strand, value.names,
                _mcols_projection(value.mcols), self._project_seqinfo(value.seqinfo),
            )
        if isinstance(value, IRanges):
            return p_iranges(value.start, value.width, value.names, _mcols_projection(value.mcols))
        if isinstance(value, Seqinfo):
            return self._project_seqinfo(value)
        if isinstance(value, Hits):
            return p_hits(value.query + 1, value.subject + 1, value.query_length, value.subject_length)
        if isinstance(value, Rle):
            return p_rle(value.values, value.lengths)
        if isinstance(value, dict) and all(isinstance(r, Rle) for r in value.values()):
            return p_rlelist(list(value), [(r.values, r.lengths) for r in value.values()])
        if isinstance(value, pl.DataFrame):
            return p_dataframe(value.height, _mcols_projection(value))
        if isinstance(value, pl.Series):
            return _column_projection(value)
        if isinstance(value, np.ndarray):
            return p_vector_from(value)
        raise ValueError(f"cannot project {type(value).__name__}")
