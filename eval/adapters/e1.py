"""Adapter for port e1 (ranges_e1, experiment e1-differential).

Port docs used: experiments/e1-differential/README.md and the public docstrings in src/ranges_e1
(IRanges, GRanges, Seqinfo, Hits, Rle).

Port index base: 0-based for *positional* results (Hits, order, nearest/precede/follow,
select=, revmap); coordinates (start/end/width) are 1-based like R. Non-positional vectors are
plain arrays; NA is a masked-array entry or None.

Translations (representation only, see e1.md):
  * R `a.b` argument -> port `a_b`; R generic -> snake_case method; `isCircular` -> `is_circular`.
  * length-1 non-NA R vector -> Python scalar; other vectors -> ndarray (masked where NA).
  * index-valued results (Hits, integer results of findOverlaps/nearest/precede/follow, order,
    revmap) are shifted +1; positive R subscripts and `new2old`-style indices are shifted -1.
  * mcols: list of (name, column) <-> polars.DataFrame; Factor <-> pl.Enum; IntegerList <-> pl.List(Int64).
  * named R vectors (seqlengths/isCircular/genome) <-> the port's dicts keyed by seqlevel.
  * coverage: Rle -> Rle projection, dict[str, Rle] -> RleList projection.
"""

from __future__ import annotations

import inspect

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

from ranges_e1 import GRanges, Hits, IRanges, Rle, Seqinfo

_RANGES = (IRanges, GRanges)
_POLARS_TYPE = {"int": pl.Int64, "dbl": pl.Float64, "lgl": pl.Boolean, "chr": pl.String}


# ---------------------------------------------------------------- R value -> port value

def _array(v: RVector) -> np.ndarray:
    """ndarray; masked when the vector holds NA (the port's NA representation)."""
    return v.np_masked() if any(x is None for x in v.values) else v.np()


def _arg(v):
    """RVector -> Python scalar (length 1, non-NA) or ndarray; everything else unchanged."""
    if isinstance(v, RVector):
        if len(v) == 1 and v.values[0] is not None:
            return v.scalar()
        return _array(v)
    return v


def _vec(v):
    """Like `_arg`, but a length-1 vector stays an array (names, where a str would be iterated)."""
    return _array(v) if isinstance(v, RVector) else v


def _named_or_plain(v):
    """Named R vector -> dict keyed by seqlevel (the port's convention); otherwise `_arg`."""
    if isinstance(v, RVector) and v.names is not None:
        if len(set(v.names)) != len(v.names):  # a dict would collapse the entries
            raise Unsupported("named vector with duplicate names")
        return dict(zip(v.names, v.list()))
    return _arg(v)


def _per_seqlevel(x, v):
    """GRanges `gaps`/`coverage` take per-seqlevel bounds as dicts keyed by seqlevel; keep R names."""
    return _named_or_plain(v) if isinstance(x, GRanges) and v is not MISSING else v


def _series(name: str, v) -> pl.Series:
    if isinstance(v, pl.Series):
        return v.alias(name)
    if isinstance(v, RVector):
        return pl.Series(name, v.list(), dtype=_POLARS_TYPE[v.type])
    if isinstance(v, np.ndarray):  # a nested call's vector: keep its dtype, masked entries -> null
        data = np.ma.getdata(v)
        dtype = pl.String if data.dtype.kind in "OUS" else pl.Series(data).dtype
        mask = np.ma.getmaskarray(v)
        return pl.Series(name, [None if m else x for x, m in zip(data.tolist(), mask.tolist())], dtype=dtype)
    raise Unsupported(f"mcols column of {type(v).__name__}")


def _frame(columns: list[tuple[str, object]]) -> pl.DataFrame | None:
    """(name, column) pairs -> polars frame; None when there are no columns."""
    if not columns:
        return None
    names = [n for n, _ in columns]
    if len(set(names)) != len(names):
        raise Unsupported("mcols: duplicate column names (polars frames cannot hold them)")
    return pl.DataFrame([_series(n, c) for n, c in columns])


def _column_from_projection(p: dict) -> pl.Series:
    kind = p["kind"]
    if kind == "vector":
        return pl.Series("", RVector(p["type"], p["values"]).list(), dtype=_POLARS_TYPE[p["type"]])
    if kind == "Factor":
        levels = p["levels"]
        vals = [None if c is None else levels[c - 1] for c in p["codes"]]
        return pl.Series("", vals, dtype=pl.Enum(levels))
    if kind == "IntegerList":
        if p["names"] is not None:
            raise Unsupported("input:IntegerList names")
        return pl.Series("", p["elements"], dtype=pl.List(pl.Int64))
    raise Unsupported(f"mcols column kind {kind}")


def _mcols_from_projection(pairs: list) -> pl.DataFrame | None:
    return _frame([(n, _column_from_projection(c)) for n, c in pairs])


_KIND_OF_DTYPE = {"b": "lgl", "i": "int", "u": "int", "f": "dbl", "O": "chr", "U": "chr", "S": "chr"}


def _subscript(i, n: int):
    """R `[` subscript (RVector or a nested call's ndarray) -> the port's 0-based Python-style subscript.

    The port's subscripts are Python-style (negatives wrap, 0 is an index), so only forms with
    the same meaning in R are passed: positive integers (shifted), full-length logicals, and
    names (the port matches them against `names`; invalid names are its own error).
    """
    if isinstance(i, RVector):
        kind, arr = i.type, _array(i)
    elif isinstance(i, np.ndarray):
        kind, arr = _KIND_OF_DTYPE.get(i.dtype.kind), i
    else:
        raise Unsupported("subscript:non-vector")
    if kind is None:
        raise Unsupported("subscript:dtype")
    if kind == "chr":
        return arr
    if np.ma.getmaskarray(arr).any():
        raise Unsupported("subscript:logical NA" if kind == "lgl" else "subscript:NA")
    arr = np.ma.getdata(arr)
    if kind == "lgl":
        if len(arr) != n:
            raise Unsupported("subscript:logical recycling")
        return arr
    if (arr != arr).any() or (arr < 1).any():
        raise Unsupported("subscript:zero/negative")
    return arr - 1


def _as_r_index(res):
    """Port 0-based positional result (Hits aside, projected separately) -> R 1-based.

    A boolean masked array is NA-only (not an index), so it is left alone.
    """
    if isinstance(res, np.ndarray) and res.dtype.kind != "b":
        return res + 1
    return res


def _revmap_to_r_base(res):
    """A `revmap` IntegerList mcol on reduce/disjoin output (port: 0-based positions) -> +1.

    These ops drop the input's mcols, so any `revmap` column is the port's own, whatever form
    the `with.revmap` flag arrived in. Two public port calls (mcols, set_mcols).
    """
    mc = res.mcols
    if mc is None or "revmap" not in mc.columns:
        return res
    shifted = mc.with_columns(pl.col("revmap").list.eval(pl.element() + 1))
    return res.set_mcols(shifted)


def _port_object(x, op: str):
    if not isinstance(x, (IRanges, GRanges, Seqinfo)):
        raise Unsupported(f"{op} on {type(x).__name__}")
    return x


def _attr(x, name: str, op: str):
    _port_object(x, op)
    if not hasattr(x, name):
        raise Unsupported(f"{type(x).__name__} has no {name}")
    return getattr(x, name)


def _send(x, method: str, *args, **kw):
    """x.method(*args, **kw) with R values converted and MISSING arguments left out."""
    fn = _attr(x, method, method)
    given = {k: v for k, v in kw.items() if v is not MISSING}
    params = inspect.signature(fn).parameters
    for k in given:  # an R argument this class's method lacks (e.g. IRanges has no ignore_strand)
        if k not in params:
            raise UnknownArgument(k)
    return fn(*[_arg(a) for a in args], **{k: _arg(v) for k, v in given.items()})


def _named_vector(kind: str, x, attr: str) -> RVector:
    """Per-seqlevel field -> named R vector: GRanges gives a dict keyed by seqlevel, Seqinfo a tuple."""
    values = _attr(x, attr, attr)
    if isinstance(values, dict):
        return RVector(kind, list(values.values()), list(values.keys()))
    return RVector(kind, list(values), list(_attr(x, "seqnames", attr)))


# ---------------------------------------------------------------- projection helpers

def _project_column(s: pl.Series) -> dict:
    dt = s.dtype
    if isinstance(dt, pl.Enum):
        levels = dt.categories.to_list()
        return p_factor(levels, s.to_physical().to_list(), base=0)
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
    raise ValueError(f"cannot project polars column of dtype {dt}")


def _mcols_pairs(df: pl.DataFrame | None) -> list:
    return [] if df is None else [(c, _project_column(df[c])) for c in df.columns]


# ---------------------------------------------------------------- adapter

class Adapter(AdapterBase):
    port = "e1"
    port_packages = ("ranges_e1",)
    unsupported_ops = {"seqlevels<-": "README: seqlevels<- (renaming/pruning) not implemented"}

    # ---- inputs ------------------------------------------------------------
    def build_IRanges(self, p: dict):
        kw = {}
        if p["names"] is not None:
            kw["names"] = p["names"]
        mc = _mcols_from_projection(p["mcols"])
        if mc is not None:
            kw["mcols"] = mc
        return IRanges(start=p["start"], width=p["width"], **kw)

    def build_GRanges(self, p: dict):
        kw = {}
        if p["names"] is not None:
            kw["names"] = p["names"]
        mc = _mcols_from_projection(p["mcols"])
        if mc is not None:
            kw["mcols"] = mc
        return GRanges(
            seqnames=p["seqnames"],
            ranges=IRanges(start=p["start"], width=p["width"]),
            strand=p["strand"],
            seqinfo=self.build_Seqinfo(p["seqinfo"]),
            **kw,
        )

    def build_Seqinfo(self, p: dict):
        return Seqinfo(p["seqnames"], p["seqlengths"], p["isCircular"], p["genome"])

    def build_DataFrame(self, p: dict):
        if p["rownames"] is not None:
            raise Unsupported("input:DataFrame rownames (port has no mcols row names)")
        df = _mcols_from_projection(p["columns"])
        if df is None:
            if p["nrow"]:
                raise Unsupported("input:DataFrame zero columns with nrow > 0")
            return pl.DataFrame()
        return df

    def build_Factor(self, p: dict):
        return _column_from_projection(p)

    def build_IntegerList(self, p: dict):
        return _column_from_projection(p)

    # ---- construct ---------------------------------------------------------
    def op_IRanges(self, *, start=MISSING, end=MISSING, width=MISSING, names=MISSING, _mcols=()):
        kw = {"start": start, "end": end, "width": width}
        if names is not MISSING:
            kw["names"] = _vec(names)
        mc = _frame(list(_mcols))
        if mc is not None:
            kw["mcols"] = mc
        return IRanges(**{k: _arg(v) for k, v in kw.items() if v is not MISSING})

    def op_GRanges(self, *, seqnames=MISSING, ranges=MISSING, strand=MISSING, seqinfo=MISSING,
                   seqlengths=MISSING, _mcols=()):
        kw = {"seqnames": seqnames, "ranges": ranges, "strand": strand, "seqinfo": seqinfo,
              "seqlengths": seqlengths}
        mc = _frame(list(_mcols))
        if mc is not None:
            kw["mcols"] = mc
        return GRanges(**{k: _arg(v) for k, v in kw.items() if v is not MISSING})

    def op_Seqinfo(self, *, seqnames=MISSING, seqlengths=MISSING, isCircular=MISSING, genome=MISSING):
        kw = {"seqnames": seqnames, "seqlengths": seqlengths, "is_circular": isCircular, "genome": genome}
        return Seqinfo(**{k: _arg(v) for k, v in kw.items() if v is not MISSING})

    # ---- accessors ---------------------------------------------------------
    def op_start(self, *, x):
        return _attr(x, "start", "start")

    def op_end(self, *, x):
        return _attr(x, "end", "end")

    def op_width(self, *, x):
        return _attr(x, "width", "width")

    def op_names(self, *, x):
        names = _attr(x, "names", "names")
        return None if names is None else RVector("chr", list(names))

    def op_length(self, *, x):
        if not isinstance(x, (IRanges, GRanges, Seqinfo, Hits, Rle, RVector)):
            raise Unsupported(f"length on {type(x).__name__}")
        return RVector("int", [len(x)])

    def op_seqnames(self, *, x):
        names = _attr(x, "seqnames", "seqnames")
        return RVector("chr", list(names)) if isinstance(names, tuple) else names  # Seqinfo: its seqlevels

    def op_strand(self, *, x):
        return _attr(x, "strand", "strand")

    def op_ranges(self, *, x):
        return _attr(x, "ranges", "ranges")

    def op_mcols(self, *, x):
        return _attr(x, "mcols", "mcols")

    # ---- metadata ----------------------------------------------------------
    def op_set_names(self, *, x, value):
        return _send(x, "set_names", _vec(value) if value is not None else None)

    def op_set_mcols(self, *, x, value):
        return _send(x, "set_mcols", value)

    # ---- seqinfo -----------------------------------------------------------
    def op_seqinfo(self, *, x):
        return _attr(x, "seqinfo", "seqinfo")

    def op_seqlevels(self, *, x):
        attr = "seqnames" if isinstance(x, Seqinfo) else "seqlevels"
        return RVector("chr", list(_attr(x, attr, "seqlevels")))

    def op_seqlengths(self, *, x):
        return _named_vector("int", x, "seqlengths")

    def op_isCircular(self, *, x):
        return _named_vector("lgl", x, "is_circular")

    def op_genome(self, *, x):
        return _named_vector("chr", x, "genome")

    def op_set_seqinfo(self, *, x, value):
        return _send(x, "set_seqinfo", value)

    def op_set_seqlengths(self, *, x, value):
        return _attr(x, "set_seqlengths", "seqlengths<-")(_named_or_plain(value))

    def op_set_isCircular(self, *, x, value):
        return _attr(x, "set_is_circular", "isCircular<-")(_named_or_plain(value))

    def op_set_genome(self, *, x, value):
        return _attr(x, "set_genome", "genome<-")(_named_or_plain(value))

    def op_merge(self, *dots, x, y=MISSING):
        if not isinstance(x, Seqinfo):
            raise Unsupported(f"merge on {type(x).__name__}")
        if y is MISSING:
            raise Unsupported("merge without y")
        return Seqinfo.merge_all([x, y, *dots])

    def op_trim(self, *, x, use_names=MISSING):
        return _send(x, "trim", use_names=use_names)

    # ---- subset / concat ---------------------------------------------------
    def op_subset(self, *, x, i=MISSING):
        if i is MISSING:
            raise Unsupported("subscript:missing")
        if not isinstance(x, _RANGES):
            raise Unsupported(f"[ on {type(x).__name__}")
        return x[_subscript(i, len(x))]

    def op_c(self, *dots, x):
        concat = getattr(type(x), "concat", None)
        if concat is None or not isinstance(x, _RANGES):
            raise Unsupported(f"c on {type(x).__name__}")
        return concat([x, *dots])

    # ---- intra-range -------------------------------------------------------
    def op_shift(self, *, x, shift=MISSING, use_names=MISSING):
        return _send(x, "shift", shift=shift, use_names=use_names)

    def op_narrow(self, *, x, start=MISSING, end=MISSING, width=MISSING, use_names=MISSING):
        return _send(x, "narrow", start=start, end=end, width=width, use_names=use_names)

    def op_resize(self, *, x, width=MISSING, fix=MISSING, use_names=MISSING, ignore_strand=MISSING):
        return _send(x, "resize", width=width, fix=fix, use_names=use_names, ignore_strand=ignore_strand)

    def op_flank(self, *, x, width=MISSING, start=MISSING, both=MISSING, use_names=MISSING,
                 ignore_strand=MISSING):
        return _send(x, "flank", width=width, start=start, both=both, use_names=use_names,
                     ignore_strand=ignore_strand)

    def op_promoters(self, *, x, upstream=MISSING, downstream=MISSING, use_names=MISSING):
        return _send(x, "promoters", upstream=upstream, downstream=downstream, use_names=use_names)

    # ---- inter-range -------------------------------------------------------
    def op_range(self, *dots, x, ignore_strand=MISSING):
        return _send(x, "range", *dots, ignore_strand=ignore_strand)

    def op_reduce(self, *, x, drop_empty_ranges=MISSING, min_gapwidth=MISSING, with_revmap=MISSING,
                  ignore_strand=MISSING):
        res = _send(x, "reduce", drop_empty_ranges=drop_empty_ranges, min_gapwidth=min_gapwidth,
                    with_revmap=with_revmap, ignore_strand=ignore_strand)
        return _revmap_to_r_base(res)

    def op_gaps(self, *, x, start=MISSING, end=MISSING, ignore_strand=MISSING):
        return _send(x, "gaps", start=_per_seqlevel(x, start), end=_per_seqlevel(x, end),
                     ignore_strand=ignore_strand)

    def op_disjoin(self, *, x, with_revmap=MISSING, ignore_strand=MISSING):
        res = _send(x, "disjoin", with_revmap=with_revmap, ignore_strand=ignore_strand)
        return _revmap_to_r_base(res)

    # ---- set operations ----------------------------------------------------
    def op_union(self, *, x, y=MISSING, ignore_strand=MISSING):
        return self._pair(x, "union", y, ignore_strand=ignore_strand)

    def op_intersect(self, *, x, y=MISSING, ignore_strand=MISSING):
        return self._pair(x, "intersect", y, ignore_strand=ignore_strand)

    def op_setdiff(self, *, x, y=MISSING, ignore_strand=MISSING):
        return self._pair(x, "setdiff", y, ignore_strand=ignore_strand)

    def op_punion(self, *, x, y=MISSING, fill_gap=MISSING, ignore_strand=MISSING):
        return self._pair(x, "punion", y, fill_gap=fill_gap, ignore_strand=ignore_strand)

    def op_pintersect(self, *, x, y=MISSING, resolve_empty=MISSING, drop_nohit_ranges=MISSING,
                      ignore_strand=MISSING, strict_strand=MISSING):
        return self._pair(x, "pintersect", y, resolve_empty=resolve_empty,
                          drop_nohit_ranges=drop_nohit_ranges, ignore_strand=ignore_strand,
                          strict_strand=strict_strand)

    def op_psetdiff(self, *, x, y=MISSING, ignore_strand=MISSING):
        return self._pair(x, "psetdiff", y, ignore_strand=ignore_strand)

    def _pair(self, x, method: str, y, **kw):
        if y is MISSING:
            raise Unsupported(f"{method} without y")
        return _send(x, method, y, **kw)

    # ---- overlaps ----------------------------------------------------------
    def op_findOverlaps(self, *, query, subject=MISSING, maxgap=MISSING, minoverlap=MISSING, type=MISSING,
                        select=MISSING, ignore_strand=MISSING, drop_self=MISSING, drop_redundant=MISSING):
        res = _send(query, "find_overlaps", subject=subject, maxgap=maxgap, minoverlap=minoverlap,
                    type=type, select=select, ignore_strand=ignore_strand, drop_self=drop_self,
                    drop_redundant=drop_redundant)
        return _as_r_index(res)

    def op_countOverlaps(self, *, query, subject=MISSING, maxgap=MISSING, minoverlap=MISSING, type=MISSING,
                         ignore_strand=MISSING):
        return _send(query, "count_overlaps", subject=subject, maxgap=maxgap, minoverlap=minoverlap,
                     type=type, ignore_strand=ignore_strand)

    def op_overlapsAny(self, *, query, subject=MISSING, maxgap=MISSING, minoverlap=MISSING, type=MISSING,
                       ignore_strand=MISSING):
        return _send(query, "overlaps_any", subject=subject, maxgap=maxgap, minoverlap=minoverlap,
                     type=type, ignore_strand=ignore_strand)

    def op_subsetByOverlaps(self, *, x, ranges=MISSING, maxgap=MISSING, minoverlap=MISSING, type=MISSING,
                            invert=MISSING, ignore_strand=MISSING):
        return _send(x, "subset_by_overlaps", subject=ranges, maxgap=maxgap, minoverlap=minoverlap,
                     type=type, invert=invert, ignore_strand=ignore_strand)

    # ---- nearest -----------------------------------------------------------
    def op_nearest(self, *, x, subject=MISSING, select=MISSING, ignore_strand=MISSING):
        return _as_r_index(_send(x, "nearest", subject=subject, select=select, ignore_strand=ignore_strand))

    def op_precede(self, *, x, subject=MISSING, select=MISSING, ignore_strand=MISSING):
        return _as_r_index(_send(x, "precede", subject=subject, select=select, ignore_strand=ignore_strand))

    def op_follow(self, *, x, subject=MISSING, select=MISSING, ignore_strand=MISSING):
        return _as_r_index(_send(x, "follow", subject=subject, select=select, ignore_strand=ignore_strand))

    def op_distance(self, *, x, y=MISSING, ignore_strand=MISSING):
        return self._pair(x, "distance", y, ignore_strand=ignore_strand)

    # ---- coverage ----------------------------------------------------------
    def op_coverage(self, *, x, shift=MISSING, width=MISSING, weight=MISSING):
        return _send(x, "coverage", shift=shift, width=_per_seqlevel(x, width), weight=weight)

    # ---- ordering ----------------------------------------------------------
    def op_order(self, *dots, decreasing=MISSING):
        if len(dots) != 1:
            raise Unsupported("order: needs exactly one key (multi-key order not implemented)")
        return _as_r_index(_send(dots[0], "order", decreasing=decreasing))

    def op_sort(self, *, x, decreasing=MISSING, ignore_strand=MISSING):
        return _send(x, "sort", decreasing=decreasing, ignore_strand=ignore_strand)

    def op_rank(self, *, x, ties_method=MISSING, ignore_strand=MISSING):
        return _send(x, "rank", ties_method=ties_method, ignore_strand=ignore_strand)

    def op_duplicated(self, *, x, fromLast=MISSING):
        return _send(x, "duplicated", from_last=fromLast)

    def op_unique(self, *, x, fromLast=MISSING):
        return _send(x, "unique", from_last=fromLast)

    # ---- projection --------------------------------------------------------
    def project(self, value) -> dict:
        if isinstance(value, GRanges):
            return p_granges(value.seqnames, value.start, value.width, value.strand, value.names,
                             _mcols_pairs(value.mcols), self.project(value.seqinfo))
        if isinstance(value, IRanges):
            return p_iranges(value.start, value.width, value.names, _mcols_pairs(value.mcols))
        if isinstance(value, Seqinfo):
            return p_seqinfo(value.seqnames, value.seqlengths, value.is_circular, value.genome)
        if isinstance(value, Hits):
            return p_hits(value.query_hits + 1, value.subject_hits + 1, value.query_length, value.subject_length)
        if isinstance(value, Rle):
            return p_rle(value.values, value.lengths)
        if isinstance(value, dict):  # GRanges coverage: seqlevel -> Rle
            return p_rlelist(list(value), [(r.values, r.lengths) for r in value.values()])
        if isinstance(value, pl.DataFrame):
            return p_dataframe(value.height, _mcols_pairs(value))
        if isinstance(value, pl.Series):
            return _project_column(value)
        if isinstance(value, np.ndarray):
            if value.dtype == object:  # seqnames / strand arrays hold str
                return p_vector("chr", value)
            return p_vector_from(value)
        raise ValueError(f"cannot project {type(value).__name__}")
