"""Adapter for port e2 (ranges_e2).

Port docs used: experiments/e2-conformance/README.md ("Python API"), and the docstrings / public names
of ranges_e2 (IRanges, GRanges, Seqinfo, Hits, Rle, Factor, NA_INDEX; mcols are polars DataFrames).
Index base: 0-based (Hits, select outputs, order). NA_INDEX = -1 marks "no hit" (select outputs,
GRanges distance). Coordinates stay 1-based closed. Method names are snake_case.

Translations (representation only):
- R `[`: positive integer subscripts are shifted by -1; logical masks and character names are passed as-is
  (logical only when its length equals len(x): the port documents no recycling). Negative / zero / NA
  subscripts, NULL, double subscripts and the `j` / `drop` arguments are Unsupported.
- Index-valued results (Hits, select integer results, order) are +1; NA_INDEX -> NA (select outputs, GRanges distance).
- R atomic arguments: RVector -> list (scalar where the port takes a scalar); NA stays None.
- mcols: extra named args / `mcols<-` values become polars DataFrames (int -> Int32, dbl -> Float64, chr -> String,
  lgl -> Boolean, Factor -> Enum, IntegerList -> List(Int32)). A DataFrame input / mcols result is a `_DF`
  (frame + nrow): the port stores "no columns" as None, so its frame cannot carry the row count.
- Any mcols list column named `revmap` is index-valued: built with shift -1 from inputs, projected with shift +1
  (symmetric, so user metadata named revmap round-trips; the port never interprets list-column values).
- Rle of a Factor (seqnames / strand) is projected as an expanded chr vector, as PROJECTION.md says.
- Port arguments an R method lacks (and vice versa) are Unsupported / UnknownArgument; nothing is filled in.
"""

from dataclasses import dataclass
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

from ranges_e2 import NA_INDEX, Factor, GRanges, Hits, IRanges, Rle, Seqinfo

_POLARS_DTYPE = {"int": pl.Int32, "dbl": pl.Float64, "chr": pl.String, "lgl": pl.Boolean}
_RVECTOR_TYPE = {"i": "int", "u": "int", "f": "dbl", "b": "lgl"}
_GR_ONLY = {"ignore_strand", "drop_nohit_ranges"}
_IR_ONLY = {"resolve_empty"}
_REVMAP = "revmap"


@dataclass(frozen=True)
class _DF:
    """A DataFrame input / mcols result: polars frame plus the row count it cannot hold when column-less."""

    frame: pl.DataFrame
    nrow: int


# ---------------------------------------------------------------- argument conversion

def _sc(v: Any) -> Any:
    """R vector -> Python scalar if length 1, else list (the port sees the vector R would); None = NULL; MISSING passes."""
    if not isinstance(v, RVector):
        return v
    return v.scalar() if len(v) == 1 else v.list()


def _ls(v: Any) -> Any:
    """R vector -> Python list (None = NA / NULL); MISSING and native objects pass through."""
    return v.list() if isinstance(v, RVector) else v


def _present(**kw: Any) -> dict:
    return {k: v for k, v in kw.items() if v is not MISSING}


def _call(x: Any, method: str, *args: Any, **kw: Any) -> Any:
    """x.<method>(*args, **kw) for an IRanges / GRanges; R arguments the class lacks are unsupported."""
    gr = isinstance(x, GRanges)
    if not gr and not isinstance(x, IRanges):
        raise Unsupported(f"{method} on {type(x).__name__}")
    kw = _present(**kw)
    for k in kw:
        if (not gr and k in _GR_ONLY) or (gr and k in _IR_ONLY):
            raise UnknownArgument(k)
    return getattr(x, method)(*args, **kw)


def _need_ranges(x: Any, what: str) -> None:
    if not isinstance(x, (IRanges, GRanges)):
        raise Unsupported(f"{what} on {type(x).__name__}")


def _need_granges(x: Any, what: str) -> None:
    if not isinstance(x, GRanges):
        raise Unsupported(f"{what} on {type(x).__name__}")


def _mapping(v: RVector) -> dict:
    """Named R vector -> {name: value}; duplicated names would silently drop entries, so they are unsupported."""
    if len(set(v.names)) != len(v.names):
        raise Unsupported("argument:duplicated names in a named vector")
    return dict(zip(v.names, v.list()))


def _seqinfo_of(x: Any, what: str) -> Seqinfo:
    if isinstance(x, GRanges):
        return x.seqinfo
    if isinstance(x, Seqinfo):
        return x
    raise Unsupported(f"{what} on {type(x).__name__}")


def _column(name: str, v: Any) -> pl.Series:
    if not isinstance(v, RVector):
        raise Unsupported(f"mcols value {type(v).__name__}")
    return pl.Series(name, v.list(), dtype=_POLARS_DTYPE[v.type])


def _frame_from(pairs: list) -> pl.DataFrame | None:
    return pl.DataFrame([_column(n, v) for n, v in pairs]) if pairs else None


def _rvector(arr: Any, na_index: bool = False, shift: int = 0) -> RVector:
    """ndarray result -> RVector. shift / na_index convert a 0-based index array (NA_INDEX -> NA) to R's base."""
    arr = np.asarray(arr)
    t = _RVECTOR_TYPE[arr.dtype.kind]
    vals = arr.tolist()
    if na_index:
        vals = [None if v == NA_INDEX else v + shift for v in vals]
    elif shift:
        vals = [v + shift for v in vals]
    return RVector(t, vals)


def _hits_or_index(res: Any) -> Any:
    return res if isinstance(res, Hits) else _rvector(res, na_index=True, shift=1)


def _subscript(x: Any, i: Any) -> Any:
    if i is MISSING:
        raise Unsupported("subscript:missing")
    if i is None:
        raise Unsupported("subscript:NULL")
    vals = i.list()
    if any(v is None for v in vals):
        raise Unsupported("subscript:NA")
    if i.type == "int":
        if any(v <= 0 for v in vals):
            raise Unsupported("subscript:negative-or-zero")
        return np.array(vals, dtype=np.int64) - 1
    if i.type == "lgl":
        if len(vals) != len(x):
            raise Unsupported("subscript:logical-recycling")
        return np.array(vals, dtype=bool)
    if i.type == "chr":
        return vals
    raise Unsupported(f"subscript:{i.type}")


class Adapter(AdapterBase):
    port = "e2"
    port_packages = ("ranges_e2",)
    unsupported_ops: dict[str, str] = {}

    # ---------------------------------------------------------------- inputs

    def _build_frame(self, pairs: list) -> pl.DataFrame | None:
        cols = []
        for name, col in pairs:
            kind = col["kind"]
            if kind == "vector":
                cols.append(pl.Series(name, col["values"] if col["type"] != "dbl" else RVector("dbl", col["values"]).list(),
                                      dtype=_POLARS_DTYPE[col["type"]]))
            elif kind == "Factor":
                lv = col["levels"]
                cols.append(pl.Series(name, [None if c is None else lv[c - 1] for c in col["codes"]], dtype=pl.Enum(lv)))
            elif kind == "IntegerList":
                # a column named revmap is index-valued: the port stores it 0-based (shift -1 here, +1 on projection)
                shift = -1 if name == _REVMAP else 0
                cols.append(pl.Series(name, [[None if v is None else v + shift for v in e] for e in col["elements"]],
                                      dtype=pl.List(pl.Int32)))
            else:
                raise Unsupported(f"input mcols column {kind}")
        return pl.DataFrame(cols) if cols else None

    def build_IRanges(self, p: dict) -> IRanges:
        return IRanges(start=p["start"], width=p["width"], names=p["names"], mcols=self._build_frame(p["mcols"]))

    def build_Seqinfo(self, p: dict) -> Seqinfo:
        return Seqinfo(p["seqnames"], p["seqlengths"], p["isCircular"], p["genome"])

    def build_GRanges(self, p: dict) -> GRanges:
        si = p["seqinfo"]
        return GRanges(
            p["seqnames"],
            IRanges(start=p["start"], width=p["width"], names=p["names"]),
            p["strand"],
            seqinfo=self.build_Seqinfo(si),
            mcols=self._build_frame(p["mcols"]),
        )

    def build_DataFrame(self, p: dict) -> _DF:
        frame = self._build_frame(p["columns"])
        return _DF(pl.DataFrame() if frame is None else frame, p["nrow"])

    def build_Factor(self, p: dict) -> Factor:
        return Factor(tuple(p["levels"]), np.array([-1 if c is None else c - 1 for c in p["codes"]], dtype=np.int64))

    # ---------------------------------------------------------------- projection

    def _columns(self, frame: pl.DataFrame) -> list:
        out = []
        for name in frame.columns:
            s = frame[name]
            dt = s.dtype
            if dt.is_integer():
                col = p_vector("int", s.to_list())
            elif dt.is_float():
                col = p_vector("dbl", s.to_list())
            elif dt == pl.String:
                col = p_vector("chr", s.to_list())
            elif dt == pl.Boolean:
                col = p_vector("lgl", s.to_list())
            elif isinstance(dt, pl.Enum):
                col = p_factor(dt.categories.to_list(), s.to_physical().to_list(), base=0)
            elif isinstance(dt, pl.List):
                # revmap is the port's documented 0-based list column
                col = p_integerlist(s.to_list(), shift=1 if name == _REVMAP else 0)
            else:
                raise Unsupported(f"project mcols column dtype {dt}")
            out.append((name, col))
        return out

    def _mcols_of(self, x: IRanges | GRanges) -> list:
        return self._columns(x.mcols) if x.has_mcols else []

    def _project_seqinfo(self, si: Seqinfo) -> dict:
        return p_seqinfo(si.seqnames, si.seqlengths, si.is_circular, si.genome)

    def _project_rle(self, r: Rle) -> dict:
        if isinstance(r.values, Factor):
            f = r.values
            return p_vector("chr", Factor(f.levels, r.to_dense()).to_list())
        return p_rle(r.values, r.lengths)

    def project(self, value: Any) -> dict:
        if isinstance(value, IRanges):
            return p_iranges(value.start, value.width, value.names, self._mcols_of(value))
        if isinstance(value, GRanges):
            return p_granges(value.seqname_strings(), value.start, value.width, value.strand_strings(),
                             value.names, self._mcols_of(value), self._project_seqinfo(value.seqinfo))
        if isinstance(value, Seqinfo):
            return self._project_seqinfo(value)
        if isinstance(value, Hits):
            return p_hits(value.query_hits + 1, value.subject_hits + 1, value.query_length, value.subject_length)
        if isinstance(value, Rle):
            return self._project_rle(value)
        if isinstance(value, dict):
            return p_rlelist(list(value), [(r.values, r.lengths) for r in value.values()])
        if isinstance(value, Factor):
            return p_factor(value.levels, value.codes, base=0)
        if isinstance(value, _DF):
            return p_dataframe(value.nrow, self._columns(value.frame))
        if isinstance(value, np.ndarray):
            return p_vector_from(value)
        raise TypeError(f"cannot project {type(value).__name__}")

    # ---------------------------------------------------------------- construct

    def op_IRanges(self, *, start=MISSING, end=MISSING, width=MISSING, names=MISSING, _mcols=()):
        kw = _present(start=_ls(start), end=_ls(end), width=_ls(width), names=_ls(names))
        return IRanges(**kw, mcols=_frame_from(list(_mcols)))

    def op_GRanges(self, *, seqnames=MISSING, ranges=MISSING, strand=MISSING, seqinfo=MISSING, seqlengths=MISSING,
                   _mcols=()):
        kw = _present(seqnames=_ls(seqnames), ranges=ranges, strand=_ls(strand), seqinfo=seqinfo)
        if seqlengths is not MISSING:
            if seqlengths is None or seqlengths.names is None:
                raise Unsupported("argument:seqlengths unnamed")  # the port takes a {seqlevel: length} mapping
            kw["seqlengths"] = _mapping(seqlengths)
        return GRanges(**kw, mcols=_frame_from(list(_mcols)))

    def op_Seqinfo(self, *, seqnames=MISSING, seqlengths=MISSING, isCircular=MISSING, genome=MISSING):
        kw = _present(seqnames=_ls(seqnames), seqlengths=_ls(seqlengths), is_circular=_ls(isCircular),
                      genome=_ls(genome))
        return Seqinfo(**kw)

    # ---------------------------------------------------------------- accessors

    def op_start(self, *, x):
        _need_ranges(x, "start")
        return _rvector(x.start)

    def op_end(self, *, x):
        _need_ranges(x, "end")
        return _rvector(x.end)

    def op_width(self, *, x):
        _need_ranges(x, "width")
        return _rvector(x.width)

    def op_names(self, *, x):
        _need_ranges(x, "names")
        return None if x.names is None else RVector("chr", x.names)

    def op_length(self, *, x):
        return RVector("int", [len(x)])

    def op_seqnames(self, *, x):
        _need_granges(x, "seqnames")
        return x.seqnames

    def op_strand(self, *, x):
        _need_granges(x, "strand")
        return x.strand

    def op_ranges(self, *, x):
        _need_granges(x, "ranges")
        return x.ranges

    def op_mcols(self, *, x):
        _need_ranges(x, "mcols")
        return _DF(x.mcols, len(x))

    # ---------------------------------------------------------------- metadata

    def op_set_names(self, *, x, value):
        _need_ranges(x, "names<-")
        return x.with_names(_ls(value))

    def op_set_mcols(self, *, x, value):
        _need_ranges(x, "mcols<-")
        return x.with_mcols(None if value is None else value.frame)

    # ---------------------------------------------------------------- seqinfo

    def op_seqinfo(self, *, x):
        _need_granges(x, "seqinfo")
        return x.seqinfo

    def op_seqlevels(self, *, x):
        return RVector("chr", list(_seqinfo_of(x, "seqlevels").seqnames))

    def op_seqlengths(self, *, x):
        return RVector("int", list(_seqinfo_of(x, "seqlengths").seqlengths))

    def op_isCircular(self, *, x):
        return RVector("lgl", list(_seqinfo_of(x, "isCircular").is_circular))

    def op_genome(self, *, x):
        return RVector("chr", list(_seqinfo_of(x, "genome").genome))

    def op_set_seqinfo(self, *, x, value):
        _need_granges(x, "seqinfo<-")
        return x.with_seqinfo(value)

    def op_set_seqlevels(self, *, x, value):
        if isinstance(x, Seqinfo):
            return x.rename_levels(value.list())
        _need_granges(x, "seqlevels<-")
        return x.with_seqlevels(value.list())

    def _set_field(self, x: Any, field: str, method: str, value: Any) -> Any:
        if value is not None and value.names is not None:
            raise Unsupported("argument:value named")  # R matches names to seqlevels; the port is positional
        if isinstance(x, Seqinfo):
            return x.with_field(field, _ls(value))
        _need_granges(x, f"{field}<-")
        return getattr(x, method)(_ls(value))

    def op_set_seqlengths(self, *, x, value):
        return self._set_field(x, "seqlengths", "with_seqlengths", value)

    def op_set_isCircular(self, *, x, value):
        return self._set_field(x, "is_circular", "with_is_circular", value)

    def op_set_genome(self, *, x, value):
        return self._set_field(x, "genome", "with_genome", value)

    def op_merge(self, *dots, x, y=MISSING):
        if dots:
            raise Unsupported("merge of more than two Seqinfo")
        if not isinstance(x, Seqinfo):
            raise Unsupported(f"merge on {type(x).__name__}")
        return x.merge(y)

    def op_trim(self, *, x):
        _need_granges(x, "trim")
        return x.trim()

    # ---------------------------------------------------------------- subset / concat

    def op_subset(self, *, x, i=MISSING):
        _need_ranges(x, "[")
        return x[_subscript(x, i)]

    def op_c(self, *dots, x):
        _need_ranges(x, "c")
        return x.concat(*dots)

    # ---------------------------------------------------------------- intra-range

    def op_shift(self, *, x, shift=MISSING):
        return _call(x, "shift", **_present(shift=_ls(shift)))

    def op_narrow(self, *, x, start=MISSING, end=MISSING, width=MISSING):
        return _call(x, "narrow", start=_ls(start), end=_ls(end), width=_ls(width))

    def op_resize(self, *, x, width=MISSING, fix=MISSING, use_names=MISSING, ignore_strand=MISSING):
        return _call(x, "resize", width=_ls(width), fix=_ls(fix), use_names=_sc(use_names),
                     ignore_strand=_sc(ignore_strand))

    def op_flank(self, *, x, width=MISSING, start=MISSING, both=MISSING, use_names=MISSING, ignore_strand=MISSING):
        return _call(x, "flank", width=_ls(width), start=_ls(start), both=_sc(both), ignore_strand=_sc(ignore_strand))

    def op_promoters(self, *, x, upstream=MISSING, downstream=MISSING):
        return _call(x, "promoters", upstream=_ls(upstream), downstream=_ls(downstream))

    # ---------------------------------------------------------------- inter-range

    def op_range(self, *dots, x, ignore_strand=MISSING):
        if dots and isinstance(x, GRanges):
            raise Unsupported("range of several GRanges")
        return _call(x, "range", *dots, ignore_strand=_sc(ignore_strand))

    def op_reduce(self, *, x, drop_empty_ranges=MISSING, min_gapwidth=MISSING, with_revmap=MISSING,
                  ignore_strand=MISSING):
        return _call(x, "reduce", drop_empty_ranges=_sc(drop_empty_ranges), min_gapwidth=_sc(min_gapwidth),
                     with_revmap=_sc(with_revmap), ignore_strand=_sc(ignore_strand))

    def op_gaps(self, *, x, start=MISSING, end=MISSING, ignore_strand=MISSING):
        def bound(v: Any) -> Any:
            if v is MISSING or v is None or not isinstance(x, GRanges):
                return _sc(v)
            # GRanges.gaps documents a {seqlevel: bound} mapping, a sequence, or a scalar
            return _mapping(v) if v.names is not None else v.list()

        return _call(x, "gaps", start=bound(start), end=bound(end), ignore_strand=_sc(ignore_strand))

    def op_disjoin(self, *, x, with_revmap=MISSING, ignore_strand=MISSING):
        return _call(x, "disjoin", with_revmap=_sc(with_revmap), ignore_strand=_sc(ignore_strand))

    # ---------------------------------------------------------------- set operations

    def op_union(self, *, x, y=MISSING, ignore_strand=MISSING):
        return _call(x, "union", y, ignore_strand=_sc(ignore_strand))

    def op_intersect(self, *, x, y=MISSING, ignore_strand=MISSING):
        return _call(x, "intersect", y, ignore_strand=_sc(ignore_strand))

    def op_setdiff(self, *, x, y=MISSING, ignore_strand=MISSING):
        return _call(x, "setdiff", y, ignore_strand=_sc(ignore_strand))

    def op_punion(self, *, x, y=MISSING, fill_gap=MISSING, ignore_strand=MISSING):
        return _call(x, "punion", y, fill_gap=_sc(fill_gap), ignore_strand=_sc(ignore_strand))

    def op_pintersect(self, *, x, y=MISSING, resolve_empty=MISSING, drop_nohit_ranges=MISSING, ignore_strand=MISSING):
        return _call(x, "pintersect", y, resolve_empty=_sc(resolve_empty), drop_nohit_ranges=_sc(drop_nohit_ranges),
                     ignore_strand=_sc(ignore_strand))

    def op_psetdiff(self, *, x, y=MISSING, ignore_strand=MISSING):
        return _call(x, "psetdiff", y, ignore_strand=_sc(ignore_strand))

    # ---------------------------------------------------------------- overlaps

    def _overlap_args(self, subject: Any, maxgap: Any, minoverlap: Any, type: Any) -> dict:
        if subject is MISSING:
            raise Unsupported("argument:subject missing (self overlaps)")
        return _present(maxgap=_sc(maxgap), minoverlap=_sc(minoverlap), type=_sc(type))

    def op_findOverlaps(self, *, query, subject=MISSING, maxgap=MISSING, minoverlap=MISSING, type=MISSING,
                        select=MISSING, ignore_strand=MISSING):
        kw = self._overlap_args(subject, maxgap, minoverlap, type)
        return _hits_or_index(_call(query, "find_overlaps", subject, select=_sc(select),
                                    ignore_strand=_sc(ignore_strand), **kw))

    def op_countOverlaps(self, *, query, subject=MISSING, maxgap=MISSING, minoverlap=MISSING, type=MISSING,
                         ignore_strand=MISSING):
        kw = self._overlap_args(subject, maxgap, minoverlap, type)
        return _rvector(_call(query, "count_overlaps", subject, ignore_strand=_sc(ignore_strand), **kw))

    def op_overlapsAny(self, *, query, subject=MISSING, maxgap=MISSING, minoverlap=MISSING, type=MISSING,
                       ignore_strand=MISSING):
        kw = self._overlap_args(subject, maxgap, minoverlap, type)
        return _rvector(_call(query, "overlaps_any", subject, ignore_strand=_sc(ignore_strand), **kw))

    def op_subsetByOverlaps(self, *, x, ranges=MISSING, maxgap=MISSING, minoverlap=MISSING, type=MISSING,
                            invert=MISSING, ignore_strand=MISSING):
        kw = self._overlap_args(ranges, maxgap, minoverlap, type)
        return _call(x, "subset_by_overlaps", ranges, invert=_sc(invert), ignore_strand=_sc(ignore_strand), **kw)

    # ---------------------------------------------------------------- nearest

    def _nearest(self, method: str, x: Any, subject: Any, select: Any, ignore_strand: Any) -> Any:
        args = () if subject is MISSING else (subject,)  # port: subject=None means self
        return _hits_or_index(_call(x, method, *args, select=_sc(select), ignore_strand=_sc(ignore_strand)))

    def op_nearest(self, *, x, subject=MISSING, select=MISSING, ignore_strand=MISSING):
        return self._nearest("nearest", x, subject, select, ignore_strand)

    def op_precede(self, *, x, subject=MISSING, select=MISSING, ignore_strand=MISSING):
        return self._nearest("precede", x, subject, select, ignore_strand)

    def op_follow(self, *, x, subject=MISSING, select=MISSING, ignore_strand=MISSING):
        return self._nearest("follow", x, subject, select, ignore_strand)

    def op_distance(self, *, x, y=MISSING, ignore_strand=MISSING):
        res = _call(x, "distance", y, ignore_strand=_sc(ignore_strand))
        # GRanges.distance documents NA_INDEX as NA; IRanges.distance never returns it
        return _rvector(res, na_index=isinstance(x, GRanges))

    # ---------------------------------------------------------------- coverage

    def op_coverage(self, *, x, shift=MISSING, width=MISSING, weight=MISSING):
        def width_arg(v: Any) -> Any:
            if v is MISSING or v is None or not isinstance(x, GRanges):
                return _sc(v)
            return _mapping(v) if v.names is not None else _sc(v)

        def weight_arg(v: Any) -> Any:
            # GRanges.coverage documents a string weight as an mcols column name
            if v is MISSING or v is None:
                return v
            return _sc(v) if v.type == "chr" else v.list()

        return _call(x, "coverage", shift=_ls(shift), width=width_arg(width), weight=weight_arg(weight))

    # ---------------------------------------------------------------- ordering

    def op_order(self, *dots, decreasing=MISSING):
        if len(dots) != 1:
            raise Unsupported("order of several keys")
        x = dots[0]
        return _rvector(_call(x, "order", decreasing=_sc(decreasing)), shift=1)

    def op_sort(self, *, x, decreasing=MISSING, ignore_strand=MISSING):
        return _call(x, "sort", decreasing=_sc(decreasing), ignore_strand=_sc(ignore_strand))

    def op_rank(self, *, x, ties_method=MISSING, ignore_strand=MISSING):
        return _rvector(_call(x, "rank", ties_method=_sc(ties_method), ignore_strand=_sc(ignore_strand)))

    def op_duplicated(self, *, x, fromLast=MISSING):
        return _rvector(_call(x, "duplicated", from_last=_sc(fromLast)))

    def op_unique(self, *, x, fromLast=MISSING):
        return _call(x, "unique", from_last=_sc(fromLast))
