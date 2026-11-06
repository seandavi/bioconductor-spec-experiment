"""Adapter for BiocPy: PyPI `iranges` 0.7.3 (IRanges) and `genomicranges` 0.9.0 (GenomicRanges, SeqInfo).

Port docs used: the IRanges and GenomicRanges READMEs (package METADATA) and the public docstrings and
signatures of `IRanges`, `GenomicRanges`, `SeqInfo`, `genomicranges.sequence_info.merge_SeqInfo`, plus the
`biocframe.BiocFrame` / `biocutils.Factor` / `biocutils.combine_sequences` containers the port is built on.
Unsupported ops and conventions in prose: eval/adapters/biocpy.md.

Index base: row indices (hits, order, nearest ..., subscripts, reverse maps) are 0-based; coordinates
(`start`, `end`) are used as given (the READMEs say 1-based with inclusive ends, as in Bioconductor).

Translations (representation only):
  * argument names R -> port (`maxgap` -> `max_gap`, `with.revmap` -> `with_reverse_map`, `type` -> `query_type`,
    `y` -> `other`, ...). An R argument the port method lacks raises UnknownArgument; a method the port lacks raises
    Unsupported (see `_invoke`).
  * roles: R `findOverlaps(query, subject)` -> `subject.find_overlaps(query)`; likewise countOverlaps, nearest,
    precede, follow. `subsetByOverlaps(x, ranges)` -> `x.subset_by_overlaps(ranges)`.
  * length-1 R vectors are passed as Python scalars (the port's documented recycling form), longer ones as
    ndarray / list. A named R vector given to a per-seqlevel setter (`seqlengths<-`, `isCircular<-`, `genome<-`) or
    to `gaps(end=)` is passed as the dict the port documents.
  * 1-based -> 0-based: positive integer `[` subscripts. 0-based -> 1-based: Hits, `order`, the integer result of
    nearest / precede / follow, the `revmap` mcols column (only when with.revmap = TRUE).
  * hits (BiocFrame `query_hits` / `self_hits`) -> Hits(from = query, to = subject, nLnode = len(query),
    nRnode = len(subject)).
  * dense coverage arrays -> Rle (runs normalised by the base helper); the `GenomicRanges.coverage` dict -> RleList.
  * seqlengths / isCircular / genome accessors return plain lists in seqlevel order; R returns them named by
    seqlevels, so the names are attached from `get_seqnames()`.
  * missing `subject` in findOverlaps / countOverlaps (self overlaps): the object is passed as its own query.
  * `range(x, ...)`: objects in `...` are first combined with the port's combine, then `range` is called.
  * `seqlengths<-` / `isCircular<-` / `genome<-` on a GRanges: `get_seqinfo()`, the SeqInfo setter, `set_seqinfo()`
    (the port has no per-GRanges setter).
  * mcols columns: numpy arrays (masked when NA is present), `biocutils.Factor`, list of int arrays (IntegerList).
"""

from __future__ import annotations

import inspect
from typing import Any, Callable

import biocutils as ut
import numpy as np
from biocframe import BiocFrame
from genomicranges import GenomicRanges, SeqInfo
from genomicranges.sequence_info import merge_SeqInfo
from iranges import IRanges

from adapter_base import (
    MISSING,
    AdapterBase,
    RVector,
    UnknownArgument,
    Unsupported,
    infer_type,
    p_dataframe,
    p_factor,
    p_granges,
    p_hits,
    p_integerlist,
    p_iranges,
    p_rle,
    p_rlelist,
    p_seqinfo,
    p_vector_from,
)


# ---------------------------------------------------------------- native wrappers (the port returns bare arrays / frames)

class _Hits:
    """A hits BiocFrame (`query_hits`, `self_hits`; 0-based) with the sizes of both sides."""

    def __init__(self, frame: BiocFrame, n_query: int, n_subject: int):
        self.frame, self.n_query, self.n_subject = frame, n_query, n_subject


class _Coverage:
    """Dense coverage: an ndarray (IRanges) or a dict seqname -> ndarray (GenomicRanges)."""

    def __init__(self, value: Any):
        self.value = value

    def __len__(self) -> int:
        return len(self.value)


class _Revmap:
    """A result computed with `with_reverse_map=True`: its `revmap` mcols column holds 0-based indices."""

    def __init__(self, obj: Any):
        self.obj = obj


# ---------------------------------------------------------------- R value -> port value

def _seq(rv: RVector) -> list:
    """R vector -> Python list (None = NA); the port documents `Sequence` for these arguments."""
    return rv.list()


def _rvector_column(rv: RVector) -> np.ndarray:
    """R vector -> typed numpy column; a masked array when it holds NA."""
    if rv.type == "chr":
        mask = [v is None for v in rv.values]
        arr = np.array(["" if m else v for v, m in zip(rv.values, mask)], dtype=str)
        return np.ma.masked_array(arr, mask=mask) if any(mask) else arr
    arr = rv.np_masked()
    return arr if np.ma.is_masked(arr) else arr.data


def _mcols_frame(columns: list, nrow: int) -> BiocFrame:
    """Extra named R arguments (`...`) -> BiocFrame; vectors become typed columns, port objects pass through."""
    return BiocFrame({name: _rvector_column(v) if isinstance(v, RVector) else v for name, v in columns},
                     number_of_rows=nrow)


def _as_rvector(value: Any) -> Any:
    """A nested call's native vector result (ndarray, masked array, non-empty list) as an R vector, NA kept."""
    if isinstance(value, np.ndarray) and value.ndim == 1:
        mask = np.ma.getmaskarray(value).tolist() if isinstance(value, np.ma.MaskedArray) else [False] * len(value)
        data = value.data if isinstance(value, np.ma.MaskedArray) else value
        vals = [None if m else v for v, m in zip(data.tolist(), mask)]
        return RVector(infer_type(data), vals)
    if isinstance(value, list) and value:
        return RVector(infer_type(value), value)
    return value


def _arr(rv: RVector) -> Any:
    """R vector -> ndarray (a list when it holds NA in an int / lgl vector), for arguments the port types as ndarray."""
    if rv.type == "chr" or (rv.type != "dbl" and any(v is None for v in rv.values)):
        return rv.list()
    return rv.np(na_as_nan=True)


def _scalar_only(kw: dict, names: tuple, what: str) -> None:
    """The port method types these arguments as a scalar `int`; a longer vector (R recycles / pairs it) is a form it lacks."""
    for name in names:
        if isinstance(kw.get(name), RVector) and len(kw[name]) > 1:
            raise Unsupported(f"{what}: {name} as a vector")


def _no_incomparables(rv: Any) -> Any:
    """R `incomparables = FALSE` (nothing is incomparable) -> the container biocutils takes; its default is `set()`."""
    if isinstance(rv, RVector) and rv.type == "lgl" and rv.values == [False]:
        return set()
    raise Unsupported("duplicated / unique: incomparables other than FALSE")


def _val(rv: Any) -> Any:
    """Length-1 R vector -> Python scalar; otherwise a list; NULL / port objects unchanged."""
    if not isinstance(rv, RVector):
        return rv
    return rv.scalar() if len(rv) == 1 else _seq(rv)


def _named_or_val(rv: Any) -> Any:
    """A named R vector -> {name: value}; otherwise `_val`."""
    if isinstance(rv, RVector) and rv.names is not None:
        return dict(zip(rv.names, rv.list()))
    return _val(rv)


def _index_vector(arr: Any) -> np.ma.MaskedArray:
    """0-based index array (None / NaN = no hit) -> 1-based int vector, no hit masked."""
    vals = list(arr)
    mask = [v is None or (isinstance(v, float) and np.isnan(v)) for v in vals]
    data = [0 if m else int(v) + 1 for v, m in zip(vals, mask)]
    return np.ma.masked_array(np.array(data, dtype=np.int64), mask=mask)


def _subscript(i: Any, n: int) -> list:
    """R `[` subscript -> port subscript. Only forms whose R and port meanings coincide are passed."""
    if not isinstance(i, RVector):
        raise Unsupported("subscript:NULL")
    vals = i.list()
    if any(v is None for v in vals):
        raise Unsupported("subscript:NA")
    if i.type == "int":
        if any(v <= 0 for v in vals):
            raise Unsupported("subscript:negative-or-zero")
        return [v - 1 for v in vals]
    if i.type == "lgl":
        if len(vals) != n:
            raise Unsupported("subscript:logical-recycling")
        return vals
    if i.type == "chr":
        return vals
    raise Unsupported(f"subscript:{i.type}")


# ---------------------------------------------------------------- projection pieces

def _project_column(col: Any) -> dict:
    if isinstance(col, ut.Factor):
        return p_factor(col.get_levels(), col.get_codes(), base=0)
    if isinstance(col, np.ndarray):
        return p_vector_from(col)
    if isinstance(col, list):
        if all(isinstance(e, (np.ndarray, list, tuple)) for e in col):  # lists hold int arrays; atomic columns are ndarrays
            return p_integerlist(col)
        return p_vector_from(col)
    raise TypeError(f"cannot project mcols column of type {type(col).__name__}")


def _project_columns(frame: BiocFrame | None) -> list:
    if frame is None:
        return []
    return [(name, _project_column(frame.get_column(name))) for name in frame.get_column_names()]


def _require(**objects: Any) -> None:
    """The port method has a required parameter for this R argument, and the case omits it."""
    for name, value in objects.items():
        if value is MISSING:
            raise Unsupported(f"argument missing: {name}")


def _names_list(names: Any) -> list | None:
    return None if names is None else list(names)


def _project_seqinfo(si: SeqInfo) -> dict:
    return p_seqinfo(si.get_seqnames(), si.get_seqlengths(), si.get_is_circular(), si.get_genome())


class Adapter(AdapterBase):
    port = "biocpy"
    port_packages = ("iranges", "genomicranges", "biocframe", "biocutils", "compressed_lists")
    unsupported_ops = {
        "overlapsAny": "no public overlapsAny equivalent (only find_overlaps / count_overlaps / subset_by_overlaps)",
        "punion": "no public parallel union",
        "psetdiff": "no public parallel set difference",
    }

    # ---------------------------------------------------------------- inputs

    def eval_arg(self, a: Any, env: dict) -> Any:
        value = super().eval_arg(a, env)
        # a nested call's native vector is an R vector to the next call
        return _as_rvector(value) if isinstance(a, dict) and "$call" in a else value

    def _column(self, p: dict) -> Any:
        kind = p["kind"]
        if kind == "Factor":
            return self.build_Factor(p)
        if kind == "IntegerList":
            return [np.ma.masked_array(np.array([0 if v is None else v for v in e], dtype=np.int32),
                                       mask=[v is None for v in e]) for e in p["elements"]]
        return _rvector_column(RVector(p["type"], p["values"]))

    def _frame(self, columns: list, nrow: int, rownames: list | None = None) -> BiocFrame:
        return BiocFrame({name: self._column(col) for name, col in columns}, number_of_rows=nrow, row_names=rownames)

    def build_IRanges(self, p: dict) -> IRanges:
        mcols = self._frame(p["mcols"], len(p["start"])) if p["mcols"] else None
        return IRanges(start=p["start"], width=p["width"], names=p["names"], mcols=mcols)

    def build_Seqinfo(self, p: dict) -> SeqInfo:
        return SeqInfo(p["seqnames"], p["seqlengths"], p["isCircular"], p["genome"])

    def build_GRanges(self, p: dict) -> GenomicRanges:
        mcols = self._frame(p["mcols"], len(p["start"])) if p["mcols"] else None
        return GenomicRanges(
            seqnames=p["seqnames"],
            ranges=IRanges(start=p["start"], width=p["width"]),
            strand=p["strand"],
            names=p["names"],
            mcols=mcols,
            seqinfo=self.build_Seqinfo(p["seqinfo"]),
        )

    def build_DataFrame(self, p: dict) -> BiocFrame:
        return self._frame(p["columns"], p["nrow"], p["rownames"])

    def build_Factor(self, p: dict) -> ut.Factor:
        return ut.Factor([-1 if c is None else c - 1 for c in p["codes"]], p["levels"])

    # ---------------------------------------------------------------- dispatch helper

    def _invoke(self, obj: Any, method: str, kw: dict, argmap: dict, *pos: Any,
                conv: dict[str, Callable] | None = None) -> Any:
        """Call `obj.method(*pos, **mapped kw)`.

        `argmap`: R argument name (dots as `_`) -> port parameter name. An R argument without a mapping, or whose
        port parameter the method lacks, is UnknownArgument; a missing method is Unsupported.
        """
        fn = getattr(obj, method, None)
        if fn is None:
            raise Unsupported(f"{type(obj).__name__} has no public {method}")
        return self._call_mapped(fn, kw, argmap, *pos, conv=conv)

    def _call_mapped(self, fn: Callable, kw: dict, argmap: dict, *pos: Any,
                     conv: dict[str, Callable] | None = None) -> Any:
        params = inspect.signature(fn).parameters
        call: dict[str, Any] = {}
        for rname, value in kw.items():
            pname = argmap.get(rname)
            if pname is None or pname not in params:
                raise UnknownArgument(rname)
            call[pname] = (conv or {}).get(rname, _val)(value)
        return fn(*pos, **call)

    def _revmap(self, result: Any, kw: dict) -> Any:
        flag = kw.get("with_revmap")
        return _Revmap(result) if isinstance(flag, RVector) and flag.values == [True] else result

    def _search(self, res: Any, n_query: int, n_subject: int, kw: dict | None = None) -> Any:
        """Result of a find_overlaps / nearest / precede / follow call."""
        if isinstance(res, BiocFrame):
            select = (kw or {}).get("select")
            if isinstance(select, RVector) and select.type == "chr" and select.values in (["first"], ["last"], ["arbitrary"]):
                return self._selected_hits(res, n_query) or _Hits(res, n_query, n_subject)
            return _Hits(res, n_query, n_subject)
        return _index_vector(res)

    @staticmethod
    def _selected_hits(frame: BiocFrame, n_query: int) -> RVector | None:
        """`select != "all"`: the port's frame holds one subject per matched query; R's result is the integer vector
        of length n_query (1-based subject, NA = no hit). None when the frame is not a list of integer pairs with one
        row per query (nothing to project to)."""
        out: list = [None] * n_query
        for q, s in zip(frame.get_column("query_hits").tolist(), frame.get_column("self_hits").tolist()):
            if not (isinstance(q, int) and isinstance(s, int)) or out[q] is not None:
                return None
            out[q] = s + 1
        return RVector("int", out)

    # ---------------------------------------------------------------- construct

    def op_IRanges(self, *, start=MISSING, end=MISSING, width=MISSING, names=MISSING, _mcols=()):
        if end is not MISSING:
            raise UnknownArgument("end")
        kw: dict[str, Any] = {}
        if start is not MISSING:
            kw["start"] = _seq(start)
        if width is not MISSING:
            kw["width"] = _seq(width)
        if names is not MISSING:
            kw["names"] = None if names is None else names.list()
        if _mcols:
            n = len(kw["start"]) if "start" in kw else 0
            kw["mcols"] = _mcols_frame(_mcols, n)
        return IRanges(**kw)

    def op_GRanges(self, *, seqnames=MISSING, ranges=MISSING, strand=MISSING, seqinfo=MISSING,
                   seqlengths=MISSING, _mcols=()):
        if seqlengths is not MISSING:
            raise UnknownArgument("seqlengths")
        if isinstance(seqinfo, RVector):
            raise Unsupported("argument:seqinfo as a character vector; the port takes a SeqInfo")
        if all(v is MISSING for v in (seqnames, ranges, strand, seqinfo)) and not _mcols:
            return GenomicRanges.empty()  # the port's documented factory for the no-argument form
        _require(seqnames=seqnames, ranges=ranges)  # GenomicRanges.__init__ has no defaults for these
        kw: dict[str, Any] = {"seqnames": seqnames.list() if isinstance(seqnames, RVector) else seqnames,
                              "ranges": ranges}
        if strand is not MISSING:
            kw["strand"] = None if strand is None else strand.list()
        if seqinfo is not MISSING:
            kw["seqinfo"] = seqinfo
        if _mcols:
            kw["mcols"] = _mcols_frame(_mcols, len(kw["seqnames"]))
        return GenomicRanges(**kw)

    def op_Seqinfo(self, *, seqnames=MISSING, seqlengths=MISSING, isCircular=MISSING, genome=MISSING):
        if all(v is MISSING for v in (seqnames, seqlengths, isCircular, genome)):
            return SeqInfo.empty()  # the port's documented factory for the no-argument form
        _require(seqnames=seqnames)  # SeqInfo.__init__ has no default for seqnames
        kw: dict[str, Any] = {"seqnames": seqnames.list()}
        if seqlengths is not MISSING:
            kw["seqlengths"] = _val(seqlengths)
        if isCircular is not MISSING:
            kw["is_circular"] = _val(isCircular)
        if genome is not MISSING:
            kw["genome"] = _val(genome)
        return SeqInfo(**kw)

    # ---------------------------------------------------------------- accessors

    def op_start(self, *, x):
        return x.get_start()

    def op_end(self, *, x):
        return x.get_end()

    def op_width(self, *, x):
        return x.get_width()

    def op_names(self, *, x):
        names = x.get_names()
        return None if names is None else RVector("chr", list(names))

    def op_length(self, *, x):
        return RVector("int", [len(x)])

    def op_seqnames(self, *, x):
        return RVector("chr", list(x.get_seqnames()))

    def op_strand(self, *, x):
        return RVector("chr", list(x.get_strand(as_type="list")))

    def op_ranges(self, *, x, **kw):
        return self._invoke(x, "get_ranges", kw, {})

    def op_mcols(self, *, x, **kw):
        return self._invoke(x, "get_mcols", kw, {})

    # ---------------------------------------------------------------- metadata

    def op_set_names(self, *, x, value):
        return x.set_names(None if value is None else value.list())

    def op_set_mcols(self, *, x, value):
        return x.set_mcols(value)

    # ---------------------------------------------------------------- seqinfo

    def _seqinfo_of(self, x: Any) -> SeqInfo:
        if isinstance(x, SeqInfo):
            return x
        if isinstance(x, GenomicRanges):
            return x.get_seqinfo()
        raise Unsupported(f"{type(x).__name__} has no seqinfo")

    def op_seqinfo(self, *, x):
        return self._seqinfo_of(x)

    def op_seqlevels(self, *, x):
        return RVector("chr", list(self._seqinfo_of(x).get_seqnames()))

    def _named_by_seqlevel(self, x: Any, type_: str, values: list) -> RVector:
        return RVector(type_, list(values), list(self._seqinfo_of(x).get_seqnames()))

    def op_seqlengths(self, *, x):
        return self._named_by_seqlevel(x, "int", self._seqinfo_of(x).get_seqlengths())

    def op_isCircular(self, *, x):
        return self._named_by_seqlevel(x, "lgl", self._seqinfo_of(x).get_is_circular())

    def op_genome(self, *, x):
        return self._named_by_seqlevel(x, "chr", self._seqinfo_of(x).get_genome())

    def op_set_seqinfo(self, *, x, value, **kw):
        for rname in kw:  # new2old, pruning.mode: the port's set_seqinfo takes only the new SeqInfo
            raise UnknownArgument(rname)
        return x.set_seqinfo(value)

    def _set_seqinfo_column(self, x: Any, setter: str, arg: Any) -> Any:
        if isinstance(x, SeqInfo):
            return getattr(x, setter)(arg)
        return x.set_seqinfo(getattr(x.get_seqinfo(), setter)(arg))

    def op_set_seqlevels(self, *, x, value, **kw):
        # only the port's positional rename (`SeqInfo.set_seqnames`); it has no match-by-name / drop / add form
        for rname in kw:  # pruning.mode
            raise UnknownArgument(rname)
        if not isinstance(value, RVector) or value.names is not None:
            raise Unsupported("seqlevels<-: only an unnamed character value maps to SeqInfo.set_seqnames")
        if len(value) != len(self._seqinfo_of(x).get_seqnames()):  # set_seqnames documents "length equal to the number of names"
            raise Unsupported("seqlevels<-: SeqInfo.set_seqnames needs as many names as current levels (no add / drop form)")
        return self._set_seqinfo_column(x, "set_seqnames", value.list())

    def op_set_seqlengths(self, *, x, value):
        return self._set_seqinfo_column(x, "set_seqlengths", None if value is None else _named_or_val(value))

    def op_set_isCircular(self, *, x, value):
        return self._set_seqinfo_column(x, "set_is_circular", None if value is None else _named_or_val(value))

    def op_set_genome(self, *, x, value):
        return self._set_seqinfo_column(x, "set_genome", None if value is None else _named_or_val(value))

    def op_merge(self, *dots, x, y=MISSING):
        objects = [x] + ([] if y is MISSING else [y]) + list(dots)
        return merge_SeqInfo(objects)

    def op_trim(self, *, x, **kw):
        return self._invoke(x, "trim", kw, {})

    # ---------------------------------------------------------------- subset / concatenate

    def op_subset(self, *, x, i=MISSING, j=MISSING, drop=MISSING):
        if j is not MISSING:
            raise UnknownArgument("j")
        if drop is not MISSING:
            raise UnknownArgument("drop")
        if i is MISSING:
            raise Unsupported("subscript:missing")
        return x[_subscript(i, len(x))]

    def _combine(self, x: Any, dots: tuple) -> Any:
        if isinstance(x, IRanges):
            return x.combine(*dots)
        return ut.combine_sequences(x, *dots)  # GenomicRanges registers with biocutils; it has no combine method

    def op_c(self, *dots, x):
        if any(isinstance(d, RVector) for d in dots):
            raise Unsupported("c: atomic operand")
        return self._combine(x, dots)

    # ---------------------------------------------------------------- intra-range

    def op_shift(self, *, x, **kw):
        return self._invoke(x, "shift", kw, {"shift": "shift"})

    def op_narrow(self, *, x, **kw):
        return self._invoke(x, "narrow", kw, {"start": "start", "end": "end", "width": "width"})

    def op_resize(self, *, x, **kw):
        return self._invoke(x, "resize", kw, {"width": "width", "fix": "fix", "ignore_strand": "ignore_strand"})

    def op_flank(self, *, x, **kw):
        return self._invoke(x, "flank", kw, {"width": "width", "start": "start", "both": "both",
                                             "ignore_strand": "ignore_strand"})

    def op_promoters(self, *, x, **kw):
        _scalar_only(kw, ("upstream", "downstream"), "promoters")
        return self._invoke(x, "promoters", kw, {"upstream": "upstream", "downstream": "downstream"})

    # ---------------------------------------------------------------- inter-range

    def op_range(self, *dots, x, **kw):
        if dots:
            x = self._combine(x, dots)
        res = self._invoke(x, "range", kw, {"with_revmap": "with_reverse_map", "ignore_strand": "ignore_strand"})
        return self._revmap(res, kw)

    def op_reduce(self, *, x, **kw):
        res = self._invoke(x, "reduce", kw, {"drop_empty_ranges": "drop_empty_ranges",
                                             "min_gapwidth": "min_gap_width",
                                             "with_revmap": "with_reverse_map",
                                             "ignore_strand": "ignore_strand"})
        return self._revmap(res, kw)

    def op_gaps(self, *, x, **kw):
        return self._invoke(x, "gaps", kw, {"start": "start", "end": "end", "ignore_strand": "ignore_strand"},
                            conv={"end": _named_or_val})

    def op_disjoin(self, *, x, **kw):
        res = self._invoke(x, "disjoin", kw, {"with_revmap": "with_reverse_map", "ignore_strand": "ignore_strand"})
        return self._revmap(res, kw)

    # ---------------------------------------------------------------- set operations

    def op_union(self, *, x, y=MISSING, **kw):
        _require(y=y)
        return self._invoke(x, "union", kw, {"ignore_strand": "ignore_strand"}, y)

    def op_intersect(self, *, x, y=MISSING, **kw):
        _require(y=y)
        return self._invoke(x, "intersect", kw, {"ignore_strand": "ignore_strand"}, y)

    def op_setdiff(self, *, x, y=MISSING, **kw):
        _require(y=y)
        return self._invoke(x, "setdiff", kw, {"ignore_strand": "ignore_strand"}, y)

    def op_pintersect(self, *, x, y=MISSING, **kw):
        _require(y=y)
        return self._invoke(x, "pintersect", kw, {"ignore_strand": "ignore_strand"}, y)

    # ---------------------------------------------------------------- overlaps

    _OVERLAP_ARGS = {"maxgap": "max_gap", "minoverlap": "min_overlap", "type": "query_type",
                     "ignore_strand": "ignore_strand"}

    @staticmethod
    def _overlap_type(kw: dict) -> None:
        """BiocPy documents query_type any / start / end / within only."""
        t = kw.get("type")
        if isinstance(t, RVector) and t.values == ["equal"]:
            raise Unsupported("type=equal")

    def op_findOverlaps(self, *, query, subject=MISSING, **kw):
        self._overlap_type(kw)
        if subject is MISSING:  # R: missing subject = self overlaps of `query`
            subject = query
        argmap = {**self._OVERLAP_ARGS, "select": "select", "nthread": "num_threads"}
        res = self._invoke(subject, "find_overlaps", kw, argmap, query)
        return self._search(res, len(query), len(subject), kw)

    def op_countOverlaps(self, *, query, subject=MISSING, **kw):
        self._overlap_type(kw)
        if subject is MISSING:  # as findOverlaps: missing subject = the query itself
            subject = query
        argmap = {**self._OVERLAP_ARGS, "nthread": "num_threads"}
        return self._invoke(subject, "count_overlaps", kw, argmap, query)

    def op_subsetByOverlaps(self, *, x, ranges=MISSING, **kw):
        self._overlap_type(kw)
        _require(ranges=ranges)
        return self._invoke(x, "subset_by_overlaps", kw, self._OVERLAP_ARGS, ranges)

    # ---------------------------------------------------------------- nearest / distance

    def _search_op(self, method: str, x: Any, subject: Any, kw: dict) -> Any:
        if subject is MISSING:
            raise Unsupported("subject missing: the port has no self-search form")
        res = self._invoke(subject, method, kw, {"select": "select", "ignore_strand": "ignore_strand"}, x)
        return self._search(res, len(x), len(subject), kw)

    def op_nearest(self, *, x, subject=MISSING, **kw):
        return self._search_op("nearest", x, subject, kw)

    def op_precede(self, *, x, subject=MISSING, **kw):
        return self._search_op("precede", x, subject, kw)

    def op_follow(self, *, x, subject=MISSING, **kw):
        return self._search_op("follow", x, subject, kw)

    def op_distance(self, *, x, y=MISSING, **kw):
        _require(y=y)
        return self._invoke(x, "distance", kw, {}, y)

    # ---------------------------------------------------------------- coverage

    def op_coverage(self, *, x, **kw):
        for name in ("weight", "shift"):
            if isinstance(kw.get(name), RVector) and kw[name].type == "chr":
                raise Unsupported(f"coverage: {name} as an mcols column name")  # the port takes numbers only
        if isinstance(x, GenomicRanges):
            _scalar_only(kw, ("shift", "width", "weight"), "GenomicRanges.coverage")
        # IRanges.coverage documents array-valued shift / weight; GenomicRanges.coverage documents scalars.
        conv = {"shift": _arr, "weight": _arr} if isinstance(x, IRanges) else {}
        res = self._invoke(x, "coverage", kw, {"shift": "shift", "width": "width", "weight": "weight",
                                               "method": "method"}, conv=conv)
        return _Coverage(res)

    # ---------------------------------------------------------------- ordering

    def op_order(self, *dots, **kw):
        if len(dots) != 1:
            raise Unsupported("order: the port orders one object")
        return _index_vector(self._invoke(dots[0], "order", kw, {"decreasing": "decreasing"}))

    def op_sort(self, *, x, **kw):
        return self._invoke(x, "sort", kw, {"decreasing": "decreasing"})

    def op_rank(self, *, x, **kw):
        return self._invoke(x, "rank", kw, {})

    def op_duplicated(self, *, x, **kw):
        # biocutils.duplicated is the port's public generic; it accepts IRanges and GenomicRanges
        return self._call_mapped(ut.duplicated, kw, {"incomparables": "incomparables", "fromLast": "from_last"}, x,
                                 conv={"incomparables": _no_incomparables})

    def op_unique(self, *, x, **kw):
        return self._call_mapped(ut.unique, kw, {"incomparables": "incomparables", "fromLast": "from_last"}, x,
                                 conv={"incomparables": _no_incomparables})

    # ---------------------------------------------------------------- projection

    def project(self, value: Any) -> dict:
        if isinstance(value, IRanges):
            return p_iranges(value.get_start(), value.get_width(), _names_list(value.get_names()),
                             _project_columns(value.get_mcols()))
        if isinstance(value, GenomicRanges):
            return p_granges(value.get_seqnames(), value.get_start(), value.get_width(),
                             value.get_strand(as_type="list"), _names_list(value.get_names()),
                             _project_columns(value.get_mcols()), _project_seqinfo(value.get_seqinfo()))
        if isinstance(value, SeqInfo):
            return _project_seqinfo(value)
        if isinstance(value, BiocFrame):
            rn = value.get_row_names()
            return p_dataframe(value.shape[0], _project_columns(value), None if rn is None else list(rn))
        if isinstance(value, ut.Factor):
            return p_factor(value.get_levels(), value.get_codes(), base=0)
        if isinstance(value, _Hits):
            frame = value.frame
            return p_hits(np.asarray(frame.get_column("query_hits")) + 1, np.asarray(frame.get_column("self_hits")) + 1,
                          value.n_query, value.n_subject)
        if isinstance(value, _Coverage):
            if isinstance(value.value, dict):
                return p_rlelist(list(value.value.keys()), list(value.value.values()))
            return p_rle(value.value)
        if isinstance(value, _Revmap):
            out = self.project(value.obj)
            for name, col in out["mcols"]:
                if name == "revmap" and col["kind"] == "IntegerList":
                    col["elements"] = [[None if v is None else v + 1 for v in e] for e in col["elements"]]
            return out
        if isinstance(value, (np.ndarray, list)):
            return p_vector_from(value)
        raise TypeError(f"cannot project {type(value).__name__}")
