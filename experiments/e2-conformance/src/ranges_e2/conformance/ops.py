"""Op table: R function name -> adapter onto the ranges_e2 API.

This is the ONLY place where R conventions (1-based indices, NA, R argument names, R subscript
rules) are translated to the Python API (0-based, -1 sentinel, snake_case, numpy subscripts).
Each adapter declares R's positional parameter order; unknown arguments raise TypeError.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import numpy as np
import polars as pl

from .. import GRanges, IRanges, Seqinfo
from ..hits import Hits
from ..rle import Factor
from .codec import DFrameOut, Vec, vec_to_series

OpFn = Callable[["Params"], Any]
Adapter = Callable[[list[Any], dict[str, Any]], Any]
OPS: dict[str, Adapter] = {}


class Params:
    """Arguments bound to R parameter names; `rest` holds `...` positionals, `dots` named extras."""

    def __init__(self, named: dict[str, Any], rest: list[Any], dots: dict[str, Any]) -> None:
        self.named, self.rest, self.dots = named, rest, dots

    def __getitem__(self, key: str) -> Any:
        if key not in self.named:
            raise TypeError(f'argument "{key}" is missing')
        return self.named[key]

    def get(self, key: str, default: Any = None) -> Any:
        return self.named.get(key, default)

    def has(self, key: str) -> bool:
        return key in self.named


def op(name: str, *params: str, varargs: bool = False, dots: bool = False) -> Callable[[OpFn], OpFn]:
    def register(fn: OpFn) -> OpFn:
        def call(args: list[Any], kwargs: dict[str, Any]) -> Any:
            named: dict[str, Any] = {}
            positional = list(args)
            for pname in params:
                if positional:
                    named[pname] = positional.pop(0)
            if positional and not varargs:
                raise TypeError(f"{name}: unused positional arguments")
            extra: dict[str, Any] = {}
            for k, v in kwargs.items():
                if k in named:
                    raise TypeError(f"{name}: formal argument {k} matched by multiple actual arguments")
                if k in params:
                    named[k] = v
                elif dots:
                    extra[k] = v
                else:
                    raise TypeError(f"{name}: unused argument {k}")
            return fn(Params(named, positional, extra))

        OPS[name] = call
        return fn

    return register


def run_op(name: str, args: list[Any], kwargs: dict[str, Any]) -> Any:
    if name not in OPS:
        raise NotImplementedError(f"op {name} is not implemented")
    return OPS[name](args, kwargs)


# ---- argument conversion ------------------------------------------------------------


def _vec(v: Any, what: str = "argument") -> Vec:
    if not isinstance(v, Vec):
        raise TypeError(f"{what} must be an atomic vector")
    return v


def scalar_int(v: Any, what: str = "argument") -> int:
    vec = _vec(v, what)
    if len(vec) != 1 or vec.has_na():
        raise ValueError(f"{what} must be a single non-NA number")
    return int(vec.ints()[0])


def opt_int(v: Any, what: str = "argument") -> int | None:
    """None or a single NA means 'not supplied'."""
    if v is None:
        return None
    vec = _vec(v, what)
    if len(vec) == 1 and vec.has_na():
        return None
    return scalar_int(vec, what)


def scalar_bool(v: Any, what: str = "argument") -> bool:
    vec = _vec(v, what)
    if len(vec) != 1 or vec.kind != "lgl" or vec.has_na():
        raise ValueError(f"{what} must be TRUE or FALSE")
    return bool(vec.values[0])


def scalar_str(v: Any, what: str = "argument") -> str:
    vec = _vec(v, what)
    if len(vec) != 1 or vec.kind != "chr" or vec.has_na():
        raise ValueError(f"{what} must be a single string")
    return str(vec.values[0])


def numbers(v: Any, what: str = "argument") -> np.ndarray | None:
    return None if v is None else _vec(v, what).numeric()


def int_array(v: Any, what: str = "argument") -> np.ndarray | None:
    return None if v is None else _vec(v, what).ints()


def bool_array(v: Any, what: str = "argument") -> np.ndarray:
    return _vec(v, what).bools()


def strings(v: Any, what: str = "argument") -> list[str | None] | None:
    return None if v is None else _vec(v, what).strings()


def names_arg(v: Any, what: str = "argument") -> list[str]:
    """Character vector without NAs."""
    vals = _vec(v, what).strings()
    if any(s is None for s in vals):
        raise ValueError(f"{what} cannot contain NAs")
    return [s for s in vals if s is not None]


def r_subscript(v: Any, n: int, names: list[str | None] | None) -> np.ndarray:
    """R `[` subscript -> 0-based positions. Zero drops out, negatives drop elements, no mixing."""
    from .._util import positions_from_subscript

    vec = _vec(v, "subscript")
    if vec.kind == "chr":
        if vec.has_na():
            raise ValueError("NA subscripts are not supported")
        return positions_from_subscript(list(vec.values), n, names)
    if vec.has_na():
        raise ValueError("NA subscripts are not supported")
    if vec.kind == "lgl":
        flags = np.array(vec.values, dtype=bool)
        if len(flags) > n:
            if flags[n:].any():
                raise IndexError("subscript is a logical vector with out-of-bounds TRUE values")
            flags = flags[:n]
        elif len(flags) < n and len(flags) > 0:
            flags = np.resize(flags, n)
        positions: np.ndarray = np.nonzero(flags)[0].astype(np.int64)
        return positions
    a = vec.ints()
    if (a < 0).any():
        if (a > 0).any():
            raise ValueError("only 0's may be mixed with negative subscripts")
        drop = -a[a < 0] - 1
        kept: np.ndarray = np.setdiff1d(np.arange(n, dtype=np.int64), drop)
        return kept
    a = a[a > 0]
    if (a > n).any():
        raise IndexError("subscript contains out-of-bounds indices")
    return a - 1


def to_r_index(a: np.ndarray) -> Vec:
    """0-based indices with -1 = no hit -> 1-based R integers with NA."""
    return Vec("int", [None if v < 0 else int(v) + 1 for v in a.tolist()])


def _select_out(res: Any) -> Any:
    return res if isinstance(res, Hits) else to_r_index(res)


def _frame_from_dframe(v: Any) -> pl.DataFrame | None:
    if v is None:
        return None
    if not isinstance(v, DFrameOut):
        raise TypeError("value must be a DataFrame or NULL")
    return v.frame


def _strand_arg(x: Any, p: Params) -> dict[str, bool]:
    """ignore.strand: GRanges only; IRanges has no such argument."""
    if not p.has("ignore.strand"):
        return {}
    if not isinstance(x, GRanges):
        raise TypeError("unused argument ignore.strand")
    return {"ignore_strand": scalar_bool(p["ignore.strand"], "ignore.strand")}


def _use_names(x: Any, p: Params) -> Any:
    return x if scalar_bool(p.get("use.names", Vec("lgl", [True]))) else x.with_names(None)


def _check_ranges(x: Any) -> None:
    if not isinstance(x, (IRanges, GRanges)):
        raise TypeError(f"expected IRanges or GRanges, got {type(x).__name__}")


# ---- constructors ----------------------------------------------------------------------


@op("IRanges", "start", "end", "width", "names")
def _iranges(p: Params) -> IRanges:
    return IRanges(
        numbers(p.get("start")), numbers(p.get("end")), numbers(p.get("width")), strings(p.get("names"))
    )


@op("Seqinfo", "seqnames", "seqlengths", "isCircular", "genome")
def _seqinfo_ctor(p: Params) -> Seqinfo:
    def values(key: str) -> list[Any] | None:
        v = p.get(key)
        return None if v is None else _vec(v, key).values

    return Seqinfo(values("seqnames") or [], values("seqlengths"), values("isCircular"), values("genome"))


@op("GRanges", "seqnames", "ranges", "strand", dots=True)
def _granges_ctor(p: Params) -> GRanges:
    dots = dict(p.dots)
    seqinfo = dots.pop("seqinfo", None)
    seqlengths = dots.pop("seqlengths", None)
    sn = p.get("seqnames")
    seqnames = sn if isinstance(sn, Factor) or sn is None else names_arg(sn, "seqnames")
    sl = None
    if seqlengths is not None:
        sv = _vec(seqlengths)
        if sv.names is None:
            raise ValueError("seqlengths must have names")
        sl = {nm: v for nm, v in zip(sv.names, sv.values) if nm is not None}
    cols = {k: vec_to_series(k, _vec(v, k)) for k, v in dots.items()}
    strand = p.get("strand")
    ctor_cols: dict[str, Any] = dict(cols)
    return GRanges(
        seqnames,
        p.get("ranges"),
        None if strand is None else _vec(strand).strings(),
        seqinfo=seqinfo,
        seqlengths=sl,
        **ctor_cols,
    )


# ---- accessors -----------------------------------------------------------------------------


@op("start", "x")
def _start(p: Params) -> Any:
    _check_ranges(p["x"])
    return p["x"].start


@op("end", "x")
def _end(p: Params) -> Any:
    _check_ranges(p["x"])
    return p["x"].end


@op("width", "x")
def _width(p: Params) -> Any:
    _check_ranges(p["x"])
    return p["x"].width


@op("names", "x")
def _names(p: Params) -> Any:
    _check_ranges(p["x"])
    return p["x"].names


@op("length", "x")
def _length(p: Params) -> int:
    return len(p["x"])


@op("names<-", "x", "value")
def _set_names(p: Params) -> Any:
    _check_ranges(p["x"])
    return p["x"].with_names(strings(p["value"]))


@op("mcols", "x")
def _mcols_get(p: Params) -> Any:
    x = p["x"]
    _check_ranges(x)
    frame = x._ranges._mcols if isinstance(x, GRanges) else x._mcols
    if isinstance(x, IRanges) and frame is None:
        return None
    return DFrameOut(frame, len(x))


@op("mcols<-", "x", "value")
def _mcols_set(p: Params) -> Any:
    _check_ranges(p["x"])
    return p["x"].with_mcols(_frame_from_dframe(p["value"]))


@op("seqnames", "x")
def _seqnames(p: Params) -> Any:
    x = p["x"]
    if isinstance(x, Seqinfo):
        return Vec("chr", list(x.seqnames))
    return x.seqnames


@op("strand", "x")
def _strand(p: Params) -> Any:
    return p["x"].strand


@op("ranges", "x")
def _ranges(p: Params) -> Any:
    return p["x"].ranges


@op("seqinfo", "x")
def _seqinfo_get(p: Params) -> Seqinfo:
    return _si(p["x"])


def _si(x: Any) -> Seqinfo:
    return x if isinstance(x, Seqinfo) else x.seqinfo


@op("seqlevels", "x")
def _seqlevels(p: Params) -> Vec:
    return Vec("chr", list(_si(p["x"]).seqnames))


@op("seqlengths", "x")
def _seqlengths(p: Params) -> Vec:
    si = _si(p["x"])
    return Vec("int", list(si.seqlengths), list(si.seqnames))


@op("isCircular", "x")
def _is_circular(p: Params) -> Vec:
    si = _si(p["x"])
    return Vec("lgl", list(si.is_circular), list(si.seqnames))


@op("genome", "x")
def _genome(p: Params) -> Vec:
    si = _si(p["x"])
    return Vec("chr", list(si.genome), list(si.seqnames))


def _field_value(si: Seqinfo, value: Any, what: str) -> list[Any]:
    vec = _vec(value, "value")
    if vec.names is not None and list(vec.names) != list(si.seqnames):
        raise ValueError(f"the names on the supplied '{what}' vector must be identical to the seqnames")
    return list(vec.values)


@op("seqlengths<-", "x", "value")
def _set_seqlengths(p: Params) -> Any:
    x = p["x"]
    vals = _field_value(_si(x), p["value"], "seqlengths")
    return _si(x).with_field("seqlengths", vals) if isinstance(x, Seqinfo) else x.with_seqlengths(vals)


@op("isCircular<-", "x", "value")
def _set_circular(p: Params) -> Any:
    x = p["x"]
    vals = _field_value(_si(x), p["value"], "isCircular")
    return _si(x).with_field("is_circular", vals) if isinstance(x, Seqinfo) else x.with_is_circular(vals)


@op("genome<-", "x", "value")
def _set_genome(p: Params) -> Any:
    x = p["x"]
    vals = _field_value(_si(x), p["value"], "genome")
    return _si(x).with_field("genome", vals) if isinstance(x, Seqinfo) else x.with_genome(vals)


@op("seqlevels<-", "x", "value")
def _set_seqlevels(p: Params) -> Any:
    x = p["x"]
    value = names_arg(p["value"], "value")
    if isinstance(x, Seqinfo):
        return x.rename_levels(value)
    return x.with_seqlevels(value)


@op("seqinfo<-", "x", "value")
def _set_seqinfo(p: Params) -> Any:
    if not isinstance(p["value"], Seqinfo):
        raise TypeError("value must be a Seqinfo object")
    return p["x"].with_seqinfo(p["value"])


@op("merge", "x", "y")
def _merge(p: Params) -> Seqinfo:
    merged: Seqinfo = p["x"].merge(p["y"])
    return merged


@op("[", "x", "i")
def _subset(p: Params) -> Any:
    x, i = p["x"], p["i"]
    if isinstance(x, Seqinfo):
        vec = _vec(i, "i")
        if vec.kind != "chr":
            raise TypeError("a Seqinfo object can be subsetted only by name")
        return x.subset(names_arg(vec, "i"))
    _check_ranges(x)
    return x[r_subscript(i, len(x), x.names)]


@op("c", varargs=True)
def _concat(p: Params) -> Any:
    if not p.rest:
        raise TypeError("nothing to combine")
    x, *others = p.rest
    _check_ranges(x)
    return x.concat(*others)


@op("trim", "x", "use.names")
def _trim(p: Params) -> Any:
    return _use_names(p["x"].trim(), p)


# ---- intra-range ---------------------------------------------------------------------------------


@op("shift", "x", "shift", "use.names")
def _shift(p: Params) -> Any:
    x = p["x"]
    _check_ranges(x)
    return _use_names(x.shift(numbers(p.get("shift", Vec("int", [0])))), p)


@op("narrow", "x", "start", "end", "width", "use.names")
def _narrow(p: Params) -> Any:
    x = p["x"]
    _check_ranges(x)
    res = x.narrow(numbers(p.get("start")), numbers(p.get("end")), numbers(p.get("width")))
    return _use_names(res, p)


@op("resize", "x", "width", "fix", "use.names", "ignore.strand")
def _resize(p: Params) -> Any:
    x = p["x"]
    _check_ranges(x)
    fix = strings(p["fix"]) if p.has("fix") else "start"
    if not p.has("width"):
        raise TypeError('argument "width" is missing')
    if isinstance(fix, list) and any(f is None for f in fix):
        raise ValueError("fix cannot contain NAs")
    return x.resize(
        numbers(p["width"]), fix, scalar_bool(p.get("use.names", Vec("lgl", [True]))), **_strand_arg(x, p)
    )


@op("flank", "x", "width", "start", "both", "use.names", "ignore.strand")
def _flank(p: Params) -> Any:
    x = p["x"]
    _check_ranges(x)
    start = bool_array(p["start"]) if p.has("start") else True
    res = x.flank(
        numbers(p["width"]),
        start,
        scalar_bool(p.get("both", Vec("lgl", [False]))),
        **_strand_arg(x, p),
    )
    return _use_names(res, p)


@op("promoters", "x", "upstream", "downstream", "use.names")
def _promoters(p: Params) -> Any:
    x = p["x"]
    _check_ranges(x)
    res = x.promoters(
        numbers(p.get("upstream", Vec("int", [2000]))), numbers(p.get("downstream", Vec("int", [200])))
    )
    return _use_names(res, p)


# ---- inter-range -----------------------------------------------------------------------------------


@op("range", "x", varargs=True, dots=True)
def _range(p: Params) -> Any:
    x = p["x"]
    _check_ranges(x)
    if set(p.dots) - {"ignore.strand"}:
        raise TypeError("unused argument")
    kw = _strand_arg(x, Params({"ignore.strand": p.dots["ignore.strand"]}, [], {})) if p.dots else {}
    if p.rest:
        if isinstance(x, GRanges):
            return x.concat(*p.rest).range(**kw)
        return x.range(*p.rest)
    return x.range(**kw)


@op("reduce", "x", "drop.empty.ranges", "min.gapwidth", "with.revmap", "with.inframe.attrib", "ignore.strand")
def _reduce(p: Params) -> Any:
    x = p["x"]
    _check_ranges(x)
    if p.has("with.inframe.attrib"):
        raise TypeError("with.inframe.attrib is not supported")
    return x.reduce(
        scalar_bool(p.get("drop.empty.ranges", Vec("lgl", [False]))),
        scalar_int(p.get("min.gapwidth", Vec("int", [1]))),
        scalar_bool(p.get("with.revmap", Vec("lgl", [False]))),
        **_strand_arg(x, p),
    )


def _per_level(v: Any, what: str) -> Any:
    """gaps start/end for GRanges: scalar, named (by seqlevel) or positional vector; NA -> None."""
    if v is None:
        return None
    vec = _vec(v, what)
    vals = [None if e is None else int(e) for e in vec.values]
    if vec.names is not None:
        return dict(zip(vec.names, vals))
    return vals[0] if len(vals) == 1 else vals


@op("gaps", "x", "start", "end", "ignore.strand")
def _gaps(p: Params) -> Any:
    x = p["x"]
    _check_ranges(x)
    if isinstance(x, GRanges):
        start = _per_level(p.get("start"), "start")
        end = _per_level(p.get("end"), "end")
        return x.gaps(start, end, **_strand_arg(x, p))
    return x.gaps(opt_int(p.get("start")), opt_int(p.get("end")))


@op("disjoin", "x", "with.revmap", "ignore.strand")
def _disjoin(p: Params) -> Any:
    x = p["x"]
    _check_ranges(x)
    return x.disjoin(scalar_bool(p.get("with.revmap", Vec("lgl", [False]))), **_strand_arg(x, p))


# ---- set operations ----------------------------------------------------------------------------------


@op("union", "x", "y", "ignore.strand")
def _union(p: Params) -> Any:
    x = p["x"]
    _check_ranges(x)
    return x.union(p["y"], **_strand_arg(x, p))


@op("intersect", "x", "y", "ignore.strand")
def _intersect(p: Params) -> Any:
    x = p["x"]
    _check_ranges(x)
    return x.intersect(p["y"], **_strand_arg(x, p))


@op("setdiff", "x", "y", "ignore.strand")
def _setdiff(p: Params) -> Any:
    x = p["x"]
    _check_ranges(x)
    return x.setdiff(p["y"], **_strand_arg(x, p))


@op("punion", "x", "y", "fill.gap", "ignore.strand")
def _punion(p: Params) -> Any:
    x = p["x"]
    _check_ranges(x)
    return x.punion(p["y"], scalar_bool(p.get("fill.gap", Vec("lgl", [False]))), **_strand_arg(x, p))


@op("pintersect", "x", "y", "resolve.empty", "drop.nohit.ranges", "ignore.strand")
def _pintersect(p: Params) -> Any:
    x = p["x"]
    _check_ranges(x)
    if isinstance(x, GRanges):
        if p.has("resolve.empty"):
            raise TypeError("unused argument resolve.empty")
        return x.pintersect(
            p["y"], scalar_bool(p.get("drop.nohit.ranges", Vec("lgl", [False]))), **_strand_arg(x, p)
        )
    if p.has("drop.nohit.ranges"):
        raise TypeError("unused argument drop.nohit.ranges")
    return x.pintersect(p["y"], scalar_str(p["resolve.empty"]) if p.has("resolve.empty") else "none")


@op("psetdiff", "x", "y", "ignore.strand")
def _psetdiff(p: Params) -> Any:
    x = p["x"]
    _check_ranges(x)
    return x.psetdiff(p["y"], **_strand_arg(x, p))


# ---- overlaps ------------------------------------------------------------------------------------------


def _overlap_kwargs(p: Params, x: Any) -> dict[str, Any]:
    kw: dict[str, Any] = {
        "maxgap": scalar_int(p.get("maxgap", Vec("int", [-1])), "maxgap"),
        "minoverlap": scalar_int(p.get("minoverlap", Vec("int", [0])), "minoverlap"),
        "type": scalar_str(p.get("type", Vec("chr", ["any"])), "type"),
    }
    kw.update(_strand_arg(x, p))
    return kw


@op("findOverlaps", "query", "subject", "maxgap", "minoverlap", "type", "select", "ignore.strand")
def _find_overlaps(p: Params) -> Any:
    q = p["query"]
    _check_ranges(q)
    select = scalar_str(p["select"]) if p.has("select") else "all"
    return _select_out(q.find_overlaps(p["subject"], select=select, **_overlap_kwargs(p, q)))


@op("countOverlaps", "query", "subject", "maxgap", "minoverlap", "type", "ignore.strand")
def _count_overlaps(p: Params) -> Any:
    q = p["query"]
    _check_ranges(q)
    counts = q.count_overlaps(p["subject"], **_overlap_kwargs(p, q))
    return Vec("int", counts.tolist(), q.names)


@op("overlapsAny", "query", "subject", "maxgap", "minoverlap", "type", "ignore.strand")
def _overlaps_any(p: Params) -> Any:
    q = p["query"]
    _check_ranges(q)
    return q.overlaps_any(p["subject"], **_overlap_kwargs(p, q))


@op("subsetByOverlaps", "x", "ranges", "maxgap", "minoverlap", "type", "invert", "ignore.strand")
def _subset_by_overlaps(p: Params) -> Any:
    x = p["x"]
    _check_ranges(x)
    return x.subset_by_overlaps(
        p["ranges"], invert=scalar_bool(p.get("invert", Vec("lgl", [False]))), **_overlap_kwargs(p, x)
    )


# ---- nearest -----------------------------------------------------------------------------------------------


def _neighbour(name: str, default_select: str) -> None:
    @op(name, "x", "subject", "select", "ignore.strand")
    def _impl(p: Params) -> Any:
        x = p["x"]
        _check_ranges(x)
        select = scalar_str(p["select"]) if p.has("select") else default_select
        fn = getattr(x, name)
        return _select_out(fn(p.get("subject"), select=select, **_strand_arg(x, p)))


_neighbour("nearest", "arbitrary")
_neighbour("precede", "first")
_neighbour("follow", "last")


@op("distance", "x", "y", "ignore.strand")
def _distance(p: Params) -> Any:
    x = p["x"]
    _check_ranges(x)
    d = x.distance(p["y"], **_strand_arg(x, p))
    if isinstance(x, GRanges):
        return Vec("int", [None if v < 0 else int(v) for v in d.tolist()])
    return d.astype(np.int32)


# ---- coverage ---------------------------------------------------------------------------------------------------


@op("coverage", "x", "shift", "width", "weight")
def _coverage(p: Params) -> Any:
    x = p["x"]
    _check_ranges(x)
    shift = numbers(p.get("shift", Vec("int", [0])))
    weight_arg = p.get("weight", Vec("int", [1]))
    if isinstance(x, GRanges):
        weight: Any = scalar_str(weight_arg) if _vec(weight_arg).kind == "chr" else numbers(weight_arg)
        wv = p.get("width")
        width: Any = None
        if wv is not None:
            vec = _vec(wv, "width")
            if vec.names is not None:
                if list(vec.names) != list(x.seqinfo.seqnames):
                    raise ValueError("when 'width' has names, they must be identical to the seqlevels")
                width = dict(zip(vec.names, vec.ints().tolist()))
            else:
                width = int(vec.ints()[0]) if len(vec) == 1 else vec.ints().tolist()
        return x.coverage(shift, width, weight)
    width = opt_int(p.get("width"))
    return x.coverage(shift, width, numbers(weight_arg))


# ---- ordering --------------------------------------------------------------------------------------------------------


@op("order", "x", varargs=True, dots=True)
def _order(p: Params) -> Any:
    x = p["x"]
    _check_ranges(x)
    from .._ordering import order_by

    decreasing = scalar_bool(p.dots["decreasing"]) if "decreasing" in p.dots else False
    if set(p.dots) - {"decreasing"}:
        raise TypeError("unused argument")
    if p.rest:
        keys = x._keys() + [k for y in p.rest for k in y._keys()]
        return order_by(keys, decreasing) + 1
    return x.order(decreasing) + 1


@op("sort", "x", "decreasing", "ignore.strand")
def _sort(p: Params) -> Any:
    x = p["x"]
    _check_ranges(x)
    return x.sort(scalar_bool(p.get("decreasing", Vec("lgl", [False]))), **_strand_arg(x, p))


@op("rank", "x", "ties.method", "ignore.strand")
def _rank(p: Params) -> Any:
    x = p["x"]
    _check_ranges(x)
    method = scalar_str(p["ties.method"]) if p.has("ties.method") else "average"
    return x.rank(method, **_strand_arg(x, p))


@op("duplicated", "x", "fromLast")
def _duplicated(p: Params) -> Any:
    x = p["x"]
    _check_ranges(x)
    return x.duplicated(scalar_bool(p.get("fromLast", Vec("lgl", [False]))))


@op("unique", "x", "fromLast")
def _unique(p: Params) -> Any:
    x = p["x"]
    _check_ranges(x)
    return x.unique(scalar_bool(p.get("fromLast", Vec("lgl", [False]))))
