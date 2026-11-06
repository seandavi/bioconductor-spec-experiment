"""Adapter for port e0 (`ranges_e0`, the E0 control: IRanges / GRanges in plain Python).

Port docs used: experiments/e0-control/README.md (API map, conventions), docstrings and public names in
src/ranges_e0/{__init__,iranges,granges,seqinfo,hits,rle}.py. Nothing else of the port was read.

Port index base: 0-based for every index-valued result (Hits, order, nearest/precede/follow and
findOverlaps `select=` results, `revmap`), with `NO_MATCH == -1` for "no match". `rank` is 1-based.
Subscripts: negative integers follow Python (not R), so only positive R subscripts are passed.

Translations (details and unsupported ops in e0.md):
- R name -> port name: `.` -> `_` (done by the base), `fromLast` -> `from_last`, second operand `y` -> `other`,
  `c(x, ...)` -> `concat`, `names<-` -> `with_names`, `mcols<-` -> `with_mcols`, `seqinfo<-` -> `with_seqinfo`.
- Atomic R vectors (`RVector`) become Python scalars (length 1), numpy arrays (int/dbl) or lists (lgl/chr).
  NA in a numeric vector becomes NaN (the port rejects it), NA in lgl/chr becomes None.
- Atomic port results (ndarray, list, dict) are wrapped back into `RVector` so nested calls see R vectors.
- Index-valued results are shifted +1; the port's `NO_MATCH` marker becomes NA (also for GRanges `distance`,
  which the port documents as using `NO_MATCH`). The `revmap` mcols column is shifted +1.
- Port mcols are `polars.DataFrame` (zero columns = none): Factor <-> `pl.Enum`, IntegerList <-> `pl.List(Int64)`.
- A zero-column mcols is returned as R `NULL` by `mcols(x)`.
- An argument the port's method does not take -> `UnknownArgument`; a call the port's signature cannot bind
  (missing required argument) or an attribute the port class lacks -> `Unsupported`.
"""

from __future__ import annotations

import inspect
import math
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
)
from ranges_e0 import NO_MATCH, GRanges, Hits, IRanges, Rle, Seqinfo, concat

PORT_OBJECTS = (IRanges, GRanges, Seqinfo, Hits, Rle)

# R argument name -> port argument name, where `.` -> `_` is not enough
RENAME = {"fromLast": "from_last", "y": "other"}

PL_DTYPE = {"int": pl.Int64, "dbl": pl.Float64, "lgl": pl.Boolean, "chr": pl.String}
R_TYPE_OF_NUMPY_KIND = {"i": "int", "u": "int", "f": "dbl", "b": "lgl", "O": "chr", "U": "chr", "S": "chr"}


# ---------------------------------------------------------------- value conversion

def num_array(rv: RVector) -> np.ndarray:
    """int/dbl/lgl vector -> ndarray; NA becomes NaN in a float array (the port rejects it)."""
    vals = rv.list()
    if any(v is None for v in vals):
        return np.array([math.nan if v is None else v for v in vals], dtype=np.float64)
    return rv.np()


def to_arg(v: Any) -> Any:
    """An R argument value -> what the port's methods take."""
    if not isinstance(v, RVector):
        return v
    if len(v) == 1:
        x = v.scalar()
        if x is None and v.type in ("int", "dbl"):
            return math.nan
        return x
    if v.type in ("int", "dbl"):
        return num_array(v)
    return v.list()


def to_rvector(arr: np.ndarray) -> RVector:
    """Atomic port result -> R vector, typed by dtype."""
    return RVector(R_TYPE_OF_NUMPY_KIND[arr.dtype.kind], arr.tolist(), None)


def named_rvector(type_: str, mapping: dict) -> RVector:
    return RVector(type_, list(mapping.values()), list(mapping.keys()))


def to_series(name: str, rv: RVector) -> pl.Series:
    return pl.Series(name, rv.list(), dtype=PL_DTYPE[rv.type])


def name_mapping(rv: RVector) -> dict:
    """Named R vector -> the port's {name: value} mapping; a mapping cannot hold repeated or NA names."""
    if len(set(rv.names)) != len(rv.names) or any(n is None for n in rv.names):
        raise Unsupported("named vector with repeated or NA names")
    return dict(zip(rv.names, rv.list()))


def mcols_dict(extra_named: list) -> dict | None:
    """`...` of IRanges() / GRanges() -> the port's `mcols=` dict of columns."""
    if not extra_named:
        return None
    names = [k for k, _ in extra_named]
    if len(set(names)) != len(names):
        raise Unsupported("mcols: duplicate column names")
    out = {}
    for k, v in extra_named:
        if not isinstance(v, RVector):
            raise Unsupported("mcols: non-atomic column")
        out[k] = to_series(k, v)
    return out


def column_series(name: str, col: dict) -> pl.Series:
    """PROJECTION mcols / DataFrame column -> polars Series."""
    kind = col["kind"]
    if kind == "vector":
        return to_series(name, RVector(col["type"], col["values"]))
    if kind == "Factor":
        levels = col["levels"]
        return pl.Series(name, [None if c is None else levels[c - 1] for c in col["codes"]], dtype=pl.Enum(levels))
    if kind == "IntegerList":
        return pl.Series(name, col["elements"], dtype=pl.List(pl.Int64))
    raise Unsupported(f"input:column:{kind}")


def mcols_frame(columns: list) -> pl.DataFrame | None:
    return pl.DataFrame([column_series(n, c) for n, c in columns]) if columns else None


def project_series(s: pl.Series) -> dict:
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
    raise ValueError(f"cannot project polars dtype {dt}")


def project_columns(df: pl.DataFrame) -> list:
    return [(c, project_series(df[c])) for c in df.columns]


def project_seqinfo(si: Seqinfo) -> dict:
    return p_seqinfo(si.seqnames, list(si.seqlengths.values()), list(si.is_circular.values()), list(si.genome.values()))


def shift_revmap_to_r_base(result: Any) -> Any:
    """`revmap` holds 0-based input indices in the port; R's are 1-based."""
    if "revmap" not in result.mcols.columns:
        return result
    shifted = result.mcols["revmap"].list.eval(pl.element() + 1)
    return result.with_mcols(result.mcols.with_columns(shifted.alias("revmap")))


# ---------------------------------------------------------------- the adapter

class Adapter(AdapterBase):
    port = "e0"
    port_packages = ("ranges_e0",)
    unsupported_ops = {
        "seqlevels<-": "no public setter: seqinfo<- renaming/pruning modes are not implemented",
        "seqlengths<-": "no public setter for a single Seqinfo column",
        "isCircular<-": "no public setter for a single Seqinfo column",
        "genome<-": "no public setter for a single Seqinfo column",
    }

    # ---- inputs ----------------------------------------------------------------
    def build_IRanges(self, p: dict) -> IRanges:
        self._require_representable_names(p["names"])
        return IRanges(start=self._ints(p["start"]), width=self._ints(p["width"]), names=p["names"],
                       mcols=mcols_frame(p["mcols"]))

    def build_GRanges(self, p: dict) -> GRanges:
        self._require_representable_names(p["names"])
        ranges = IRanges(start=self._ints(p["start"]), width=self._ints(p["width"]), names=p["names"])
        return GRanges(p["seqnames"], ranges, p["strand"], mcols=mcols_frame(p["mcols"]),
                       seqinfo=self.build_Seqinfo(p["seqinfo"]))

    def build_Seqinfo(self, p: dict) -> Seqinfo:
        return Seqinfo(p["seqnames"], p["seqlengths"], p["isCircular"], p["genome"])

    def build_DataFrame(self, p: dict) -> pl.DataFrame:
        if p["rownames"] is not None:
            raise Unsupported("input:DataFrame rownames")
        if not p["columns"] and p["nrow"] > 0:
            raise Unsupported("input:DataFrame zero columns with rows")
        return mcols_frame(p["columns"]) if p["columns"] else pl.DataFrame()

    def build_Factor(self, p: dict) -> pl.Series:
        return column_series("", p)

    @staticmethod
    def _require_representable_names(names: list | None) -> None:
        """The port stores an NA name as ""."""
        if names is not None and any(n is None for n in names):
            raise Unsupported("input:names with NA")

    @staticmethod
    def _ints(values: list) -> np.ndarray:
        """Integer list -> ndarray; NA becomes NaN in a float array (the port rejects it)."""
        if any(v is None for v in values):
            return np.array([math.nan if v is None else v for v in values], dtype=np.float64)
        return np.array(values, dtype=np.int64)

    # ---- calling the port ---------------------------------------------------------
    def _method(self, name: str, obj: Any, *args: Any, **kw: Any) -> Any:
        """obj.name(*args, **kw) with R argument values converted; kw keys are R names (`.` -> `_`)."""
        if not isinstance(obj, PORT_OBJECTS):
            raise Unsupported(f"{name}: not a port object")
        fn = getattr(obj, name, None)
        if fn is None:
            raise Unsupported(f"{type(obj).__name__}.{name}")
        sig = inspect.signature(fn)
        call_kw = {}
        for k, v in kw.items():
            pk = RENAME.get(k, k)
            if pk not in sig.parameters:
                raise UnknownArgument(k)
            call_kw[pk] = to_arg(v)
        try:
            sig.bind(*args, **call_kw)
        except TypeError as e:
            raise Unsupported(f"{type(obj).__name__}.{name}: {e}")
        return fn(*args, **call_kw)

    def _attr(self, obj: Any, name: str) -> Any:
        if not isinstance(obj, PORT_OBJECTS):
            raise Unsupported(f"{name}: not a port object")
        if not hasattr(obj, name):
            raise Unsupported(f"{type(obj).__name__}.{name}")
        return getattr(obj, name)

    def _construct(self, cls: type, **kw: Any) -> Any:
        sig = inspect.signature(cls)
        try:
            sig.bind(**kw)
        except TypeError as e:
            raise Unsupported(f"{cls.__name__}(): {e}")
        return cls(**kw)

    @staticmethod
    def _index_vector(arr: np.ndarray, no_match_is_na: bool = True) -> RVector:
        """0-based index array -> 1-based R integer vector; the port's NO_MATCH marker becomes NA."""
        return RVector("int", [None if (no_match_is_na and v == NO_MATCH) else v + 1 for v in arr.tolist()], None)

    def _index_result(self, res: Any) -> Any:
        """Hits stay native (shifted when projected); `select=` results are 0-based index arrays."""
        return res if isinstance(res, Hits) else self._index_vector(res)

    # ---- constructors -----------------------------------------------------------------
    def op_IRanges(self, *, start=MISSING, end=MISSING, width=MISSING, names=MISSING, _mcols=()):
        kw = {k: to_arg(v) for k, v in (("start", start), ("end", end), ("width", width)) if v is not MISSING}
        if names is not MISSING:
            kw["names"] = names.list() if isinstance(names, RVector) else names
        mcols = mcols_dict(list(_mcols))
        if mcols is not None:
            kw["mcols"] = mcols
        return self._construct(IRanges, **kw)

    def op_GRanges(self, *, seqnames=MISSING, ranges=MISSING, strand=MISSING, seqinfo=MISSING,
                   seqlengths=MISSING, _mcols=()):
        kw = {k: to_arg(v) for k, v in (("ranges", ranges), ("seqinfo", seqinfo)) if v is not MISSING}
        # atomic seqnames / strand go as lists, so that NA stays distinct from Python None (R NULL)
        for k, v in (("seqnames", seqnames), ("strand", strand)):
            if v is not MISSING:
                kw[k] = v.list() if isinstance(v, RVector) else v
        if seqlengths is not MISSING:
            # the port takes {seqname: length}; an unnamed vector is handed over as is and rejected by the port
            if isinstance(seqlengths, RVector) and seqlengths.names is not None:
                kw["seqlengths"] = name_mapping(seqlengths)
            else:
                kw["seqlengths"] = to_arg(seqlengths)
        mcols = mcols_dict(list(_mcols))
        if mcols is not None:
            kw["mcols"] = mcols
        return self._construct(GRanges, **kw)

    def op_Seqinfo(self, *, seqnames=MISSING, seqlengths=MISSING, isCircular=MISSING, genome=MISSING):
        kw = {}
        for pk, v in (("seqnames", seqnames), ("seqlengths", seqlengths), ("is_circular", isCircular)):
            if v is not MISSING:
                kw[pk] = v.list() if isinstance(v, RVector) else v
        if genome is not MISSING:
            # the port takes a str for a single genome shared by every sequence, else one entry per sequence
            kw["genome"] = to_arg(genome) if isinstance(genome, RVector) and len(genome) == 1 else (
                genome.list() if isinstance(genome, RVector) else genome)
        return self._construct(Seqinfo, **kw)

    # ---- accessors ------------------------------------------------------------------------
    def op_start(self, *, x):
        return to_rvector(self._attr(x, "start"))

    def op_end(self, *, x):
        return to_rvector(self._attr(x, "end"))

    def op_width(self, *, x):
        return to_rvector(self._attr(x, "width"))

    def op_names(self, *, x):
        names = self._attr(x, "names")
        return None if names is None else RVector("chr", names, None)

    def op_length(self, *, x):
        if not isinstance(x, PORT_OBJECTS):
            raise Unsupported("length: not a port object")
        return RVector("int", [len(x)], None)

    def op_seqnames(self, *, x):
        names = self._attr(x, "seqnames")
        # Seqinfo.seqnames is a tuple of level names, GRanges.seqnames a per-range array
        return RVector("chr", list(names), None) if isinstance(names, tuple) else to_rvector(names)

    def op_strand(self, *, x):
        return to_rvector(self._attr(x, "strand"))

    def op_ranges(self, *, x):
        return self._attr(x, "ranges")

    def op_mcols(self, *, x):
        df = self._attr(x, "mcols")
        return df if df.width else None  # the port represents "no metadata" as a zero-column frame

    # ---- metadata ------------------------------------------------------------------------------
    def op_set_names(self, *, x, value):
        return self._method("with_names", x, value.list() if isinstance(value, RVector) else value)

    def op_set_mcols(self, *, x, value):
        if isinstance(value, RVector):
            raise Unsupported("mcols<-: atomic value")
        return self._method("with_mcols", x, value)

    # ---- seqinfo -----------------------------------------------------------------------------------
    def op_seqinfo(self, *, x):
        return self._attr(x, "seqinfo")

    def op_seqlevels(self, *, x):
        return RVector("chr", list(self._attr(x, "seqlevels")), None)

    def op_seqlengths(self, *, x):
        return named_rvector("int", self._attr(x, "seqlengths"))

    def op_isCircular(self, *, x):
        return named_rvector("lgl", self._attr(x, "is_circular"))

    def op_genome(self, *, x):
        return named_rvector("chr", self._attr(x, "genome"))

    def op_set_seqinfo(self, *, x, value):
        return self._method("with_seqinfo", x, value)

    def op_merge(self, *, x, y):
        return self._method("merge", x, y)

    def op_trim(self, *, x):
        return self._method("trim", x)

    # ---- subset / concat ---------------------------------------------------------------------------------
    def op_subset(self, *, x, i=MISSING):
        if not isinstance(x, PORT_OBJECTS) or not hasattr(x, "__getitem__"):
            raise Unsupported("[: not a subsettable port object")
        if i is MISSING:
            raise Unsupported("subscript:missing")
        return x[self._subscript(x, i)]

    @staticmethod
    def _subscript(x: Any, i: Any) -> Any:
        """R subscript -> port subscript. The port documents Python semantics for negative integers and
        no recycling for logicals, so only the forms whose meaning is the same are passed."""
        if not isinstance(i, RVector):
            raise Unsupported(f"subscript:{type(i).__name__}")
        vals = i.list()
        if any(v is None for v in vals):
            raise Unsupported("subscript:NA")
        if i.type in ("int", "dbl"):
            if any(v != int(v) or v <= 0 for v in vals):
                raise Unsupported("subscript:non-positive or fractional")
            return np.array([int(v) - 1 for v in vals], dtype=np.int64)  # R is 1-based, the port 0-based
        if i.type == "lgl":
            if len(vals) != len(x):
                raise Unsupported("subscript:logical recycling")
            return np.array(vals, dtype=bool)
        return vals

    def op_c(self, *dots, x):
        if not isinstance(x, PORT_OBJECTS):
            raise Unsupported("c: not a port object")
        return concat(x, *dots)

    # ---- intra-range -----------------------------------------------------------------------------------------
    def op_shift(self, *, x, **kw):
        return self._method("shift", x, **kw)

    def op_narrow(self, *, x, **kw):
        return self._method("narrow", x, **kw)

    def op_resize(self, *, x, **kw):
        return self._method("resize", x, **kw)

    def op_flank(self, *, x, **kw):
        return self._method("flank", x, **kw)

    def op_promoters(self, *, x, **kw):
        return self._method("promoters", x, **kw)

    # ---- inter-range ------------------------------------------------------------------------------------------------
    def op_range(self, *dots, x, **kw):
        return self._method("range", x, *dots, **kw)

    def op_reduce(self, *, x, **kw):
        return shift_revmap_to_r_base(self._method("reduce", x, **kw))

    def op_gaps(self, *, x, **kw):
        return self._method("gaps", x, **kw)

    def op_disjoin(self, *, x, **kw):
        return shift_revmap_to_r_base(self._method("disjoin", x, **kw))

    # ---- set operations ------------------------------------------------------------------------------------------------
    def op_union(self, *, x, y, **kw):
        return self._method("union", x, y, **kw)

    def op_intersect(self, *, x, y, **kw):
        return self._method("intersect", x, y, **kw)

    def op_setdiff(self, *, x, y, **kw):
        return self._method("setdiff", x, y, **kw)

    def op_punion(self, *, x, y, **kw):
        return self._method("punion", x, y, **kw)

    def op_pintersect(self, *, x, y, **kw):
        return self._method("pintersect", x, y, **kw)

    def op_psetdiff(self, *, x, y, **kw):
        return self._method("psetdiff", x, y, **kw)

    # ---- overlaps ----------------------------------------------------------------------------------------------------------
    def op_findOverlaps(self, *, query, **kw):
        return self._index_result(self._method("find_overlaps", query, **kw))

    def op_countOverlaps(self, *, query, **kw):
        return to_rvector(self._method("count_overlaps", query, **kw))

    def op_overlapsAny(self, *, query, **kw):
        return to_rvector(self._method("overlaps_any", query, **kw))

    def op_subsetByOverlaps(self, *, x, **kw):
        return self._method("subset_by_overlaps", x, **kw)

    # ---- nearest -------------------------------------------------------------------------------------------------------------------
    def op_nearest(self, *, x, **kw):
        return self._index_result(self._method("nearest", x, **kw))

    def op_precede(self, *, x, **kw):
        return self._index_result(self._method("precede", x, **kw))

    def op_follow(self, *, x, **kw):
        return self._index_result(self._method("follow", x, **kw))

    def op_distance(self, *, x, y, **kw):
        # the port marks incomparable GRanges pairs with NO_MATCH where R has NA
        d = self._method("distance", x, y, **kw)
        return RVector("int", [None if v == NO_MATCH else v for v in d.tolist()], None)

    # ---- coverage ------------------------------------------------------------------------------------------------------------------------
    def op_coverage(self, *, x, **kw):
        width = kw.get("width")
        if isinstance(width, RVector) and (len(width) > 1 or width.names is not None):
            if width.names is None:
                raise Unsupported("argument:width unnamed vector")
            kw["width"] = name_mapping(width)  # the port takes {seqname: length}, also for one name
        return self._method("coverage", x, **kw)

    # ---- ordering -------------------------------------------------------------------------------------------------------------------------------
    def op_order(self, *dots, **kw):
        if len(dots) != 1:
            raise Unsupported("order: the port orders exactly one object")
        return self._index_vector(self._method("order", dots[0], **kw), no_match_is_na=False)

    def op_sort(self, *, x, **kw):
        return self._method("sort", x, **kw)

    def op_rank(self, *, x, **kw):
        return to_rvector(self._method("rank", x, **kw))

    def op_duplicated(self, *, x, **kw):
        return to_rvector(self._method("duplicated", x, **kw))

    def op_unique(self, *, x, **kw):
        return self._method("unique", x, **kw)

    # ---- projection ---------------------------------------------------------------------------------------------------------------------------------------
    def project(self, value: Any) -> dict:
        if isinstance(value, IRanges):
            return p_iranges(value.start, value.width, value.names, project_columns(value.mcols))
        if isinstance(value, GRanges):
            return p_granges(value.seqnames, value.start, value.width, value.strand, value.names,
                             project_columns(value.mcols), project_seqinfo(value.seqinfo))
        if isinstance(value, Seqinfo):
            return project_seqinfo(value)
        if isinstance(value, pl.DataFrame):
            return p_dataframe(value.height, project_columns(value))
        if isinstance(value, pl.Series):
            return project_series(value)
        if isinstance(value, Hits):
            return p_hits(value.query_hits + 1, value.subject_hits + 1, value.query_length, value.subject_length)
        if isinstance(value, Rle):
            return p_rle(value.values, value.lengths)
        if isinstance(value, dict) and all(isinstance(r, Rle) for r in value.values()):
            return p_rlelist(list(value.keys()), [(r.values, r.lengths) for r in value.values()])
        raise ValueError(f"cannot project {type(value).__name__}")
