"""R operation name -> Python callable. Explicit on purpose (PLAN 4.4).

Naming: R camelCase -> snake_case methods; R ``.`` in argument names -> ``_``.
``index`` marks results that are 0-based positions in Python (projected +1);
``dec_args`` marks positional int inputs that are R 1-based subscripts (built -1).
"""

from __future__ import annotations

import re
from collections.abc import Callable
from dataclasses import dataclass
from typing import Any

import numpy as np


@dataclass(frozen=True)
class Op:
    fn: Callable[..., Any]
    index: Callable[[dict[str, Any]], bool] = lambda kw: False
    dec_args: tuple[int, ...] = ()
    kind: str | None = None  # expected result kind when it cannot be inferred from the value (empty dicts)


def _prop(name: str) -> Op:
    return Op(lambda x: getattr(x, name))


def _method(name: str, index: Callable[[dict[str, Any]], bool] = lambda kw: False) -> Op:
    return Op(lambda x, *a, **kw: getattr(x, name)(*a, **kw), index)


def _select_index(default: str) -> Callable[[dict[str, Any]], bool]:
    return lambda kw: kw.get("select", default) != "all"


OPS: dict[str, Op] = {
    "identity": Op(lambda x: x),
    "length": Op(len),
    "[": Op(lambda x, i: x[i], dec_args=(1,)),
    "c": Op(lambda x, *rest: type(x).concat(x, *rest)),
    # accessors
    **{r: _prop(p) for r, p in {
        "start": "start", "end": "end", "width": "width", "names": "names",
        "seqnames": "seqnames", "strand": "strand", "ranges": "ranges",
        "seqlevels": "seqlevels", "seqlengths": "seqlengths",
        "isCircular": "is_circular", "genome": "genome", "seqinfo": "seqinfo",
    }.items()},
    # intra-range, inter-range, set ops, ordering
    **{r: _method(r) for r in (
        "shift", "narrow", "resize", "flank", "promoters", "range", "reduce",
        "gaps", "disjoin", "union", "intersect", "setdiff", "punion",
        "pintersect", "psetdiff", "unique", "duplicated", "sort", "rank",
        "distance", "trim",
    )},
    "coverage": Op(lambda x, *a, **kw: x.coverage(*a, **kw), kind="coverage"),
    "order": _method("order", lambda kw: True),
    "findOverlaps": _method("find_overlaps", _select_index("all")),
    "countOverlaps": _method("count_overlaps"),
    "overlapsAny": _method("overlaps_any"),
    "subsetByOverlaps": _method("subset_by_overlaps"),
    "nearest": _method("nearest", _select_index("arbitrary")),
    "precede": _method("precede", _select_index("first")),
    "follow": _method("follow", _select_index("first")),
}


def _snake(r_name: str) -> str:
    """R argument name -> Python: `.` -> `_`, camelCase -> snake_case."""
    return re.sub(r"(?<=[a-z])([A-Z])", lambda m: "_" + m.group(1).lower(), r_name.replace(".", "_"))


def result_kind(call: dict[str, Any]) -> str | None:
    return OPS[call["op"]].kind


def run_op(call: dict[str, Any], inputs: dict[str, Any]) -> tuple[Any, bool]:
    """Run a case's call; return (result, result_is_index_vector)."""
    op = OPS[call["op"]]
    args = [inputs[n] for n in call.get("args", [])]
    for pos in op.dec_args:
        if isinstance(args[pos], np.ndarray) and args[pos].dtype.kind == "i":
            args[pos] = args[pos] - 1
    kwargs: dict[str, Any] = {}
    for k, v in (call.get("kwargs") or {}).items():
        if isinstance(v, dict) and set(v) == {"ref"}:
            v = inputs[v["ref"]]
        kwargs[_snake(k)] = v
    return op.fn(*args, **kwargs), op.index(kwargs)
