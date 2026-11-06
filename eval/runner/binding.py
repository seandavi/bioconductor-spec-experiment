"""Bind a case's call to R argument names via ops.json (stdlib only; shared by oracle and port runner)."""

from __future__ import annotations


class BindError(Exception):
    """R itself would reject this argument list (unused / duplicated argument)."""


def bind(entry: dict, args: list, kwargs: dict) -> tuple[dict, list, list]:
    """Return (bound, extra_positional, extra_named).

    bound: R argument name -> value, in call order. extra_positional: values for `...` ('objects' ops).
    extra_named: [(R name, value)] for `...` ('mcols' ops).
    """
    positional = entry["positional"]
    named = set(entry["named"])
    dots = entry["dots"]
    bound: dict = {}
    extra_pos: list = []
    extra_named: list = []
    for i, a in enumerate(args):
        if i < len(positional):
            bound[positional[i]] = a
        else:
            extra_pos.append(a)
    for k, v in kwargs.items():
        if k in named:
            if k in bound:
                raise BindError(f"formal argument '{k}' matched by multiple actual arguments")
            bound[k] = v
        elif dots == "mcols":
            extra_named.append((k, v))
        else:
            raise BindError(f"unused argument ({k})")
    if extra_pos and dots != "objects":
        raise BindError(f"unused argument ({len(extra_pos)} extra positional)")
    return bound, extra_pos, extra_named
