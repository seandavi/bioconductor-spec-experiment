"""Brute-force helpers and hypothesis strategies shared by tests."""

from hypothesis import strategies as st

from ranges_e0b import IRanges


def se(x: IRanges) -> list[tuple[int, int]]:
    return list(zip(x.start.tolist(), x.end.tolist()))


def positions(x: IRanges) -> set[int]:
    out: set[int] = set()
    for s, w in zip(x.start.tolist(), x.width.tolist()):
        out.update(range(s, s + w))
    return out


def runs(pos: set[int]) -> list[tuple[int, int]]:
    """Maximal runs of consecutive integers as (start, end)."""
    out: list[tuple[int, int]] = []
    for p in sorted(pos):
        if out and out[-1][1] == p - 1:
            out[-1] = (out[-1][0], p)
        else:
            out.append((p, p))
    return out


@st.composite
def iranges_st(draw: st.DrawFn, max_size: int = 12, allow_empty: bool = True) -> IRanges:
    n = draw(st.integers(0, max_size))
    starts = [draw(st.integers(-5, 25)) for _ in range(n)]
    widths = [draw(st.integers(0 if allow_empty else 1, 8)) for _ in range(n)]
    return IRanges(start=starts, width=widths) if n else IRanges()
