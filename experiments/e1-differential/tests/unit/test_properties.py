"""Python-only invariants (no R): algebraic identities the kernels must satisfy on random input."""

import numpy as np
import pytest

from ranges_e1 import GRanges, IRanges


def rand_ir(rng: np.random.Generator, n: int) -> IRanges:
    return IRanges(start=rng.integers(-5, 40, n), width=rng.choice([0, 1, 2, 5, 9], n))


def covered(x: IRanges, lo: int = -10, hi: int = 60) -> np.ndarray:
    """Boolean membership of positions lo..hi (empty ranges cover nothing)."""
    pos = np.arange(lo, hi + 1)
    out = np.zeros(len(pos), dtype=bool)
    for s, w in zip(x.start, x.width):
        out |= (pos >= s) & (pos < s + w)
    return out


SEEDS = range(25)


@pytest.mark.parametrize("seed", SEEDS)
def test_reduce_is_idempotent_and_disjoint(seed: int) -> None:
    x = rand_ir(np.random.default_rng(seed), 30)
    r = x.reduce(drop_empty_ranges=True)
    assert r.reduce(drop_empty_ranges=True).start.tolist() == r.start.tolist()
    assert (r.start[1:] > r.end[:-1] + 1).all()  # sorted with a gap between neighbours
    assert (covered(r) == covered(x)).all()


@pytest.mark.parametrize("seed", SEEDS)
def test_gaps_complement_coverage(seed: int) -> None:
    x = rand_ir(np.random.default_rng(seed), 20)
    g = x.gaps(start=-10, end=60)
    assert (covered(g) == ~covered(x)).all()


@pytest.mark.parametrize("seed", SEEDS)
def test_set_operations_match_membership(seed: int) -> None:
    rng = np.random.default_rng(seed)
    x, y = rand_ir(rng, 15), rand_ir(rng, 15)
    cx, cy = covered(x), covered(y)
    assert (covered(x.union(y)) == (cx | cy)).all()
    assert (covered(x.intersect(y)) == (cx & cy)).all()
    assert (covered(x.setdiff(y)) == (cx & ~cy)).all()


@pytest.mark.parametrize("seed", SEEDS)
def test_coverage_sums_to_total_width(seed: int) -> None:
    x = rand_ir(np.random.default_rng(seed), 20)
    cov = x.coverage(width=100)
    clipped = np.clip(x.end, 0, 100) - np.clip(x.start - 1, 0, 100)
    assert (cov.decode() * 1).sum() == clipped[x.width > 0].sum()
    assert len(cov) == 100


@pytest.mark.parametrize("seed", SEEDS)
def test_find_overlaps_agrees_with_brute_force(seed: int) -> None:
    rng = np.random.default_rng(seed)
    q, s = rand_ir(rng, 12), rand_ir(rng, 14)
    hits = q.find_overlaps(s)
    got = set(zip(hits.query_hits.tolist(), hits.subject_hits.tolist()))
    want = {
        (i, j)
        for i in range(len(q))
        for j in range(len(s))
        if q.start[i] <= s.end[j] and s.start[j] <= q.end[i]
    }
    assert got == want


@pytest.mark.parametrize("seed", SEEDS)
def test_order_sorts_and_rank_is_consistent(seed: int) -> None:
    x = rand_ir(np.random.default_rng(seed), 25)
    o = x.order()
    assert sorted(o.tolist()) == list(range(25))
    key = list(zip(x.start[o].tolist(), x.width[o].tolist()))
    assert key == sorted(key)
    r = x.rank(ties_method="first")
    assert sorted(r.tolist()) == list(range(1, 26))
    assert not x.unique().duplicated().any()


@pytest.mark.parametrize("seed", SEEDS)
def test_granges_ops_respect_strand_groups(seed: int) -> None:
    rng = np.random.default_rng(seed)
    n = 20
    g = GRanges(
        rng.choice(["a", "b"], n), start=rng.integers(1, 40, n), width=rng.choice([0, 1, 4], n),
        strand=rng.choice(["+", "-", "*"], n),
    )
    r = g.reduce()
    for sn in ("a", "b"):
        for st in ("+", "-", "*"):
            sub = g[(g.seqnames == sn) & (g.strand == st)]
            out = r[(r.seqnames == sn) & (r.strand == st)]
            same = covered(IRanges(out.start, width=out.width)) == covered(IRanges(sub.start, width=sub.width))
            assert same.all()
            assert len(out) <= len(sub)
    assert len(g.reduce(ignore_strand=True)) <= len(r)
    assert set(g.sort().strand) <= {"+", "-", "*"}
