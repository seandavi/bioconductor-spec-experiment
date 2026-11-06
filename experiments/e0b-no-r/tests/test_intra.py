import pytest

from ranges_e0b import IRanges


def se(x: IRanges) -> list[tuple[int, int]]:
    return list(zip(x.start.tolist(), x.end.tolist()))


def test_shift() -> None:
    # R: shift(IRanges(c(1,5), c(3,9)), 10)
    x = IRanges(start=[1, 5], end=[3, 9], names=["a", "b"], mcols={"v": [1, 2]})
    y = x.shift(10)
    assert se(y) == [(11, 13), (15, 19)]
    assert y.names == ["a", "b"] and y.mcols["v"].to_list() == [1, 2]
    assert se(x.shift([1, -1])) == [(2, 4), (4, 8)]
    assert se(x.shift(-5)) == [(-4, -2), (0, 4)]


def test_narrow_relative() -> None:
    x = IRanges(start=1, end=10)
    # R: narrow(IRanges(1,10), start=3)
    assert se(x.narrow(start=3)) == [(3, 10)]
    assert se(x.narrow(end=4)) == [(1, 4)]
    assert se(x.narrow(start=2, end=5)) == [(2, 5)]
    assert se(x.narrow(start=2, width=3)) == [(2, 4)]
    assert se(x.narrow(end=5, width=2)) == [(4, 5)]


def test_narrow_negative() -> None:
    x = IRanges(start=11, end=20)
    # R: narrow(IRanges(11,20), start=-3)  # last 3 positions
    assert se(x.narrow(start=-3)) == [(18, 20)]
    assert se(x.narrow(end=-2)) == [(11, 19)]
    assert se(x.narrow(start=-5, end=-3)) == [(16, 18)]


def test_narrow_width_only_keeps_leading_positions() -> None:
    # UNVERIFIED: R behaviour for width alone
    assert se(IRanges(start=11, end=20).narrow(width=3)) == [(11, 13)]


def test_narrow_out_of_range_errors() -> None:
    x = IRanges(start=1, end=10)
    assert x.narrow(start=11).width.tolist() == [0]  # empty at end + 1
    with pytest.raises(ValueError):
        x.narrow(start=12)
    with pytest.raises(ValueError):
        x.narrow(end=11)
    with pytest.raises(ValueError):
        x.narrow(start=0)
    with pytest.raises(ValueError):
        x.narrow(start=1, end=5, width=5)


def test_narrow_to_empty() -> None:
    # narrowing to nothing: start = end + 1
    x = IRanges(start=1, end=10).narrow(start=4, end=3)
    assert x.width.tolist() == [0] and x.start.tolist() == [4]


def test_narrow_noop() -> None:
    x = IRanges(start=[1, 5], end=[3, 9])
    assert x.narrow().equals(x)


def test_resize_fix() -> None:
    x = IRanges(start=5, end=10)
    # R: resize(IRanges(5,10), 3, fix="start")
    assert se(x.resize(3)) == [(5, 7)]
    assert se(x.resize(3, fix="end")) == [(8, 10)]
    assert se(x.resize(10)) == [(5, 14)]
    assert se(x.resize(10, fix="end")) == [(1, 10)]


def test_resize_center() -> None:
    x = IRanges(start=5, end=10)
    assert se(x.resize(2, fix="center")) == [(7, 8)]
    assert se(x.resize(10, fix="center")) == [(3, 12)]
    # UNVERIFIED: odd difference rounds start down (floor)
    assert se(x.resize(5, fix="center")) == [(5, 9)]
    assert se(x.resize(7, fix="center")) == [(4, 10)]


def test_resize_per_range_fix_and_zero() -> None:
    x = IRanges(start=[5, 5], end=[10, 10])
    assert se(x.resize(2, fix=["start", "end"])) == [(5, 6), (9, 10)]
    z = x.resize(0)
    assert z.width.tolist() == [0, 0]
    with pytest.raises(ValueError):
        x.resize(-1)
    with pytest.raises(ValueError):
        x.resize(2, fix="middle")


def test_flank_start() -> None:
    x = IRanges(start=10, end=20)
    # R: flank(IRanges(10,20), 5)
    assert se(x.flank(5)) == [(5, 9)]
    # R: flank(IRanges(10,20), 5, start=FALSE)
    assert se(x.flank(5, start=False)) == [(21, 25)]


def test_flank_both() -> None:
    x = IRanges(start=10, end=20)
    # R: flank(IRanges(10,20), 5, both=TRUE)
    assert se(x.flank(5, both=True)) == [(5, 14)]
    assert se(x.flank(5, start=False, both=True)) == [(16, 25)]


def test_flank_negative_is_inward() -> None:
    x = IRanges(start=10, end=20)
    # R: flank(IRanges(10,20), -3)
    assert se(x.flank(-3)) == [(10, 12)]
    assert se(x.flank(-3, start=False)) == [(18, 20)]


def test_flank_zero() -> None:
    assert IRanges(start=10, end=20).flank(0).width.tolist() == [0]


def test_flank_per_range_start() -> None:
    x = IRanges(start=[10, 10], end=[20, 20])
    assert se(x.flank(2, start=[True, False])) == [(8, 9), (21, 22)]


def test_promoters_defaults() -> None:
    # R: promoters(IRanges(5000, 5100))
    assert se(IRanges(start=5000, end=5100).promoters()) == [(3000, 5199)]
    assert se(IRanges(start=50, end=60).promoters(upstream=10, downstream=5)) == [(40, 54)]
    assert se(IRanges(start=50, end=60).promoters(0, 0)) == [(50, 49)]
    with pytest.raises(ValueError):
        IRanges(start=50, end=60).promoters(-1, 0)
