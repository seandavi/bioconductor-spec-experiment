import pytest

from ranges_e0 import IRanges

X = IRanges(10, 20)


def tup(x):
    return x.to_tuples()


def test_shift():
    assert tup(X.shift(-3)) == [(7, 17)]  # R: shift(IRanges(10,20), -3) -> 7-17
    assert tup(IRanges([1, 5], [2, 6]).shift([1, 2])) == [(2, 3), (7, 8)]


def test_shift_keeps_names_and_mcols():
    x = IRanges([10, 20], width=5, names=["a", "b"], mcols={"s": [1, 2]})
    y = x.shift(1)
    assert y.names == ["a", "b"] and y.mcols["s"].to_list() == [1, 2]
    assert x.shift(1, use_names=False).names is None


def test_narrow_start_only():
    assert tup(X.narrow(start=2)) == [(11, 20)]  # R: narrow(IRanges(10,20), start=2)
    assert tup(X.narrow(start=-3)) == [(18, 20)]  # negative = from end


def test_narrow_end_only():
    assert tup(X.narrow(end=-2)) == [(10, 19)]


def test_narrow_start_end():
    assert tup(X.narrow(start=-3, end=-1)) == [(18, 20)]
    assert tup(X.narrow(start=12, end=11)) == [(21, 20)]  # zero-width at the end


def test_narrow_with_width():
    assert tup(X.narrow(start=3, width=2)) == [(12, 13)]
    assert tup(X.narrow(end=-3, width=2)) == [(17, 18)]


def test_narrow_errors():
    # R: narrow() never extends a range, and width needs exactly one of start/end
    for kw in (
        dict(width=3),
        dict(end=15, width=2),
        dict(start=12, end=30),
        dict(start=0),
        dict(start=13),
        dict(start=1, end=2, width=2),
    ):
        with pytest.raises(ValueError):
            X.narrow(**kw)


def test_resize_start_end():
    assert tup(X.resize(5)) == [(10, 14)]
    assert tup(X.resize(5, fix="end")) == [(16, 20)]
    assert tup(X.resize(0)) == [(10, 9)]


def test_resize_center():
    # R: resize(IRanges(10,20), 4, fix="center") -> 13-16 ; width 5 -> 13-17
    assert tup(X.resize(4, fix="center")) == [(13, 16)]
    assert tup(X.resize(5, fix="center")) == [(13, 17)]
    # even-width input: IRanges(10,19) -> 5: 12-16, 4: 13-16
    assert tup(IRanges(10, 19).resize(5, fix="center")) == [(12, 16)]
    assert tup(IRanges(10, 19).resize(4, fix="center")) == [(13, 16)]
    # zero-width (5,4) -> 3 centered: 3-5
    assert tup(IRanges(5, 4).resize(3, fix="center")) == [(3, 5)]


def test_resize_vector_fix_and_errors():
    x = IRanges([10, 10], [20, 20])
    assert tup(x.resize(3, fix=["start", "end"])) == [(10, 12), (18, 20)]
    with pytest.raises(ValueError):
        X.resize(-1)  # R: 'width' values must be non-negative
    with pytest.raises(ValueError):
        X.resize(3, fix="middle")


def test_flank_start_end():
    assert tup(X.flank(3)) == [(7, 9)]
    assert tup(X.flank(3, start=False)) == [(21, 23)]


def test_flank_negative_goes_inside():
    assert tup(X.flank(-3)) == [(10, 12)]
    assert tup(X.flank(-3, start=False)) == [(18, 20)]


def test_flank_both():
    assert tup(X.flank(3, both=True)) == [(7, 12)]
    assert tup(X.flank(3, start=False, both=True)) == [(18, 23)]
    assert tup(X.flank(-3, both=True)) == [(7, 12)]  # both uses abs(width)


def test_flank_zero_width():
    assert tup(IRanges(5, 4).flank(2)) == [(3, 4)]


def test_promoters():
    assert tup(X.promoters()) == [(-1990, 209)]  # R defaults 2000/200
    assert tup(X.promoters(upstream=5, downstream=3)) == [(5, 12)]
    assert tup(IRanges(5, 4).promoters()) == [(-1995, 204)]
    assert tup(X.promoters(0, 0)) == [(10, 9)]
    with pytest.raises(ValueError):
        X.promoters(-1, 0)


def test_empty_input():
    e = IRanges.empty()
    assert len(e.shift(1)) == len(e.narrow(start=1)) == len(e.resize(1)) == len(e.flank(1)) == len(e.promoters()) == 0
