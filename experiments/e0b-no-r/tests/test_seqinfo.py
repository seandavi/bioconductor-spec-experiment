import pytest

from ranges_e0b.seqinfo import Seqinfo


def test_construct_and_broadcast() -> None:
    si = Seqinfo(["chr1", "chr2"], seqlengths=[100, None], is_circular=[False, None], genome="hg38")
    assert si.seqnames == ("chr1", "chr2")
    assert si.seqlengths == (100, None)
    assert si.is_circular == (False, None)
    assert si.genome == ("hg38", "hg38")
    assert len(si) == 2


def test_defaults_are_na() -> None:
    si = Seqinfo(["a"])
    assert si.seqlengths == (None,) and si.is_circular == (None,) and si.genome == (None,)


def test_validation() -> None:
    with pytest.raises(ValueError):
        Seqinfo(["a", "a"])
    with pytest.raises(ValueError):
        Seqinfo(["a", "b"], seqlengths=[1])
    with pytest.raises(ValueError):
        Seqinfo(["a"], seqlengths=[-1])
    with pytest.raises(ValueError):
        Seqinfo(["a"], is_circular=["yes"])


def test_equality_hash_index() -> None:
    a = Seqinfo(["x", "y"], [1, 2])
    assert a == Seqinfo(["x", "y"], [1, 2]) and hash(a) == hash(Seqinfo(["x", "y"], [1, 2]))
    assert a != Seqinfo(["y", "x"], [2, 1])
    assert a.index("y") == 1
    with pytest.raises(KeyError):
        a.index("z")


def test_merge_union_and_order() -> None:
    a = Seqinfo(["chr1", "chr2"], [100, None])
    b = Seqinfo(["chr2", "chr3"], [200, 300])
    m = a.merge(b)
    assert m.seqnames == ("chr1", "chr2", "chr3")
    assert m.seqlengths == (100, 200, 300)  # None yields to a value


def test_merge_conflict() -> None:
    with pytest.raises(ValueError, match="seqlengths"):
        Seqinfo(["a"], [1]).merge(Seqinfo(["a"], [2]))
    with pytest.raises(ValueError, match="genome"):
        Seqinfo(["a"], genome="hg19").merge(Seqinfo(["a"], genome="hg38"))
    with pytest.raises(ValueError, match="is_circular"):
        Seqinfo(["a"], is_circular=True).merge(Seqinfo(["a"], is_circular=False))


def test_merge_genome_na_yields() -> None:
    m = Seqinfo(["a"]).merge(Seqinfo(["a"], genome="hg38", is_circular=False))
    assert m.genome == ("hg38",) and m.is_circular == (False,)


def test_merge_no_common_names_warns() -> None:
    with pytest.warns(UserWarning, match="no sequence names in common"):
        Seqinfo(["a"]).merge(Seqinfo(["b"]))
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        Seqinfo([]).merge(Seqinfo(["b"]))  # empty side: no warning
