import warnings

import pytest

from ranges_e0 import Seqinfo


def make() -> Seqinfo:
    return Seqinfo(["chr1", "chr2", "chrM"], [100, 200, None], [False, False, True], "hg19")


def test_accessors():
    si = make()
    assert si.seqnames == si.seqlevels == ("chr1", "chr2", "chrM")
    assert si.seqlengths == {"chr1": 100, "chr2": 200, "chrM": None}
    assert si.is_circular == {"chr1": False, "chr2": False, "chrM": True}
    assert si.genome == {"chr1": "hg19", "chr2": "hg19", "chrM": "hg19"}
    assert len(si) == 3 and "chr2" in si and si.index("chrM") == 2


def test_defaults_are_na():
    si = Seqinfo(["a", "b"])
    assert si.seqlengths == {"a": None, "b": None}
    assert si.is_circular == {"a": None, "b": None}
    assert si.genome == {"a": None, "b": None}


def test_validation():
    with pytest.raises(ValueError):
        Seqinfo(["a", "a"])
    with pytest.raises(ValueError):
        Seqinfo(["a"], [1, 2])
    with pytest.raises(ValueError):
        Seqinfo(["a"], [-1])


def test_merge_order_and_na_fill():
    # R: merge(si, Seqinfo(c("chr2","chr3"), c(200,300))) -> chr1 chr2 chrM chr3 ; chr3 genome NA
    a = make()
    b = Seqinfo(["chr2", "chr3"], [200, 300])
    with pytest.warns(UserWarning, match="not in the other"):
        m = a.merge(b)
    assert m.seqnames == ("chr1", "chr2", "chrM", "chr3")
    assert m.seqlengths == {"chr1": 100, "chr2": 200, "chrM": None, "chr3": 300}
    assert m.genome["chr3"] is None and m.genome["chr1"] == "hg19"
    assert m.is_circular["chr3"] is None


def test_merge_fills_na_from_other_side_without_warning():
    a = Seqinfo(["chr1", "chr2"], [None, 5])
    b = Seqinfo(["chr2", "chr1"], [5, 100], genome="g")
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # only one side has extra levels: no warning (none here)
        m = a.merge(b)
    assert m.seqnames == ("chr1", "chr2")
    assert m.seqlengths == {"chr1": 100, "chr2": 5}
    assert m.genome == {"chr1": "g", "chr2": "g"}


def test_merge_conflicts_raise():
    # R: "sequence chr1 has incompatible seqlengths" / genomes / circularity flags
    a = make()
    with pytest.raises(ValueError, match="chr1.*seqlengths"):
        a.merge(Seqinfo(["chr1"], [50]))
    with pytest.raises(ValueError, match="genomes"):
        a.merge(Seqinfo(["chr1"], [100], genome="hg38"))
    with pytest.raises(ValueError, match="circularity"):
        a.merge(Seqinfo(["chrM"], None, [False]))


def test_merge_na_vs_value_is_compatible():
    a = make()
    m = a.merge(Seqinfo(["chr1"], None, None, "hg19"))
    assert m == a


def test_equality_and_from_mapping():
    assert make() == make()
    assert Seqinfo.from_mapping({"a": 1, "b": None}).seqlengths == {"a": 1, "b": None}
    assert Seqinfo(["a"]) != Seqinfo(["a"], [3])
