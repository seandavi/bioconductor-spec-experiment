"""Random range / seqinfo data for the heldout generator (stdlib only).

Scope comes from TARGET.md only. Nothing here reads a suite, a port or R output.
"""

from __future__ import annotations

import random
from dataclasses import dataclass, field

INT_MAX = 2147483647
CHRS = ["chr1", "chr2", "chr3", "chrX", "chrM", "chrY", "scaf_1"]
NAMES_POOL = ["a", "b", "c", "d", "e", "f", "g", "h"]
GENOMES = ["hg19", "mm10", "toyGenome"]
MCOL_NAMES = ["m1", "score", "tag", "flag", "m2"]


def ivec(vals, names=None) -> dict:
    d = {"int": list(vals)}
    if names is not None:
        d["names"] = list(names)
    return d


def dvec(vals, names=None) -> dict:
    d = {"dbl": list(vals)}
    if names is not None:
        d["names"] = list(names)
    return d


def cvec(vals, names=None) -> dict:
    d = {"chr": list(vals)}
    if names is not None:
        d["names"] = list(names)
    return d


def lvec(vals, names=None) -> dict:
    d = {"lgl": list(vals)}
    if names is not None:
        d["names"] = list(names)
    return d


def scalar(v):
    """Case-schema scalar: integral doubles must be written as {"dbl": [v]}."""
    if isinstance(v, float) and v.is_integer():
        return {"dbl": [v]}
    return v


@dataclass
class SeqinfoSpec:
    levels: list[str]
    lengths: list[int | None] | None = None
    circular: list[bool | None] | None = None
    genome: list[str | None] | None = None

    def body(self) -> dict:
        b: dict = {"seqnames": cvec(self.levels)}
        if self.lengths is not None:
            b["seqlengths"] = ivec(self.lengths)
        if self.circular is not None:
            b["isCircular"] = lvec(self.circular)
        if self.genome is not None:
            b["genome"] = cvec(self.genome)
        return b

    def feature_tags(self) -> list[str]:
        tags = []
        if self.lengths is not None and any(v is not None for v in self.lengths):
            tags.append("seqinfo_lengths")
        if self.circular is not None and any(self.circular):
            tags.append("seqinfo_circular")
        if self.genome is not None and any(self.genome):
            tags.append("seqinfo_genome")
        if len(self.levels) > 1:
            tags.append("seqinfo_multi")
        if self.levels != sorted(self.levels):
            tags.append("seqinfo_unsorted")
        return tags


def sample_seqinfo(rng: random.Random, *, min_levels: int = 1, small_lengths: bool = False) -> SeqinfoSpec:
    nl = max(min_levels, rng.choice([1, 2, 2, 3, 3, 4]))
    levels = rng.sample(CHRS, nl)  # random order => usually not sorted
    mode = rng.random()
    lengths = None
    if mode < 0.30 or small_lengths:
        lo, hi = (15, 45) if small_lengths else (100, 400)
        lengths = [rng.randint(lo, hi) for _ in levels]
        if rng.random() < 0.3:
            lengths[rng.randrange(nl)] = None
    elif mode < 0.40:
        lengths = [None] * nl
    circular = None
    if rng.random() < 0.35:
        circular = [rng.random() < 0.3 for _ in levels]
        if "chrM" in levels and rng.random() < 0.7:
            circular[levels.index("chrM")] = True
        if rng.random() < 0.15:
            circular[rng.randrange(nl)] = None
    genome = None
    if rng.random() < 0.4:
        g = rng.choice(GENOMES)
        genome = [g] * nl
        if nl > 1 and rng.random() < 0.15:
            genome[rng.randrange(nl)] = None
    return SeqinfoSpec(levels, lengths, circular, genome)


def variant_seqinfo(rng: random.Random, base: SeqinfoSpec) -> SeqinfoSpec:
    """A seqinfo overlapping `base`: other order / subset / extra level / maybe conflicting lengths."""
    levels = list(base.levels)
    rng.shuffle(levels)
    if len(levels) > 1 and rng.random() < 0.4:
        levels = levels[:-1]
    if rng.random() < 0.5:
        extra = [c for c in CHRS if c not in levels]
        if extra:
            levels.append(rng.choice(extra))
    lengths = None
    if rng.random() < 0.6:
        lengths = []
        for lv in levels:
            if lv in base.levels and base.lengths and base.lengths[base.levels.index(lv)] is not None:
                v = base.lengths[base.levels.index(lv)]
                if rng.random() < 0.15:
                    v += rng.randint(1, 5)  # conflicting length
                lengths.append(v if rng.random() < 0.8 else None)
            else:
                lengths.append(rng.choice([None, rng.randint(100, 400)]))
    genome = None
    if rng.random() < 0.4:
        gname = rng.choice(GENOMES) if rng.random() < 0.3 or not base.genome else (base.genome[0] or "hg19")
        genome = [gname] * len(levels)
    circ = None
    if rng.random() < 0.3:
        circ = [rng.random() < 0.25 for _ in levels]
    return SeqinfoSpec(levels, lengths, circ, genome)


@dataclass
class Ranges:
    """A python-side range set (IRanges or GRanges) that serializes to a case literal."""

    start: list[int]
    width: list[int]
    names: list[str] | None = None
    mcols: list[tuple[str, dict]] = field(default_factory=list)
    seqnames: list[str] | None = None  # GRanges iff not None
    strand: list[str] | None = None
    seqinfo: SeqinfoSpec | None = None
    seqinfo_ref: str | None = None
    form: str = "sw"  # sw | se | ew

    @property
    def n(self) -> int:
        return len(self.start)

    @property
    def is_gr(self) -> bool:
        return self.seqnames is not None

    def end(self) -> list[int]:
        return [s + w - 1 for s, w in zip(self.start, self.width)]

    def literal(self) -> dict:
        body: dict = {}
        if self.is_gr:
            body["seqnames"] = cvec(self.seqnames)
        if self.form == "sw":
            body["start"], body["width"] = ivec(self.start), ivec(self.width)
        elif self.form == "se":
            body["start"], body["end"] = ivec(self.start), ivec(self.end())
        else:
            body["end"], body["width"] = ivec(self.end()), ivec(self.width)
        if self.names is not None:
            body["names"] = cvec(self.names)
        if self.is_gr:
            if self.strand is not None:
                body["strand"] = cvec(self.strand)
            if self.seqinfo_ref is not None:
                body["seqinfo"] = {"$ref": self.seqinfo_ref}
            elif self.seqinfo is not None:
                body["seqinfo"] = self.seqinfo.body()
        if self.mcols:
            body["mcols"] = [[k, v] for k, v in self.mcols]
        return {"GRanges" if self.is_gr else "IRanges": body}

    def effective_levels(self) -> list[str]:
        if self.seqinfo is not None:
            return self.seqinfo.levels
        seen: list[str] = []
        for s in self.seqnames or []:
            if s not in seen:
                seen.append(s)
        return seen
