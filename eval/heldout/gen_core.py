"""Case-building core for the heldout generator: seeds, G (per-case state), call assembly."""

from __future__ import annotations

import hashlib
import random

from gen_data import (
    INT_MAX, MCOL_NAMES, NAMES_POOL, Ranges, SeqinfoSpec, cvec, dvec, ivec, lvec,
    sample_seqinfo, variant_seqinfo,
)

NAMESPACE = "bioc-spec-eval/heldout/v1"


def seed_for(op: str, k: int) -> int:
    return int.from_bytes(hashlib.sha256(f"{NAMESPACE}/{op}/{k}".encode()).digest()[:8], "big")


def mk_call(rng: random.Random, entry: dict, chosen: dict, dots: list | None = None,
            mcols_kw: list | None = None) -> dict:
    """Assemble {op, args, kwargs}. `chosen` maps R argument name -> arg value (ordered).

    A random contiguous prefix of the op's positional formals is passed positionally; the
    rest go by name. Extra objects (`dots`) force every positional formal to be positional.
    """
    pos_names = entry["positional"]
    prefix = []
    for nm in pos_names:
        if nm not in chosen:
            break
        prefix.append(nm)
    if dots:
        m = len(prefix)
    elif rng.random() < 0.4:
        m = len(prefix)
    else:
        m = rng.randint(1 if prefix else 0, len(prefix))
    args = [chosen[nm] for nm in prefix[:m]]
    args += list(dots or [])
    kwargs = {nm: v for nm, v in chosen.items() if nm not in prefix[:m]}
    for nm, v in mcols_kw or []:
        kwargs[nm] = v
    return {"op": entry["op"], "args": args, "kwargs": kwargs}


class G:
    """Per-case state: rng, inputs, tags, invalid-injection slot."""

    def __init__(self, rng: random.Random, entry: dict, cls: str, k: int):
        self.rng = rng
        self.entry = entry
        self.cls = cls
        self.k = k
        self.inputs: dict[str, dict] = {}
        self.tags: list[str] = []
        self.bad_slot: str | None = None
        self.extreme = rng.random() < 0.02
        self.p_names = rng.choice([0.0, 0.15, 0.4])
        self.p_mcols = rng.choice([0.0, 0.15, 0.4])
        self.p_zero = rng.choice([0.0, 0.1, 0.25])
        self.p_empty = 0.07
        self.counter = 0

    # ---- tags / invalid injection -------------------------------------------------
    def tag(self, *tags: str) -> None:
        for t in tags:
            if t not in self.tags:
                self.tags.append(t)

    def begin(self, slots: list[str], p_bad: float = 0.22) -> None:
        """Maybe pick one argument slot to receive a deliberately invalid value."""
        if slots and self.rng.random() < p_bad:
            self.bad_slot = self.rng.choice(slots)

    def arg(self, slot: str, good, bads: list):
        """Return `good` unless this slot was picked for invalid injection."""
        if self.bad_slot == slot:
            self.tag("invalid")
            return self.rng.choice(bads)
        return good

    @property
    def is_bad(self) -> bool:
        return self.bad_slot is not None

    # ---- ranges -------------------------------------------------------------------
    def size(self) -> int:
        r = self.rng
        if r.random() < self.p_empty:
            return 0
        return r.choice([1, 2, 3, 3, 4, 5, 5, 6, 7, 8, 9, 10, 12])

    def _coords(self, n: int, window=(-5, 60), maxw: int = 12, minw: int = 0) -> tuple[list[int], list[int]]:
        r = self.rng
        if self.extreme:
            top = r.random() < 0.7
            base = INT_MAX - 70 if top else -INT_MAX + 5
            starts = [base + r.randint(0, 40) for _ in range(n)]
            widths = [r.randint(0, 8) for _ in range(n)]
            self.tag("extreme")
            return starts, widths
        starts, widths = [], []
        for _ in range(n):
            starts.append(r.randint(*window))
            if r.random() < self.p_zero:
                widths.append(0)
            else:
                widths.append(max(minw, r.choice([1, 1, 2, 3, 4, 5, 6, 8, maxw])))
        return starts, widths

    def mcol(self, n: int, complex_ok: bool = False) -> dict:
        r = self.rng
        kinds = ["int", "dbl", "chr", "lgl"] + (["factor", "ilist"] if complex_ok and r.random() < 0.4 else [])
        kind = r.choice(kinds)
        na = r.random() < 0.4
        def maybe(v):
            return None if na and r.random() < 0.25 else v
        if kind == "int":
            return ivec([maybe(r.randint(-5, 20)) for _ in range(n)])
        if kind == "dbl":
            return dvec([maybe(r.randint(-8, 40) / 4) for _ in range(n)])
        if kind == "chr":
            return cvec([maybe(r.choice(["x", "y", "zz", "", "x"])) for _ in range(n)])
        if kind == "lgl":
            return lvec([maybe(r.random() < 0.5) for _ in range(n)])
        if kind == "factor":
            levels = r.sample(["lo", "mid", "hi", "zz"], r.randint(1, 3))
            return {"Factor": {"levels": levels, "codes": [maybe(r.randint(1, len(levels))) for _ in range(n)]}}
        elems = [[r.randint(-3, 9) for _ in range(r.randint(0, 3))] for _ in range(n)]
        return {"IntegerList": {"elements": elems}}

    def mcols_for(self, n: int, force: bool | None = None, complex_ok: bool = False) -> list:
        r = self.rng
        if not (force if force is not None else r.random() < self.p_mcols):
            return []
        names = r.sample(MCOL_NAMES, r.choice([1, 1, 2]))
        return [(nm, self.mcol(n, complex_ok)) for nm in names]

    def names_for(self, n: int, force: bool | None = None) -> list[str] | None:
        r = self.rng
        if n == 0 or not (force if force is not None else r.random() < self.p_names):
            return None
        pool = NAMES_POOL[: r.choice([3, 8])]  # small pool => duplicates
        if r.random() < 0.5 and n <= len(NAMES_POOL):
            return r.sample(NAMES_POOL, n)
        return [r.choice(pool) for _ in range(n)]

    def strand_for(self, n: int, mode: str | None = None) -> list[str] | None:
        r = self.rng
        mode = mode or r.choices(["mixed", "absent", "plus", "star"], [0.55, 0.2, 0.15, 0.1])[0]
        if mode == "absent":
            return None
        if mode == "plus":
            return ["+"] * n
        if mode == "star":
            return ["*"] * n
        return [r.choice("+-*") for _ in range(n)]

    def ranges(self, n: int | None = None, *, gr: bool | None = None, si: SeqinfoSpec | None = None,
               si_ref: str | None = None, names: bool | None = None, mcols: bool | None = None,
               strand: str | None = None, window=(-5, 60), maxw: int = 12, minw: int = 0,
               complex_mcols: bool = False, form: str | None = None) -> Ranges:
        r = self.rng
        gr = (self.cls == "GRanges") if gr is None else gr
        n = self.size() if n is None else n
        start, width = self._coords(n, window, maxw, minw)
        rs = Ranges(start, width, self.names_for(n, names), self.mcols_for(n, mcols, complex_mcols),
                    form=form or r.choice(["sw", "sw", "se", "ew"]))
        if gr:
            if si is None and si_ref is None:
                si = sample_seqinfo(r) if r.random() < 0.75 else None
            levels_pool = si.levels if si is not None else r.sample(["chr1", "chr2", "chr3", "chrX"], r.randint(1, 3))
            used = r.sample(levels_pool, r.randint(1, len(levels_pool)))
            rs.seqnames = [r.choice(used) for _ in range(n)]
            rs.strand = self.strand_for(n, strand)
            rs.seqinfo = si
            rs.seqinfo_ref = si_ref
        return rs

    def seqinfo_input(self, name: str = "si", spec: SeqinfoSpec | None = None, **kw) -> SeqinfoSpec:
        spec = spec or sample_seqinfo(self.rng, **kw)
        self.inputs[name] = {"Seqinfo": spec.body()}
        return spec

    def add(self, name: str, rs: Ranges) -> Ranges:
        self.inputs[name] = rs.literal()
        return rs

    def shared_seqinfo(self, name: str = "si", **kw) -> tuple[SeqinfoSpec, str | None]:
        """A seqinfo for several GRanges: inline in each, or a shared Seqinfo input by $ref."""
        spec = sample_seqinfo(self.rng, **kw)
        if self.rng.random() < 0.4:
            self.seqinfo_input(name, spec)
            return spec, name
        return spec, None

    def pair(self, *, derived: float = 0.4, same_si: float = 0.75, n_a: int | None = None,
             n_b: int | None = None, **kw) -> tuple[Ranges, Ranges]:
        """Inputs x, y of the case class; y is independent or derived from x."""
        r = self.rng
        a_si = b_si = None
        a_ref = b_ref = None
        if self.cls == "GRanges":
            spec, ref = self.shared_seqinfo(**{k: v for k, v in kw.items() if k in ("small_lengths", "min_levels")})
            a_si, a_ref = spec, ref
            if r.random() < same_si:
                b_si, b_ref = a_si, a_ref
            else:
                b_si, b_ref = variant_seqinfo(r, spec), None
        rkw = {k: v for k, v in kw.items() if k not in ("small_lengths", "min_levels")}
        a = self.ranges(n_a, si=a_si, si_ref=a_ref, **rkw)
        if r.random() < derived and a.n > 0 and n_b is None:
            b = self.derive(a, b_si, b_ref)
        else:
            b = self.ranges(n_b, si=b_si, si_ref=b_ref, **rkw)
            if self.cls == "GRanges" and b_si is a_si and a.seqnames:
                pass
        self.add("x", a)
        self.add("y", b)
        return a, b

    def parallel_pair(self, n: int | None = None, m: int | None = None, **kw) -> tuple[Ranges, Ranges]:
        """x of length n and y of length m (default n) with y[i] a nudged copy of x[i % n], so pairs are compatible."""
        r = self.rng
        si = ref = None
        if self.cls == "GRanges":
            si, ref = self.shared_seqinfo(**{k: v for k, v in kw.items() if k in ("small_lengths", "min_levels")})
        rkw = {k: v for k, v in kw.items() if k not in ("small_lengths", "min_levels")}
        a = self.ranges(n, si=si, si_ref=ref, **rkw)
        m = a.n if m is None else m
        if a.n == 0:
            b = self.ranges(m, si=si, si_ref=ref, **rkw)
        else:
            idx = [i % a.n for i in range(m)]
            b = Ranges([a.start[i] + r.choice([-2, -1, 0, 0, 1, 2]) for i in idx],
                       [max(0, a.width[i] + r.choice([-2, 0, 0, 1, 3])) for i in idx], None, [],
                       form=r.choice(["sw", "se", "ew"]))
            if r.random() < self.p_names:
                b.names = self.names_for(m, True)
            if a.is_gr:
                b.seqnames = [a.seqnames[i] if r.random() < 0.9 else r.choice(a.effective_levels()) for i in idx]
                st = a.strand
                b.strand = None if st is None else [st[i] if r.random() < 0.8 else r.choice("+-*") for i in idx]
                b.seqinfo, b.seqinfo_ref = si, ref
        self.add("x", a)
        self.add("y", b)
        return a, b

    def derive(self, a: Ranges, si, ref) -> Ranges:
        """A range set resembling `a`: some copied, some nudged, some dropped or duplicated."""
        r = self.rng
        idx = [i for i in range(a.n) if r.random() < 0.8] or [0]
        if r.random() < 0.3:
            idx += [r.choice(idx) for _ in range(r.randint(1, 2))]
        start = [a.start[i] + r.choice([0, 0, 0, 1, -1, 2]) for i in idx]
        width = [max(0, a.width[i] + r.choice([0, 0, 0, 1, -1, 3])) for i in idx]
        rs = Ranges(start, width, None, [], form=r.choice(["sw", "se", "ew"]))
        if r.random() < self.p_names:
            rs.names = self.names_for(len(idx), True)
        if r.random() < self.p_mcols:
            rs.mcols = self.mcols_for(len(idx), True)
        if a.is_gr:
            rs.seqnames = [a.seqnames[i] for i in idx]
            st = a.strand
            rs.strand = None if st is None else [x if r.random() < 0.8 else r.choice("+-*") for x in (st[i] for i in idx)]
            rs.seqinfo, rs.seqinfo_ref = si, ref
            if si is not None:
                lv = si.levels
                rs.seqnames = [s if s in lv else lv[0] for s in rs.seqnames]
        return rs

    # ---- argument samplers -----------------------------------------------------------
    def int_arg(self, n: int, lo: int, hi: int, *, vec_p: float = 0.3, na_p: float = 0.0, dbl_p: float = 0.0):
        """Scalar or per-range int argument (bare scalar / int vector / dbl vector / NA)."""
        r = self.rng
        if r.random() < na_p:
            self.tag("na_arg")
            return {"int": [None]} if r.random() < 0.7 or n == 0 else ivec([None if i == 0 else r.randint(lo, hi) for i in range(n)])
        if r.random() < dbl_p:
            if n > 0 and r.random() < vec_p:
                return dvec([float(r.randint(lo, hi)) for _ in range(n)])
            v = r.choice([float(r.randint(lo, hi)), r.randint(lo, hi) + 0.5])
            return {"dbl": [v]}
        if n > 0 and r.random() < vec_p:
            return ivec([r.randint(lo, hi) for _ in range(n)])
        return r.randint(lo, hi)

    def flag(self, p_na: float = 0.0):
        if self.rng.random() < p_na:
            self.tag("na_arg")
            return {"lgl": [None]}
        return self.rng.random() < 0.5
