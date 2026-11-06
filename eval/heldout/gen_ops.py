"""Per-op heldout case builders. Each takes a G and returns a call dict; inputs live on G."""

from __future__ import annotations

from gen_core import G, mk_call
from gen_data import (
    CHRS, INT_MAX, Ranges, SeqinfoSpec, cvec, dvec, ivec, lvec, sample_seqinfo, variant_seqinfo,
)

# ----------------------------------------------------------------------------- helpers

def unary(g: G, extra=None, *, gr_kw=None, mcols: bool | None = None, **rkw):
    """Input x only, plus optional pre-sampled chosen args."""
    x = g.ranges(mcols=mcols, **rkw)
    g.add("x", x)
    chosen = {"x": {"$ref": "x"}}
    chosen.update(extra or {})
    return x, chosen


def ref(name: str) -> dict:
    return {"$ref": name}


def ignore_strand(g: G, chosen: dict, p: float = 0.3, p_na: float = 0.08) -> None:
    if g.cls == "GRanges" and g.rng.random() < p:
        chosen["ignore.strand"] = g.flag(p_na)
        g.tag("strand")


def bad_enum(slot_vals):
    return slot_vals


# ----------------------------------------------------------------------------- construct

def op_IRanges(g: G):
    r = g.rng
    g.begin(["width", "end", "names"])
    n = g.size()
    s, w = g._coords(n)
    form = r.choice(["sw", "se", "ew", "sw", "se", "ew", "sw", "se", "s", "e", "w", "sew", "none"])
    if n == 0:
        g.tag("empty")
    if 0 in w:
        g.tag("zero_width")
    end = [a + b - 1 for a, b in zip(s, w)]
    ch: dict = {}
    if form in ("sw", "se", "s", "sew"):
        ch["start"] = ivec(s) if r.random() < 0.9 else dvec([float(v) for v in s])
    if form in ("se", "e", "sew", "ew"):
        ch["end"] = ivec(g.arg("end", end, [[a - 4 for a in s]]))
    if form in ("sw", "w", "sew", "ew"):
        ch["width"] = ivec(g.arg("width", w, [[-1] + w[1:] if w else [-1], [None] * len(w)]))
    if form == "sew" and r.random() < 0.3:
        ch["width"] = ivec([x + 1 for x in w])  # inconsistent triple
        g.tag("invalid")
    names = g.names_for(n)
    if names is not None:
        g.tag("names")
        ch["names"] = cvec(g.arg("names", names, [names + ["extra"], names[:-1]]))
    mc = []
    if r.random() < g.p_mcols and n > 0:
        mc = g.mcols_for(n, True)
        g.tag("mcols")
    return mk_call(r, g.entry, ch, mcols_kw=mc)


def op_GRanges(g: G):
    r = g.rng
    g.begin(["seqnames", "strand", "seqinfo", "ranges"])
    n = g.size()
    rs = g.ranges(n, gr=False)
    spec = sample_seqinfo(r)
    use_si = r.random() < 0.7
    if use_si:
        g.seqinfo_input("si", spec)
        levels = spec.levels
    else:
        levels = r.sample(CHRS, r.randint(1, 3))
    sn = [r.choice(levels) for _ in range(n)]
    if n == 0:
        g.tag("empty")
    if 0 in rs.width:
        g.tag("zero_width")
    if rs.names:
        g.tag("names")
    if r.random() < 0.8:
        g.add("r", rs)
        rng_arg = ref("r")
    else:
        rng_arg = {"$call": {"op": "IRanges", "args": [ivec(rs.start)], "kwargs": {"width": ivec(rs.width)}}}
    ch: dict = {}
    sn_arg = cvec(sn)
    if n > 0 and r.random() < 0.15:
        sn_arg = cvec([sn[0]])  # recycled seqnames
    ch["seqnames"] = g.arg("seqnames", sn_arg, [cvec(sn + ["chr1"]), cvec(["not_a_level"] * n if n else ["q"])] if use_si else [cvec(sn + ["chr1"])])
    ch["ranges"] = g.arg("ranges", rng_arg, [ivec([1, 2]), None])
    if r.random() < 0.75:
        st = g.strand_for(n, "mixed" if r.random() < 0.7 else None) or ["*"] * n
        if n > 0 and r.random() < 0.15:
            st = [st[0]]
        ch["strand"] = cvec(g.arg("strand", st, [["x"] * max(n, 1), ["+"] * (n + 1)]))
        if len(set(st)) > 1:
            g.tag("strand")
    if use_si:
        ch["seqinfo"] = g.arg("seqinfo", ref("si"), [cvec(["chr1"])])
        g.tag("seqinfo", *spec.feature_tags())
    elif r.random() < 0.2:
        ch["seqlengths"] = ivec([r.randint(100, 300) for _ in levels], names=levels)
        g.tag("seqinfo", "seqinfo_lengths")
    mc = []
    if r.random() < max(g.p_mcols, 0.2) and n > 0:
        mc = g.mcols_for(n, True)
        g.tag("mcols")
    return mk_call(r, g.entry, ch, mcols_kw=mc)


def op_Seqinfo(g: G):
    r = g.rng
    g.begin(["seqnames", "seqlengths", "isCircular", "genome"])
    spec = sample_seqinfo(r)
    if r.random() < 0.1:
        spec = SeqinfoSpec([])
        g.tag("empty")
    nl = len(spec.levels)
    ch = {"seqnames": cvec(g.arg("seqnames", spec.levels, [spec.levels + spec.levels[:1], [None] + spec.levels[1:]] if nl else [["a", "a"]]))}
    if spec.lengths is not None or r.random() < 0.2:
        sl = spec.lengths if spec.lengths is not None else [r.randint(10, 300) for _ in spec.levels]
        ch["seqlengths"] = ivec(g.arg("seqlengths", sl, [sl + [5], [-3] + sl[1:] if sl else [-3]]))
        g.tag("seqinfo_lengths")
    if spec.circular is not None:
        ch["isCircular"] = lvec(g.arg("isCircular", spec.circular, [spec.circular + [True]]))
        g.tag("seqinfo_circular")
    if spec.genome is not None:
        ch["genome"] = cvec(g.arg("genome", spec.genome, [spec.genome + ["x", "y"]]))
        g.tag("seqinfo_genome")
    g.tag("seqinfo", *spec.feature_tags())
    return mk_call(r, g.entry, ch)


# ----------------------------------------------------------------------------- accessors

def op_plain_accessor(g: G):
    x, ch = unary(g)
    if x.n == 0:
        g.tag("empty")
    return mk_call(g.rng, g.entry, ch)


def op_ranges(g: G):
    x, ch = unary(g, mcols=g.rng.random() < 0.5)
    if g.rng.random() < 0.5:
        ch["use.names"] = g.flag(0.1)
    if g.rng.random() < 0.4:
        ch["use.mcols"] = g.flag(0.1)
    return mk_call(g.rng, g.entry, ch)


def op_mcols(g: G):
    x, ch = unary(g, mcols=g.rng.random() < 0.8, complex_mcols=True)
    if g.rng.random() < 0.4:
        ch["use.names"] = g.flag(0.15)
    return mk_call(g.rng, g.entry, ch)


# ----------------------------------------------------------------------------- metadata setters

def op_set_names(g: G):
    r = g.rng
    g.begin(["value"])
    x, ch = unary(g, names=r.random() < 0.4)
    n = x.n
    mode = r.random()
    if mode < 0.15 or n == 0 and mode < 0.5:
        val = None
    else:
        val = cvec(g.names_for(n, True) or [])
        if n > 0 and r.random() < 0.15:
            nm = list(val["chr"])
            nm[r.randrange(n)] = None
            val = cvec(nm)
            g.tag("na_arg")
    val = g.arg("value", val, [cvec(["a"] * (n + 2)), ivec(list(range(n))), cvec(["a"] * max(n - 1, 0))])
    ch["value"] = val
    return mk_call(r, g.entry, ch)


def op_set_mcols(g: G):
    r = g.rng
    g.begin(["value"])
    x, ch = unary(g, mcols=r.random() < 0.4, complex_mcols=True)
    n = x.n
    if r.random() < 0.2:
        val = None
    else:
        cols = g.mcols_for(n, True, complex_ok=True)
        g.inputs["df"] = {"DataFrame": {"columns": [[k, v] for k, v in cols]}} if cols or n else {"DataFrame": {"columns": [], "nrow": n}}
        if not cols:
            g.inputs["df"] = {"DataFrame": {"columns": [], "nrow": n}}
        val = ref("df")
        if g.bad_slot == "value":
            g.tag("invalid")
            cols2 = g.mcols_for(n + 2, True)
            g.inputs["df"] = {"DataFrame": {"columns": [[k, v] for k, v in cols2]}} if cols2 else {"DataFrame": {"columns": [], "nrow": n + 2}}
    ch["value"] = val
    g.tag("mcols")
    return mk_call(r, g.entry, ch)


# ----------------------------------------------------------------------------- seqinfo

def _gr_or_seqinfo(g: G, allow_seqinfo: bool = True):
    """Input x: GRanges (cls) or a bare Seqinfo (cls == 'Seqinfo')."""
    if g.cls == "Seqinfo":
        spec = g.seqinfo_input("x")
        g.tag("seqinfo", *spec.feature_tags())
        return None, spec
    si = sample_seqinfo(g.rng, min_levels=2 if g.rng.random() < 0.5 else 1)
    x = g.ranges(si=si, mcols=g.rng.random() < 0.2)
    g.add("x", x)
    g.tag("seqinfo", *si.feature_tags())
    return x, si


def op_seqinfo_getter(g: G):
    x, si = _gr_or_seqinfo(g)
    return mk_call(g.rng, g.entry, {"x": ref("x")})


def op_set_seqlevels(g: G):
    r = g.rng
    g.begin(["pruning.mode", "value"])
    x, si = _gr_or_seqinfo(g)
    levels = list(si.levels)
    used = set(x.seqnames) if x is not None else set()
    mode = r.choice(["reorder", "rename", "rename", "add", "drop", "named", "same"])
    new = list(levels)
    if mode == "reorder":
        r.shuffle(new)
    elif mode == "rename":
        new = [lv + "_r" if r.random() < 0.6 else lv for lv in levels]
    elif mode == "add":
        extra = [c for c in CHRS if c not in levels]
        new = levels + ([r.choice(extra)] if extra else ["extra"])
        r.shuffle(new)
    elif mode == "drop" and len(levels) > 1:
        new = levels[:-1] if r.random() < 0.5 else levels[1:]
    val = cvec(new)
    if mode == "named":
        # old -> new by name, as in `seqlevels(x) <- c(chr1 = "A")`
        val = cvec([lv.upper() for lv in levels], names=levels)
    elif mode == "drop":
        pass
    ch: dict = {"x": ref("x")}
    if r.random() < 0.6:
        ch["pruning.mode"] = g.arg("pruning.mode", r.choice(["error", "coarse", "fine", "tidy"]), ["bogus", "COARSE"])
    ch["value"] = g.arg("value", val, [cvec(levels + levels[:1]), cvec([None] * len(levels)), ivec([1, 2])])
    return mk_call(r, g.entry, ch)


def op_set_seqinfo(g: G):
    r = g.rng
    g.begin(["pruning.mode", "new2old"])
    x, si = _gr_or_seqinfo(g)
    levels = list(si.levels)
    target = variant_seqinfo(r, si) if r.random() < 0.8 else SeqinfoSpec(list(levels), si.lengths, si.circular, si.genome)
    g.seqinfo_input("val", target)
    ch: dict = {"x": ref("x")}
    if r.random() < 0.85:
        new2old = [(levels.index(lv) + 1) if lv in levels else None for lv in target.levels]
        ch["new2old"] = ivec(g.arg("new2old", new2old, [[0] * len(new2old), [99] * len(new2old)]))
    if r.random() < 0.6:
        ch["pruning.mode"] = g.arg("pruning.mode", r.choice(["error", "coarse", "tidy"]), ["bogus"])
    ch["value"] = ref("val")
    return mk_call(r, g.entry, ch)


def _named_or_plain(g: G, levels: list[str], vals: list, mk, p_named: float = 0.4):
    r = g.rng
    if r.random() < p_named:
        return mk(vals, names=levels)
    return mk(vals)


def op_set_seqlengths(g: G):
    r = g.rng
    g.begin(["value"])
    x, si = _gr_or_seqinfo(g)
    L = si.levels
    vals = [r.randint(80, 400) if r.random() < 0.8 else None for _ in L]
    if any(v is None for v in vals):
        g.tag("na_arg")
    v = _named_or_plain(g, L, vals, ivec)
    ch = {"x": ref("x"), "value": g.arg("value", v, [ivec(vals + [5]), ivec([-4] * len(L)), cvec(["a"] * len(L))])}
    return mk_call(r, g.entry, ch)


def op_set_isCircular(g: G):
    r = g.rng
    g.begin(["value"])
    x, si = _gr_or_seqinfo(g)
    L = si.levels
    vals = [r.choice([True, False, False, None]) for _ in L]
    if None in vals:
        g.tag("na_arg")
    v = _named_or_plain(g, L, vals, lvec)
    ch = {"x": ref("x"), "value": g.arg("value", v, [lvec(vals + [True]), ivec([1] * len(L))])}
    return mk_call(r, g.entry, ch)


def op_set_genome(g: G):
    r = g.rng
    g.begin(["value"])
    x, si = _gr_or_seqinfo(g)
    L = si.levels
    gname = r.choice(["hg38", "mm39", "tmp"])
    vals = [gname if r.random() < 0.85 else None for _ in L]
    if None in vals:
        g.tag("na_arg")
    if r.random() < 0.3:
        v = cvec([gname])
    else:
        v = _named_or_plain(g, L, vals, cvec)
    ch = {"x": ref("x"), "value": g.arg("value", v, [cvec(vals + ["a", "b"]), ivec([1])])}
    return mk_call(r, g.entry, ch)


def op_merge(g: G):
    r = g.rng
    a = sample_seqinfo(r)
    b = variant_seqinfo(r, a) if r.random() < 0.8 else sample_seqinfo(r)
    g.seqinfo_input("x", a)
    g.seqinfo_input("y", b)
    g.tag("seqinfo", "c_merge_seqinfo", *a.feature_tags())
    dots = []
    if r.random() < 0.2:
        c = variant_seqinfo(r, a)
        g.seqinfo_input("z", c)
        dots = [ref("z")]
    return mk_call(r, g.entry, {"x": ref("x"), "y": ref("y")}, dots=dots)


def op_trim(g: G):
    r = g.rng
    si = sample_seqinfo(r, small_lengths=r.random() < 0.8, min_levels=1)
    g.tag("seqinfo", "trim_seqlengths", *si.feature_tags())
    x = g.ranges(si=si, window=(-8, 50), maxw=20)
    g.add("x", x)
    ch: dict = {"x": ref("x")}
    if r.random() < 0.4:
        ch["use.names"] = g.flag(0.1)
    return mk_call(r, g.entry, ch)


# ----------------------------------------------------------------------------- subset / c

def _subscript(g: G, n: int, names: list[str] | None):
    r = g.rng
    kinds = ["pos", "pos", "neg", "lgl", "lgl_short", "zero", "na", "na", "oob", "dbl", "names", "emptyint", "all", "ok_rep"]
    kind = r.choice(kinds)
    if kind == "pos" or (kind in ("neg", "ok_rep", "all") and n == 0):
        return ivec([r.randint(1, n) for _ in range(r.randint(0, 5))]) if n else ivec([])
    if kind == "all":
        return ivec(list(range(1, n + 1)))
    if kind == "ok_rep":
        return ivec([r.randint(1, n) for _ in range(r.randint(3, 9))])
    if kind == "neg":
        return ivec([-v for v in sorted(r.sample(range(1, n + 1), r.randint(1, min(n, 3))))])
    if kind == "lgl":
        return lvec([r.random() < 0.5 for _ in range(n)])
    if kind == "lgl_short":
        g.tag("invalid_possible")
        return lvec([r.random() < 0.5 for _ in range(max(1, n // 2))])
    if kind == "zero":
        return ivec([0] + [r.randint(1, max(n, 1)) for _ in range(2)] if n else [0])
    if kind == "na":
        g.tag("na_arg")
        return ivec([None, r.randint(1, max(n, 1))] if n else [None])
    if kind == "oob":
        return ivec([n + 1, 1] if n else [1])
    if kind == "dbl":
        return dvec([float(r.randint(1, max(n, 1)))])
    if kind == "emptyint":
        return ivec([])
    if names:
        pick = [r.choice(names + ["nope"]) for _ in range(r.randint(1, 3))]
        return cvec(pick)
    return cvec(["a", "b"])


def op_subset(g: G):
    r = g.rng
    x, ch = unary(g, names=r.random() < 0.5, mcols=r.random() < 0.4, complex_mcols=True)
    ch["i"] = _subscript(g, x.n, x.names)
    if g.cls == "GRanges" and r.random() < 0.25 and x.mcols:
        ch["j"] = r.choice([ivec([1]), cvec([x.mcols[0][0]]), lvec([True] * len(x.mcols)), cvec(["missing_col"])])
    return mk_call(r, g.entry, ch)


def op_c(g: G):
    r = g.rng
    g.begin(["y"])
    mode = r.random()
    if g.cls == "GRanges":
        si = sample_seqinfo(r)
        differ = mode < 0.6
        g.tag("seqinfo", *si.feature_tags())
        parts = []
        k = r.choice([2, 2, 2, 3])
        has_mc = r.random() < 0.3
        mc_sets = None
        for i in range(k):
            sp = si if (not differ or i == 0) else variant_seqinfo(r, si)
            n = g.size()
            rs = g.ranges(n, si=sp)
            if has_mc:
                rs.mcols = [("m1", g.mcol(n))] if r.random() < 0.85 else g.mcols_for(n, True)
            parts.append(rs)
        if differ:
            g.tag("c_merge_seqinfo")
    else:
        k = r.choice([2, 2, 2, 3])
        has_mc = r.random() < 0.3
        parts = []
        for i in range(k):
            n = g.size()
            rs = g.ranges(n)
            if has_mc:
                rs.mcols = [("m1", g.mcol(n))] if r.random() < 0.85 else g.mcols_for(n, True)
            parts.append(rs)
    names = ["x", "y", "z"]
    for nm, rs in zip(names, parts):
        g.add(nm, rs)
    ch: dict = {"x": ref("x")}
    dots = [ref(nm) for nm in names[1:len(parts)]]
    if r.random() < 0.2 and g.bad_slot is None:
        # a non-range argument among the objects
        dots.append(ivec([1]))
        g.tag("invalid")
    if r.random() < 0.25:
        ch["ignore.mcols"] = g.flag(0.1)
    return mk_call(r, g.entry, ch, dots=dots)


# ----------------------------------------------------------------------------- intra-range

def op_shift(g: G):
    r = g.rng
    g.begin(["shift"])
    x, ch = unary(g)
    big = [INT_MAX - 10, -INT_MAX + 3]
    sh = g.int_arg(x.n, -20, 20, vec_p=0.35, na_p=0.1, dbl_p=0.08)
    if x.n and r.random() < 0.05:
        sh = r.choice(big)
    ch["shift"] = g.arg("shift", sh, [ivec([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13][: x.n + 2]), cvec(["a"]), ivec([])])
    if r.random() < 0.3:
        ch["use.names"] = g.flag()
    return mk_call(r, g.entry, ch)


def op_narrow(g: G):
    r = g.rng
    g.begin(["start", "end", "width"])
    g.p_zero = 0.0 if r.random() < 0.85 else g.p_zero  # mostly legal narrowing; zero-width inputs stay possible
    x, ch = unary(g, minw=7, maxw=14)
    combo = r.choice(["start", "end", "width", "se", "sw", "ew", "se", "sw", "all3"])
    n = x.n
    if combo in ("start", "se", "sw", "all3"):
        ch["start"] = g.arg("start", g.int_arg(n, 1, 3, na_p=0.05) if r.random() < 0.9 else g.int_arg(n, -4, 5), [ivec([100]), ivec([None])])
    if combo in ("end", "se", "ew", "all3"):
        ch["end"] = g.arg("end", g.int_arg(n, 4, 6, na_p=0.05) if r.random() < 0.6 else g.int_arg(n, -3, -1), [ivec([-100])])
    if combo in ("width", "sw", "ew", "all3"):
        ch["width"] = g.arg("width", g.int_arg(n, 0, 4, na_p=0.05), [ivec([-1]), ivec([99])])
    if combo == "all3":
        g.tag("invalid")
    out = {k: ch[k] for k in ("x", "start", "end", "width") if k in ch}
    if r.random() < 0.25:
        out["use.names"] = g.flag()
    return mk_call(r, g.entry, out)


def op_resize(g: G):
    r = g.rng
    g.begin(["width", "fix"])
    x, ch = unary(g)
    n = x.n
    ch["width"] = g.arg("width", g.int_arg(n, 0, 12, na_p=0.06, dbl_p=0.05), [ivec([-3]), cvec(["a"]), ivec([1] * (n + 3))])
    if r.random() < 0.8:
        fixes = ["start", "end", "center"]
        fx = cvec([r.choice(fixes)]) if r.random() < 0.6 or n == 0 else cvec([r.choice(fixes) for _ in range(n)])
        if fx["chr"] == ["start"] and r.random() < 0.4:
            fx = "start"
        ch["fix"] = g.arg("fix", fx, ["middle", cvec(["start", "bogus"]), ivec([1])])
    if r.random() < 0.3:
        ch["use.names"] = g.flag()
    ignore_strand(g, ch, 0.4)
    return mk_call(r, g.entry, ch)


def op_flank(g: G):
    r = g.rng
    g.begin(["width", "start", "both"])
    x, ch = unary(g)
    n = x.n
    ch["width"] = g.arg("width", g.int_arg(n, -8, 10, na_p=0.06, dbl_p=0.05), [cvec(["a"]), ivec([1] * (n + 3)), ivec([])])
    if r.random() < 0.7:
        st = g.flag(0.04) if r.random() < 0.7 or n == 0 else lvec([r.random() < 0.5 for _ in range(n)])
        ch["start"] = g.arg("start", st, [cvec(["start"]), ivec([2])])
    if r.random() < 0.5:
        ch["both"] = g.arg("both", g.flag(0.04), [cvec(["yes"])])
    if r.random() < 0.2:
        ch["use.names"] = g.flag()
    ignore_strand(g, ch, 0.4)
    return mk_call(r, g.entry, ch)


def op_promoters(g: G):
    r = g.rng
    g.begin(["upstream", "downstream"])
    x, ch = unary(g)
    n = x.n
    if r.random() < 0.8:
        ch["upstream"] = g.arg("upstream", g.int_arg(n, 0, 15, na_p=0.05), [ivec([-5]), cvec(["a"])])
    if r.random() < 0.8:
        ch["downstream"] = g.arg("downstream", g.int_arg(n, 0, 10, na_p=0.05), [ivec([-1]), cvec(["a"])])
    if r.random() < 0.2:
        ch["use.names"] = g.flag()
    return mk_call(r, g.entry, ch)


# ----------------------------------------------------------------------------- inter-range

def op_range(g: G):
    r = g.rng
    x, ch = unary(g)
    dots = []
    if r.random() < 0.35:
        g.cls == "GRanges"
        y = g.ranges(si=x.seqinfo, si_ref=x.seqinfo_ref) if g.cls == "GRanges" else g.ranges()
        g.add("y", y)
        dots = [ref("y")]
    if r.random() < 0.3:
        ch["na.rm"] = g.flag(0.05)
    if r.random() < 0.5:
        ch["with.revmap"] = g.flag(0.03)
        g.tag("revmap")
    ignore_strand(g, ch, 0.4)
    return mk_call(r, g.entry, ch, dots=dots)


def op_reduce(g: G):
    r = g.rng
    g.begin(["min.gapwidth", "drop.empty.ranges", "with.revmap"])
    x, ch = unary(g, maxw=6)
    if r.random() < 0.5:
        ch["drop.empty.ranges"] = g.arg("drop.empty.ranges", g.flag(0.03), [cvec(["maybe"])])
    if r.random() < 0.6:
        ch["min.gapwidth"] = g.arg("min.gapwidth", g.int_arg(0, 0, 6, na_p=0.04, dbl_p=0.05), [-1, ivec([1, 2]), cvec(["a"])])
    if r.random() < 0.5:
        ch["with.revmap"] = g.arg("with.revmap", g.flag(0.03), [cvec(["x"])])
        g.tag("revmap")
    ignore_strand(g, ch, 0.4)
    return mk_call(r, g.entry, ch)


def op_gaps(g: G):
    r = g.rng
    g.begin(["start", "end"])
    x, ch = unary(g)
    if r.random() < 0.6:
        ch["start"] = g.arg("start", g.int_arg(0, -3, 10, na_p=0.08), [cvec(["a"]), ivec([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15])])
    if r.random() < 0.6:
        ch["end"] = g.arg("end", g.int_arg(0, 10, 70, na_p=0.08), [cvec(["a"])])
    ignore_strand(g, ch, 0.4)
    return mk_call(r, g.entry, ch)


def op_disjoin(g: G):
    r = g.rng
    x, ch = unary(g, maxw=8)
    if r.random() < 0.5:
        ch["with.revmap"] = g.flag(0.03)
        g.tag("revmap")
    ignore_strand(g, ch, 0.4)
    return mk_call(r, g.entry, ch)


# ----------------------------------------------------------------------------- set ops

def op_setop(g: G):
    r = g.rng
    g.begin(["y"], p_bad=0.04)
    x, y = g.pair(maxw=8)
    ch = {"x": ref("x"), "y": ref("y")}
    ignore_strand(g, ch, 0.4)
    return mk_call(r, g.entry, ch)


def op_parallel(g: G):
    """punion / pintersect / psetdiff: x and y must be equal length (or one of length 1)."""
    r = g.rng
    g.begin(["y", "fill.gap", "resolve.empty"], p_bad=0.08)
    n = g.size()
    mode = r.random()
    nb = n if mode < 0.85 else (1 if mode < 0.93 else max(n + r.randint(1, 3), 1))
    x, y = g.parallel_pair(n, nb, maxw=8)
    if mode >= 0.93:
        g.tag("invalid")
    ch: dict = {"x": ref("x"), "y": ref("y")}
    op = g.entry["op"]
    if op == "punion" and r.random() < 0.6:
        ch["fill.gap"] = g.arg("fill.gap", g.flag(0.03), [cvec(["a"])])
    if op == "pintersect" and g.cls == "IRanges":
        if r.random() < 0.5:
            ch["resolve.empty"] = g.arg("resolve.empty", r.choice(["none", "max.start", "start.x"]), ["bogus"])
        if r.random() < 0.4:
            ch["drop.nohit.ranges"] = g.flag(0.03)
    if op == "pintersect" and g.cls == "GRanges" and r.random() < 0.3:
        ch["strict.strand"] = g.flag()
    ignore_strand(g, ch, 0.4)
    return mk_call(r, g.entry, ch)


# ----------------------------------------------------------------------------- overlaps / nearest

def _overlap_common(g: G, ch: dict, *, with_type=True, with_gap=True, selects=False):
    r = g.rng
    typ = None
    if with_type and r.random() < 0.7:
        typ = g.arg("type", r.choice(["any", "start", "end", "within", "equal", "any", "within"]), ["bogus", "ANY", ivec([1])])
    only_one = typ in (None, "any") and g.bad_slot not in ("maxgap", "minoverlap")  # R: type "any" takes one of them
    use_gap = with_gap and r.random() < 0.5
    use_min = with_gap and r.random() < 0.5
    if only_one and use_gap and use_min:
        use_gap, use_min = (True, False) if r.random() < 0.5 else (False, True)
    if use_gap:
        ch["maxgap"] = g.arg("maxgap", g.int_arg(0, -1, 12, na_p=0.04, dbl_p=0.05), [cvec(["a"]), ivec([1, 2])])
    if use_min:
        ch["minoverlap"] = g.arg("minoverlap", g.int_arg(0, 0, 6, na_p=0.04, dbl_p=0.05), [-2, cvec(["a"])])
    if typ is not None:
        ch["type"] = typ
    if selects and r.random() < 0.8:
        ch["select"] = g.arg("select", r.choice(["all", "first", "last", "arbitrary", "all", "arbitrary"]), ["bogus", "ALL"])


def op_findOverlaps(g: G):
    r = g.rng
    g.begin(["type", "select", "maxgap", "minoverlap"])
    x, y = g.pair(derived=0.5, maxw=8)
    ch: dict = {"query": ref("x")}
    self_ov = r.random() < 0.12
    if self_ov:
        g.inputs.pop("y")
        g.tag("self")
    else:
        ch["subject"] = ref("y")
    _overlap_common(g, ch, selects=True)
    if self_ov and r.random() < 0.6:
        ch["drop.self"] = g.flag()
        ch["drop.redundant"] = g.flag()
    ignore_strand(g, ch, 0.35)
    return mk_call(r, g.entry, ch)


def op_countOverlaps(g: G):
    r = g.rng
    g.begin(["type", "maxgap", "minoverlap"])
    x, y = g.pair(derived=0.4, maxw=8)
    ch: dict = {"query": ref("x"), "subject": ref("y")}
    _overlap_common(g, ch)
    ignore_strand(g, ch, 0.35)
    return mk_call(r, g.entry, ch)


def op_overlapsAny(g: G):
    r = g.rng
    g.begin(["type", "maxgap", "minoverlap"])
    x, y = g.pair(derived=0.4, maxw=8)
    ch: dict = {"query": ref("x"), "subject": ref("y")}
    _overlap_common(g, ch)
    ignore_strand(g, ch, 0.35)
    return mk_call(r, g.entry, ch)


def op_subsetByOverlaps(g: G):
    r = g.rng
    g.begin(["type", "maxgap", "minoverlap"])
    x, y = g.pair(derived=0.4, maxw=8)
    ch: dict = {"x": ref("x"), "ranges": ref("y")}
    _overlap_common(g, ch)
    if r.random() < 0.4:
        ch["invert"] = g.flag(0.03)
    ignore_strand(g, ch, 0.35)
    return mk_call(r, g.entry, ch)


def op_nearest_like(g: G):
    r = g.rng
    g.begin(["select"])
    x, y = g.pair(derived=0.3, maxw=6)
    ch: dict = {"x": ref("x")}
    if r.random() < 0.1:
        g.inputs.pop("y")
        g.tag("self")
    else:
        ch["subject"] = ref("y")
    if r.random() < 0.7:
        legal = {"nearest": ["arbitrary", "all"], "precede": ["first", "all"], "follow": ["last", "all"]}[g.entry["op"]]
        ch["select"] = g.arg("select", r.choice(legal), ["bogus", "ALL", "first" if g.entry["op"] == "follow" else "last"])
    ignore_strand(g, ch, 0.7, 0.15)
    return mk_call(r, g.entry, ch)


def op_distance(g: G):
    r = g.rng
    g.begin(["y"], p_bad=0.06)
    n = g.size()
    mode = r.random()
    nb = n if mode < 0.85 else (1 if mode < 0.93 else n + 2)
    x, y = g.parallel_pair(n, nb, maxw=6)
    if mode >= 0.93:
        g.tag("invalid_possible")
    ch = {"x": ref("x"), "y": ref("y")}
    ignore_strand(g, ch, 0.6)
    return mk_call(r, g.entry, ch)


# ----------------------------------------------------------------------------- coverage

def op_coverage(g: G):
    r = g.rng
    g.begin(["shift", "width", "weight", "method"])
    x, ch = unary(g, maxw=10, window=(1, 50))
    n = x.n
    gr = g.cls == "GRanges"
    if r.random() < 0.4:
        if gr:
            ch["shift"] = g.arg("shift", g.int_arg(0, -3, 8, na_p=0.04), [cvec(["a"])])
        else:
            ch["shift"] = g.arg("shift", g.int_arg(n, -3, 8, vec_p=0.3, na_p=0.05), [cvec(["a"]), ivec([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14])])
    if r.random() < 0.5:
        ch["width"] = g.arg("width", g.int_arg(0, 0, 80, na_p=0.05), [-4, cvec(["a"]), ivec([1, 2, 3])])
    if r.random() < 0.4:
        if r.random() < 0.5:
            w = ivec([r.randint(0, 5) for _ in range(n)])
        else:
            w = dvec([r.randint(-8, 20) / 4 for _ in range(n)])
        if n and r.random() < 0.08:
            vals = list(w[next(iter(w))])
            vals[0] = None
            w = {next(iter(w)): vals}
            g.tag("na_arg")
        ch["weight"] = g.arg("weight", w, [cvec(["a"] * n), ivec([1] * (n + 2))])
    if r.random() < 0.3:
        ch["method"] = g.arg("method", r.choice(["auto", "sort", "hash"]), ["bogus"])
    return mk_call(r, g.entry, ch)


# ----------------------------------------------------------------------------- ordering

def _order_input(g: G, n: int | None = None, si=None, **kw):
    r = g.rng
    if g.cls == "GRanges":
        si = si or sample_seqinfo(r, min_levels=2 if r.random() < 0.7 else 1)
        g.tag("seqinfo", *si.feature_tags())
        x = g.ranges(n, si=si, strand="mixed" if r.random() < 0.7 else None, window=(0, 20), maxw=5, **kw)
        if x.n and r.random() < 0.5:  # ties
            k = r.randrange(x.n)
            x.start += [x.start[k]]
            x.width += [x.width[k]]
            x.seqnames.append(x.seqnames[k])
            if x.strand is not None:
                x.strand.append(r.choice("+-*"))
            if x.names:
                x.names.append("a")
            x.mcols = [(nm, _extend(g, col)) for nm, col in x.mcols]
        if len({s for s in x.seqnames}) > 1 and x.strand and len(set(x.strand)) > 1:
            g.tag("order_cross_seqlevel_strand")
        return x
    x = g.ranges(n, window=(0, 20), maxw=5, **kw)
    if x.n and r.random() < 0.5:
        k = r.randrange(x.n)
        x.start.append(x.start[k])
        x.width.append(x.width[k])
        if x.names:
            x.names.append("a")
        x.mcols = [(nm, _extend(g, col)) for nm, col in x.mcols]
    return x


def _extend(g: G, col: dict) -> dict:
    """Append one value to a simple mcols column (keeps lengths consistent after adding a tie)."""
    for key in ("int", "dbl", "chr", "lgl"):
        if key in col:
            return {key: col[key] + [col[key][0] if col[key] else None]} if col[key] else {key: [None]}
    if "Factor" in col:
        f = col["Factor"]
        return {"Factor": {"levels": f["levels"], "codes": f["codes"] + [None]}}
    return {"IntegerList": {"elements": col["IntegerList"]["elements"] + [[]]}}


def op_order(g: G):
    r = g.rng
    g.begin(["decreasing", "na.last", "method"])
    x = _order_input(g)
    g.add("x", x)
    dots = []
    if r.random() < 0.2:
        if g.cls == "GRanges":
            y = _order_input(g, x.n, si=x.seqinfo)
        else:
            y = _order_input(g, x.n)
        if y.n == x.n:
            g.add("y", y)
            dots = [ref("y")]
    ch: dict = {}
    if r.random() < 0.4:
        ch["decreasing"] = g.arg("decreasing", g.flag(0.03), [cvec(["up"])])
    if r.random() < 0.2:
        ch["na.last"] = g.flag(0.3)
    if r.random() < 0.15:
        ch["method"] = g.arg("method", r.choice(["auto", "radix", "shell"]), ["bogus"])
    return mk_call(r, g.entry, {}, dots=[ref("x")] + dots, mcols_kw=list(ch.items()))


def op_sort(g: G):
    r = g.rng
    g.begin(["decreasing"])
    x = _order_input(g)
    g.add("x", x)
    ch: dict = {"x": ref("x")}
    if r.random() < 0.5:
        ch["decreasing"] = g.arg("decreasing", g.flag(0.03), [cvec(["up"])])
    if r.random() < 0.2:
        ch["na.last"] = g.flag(0.3)
    ignore_strand(g, ch, 0.4)
    return mk_call(r, g.entry, ch)


def op_rank(g: G):
    r = g.rng
    g.begin(["ties.method"])
    x = _order_input(g)
    g.add("x", x)
    ch: dict = {"x": ref("x")}
    if r.random() < 0.3:
        ch["na.last"] = g.flag(0.2)
    if r.random() < 0.88:
        ch["ties.method"] = g.arg("ties.method", r.choice(["average", "first", "last", "min", "max"]), ["bogus"])
    ignore_strand(g, ch, 0.4)
    return mk_call(r, g.entry, ch)


def op_duplicated(g: G):
    r = g.rng
    x = _order_input(g)
    g.add("x", x)
    ch: dict = {"x": ref("x")}
    if r.random() < 0.3:
        ch["incomparables"] = False
    if r.random() < 0.5:
        ch["fromLast"] = g.flag(0.03)
    if r.random() < 0.1:
        ch["method"] = r.choice(["auto", "quick"])
    return mk_call(r, g.entry, ch)


def op_unique(g: G):
    r = g.rng
    x = _order_input(g)
    g.add("x", x)
    ch: dict = {"x": ref("x")}
    if r.random() < 0.3:
        ch["incomparables"] = False
    if r.random() < 0.5:
        ch["fromLast"] = g.flag(0.03)
    return mk_call(r, g.entry, ch)
