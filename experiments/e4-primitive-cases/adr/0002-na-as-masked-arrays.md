# ADR-0002: NA in integer results is a masked array

- Status: accepted.
- `select="first"` without a hit, `nearest`/`precede`/`follow` without a
  subject, `distance` across seqnames: `np.ma.MaskedArray`, mask = NA.
- Projection maps masked -> `null`, matching golden `null`.
- Why: numpy ints have no NA; sentinel values (-1) collide with valid
  negative indices.
