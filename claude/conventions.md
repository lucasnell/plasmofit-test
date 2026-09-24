# Writing conventions for these notes

*Split out of `CLAUDE.md`, which is now an index. See `CLAUDE.md` for
orientation and the current state of play.*

**Every table must say what is in its cells.** This applies to every file in
`claude/`. A number here is almost never
self-describing: `8.89` could be a posterior mean, a median, a single group's
value, or a mean over groups, and `-40%` could be relative to a truth, to
another arm, or to a prior. State it once, immediately above the table, in a
line or two. Cover, as applicable:

- **what the number is** -- posterior mean, median, mean over groups, sum,
  difference, ratio
- **its units** -- hours, log10 units, dimensionless, percent of what
- **what it is relative to**, if it is a comparison -- which baseline, and
  which direction is "better" or "less biased"
- **what it is aggregated over**, if anything -- groups, trials, replicates,
  draws -- and how many
- **whether it is paired**, when replicates share a dataset, since paired and
  unpaired numbers of the same quantity differ here by more than the effects
  being measured

This applies to any table of **numbers**. Purely descriptive tables -- the
script list, the `_data/` naming key, a list of which fit came from which code
and seed -- are exempt, since their cells are prose.

Do not rely on the column header alone. A header names the quantity; it does
not say how it was computed, and these notes' whole value is that a number in
it can be re-derived a month later.

The same applies to a number quoted in prose. Write "posterior mean over 14
`grp_init` groups" rather than "the estimate", and give the script and the
saved output it came from.
