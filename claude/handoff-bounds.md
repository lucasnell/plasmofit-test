# Handoff: ladder extension and bound geometry, 2026-10-01

*Transient. Delete once the results are in `findings.md` and thread 2 is
updated.*

## What is running

| job | what | tasks | outputs |
|---|---|---|---|
| **29493** | `b_shape` ladder extension, `wockner-fit.R` 19-21 | 3 | `_data/wock-fit-np_bs{400,600}.rds`, `np_bs250_ms1000`, + `-LOO-` |
| **29494** | bound geometry, two simulation arms | 6 (36-38, 41-43) | `_data/wock-schedsim-{fit,RES}-cl_{move,wide_move}-rep{1,2,3}.rds` |

Descriptive table, no numeric cells. Both mail `lan68@cornell.edu` on
`END,FAIL`. Real-data fits ~2 h; the simulation arms ~2 h, `cl_wide_move`
longer because a lower `min_cl` widens the Erlang window.

## 29493: where does the ladder stop climbing?

elpd rose monotonically across the whole ladder and was still rising at the
250 rung (+36.2 against a free `b_shape`, +2.26 against the 100 rung), so 250
was the top of what was run, not an optimum. 400 and 600 bracket it.

**`np_bs250_ms1000` is the control and must be read first.** 400 and 600 need
`max_shape` above them, so they use 1000 where rungs 14-18 used 400 --
changing two things at once. `max_shape` also sets the scale of Stan's
bounded transform, so it can move divergences and step size by itself, and
*divergences rising with the pin* is one of the ladder's findings. That rung
is `np_bs250` with `max_shape` alone changed:

- `np_bs250_ms1000` vs `np_bs250` = what `max_shape` is worth on its own.
- Anything in 400/600 beyond that is the pin.

If the control moves elpd or divergences appreciably, the extended rungs
cannot be compared with 14-18 without subtracting it.

### How to read it

Re-run `_scripts/wockner-bshape-ladder.R` after adding the three labels to
`LADDER`. Expect from the existing trend, as a prediction rather than a
result: `cycle_length` keeps falling (43.6 h at 250), divergences keep
rising (4.58% at 250), and elpd keeps rising but by less each rung. **The
question is whether elpd plateaus.** An age range of 5.2 h at 250 falls to
4.1 h at 400 and 3.4 h at 600, against 0.47 h per stage cell, so this is not
yet asking the model to resolve below its own discretisation -- if elpd is
still climbing at 600, that is a real preference for near-perfect synchrony
and not a discretisation artefact.

## 29494: bound geometry, thread 2's last cheap suspect

The ~1.5 h of simulated cycle-length bias that survives correcting both
nuisance priors is unexplained. The default window is `[35, 50]` against a
truth of 45.012 h -- 5 h above, 10 h below -- and the `cycle_length` prior is
normal on the **logit** scale between those bounds, so an asymmetric window
is also an asymmetric prior in hours.

| arm | bounds | what it isolates |
|---|---|---|
| `cl_move` | `[37.5, 52.5]` | asymmetry, at the same 15 h width |
| `cl_wide_move` | `[30, 60]` | distance from any bound (moved *and* widened) |

Both **move** the bounds rather than only widening them: `max_cl = 55` is
what reintroduced the boundary mode once already.

**Read them together.** Only `cl_wide_move` moving the bias means it is
distance from the bounds. Both moving it means it is the asymmetry. Neither
moving it means bound geometry is not the remaining explanation, the ~1.5 h
is somewhere else, and thread 2 is out of cheap suspects.

**The confound to state either way**: `cl_prior_center` stays at 48 h so the
prior's *location* in hours is unchanged, but its *width* in hours is not --
the same `sd_logit_cl` spans more hours in a wider window. Bounded by the
cycle-length prior's measured weight of ~0.15 h, so it cannot account for a
large move, but it is not zero and `cl_wide_move` carries more of it than
`cl_move`.

**These arms use `build`, not `data`.** `min_cl`/`max_cl` feed
`mean_logit_cl` and the Erlang window, so a post-hoc override of the built
list would leave every derived field describing the old bounds and produce a
wrong fit with no error. `wockner-schedule-sim.R` now builds a second
`archer_stan_data()` for any arm with `build`, with a `stopifnot` that the
design is unchanged and the truth sits inside the fitted bounds.

### How to read it

`_scripts/wockner-schedule-sim-analyze.R` globs the RES files, so both arms
appear with no argument. They are reps 1-3, which are the replicates with a
converged `default` on the same simulated dataset, so each pairs against it.
The number is the paired change in cycle-length bias, as for every other arm
in the budget. Note the budget is at **n=4** for the existing arms and these
will be **n=3**.

## Outstanding, not started

- **`default-rep4`** still misses the R-hat gate at 1.06 on a second seed --
  a near-miss (divergences 478 -> 39, ESS 13 -> 117), not a hard dataset. A
  third seed or a longer warmup would likely clear it and take the paired
  budget to n=5. Deliberately skipped this round.
- **Thread 4, `hold_out` masking then Design A.** Still the only honest route
  to the hierarchy answer, and still the real engineering job: observation
  level `loo` is the weak comparison and trial-level `loo` is broken (Pareto
  k > 0.7 for all 13 units, both models). All-zero `hold_out` reproducing
  current fits exactly is the regression test.

## When these land

1. `sacct -j 29493,29494 --format=JobID,State,ExitCode -P`.
2. Ladder: add the three labels to `LADDER` in
   `_scripts/wockner-bshape-ladder.R` and re-run; read the `ms1000` control
   before the new rungs.
3. Bounds: re-run `_scripts/wockner-schedule-sim-analyze.R`.
4. Write into `findings.md` with `Cells:` per `conventions.md`, update
   `CLAUDE.md` and `threads.md`, and delete this file.
