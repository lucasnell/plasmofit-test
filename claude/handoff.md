# Handoff: threads 2 and 3 on real data, 2026-09-24

*Transient. Written 16:15 while SLURM 28931 was running. Once its results are
written into `findings.md` and threads 2 and 3 in `threads.md` are updated,
**delete this file** -- it describes a run, not a conclusion.*

## What is running

SLURM **28931**, tasks 10-13, started 15:13. ETAs from per-chain iteration
rates measured at 16:03 -- the slowest of 4 chains governs each task:

| task | config | slowest chain at 16:03 | s/iter | ETA |
|---|---|---|---|---|
| 12 | `pl_wide_total0` | 1040 / 1700 | 2.88 | ~16:34 |
| 13 | `pl_wide_both` | 701 / 1700 | 4.28 | ~17:14 |
| 10 | `np_wide_bshape` | 510 / 1700 | 5.88 | ~17:59 |
| 11 | `np_wide_both` | 510 / 1700 | 5.88 | ~17:59 |

Cells: `slowest chain` is the lowest iteration reached by any of the 4 chains
out of 1700 (700 warmup + 1000 sampling); `s/iter` is wall clock since 15:13
divided by that; the ETA extrapolates the remainder at that rate. Linear
extrapolation is defensible because on task 9 warmup and sampling cost nearly
the same per iteration (3.49 vs 3.65 s). A post-adaptation speed change of
+-20% moves 17:59 by about +-20 min.

Panel therefore ~18:05. `_data/wock-fit-<config>.rds` is written before the
summary, so a job that dies late has still saved its fit.

**A prediction, not evidence.** Both `b_shape`-widened tasks run ~68% slower
per iteration than task 9 (`np_wide_total0`, ~3.5 s/iter), and `pl_wide_both`
sits between. Cost per iteration tracks leapfrog steps, which tracks
curvature, so widening `b_shape`'s prior appears to make the geometry harder
where widening `log10_total0`'s did not. That is consistent with thread 2's
simulation result, but it is a weak, indirect signal confounded with the
prior simply being wider. **It must not count toward the conclusion.** Check
it against the leapfrog and divergence counts in section 1 and leave it
there.

## The analysis is armed -- do not re-run it blindly

`_data/panel-watch.sh` runs detached (`setsid`, PPID 1, so it survives a
cleared session, a killed shell, or a dropped SSH). It waits for 28931 to
drain, logs each task's exit state and which of the eight fits are
`PRESENT`/`MISSING`, then runs `_scripts/wockner-prior-panel.R` into
`_data/prior-panel.log`.

**Check `_data/panel-watch.log` first:**

- ends at `waiting for SLURM 28931` -> jobs still running, nothing to do
- contains `panel exited 0` -> **read `_data/prior-panel.log`**; the analysis
  is done, do not run it again
- contains `panel exited` with any other code -> the panel itself failed and
  its error is in `_data/prior-panel.log`; re-run by hand
- log ends mid-way and `pgrep -f panel-watch.sh` finds nothing -> the watcher
  was killed; re-run by hand

By hand:

```
squeue -u lan68                      # empty means all four are done
srun -N 1 -n 1 -c 4 --mem=48G Rscript --vanilla \
    _scripts/wockner-prior-panel.R 2>&1 | tee _data/prior-panel.log
```

`_data/panel-watch.sh` is scaffolding: delete it once its log reads clean.
Keep the log.

## The eight-fit panel

`_scripts/wockner-fit.R` `CONFIGS` index -> what it is. Tasks 1, 2, 8, and 9
were already on disk before 28931; 10-13 are what it is fitting.

| task | config | model | prior override |
|---|---|---|---|
| 1 | `no_pool` | `no_pool` | none -- the original baseline |
| 2 | `pooled_cl` | `pooled_cl` | none |
| 8 | `np_anchor` | `no_pool` | `inoc_size` -- inoculum-anchored |
| 9 | `np_wide_total0` | `no_pool` | `sd_log10_total0 = 1` |
| 10 | `np_wide_bshape` | `no_pool` | `sd_log_b_shape = 1.5` |
| 11 | `np_wide_both` | `no_pool` | both widened |
| 12 | `pl_wide_total0` | `pooled_cl` | `sd_log10_total0 = 1` |
| 13 | `pl_wide_both` | `pooled_cl` | both widened |

Thread 2 reads `no_pool` -> `np_wide_bshape` (the `b_shape` prior in
isolation) and `np_wide_total0` -> `np_wide_both` (the same on top of a
corrected `log10_total0` prior). Thread 3 reads `no_pool` vs `pooled_cl` at
each of the three prior settings.

## How to read the output

In the header comment of `_scripts/wockner-prior-panel.R`, which is kept next
to the code so it cannot drift from it. In short: section 1 is sampler health
and voids any row that fails it, section 4 is thread 2, section 5 is thread 3.

The script's own regression check: it reproduces the pooling offset
(-0.304 h, paired mean z -3.79) and the original `no_pool`/`pooled_cl` elpd
difference (-1.38, se 1.03) from the fits that already existed. A missing fit
produces `NA` rows rather than an error, so `NA` in section 4 or 5 means a
fit is absent, not that the code is broken.

**The one interpretive trap.** On real data there is no truth, so a shift in
`cycle_length` when a prior is widened is a **shift, not a bias**. It bounds
from below how much of the reported cycle length the prior is responsible
for. It cannot show that the remainder is biological.

## Then write it down, or it is lost

The panel prints to a log and nothing else. To finish the job:

1. A new section in `claude/findings.md` -- the tables from sections 4 and 5,
   each with a `Cells:` definition per `conventions.md`, and the sampler
   health that licenses them.
2. One-line verdicts under "Settled" in `CLAUDE.md`, and the stale thread 2
   and 3 entries in `claude/threads.md` updated or struck through.
3. Delete this file and `_data/panel-watch.sh`.

**If a fit fails the R-hat < 1.05 gate**, section 1 drops it and everything
about it is void. Thread 9 found both earlier non-converged replicates were
bad *fit seeds*, not hard datasets, so try another seed before concluding
anything about the posterior. That is a ~2 h refit of that one task, so it is
a decision to take deliberately, not an automatic step.
