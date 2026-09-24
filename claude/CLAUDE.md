# plasmofit-test — scripts and fits

Driver scripts and saved fits for the `plasmofit` package; not itself a
package. **Working on the cluster**, primary directory
`/home2/lan68/plasmofit/plasmofit-test`, package source
`/home2/lan68/plasmofit/plasmofit`. Those are **two separate git
repositories** and both need committing separately; package-side notes live
in `plasmofit/claude/` (`chat3/CLAUDE.md` for model internals,
`chat5/CLAUDE.md` for the inoculum-anchored prior). These scripts used to
live in `plasmofit/_testing/`, so old commits still name that path.

This file is an index and the current state. **Read the file you need:**

| file | when to read it |
|---|---|
| `claude/scripts.md` | what each script does, the `CONFIGS`/arm mechanics, `_data/` naming, how to run on the cluster |
| `claude/gotchas.md` | **before running or editing anything** — the traps that have cost hours each |
| `claude/findings.md` | the modelling results: hierarchy, pooling offset, schedule-bias simulation, nuisance priors |
| `claude/threads.md` | open threads 1–10, in priority order |
| `claude/conventions.md` | how to write in these files; **every numeric table must define its cells** |

## In flight as of 2026-09-24 15:20 — read this first on coming back

SLURM **28931**, tasks 10–13, started 15:13, four real-data fits at ~1.7 h
each, so expect them around **17:00**. Logs `_data/wock-fit-1{0,1,2,3}.out`;
`_data/wock-fit-<config>.rds` is written before the summary, so a job that
dies late has still saved its fit.

| task | config | model | prior override |
|---|---|---|---|
| 10 | `np_wide_bshape` | `no_pool` | `sd_log_b_shape = 1.5` |
| 11 | `np_wide_both` | `no_pool` | `sd_log_b_shape = 1.5`, `sd_log10_total0 = 1` |
| 12 | `pl_wide_total0` | `pooled_cl` | `sd_log10_total0 = 1` |
| 13 | `pl_wide_both` | `pooled_cl` | `sd_log_b_shape = 1.5`, `sd_log10_total0 = 1` |

**The one command to run when they land:**

```
squeue -u lan68                      # empty means all four are done
srun -N 1 -n 1 -c 4 --mem=48G Rscript --vanilla \
    _scripts/wockner-prior-panel.R 2>&1 | tee _data/prior-panel.log
```

`wockner-prior-panel.R` answers threads 2 and 3 off the resulting eight-fit
panel and needs no editing. It was run end to end against the four fits that
existed at 15:20, so it is known to work; a missing fit becomes `NA` rows
rather than an error, which is also how to see at a glance which of the four
failed. It reproduces the pooling offset (−0.304 h, paired mean z −3.79) and
the original `no_pool`/`pooled_cl` elpd difference (−1.38, se 1.03), which is
the regression check on it.

Read the output in order and stop at the first failure:

1. **§1 sampler health.** Any row with max R-hat ≥ 1.05 is dropped and
   everything about it is void. Thread 9's lesson: that is more likely a bad
   fit seed than a hard posterior, so re-run the task before concluding.
2. **§4 thread 2.** `no_pool → np_wide_bshape` is `b_shape`'s prior in
   isolation; `np_wide_total0 → np_wide_both` is the same on top of the
   corrected `log10_total0` prior. In simulation, widening `b_shape` moves
   `cycle_length` −0.50 h and widening `log10_total0` moves nothing. On real
   data there is no truth, so that row is a **shift, not a bias**: it bounds
   the prior's contribution from below and cannot show the remainder is
   biological.
3. **§5 thread 3.** `no_pool` vs `pooled_cl` at three prior settings. The
   last table is the number the hierarchy conclusion rests on — the
   cycle-length pooling offset in hours. If it moves, every hierarchy
   conclusion was drawn under a prior costing 104.7 elpd and needs restating.

**Do not edit `_scripts/wockner-fit.R` while 28931 is running.** `Rscript`
reads source incrementally; an edit mid-run killed two 1.5 h jobs today.

## Settled, 2026-09-24

- **The `normal(1, 0.25)` prior on `log10_total0` is misspecified.** Relaxing
  it gains **104.7 elpd (se 12.8)** and moves `R` from 6.2 into the 15–18
  range burst size implies. The inoculum anchor adds nothing over simply
  widening it (2.13 ± 2.91). **`R` is not identified by these data.**
- **The `b_shape` prior owns ~0.50 h of the simulated cycle-length bias and
  the `log10_total0` prior owns none**, so the `log10_total0`/`R` ridge is
  not the mechanism. ~1.07 h of ~1.97 h is still unexplained.
- The simulated designs are **not** less informative than the real one (0.97
  of its information about the oscillation). Thread 3 as originally posed is
  dead.
- Truth rebuilt from posterior **draws** rather than the mean vector does not
  rescue nuisance recovery, and the bias survives. That explanation is dead.
- The inoculum anchor switched off targets the same posterior as the
  pre-change code: width inflation 1.01×, 0 of 208 means beyond 3 MCSE.
- Both non-converged replicates were bad **fit seeds**, not hard datasets.

Detail and the tables behind each of these are in `claude/findings.md`.

## Known thin spots

- Thread 2's ~0.50 h `b_shape` effect is measured **in simulation**, three
  paired replicates. Whether it carries to real data is what 28931 asks.
- The real-data prior comparisons are n=1 by construction — one dataset, one
  fit per setting — so shifts get a Monte Carlo z, not a confidence interval.
- Three noise realizations, and three posterior-draw replicates (all usable
  only after a refit on a second seed).
- No cycle-length number should yet be reported as biological.
