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

## In flight as of 2026-09-24 15:40 — read this first on coming back

SLURM **28931**, tasks 10–13, started 15:13, four real-data fits at ~1.7 h
each, so expect them around **17:00**. Logs `_data/wock-fit-1{0,1,2,3}.out`;
`_data/wock-fit-<config>.rds` is written before the summary, so a job that
dies late has still saved its fit.

**The analysis is already armed — do not re-run it blindly.**
`_data/panel-watch.sh` runs detached (`setsid`, PPID 1, so it survives a
cleared session or a dropped SSH), waits for 28931 to drain, then runs
`_scripts/wockner-prior-panel.R`. **Check `_data/panel-watch.log` first:**

- ends at `waiting for SLURM 28931` → jobs still running, nothing to do
- contains `panel exited 0` → **read `_data/prior-panel.log`**, the analysis
  is done; the watcher also logged each task's exit state and which of the
  eight fits are `PRESENT`/`MISSING`
- contains `panel exited` with anything else → the panel itself failed; its
  error is in `_data/prior-panel.log`, re-run the command below by hand
- the watcher process is gone and the log ends mid-way → it was killed;
  re-run the command below by hand

`_data/panel-watch.sh` is scaffolding: delete it once its log reads clean.
Its log is kept.

| task | config | model | prior override |
|---|---|---|---|
| 10 | `np_wide_bshape` | `no_pool` | `sd_log_b_shape = 1.5` |
| 11 | `np_wide_both` | `no_pool` | `sd_log_b_shape = 1.5`, `sd_log10_total0 = 1` |
| 12 | `pl_wide_total0` | `pooled_cl` | `sd_log10_total0 = 1` |
| 13 | `pl_wide_both` | `pooled_cl` | `sd_log_b_shape = 1.5`, `sd_log10_total0 = 1` |

**The command, if you do need to run it by hand:**

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

**How to read the output** — in order, stopping at the first failure — is
in the header comment of `_scripts/wockner-prior-panel.R`. In short: §1 is
sampler health and voids any row it fails, §4 is thread 2, §5 is thread 3.

**Do not edit `_scripts/wockner-fit.R` while 28931 is running.** `Rscript`
reads source incrementally; an edit mid-run killed two 1.5 h jobs today.

## Settled, 2026-09-24

Tables and derivations for all of these are in `claude/findings.md`.

- **The `normal(1, 0.25)` prior on `log10_total0` is misspecified.** Relaxing
  it gains **104.7 elpd (se 12.8)** and moves `R` from 6.2 into the 15–18
  range burst size implies. The inoculum anchor adds nothing over simply
  widening it (2.13 ± 2.91). **`R` is not identified by these data.**
- **The `b_shape` prior owns ~0.50 h of the simulated cycle-length bias and
  the `log10_total0` prior owns none**, so the `log10_total0`/`R` ridge is
  not the mechanism. ~1.07 h of ~1.97 h is still unexplained.
- Dead explanations: the simulated designs are *not* less informative (0.97
  of the real design's information about the oscillation), and truth rebuilt
  from posterior **draws** rather than the mean vector rescues nothing.
- The anchor switched off reproduces the pre-change posterior (width
  inflation 1.01×); both non-converged replicates were bad **fit seeds**.

## Known thin spots

- Thread 2's ~0.50 h `b_shape` effect is measured **in simulation**, three
  paired replicates. Whether it carries to real data is what 28931 asks.
- The real-data prior comparisons are n=1 by construction — one dataset, one
  fit per setting — so shifts get a Monte Carlo z, not a confidence interval.
- Three noise realizations; three posterior-draw replicates, usable only
  after a refit on a second seed. No cycle-length number is yet biological.
