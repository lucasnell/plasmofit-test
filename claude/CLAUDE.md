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
| `claude/handoff.md` | **transient** — the run in flight right now and how to finish it; delete when done |
| `claude/scripts.md` | what each script does, the `CONFIGS`/arm mechanics, `_data/` naming, how to run on the cluster |
| `claude/gotchas.md` | **before running or editing anything** — the traps that have cost hours each |
| `claude/findings.md` | the modelling results: hierarchy, pooling offset, schedule-bias simulation, nuisance priors |
| `claude/threads.md` | open threads 1–10, in priority order |
| `claude/conventions.md` | how to write in these files; **every numeric table must define its cells** |

## In flight as of 2026-09-24 16:15 — read `claude/handoff.md` first

SLURM **28931**, tasks 10–13, is fitting the last four of an eight-fit
real-data panel that answers open threads 2 and 3. ETAs: task 12 ~16:34,
task 13 ~17:14, tasks 10 and 11 ~17:59.

**The analysis is armed. Do not re-run it blindly.** `_data/panel-watch.sh`
runs detached and will run `_scripts/wockner-prior-panel.R` into
`_data/prior-panel.log` when the queue drains, ~18:05. Check
`_data/panel-watch.log` before doing anything.

**`claude/handoff.md` has all of it**: what each task is, how to tell from
the watcher log whether the analysis already ran, how to read the output,
what to do if a fit fails convergence, and — the step most easily lost —
where to write the results down. Delete it once that is done.

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
