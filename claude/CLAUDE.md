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

## Nothing in flight as of 2026-09-25

The queue is empty. SLURM 28931 (tasks 10–13) completed and the eight-fit
prior panel has been read and written up; see "the eight-fit prior panel" in
`claude/findings.md`.

**Do not edit `_scripts/wockner-fit.R` while a fit is running.** `Rscript`
reads source incrementally; an edit mid-run killed two 1.5 h jobs.

## Settled, 2026-09-25 — the eight-fit prior panel

`_scripts/wockner-prior-panel.R`, output `_data/prior-panel.log`. All eight
fits clear R-hat < 1.05; divergences run 1.05–4.0% and are the standing
caveat, worst in the two fits with both priors widened.

- **Thread 2 answered: the `b_shape` prior carries cycle length on real data
  too.** Widening it moves `cycle_length` **−0.47 h** in isolation (against
  −0.501 h in simulation) and **−1.29 h** on top of a corrected
  `log10_total0` prior, and gains **+24.8 elpd (se 2.2)**. The two nuisance
  priors are additive in elpd but **super-additive on `cycle_length`**:
  jointly −1.1 h against −0.27 h from the one-at-a-time contrasts. On real
  data this is a **shift, not a bias** — it bounds the priors' share of the
  reported 45.3 h from below and says nothing about the remainder.
- **`b_shape` is not identified in location.** Under `lognormal(2, 0.5)` its
  posterior is barely narrower than its prior (CV 0.433 against 0.533), so
  **the 14.9 reported throughout this project is a prior artefact**; widened,
  it reads 65.4 with a posterior sd of 50. Same standing as `R`.
- **Thread 3 answered: the hierarchy conclusion survives, its number does
  not.** `pooled_cl` fails to beat `no_pool` at all three prior settings
  (−1.38, −1.67, −1.68 elpd, se ~1), so correcting the misspecified prior
  does not rescue the hierarchy. But the pooling offset runs **−0.304,
  −0.110, −0.509 h** across those settings — a factor of 4.6, and
  non-monotone. **Quote −0.509 h**, the best-fitting setting; −0.304 h was
  drawn under a prior costing 104.7 elpd.

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

- The real-data prior comparisons are n=1 by construction — one dataset, one
  fit per setting — so shifts get a Monte Carlo z, not a confidence interval.
  The panel's `mean_z` says a shift beats MCMC noise and nothing more.
- Divergences are non-zero in all eight panel fits and reach 4.0% in
  `pl_wide_both`, so every posterior mean in the panel is a mean over an
  imperfectly explored posterior.
- The pooling offset (−0.509 h) is ~0.4 of the within-fit posterior sd of one
  trial's `cycle_length` (1.28 h under `no_pool` with both priors corrected).
- Three noise realizations; three posterior-draw replicates, usable only
  after a refit on a second seed. **No cycle-length number is yet
  biological** — the panel bounds the priors' share from below and stops
  there.
