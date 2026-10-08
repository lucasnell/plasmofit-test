# Handoff — 2026-10-08

## State of play right now — read this first

**Nothing is running.** SLURM 29684, the decay-law test, finished and has been
read. The queue is empty and nothing is pending submission.

**One job should be submitted next**, and until it lands no reported
cycle-length number should change. That is the single most important standing
constraint. See "The one thing to run next" below.

## What 29684 decided: the test cannot tell, and it says why

All 14 tasks COMPLETED, every `.err` empty, 9:07 to 21:04 elapsed each. The
6-9 h estimate was **2.3x low** — model B at `n_c` = 384 was extrapolated from
model A's cost, and extrapolating one model's cost from another's does not
work, because the extra parameter changes the *number* of optimiser
evaluations as well as the cost of each. Recorded in `claude/gotchas.md`.

Full numbers: `claude/findings.md`, section "The decay law: the test cannot
tell, and it says why". Saved output `_data/decay-law-read-2026-10-08.txt`,
regenerate with `/programs/R-4.6.1/bin/Rscript --vanilla
_scripts/decay-law-read.R`.

The question was whether synchrony decays as **√(cycles)** — what the Erlang
chain and a 1-D Gaussian-kernel IPM both give — or **linearly in cycles**,
which is what between-parasite heterogeneity in cycle duration gives.

- Total `d_ll` (model B minus model A) is **+1.53** over 14 units, B ahead in
  8 and A ahead in 6. **One unit supplies +1.449 of that**; the other 13
  together give +0.08.
- **The six negatives are not evidence for A.** B contains A at sigma = 0, and
  A's maximum is never at `n_c` = 96, so both maximise over the same
  {192, 384} and `ll_B >= ll_A` must hold. Those six are Nelder-Mead stopping
  short. Their worst, **−0.521**, is the noise floor, and **1 of 14 units
  clears it**.
- **The power bound is the real result.** A maximised gain is at least the
  gain at any fixed sigma, so each unit's `d_ll` bounds from above what a
  linear component of the pre-registered size — sigma = 0.033 — could have
  bought. Outside the one unit above the floor that bound is **at most +0.396,
  median +0.0000**, over 1130 observations. Five units' optimisers landed
  within 0.01 of 0.033 and gained +1.449, +0.097, +0.069, +0.069, and +0.053.

So a linear component of exactly the size the project had reason to expect is
nearly free in likelihood terms. **The data are indifferent between the two
laws.** They do not prefer √; √ was simply not beaten.

**Consequence for the IPM (thread 15).** The branch that would have killed a
Gaussian-kernel IPM did not fire, and neither did the branch that would have
endorsed one. An IPM **cannot be justified on the grounds that √ is
established**, and a 1-D IPM with a Gaussian development kernel would
reproduce √ by construction anyway. The case for it rests entirely on the
structural argument — `n_c` is simultaneously the numerical mesh and the
biological desynchronisation rate — which is a biology-and-design judgement.
That decision is Lucas's, and it is weighed against discarding the validated
Erlang-window series, its `matrix_exp` cross-check, and the comparability of
every existing fit.

**A rerun would not fix this.** Warm-starting B from A's solution (see below)
would remove the noise floor, but it cannot create power the likelihood
surface does not have, and it costs the full ~21 h wall again. Sharpening the
design — a longer observation window, units with more cycles — would do more.

## The by-product, which is worth more than the test

Model A alone is a **likelihood profile over `n_c` with no priors anywhere**.
Cells are differences in maximised profile log-likelihood within a unit,
summed over the 14 `grp_init` units that hold all 1130 observations.

| comparison | summed | units favouring the higher `n_c` |
|---|---|---|
| 192 over 96 | **+72.4** | **14 of 14** |
| 384 over 96 | +82.8 | 14 of 14 |
| 384 over 192 | +10.4 | 10 of 14 |

The Bayesian comparison on the same observations put 96 **71.8 elpd** (se
13.2) behind 192, with 192 and 384 tied. Same ordering, same magnitude, same
plateau, from a different inferential machine with no priors.

Two things this is **not**. It is not independent confirmation of the number —
an in-sample maximised likelihood and an out-of-sample elpd difference are
different quantities, and agreeing to 0.6 units is coincidence. And the gains
are concentrated: `MMV048_PIB|1800` and `Piperaquine|1800` give +10.2 and
+10.8 of the +72.4 while `MMV048_PartB|2800` and `OZ439|1800` give +0.30 and
+0.33. What it does rule out is a **prior artefact**: `n_c` = 96 being
structurally wrong no longer depends on any part of the Bayesian
specification.

## The one thing to run next

Resubmit the production `n_c` ladder. 29635 failed the convergence gate (max
R-hat 6.13 and 8.34, one chain stuck 96 and 70 lp units below the others) and
was not reported.

**Submit a reseed-only version of entry 41 (`n_c` = 192) first, ~9 h.**
`_scripts/wockner-fit-nc-bs400-retry.sh` is written and **deliberately not
submitted**: it adds `adapt_delta` 0.95 and `max_treedepth` 12, which costs
2.5-5x because the sampler already saturated treedepth in 46-61% of
transitions. Entry 41 would run 21-42 h and entry 42 would likely exceed its
own 3-day walltime. Reseeding alone targets the observed failure — one stuck
chain — at a fifth of the cost.

**What it decides.** Whether the `n_c` effect on `cycle_length` survives a
pinned `b_shape`. The deterministic tests *predict it should persist and be
larger*: 0.98 h at `b_shape` 400 against 0.72 h at 15, in all twelve cells of
the least-squares table. **If it converges and the effect is gone, the whole
mechanistic account is wrong** and should be revisited rather than patched.

## Key decisions this session

- **Do not treat the √ decay law as established.** 29684 could not distinguish
  it from linear at the amplitude that matters.
- **A nested model scoring worse is an optimiser failure, not evidence.** This
  corrected the pre-registered reading rule, which had counted six optimiser
  failures as wins for model A. The nesting was stated in
  `decay-law-read.R`'s own header from the start and the rule did not use it.
  The reader now asserts it.
- **`n_c` = 96 is wrong without any prior**, from the ML profile.

## Context for the next session

**What did not work, so it is not retried.**

- *Sequestration-grid discretisation* as the cause of the `n_c` shift is
  **wrong**: the duty cycle *rises* with `n_c` (0.39765 → 0.40242), implying
  `cycle_length` should rise by +0.284 h against an observed −2.17 h. Wrong
  sign, eight times too small.
- *The pre-registered reading rule for `_scripts/nc-mechanism.R`* — separate
  the channels by whether the shift is constant or growing with window length
  — **did not work**. The shift *shrinks* in every cell, because a longer
  window pins the period harder as well as accumulating more spread.
  `_scripts/nc-period-check.R` is the clean instrument for the period channel.
- *A per-leapfrog cost probe* under-predicted the production fit by 2x, and
  *extrapolating model B's cost from model A's* under-predicted 29684 by 2.3x.
  Both have the same cause: a cost model that holds the evaluation count fixed
  measures only one of the two terms.

**If the decay-law test is ever rerun**, `fit_unit()` in
`_scripts/decay-law-test.R` starts model B from two fixed points and never
warm-starts it from model A's solution. A third start at A's optimum with a
small sigma makes the nesting violation impossible, at the cost of one extra
optimisation. This is also why the present result is **biased toward A in
magnitude**, by up to the noise floor.

**A fact worth carrying.** Observations begin at 72 h — **1.6 cycles in** —
and run to 4.8 cycles. Nothing is observed in the first cycle and a half. That
is why `b_shape`, which describes t = 0, is unidentified: it is pure backward
extrapolation. A dispersion rate would govern change *within* the window and
should be better identified.

**Corrections standing from earlier sessions.** The forward map accounts for
**half to two thirds** of the real `n_c` shift, not the quarter first
estimated from the period channel alone; both figures are in
`claude/findings.md` and the later one supersedes. And `mmcm.pdf` is
**Greischar & Childs (2023)**, *Trends in Parasitology* 39(8), doi
10.1016/j.pt.2023.05.006 — an earlier note transposed the authors from the
adjacent 2019 entry.

**The two PDFs at the repo root are untracked**, deliberately: they are
published articles and this repo is public. A fresh clone will not have them.

**Where things live.** Root `CLAUDE.md` (stable context and settled
decisions), `PROJECT_INDEX.md` (status, workstreams, decision log), `TODO.md`
(actionable layer), this file. The detail is in `claude/`: `findings.md`
(results, every table's cells defined), `gotchas.md` (**read before running
anything**), `scripts.md`, `threads.md` (long-form behind `TODO.md`), and
`references.md`.
