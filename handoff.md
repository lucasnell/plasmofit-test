# Handoff — 2026-10-08

## State of play right now — read this first

**Two jobs running**, both submitted 2026-10-08 after 29684 read out.

```bash
cd /home2/lan68/plasmofit/plasmofit-test
squeue -u lan68
```

**SLURM 30524 — the `n_c` = 768 rung, 14 tasks, ~2–7 h.** This is the last
cheap evidence bearing on the IPM decision. Model A only.

```bash
ls _data/decay-law-a768-unit*.rds | wc -l          # expect 14
/programs/R-4.6.1/bin/Rscript --vanilla _scripts/nc-768-read.R
```

**SLURM 30525 — reseed of the production `n_c` ladder, entry 39, ~8.5 h.**
Writes with `WOCKFIT_SUFFIX=-seed2`, so it lands beside the failed fit rather
than overwriting it. **Until this lands, no reported cycle-length number
should change.** That is the single most important standing constraint.

Neither script may be edited while its job runs.

**Both jobs are inside the half-node budget**, checked 2026-10-08 against the
rule now in `CLAUDE.md`: 30524 is 14 x 1 CPU x 16G and 30525 is 4 CPUs x 24G,
so together **18 of 128 CPUs and 253,952 of 515,670 MB** — 49% of the memory
budget. Neither needs an `--array=1-N%M` throttle. Both are `PD (Priority)`
because someone else holds 250 of the node's 256 CPUs, which is a queueing
fact, not a budget violation.

### What 30524 decides, and the rule for reading it

**The full case for and against is already written**, in
`claude/ipm-decision.md`, drafted while 30524 was still running so the
reading rule could not be reverse-engineered from the answer. It covers the
floor argument below, the two concrete designs, the numerical-diffusion
hazard that makes a spectral scheme the right default for an IPM, and two
things that hold whichever way 30524 reads.

**The IPM question reduces to one inequality.** In the Erlang chain, transit
time has mean `cycle_length` and coefficient of variation **`1/√n_c`**. For a
fixed stage count the Erlang is the **minimum-variance** case — minimising
Σ1/λᵢ² subject to Σ1/λᵢ = mean puts all rates equal — so `1/√n_c` is a
**floor the family cannot go under**. Everything turns on whether the data's
preferred dispersion sits at that floor or below it.

- **At the floor** → the mesh/rate coupling is not distorting the fit. Pick
  `n_c` by elpd, and free the stage rates only if dispersion *above* the floor
  is wanted. Non-uniform rates give a hypoexponential, which changes the
  generator's diagonal and keeps `mat_exp_series` and the `matrix_exp`
  cross-check.
- **Below the floor** → the data want transit closer to deterministic than the
  chain reaches at any affordable `n_c`. The family is fighting them, and an
  IPM's free dispersion width is justified on evidence rather than taste.

**Pre-registered rule, fixed before submission** and repeated in
`_scripts/decay-law-768.sh` and in the reader. Gains so far are **+72.4**
(96→192) and **+10.4** (192→384), ratio 0.14, so geometric decay predicts
**+1.5** for 384→768.

| summed gain 384→768 | units positive | verdict |
|---|---|---|
| < +3 | < 9/14 | plateau — no rewrite justified |
| ≥ +8 | ≥ 10/14 | still climbing — IPM justified |
| anything else | | ambiguous — decide on biology |

**Why this is needed at all.** The two existing pieces of evidence disagree.
Bayesian elpd put 192 and 384 tied (−3.2, se 5.0), implying a plateau; the ML
profile from 29684 put 384 ahead by +10.4 summed in 10/14 units, implying the
preference is still drifting. Out-of-sample is the better criterion, but it is
one comparison with se 5.0. The ML profile is a fairer instrument here than
in-sample comparisons usually are, because `n_c` is a fixed structural choice
and **every rung has the same parameter count**.

### What 30525 decides

Whether the `n_c` effect on `cycle_length` survives a pinned `b_shape`. 29635
failed the gate: max R-hat 6.13 with **one** chain stuck at lp__ −929.0
against −832.3, −833.2, −836.7, while the majority agreed (`cycle_length[1]`
42.18, 42.28, 42.10, outlier 43.09). A single stuck chain is what a different
seed fixes, and reseeding changes one thing rather than three.

The deterministic tests **predict the effect persists and grows**: 0.98 h at
`b_shape` 400 against 0.72 h at 15, in all twelve cells of the least-squares
table. **If it converges and the effect is gone, the whole mechanistic account
is wrong** and should be revisited rather than patched.

**If a chain sticks again** at a similar lp gap, the mode is real: report it as
multimodality, show both modes, and move to entry 41 (`adapt_delta` 0.95,
`max_treedepth` 12, `_scripts/wockner-fit-nc-bs400-retry.sh`, still
unsubmitted) rather than reseeding a third time. That retry costs 2.5–5×
because the sampler already saturated treedepth in 46–61% of transitions.

## What 29684 decided: the decay-law test cannot tell

All 14 tasks COMPLETED, every `.err` empty, 9:07 to 21:04 elapsed. The 6–9 h
estimate was 2.3× low. Full numbers in `claude/findings.md`, "The decay law:
the test cannot tell, and it says why"; saved output
`_data/decay-law-read-2026-10-08.txt`.

- Total `d_ll` (model B minus model A) is **+1.53** over 14 units, 8–6 on
  sign. **One unit supplies +1.449**; the other 13 give +0.08.
- **The six negatives are not evidence for A.** B contains A at sigma = 0, and
  A's maximum is never at `n_c` = 96, so `ll_B >= ll_A` must hold. They are
  Nelder-Mead stopping short. Their worst, **−0.521**, is the noise floor, and
  **1 of 14 units clears it**.
- **The power bound is the result.** A maximised gain is at least the gain at
  any fixed sigma, so each unit's `d_ll` bounds from above what a linear
  component at the pre-registered sigma = 0.033 could have bought: **at most
  +0.396, median +0.0000** over 1130 observations.

So a linear component of exactly the predicted size is nearly free. The data
are **indifferent** between the laws — √ was not confirmed, only unbeaten.
An IPM therefore cannot be justified by appeal to the decay law, and a 1-D
Gaussian-kernel IPM would reproduce √ by construction. That is why the
decision moved to the variance-floor question above.

**One caveat.** The power bound needs B's optimiser to have found its maximum,
which in 6 of 14 units it did not. `fit_unit()` starts B from two fixed points
and never warm-starts it from A's solution; a third start at A's optimum with
a small sigma would make the violation impossible. The present result is
**biased toward A in magnitude**, by up to the floor.

## The by-product, which is worth more than the test

Model A alone is a **likelihood profile over `n_c` with no priors anywhere**.
Cells are differences in maximised profile log-likelihood within a unit,
summed over the 14 `grp_init` units holding all 1130 observations.

| comparison | summed | units favouring the higher `n_c` |
|---|---|---|
| 192 over 96 | **+72.4** | **14 of 14** |
| 384 over 96 | +82.8 | 14 of 14 |
| 384 over 192 | +10.4 | 10 of 14 |

The Bayesian comparison put 96 **71.8 elpd** (se 13.2) behind 192. Same
ordering, same magnitude, same plateau, from a different machine.

Two things it is **not**: not independent confirmation of the number, since an
in-sample maximised likelihood and an out-of-sample elpd difference are
different quantities and matching to 0.6 units is coincidence; and not a
uniform effect, since `MMV048_PIB|1800` and `Piperaquine|1800` give +10.2 and
+10.8 of the +72.4 while `MMV048_PartB|2800` and `OZ439|1800` give +0.30 and
+0.33. What it rules out is a **prior artefact**.

## Key decisions this session

- **Do not treat the √ decay law as established**, only unbeaten.
- **A nested model scoring worse is optimiser failure, not evidence.** This
  corrected the pre-registered reading rule, which counted six optimiser
  failures as wins for model A. The reader now asserts the inequality.
- **`n_c` = 96 is wrong without any prior.**
- **Settle the IPM question on the Erlang variance floor**, not on the decay
  law. That is what 30524 measures.
- **Reseed entry 39, not entry 41.** Entries 41-42 already carry the control
  changes, so "reseed-only run of entry 41" was a contradiction in earlier
  notes.

## Context for the next session

**What did not work, so it is not retried.**

- *Sequestration-grid discretisation* as the cause of the `n_c` shift is
  **wrong**: the duty cycle *rises* with `n_c` (0.39765 → 0.40242), implying
  `cycle_length` should rise by +0.284 h against an observed −2.17 h. Wrong
  sign, eight times too small.
- *The pre-registered reading rule for `_scripts/nc-mechanism.R`* **did not
  work**. The shift *shrinks* in every cell, because a longer window pins the
  period harder as well as accumulating more spread.
  `_scripts/nc-period-check.R` is the clean instrument for the period channel.
- *Cost extrapolation has now failed twice.* A per-leapfrog probe
  under-predicted a production fit by 2×, and extrapolating model B's cost
  from model A's under-predicted 29684 by 2.3×. Same cause: a cost model that
  holds the evaluation count fixed measures only one of two terms. The 768
  sizing avoids it by measuring the same model at both rungs — 1.925 s at 384
  against 14.894 s at 768, ratio 7.74.

**A fact worth carrying.** Observations begin at 72 h — **1.6 cycles in** —
and run to 4.8 cycles. Nothing is observed in the first cycle and a half. That
is why `b_shape`, which describes t = 0, is unidentified: it is pure backward
extrapolation. A dispersion rate would govern change *within* the window and
should be better identified. It is also why the decay-law test had no power:
√k and k barely differ over about three cycles.

**Corrections standing from earlier sessions.** The forward map accounts for
**half to two thirds** of the real `n_c` shift, not the quarter first
estimated; both are in `claude/findings.md` and the later supersedes. And
`mmcm.pdf` is **Greischar & Childs (2023)**, *Trends in Parasitology* 39(8),
doi 10.1016/j.pt.2023.05.006 — an earlier note transposed the authors from the
adjacent 2019 entry.

**The two PDFs at the repo root are untracked**, deliberately: they are
published articles and this repo is public. A fresh clone will not have them.

**Almost everything in `_data/` is gitignored**, so the fits — 726 files,
5.2 GB — exist only on the cluster filesystem. The exception, added
2026-10-08: `decay-law-read-*.txt` and `nc-768-read-*.txt` are negated in
`_data/.gitignore` and tracked, because numbers from them are quoted in the
notes and have to be re-derivable from the repo alone.
**Where things live.** Root `CLAUDE.md` (stable context and settled
decisions), `PROJECT_INDEX.md` (status, workstreams, decision log), `TODO.md`
(actionable layer), this file. The detail is in `claude/`: `findings.md`
(results, every table's cells defined), `gotchas.md` (**read before running
anything**), `scripts.md`, `threads.md` (long-form behind `TODO.md`), and
`references.md`.
