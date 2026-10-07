# plasmofit — malaria within-host dynamics — index

**Status as of 2026-10-07:** The modelling is in a late diagnostic phase, not
a writing phase. Three large questions are settled — the `log10_total0` prior
was misspecified, `b_shape` must be fixed rather than estimated, and the
hierarchy earns nothing either way — but **no cycle-length number is
reportable yet**. Two live issues block it: the production `n_c` ladder failed
to converge, and the age-structure model conflates a numerical mesh with a
biological rate. The decay-law test (SLURM 29684) is running and decides
whether the structural fix is worth building.

## The question

How long is the asexual blood-stage replication cycle of *P. falciparum* in
human volunteer infection studies, and — the part that has taken the work —
how much of any estimate is an artefact of modelling choices nobody checked?

The raw estimate has moved between 38.7 h and 45.5 h across settings that were
chosen for numerical convenience or left at defaults. Establishing which
number to believe, and with what qualifications, is the deliverable.

## Approach

A hierarchical Bayesian ODE model (Stan, via the `plasmofit` package) of
stage-structured parasite dynamics. Parasites advance through `n_c` sequential
exponential compartments, so transit over one cycle is Erlang; late stages
sequester and leave circulation, and only circulating parasites are observed
as iRBC/mL. Four model variants differ in what is pooled across trials
(`no_pool`, `pooled_R`, `pooled_cl`, `pooled_both`).

Two instruments run alongside the real-data fits:

1. **A schedule simulation** that generates data from the fitted model and
   refits, measuring estimator bias against a known truth. This is what showed
   the +1.97 h cycle-length bias and attributed it.
2. **Deterministic forward-map tests**, noiseless and prior-free, which
   isolate what is built into the model's arithmetic from what comes from
   estimation.

Model comparison is by PSIS-LOO on observations, with held-out elpd where a
mask is appropriate. The convergence gate is max R-hat < 1.05.

## Current state

| Workstream | State | Next |
|---|---|---|
| Cycle-length bias attribution | **Attributed** to `b_shape` (−0.92 h of +1.97 h); `fix_bshape` stands at n_pair = 1 | Nothing. Leans on `wide_bshape`, 6 clean pairs |
| `n_c` as a biological assumption | **Established** that 96 is wrong by 71.8 elpd, structurally; mechanism explained | Production ladder retry, written and unsubmitted |
| Decay law (√ vs linear) | **Running**, SLURM 29684, 14 tasks, 6-9 h | Read with `_scripts/decay-law-read.R` |
| Age-structure rewrite (IPM) | **Gated** on the decay law | Do not start before it reads out |
| Hierarchy / model comparison | **Closed** — no robust evidence either way | Nothing unless a design argument changes |
| Per-individual initial density | **Scoped**, not estimable as a free effect | Blocked on subject weights |
| Package | Pushed, `8dde0c1`; generated bindings now match the Stan sources | — |

## Key links

- Package repo: `github.com/lucasnell/plasmofit` (`/home2/lan68/plasmofit/plasmofit`)
- Analysis repo: `github.com/lucasnell/plasmofit-test` (this folder)
- Data: `_data/wockner-cleaned.csv`
- Results in full: `claude/findings.md`
- Prior art on the same confound: Greischar & Childs, `mmcm.pdf` / `mmc1.pdf`,
  and `claude/references.md`

## Decision log

Append-only. Supersede a line with a new one rather than editing it.

| Date | Decision | Why |
|---|---|---|
| 2026-09-24 | Relax the `log10_total0` prior to sd 1 | `normal(1, 0.25)` is misspecified; relaxing gains 104.7 elpd (se 12.8) and moves `R` into the range burst size implies |
| 2026-09-24 | Drop the inoculum anchor as a route to `R` | It adds nothing over simply widening the prior (2.13 ± 2.91); `R` is not identified by these data |
| 2026-10-01 | Fix `b_shape` rather than estimate it | +28 to +36 elpd at every rung, 10–12 standard errors |
| 2026-10-02 | Pin `b_shape` at 400 | The ladder plateaus there; `cycle_length` is flat from 250 upward |
| 2026-10-02 | Report any cycle length with the `b_shape` ladder | The estimate is conditional on the pin: 44.0 h at 50, 43.6 h at 250 |
| 2026-10-04 | Withdraw the claim that the hierarchy earns its keep | Design A is mask-dependent; the signal sat on training-set points |
| 2026-10-04 | Masks use `rank(time)`, never `arrange()` | Reordering rows changes derived group factor levels and invalidates comparison with every earlier fit |
| 2026-10-05 | Close the horizon ladder as flat | Its one signal came from a fit that failed to converge; reseeded, no rung reaches \|z\| = 2 |
| 2026-10-05 | Park SBC | It was the branch for `fix_bshape` returning near zero, which it did not; it checks code that has not changed |
| 2026-10-06 | Do not reseed `fix_bshape` a third time | Two seeds reproduced the same biases and both stalled just above the gate; pre-registered stopping rule |
| 2026-10-06 | Treat `n_c` as a biological assumption, not a numerical setting | 96 predicts 71.8 elpd worse than 192, and `max_rel_diff` of 8e−13 rules out an arithmetic cause |
| 2026-10-07 | Gate any IPM rewrite on the decay-law test | A 1-D IPM with a Gaussian kernel reproduces the √ law it would replace, so the law must be settled first |
| 2026-10-07 | Hold all cycle-length numbers until the production `n_c` ladder converges | The failed run leaves open whether the effect survives a pinned `b_shape` |
