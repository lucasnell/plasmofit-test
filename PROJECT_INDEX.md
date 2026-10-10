# plasmofit — malaria within-host dynamics — index

**Status as of 2026-10-09:** The age-structure question is **answered**, and
the answer is that the data cannot separate initial synchrony from accumulated
desynchronisation — they constrain only the total spread at the end of the
window. That closes the IPM question (do not build it), explains why the decay
law was not learnable, and leaves `cycle_length` **stable at 40.93–41.15 h once
`n_c` >= 384**. The first converged production fit with `b_shape` pinned gives
**41.20 h** (13 trials, `n_c` = 192, SLURM 30525), but at 5.8% divergences and
with **no second rung**: entry 40 (`n_c` = 384) failed again on its reseed
(SLURM 30829, 2026-10-10, max R-hat 7.74, all four chains at different `lp__`
levels). The ladder is still open and no headline number is final. **SLURM
30830 is running**: the tuned retry of that rung (entry 42) and a `b_shape`
ladder at `n_c` = 384 (entries 43–48), since the 400 pin was chosen at
`n_c` = 96 and may not carry.

The forward-map question is also closed: the convolution was added to the
package and removed again, slower and less accurate than the series it would
have replaced. Nothing in the package changed today.

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
| Production `n_c` ladder | **One rung converged** (30525, `n_c` = 192, 41.20 h, max R-hat 1.043) at 5.8% divergences | Entry 40 failed twice (29635; reseed 30829, 2026-10-10, general mixing failure, 34% treedepth saturation). **SLURM 30830 running**: tuned retry (42) plus `b_shape` ladder at 384 (43–48) |
| `n_c` as a biological assumption | **Resolved**: 96 and 192 are wrong, and `cycle_length` is stable at 40.93–41.15 h once `n_c` >= 384. The `b_shape`/rung trade-off is the mechanism (Spearman −1.00 in 14/14 units) | Second production rung at `b_shape` 400 failed twice; the 400 pin was chosen at `n_c` = 96 and may not carry to 384 |
| Decay law (√ vs linear) | **Indecisive, and now explained**: the total spread is nearly fixed across rungs, so its growth law is barely constrained | Closed |
| Age-structure rewrite (IPM) | **Closed — do not build it.** Both profiles saturate, so an IPM returns `sigma_d` at a boundary rather than a rate | Nothing |
| Forward map | **Closed.** The convolution was added to the package and removed again (`bba580f`): slower and less accurate than the Erlang-window series at every usable `n_c` | Nothing; do not re-open without new evidence |
| Hierarchy / model comparison | **Closed** — no robust evidence either way | Nothing unless a design argument changes |
| Per-individual initial density | **Scoped**, not estimable as a free effect | Blocked on subject weights |
| Package | Pushed, `8dde0c1`; generated bindings now match the Stan sources | — |

## Key links

- Package repo: `github.com/lucasnell/plasmofit` (`/home2/lan68/plasmofit/plasmofit`)
- Analysis repo: `github.com/lucasnell/plasmofit-test` (this folder)
- Data: `_data/wockner-cleaned.csv`
- Results in full: `claude/findings.md`
- Age-structure design decision and the validated prototype: `claude/ipm-decision.md`
- Prior art on the same confound: Greischar & Childs (2023) *Trends
  Parasitol* 39(8). PDFs at the repo root as `mmcm.pdf` / `mmc1.pdf`,
  **untracked** (published, and this repo is public), so a clone will not
  have them. Summary in `claude/references.md`.

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
| 2026-10-08 | Do not treat the √ decay law as established | 29684 could not distinguish it from linear: at the predicted sigma = 0.033 a linear component buys at most +0.396 and typically +0.0000 log-likelihood units over 1130 observations |
| 2026-10-08 | A nested model scoring worse counts as optimiser failure, not evidence | B contains A at sigma = 0, so six negative `d_ll` values are impossible as evidence; their worst, 0.521, is the noise floor |
| 2026-10-08 | Settle the IPM question on the Erlang variance floor, not on the decay law | Transit CV is `1/√n_c` and Erlang is the minimum-variance case at fixed stage count, so the only thing a chain cannot do is go below that floor |
| 2026-10-08 | Reseed entry 39, not entry 41 | 41-42 already carry `adapt_delta` 0.95 and `max_treedepth` 12, so a reseed-only run means 39; earlier notes said 41 and were wrong |
| 2026-10-08 | Never take more than half the node | `cbsugreischar` is one shared node (256 CPUs, 1,031,340 MB); arrays that would exceed 128 CPUs or 515,670 MB get `--array=1-N%M`. Rule and the formula for M are in `CLAUDE.md` |
| 2026-10-08 | Judge the dispersion profile on whether it is bounded in the fine direction, not on an interior maximum | A saturating profile never has an interior maximum, and saturation means the data do not exclude zero dispersion — which argues against an IPM rather than for it |
| 2026-10-08 | The exact chain is a convolution, not only a matrix exponential | In absolute developmental age, transport is pure-birth, growth and sequestration are weights; reproduces `mat_exp_series` to 1e-12 and runs 38x faster at `n_c` = 384, so cost is no longer a reason to rewrite |
| 2026-10-08 | An IPM is a different model, not a reparameterisation | The stage at time t is exactly Poisson, a lattice distribution; any continuous kernel differs in the tails, worth 0.06 log10 units at the deepest trough and growing with `n_c` |
| 2026-10-08 | Keep 30527 running although 30576 supersedes it | 30527 is the production `mat_exp_series` path and the convolution's cross-check; if the two ever disagree, believe `mat_exp_series` |
| 2026-10-09 | Do not build the IPM | Both profiles saturate, so the data put an upper bound on dispersion and no lower bound; `sigma_d` would come back at a boundary rather than as a measured rate |
| 2026-10-09 | The data constrain only the TOTAL stage spread | Accumulated spread varies 2.0–4.4x across rungs while the total varies 12.7–23.8%; `b_shape` falls monotonically with the rung in 14 of 14 units |
| 2026-10-09 | Arbitrate optimiser disagreements by the higher log-likelihood | Two-start `mat_exp_series` said TURNS OVER, eight-start convolution said SATURATES; forward maps agree to 1e-11 and the two-start path was 4.10 units short at `n_c` = 768 |
| 2026-10-09 | Keep the Erlang-window series; remove the convolution | Measured through `grad_log_prob` on the real data it is 5x slower at `n_c` = 96 and 2.8x at 768, and less accurate at high `n_c`; the 394x figure was against `mat_exp_series`, which the likelihood never calls |
| 2026-10-09 | Benchmark through the interface production uses | Two cost claims this session were wrong from measuring the wrong function: the 394x speedup, and ~540 h for an `n_c` = 768 fit that is really ~3 days |
| 2026-10-10 | Run the tuned retry of entry 40 AND a `b_shape` ladder at `n_c` = 384 (SLURM 30830) | Entry 40's reseed failed by general mixing failure, not one stuck chain, with 34% treedepth saturation; the 400 pin was chosen at `n_c` = 96 and `b_shape` falls as `n_c` rises, so the pin itself is in question at 384 |
