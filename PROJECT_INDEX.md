# plasmofit — malaria within-host dynamics — index

**Status as of 2026-10-15 (body last revised 2026-10-13; see the workstream table for later results):** The age-structure question is **answered**, and
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
`n_c` = 96 and may not carry. **First read 2026-10-13:** at `n_c` = 384 only
`b_shape` pinned at 15 and `b_shape` free converge; every pin from 50 up fails
like 400 did. Free `b_shape` is 12.2, pin 15 ties it on elpd, and
`cycle_length` is 41.34 h and 41.82 h respectively (means over 13 trials).
The tuned high pins (42, 49–53) were cancelled as projected to overrun.
**Runtime:** a benchmark found production fits run ~5x slower per gradient
than they need to, because chains share L3 caches and hyperthreads. The pinned
one-chain-per-L3-group driver refit entry 43 in 3.2 h against 38.5 h (12x),
but that refit fails the gate: a second region ~3 h lower for MMV048_PIB,
which the original run never visited, so the pin-15 result is less secure
than it looked.

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
| Cycle-length bias attribution | **Attributed** to `b_shape` (−0.92 h of +1.97 h); `fix_bshape` stands at n_pair = 1 | Simulation study read 2026-10-20: 9 of 48 cells readable, biases −0.60 to +1.19 h; fitting `n_c` 384 instead of 192 shortens `cycle_length` in 12 of 12 paired replicates (mean −0.55 h); the misspecified gen 384 / fit 192 cells are unreadable. Reported as it stands, no second stage (Lucas, 2026-10-21): `reports/simstudy-2026-10-21.html` |
| Production `n_c` ladder | **One rung converged** (30525, `n_c` = 192, 41.20 h, max R-hat 1.043) at 5.8% divergences | Entry 40 failed twice (29635; reseed 30829, 2026-10-10, general mixing failure, 34% treedepth saturation). `n_c` = 384 ladder partly read: only pin 15 and free converge (41.34, 41.82 h). 42 and 49–53 cancelled. Pinned driver 12x faster (30856) and adopted; its refit of pin 15 failed the gate via a second region for MMV048_PIB. **SLURM 30858**: pin 15 passes the gate in 1 of 5 runs (failures in four different trials; two regions ~3 h apart in both MMV048 trials); pooled over 5 runs the mean over trials spans only 41.29–41.41 h. **Profile (30869, 2026-10-14):** MMV048_PIB at pin 15 has two likelihood peaks, 38.5 and 43.0 h, within 0.65 ll units; MMV048_PartB's free profile is flat over 38.5–43.5 h. **Free seeds (30858, 2026-10-15):** 2 of 4 runs pass, giving 41.82 and 41.79 h; both failures have one chain at higher `lp__`. **Simulation study of `cycle_length` bias running** (SLURM 30875, 2026-10-15; `TODO.md`) |
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
- Shareable account of the cycle-length bias (2026-10-18):
  `reports/cycle-length-bias-2026-10-18.md` (markdown for Google Docs
  import) and `.html`, published privately at
  https://claude.ai/artifact/4z8jeubKob6e341xXayFcK (version 1)
- Simulation study report (2026-10-21, for Lucas, statistics explained):
  `reports/simstudy-2026-10-21.html`, published privately at
  https://claude.ai/artifact/7seQba92rjtufpQgHxWn3H
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
| 2026-10-13 | Do not read an elpd plateau at `n_c` = 384 from the default-settings ladder | Only pin 15 and free converged; every pin from 50 up failed the gate, so there is no ladder of converged rungs to plateau |
| 2026-10-13 | Cancel 30837 (49–53) | Projected 250–370 h against a 168 h walltime, and 49 duplicated converged entry 43 |
| 2026-10-13 | Treat fit runtime as a layout problem before a model problem | Same model ran 89–345 ms per gradient by node load; benchmark: sharing an L3 group halves speed, best layout ~16 ms against 83 ms in production |
| 2026-10-13 | Cancel entry 42 | Projected ~180 h against a 168 h walltime; a timed-out task writes nothing |
| 2026-10-13 | Run chains as separate pinned processes, one per L3 cache group | Thread count follows the CPUs a process may use, and only pinning sets it; benchmark puts this layout at ~16 ms per gradient against 83 ms |
| 2026-10-13 | Adopt the pinned driver for all fits | Refit entry 43 in 3.2 h against 38.5 h with `log_prob` unchanged by layout; posterior agreed in every trial apart from a region the original never visited |
| 2026-10-13 | Probe MMV048_PIB's second region with three more seeds each of pin 15 and free | One of eight pin-15 chains found a region ~3 h lower at equal `lp__`; at ~3-13 h per fit, seeds are now cheap enough to measure how often |
| 2026-10-15 | Run the simulation study of `cycle_length` bias: 2 truths × `n_c` 2x2 {192, 384} × true `b_shape` ×{0.5, 1, 3} × {estimated under the wide prior, fixed at truth} × 5 reps, 8 chains, simulations only | Measure the bias at production settings, never measured at `n_c` 384 or with `b_shape` fixed; Lucas chose the design, no pilot, and a 75% node allowance for this study only |
| 2026-10-18 | Cancel every remaining fixed-`b_shape` `n_c` 384 chain task (reps 3-5) | 23 of 24 such fits failed the gate in reps 1-2, 20 severely, even at the true value; they were the slowest tasks (median 18 h) and their bias could not become readable |
| 2026-10-18 | Move the simulation study's remaining fit `n_c` 384 chains (192 pending, 16 just started) to the pinned layout; always pin from now on | Packed single-threaded chains ran at ~260 ms per gradient against 28 ms pinned; most of a 13 h chain was layout (`findings.md`, "Where the simulation study's time goes"). Lucas, 2026-10-18 |
| 2026-10-21 | Report the simulation study as it stands; no 1000-draw second stage | 39 of 48 cells unreadable, but a refit would take ~3.7-7 days at the half-node budget, and refitting only the failures gives them a second try the passing fits do not get. Lucas |
| 2026-10-21 | Second simulation work: one truth (posterior means of `np_wide_total0_nc384_t`), generate at `n_c` 384 and 1024, always fit at 384, `b_shape` prior lognormal(2, 0.5), one fitting procedure (that real-data fit's) | Fewer fits, each sampled to converge; the ML profile puts the best rung at 1024 (2-unit interval {512, 1024}), so the plausible misspecification is a true `n_c` above the fit's. Lucas |
| 2026-10-21 | Run a single-unit study first (`_scripts/unitsim.R`); the full-model study (`_scripts/simstudy2.R`) is on hold and runs only on Lucas's explicit say-so | If the concern is structural bias, start with the smallest problem that has it: one trajectory, no hierarchy, cheap enough for many replicates. Lucas |
| 2026-10-21 | Single-unit study: the shortest-window unit (MMV048_PIB.1800) and the longest (Mefloquine.1800); `b_shape` both estimated and fixed at truth, the fixed arm's sampler settings to be chosen from a convergence probe | Units chosen on design alone; the fixed arm tests whether the full model's fixed-`b_shape` failures (1 of 25 passed) occur within one trajectory. Lucas |
| 2026-10-21 | Single-unit study: n = 60 replicates per cell, fixed in advance; benchmark two chains per L3 group rather than change the package's threading | At 5-15 min per chain cost no longer binds, so n is set for precision (se ~0.25 h if replicate sd ~2 h). Lucas |
| 2026-10-22 | Single-unit study: fixed `b_shape` arm at 2000 warmup + 2000 draws for every fit, estimated arm at production length; no retries; a cell readable at >= 45 of 60 passing; refit replicate 1 | The fixed arm passed 3 of 4 probe fits at production length, the failure marginal; one setting per arm, chosen before the study, avoids giving failures a second try. Lucas |
| 2026-10-22 | Benchmark two chains per L3 group against one before the study; Lucas decides the layout from the numbers | Single-unit fits keep 2 of 4 pinned cores busy; a package change to thread over time points was weighed and not chosen. Lucas |
| 2026-10-18 | Correct the `n_c` elpd figure: 96 is 68.6 (se 9.7) behind 192 and 71.8 (se 13.2) behind 384 | The 2026-10-06 row and `CLAUDE.md` gave 71.8 against 192; recomputed from the saved LOO objects (`_scripts/nc-ladder-loo.R`, `_data/nc-ladder-loo-2026-10-18.txt`). The conclusion is unchanged |
