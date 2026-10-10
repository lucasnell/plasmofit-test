# plasmofit — malaria within-host dynamics — to do

**Key dates**

- No external deadline recorded. Add submission or meeting dates here.

<!-- Long-form reasoning for each item lives in claude/threads.md, numbered to
match where a thread exists. This file is the actionable layer. -->

## Age structure and `n_c` (threads 5, 15)

### In progress

- [ ] **SLURM 30830**, submitted 2026-10-10, `_scripts/wockner-fit-nc384-ladder.sh`.
      Entry 42 is the tuned retry of the 400 rung at `n_c` = 384 (51–102 h
      expected, 7-day walltime); entries 43–48 are a `b_shape` ladder at
      `n_c` = 384 (15, 50, 100, 250, 600 pinned, and free; default sampler
      settings; ~20 h each). Chosen by Lucas after entry 40 failed twice.
      Prediction in `wockner-fit.R` above the configs. **Do not edit
      `wockner-fit.R` while it runs.** Read `lp__` per chain, then the gate,
      then paired PSIS-LOO across converged rungs; report `cycle_length` as its
      spread across converged rungs.
- [ ] **SLURM 30837 and 30842**, submitted 2026-10-10: the insurance runs at
      entry 42's tuned settings. 30837 is entries 49–53 (the pinned rungs,
      `_scripts/wockner-fit-nc384t-ladder.sh`); 30842 is entry 54 (the free
      rung, `_scripts/wockner-fit-nc384t48.sh`). **Cancel both if 43–48
      converge.** They run copies (`wockner-fit-nc384t.R`,
      `wockner-fit-nc384t48.R`); merge 49–54 into `wockner-fit.R` and delete
      the copies once nothing is running.

### Next

- [ ] **Decide the `n_c` default.** The ML profile says `cycle_length` is
      stable at **40.93–41.15 h once `n_c` >= 384** and the gains beyond 512
      are under 1.3 ll units over 1130 observations, so the choice is between
      384 and 512 on cost. 96 and 192 are out.
- [x] **Nothing to install.** `conv_series()` was added, measured, and removed
      again (`bba580f`); the package is byte-identical to `8dde0c1`, which is
      what is already in the live library. No reinstall is needed and the
      scratch library `/home2/lan68/plasmofit/.Rlib-dev` has been deleted.
- [x] **Decided: do NOT use `conv_series()` for fitting.** It is 5x slower
      than the Erlang-window series at `n_c` = 96 and 2.8x at 768 (the
      crossover is near `n_c` ~ 14,000), and it is also LESS accurate at high
      `n_c` — 5e-10 at 768 against `mat_exp_series`, which is itself stable to
      3e-16 across step sizes. Keep `use_conv = 0L`. The switch stays as a
      cross-check at `n_c` <= 384, where the two agree to 1e-13 on `log_prob`
      and its gradient.
- [ ] Nothing.

### Done

- [x] **2026-10-10 — entry 40 reseed (SLURM 30829) read: FAIL**, max R-hat
      7.744. Reader `_scripts/nc-bs400-ladder-read.R`, output
      `_data/nc-bs400-ladder-read-2026-10-10.txt`.
- [x] **2026-10-09 — the age-structure question is answered.** The data
      constrain only the **total** stage spread at the end of the window, not
      how it splits between initial synchrony and accumulated
      desynchronisation. Accumulated spread varies 2.0x (chain) and 4.4x
      (gamma IPM) across rungs while the total varies 12.7% and 23.8%, and
      `b_shape` falls monotonically with the rung in **14 of 14 units**
      (Spearman −1.00). That explains the decay law's unlearnability, the
      profile's saturation, and `b_shape`'s non-identification at once.
- [x] **2026-10-09 — do not build the IPM.** Both profiles SATURATE: upper
      bound on dispersion, no lower bound, so `sigma_d` would come back at a
      boundary rather than as a rate. Non-uniform stage rates do not help
      either — they only ADD dispersion above the Erlang floor, and the data
      want less. The cost argument had already gone when the exact chain
      turned out to be a convolution.
- [x] **2026-10-08/09 — the production ladder converged on a reseed.** SLURM
      30525, `n_c` = 192, `b_shape` 400: max R-hat 1.0432, `lp__` spread 8.3
      across chains against 29635's ~96. **41.203 h**, mean over 13 trials,
      at 5.8% divergences.
- [x] **2026-10-09 — the two-start/eight-start conflict is resolved.**
      `mat_exp_series` said TURNS OVER at 512, the eight-start convolution
      said SATURATES; the forward maps agree to 1e-11, every rung-wise
      difference was >= 0, and the two-start path was 4.10 units short at
      `n_c` = 768. The higher log-likelihood wins.
- [x] Decay-law test (29684): cannot distinguish √ from linear — and now
      explained, since the total spread is nearly fixed across rungs.
- [x] `n_c` = 96 is wrong without any prior: +72.4 ll units in 14 of 14 units.
- [x] Simulation 2x2 (29637): misspecifying `n_c` costs ~1–1.45 h either way.
- [x] The observable period is not the `cycle_length` parameter.
- [x] IPM prototype built and validated: the exact chain is a convolution,
      1e-12 against `mat_exp_series` and 38x faster per trajectory.

## Cycle-length bias (thread 2)

### Next

- [ ] **Decide what cycle-length number to report, and with what
      qualifications.** It is conditional on the `b_shape` pin *and* now on
      `n_c`. Bias-corrected estimates sit in the low 40s under either
      hypothesis about the true `n_c`.

### Blocked

- [ ] **All reported cycle lengths** — **blocked on:** the production `n_c`
      ladder, since 2026-10-06.

### Done

- [x] Attributed the bias to `b_shape` (−0.92 h pinned at truth, −0.51 h
      widened; `sd_iRBC` and `log10_total0` carry none).
- [x] Localised the rest: +0.41 h sits in the likelihood with all nuisances
      known, ~0.44 h in prior location.
- [x] Refuted phase-volume, phase-period compensation, bound geometry,
      median-vs-mean, `sd_iRBC`, and MAP as explanations.

## Model structure and misspecification (threads 14, plus censoring)

### Blocked

- [ ] **Per-individual initial density** — **blocked on:** subject-level
      weights from the source trials, since 2026-10-06.
      `wockner-cleaned.csv` has none. The mechanism is real (`para` is a
      concentration, `inoc_size` a count, blood volume differs) but as a free
      random effect it is not estimable: 0.05–0.09 log10 against a per-series
      standard error of 0.217. The version worth doing *computes* it from
      weight via Nadler's equation, at zero free parameters.

### Next

- [ ] **Check the low-end censoring.** `para` has no zeros and a minimum of
      exactly 1.0 across all 1130 observations, with 8.1% below 10 iRBC/mL.
      That looks like a quantitation floor the likelihood does not model.

## Validation and infrastructure (threads 7, 9, 10, 12, 13)

### Next

- [ ] More schedule-simulation replicates — the cheapest way to tighten every
      paired number, including `fix_bshape`'s n_pair = 1.
- [ ] Run `_scripts/test-archer-fit.R` to completion with a full-length fit;
      it has only ever been smoke-tested.

### Blocked

- [ ] **SBC** — **parked** 2026-10-05, not blocked. Unpark only if the
      `fix_bshape` attribution collapses, if a cycle-length number goes into a
      manuscript, or if the model changes structurally. Do not run it merely
      because it is cheap.

### Done

- [x] Build-to-build drift: no detectable drift; the earlier caveat withdrawn.
- [x] Reseeded the two non-converged schedule-simulation replicates.
- [x] Package committed and pushed (`8dde0c1`); generated Rcpp bindings had
      been left at `ocombo_len` while the Stan sources used `ocombo_len_fit`.
- [x] Deleted `.prechange/` scratch (494 MB) after verifying its worktree HEAD
      was an ancestor of `main`.

## Hierarchy and model comparison (threads 6, 8)

### Next

- [ ] **`cl_prior_center`** — decide whether the default moves off 48 h. Now
      entangled with `n_c`: the data pull downward at higher `n_c`, so a prior
      centred at 48 is further off than it looked.

### Done

- [x] Design A: mask-dependent, so the hierarchy claim is withdrawn.
- [x] Horizon ladder: flat once the non-converged rung is reseeded.
- [x] Design B: argued against as structurally near-rigged; not run.

## Forward map, opened and closed 2026-10-08/09 (thread 15)

### Done

- [x] **The Poisson convolution was ported into the package and removed
      again.** Added as `conv_series()` (`815f3a0`), removed by revert
      (`bba580f`); the package is byte-identical to `8dde0c1`. **Do not
      re-open this without new evidence.** It was proposed on a 394x speedup
      over `mat_exp_series`, which the likelihood never calls — it calls
      `ew_poly_eval`. Through `grad_log_prob` on the real data the convolution
      is 5x SLOWER at `n_c` = 96 and 2.8x at 768, crossover near `n_c` ~
      14,000, and LESS accurate at high `n_c` (5.1e-10 at 768 against a
      `mat_exp_series` stable to 3e-16 across step sizes). Two earlier notes
      here were wrong and are corrected: Stan **does** have an FFT (2.39.0),
      and an `n_c` = 768 production fit is ~3 days, not ~540 h.
- [x] **What it did establish, and which survives the removal.** Two
      independently written forward maps agree to 6.4e-14 on `log_prob` and
      1.7e-13 on its gradient on the real data, and `mat_exp_series` is stable
      to 3e-16 across step sizes. Neither was known before. Scripts in
      `_scripts/` need a package built from `05c9c1f`; outputs are tracked.
