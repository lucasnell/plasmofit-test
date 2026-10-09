# plasmofit — malaria within-host dynamics — to do

**Key dates**

- No external deadline recorded. Add submission or meeting dates here.

<!-- Long-form reasoning for each item lives in claude/threads.md, numbered to
match where a thread exists. This file is the actionable layer. -->

## Age structure and `n_c` (threads 5, 15)

### In progress

- [ ] **Entry 40 reseed — SLURM 30829**, submitted 2026-10-09, **~24.5 h**
      (29635_40 ran 24:34:15; walltime 3 days). `np_bs400_nc384`,
      `WOCKFIT_SEED=20261009`, `WOCKFIT_SUFFIX=-seed2`. **The only thing
      blocking a reportable cycle length.** Check `lp__` per chain before the
      R-hat gate, then compare against entry 39's 41.203 h: the deterministic
      tests predict the `n_c` effect PERSISTS AND GROWS with `b_shape` pinned
      (0.98 h at 400 against 0.72 h at 15, all twelve cells). If it converges
      and the effect is gone, the mechanistic account is wrong and should be
      revisited, not patched. If a chain sticks again at a similar lp gap the
      mode is real -- report multimodality and do not reseed a third time.

### Next

- [ ] **Decide the `n_c` default.** The ML profile says `cycle_length` is
      stable at **40.93–41.15 h once `n_c` >= 384** and the gains beyond 512
      are under 1.3 ll units over 1130 observations, so the choice is between
      384 and 512 on cost. 96 and 192 are out.
- [ ] **Decide whether to port the Poisson convolution into the package** —
      see the forward-map workstream below. Independent of everything here.

### Blocked

- [ ] Nothing.

### Done

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

## Forward-map speedup, opened 2026-10-08 (thread 15)

### Next

- [ ] **Decide whether to port the Poisson convolution into the package.** It
      is the *same* model as `mat_exp_series`, agreeing to 1e-12, and runs 38x
      faster at `n_c` = 384 and 14x at 192 (`_scripts/ipm-prototype.R`). That
      would make the `n_c` ladder cheap enough to extend, and a production
      `n_c` = 768 fit feasible where it is currently ~540 h. **Independent of
      the IPM decision** — it changes no model and invalidates no fit. Needs:
      a Stan implementation of the convolution (FFT is not available in Stan,
      so a direct lattice convolution, which is O(n_c * window) rather than
      O(n_c^3)), gradient checks, and the same `max_rel_diff` cross-check
      against `mat_exp_series` that the current path has.
