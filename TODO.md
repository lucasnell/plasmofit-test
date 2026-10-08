# plasmofit — malaria within-host dynamics — to do

**Key dates**

- No external deadline recorded. Add submission or meeting dates here.

<!-- Long-form reasoning for each item lives in claude/threads.md, numbered to
match where a thread exists. This file is the actionable layer. -->

## Age structure and `n_c` (threads 5, 15)

### In progress

- [ ] **Resubmit the production `n_c` ladder.** This is now the only thing
      standing between the project and a reportable cycle length. 29635 failed
      (max R-hat 6.13 and 8.34, one chain stuck 96 and 70 lp units below the
      others). `_scripts/wockner-fit-nc-bs400-retry.sh` is written but **the
      costing says reseed-only first**: the sampler already saturated
      `max_treedepth` in 46–61% of transitions, so `adapt_delta` 0.95 plus
      `max_treedepth` 12 costs 2.5–5× (entry 41 ~21–42 h, entry 42 ~60–120 h,
      which exceeds the script's 3-day walltime). Reseeding alone targets the
      observed failure and costs ~9 h at `n_c` = 192.

### Next

- [ ] **Decide the `n_c` default** once the above lands. 96 is wrong — by 71.8
      elpd under the Bayesian comparison and by +72.4 log-likelihood units in
      14 of 14 units under a priors-free maximum-likelihood profile. 192 and
      384 are tied on elpd with `b_shape` free; the ML profile puts 384 ahead
      of 192 by +10.4 summed, in 10 of 14 units.
- [ ] **Give the desynchronisation rate its own parameter.** The real defect:
      `n_c` is both the numerical mesh and the biological rate, and raising it
      to fix the numerics removes a biological process. **The decay-law test
      did not determine the form** — see below — so the form is now a biology
      judgement, not something these data settle.

### Blocked

- [ ] **IPM / transport-with-dispersion rewrite** — **no longer blocked on
      evidence, because the evidence came back indecisive.** The decay-law test
      could not distinguish √ from linear: at the predicted sigma = 0.033 a
      linear component buys at most +0.396 and typically +0.0000
      log-likelihood units over 1130 observations. So an IPM cannot be
      justified by "√ is established", and a 1-D Gaussian-kernel IPM would
      reproduce √ by construction. **Deciding this needs a judgement call from
      Lucas**, weighed against the cost: discarding the validated Erlang-window
      series and its `matrix_exp` cross-check, and making every existing fit
      incomparable.

### Done

- [x] **Decay-law test (29684), 2026-10-08: cannot tell, and says why.** All 14
      tasks COMPLETED, 9–21 h each. Total `d_ll` +1.53, 8–6 on sign, +1.449 of
      it from one unit. B nests A, so the six negatives are optimiser failure,
      not evidence; the worst, −0.521, is the noise floor and only 1 of 14
      units clears it. Thread 15 and `claude/findings.md` have the numbers.
- [x] **`n_c` = 96 is wrong without any prior.** Model A alone is a
      priors-free likelihood profile: 192 beats 96 in 14 of 14 units, +72.4
      summed, against the Bayesian 71.8 elpd on the same 1130 observations.
- [x] Real-data `n_c` ladder with `b_shape` free (29633): `cycle_length` moves
      −3.20 h, `b_shape` falls monotonically, 96 loses 71.8 elpd.
- [x] Ruled out a numerical cause: `max_rel_diff` is 8e−13 at `n_c` = 96 and
      *grows* with `n_c`.
- [x] Simulation 2×2 (29637): misspecifying `n_c` costs ~1–1.45 h in either
      direction; correcting both real rungs favours true `n_c` = 192.
- [x] Explained the mechanism: the observable period is not the
      `cycle_length` parameter, because a one-sided sequestration window meets
      a broadening skewed age distribution.

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
