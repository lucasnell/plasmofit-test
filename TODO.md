# plasmofit — malaria within-host dynamics — to do

**Key dates**

- No external deadline recorded. Add submission or meeting dates here.

<!-- Long-form reasoning for each item lives in claude/threads.md, numbered to
match where a thread exists. This file is the actionable layer. -->

## Age structure and `n_c` (threads 5, 15)

### In progress

- [ ] **SLURM 30875 — simulation study of `cycle_length` bias**, submitted
      2026-10-15. 48 cells × 5 replicates = 240 fits × 8 chains = 1920
      chain tasks, `--array=1-1920%96`, 2 CPUs and 4G each (192 CPUs: the
      75% allowance Lucas gave for this study only). Design and reading rule
      are in the `_scripts/simstudy.R` header: truths A (`pooled_cl`, n_c 96,
      45.012 h) and B (entry 54, n_c 384, per-trial truths); generate × fit
      `n_c` in {192, 384}; true `b_shape` × {0.5, 1, 3}; `b_shape` estimated
      under sd_log_b_shape 1.5 or fixed at the truth. **SLURM 30876** merges
      and runs `simstudy-read.R` once the array ends (`afterany`), writing
      `_data/simstudy-read-<date>.txt`. Estimated 1-3 days at submission, which was wrong: chains fitting at `n_c` 384 average 12.9-13.6 h, at 192 1.8-1.9 h; the 2026-10-18 replay of the remaining queue (`_scripts/simstudy-eta.R`, `_data/simstudy-eta-2026-10-18.txt`) puts the chain array's end at Tue 2026-10-20 ~18:00 (10-90%: 15:00-22:00), then ~20 min of merging. Tasks run the
      frozen copy `_data/simstudy/simstudy-run.R` (md5 526f217f…, =
      commit `4ab83b9`), so `_scripts/simstudy.R` may be edited; **do not
      reinstall the package until the array ends.** Partial read at any
      time: `SIMSTUDY_MODE=merge Rscript --vanilla _scripts/simstudy.R`, then
      the reader.
      **Partial read 2026-10-18** (96 of 240 fits, reps 1-2;
      `findings.md`, "Simulation study, partial read"): 24 of 96 pass; no
      cell readable yet. Fixed `b_shape` at fit `n_c` 384 fails in 23 of
      24 fits, 20 severely, even at the truth. **Lucas chose to cancel the
      still-pending fixed-`b_shape` `n_c` 384 chain tasks** (reps 3-5):
      207 cancelled 2026-10-18, listed in
      `_data/simstudy/fixed384-cancelled-2026-10-18.txt`; the 46 running
      ones were cancelled too, on Lucas's say-so the same day
      (`_data/simstudy/fixed384-cancelled-running-2026-10-18.txt`). Task 1728 escaped both (a parsing slip) and was cancelled separately; it is appended to the first list. Those 12 cells will hold 2 replicates
      plus any rep-3+ fit whose 8 chains had all finished before the cancellation; merge skips fits
      with missing chains, and 30876 still runs (`afterany` counts
      cancelled tasks as ended).
- [ ] **Open, 2026-10-18: the study's `n_c` 384 chains run ~7x slower per
      gradient than with few chains running** (260 against 38 ms; ~1.8x at 192),
      so most of their 13 h is layout, not model (`findings.md`, "Where the
      simulation study's time goes"; `_data/simstudy-cost-2026-10-18.txt`).
      Lucas's question, a low `n_c` plus a desynchronisation-rate parameter
      to cut fit time, was answered against the existing record
      (`claude/ipm-decision.md`): inside the chain it cannot reach the
      dispersion the data want (Erlang floor); outside it, the only forward
      map we have that separates mesh from dispersion is slower per
      gradient than the chain at 384. **Decided (Lucas): re-run the 192
      pending fit `n_c` 384 chain tasks under the pinned layout.** Built:
      `_scripts/simstudy-pinned.sh` (workers of 4 physical cores claim task
      IDs from `_data/simstudy/pinned-tasks.txt` via `mkdir` locks in
      `_data/simstudy/claims/`) and `_scripts/simstudy-relayout.sh` (drops
      30875's throttle to 48, cancels those tasks, submits a 64-CPU and a
      32-CPU pinned job that wait on enough running 30875 tasks to stay
      within 192 CPUs, and adds both to 30876's dependency; records to
      `_data/simstudy/relayout-<date>.txt`). **Not yet run**: the
      `scontrol`/`scancel`/`sbatch` calls were blocked by the permission
      classifier, so Lucas runs it. After it runs, check one worker's
      first-chain log line ("threads, % cpu after 10 min") to confirm 4
      threads.

### Next

- [ ] **Examine the higher-density chains in the two failing free runs**
      (seed 20261013 chain 2, +5.9 to +8.1 `lp__`; seed 20261015 chain 4,
      +3.8 to +5.4). Which parameters differ, and are the passing free runs
      (41.82, 41.79 h) missing a better region? Cheap: reads saved fits.
- [ ] **Extend the `cycle_length` profile down to 35 h** for the three
      profiles with a maximum at the 36.0 h grid edge (KAE609 fixed 15 and
      estimated, MMV048_PartB fixed 15). Question: is the short-period
      region with `b_shape` at its cap a real peak, or does it rise to the
      35 h bound? Cheap: under an hour per task on 1 CPU.
      `cl-profile.R` refuses to overwrite, so the extension needs new output
      names. Source: `_data/cl-profile-read-2026-10-14.txt`.
- [ ] **Do the priors or the hierarchy put MMV048_PIB in its upper
      region?** The likelihood has two peaks of near-equal height (38.5 and
      43.0 h, 0.65 ll units apart, the lower one marginally higher), but the
      pooled pin-15 MCMC puts the trial at 43.01 h. Test, for example, by
      refitting with a different `cl_prior_center` and checking which region
      the chains choose. Ties into the `cl_prior_center` item below.
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

### Done

- [x] **2026-10-15 — seed probe (SLURM 30858) complete.** Free `b_shape`
      passes the gate in 2 of 4 runs (pin 15: 1 of 5); the two passing runs
      give 41.82 and 41.79 h (mean over 13 trials) and agree within 0.17 h
      per trial. Both free failures are narrow (1-2 quantities) and have one
      chain at HIGHER `lp__` (+4 to +8), not a stuck chain; not yet examined.
      `_data/seed-probe-read-2026-10-15.txt`; `findings.md`, "Seed probe,
      complete".
- [x] **2026-10-14 — `cycle_length` profile (SLURM 30869) read.**
      MMV048_PIB at pin 15: TWO PEAKS, 38.5 and 43.0 h, heights within 0.65
      ll units (lower marginally higher), 5.4-unit dip. Free `b_shape`
      merges them into one peak at 40.0 h with `b_shape` at its 5000 cap.
      MMV048_PartB: free profile flat within 0.5 units over 38.5-43.5 h.
      Control KAE609 free: TWO PEAKS by the rule, via the 36.0 h grid edge.
      `_data/cl-profile-read-2026-10-14.txt`; `findings.md`, "The
      `cycle_length` profile".
- [x] **2026-10-14 — pooled the five pin-15 runs**: mean over trials
      41.29-41.41 h across runs, pooled 41.335; R-hat 1.039 over 20 chains,
      which is dilution, not convergence. `_data/pin15-pooled-2026-10-14.txt`.
- [x] **2026-10-13 — pinned driver adopted** as the standard way to run fits
      (`CLAUDE.md`, "How this project works").
- [x] **2026-10-13 — pinned validation (30856) read**: 12.0x faster
      (3.20 h against 38.52 h); per-trial `cycle_length` |z| < 3 and elpd
      within 2 se, but the gate fails through MMV048_PIB.
      `_data/pinned-validate-2026-10-13.txt`.
- [x] **2026-10-13 — pinned driver built**: `WOCKFIT_MODE` chain/merge in
      `wockner-fit.R`, `_scripts/wockner-fit-pinned.sh`. Entries 49-54 folded
      back into `wockner-fit.R`, copies deleted, threading comment corrected.
- [x] **2026-10-13 — cancelled entry 42** after 3 d 4 h (projected ~180 h
      against a 168 h walltime).
- [x] **2026-10-13 — thread-layout benchmark** (30843-30853): sharing an L3
      group halves speed; best layout ~16 ms per gradient against 83 ms for
      production. `_data/thread-bench-read-2026-10-13.txt`.
- [x] **2026-10-13 — cancelled 30837 (49-53)**: projected to overrun the
      7-day walltime; 49 duplicated converged 43.
- [x] **2026-10-13 — `n_c` = 384 ladder, first read** (43-48, 54): only pin
      15 and free (tuned) converge; free `b_shape` 12.2; pin 15 ties free
      (−0.49, se 4.06); `cycle_length` 41.34 and 41.82 h. Output
      `_data/nc384-ladder-read-2026-10-13.txt`.
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
