# Open threads

*Split out of `CLAUDE.md`, which is now an index. See `CLAUDE.md` for
orientation and the current state of play.*

Roughly in priority order.

1. **Finish the inoculum-anchor regression test and take the first anchored
   fit.** The test is one fit short of a verdict -- see `findings.md`,
   "Regression test for the anchor". Posterior means already agree within
   Monte Carlo error and the parameter sets match exactly; what is missing
   is the null run (SLURM
   `28892`) that posterior widths get judged against. Re-run
   `_scripts/wockner-anchor-regression.R` once it lands. **Do not** weaken
   the thresholds in that script to get a PASS: they were fixed before the
   null run existed, and the whole point of the null is that it is the only
   thing entitled to move them. Then fit the real data with
   `inoc_size = "inoc_size"` and read `delta_total0` against the predicted
   +0.811, and `R` against 9.96.
2. **Attribute the cycle-length bias.** Now the live scientific question,
   since the truth-construction explanation is dead (see `findings.md`,
   "Rebuilt from posterior draws"). Paired measurement puts ~0.25 h on the
   hierarchy and ~0.15 h on the cycle-length prior, against a total of
   +0.86 to +2.58 h.
   The remaining suspect is the nuisance priors dragging `cycle_length`
   along the `log10_total0`/`R` ridge; the maximum-likelihood/posterior gap
   is 1.3-2.1 h on the same data. **The direct test**: refit simulated data
   with `mean_log_b_shape`/`sd_log_b_shape` and
   `mean_log10_total0`/`sd_log10_total0` widened, and see how much of the
   bias moves with them. Both are already `archer_stan_data()` arguments, so
   this needs no package change. Note the inoculum anchor is a second,
   independent route at the same target: it replaces the `log10_total0`
   prior with a data-derived one, so if that prior is carrying the bias, the
   anchored fit should move `cycle_length`.
   Checked and closed:
   - ~~Prior pull on `cycle_length` (explanation 1).~~ Quadrupling the prior
     variance moves the estimate 0.15 h.
   - ~~The hierarchy.~~ `no_hier` (`pooled_cl`) still biases +0.93 and
     +1.97 h on replicates 2 and 3.
   - ~~Posterior mean versus mode on a bounded parameter.~~ Differ by
     0.044 h, posterior symmetric (skew -0.016), via
     `wockner-schedule-sim-mode.R`.
   - ~~Incoherent (mean-vector) truth.~~ Rebuilt from posterior draws;
     recovery is no better and the bias survives.
   Still open within this: **bound geometry**. `[35, 50]` leaves 5 h above a
   truth of 45 and 10 h below. Re-simulate with the bounds *moved*, e.g.
   `[30, 60]`, not merely widened -- `max_cl = 55` reintroduces the boundary
   mode.
   **Now running on real data**, SLURM 28931 tasks 10-13: the simulation
   result says `b_shape`'s prior owns ~0.50 h and `log10_total0`'s owns
   none, so the real-data question is whether the same dissociation holds
   where there is no truth to recover. Read with
   `_scripts/wockner-prior-panel.R`, section 4. Note what a real-data answer
   can and cannot be: a shift in `cycle_length` when the prior is widened
   bounds the prior's contribution from below, but it cannot show the
   remainder is biological.
   Until one of these lands, no cycle-length number should be reported as an
   estimate of anything biological.
3. ~~**Why are the simulated data less informative than the real data?**~~
   **Answered and closed**: they are not. Information about the oscillation
   is 0.97 of the real design's, comparing like with like (see
   `findings.md`, "Information content"). What replaces it: **is the
   estimator biased at this design, and
   is the real fit subject to the same bias?** If simulating from `b_shape`
   = 14.9 returns 8.90, the real fit's 18.9 is an estimate from the same
   biased estimator rather than evidence the real data pin it down. The test
   is the one already running for thread 2 -- if widening the nuisance priors
   removes the recovery bias, the estimator is prior-driven; if it does not,
   the bias is structural and every nuisance number in this project,
   including the real ones, needs re-reading. Separately, the real data's
   within-series time trends are steeper than the fitted trajectories (median
   range 2.56 against 1.95-2.21 log10 units), which is an unexplored
   lack-of-fit signal.
   **Also now running**, SLURM 28931 tasks 12-13: every `no_pool` vs
   `pooled_cl` conclusion was drawn under the `normal(1, 0.25)` prior on
   `log10_total0` that costs 104.7 elpd, so the comparison is being re-run
   with that prior corrected and with both nuisance priors corrected. Read
   with `_scripts/wockner-prior-panel.R`, section 5. The number to watch is
   the cycle-length pooling offset, -0.304 h as originally run.
4. **`hold_out` masking in `plasmofit`, then Design A.** The enabling change
   for cross-validation that does not rely on PSIS: a per-observation 0/1
   `hold_out` in `data`, with `transformed data` ordering each combo's kept
   observations first and storing a second length, so the model block's
   `reduce_sum` uses the fitting length while `generated quantities` keeps
   the full range. `log_lik` there (archer_fit.stan:492-515) is already
   computed independently of the model block, which is what makes this work.
   No change to `traj_combo_partial_sum` or the dedup; all-zero `hold_out`
   must reproduce current fits exactly, which is the regression test. Apply
   to all four Stan programs, add `archer_stan_data(hold_out = )`.
   Then **Design A**: mask the later portion of every series and compare
   held-out elpd across all four model variants (~4 fits). Each trial's
   initial conditions, error scale and `eta_cl[j]` stay informed, so
   `no_pool` can adapt per trial while `pooled_cl` cannot, and cycle-length
   error shows up as accumulated phase drift exactly in the held-out window.
5. **Design B, only if A is ambiguous.** True leave-one-trial-out K-fold,
   13 folds x 2-4 models. Note it is structurally near-rigged against the
   hierarchy: for a never-seen trial, `no_pool`'s point prediction collapses
   to the population mean, the same location `pooled_cl` gives, so it can
   only win on calibration. That likely explains why trial-level `loo` put
   `pooled_cl` marginally ahead.
6. **`cl_prior_center` decision.** Decide whether the default should move
   off 48 h. Demoted: the simulation showed the prior carries less of the
   error than thought (weight 0.113), so this mostly does not fix anything.
   Matters for reporting a cycle-length number; mostly cancels for model
   comparison.
7. `_scripts/test-archer-fit.R` (the driver script, not the package's
   `tests/testthat/test-archer-fit.R`) has not been run to completion with a
   full-length fit; it has only been smoke-tested with a short one, where
   recovery was good (23/23 parameters inside their 95% intervals, max |z|
   0.76). Worth revisiting in light of the schedule-bias finding: that smoke
   test used a denser design (8 observations per series, against 4-8 in
   Wockner).

### Lower priority

8. **More schedule-simulation replicates**, and more posterior-draw
   replicates specifically -- there are only two usable ones. The
   correlation question is limited by noise realizations, not by compute.
   Add seeds to `REP_SEEDS` in `wockner-schedule-sim.R` and widen the array,
   or submit more `SCHEDSIM_TRUTH_DRAW` values; ~2 h wall clock each.
9. ~~**Re-run the two non-converged replicates with a different fit seed.**~~
   **Done**, `SCHEDSIM_FIT_SEED=415926535`. Both converge on the new seed, so
   neither was a hard dataset -- it was the chain initialization.
   `wide-rep1`: R-hat 1.23 -> 1.0255, divergences 11.5% -> 2.9%, ESS 13 ->
   180, bias **+2.83 h**. `default-rep2-draw1500`: R-hat 1.107 -> 1.0083,
   divergences 7.5% -> 0.6%, ESS 26 -> 239, bias **+1.82 h**. The
   posterior-draw result is therefore n=3, not n=2: +1.68, +1.82, +0.86,
   mean +1.45. Worth noting for any future replicate that fails: try another
   fit seed before concluding anything about the dataset.
10. **Submit the cycle-length prior sensitivity runs.** `wockner-fit.R`
    entries 3-7, `--array=3-7`. Written but never submitted. Demoted: these
    were to discriminate explanation 1 from 2-3 for the sampling-density
    correlation. The simulation has since killed 1 and found against 2, so
    these would now be confirming on real data that the prior is not the
    story -- worth something, no longer decisive.
