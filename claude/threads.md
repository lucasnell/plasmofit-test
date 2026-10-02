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
   **Budget updated 2026-10-01 at n=4**: the `b_shape` prior's share is
   **0.43 h**, not 0.50 h, and its spread across replicates is four times
   wider than three replicates suggested (-0.228 to -0.531). The old
   precision was luck. Separately, fixing `b_shape` rather than estimating
   it gains **+28 to +36 elpd** on real data, and `cycle_length` then depends
   on where it is pinned -- 44.0 h at 50 down to 43.6 h at 250 (see
   `findings.md`, "The `b_shape` ladder"). **Any cycle-length number is now
   conditional on that choice and must be reported with the ladder.**
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
   ~~Still open within this: **bound geometry**.~~ **Closed 2026-10-02 and it
   is not the answer.** Both `[37.5, 52.5]` (same width, centred on the
   truth, isolating asymmetry) and `[30, 60]` (moved and widened) change the
   paired cycle-length bias by **+0.059 h** and **+0.03 h** -- against the
   ~-1.5 h that would be needed, with the sign wrong. Three of six replicates
   failed the gate, so this is n=2 and n=1: adequate for "nowhere near the
   size needed", not for an estimate. **This thread is now out of cheap
   suspects**; what is left is structural, and thread 4's `hold_out` masking
   is the next real instrument.
   ~~**Now running on real data**, SLURM 28931 tasks 10-13.~~ **Done, and
   the dissociation reproduces** -- see `findings.md`, "the eight-fit prior
   panel". Widening `b_shape`'s prior moves real-data `cycle_length`
   -0.47 h in isolation (-0.501 h in simulation) and -1.29 h on top of a
   corrected `log10_total0` prior, leaving `log10_total0` and `R` untouched,
   and gains +24.8 elpd (se 2.23). The two nuisance priors turn out to be
   **super-additive on `cycle_length`** -- jointly -1.1 h against -0.27 h
   from the one-at-a-time contrasts -- so ~0.50 h is a floor for what the
   pair does together, not the pair's total. This is a **shift, not a
   bias**: it bounds the priors' share of the reported 45.3 h from below and
   cannot show the remainder is biological. The ~1.07 h unexplained in the
   simulation budget is untouched by it.
   **Answered 2026-09-25, and it does not shrink the remainder.** Simulating
   from a truth whose `b_shape` is 65 and fitting with both nuisance priors
   widened still leaves **+1.83 h** of cycle-length bias; the `pooled_cl`
   truth with the same widened priors leaves +1.63 h against `default`'s
   +1.97 h. Widening the nuisance priors buys ~0.5 h of ~2 h. **A shift in
   `cycle_length` when a prior is widened is not bias removal.** ~1.5 h is
   unaccounted for and **bound geometry is now the one cheap untested
   suspect left** -- re-run with `[35, 50]` *moved*, e.g. `[30, 60]`, not
   merely widened.
   SLURM 28940 added replicates 4-5 of every budget arm, but `default-rep4`
   and `default-rep5` both failed the convergence gate, so the **paired**
   budget is still n=3 and every number in it is unchanged. Refit those two
   on a new `SCHEDSIM_FIT_SEED` before reading anything more into the
   decomposition.
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
   **Answered 2026-09-25, both halves.** The estimator is biased at this
   design and the real fit is subject to it: under identical widened priors a
   true `b_shape` of 14.9 returns 19.5-27.9 and a true 65.0 returns
   32.2-37.4, so a factor of 4.36 in the truth becomes 1.53 in the estimate.
   `b_shape` is unidentified over at least 15-65 and no point estimate of it
   is defensible. Note the direction: the real data report 65.4, *above* what
   the design returns when the truth really is 65, so the attenuation runs
   the opposite way from "the real 65 is an artefact". Do not invert the
   calibration to get a number.
   ~~The lack-of-fit signal in the within-series time trends.~~ **Explained**:
   it was the misspecified `log10_total0` prior flattening the trajectories.
   Correcting it takes the median within-series range from 1.66 to 2.29
   against an observed 2.563, and the posterior predictive p-value from 0 to
   0.155. Widening `b_shape` does nothing for it. A residual flatness remains
   (`ppp` 0.155-0.17) but it is no longer a clear misfit.
   ~~**Also now running**, SLURM 28931 tasks 12-13.~~ **Done**: the
   `no_pool` vs `pooled_cl` comparison was re-run with `log10_total0`
   corrected and with both nuisance priors corrected. `pooled_cl` fails to
   beat `no_pool` at all three settings (-1.38, -1.67, -1.68 elpd, se ~1),
   so **the hierarchy conclusion survives a defensible prior**. The pooling
   offset does not survive as a number: -0.304, -0.110, -0.509 h across the
   three settings, a factor of 4.6 and non-monotone. Quote **-0.509 h**, the
   best-fitting setting.
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
   **Done**, `SCHEDSIM_FIT_SEED=415926535`. **Qualified 2026-10-01**: a new
   seed is not always enough. Of three refits at
   `SCHEDSIM_FIT_SEED=1618033989`, two converged and `default-rep4` did not
   -- max R-hat 1.06 against a 1.05 gate, though divergences fell 478 -> 39
   and ESS rose 13 -> 117, so the seed did most of the work and the gate is
   close. Try a third seed or a longer warmup there; do not assume one reseed
   settles it. Both converge on the new seed, so
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
