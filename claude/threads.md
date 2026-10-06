# Open threads

*Split out of `CLAUDE.md`, which is now an index. See `CLAUDE.md` for
orientation and the current state of play.*

Roughly in priority order.

1. ~~**Finish the inoculum-anchor regression test and take the first anchored
   fit.**~~ **CLOSED 2026-10-02, PASS.** The null run was never missing:
   SLURM 28892 died at the summary stage after writing its fit, and
   `wock-schedsim-fit-default-rep2-seed271828183.rds` had been on disk since
   2026-09-24. Posteriors agree better across the change than two runs of
   identical code differ (|z| median 0.73 against the null's 1.12; width
   inflation 1.01x; nothing beyond 3 MCSE in 208 entries). See
   `findings.md`, "Thread 1's anchor regression". The anchored fit itself was
   taken as `np_anchor` and is written up under "The `log10_total0` prior was
   misspecified". **Lesson worth keeping: a FAILED SLURM state does not mean
   no output** -- the fit is written before the summary.
   ~~The test is one fit short of a verdict -- see `findings.md`,
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
   **Attributed 2026-10-05**: `b_shape` is the nuisance that carries it.
   Pinning it at the simulated truth (`fix_bshape`) removes **−0.92 h** of
   +1.97 h, and merely widening its prior (`wide_bshape`, 6 converged pairs,
   every one negative) removes **−0.51 h**. `fix_sd` (+0.104) and
   `wide_total0` (+0.147) remove nothing, and `wide_nuis` — all five widened
   — gains no more than `wide_bshape` alone. **Of the nuisances only
   `b_shape` moves `cycle_length`.** With the likelihood's own +0.41 h and
   prior location's ~0.44 h, the budget is of the right order to close, but
   the arms overlap and must not be summed. See `findings.md`, "`b_shape` is
   the nuisance whose estimation carries the bias".
   **Open**: `fix_bshape` has only one converged pair of three; reps 2 and 3
   are refitting on a second seed. The conclusion rests on `wide_bshape`
   until they land.
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
   **2026-10-05: the last named component is eliminated and this thread is
   now IDEA-LIMITED.** The MLE-to-posterior gap -- per-trial MLE −0.195 h,
   pooled MLE +0.544 h, fitted posterior +1.56 to +1.71 h -- had two named
   untested components left. The MAP check could not test marginalisation
   (the joint surface has no usable mode: modes span ~2000 nats and 7-13 h,
   the best never found twice in ~200 starts). `fix_sd` tested the other and
   it is **+0.104 h**, the wrong sign and an order of magnitude short. There
   is nothing named left to try, and more fits will not help; what is needed
   is a different idea about where a correctly-specified posterior mean can
   sit 1.5 h from the truth.
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
4. ~~**`hold_out` masking in `plasmofit`, then Design A.**~~ **DONE
   2026-10-03, and it reversed the standing answer.** The masking shipped in
   the package (per-observation 0/1, kept observations ordered first within
   each trajectory combo, `generated quantities` scoring every observation
   whether fitted or not; all-zero reproduces an unmasked fit bit-for-bit,
   asserted in `tests/testthat/test-holdout.R`). Design A then masked the
   last third of every series and scored the held-out window.
   **The answer is mask-dependent and so the hierarchy question is NOT
   settled.** At the last third, collapsing `cycle_length` costs 3.27 log
   units (z −4.5); at the last quarter, −0.04 (z −0.1). Both fractions were
   fixed before any fit ran. Decomposing the third-mask fits shows 89% of the
   deficit sits on the 218 points the quarter mask *also* holds out, so it is
   not a horizon effect -- the difference is that the quarter-mask fits were
   TRAINED on the 88 deepest points and the third-mask fits were not. See
   `findings.md`, "Design A's verdict is mask-dependent".
   **This is what running both fractions was for.** A single mask would have
   reported a reversal that does not hold. And a masked fit scores a held-out WINDOW, not a held-out
   TRIAL -- so this answers whether per-trial cycle lengths predict better,
   not whether a never-seen trial would be predicted better, which is thread
   5's question and remains structurally near-rigged.
5. **`n_c` is a biological assumption, not a numerical setting, and 96 is
   wrong.** Raised, parked, and unparked 2026-10-05; **answered 2026-10-06**
   by SLURM 29633 (entries 37-38). See `findings.md`, "`n_c` is a biological
   assumption, not a numerical setting".
   `b_shape` falls monotonically along the ladder, 64.45 -> 61.57 -> 25.94,
   which was the predicted branch. `cycle_length` moves **-3.20 h**
   (44.25 -> 42.07 -> 41.04), larger than the whole +1.97 h simulated bias.
   **`n_c = 96` is 71.8 elpd worse than 384** (se 13.2); 192 vs 384 is -3.2
   (se 5.0), so the ladder **plateaus at 192**.
   **The effect is structural, not numerical**: `max_rel_diff` at `n_c = 96`
   is 8e-13, so the series solution is faithful to one part in 10^12 and
   `check_erlang_window()` was doing its job. No earlier fit is suspect
   arithmetic -- the 96 model is computed correctly and fits worse.
   **STILL OPEN, and it gates every reported cycle length.** 29633 ran with
   `b_shape` estimated, which is right for testing the confound but is not
   the production configuration (`b_shape` pinned at 400). Entries 39-40
   (`np_bs400_nc192`, `np_bs400_nc384`), script
   `_scripts/wockner-fit-nc-bs400.sh`, **written and not yet submitted**, ask
   whether the effect survives a pinned `b_shape`. Outcomes recorded in
   advance there and in `wockner-fit.R`.
   **IN FLIGHT, SLURM 29635** (entries 39-40, `np_bs400_nc192`,
   `np_bs400_nc384`): does the effect survive a **pinned** `b_shape`? This
   gates every reported cycle length.
   **IN FLIGHT, SLURM 29637**: the simulation **2x2**. Generating `n_c` and
   fitting `n_c` are now decoupled in `wockner-schedule-sim.R` via
   `arm$sim_n_c` and `arm$build$n_c`, so the four cells separate ESTIMATOR
   bias (the correctly specified diagonal) from MISSPECIFICATION bias (the
   off-diagonal). Arms `sim192_fit96`, `sim192_fit192`, `sim96_fit192`
   against the existing `default` (96/96). Read with
   `_scripts/nc-2x2-read.R`; see `findings.md`, "The `n_c` 2x2".
   **`sim96_fit192` is the cell that could overturn the real-data result**: if
   fitting ABOVE the true `n_c` drags `cycle_length` down, the -3.20 h is an
   artefact of over-large `n_c` rather than a correction.
   **Pairing, deliberately asymmetric**: `sim96_fit192` shares `default`'s
   simulated data (verified, y_sim mean log10 2.7590 both) and is in
   `PAIRED_ARMS`; the gen-192 arms share data with each other (2.6934), NOT
   with `default`, so they are excluded from it and read against the true
   `cycle_length` instead.
   `build_A` uses `lambda = n_c / cycle_length` over `n_c` sequential
   exponential compartments, so transit over one cycle is Erlang(`n_c`,
   lambda) and the stage distribution's sd after k cycles is
   `sqrt(k / n_c)` **cycles**. Over the Wockner window (k = 4.80):

   | `n_c` | sd after 1 cycle | sd after 4.8 cycles |
   |---|---|---|
   | 48 | 0.144 cyc | 0.316 cyc = 14.2 h |
   | **96 (default)** | **0.102 cyc** | **0.224 cyc = 10.1 h** |
   | 192 | 0.072 cyc | 0.158 cyc = 7.1 h |
   | 384 | 0.051 cyc | 0.112 cyc = 5.0 h |

   Cells: sd of the stage distribution, in cycles and in hours at a 45.012 h
   period, from the Erlang transit time. Deterministic, no data.

   So **the model already desynchronises deterministically**, and the rate is
   fixed entirely by `n_c` — which `check_erlang_window()` sizes for
   numerical accuracy ("roughly `n_c >= 6 * max(time) / min_cl`"), not for
   biology. The two uses are coupled: raising `n_c` for a longer window or a
   shorter `min_cl` also makes the parasites desynchronise more slowly.
   Nothing in this project has ever checked that rate against data.
   **There is no demographic stochasticity anywhere in the model** — the
   trajectory is deterministic given parameters and all noise is
   observational — so the stochastic half of the decay has no representation
   at all.
   **This cannot explain the simulated cycle-length bias**: the simulation
   generates from the same model with the same `n_c`, so the decay rate
   matches by construction. It is a real-data misspecification risk only.
   **Upgraded 2026-10-05, and the sentence above still stands.** `fix_bshape`
   showed that errors landing in `b_shape` propagate into `cycle_length`:
   pinning it at the truth removes −0.92 h of the +1.97 h bias, and of the
   five nuisances it is the only one that moves `cycle_length` at all.
   That does **not** make the decay rate the cause of the simulated bias —
   `n_c` is identical in the generating and fitted model, so the rate is
   right by construction, as above. What it changes is the stakes on real
   data: `b_shape` is the only free synchrony knob, so a wrong fixed decay
   rate would land there, and there is now a **measured channel** from
   `b_shape` to the headline estimate. This moves from "misspecification risk
   with no known consequence" to "misspecification risk with a quantified
   route to the number we report".
   **The connection worth testing if it is ever pursued**: initial synchrony
   and decay rate trade off against the observed late-time oscillation
   amplitude. `b_shape` 14.9 gives an initial sd of 0.090 cycles, 65 gives
   0.044, 400 gives 0.018 — all smaller than the 0.224 cycles the chain adds
   by the end of the window. So at late times the amplitude is dominated by
   the decay, which is fixed, and `b_shape` is the only free knob. A `b_shape`
   that runs to 65+ when its prior is relaxed may be absorbing a
   desynchronisation rate that is wrong, which would tie this to the
   `b_shape` non-identification rather than to the cycle-length bias.
6. **Design B, only if A is ambiguous.** True leave-one-trial-out K-fold,
   13 folds x 2-4 models. Note it is structurally near-rigged against the
   hierarchy: for a never-seen trial, `no_pool`'s point prediction collapses
   to the population mean, the same location `pooled_cl` gives, so it can
   only win on calibration. That likely explains why trial-level `loo` put
   `pooled_cl` marginally ahead.
7. ~~**Build-to-build drift in the weakly identified directions.**~~
   **CLOSED 2026-10-03, and the original claim was wrong.** The pre-change
   source was rebuilt into its own library and refitted at the same seed, so
   the comparison is identical source across two builds with nothing else
   mixed in. Drift is **mean |z| 0.79** (`b_shape` free) and **0.76**
   (pinned), against a seed-only null of 0.85 -- at or BELOW it. **There is
   no detectable rebuild drift**, and the caveat that `R` carries
   irreproducibility beyond its MCSE is **withdrawn**.
   The 2.5-3.5x that prompted this thread came from comparisons that crossed
   the `b_shape` code change, not from rebuilding. Of those, only
   `np_wide_total0` (1.23) exceeds the null; `np_bs400` (0.92) and the new
   data path (0.95) sit inside it. That one figure is unexplained, and with
   the null itself measured once there is nothing to say how far a single
   comparison should scatter. **If it matters, the fix is several null runs,
   not more rebuilds** -- `WOCKFIT_SEED` makes them cheap.
   Beware: the first attempt at this read mean |z| 25.5 because the two fits
   had different group level ORDERS, from an edit made while the job was
   running. See `gotchas.md`.
8. **`cl_prior_center` decision.** Decide whether the default should move
   off 48 h. Demoted: the simulation showed the prior carries less of the
   error than thought (weight 0.113), so this mostly does not fix anything.
   Matters for reporting a cycle-length number; mostly cancels for model
   comparison.
9. `_scripts/test-archer-fit.R` (the driver script, not the package's
   `tests/testthat/test-archer-fit.R`) has not been run to completion with a
   full-length fit; it has only been smoke-tested with a short one, where
   recovery was good (23/23 parameters inside their 95% intervals, max |z|
   0.76). Worth revisiting in light of the schedule-bias finding: that smoke
   test used a denser design (8 observations per series, against 4-8 in
   Wockner).

### Lower priority

10. **More schedule-simulation replicates**, and more posterior-draw
   replicates specifically -- there are only two usable ones. The
   correlation question is limited by noise realizations, not by compute.
   Add seeds to `REP_SEEDS` in `wockner-schedule-sim.R` and widen the array,
   or submit more `SCHEDSIM_TRUTH_DRAW` values; ~2 h wall clock each.
11. ~~**Re-run the two non-converged replicates with a different fit seed.**~~
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
12. **Submit the cycle-length prior sensitivity runs.** `wockner-fit.R`
    entries 3-7, `--array=3-7`. Written but never submitted. Demoted: these
    were to discriminate explanation 1 from 2-3 for the sampling-density
    correlation. The simulation has since killed 1 and found against 2, so
    these would now be confirming on real data that the prior is not the
    story -- worth something, no longer decisive.
13. **SBC (simulation-based calibration), parked 2026-10-05 with a condition
    for unparking.** Sized but never run: ~38 core-hours for N = 100 with
    short chains (L = 100 draws), ~1.5 h wall clock at 25 concurrent.
    **Why it is parked.** It was the pre-registered branch for `fix_bshape`
    coming back near zero, which would have meant no single nuisance carried
    the cycle-length bias and the cost was joint -- a diffuse result that
    needs a calibration instrument rather than another arm. `fix_bshape` came
    back at **-0.92 h**, so the bias is attributed to `b_shape` and the arms
    are still working. Running SBC now would answer a question that is no
    longer the blocking one.
    **What SBC would and would not tell us.** It checks that the posterior is
    calibrated under the model's own prior -- that is, over datasets drawn
    from the prior predictive. It therefore speaks to whether the *sampler
    and model code* are correct, not to whether the cycle-length bias is real
    on Wockner-like designs. The schedule simulation already answers the
    latter directly and at the design we care about, which is why it has been
    the instrument throughout.
    **Unpark it if** any of these happen:
    - the `fix_bshape` reseeds come back and the attribution collapses toward
      zero, putting us back in the joint-cost branch;
    - a cycle-length number is going into a manuscript and we need a
      calibration statement to accompany it;
    - thread 5's `n_c` work changes the model (not just its settings), since
      a structural change reopens whether the code is right.
    **Do not** run it merely because it is cheap and sitting there. It is a
    correctness check on code that has not changed, and the open questions
    are about identification, not implementation.
14. **Initial density should vary between individuals, and the model has no
    term for it.** Raised 2026-10-06, **not to be addressed now**. Sizing in
    `_scripts/total0-bloodvol-scale.R`.
    **The mechanism.** `para` is a concentration (iRBC/mL) but `inoc_size` is
    a COUNT -- 1800, 2300 or 2800, identical for everyone in a cohort. The
    conversion is a division by blood volume, which differs between people.
    Two people given the same number of parasites start at different
    *densities*, and density is what is observed.
    **What the model does now**, verified rather than assumed:
    - `archer_stan_data(blood_volume_ml = 5000)` -- one 5 L adult for
      everyone.
    - `log10_total0` is `vector[n_grp_init]`, so **14 values over 177
      series**. There is no per-individual initial-density parameter at all.
    - `sigma_total0` is the between-GROUP sd around the anchor, not
      between-individual, so it does not cover this either.
    **The package's own claim is half right.** The docs say an error in
    `blood_volume_ml` "is absorbed by the offset, which it is not separately
    identifiable from". True for the MEAN -- a constant error is exactly what
    `delta_total0` absorbs. False for the VARIANCE: one offset cannot absorb a
    spread. The claim should be qualified when this is written up.
    **But the spread is small.** Blood volume in screened healthy adults has
    CV roughly 12-20% (Nadler; volume scales sublinearly with weight, and the
    cohort is weight-screened), so sd on the log10 scale is CV/ln(10):

    | CV | sd (log10) | vs per-series SE | share of obs variance |
    |---|---|---|---|
    | 12% | 0.052 | 0.24x | 1.0% |
    | 15% | 0.065 | 0.30x | 1.5% |
    | 20% | 0.087 | 0.40x | 2.7% |

    Cells: implied between-individual sd of `log10_total0`; the per-series
    standard error is mean `sd_iRBC` 0.533 over a median 6 observations,
    = 0.217 log10.
    **Verdict: the mechanism is real and correctly identified, and as a free
    random effect it is not estimable from these data.** It is a quarter to
    two fifths of the noise on a single series' mean level. Adding a per-series
    random effect on `log10_total0` would shuffle variance between it and
    `sd_iRBC` without changing fit, and would add 177 weakly informed
    parameters. It does **not** threaten the anchor conclusion either: group
    means average over 13-24 individuals, so the blood-volume component shrinks
    to about 0.065/sqrt(17) = 0.016 log10 at group level.
    **The fix worth doing, if the covariate can be got.** Do not estimate it --
    COMPUTE it. The source volunteer-infection trials record weight (and often
    height and sex); Nadler's equation turns those into a per-individual blood
    volume, which makes `log10(inoc_size / blood_volume_i)` a **known**
    per-series anchor. That converts an unidentifiable variance component into
    a known offset at the cost of **zero** free parameters, keeping the one
    shared `delta_total0`. It needs the anchor and `log10_total0` to become
    per-series, which is a package change.
    **Blocker**: `wockner-cleaned.csv` has only id, trial, cohort, inoc_size,
    subject, time, para. No weight, height or sex. This cannot start until
    subject-level covariates are recovered from the source trials.
    **What would change the verdict**: a cohort with a wider weight range; a
    design with more observations per series, which shrinks the per-series SE
    the effect has to clear; or evidence that the low end is censored, since a
    persistent offset has the most leverage where the response is floored --
    see the note below.
    **Related and separate**: `para` has no zeros and a minimum of exactly
    1.0 across all 1130 observations, with 8.1% below 10 iRBC/mL. That looks
    like a floor at the quantitation limit, and the likelihood treats
    `log10(y + 1)` as normal with no censoring term. Worth a thread of its own
    eventually; it is where a small persistent per-individual offset would
    matter most.
