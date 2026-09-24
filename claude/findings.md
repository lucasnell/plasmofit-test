# Findings: the hierarchy, the schedule-bias simulation, and the priors

*Split out of `CLAUDE.md`, which is now an index. See `CLAUDE.md` for
orientation and the current state of play.*

Every numeric table below states what is in its cells; see
`conventions.md`. Roughly chronological, so a later section can overturn
an earlier one -- the state of play in `CLAUDE.md` says which conclusions
are current.

## Cycle-length hierarchy: no_pool vs pooled_cl

Run via `wockner-fit.R`'s `no_pool`/`pooled_cl` configs plus
`wockner-fit-kfold.R`'s leave-one-trial-out refits (26 folds), both with the
settled priors and default `adapt_delta`/`max_treedepth` for both models (see
git history for the reasoning -- a small simulated-data sampling issue with
`pooled_cl` from `plasmofit`'s `chat4/CLAUDE.md`, with an unresolved
discrepancy about that simulation's true cycle_length, pooled vs per-series).
That issue **did not reproduce** at Wockner scale: both configs sampled
cleanly (~1% divergences, max R-hat 1.01, chains agreeing on `lp__` within
~0.5 units, no sign of the mode-splitting seen with `max_cl = 55`).

**Result: leans toward the hierarchy not earning its keep, but the sharpest
comparison is uninterpretable.**

- **Observation-level `loo`** (the weak comparison -- these are short series,
  a held-out point is pinned down by its neighbours regardless of model):
  `no_pool` elpd -1053, `pooled_cl` elpd -1055, `elpd_diff` -1.4 (se 1.0).
  Not distinguishable; matches `archer_log_lik()`'s own docs on what to
  expect from this grouping.
- **Trial-level `loo`** (the comparison that's actually supposed to answer
  this, per `archer_log_lik()`'s docs) **is broken**: Pareto k > 0.7 for all
  13 trial-units, in both models. Confirms the earlier "Leave-one-trial-out
  strains `loo`" gotcha applies here too, not just to a single model's
  self-assessment -- a trustworthy answer needs the actual K-fold refitting
  this section's robustness check stops short of (see below).
- **`div_by_sigma`**: `no_pool` still shows a mild funnel (divergences
  concentrate at lower `sigma_logit_cl`, 0.40 vs 0.50 for non-divergent
  draws -- down from the historical ~11x gap but not gone at ~1.25x).
  `pooled_cl` shows nothing analogous for `sigma_logit_R` (0.206 vs 0.213),
  as expected since it has no `sigma_logit_cl` to have a funnel in.
- **`wockner-fit-kfold.R`'s robustness check** (13 leave-one-trial-out
  refits per model -- **not** full predictive K-fold: scoring a held-out
  trial under `no_pool`'s posterior would need a way to predict a brand-new
  hierarchy group, which `archer_fit.stan` has no machinery for; this instead
  checks whether the population-level estimates are stable, see the script's
  header) backs up the full-data read rather than overturning it:
  - `sigma_logit_cl` stays close to its prior scale in every fold (mean 0.49
    across folds, range 0.34-0.61) -- the "data can't rule out zero
    between-trial variation" finding isn't a one-off of the full dataset.
  - The population cycle-length estimate is stable across which trial is
    held out: `no_pool` 44.9-45.9 h (sd 0.33), `pooled_cl` 44.6-45.3 h
    (sd 0.25). No trial is an outlier that's driving the result.
  - `no_pool`'s population estimate came out consistently ~0.3-0.7 h *higher*
    than `pooled_cl`'s in every one of the 13 folds. Since decomposed --
    see "The no_pool / pooled_cl cycle-length offset" below. Roughly half of
    it was an artifact of which summary got reported.

**Bottom line:** every comparison that actually ran points the same
direction (no detectable benefit from the hierarchy), but the one designed to
settle it decisively (trial-level `loo`) failed its own diagnostics in both
models. Proper K-fold (Stan-side masking of a held-out trial's likelihood, so
`log_lik` scores it honestly against the rest of the hierarchy -- deferred
when this script was written, see `wockner-fit-kfold.R`'s header) is still
the way to actually close this out.

## The no_pool / pooled_cl cycle-length offset

Decomposed by `wockner-pooling-offset.R` (post-hoc on the saved full-data
fits, no refitting). Of the 0.526 h gap as originally reported:

Cells: `size` is hours of the `no_pool` minus `pooled_cl` difference in
population `cycle_length` attributable to that component, from one pair of
saved full-data fits (no replicates).

| component | size | what it is |
|---|---|---|
| summary statistic | 0.222 h | artifact, not a model difference |
| differential prior pull | ~0.118 h | real |
| trial- vs observation-weighting | ~0.129 h | real |

(The two real components overlap slightly; they don't sum exactly to the
0.304 h that survives removing the artifact.)

- **Nearly half was the statistic, not the models.** `wockner-fit-kfold.R`
  reported `transform(mu_logit_cl)` for `no_pool` -- a logit-scale hyper-mean
  pushed through a nonlinear transform -- against `pooled_cl`'s single value
  in hours. `inv_logit` is concave at p ~ 0.67, so by Jensen the hours-scale
  average across trials sits *below* the transform of the mean: 45.316 h vs
  45.538 h (a 2nd-order expansion predicts -0.162 h, observed -0.222 h).
  **When comparing a hierarchical population location against a pooled value,
  compare the trial-averaged `cycle_length`, not `transform(mu_logit_cl)`.**
- **Larger trials have shorter cycles** (`corr(cl, n_obs)` = -0.33), so the
  observation-weighted average (45.187 h) sits below the unweighted one
  (45.316 h) and much closer to `pooled_cl`'s 45.012 h. The hierarchical mean
  weights trials ~equally; the pooled value is observation-weighted.
- **Nothing else moved.** `R`, `sd_iRBC`, `log10_total0` and `mu_logit_R` are
  effectively identical between the two fits. The only other parameter that
  shifted is `b_offset` (phase, mean diff 0.089, max 0.323) -- forcing one
  cycle length pushes per-trial timing differences into phase, which is
  where you'd expect them to go.
- Keep the magnitude in perspective: the surviving 0.304 h is well inside one
  posterior sd of either estimate (~0.6-0.9 h). It is systematic, not large.

### Per-trial cycle_length is almost entirely predicted by sampling density

From `wockner-cl-clustering.R`. Two findings, and the second matters much more
than the first.

**The apparent clustering is not real.** The `no_pool` per-trial posterior
means sort into what looks like two groups (seven at 44.1-45.0 h, a 1.06 h
gap, six at 46.0-46.6 h). It doesn't survive either check:

- Spacing structure is exactly what the fitted unimodal hierarchy produces.
  Comparing, within each draw, the 13 fitted values against 13 fresh values
  simulated from that draw's own `normal(mu_logit_cl, sigma_logit_cl)`:
  P(fitted more clumped than simulated) = **0.43**. Slightly *less* clumped
  than the model expects, not more.
- The gap is smaller than one per-trial posterior sd (1.06 h vs ~1.36 h). No
  trial sits confidently on either side (`P(cl > 45.5)` spans only 0.14-0.74),
  and the posterior for how many trials exceed the split is broad (modal value
  7 carries just 0.16). The split is an artifact of reading point estimates.

**But: `corr(per-trial cycle_length, observations per series) = -0.93.**
That is the correlation of posterior *means*, with n = 13. Computed instead
within each posterior draw, so that each trial's 1.15-1.51 h uncertainty is
carried rather than conditioned away, it is **-0.52 (95% -0.85 to +0.14)** --
the same sign, but an interval crossing zero. Quote both: the point-estimate
version describes the systematic component and overstates how sharp the
relationship is. The sparsest-sampled trial (OZ439, 4.8 obs/series)
has the highest estimate (46.6 h); the densest (MMV048_PartB, 7.7) has the
lowest (44.1 h). Related measures are far weaker -- `span_h` and `n_cycles`
-0.40, `n_obs` -0.33, `n_series` -0.03 -- so it is specifically **how densely
each series is sampled**, not how much data the trial has overall. That fits
the mechanism: cycle length is identified from the phase of the oscillation
*within* a series, and adding more sparsely-sampled series doesn't add phase
resolution the way adding observations within a series does.

Three candidate explanations. **1 is dead; 2 was tested directly and does not
hold in the form stated here** -- see "Schedule-bias simulation" below. The
list is kept because the reasoning for each is still what the tests were
built against.

1. **Prior pull.** Weakly identified trials drift toward the 48 h prior,
   which sits above every trial. Argues against it: per-trial posterior sds
   are nearly uniform (1.15-1.51 h) and correlate only +0.31 with the
   estimate, so the trials don't differ much in *uncertainty* even as they
   differ systematically in *location*.
2. **Likelihood-side bias from the observation schedule.** Sparse sampling
   aliases the oscillation, biasing the recovered phase and hence cycle
   length. This would not be fixed by changing the prior, which makes it the
   more worrying possibility.
3. **Confounding with study protocol.** The 72 h-span studies are different
   trials run differently; sampling density may proxy something real.

**This bears directly on the hierarchy question.** If between-trial variation
in cycle length is substantially an artifact of differing observation
schedules rather than biology, then the hierarchy is modelling measurement
design, and the earn-its-keep comparison is answering a different question
than intended.

The decisive test for 2 was run; see the next section. The prior-sensitivity
configs (entries 3-7 in `wockner-fit.R`) are still unsubmitted and would now
be confirming on real data what the simulation already showed on simulated.

### `cl_prior_center = 48` is doing real work, and probably shouldn't be

The bigger finding here. Both models put the *same* prior on their population
parameter -- `normal(mean_logit_cl, sd_logit_cl)`, i.e. centred at
`cl_prior_center` = **48 h** with `sd_logit_cl = 1` -- and the data sit at
~45 h. Under a normal-normal approximation (prior weight =
(posterior sd / prior sd)^2) that prior carries:

- `no_pool`'s `mu_logit_cl`: weight 0.079, pulling the estimate **+0.271 h** up
- `pooled_cl`'s `logit_cl`: weight 0.038, pulling it **+0.153 h** up

`mu_logit_cl` never enters the likelihood directly (it reaches the data only
through `eta_cl ~ normal(mu_logit_cl, sigma_logit_cl)`), so the same prior has
about twice the leverage there -- that asymmetry *is* the differential-pull
component above. But the part worth acting on is that **both** estimates are
being pulled up by a prior centred 3 h above where the data are, and
`sd_logit_cl = 1` is not weak enough for that to be ignorable. Nothing in the
settled-configuration history (`max_cl`, `sd_bs_cl`, `center_cl`) examined
`cl_prior_center`. Worth a sensitivity check -- refit with `sd_logit_cl` wide
or `cl_prior_center` at ~45 -- before any cycle-length estimate is reported
as a number rather than used for model comparison.

**Revised by the schedule-bias simulation.** On simulated data where the true
cycle length is known, widening `sd_logit_cl` from 1 to 2 moved the estimate
only 0.145 h, and the implied prior weight was 0.113 -- the same order as the
0.079 above, and confirming the prior is not the main problem. The estimate
was still +1.5 h high with the prior's influence removed. So moving
`cl_prior_center` off 48 h would buy less than this section implies, and
would not touch the larger, likelihood-side bias. Decide it on its own
merits, not as a fix for the cycle-length estimate.

### Schedule-bias simulation: explanation 1 is dead, 2 is not the answer either

`wockner-schedule-sim.R` (+ `.sh`), analysed by
`wockner-schedule-sim-analyze.R`, verified by `wockner-schedule-sim-check.R`.
Six fits: two prior arms (`sd_logit_cl` 1 and 2) x three noise replicates,
sharing noise seeds across arms so the prior comparison is paired. Outputs
are `_data/wock-schedsim-{fit,RES}-<arm>-rep<n>.rds`.

Design: simulate **every trial from one true `cycle_length`**, at the real
observation times and group structure, then fit `no_pool` and look for the
sampling-density correlation. There is no between-trial variation in the
truth, so anything that appears is manufactured by the design and the prior.
Generative parameters are the `pooled_cl` fit's posterior means, **all of
them** -- taking `no_pool`'s parameters and overwriting `cycle_length` would
pair per-trial phases with a pooled period, and `b_offset` is exactly where
forcing one cycle length pushes the timing (see the section above). True
`cycle_length` = 45.012 h.

**The headline is not the correlation. It is that the model does not recover
a cycle length it generated from.**

Cells: per-trial posterior mean `cycle_length` averaged over 13 trials, in
hours, with `(bias vs the true 45.012 h)` in brackets. The two columns are
the same simulated dataset fitted under two cycle-length prior widths, so
they are paired within a row.

| replicate | `sd_logit_cl` = 1 | `sd_logit_cl` = 2 |
|---|---|---|
| rep1 | 47.59 (+2.58) | did not converge |
| rep2 | 46.23 (+1.22) | 46.06 (+1.04) |
| rep3 | 47.21 (+2.19) | 47.09 (+2.08) |

- **Not a solver mismatch.** The simulation generates with `mat_exp_series`
  and the likelihood evaluates with `ew_poly_series`. `run_check = 1` puts
  `max_rel_diff` at 5e-13 over 150 draws, so they share a forward model. This
  check is the reason the bias can be read as a property of the design.
- **Not the prior.** Quadrupling the prior variance moved the estimate by
  0.145 h. That part is a direct measurement and stands.
  - An earlier version of this section went further, extrapolating the two
    arms through a normal-normal approximation to an "implied likelihood-only
    estimate" of 46.52 h, i.e. +1.5 h. **That number was wrong and is
    withdrawn.** `wockner-schedule-bias-profile.R` measured the likelihood
    directly instead of extrapolating to it, and found no such bias (below).
    The lesson is the obvious one: a two-point extrapolation through a
    nonlinear transform, on two replicates, was not evidence, and labelling
    it "read the sign, not the third digit" did not make it safe to report.
- **The noise draw matters far more than the prior.** Within a replicate the
  two arms agree closely; across replicates the population estimate moves
  ~1 h. Any single simulated dataset is a weak read.

**The correlation does reappear, but not demonstrably at full strength.**

Cells: Pearson correlation across 13 trials between each trial's posterior
mean `cycle_length` and its observations per series. Dimensionless; the
simulated row lists one value per replicate, not a range.

| | correlation of posterior means |
|---|---|
| real data | -0.933 |
| simulated, `sd_logit_cl` = 1 | -0.874, -0.354, -0.186 |
| simulated, `sd_logit_cl` = 2 | -0.893, -0.338 |

Zero of five replicates reached -0.933; the closest was -0.893. On the
within-draw statistic the real value sits at P ~ 0.16-0.19 against the
simulated distributions. So the schedules manufacture a correlation of the
right sign that sometimes gets close, and three noise realizations cannot
say whether they routinely reach -0.93.

**Do not mix the two correlation statistics.** Comparing the real
*point-estimate* correlation against a simulated *per-draw* distribution
guarantees an extreme-looking answer and means nothing; the per-draw version
is attenuated by posterior noise. An earlier version of the analysis script
did exactly that and reported P = 0.0000. It now reports the two separately
and refuses to cross them.

**What this changes.** The hierarchical fits are biased upward by +1.0 to
+2.6 h on data generated from a known cycle length, and widening the prior
barely touches it. Any reported cycle-length number inherits that. But
fixing it is not a matter of choosing a better prior, and it is not the
observation schedules aliasing the likelihood either -- see the next section.

### Where the bias is, by elimination

`wockner-schedule-bias-profile.R`, output `_data/wock-schedbias-profile.rds`.
Maximum likelihood at the real observation schedules: no prior, no hierarchy,
no MCMC. Simulate from the known truth, maximize, see where the maximum
lands. Cheap because all series in a trial share one trajectory (the same
dedup `archer_fit.stan` does), so one full-grid solve per likelihood
evaluation.

**The control is what licenses the rest**: with noiseless data the maximum
must land exactly on the truth, and it does, to 0.00e+00 h for all 13 trials.
The first version of this script failed that control on one trial, because
`n_full` was set to `max(ts)/dt_full` instead of
`round(max(ts)/dt_full) + 1` (`archer_fit.stan:179`), so every trial observed
at t = 216 indexed past the end of its trajectory and scored `1e12`
everywhere. Without the control that would have been a fabricated row in a
results table.

Cells: `bias` is the estimated `cycle_length` minus the true 45.012 h, in
hours, averaged over replicates; **positive means overestimation**. These are
maximum-likelihood point estimates with no prior and no MCMC, so `precision`
is the spread across independent simulated datasets, not a posterior width.

| estimator | bias vs truth 45.012 h | precision |
|---|---|---|
| per-trial MLE, averaged over 13 trials | **-0.195 h** | 16 replicates, se ~0.22 |
| pooled MLE, one `cycle_length` for all trials | **+0.544 h** | 8 replicates, 95% CI -0.26 to +1.34 |
| hierarchical posterior, `sd_logit_cl = 2` | +1.56 h | paired replicates 2-3 |
| hierarchical posterior, `sd_logit_cl = 1` | +1.71 h | paired replicates 2-3 |

- **The likelihood is not biased upward at these schedules.** Per trial the
  mean bias is indistinguishable from zero, and the pooled maximum is
  +0.54 h with a confidence interval spanning zero (p = 0.15, n = 8). There
  may be a small positive bias; there is certainly not a +1.5 h one.
- **The schedules do not manufacture the correlation at the likelihood
  level either.** `corr(per-trial MLE, obs_per_series)` is **+0.126** on the
  means and **+0.031 (95% -0.317 to +0.477)** within a replicate. The wrong
  sign, centred on zero. This is the cleanest test of explanation 2 run so
  far and it comes back negative.
- **A single trial barely identifies `cycle_length` at all.** Per-trial MLE
  sds are 1.9-4.6 h, against per-trial *posterior* sds of 1.15-1.51 h in the
  hierarchical fit. Most of what pins down a trial's cycle length in the
  fitted model comes from the other trials, not from that trial's own data.
- **The hierarchy is not where it lives either.** An earlier version of this
  section put ~1.1 h on the hierarchical structure by subtracting an
  unpaired likelihood figure from an unpaired posterior one. That arithmetic
  was invalid -- replicate-to-replicate spread is ~1 h, larger than the
  differences being separated, so only within-replicate comparisons can be
  subtracted. Fitting `pooled_cl` to the same simulated datasets (arm
  `no_hier`) gives the paired answer:

  Cells: cycle-length **bias** in hours -- posterior mean `cycle_length`
  minus the simulated truth of 45 h -- one row per noise replicate, so the
  four numbers in a row are paired on the same simulated dataset. `sd 1` and
  `sd 2` are the `sd_logit_cl` prior scale. `hierarchy` is `no_pool` sd 1
  minus `pooled_cl`, the part of the bias the hierarchical structure
  contributes; positive means biased high.

  | replicate | `no_pool` sd 1 | `no_pool` sd 2 | `pooled_cl` | hierarchy |
  |---|---|---|---|---|
  | rep2 | +1.218 | +1.044 | +0.928 | +0.290 |
  | rep3 | +2.194 | +2.077 | +1.968 | +0.226 |

  The hierarchy contributes ~0.25 h and the prior (sd 1 to 2) ~0.15 h.
  Removing the hierarchy entirely leaves +0.93 h and +1.97 h. So the bulk of
  the bias is in the single-cycle-length fit itself, and the replicate-to-
  replicate swing (+0.93 vs +1.97 on the same design, differing only in
  noise) is larger than every structural effect measured.

  This does not contradict the unbiased per-trial MLEs above: a pooled
  maximum weights trials by information rather than averaging them, and its
  own replicate spread is ~0.96 h. The paired pooled MLE on these exact two
  datasets is what splits the remainder into likelihood and prior;
  `wockner-schedule-bias-profile.R`'s last stage computes it.

  A third arm, `tight_sigma` (`no_pool` with `sd_bs_cl = 0.001`), was meant
  as the within-model control, since `no_hier` answers the question by
  changing the model. It was **cancelled**: it ran ~5x slower than every
  other arm -- a centred parameterization pinned at a near-zero scale is the
  geometry HMC handles worst -- and was heading for 4-6 h and probably a
  non-converged fit. The control it would have provided is instead a reading
  of the two Stan programs: `archer_fit_pooled_cl.stan` differs from
  `archer_fit.stan` only in the cycle-length block (a scalar `logit_cl`
  broadcast by `rep_vector` in place of `mu_logit_cl`/`sigma_logit_cl`/
  `eta_cl`, and the `n_grp_cl >= 2` reject dropped). Likelihood,
  `reduce_sum`, trajectory dedup, Erlang window and generated quantities are
  identical, and both put the same `normal(mean_logit_cl, sd_logit_cl)` prior
  on their population location. The arm is still defined in
  `wockner-schedule-sim.R` if it is ever wanted; note it needs
  `center_cl = 0` to stand a chance of sampling.

**Not a posterior-summary artifact.** `wockner-schedule-sim-mode.R`, output
`_data/wock-schedsim-mode.rds`. `cycle_length` is bounded and nonlinearly
transformed, so a skewed posterior summarized by its mean would manufacture
part of the gap. It does not: across the five converged fits the posterior
mean and mode of the population `cycle_length` differ by **-0.044 h**, mean
posterior skew is **-0.016**, and the bias is +1.82 h by the mean against
+1.87 h by the mode. The posterior is symmetric and the bias is in it, not in
the choice of summary. Same conclusion applying the transform after
summarizing `mu_logit_cl` rather than before.

### The simulation does not recover its own nuisance parameters

The finding that matters most here, and it was turned up by chasing the
cycle-length bias rather than looked for. `wockner-schedule-sim-analyze.R`'s
nuisance-recovery stage compares each simulation fit's posterior means against
the values the data were generated from:

Cells: `truth` and `estimated` are posterior means averaged over the
parameter's groups (14 `grp_init`, 13 `grp_R`, 27 `grp_sd`), on each
parameter's own scale; `difference` is `estimated - truth` in those units;
`relative` is that divided by `truth`, so **0% is perfect recovery**.
Averaged over 7 replicates that share the mean-vector truth.

| parameter | truth (mean) | estimated | difference | relative | prior |
|---|---|---|---|---|---|
| `b_shape` | 14.9 | 8.90 | **-5.95** | -40% | `lognormal(2, 0.5)`, median 7.39 |
| `log10_total0` | 0.425 | 0.854 | **+0.429** | +101% | `normal(1, 0.25)` |
| `b_offset` | 0.254 | 0.369 | +0.114 | +45% | -- |
| `R` | 6.24 | 4.73 | -1.50 | -24% | median `max_R * inv_logit(-2)` = 5.98 |
| `sd_iRBC` | 0.598 | 0.607 | +0.009 | +1% | -- |

`b_shape` lands almost exactly on its prior median instead of on the truth,
and `log10_total0` moves halfway to its prior mean. Only `sd_iRBC` is
recovered.

An earlier version of this section explained that as a three-way
`(log10_total0, R, b_shape)` ridge. **That was wrong**: `wockner-ridge.R`
finds `b_shape` essentially uncorrelated with the others in the real
posterior (median +0.026 against `log10_total0`, +0.004 against `R`, over all
14 groups). The ridge is `log10_total0` against `R` alone -- see the next
section.

This is also where the maximum-likelihood/posterior gap comes from. The paired
ladder, all on the same simulated datasets:

Cells: bias in hours -- estimated population `cycle_length` minus the true
45.012 h, **positive means overestimation**. All three columns are fitted to
the same simulated dataset within a row, so differences across a row are
paired and differences down a column are not.

| replicate | pooled MLE | `pooled_cl` posterior | `no_pool` posterior |
|---|---|---|---|
| rep1 | +1.01 | -- | +2.58 |
| rep2 | **-1.18** | +0.928 | +1.218 |
| rep3 | +0.711 | +1.968 | +2.194 |

On rep2 the data alone say 43.8 h and the fitted model says 45.9 h. That
2.1 h gap is not the cycle-length prior, which moves the estimate 0.15 h when
widened; it is the nuisance priors dragging `cycle_length` along the ridge,
plus the MLE fixing `sd_iRBC` at truth where the fit estimates it.

The truth used here was the `pooled_cl` fit's posterior **mean vector**, and
with correlated parameters a mean vector is not a coherent parameter set: it
can sit where no posterior draw sits. That was the leading suspect for the
mis-recovery, and the simulation was rebuilt from single posterior draws to
test it. **It is not the explanation** -- recovery is no better from a draw.
See "Rebuilt from posterior draws" below; the numbers in this table stand as
a property of the model at this design, not as an artifact of the truth
construction.

### The ridge is real, but the simulation's weak identification is not

`wockner-ridge.R`, output `_data/wock-ridge.rds`. Post-hoc on the saved real
`no_pool` fit, no refitting. Asks whether the nuisance mis-recovery above is a
property of the model on real data or an artifact of how the simulation's
truth was built.

**Posterior correlations, over all 14 `grp_init` groups:**

Cells: within-draw Pearson correlation between the two parameters, computed
separately in each of the 14 `grp_init` groups of one saved real-data
`no_pool` fit; `median` and `range` are over those 14 groups. Dimensionless.

| pair | median | range |
|---|---|---|
| `log10_total0` vs `R` | **-0.886** | -0.958 to -0.811 |
| `b_offset` vs `cycle_length` | +0.517 | **-0.932 to +0.667** |
| `b_shape` vs `log10_total0` | +0.026 | -0.050 to +0.064 |
| `b_shape` vs `R` | +0.004 | -0.030 to +0.066 |
| `b_shape` vs `cycle_length` | -0.094 | -0.267 to +0.083 |

- **`log10_total0` against `R` is a genuine ridge**, in every group, |r| ~0.89:
  a larger starting population with slower growth gives a similar trajectory.
- **`b_shape` is not in it.** It is independent of everything else, so no
  trade-off explains its mis-recovery in the simulation.
- **`b_offset` against `cycle_length` does not generalize.** Group 1 gives
  -0.932, which is where a single-group read would have stopped; the median
  over all groups is +0.52 and the sign flips. Check every group before
  quoting this pair.

**The real data are informative where the simulated data were not**, compared
on the scale each prior is written on (a lognormal's natural-scale sd grows
with its location, so `b_shape` must be compared on the log scale):

Cells: `real posterior` is the posterior mean over groups on the parameter's
natural scale; the last column is the posterior sd divided by the prior sd,
**each computed on the scale that prior is written on** (log scale for the
lognormal), so **a ratio near 1 means the data added nothing** and near 0
means the likelihood dominates.

| parameter | prior | real posterior | posterior sd / prior sd |
|---|---|---|---|
| `b_shape` | `lognormal(2, 0.5)`, median 7.39 | **18.9** | 0.759 (log scale) |
| `log10_total0` | `normal(1, 0.25)` | **0.312** | 0.697 |
| `R` | median `max_R * inv_logit(-2)` = 5.96 | 5.88 | 0.106 |

On real data the likelihood moves `b_shape` from 7.39 to 18.9 and
`log10_total0` from 1 to 0.312, well away from both priors. In the simulation,
generated from `b_shape` = 14.9, the fit returned 8.90 -- back at the prior
median. The simulated data are therefore less informative than the real data.

Two cautions on this comparison. First, it is across fits: 18.9 is the real
**`no_pool`** posterior, while the simulation's truth came from **`pooled_cl`**
(14.9). Second, the obvious explanation -- that the mean-vector truth broke
the parameter correlations -- has since been tested directly and **rejected**;
see the next section. Whatever makes the simulated data less informative than
the real data at the same observation times and the same sample size is not
yet identified.

Note what this does *not* say. The `log10_total0`/`R` ridge is real and
pervasive, and `R` sitting near its prior median is a coincidence of location,
not evidence the prior is driving it -- its posterior is a tenth of the prior
width.

### Rebuilt from posterior draws: the mean-vector explanation is dead

`wockner-schedule-sim.R` with `SCHEDSIM_TRUTH_DRAW=<i>` in the environment
takes draw `i` of `wock-fit-pooled_cl.rds` as the truth instead of the
posterior mean vector, writes `truth_source` and `truth_nuisance` into the
summary, and names the output `default-rep2-draw<i>`. Three draws were run
(500, 1500, 2500), all on replicate 2's noise seed so the only thing that
changes against `default-rep2` is the truth. Output
`_data/schedsim-analyze.log`.

`draw1500` is **excluded**: R-hat 1.107, 301 divergences (7.5%), min ESS 26.
That leaves two usable draw-based replicates, which is thin, and every number
below should be read with that in mind.

**Nuisance recovery does not improve.** Relative difference, posterior mean
against the truth actually simulated from:

Cells: `(estimate - truth) / truth` as a percentage, where the estimate is
the posterior mean averaged over the parameter's groups and the truth is what
that replicate was simulated from; **0% is perfect recovery**. The mean-vector
column averages 7 replicates, the draw columns are single replicates.

| parameter | mean vector (n=7) | draw 500 | draw 2500 |
|---|---|---|---|
| `b_shape` | -40% | -26% | -22% |
| `log10_total0` | +101% | **+110%** | **+133%** |
| `R` | -24% | -24% | -27% |
| `b_offset` | +45% | +199% | +45% |
| `sd_iRBC` | +1% | -0.2% | -0.7% |

`b_shape` recovers somewhat better and `log10_total0` recovers *worse*. Only
`sd_iRBC` is recovered under either truth. A coherent parameter vector was the
leading explanation for the mis-recovery and it does not survive contact with
the test.

**The cycle-length bias survives:**

Cells: `truth` is the single `cycle_length` all 13 trials were simulated
from, in hours; `fitted mean` is the per-trial posterior mean averaged over
13 trials, in hours; `bias` is `fitted mean - truth`, so **positive means the
fit overestimates**. One simulated dataset per row.

| replicate | truth | fitted mean | bias |
|---|---|---|---|
| `default-rep1` | mean vector, 45.012 | 47.59 | +2.58 |
| `default-rep2` | mean vector, 45.012 | 46.23 | +1.22 |
| `default-rep3` | mean vector, 45.012 | 47.20 | +2.19 |
| `default-rep2-draw500` | draw, 45.443 | 47.12 | **+1.68** |
| `default-rep2-draw2500` | draw, 45.426 | 46.28 | **+0.86** |

Both draw-based values sit inside the replicate-to-replicate spread of the
mean-vector runs (+1.22 to +2.58), so the truth construction is not
distinguishable as a source of bias at n=2. What is now established is the
negative: the bias is **not** an artifact of simulating from an incoherent
parameter vector.

**Sampling-density correlation, per-trial posterior means against observations
per series:**

Cells: Pearson correlation across the 13 trials between each trial's
posterior mean `cycle_length` and its observations per series. Dimensionless,
one value per fit; **more negative means sparser-sampled trials got longer
cycle lengths**.

| source | r |
|---|---|
| real data | **-0.933** |
| `default`, mean vector (n=3) | -0.186, -0.354, -0.874 |
| `wide`, mean vector (n=2) | -0.338, -0.893 |
| `default-rep2-draw500` | -0.634 |
| `default-rep2-draw2500` | -0.557 |

0 of 7 simulated replicates reach the real value. Read alongside the per-draw
statistic, which carries each trial's uncertainty: real -0.515 (95% -0.852 to
+0.137) against simulated -0.13 to -0.20 depending on arm. The two statistics
are not interchangeable and are reported separately for that reason -- see the
comment block in `wockner-schedule-sim-analyze.R`.

So the schedules plus the model manufacture a correlation of roughly -0.2 to
-0.6 out of nothing, which is a large share of -0.93 but does not account for
it.

### Inoculum-anchored prior for `log10_total0` (package change)

Idea: `y` is in infected RBC per mL and the inoculation sizes are known in
viable parasites, so `log10(inoculum / 5000 mL)` is a per-`grp_init`
prediction of `log10_total0` rather than a guess. It replaces
`normal(mean_log10_total0, sd_log10_total0)` with a hierarchy centred on the
anchor:

```stan
delta_total0[1] ~ normal(mean_delta_total0, sd_delta_total0);
sigma_total0[1] ~ normal(0, sd_bs_total0) T[0,];
log10_total0 ~ normal(anchor_log10_total0 + delta_total0[1], sigma_total0[1]);
```

`delta_total0` is a single shared offset on the log10 scale, estimated, not
fixed. It is unbounded by design: a fraction bounded at 1 would be violated by
anyone whose blood volume is under 5 L, and the establishment fraction has no
reliable published value to fix it to. `delta_total0` and `sigma_total0` are
declared `array[total0_anchor]`, so they are **zero-sized and cost nothing**
when the anchor is off, which is the default.

What `wockner-inoc-prior.R` (output `_data/wock-inoc-prior.rds`) found on the
real fit, and what the anchored fit has to reproduce or contradict:

- The fit wants **+0.811 log10 more** parasites at t=0 than were inoculated,
  a factor of 6.47. That is the value `delta_total0` should land on.
- Blood volume cannot absorb it. +-10% moves the anchor 0.041 log10; closing
  0.811 would need a blood volume of 772 mL.
- Between-group sd of `log10_total0` is 0.165 against a within-group posterior
  sd of 0.181, ratio 0.91. One shared offset is compatible with the data;
  a per-group offset is not required.
- Implied `R` at the anchor is **9.96**, from regressing `log10(R)` on
  `log10_total0` across draws (slope -0.25 along the ridge). An earlier
  arithmetic estimate of ~20 ignored sequestration and is wrong. 9.96 sits
  between *in vitro* 3D7 (~8) and the burst-size ceiling (median 15-18,
  max 32), and well under Wockner's 28.7-35.4.

`max_R` stays at 50. It is deliberately permissive so the model does not
structurally exclude the Wockner range, which would make the two sets of
estimates incomparable.

Package state: `plasmofit` **0.0.0.9008**, installed and tested (171 passing,
0 failing). Documentation regenerated at 0.0.0.9009 with roxygen2 8.1.0;
`tools::checkDocFiles()` and `tools::undoc()` are both clean. The change
touches all four Stan programs, `R/archer-fit.R`, and
`tests/testthat/test-archer-stan-data.R`. Validation run: anchor off by
default with an all-zero anchor vector; `inoc_size = "inoc_size"` reproduces
`wockner-inoc-prior.R`'s independent per-group prediction; `blood_volume_ml`
rescales by exactly log10 of the ratio; all four bad-input paths error.

### Regression test for the anchor

`_scripts/wockner-anchor-regression.R`, output `_data/anchor-regression.log`.
The question is whether `total0_anchor = 0` under 0.0.0.9008 still targets the
posterior that the pre-anchor code did. It matters because the `else` branch
is supposed to be the original prior statement verbatim and the two new
parameters are supposed to be zero-sized.

**Do not test this by comparing draws.** The first version of this script
demanded bit-identical draw matrices and reported FAIL. That bar is
unreachable: the reference was sampled by a **different compiled DSO**, and
recompiling can reorder floating-point operations, which HMC amplifies into a
completely different trajectory within a few leapfrog steps. A second version
compared posterior sds against a flat 0.9-1.1 band and also reported FAIL;
that band is arbitrary in the other direction, because the Monte Carlo error
on a posterior sd scales with effective sample size, and the parameters that
tripped it (`b_off_vec`, ESS 500-800) carry several percent of it per fit.
Neither FAIL was evidence about the code.

**The design that does work is three fits:**

| fit | code | seed | file |
|---|---|---|---|
| reference | pre-change | A | `_data/regress-ref-fit-default-rep2.rds` |
| rebuilt | post-change | A | `_data/wock-schedsim-fit-default-rep2.rds` |
| null run | post-change | B | `_data/wock-schedsim-fit-default-rep2-seed*.rds` |

reference-vs-rebuilt is the test. rebuilt-vs-null-run is the null: identical
code and identical data, different sampler seed, so everything it shows is
run-to-run variation. Splitting one fit's chains is **not** a substitute for
it -- that holds the step-size and mass-matrix adaptation constant, and those
adapt separately in every run (0.0210 against 0.0203 here, with 1.04M against
0.90M leapfrog steps and 49 against 25 divergences).

`SCHEDSIM_FIT_SEED=<n>` in `wockner-schedule-sim.R` refits the same simulated
data with a different sampler seed and names the output `...-seed<n>`:

```
sbatch --array=2 --job-name=wock-seednull \
    --output=_data/wock-seednull.out --error=_data/wock-seednull.err \
    --export=ALL,SCHEDSIM_FIT_SEED=271828183 \
    _scripts/wockner-schedule-sim.sh
```

**Result so far, means (the substantive half, and it passes):**

- parameter name sets identical, 208 each -- so the zero-sized
  `delta_total0` and `sigma_total0` really do contribute no columns
- `|z|` on posterior means, scaled by each fit's MCSE: median 0.73,
  90th percentile 1.60, max 2.67; **0 of 208 above 3**
- `log10_total0`, the one parameter whose prior statement the change touches:
  max `|z|` 1.79, sd ratio 0.964-1.049
- widest movers are `sd_iRBC[2]` (z 2.67) and `b_shape[4]` (z -2.34), neither
  connected to the anchor

**Widths, against the null run (seed 271828183), and the verdict:**

Cells: for each of 208 parameters, `|log(sd_A / sd_B)|` where sd is the
posterior sd; the columns are the median and 95th percentile of that over the
208. Dimensionless; **0 means the two fits give identical posterior widths**.

| | median \|log sd ratio\| | 95th percentile |
|---|---|---|
| test, pre-change vs post-change | 0.0262 | 0.1248 |
| null, same code, different seed | 0.0368 | 0.1241 |

Inflation at the 95th percentile: **1.01x**, against a threshold of 1.5x fixed
before the null run existed. Null run health: max R-hat 1.0287, 39
divergences, so it converged and is entitled to serve as the null. Means on
the null pair behave like the test pair (median \|z\| 1.12 against 0.73, 0 of
208 above 3 in both).

**Verdict: PASS.** The anchor switched off targets the same posterior as the
pre-change code.

Note how badly the two earlier bars misled. Bit-identity said FAIL. The flat
0.9-1.1 band on sd ratios said FAIL. The within-fit split-half null put the
inflation at 1.41x, which also would have said FAIL at any sensible
threshold — and the reason it is wrong is visible in the numbers: two halves
of one fit share an adapted step size and mass matrix, so they agree more
closely than two independent runs do. Only the between-run null is the right
comparison, and it gives 1.01x.

The script also refuses to pass on a null run with R-hat above 1.05. A badly
mixed null has inflated widths and would make **any** test pair look
acceptable, so a pass obtained against one means nothing. In that case it
reports INCONCLUSIVE and asks for another `SCHEDSIM_FIT_SEED`.

### Information content: the simulated designs are not the poorer ones

`wockner-sim-information.R`, output `_data/sim-information.log` and
`_data/wock-sim-information.rds`. Open thread 3, answered and closed.

The lead was that simulated `log10(y + 1)` has sd 0.85-0.89 against a real
1.01, which looked like it might explain why the simulation recovers
`b_shape` low and `log10_total0` high. **It does not.** Raw spread is the
wrong quantity: information about a periodic parameter comes from the
oscillation, not from the trend or the overall scale. Detrending each series
on time and dividing the remaining signal by that series' `sd_iRBC` gives

Cells: `signal` is the median over 177 series of the sd of the **noiseless**
trajectory after removing a within-series linear trend in time, in log10
units; `noise` is the median over series of that series' `sd_iRBC`, same
units; the fourth column sums `(signal/noise)^2` weighted by each series'
observation count over all 1130 observations, so it is dimensionless and
**larger means more information about the oscillation**; `vs real` is that
sum divided by the real data's.

| dataset | signal (wiggle) | noise | sum (signal/noise)^2 | vs real |
|---|---|---|---|---|
| real | 0.204 | 0.585 | 149 | -- |
| `sim-rep1/2/3` | 0.207 | 0.590 | 145 | **0.97** |
| `sim-draw500` | 0.229 | 0.608 | 195 | 1.31 |
| `sim-draw1500` | 0.230 | 0.606 | 197 | 1.32 |
| `sim-draw2500` | 0.218 | 0.624 | 174 | 1.17 |

The apples-to-apples comparison is `sim-rep1/2/3` against real, because both
build their trajectory from a posterior **mean**: ratio **0.97**, which is no
deficit at all. The draw-based replicates come out 1.17-1.32 precisely
because a draw is less shrunk than a mean, which is a useful internal check
that the machinery is measuring what it claims to.

So the mis-recovery is not an information deficit and needs another
explanation. **The natural successor hypothesis: the estimator is biased at
this design, and the real fit is subject to the same bias.** If simulating
from `b_shape` = 14.9 returns 8.90, a 40% shortfall, then the real fit's 18.9
is not evidence that the real data pin `b_shape` down — it is an estimate
from the same biased estimator, and the true value would be higher still.
That reframes the ridge section's "the real data are informative where the
simulated data were not", which compared a real estimate against a simulated
truth as though the former were unbiased.

**A separate finding, not what was being looked for.** The real data's extra
spread is entirely within series (ratio 0.80-0.86) and not between them
(0.83-1.01), while the *detrended* within-series scatter is if anything
larger in simulation (0.53-0.57 against a real 0.508). All of the real
excess therefore sits in the within-series **linear trend**: real series rise
and fall more steeply over time than the fitted trajectories do, with a
median within-series range of 2.56 log10 units against 1.95-2.21 simulated.
That is a lack-of-fit signal in the model's time course, independent of
everything above, and nobody has looked at it.

### The `log10_total0` prior was misspecified, and that is most of the story

`wockner-anchor-check.R`, output `_data/anchor-check.log` and
`_data/wock-anchor-check.rds`. Three fits to the **real** data.

Cells: each parameter column is the **posterior mean, then averaged over that
parameter's groups** -- 14 `grp_init` for `log10_total0`, 13 `grp_R` for `R`,
13 `grp_cl` for `cycle_length`, 27 `grp_sd` for `sd_iRBC`. `log10_total0` is
in log10 iRBC/mL, `cycle_length` in hours, `R` and `sd_iRBC` dimensionless.
`loo elpd` is expected log pointwise predictive density over all 1130
observations; **higher (less negative) is better**. One fit per row, no
replicates.

| config | prior on `log10_total0` | `log10_total0` | `R` | `sd_iRBC` | `cycle_length` | loo elpd |
|---|---|---|---|---|---|---|
| `no_pool` | `normal(1, 0.25)` | +0.430 | 6.23 | 0.597 | 45.32 | -1053.2 |
| `np_wide_total0` | `normal(1, 1)` | -0.913 | 14.9 | 0.544 | 45.52 | **-950.7** |
| `np_anchor` | anchored on the inoculum | -1.18 | 17.9 | 0.543 | 46.14 | **-948.5** |

`loo_compare`, over the same 1130 observations. Cells: `elpd_diff` is each
row's elpd minus the best row's, so it is 0 for the best and negative for the
rest; `se_diff` is the standard error **of that paired difference**, which is
much smaller than the se of either elpd because the two share observations.

| | elpd_diff | se_diff |
|---|---|---|
| `np_anchor` | 0.00 | -- |
| `np_wide_total0` | -2.13 | 2.91 |
| `no_pool` | **-104.70** | **12.83** |

**Read this the right way round.** The anchored fit beats the original by 8.2
standard errors, but it does **not** beat the control that merely widens the
old prior without using the inoculum at all (2.13 +- 2.91). So the entire
predictive gain comes from `normal(1, 0.25)` being **wrong**, and the
inoculum information adds nothing detectable. The anchor was worth building
because it is what exposed this, not because it wins.

Two things independent of elpd say the new location is the right one:

- The original fit implied **6.47x more parasites at t = 0 than were
  inoculated** -- an establishment fraction of 647%, which is impossible. The
  anchored fit gives `delta_total0` = **-0.798**, a factor of 0.16, i.e. a
  16% establishment fraction, which is unremarkable.
- `R` moves from 6.23 into the **15-18** range, which is where burst size
  puts it (median 15-18 per schizont, max 32). The old 6.23 was near the
  *in vitro* 3D7 figure of ~8.

**The prediction of `delta_total0` = +0.811 failed, and the failure is the
finding.** It assumed the anchored fit would leave `log10_total0` where the
unanchored fit put it. It did not: `log10_total0` moved -1.61. The ridge
itself held up exactly -- applying the measured slope
(`d log10(R) / d log10_total0` = **-0.280**) to the distance actually
travelled predicts `R` = **17.63** against an observed **17.95**. What was
wrong was the assumed landing point, which took a prior-dominated posterior
for a likelihood-dominated one.

**Consequence: `R` is not identified by these data.** It reads 6.2, 14.9, or
17.9 depending on the prior, and the last two are predictively
indistinguishable. Any `R` reported from this model is a statement about the
prior unless the prior is defended. The same caution applies, more weakly, to
`cycle_length`, which spans 45.3 to 46.1 across the three.

This also retires the ridge section's claim that the real data move
`log10_total0` "well away from both priors". They do not: relax the prior and
it moves another 1.3 decades.

### Thread 2: it is the `b_shape` prior, and the ridge is not the mechanism

From `wockner-schedule-sim-analyze.R`, paired within replicate because
replicate-to-replicate spread (~1 h) exceeds the effects being measured. All
arms run on the same simulated datasets as `default`.

**Change in cycle-length bias when a prior is widened.** Cells: bias is the
per-trial posterior mean `cycle_length` averaged over 13 trials, minus the
known true `cycle_length`, in **hours**. Each entry is that arm's bias minus
`default`'s bias **on the same simulated dataset** (noise seeds are shared),
so it is paired; **negative means widening reduced the bias**. `mean change`
averages the three replicates; `per replicate` lists them in order 1, 2, 3.

| arm | widened | mean change | per replicate |
|---|---|---|---|
| `wide_bshape` | `sd_log_b_shape` 0.5 -> 1.5 | **-0.501 h** | -0.531, -0.452, -0.519 |
| `wide_nuis` | both | -0.540 h | -0.142, -1.08, -0.398 |
| `wide_total0` | `sd_log10_total0` 0.25 -> 1 | +0.129 h | +0.429, -0.156, +0.113 |

**Recovery of the nuisance parameters themselves, by arm.** Cells: posterior
mean averaged over the parameter's groups (14 `grp_init` for `b_shape` and
`log10_total0`, 13 `grp_R` for `R`), then averaged over the 3 replicates of
that arm; the bracketed percentage is `(estimate - truth) / truth`, so **0%
is perfect recovery** and the sign says which way it misses. `truth` is the
value the data were simulated from, identical across arms.

| parameter | truth | `default` | `wide_bshape` | `wide_total0` |
|---|---|---|---|---|
| `b_shape` | 14.9 | 8.89 (-40%) | **23.9 (+61%)** | 8.68 (-42%) |
| `log10_total0` | 0.425 | 0.858 (+102%) | 0.858 (+102%) | **0.533 (+25%)** |
| `R` | 6.24 | 4.75 (-24%) | 4.71 (-24%) | **5.92 (-5%)** |

The dissociation is clean and it is the opposite of what was expected:

- Widening the **`log10_total0`** prior fixes `log10_total0` (+102% to +25%)
  and `R` (-24% to -5%) -- and does **nothing** to the cycle-length bias
  (+0.13 h, sign inconsistent across replicates).
- Widening the **`b_shape`** prior does nothing for `log10_total0` or `R`,
  overshoots `b_shape` itself (-40% to +61%), and removes **0.50 h** of
  cycle-length bias, with a spread of only 0.08 h across three replicates.

**So the `log10_total0`/`R` ridge is not what drags `cycle_length`.** The
`b_shape` prior is, and `b_shape` is off the ridge entirely -- median
posterior correlation +0.026 with `log10_total0` and +0.004 with `R` (see
"The ridge is real"). The standing hypothesis in this file, that the bias came
from nuisance priors dragging `cycle_length` **along the ridge**, is wrong.

The real-data fits agree on the part they can speak to: relaxing the
`log10_total0` prior there moved `cycle_length` +0.20 h, against +0.13 h in
the simulation.

**Budget for the cycle-length bias** at the `default` arm, whose mean bias is
+1.97 h over replicates 1-3. Cells: `size` is hours of cycle-length bias
attributable to that source, each from a paired within-replicate contrast
against `default` on the same simulated datasets. These are not guaranteed to
sum to the total, and do not.

| source | size | how measured |
|---|---|---|
| `b_shape` prior | ~0.50 h | paired, `wide_bshape` vs `default` |
| hierarchy | ~0.25 h | paired, `no_hier` vs `default` |
| cycle-length prior | ~0.15 h | paired, `wide` vs `default` |
| `log10_total0` prior | ~0 | paired, `wide_total0` vs `default` |
| **unexplained** | **~1.07 h** | the remainder |

Over half is still unaccounted for.
