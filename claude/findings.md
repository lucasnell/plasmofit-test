# Findings: the hierarchy, the schedule-bias simulation, and the priors

*Detail file. Orientation and settled decisions are in the repo-root
`CLAUDE.md`; current status in `PROJECT_INDEX.md`; live work in `TODO.md`;
the last session in `handoff.md`. (These four replaced `claude/CLAUDE.md`
as the top layer on 2026-10-07.)*

Every numeric table below states what is in its cells; the rule is in the
repo-root `CLAUDE.md`. Roughly chronological, so a later section can overturn
an earlier one -- `PROJECT_INDEX.md` says which conclusions
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

Cells: `real posterior` is the posterior mean **in `grp_init` group 1 only**,
on the parameter's natural scale -- `wockner-ridge.R:49` fixes `gi` at group
1, which is the group `b_shape`, `b_offset` and `log10_total0` are indexed
by. It is **not** a mean over the 14 groups, as an earlier version of this
caption said; over all 14, `b_shape` is **14.91** (range 10.99-19.23), from
`_data/wock-prior-panel.rds`. The last column is the posterior sd divided by
the prior sd, group 1 only, **each computed on the scale that prior is
written on** (log scale for the lognormal), so **a ratio near 1 means the
data added nothing** and near 0 means the likelihood dominates.

| parameter | prior | real posterior | posterior sd / prior sd |
|---|---|---|---|
| `b_shape` | `lognormal(2, 0.5)`, median 7.39 | **18.9** | 0.759 (log scale) |
| `log10_total0` | `normal(1, 0.25)` | **0.312** | 0.697 |
| `R` | median `max_R * inv_logit(-2)` = 5.96 | 5.88 | 0.106 |

On real data the likelihood moves `b_shape` from 7.39 to 18.9 and
`log10_total0` from 1 to 0.312, well away from both priors. In the simulation,
generated from `b_shape` = 14.9, the fit returned 8.90 -- back at the prior
median. The simulated data are therefore less informative than the real data.

**Both halves of "well away from both priors" have since been retired.** For
`log10_total0`, by the misspecification section below. For `b_shape`, by the
eight-fit panel: its posterior CV over all 14 groups is 0.433 against a prior
CV of 0.533, so the default prior is doing nearly all the work, and widening
it moves the posterior mean from 14.9 to 65.4. This section's group-1 read of
18.9 is the **highest** of the 14 groups, which is what made the move look
larger than it is -- the single-group trap this section warns about for
`b_offset`, in its own table.

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
| `b_shape` prior | ~0.50 h **(superseded: 0.43 h)** | paired, `wide_bshape` vs `default` |
| hierarchy | ~0.25 h | paired, `no_hier` vs `default` |
| cycle-length prior | ~0.15 h | paired, `wide` vs `default` |
| `log10_total0` prior | ~0 | paired, `wide_total0` vs `default` |
| **unexplained** | **~1.07 h** | the remainder |

Over half is still unaccounted for.

**Superseded in two places by later sections.** The `b_shape` prior's share
is **0.43 h** at n=4, not 0.50 h at n=3, with a spread four times wider (see
"Refits: the budget reaches n=4"). And the one-at-a-time decomposition is the
wrong shape: the two nuisance priors are super-additive on `cycle_length`
(see "the eight-fit prior panel"). The unexplained remainder is unchanged and
is now ~1.5 h on the real-data scale (see "`b_shape` is not identified").

### Threads 2 and 3 on real data: the eight-fit prior panel

`_scripts/wockner-prior-panel.R`, output `_data/prior-panel.log`, summaries
`_data/wock-prior-panel.rds`. Fits are `_scripts/wockner-fit.R` tasks 1, 2,
8-13; tasks 10-13 are SLURM 28931, all four `COMPLETED`. The supplementary
spread check below is `_scripts/panel-bshape-spread.R`, output
`_data/panel-bshape-spread.log`.

**What licenses the rest.** Cells: one row per fit; `div` is divergent
transitions after warmup summed over 4 chains and `div_pct` the same as a
percent of the 4000 post-warmup transitions, `max_rhat` the largest split
R-hat over all parameters, `min_ess` the smallest bulk n_eff over all
parameters. Absolute, not a comparison.

| fit | div | div_pct | max_rhat | min_ess |
|---|---|---|---|---|
| `no_pool` | 42 | 1.05 | 1.01 | 246 |
| `np_wide_total0` | 82 | 2.05 | 1.02 | 160 |
| `np_wide_bshape` | 45 | 1.12 | 1.01 | 324 |
| `np_wide_both` | 137 | 3.42 | 1.02 | 214 |
| `np_anchor` | 74 | 1.85 | 1.03 | 186 |
| `pooled_cl` | 45 | 1.12 | 1.01 | 569 |
| `pl_wide_total0` | 78 | 1.95 | 1.01 | 205 |
| `pl_wide_both` | 160 | 4.00 | 1.02 | 230 |

All eight clear the R-hat < 1.05 gate, so nothing is dropped. **The
divergences are the standing caveat**: they are non-zero everywhere and rise
with each widening, to 3.4% and 4.0% in the two fits with both priors
widened. Every posterior mean below is a mean over an imperfectly explored
posterior, and the two `wide_both` fits are the ones to distrust first.

The panel reproduces both numbers it inherits -- the pooling offset (-0.304
h, paired mean z -3.79) and the original `no_pool`/`pooled_cl` elpd
difference (-1.38, se 1.03) -- from the fits that already existed.

#### Thread 2: the `b_shape` prior carries cycle length on real data too

Two contrasts, because the isolated effect and the effect on top of a
corrected `log10_total0` prior are different questions. Cells: `from` and
`to` are posterior means averaged over the parameter's groups (14 `grp_init`
for `log10_total0` and `b_shape`, 13 for `R` and `cycle_length`, 27 for
`sd_iRBC`); `shift` is `to` minus `from` averaged over groups and **paired by
group**; `mean_z` is that paired shift over its Monte Carlo standard error
averaged over groups, so it says the shift beats MCMC noise and **nothing
more** -- it is not an inferential z. Units: `cycle_length` in hours,
`log10_total0` and `sd_iRBC` in log10 units, `R` and `b_shape`
dimensionless.

| contrast | parameter | from | to | shift | mean_z |
|---|---|---|---|---|---|
| `no_pool` -> `np_wide_bshape` | `cycle_length` | 45.3 | 44.8 | **-0.468** | -8.41 |
| | `b_shape` | 14.9 | 65.4 | +50.5 | 40.0 |
| | `log10_total0` | +0.430 | +0.431 | +0.0007 | -0.12 |
| | `R` | 6.23 | 6.24 | +0.0095 | 0.77 |
| | `sd_iRBC` | 0.597 | 0.584 | -0.0135 | -7.33 |
| `np_wide_total0` -> `np_wide_both` | `cycle_length` | 45.5 | 44.2 | **-1.29** | -18.2 |
| | `b_shape` | 14.4 | 64.4 | +50.0 | 42.1 |
| | `log10_total0` | -0.913 | -0.881 | +0.0316 | 2.93 |
| | `R` | 14.9 | 14.4 | -0.447 | -5.03 |
| | `sd_iRBC` | 0.544 | 0.533 | -0.0115 | -7.47 |

Paired PSIS-LOO, `to` minus `from`, summed over the 1130 observations both
fits share; positive favours `to`; `se_diff` is the standard error of the
paired difference.

| contrast | elpd_diff | se_diff |
|---|---|---|
| `no_pool` -> `np_wide_bshape` | **+24.8** | 2.23 |
| `np_wide_total0` -> `np_wide_both` | **+24.4** | 2.06 |

**The simulation's dissociation reproduces.** Widening `b_shape`'s prior
moves `cycle_length` by -0.47 h in isolation, against -0.501 h in
simulation, and leaves `log10_total0` and `R` untouched (z -0.12 and 0.77) --
the same off-the-ridge behaviour the simulation found. Widening it is also a
genuine predictive gain, +24.8 elpd at 11 standard errors, which the
simulation could not have shown.

**The effect is larger once `log10_total0` is fixed**, -1.29 h against -0.47
h. The two priors are additive in elpd -- widening `log10_total0` alone gains
102.6, `b_shape` alone 24.8, both 127.0 -- but **super-additive on
`cycle_length`**: singly +0.20 and -0.47 h, jointly -1.1 h. So the joint
effect is not recoverable from the one-at-a-time contrasts, and the ~0.50 h
the simulation budget assigns to `b_shape` is a floor for what the pair does
together.

**This is a shift, not a bias.** There is no truth on real data. It bounds
from below how much of the reported 45.3 h is the nuisance priors' doing --
at least 1.1 h of it -- and says nothing about whether the remainder is
biological.

**`b_shape` is not identified in location, only bounded.** From
`_data/panel-bshape-spread.log`. Cells: `mean` and `sd` are the posterior
mean and posterior sd of `b_shape` averaged over its 14 `grp_init` groups;
`CV` is the mean over groups of each group's own sd/mean, comparable to the
prior CV, which for `lognormal(mu, sigma)` is `sqrt(exp(sigma^2) - 1)` and
does not depend on `mu`. Dimensionless.

| `b_shape` prior | prior CV | posterior mean | posterior sd | posterior CV |
|---|---|---|---|---|
| `lognormal(2, 0.5)` (default) | 0.533 | 14.9 | 6.5 | 0.433 |
| `lognormal(2, 1.5)` (widened) | 2.913 | 65.4 | 50.7 | 0.776 |

Under the default prior the posterior is barely narrower than the prior
(0.433 against 0.533) and sits at its upper end -- the data add almost
nothing, and **the 14.9 reported throughout this project is a prior
artefact**. Widening it lets the likelihood speak, and it does: the CV falls
from 2.913 to 0.776, so the data genuinely rule out the wide prior's tails.
But the location moves by a factor of 4.4, and a posterior sd of 50 on a mean
of 65 means that mean is not a usable point estimate. `b_shape` is in the
same position as `R`: **a number reported from it is a statement about the
prior.**

**Ruled out as evidence: the runtime signal.** Both `b_shape`-widened tasks
ran ~68% slower per iteration than `np_wide_total0`, which would be
consistent with widening `b_shape` making the geometry harder. Section 1 does
not support it. In isolation (`no_pool` -> `np_wide_bshape`) divergences are
flat (1.05% -> 1.12%), max treedepth hits fall 1 -> 0, and the adapted step
size *rises* 0.0209 -> 0.0219, all of which say the geometry got no harder.
The quantity that would settle it, mean leapfrog steps per iteration, is not
recorded by the panel. **It does not count toward anything.**

#### Thread 3: the hierarchy conclusion survives, its number does not

Cells: `no_pool` and `pooled_cl` are posterior mean `cycle_length` averaged
over 13 trials, in hours; `offset_h` is `pooled_cl` minus `no_pool` averaged
over trials and **paired by trial**; `mean_z` is that paired offset over its
Monte Carlo standard error, again a noise check and not an inferential z.
`elpd_diff` is paired PSIS-LOO, `pooled_cl` minus `no_pool`, summed over
1130 observations, **positive favouring `pooled_cl`**, with the standard
error of the paired difference.

| prior setting | `no_pool` | `pooled_cl` | offset_h | mean_z | elpd_diff | se_diff |
|---|---|---|---|---|---|---|
| as originally run | 45.3 | 45.0 | -0.304 | -3.79 | -1.38 | 1.03 |
| `log10_total0` corrected | 45.5 | 45.4 | **-0.110** | -2.16 | -1.67 | 0.86 |
| both corrected | 44.2 | 43.7 | **-0.509** | -7.19 | -1.68 | 1.21 |

**The predictive conclusion is robust.** `pooled_cl` fails to beat `no_pool`
at all three prior settings, by -1.4 to -1.7 elpd, never more than about 2
standard errors and always in the same direction. Correcting the prior that
costs 104.7 elpd does not rescue the hierarchy, which is the thing thread 3
asked.

**The offset's magnitude is not robust.** Its sign is -- pooling shortens
cycle length at every setting -- but the size runs -0.304, -0.110, -0.509 h,
a factor of 4.6, and **non-monotonically**: correcting `log10_total0` alone
nearly removes it, correcting both restores it larger than it started.
Quoting "-0.304 h" as the pooling offset is quoting a prior-dependent number
drawn under a prior that costs 104.7 elpd. On the best-fitting setting
available the offset is **-0.509 h**, and that is the number to use.

**Keep the offset next to a posterior width.** From
`_data/panel-bshape-spread.log`, the within-fit posterior sd of one trial's
`cycle_length`, averaged over 13 trials, is 1.36 h under `no_pool` and 0.63 h
under `pooled_cl` (1.28 and 0.56 h with both priors corrected). The -0.509 h
offset is ~0.4 of one trial's posterior sd. Its `mean_z` of -7.19 says only
that it is not MCMC noise; with one dataset and one fit per setting there is
no interval on it.

### The trend misfit is the `log10_total0` prior, not `b_shape`

`_scripts/wockner-ppc-trend.R`, output `_data/ppc-trend.log` and
`_data/wock-ppc-trend.rds`. SLURM 28957, COMPLETED. Post-hoc on four saved
real-data fits, no refitting.

The hypothesis under test was that the panel's `b_shape` ~ 65 is the model
absorbing structural lack of fit -- real series rise and fall more steeply
than the fitted trajectories, and a sharper stage window is what would let
the model chase a steeper decline. **It is wrong.**

Cells: one row per fit. `b_shape` is the posterior mean averaged over 14
`grp_init` groups. `traj` is the median-over-177-series within-series range
of `log10(y_hat + 1)` for the **noiseless trajectory**, `rep` the same for a
full posterior predictive replicate (trajectory plus observation noise), each
the median over 200 posterior draws with a 95% interval, in log10 units.
`ppp` is `P(T(y_rep) >= T(y_obs))` over those draws: ~0.5 reproduces the
statistic, near 0 means the model generates series flatter than the real
ones. **Observed `T(y_obs)` = 2.563.**

| fit | `b_shape` | traj (95%) | rep (95%) | `ppp` |
|---|---|---|---|---|
| `no_pool` | 15.0 | 1.66 (1.51-1.82) | 2.16 (2.02-2.30) | **0** |
| `np_wide_bshape` | 66.6 | 1.73 (1.53-1.90) | 2.19 (2.04-2.36) | **0** |
| `np_wide_total0` | 14.4 | 2.29 (2.06-2.44) | 2.48 (2.36-2.66) | 0.155 |
| `np_wide_both` | 65.3 | 2.31 (2.11-2.48) | 2.50 (2.34-2.64) | 0.17 |

**Widening `b_shape` does essentially nothing to the trend.** At an unchanged
`log10_total0` prior it moves the trajectory range 1.66 -> 1.73 and the
replicate range 2.16 -> 2.19, and `ppp` stays at 0. Widening `log10_total0`
moves the trajectory range 1.66 -> 2.29 and takes `ppp` from 0 to 0.155. The
two rows with the corrected `log10_total0` prior differ by a factor of 4.4 in
`b_shape` and are indistinguishable on this statistic.

So the lack-of-fit signal recorded under "Information content" was the
misspecified `log10_total0` prior flattening the trajectories, and correcting
it largely closes the gap. `ppp` of 0.155-0.17 is still on the low side --
the model runs slightly flatter than the real series -- but it is no longer a
clear misfit.

**What this does and does not settle.** It rules out the absorption
explanation for `b_shape` on this statistic: whatever the 24.8 elpd is
buying, it is not a steeper time course. It does not show that `b_shape` ~ 65
is a real measurement; that is the next section.

### `b_shape` is not identified: the design cannot separate 15 from 65

`_scripts/wockner-schedule-sim.R` with `SCHEDSIM_TRUTH_FIT=pl_wide_both`
(SLURM 28941, 3 tasks, all COMPLETED), read by
`_scripts/schedsim-truthfit-check.R`, output
`_data/schedsim-truthfit-check.log`. Simulates from `pl_wide_both`, whose
truth has `b_shape` = 65.0 and `cycle_length` = 43.7361 h, and fits under the
same widened priors the truth was generated under.

Cells: one row per replicate, all of arm `wide_nuis`, so
`sd_log_b_shape = 1.5` and `sd_log10_total0 = 1` in **every** row and the
only thing that differs between the blocks is which fit the truth came from.
`b_shape` true and est are that parameter averaged over 14 `grp_init` groups
-- the value simulated from, and the posterior mean recovered; `rel` is
`(est - true) / true`. `cl` columns are the population `cycle_length` in
hours and `cl_bias` is est minus true. One fit per row, nothing averaged.

| replicate | truth fit | `b_shape` true | est | `rel` | `cl` true | est | `cl_bias` |
|---|---|---|---|---|---|---|---|
| rep1 | `pl_wide_both` | 65.0 | 32.2 | -0.505 | 43.74 | 45.8 | +2.05 |
| rep2 | `pl_wide_both` | 65.0 | 35.6 | -0.452 | 43.74 | 45.5 | +1.74 |
| rep3 | `pl_wide_both` | 65.0 | 37.4 | -0.424 | 43.74 | 45.4 | +1.69 |
| rep1 | `pooled_cl` | 14.9 | 21.9 | +0.475 | 45.01 | 47.4 | +2.44 |
| rep2 | `pooled_cl` | 14.9 | 23.0 | +0.546 | 45.01 | 45.1 | +0.06 |
| rep3 | `pooled_cl` | 14.9 | 22.9 | +0.540 | 45.01 | 46.8 | +1.80 |
| rep4 | `pooled_cl` | 14.9 | 27.9 | +0.877 | 45.01 | 47.0 | +2.02 |
| rep5 | `pooled_cl` | 14.9 | 19.5 | +0.313 | 45.01 | 46.8 | +1.82 |

**The estimate is compressed toward ~25 whatever the truth is.** A truth of
14.9 returns 19.5-27.9; a truth of 65.0 returns 32.2-37.4. The truth changes
by a factor of 4.36 and the estimate by a factor of 1.53 -- on the log scale
the design transmits **0.287** of a change in `b_shape` onto the estimate.
The two ranges do not overlap, so the design carries *some* information, but
it never gets near either truth: it overshoots a small `b_shape` by ~50% and
undershoots a large one by ~46%.

**This runs the opposite way from the artefact explanation.** The real data
under the same widened prior report 65.4, which is *above* what this design
returns even when the truth is genuinely 65. So the real-data 65.4 is not an
estimator wandering up from a modest true value -- the attenuation pushes the
other way. **Do not inverse-calibrate it.** Reading 65.4 back through the
0.287 slope extrapolates far outside the two points that define it and gives
a number this design cannot support. The defensible statement is
directional: `b_shape` is unidentified over at least 15-65, and nothing in
these data licenses a point estimate for it.

**The cycle-length bias survives the corrected priors.** From the
`pl_wide_both` truth, fitted with both nuisance priors widened, bias is
+1.69 to +2.05 h (mean +1.83). From the `pooled_cl` truth with the same
widened priors it is +1.63 h over 5 replicates. Against the `default` arm's
+1.97 h, widening buys ~0.5 h of ~2 h and leaves ~1.5 h. **A shift in
`cycle_length` when a prior is widened is not bias removal**, on real data or
in simulation.

### Replicates 4-5: the arms gained, the paired budget did not

`_scripts/wockner-schedule-sim.sh`, SLURM 28940, 13 tasks, all COMPLETED
(exit 0), read by `_scripts/wockner-schedule-sim-analyze.R`, output
`_data/schedsim-analyze.log`.

**Three of the new replicates failed the convergence gate** and are excluded:
`default-rep4` (max R-hat 1.22, 478 divergences, min ESS 13),
`default-rep5` (1.07, ESS 57), and `wide_total0-rep4` (1.06, ESS 63). Per
thread 9 both earlier failures of this kind were bad **fit seeds** rather
than hard datasets, so these are candidates for a refit on a new
`SCHEDSIM_FIT_SEED`, not evidence about the posterior.

Cells: `n_arm` is how many replicates of that arm converged; `n_pair` how
many have a converged `default` replicate on the **same simulated dataset**,
which is what the paired contrast can use. `mean_delta` is the paired change
in cycle-length bias in hours, that arm minus `default`, **negative means
widening reduced the bias**; `min`/`max` are over the paired replicates.

| arm | `n_arm` | `n_pair` | `mean_delta` | min | max |
|---|---|---|---|---|---|
| `wide_bshape` | 5 | 3 | -0.501 | -0.531 | -0.452 |
| `wide_nuis` | 5 | 3 | -0.540 | -1.08 | -0.142 |
| `wide_total0` | 4 | 3 | +0.129 | -0.156 | +0.429 |

**`n_pair` is still 3.** Both new `default` replicates failed, and `default`
is the baseline every paired contrast subtracts, so the batch added no paired
replicates at all and every number in the budget is unchanged. The unpaired
arms did gain: `no_hier` is now n=5 (+2.56, +0.928, +1.97, +2.17, +1.57) and
`wide`, `wide_bshape` and `wide_nuis` are n=5.

Nuisance recovery moved very little with the extra replicates, which is worth
having: `wide_bshape`'s `b_shape` goes 23.9 (n=3) to 24.3 (n=5),
`wide_nuis`'s `log10_total0` +25.6% to +20.0%. The recovery story is not an
artefact of three noise draws.

**Until `default-rep4` and `default-rep5` are refit on a new seed, option 3
has not delivered what it was run for.**

### What `b_shape` means biologically, and where the data sit

`_scripts/bshape-upper-range.R`, output `_data/bshape-upper-range.log`.
Post-hoc on saved fits plus a deterministic mapping; no refitting.

`b_shape` is the shape of a **symmetric Beta(`b_shape`, `b_shape`)** over
position in the cycle (`inst/stan/functions/plasmofit.stan`, `beta_starts`),
discretised into `n_c` = 96 stage cells and declared
`<lower=2, upper=max_shape>` with `max_shape` = 250. So it is a synchrony
parameter, and it converts into an age range.

Cells: `age_range_h` is the width of the central 99% of the starting-stage
distribution -- the 0.5% and 99.5% quantiles of Beta(`b_shape`, `b_shape`)
scaled by a 45.3 h cycle, in hours. Exact Beta quantiles, not a normal
approximation. Deterministic; no data enter.

| `b_shape` | age range (h) | where it comes from |
|---|---|---|
| 14.9 | **20.4** | posterior mean under the default `lognormal(2, 0.5)` prior |
| 48.3 | 11.7 | posterior **median**, widened prior |
| 65.4 | 10.1 | posterior **mean**, widened prior |
| 84 | 8.93 | the 9 h initial age range assumed in `mmcm.pdf`, Figure 1 |
| 100 | 8.19 | -- |
| 250 | 5.2 | the `max_shape` bound |

**The default prior implies a starting-stage spread of 20.4 h on a ~45 h
cycle.** That is close to no synchrony at all, which is not what an
inoculation is. The misspecification of the `b_shape` prior is therefore
visible on biological grounds and not only through elpd, and the direction
the data move when it is relaxed -- toward 10-12 h -- is *toward* the
literature's assumption, not away from it.

Note the anchor's status: the 9 h figure is a **modelling assumption** in
`mmcm.pdf`'s Figure 1 for controlled human infection trials, not a
measurement of any inoculum. It fixes an order of magnitude, nothing finer.

**The widened-prior posterior is strongly right-skewed, and the panel
reported its mean.** Cells: quantiles and tail probabilities of `b_shape`
pooled over 14 `grp_init` groups and all post-warmup draws of one real-data
fit. Dimensionless.

| fit | mean | median | q90 | q97.5 | P(>100) | P(>240) |
|---|---|---|---|---|---|---|
| `no_pool` | 14.9 | 13.6 | 23.9 | 32.2 | 0 | 0 |
| `np_wide_bshape` | 65.4 | **48.3** | 144 | 210 | 0.208 | 0.005 |
| `np_wide_both` | 64.4 | 47.5 | 143 | 208 | 0.203 | 0.004 |

Two consequences. **The bound is not binding** -- 0.4-0.5% of the mass lies
within 4% of `max_shape` = 250, so the 65.4 is not an artefact of where the
model was cut off. And **65.4 overstates the centre**: the median is 48.3, so
"the data say 65" is a mean of a long-tailed posterior. Quote the median, or
the age range, not the mean.

**The effect of `b_shape` collapses at high values.** The starting-stage sd
is `cycle_length / (2 * sqrt(2 * b_shape + 1))`, so a further +10 in
`b_shape` is worth 0.90 h at `b_shape` = 15, 0.136 h at 65, 0.074 h at 100,
and 0.020 h at 250. Anything above ~100 is nearly unobservable at this
design, which is the mechanism behind the attenuation measured in "the design
cannot separate 15 from 65".

### The `b_shape` ladder: fixing synchrony is a large predictive gain

`_scripts/wockner-bshape-ladder.R`, output `_data/bshape-ladder.log` and
`_data/wock-bshape-ladder.rds`. Fits are `wockner-fit.R` tasks 14-18, SLURM
28963, all five COMPLETED. `b_shape` pinned by a tight prior
(`sd_log_b_shape` = 0.05, ±10% at 95%) at five centres, all carrying the
corrected `log10_total0` prior, so the ladder reads against `np_wide_total0`
(same prior, `b_shape` free) and `np_wide_both` (same prior, `b_shape` prior
widened). `max_shape` is 400 for the five rungs so the 250 rung is not on its
own boundary; it is 250 for the two references, and non-binding in all seven.

Cells: `pinned` is where the tight prior was centred, blank where `b_shape`
was free. `b_shape` and `cycle_length` are posterior means averaged over 14
`grp_init` and 13 `grp_cl` groups; `age_range_h` is the central 99% of the
starting-stage distribution implied by the fitted `b_shape` at a 45.5 h
cycle, in hours; `d_cl` is `cycle_length` minus `np_wide_total0`'s, in hours.
`div_pct` is divergent transitions as a percent of 4000 post-warmup
transitions. `vs_free` is the paired PSIS-LOO difference against
`np_wide_total0` and `vs_bs100` against the `np_bs100` rung, summed over the
same 1130 observations, **positive means better**, with the standard error of
the paired difference.

| fit | pinned | `b_shape` | age range | `cycle_length` | `d_cl` | `div_pct` | max R-hat | `vs_free` (se) | `vs_bs100` (se) |
|---|---|---|---|---|---|---|---|---|---|
| `np_wide_total0` | free | 14.4 | 20.8 | 45.5 | 0 | 2.05 | 1.02 | 0 | −34.0 (2.91) |
| `np_bs50` | 50 | 50.2 | 11.5 | 44.0 | −1.53 | 2.33 | 1.03 | **+28.5** (2.46) | −5.44 (0.58) |
| `np_bs84` | 84 | 84.2 | 8.96 | 43.9 | −1.66 | 3.35 | 1.03 | **+32.6** (2.85) | −1.34 (0.32) |
| `np_bs100` | 100 | 100 | 8.22 | 43.8 | −1.74 | 3.58 | 1.02 | **+34.0** (2.91) | 0 |
| `np_bs150` | 150 | 150 | 6.73 | 43.8 | −1.77 | 3.48 | **1.05** | **+34.8** (3.12) | +0.84 (0.51) |
| `np_bs250` | 250 | 250 | 5.22 | 43.6 | −1.92 | 4.58 | 1.04 | **+36.2** (3.20) | +2.26 (0.50) |
| `np_wide_both` | free (wide) | 64.4 | 10.2 | 44.2 | −1.29 | 3.42 | 1.02 | +24.4 (2.06) | −9.51 (1.06) |

**Fixing `b_shape` is worth +28 to +36 elpd over leaving it free, at every
rung tested**, and every rung also beats `np_wide_both`, which merely widens
its prior. That is a large gain, 10-12 standard errors, and it is the main
result: the argument for fixing this parameter does not rest on the biology
alone.

**Where it is fixed matters far less, but the data are not indifferent.**
Across 84-250 the spread is 3.6 elpd, against 28-36 for fixing at all.
Within that, higher is consistently better and detectably so: `np_bs250`
beats `np_bs100` by 2.26 (se 0.50) and `np_bs100` beats `np_bs84` by 1.34
(se 0.32). `np_bs50` is clearly the worst rung, −5.44 (se 0.58).

**The ladder is still climbing at its top rung.** elpd rises monotonically to
250 with no plateau, so the data have not been shown to prefer any value --
only "higher than 50". 250 is the top of what was run, not an optimum, and
`max_shape` would have to be raised again to find one.

**`cycle_length` does not plateau either, and this is the cost.** It falls
monotonically 44.0 -> 43.6 h across the ladder, and sits 1.53-1.92 h below
the free-`b_shape` fit. So a reported cycle length is conditional on where
`b_shape` is pinned, to about 0.4 h across plausible rungs and about 1.7 h
against leaving it free. **This is the sensitivity, and it is the thing to
report** -- not the value at any one rung.

Note what does *not* move: `R` (14.2-14.3), `log10_total0` (−0.857 to
−0.867), and `sd_iRBC` (0.528-0.532) are flat across the whole ladder. Only
`cycle_length` tracks the pin.

**Pinning samples worse, not better.** Divergences rise monotonically with
the pin, 2.05% free to 4.58% at 250, and the adapted step size falls from
0.0196 to 0.0115-0.0162; `np_bs150` reaches max R-hat 1.05, at the gate. A
narrow starting distribution makes the trajectory sharper and the geometry
stiffer. **The best-fitting rung is also the worst-sampling one.**

**Two predictions recorded before the run were wrong**, and the record is
kept because it bounds how far the heuristic reasoning can be trusted. The
first was stated against the wrong reference: `np_wide_total0` leaves
`b_shape` free at 14.4, so moving to a pinned 100 is a change of 20.8 -> 8.2
h in the age range and could never have been a 0.1-0.3 h move. Against the
right reference, `np_wide_both` at `b_shape` 64.4, pinning at 100 moves
`cycle_length` −0.45 h, still above the predicted 0.1-0.3 h. The second was
that pinning a bounded scalar would sample no worse than widening its prior;
it samples worse at every rung. **The argument that `b_shape`'s effect
collapses at high values is true on the `b_shape` scale and misleading on the
age-range scale**, which is what `cycle_length` tracks: per hour of age-range
narrowing the effect on `cycle_length` does not fall off at all, and between
the 150 and 250 rungs it is larger than between 50 and 84.

### Refits: the budget reaches n=4, and `default-rep4` is a near-miss

SLURM 28964, 3 tasks, all COMPLETED, read by
`_scripts/wockner-schedule-sim-analyze.R`, output `_data/schedsim-analyze.log`.
`SCHEDSIM_FIT_SEED` = 1618033989 and nothing else changed.

- `default-rep5` converged (max R-hat 1.02, 46 divergences, min ESS 144).
- `wide_total0-rep4` converged (1.01, 70, 203).
- **`default-rep4` did not** -- max R-hat **1.06** against a 1.05 gate. But
  it is a near-miss, not a hard dataset: divergences fell 478 -> 39 and min
  ESS rose 13 -> 117 against the original seed. The seed did most of the
  work and the gate is close; a third seed or a longer warmup would likely
  clear it. This is weaker than thread 9's two cases, where a new seed
  settled it outright, and stronger than "the dataset is hard".

**The paired budget is now n=4, and the headline number moved.** Cells: as in
"Thread 2: it is the `b_shape` prior" -- paired change in cycle-length bias
in hours, arm minus `default` on the same simulated dataset, negative means
widening reduced the bias; `n_pair` is replicates where both converged.

| arm | `n_pair` | mean change | range | was (n=3) |
|---|---|---|---|---|
| `wide_bshape` | 4 | **−0.433 h** | −0.531 to −0.228 | −0.501 h |
| `wide_nuis` | 4 | −0.423 h | −1.08 to −0.071 | −0.540 h |
| `wide_total0` | 4 | +0.099 h | −0.156 to +0.429 | +0.129 h |

**The `b_shape` prior's share of the simulated cycle-length bias is ~0.43 h,
not ~0.50 h, and the spread is four times wider than three replicates
suggested** -- −0.228 to −0.531 against −0.452 to −0.531. The old figure's
apparent precision (0.08 h across three replicates) was luck. This is what
option 3 was run for and it is the one number in the budget that the extra
replicates actually changed.

Nuisance recovery barely moved at n=4-5: `default` `b_shape` 8.89 -> 8.77,
`log10_total0` +102% -> +101%, `R` −23.8% -> −23.6%.

### The ladder plateaus at ~400, and bound geometry is not the bias

`_scripts/wockner-bshape-ladder.R` (output `_data/bshape-ladder.log`,
`_data/wock-bshape-ladder.rds`) and
`_scripts/wockner-schedule-sim-analyze.R` (output
`_data/schedsim-analyze.log`). SLURM 29493 (3 real-data fits) and 29494 (6
simulation fits), all COMPLETED.

**Correction to the previous ladder section.** `np_bs150`'s max R-hat is
**1.0513**, which fails the R-hat < 1.05 gate; it was described there as "at
the gate" and should have been dropped. Nothing rests on it -- elpd is
monotone through it either way -- but it is not a valid rung.

#### Where the ladder stops climbing

Cells: `age_range_h` is the central 99% of the starting-stage distribution at
a 45.5 h cycle, in hours. `d_cl` is `cycle_length` minus `np_wide_total0`'s,
in hours. `div_pct` is divergent transitions as a percent of 4000 post-warmup
transitions. `vs_bs100` is the paired PSIS-LOO difference against the
`np_bs100` rung over the same 1130 observations, **positive is better**, with
the standard error of the paired difference.

| rung | age range | `d_cl` | `div_pct` | max R-hat | `vs_bs100` (se) |
|---|---|---|---|---|---|
| 100 | 8.22 | −1.74 | 3.58 | 1.021 | 0 |
| 150 | 6.73 | −1.77 | 3.48 | **1.051 (fails)** | +0.84 (0.51) |
| 250 | 5.22 | −1.92 | 4.58 | 1.043 | +2.26 (0.50) |
| 250, `max_shape` 1000 | 5.22 | −1.83 | 4.68 | **1.055 (fails)** | +2.26 (0.42) |
| 400 | 4.14 | −1.89 | 4.32 | 1.029 | **+3.22** (0.54) |
| 600 | 3.38 | −1.95 | 5.92 | 1.024 | +3.11 (0.51) |

Direct paired comparisons at the top of the ladder, same cells:

| contrast | elpd_diff (se) | z |
|---|---|---|
| 400 vs 250 | **+0.96** (0.29) | +3.3 |
| 600 vs 400 | −0.11 (0.35) | −0.3 |
| 250 `max_shape` 1000 vs 250 | **+0.00** (0.38) | 0.0 |

**elpd plateaus at 400.** It is still detectably climbing from 250 to 400
(z 3.3) and flat from 400 to 600 (z −0.3). So the data prefer a
starting-stage age range of about **4 h**, and asking for narrower buys
nothing. That is the answer the extension was run for.

**`max_shape` on its own is worth nothing predictively** -- +0.00 elpd
(se 0.38) between two fits pinned identically at 250 -- so the extended rungs
are comparable with 14-18 despite using 1000 against 400. Note the control
fit itself fails the R-hat gate at 1.0550, so by this project's own rule it
is void. The conclusion does not depend on it: `np_bs400` beats `np_bs250`
by +0.96 (se 0.29) and `np_bs250_ms1000` by +0.96 (se 0.31), identically, so
which of the two is used as reference does not matter.

**`cycle_length` is flat across the top of the ladder** -- 43.6 to 43.7 h
from 250 to 600 -- against 44.0 h at 50 and 45.5 h free. So the sensitivity
reported earlier does not widen: anywhere from 250 upward gives the same
answer to within 0.1 h.

**Pushing past the plateau costs sampling.** Divergences are worst at 600
(5.92%, against 2.05% with `b_shape` free). There is no reason to pin above
~400 and a reason not to.

**Against the biology**, 4 h is *more* synchronous than the 9 h initial age
range assumed for controlled human infection trials in `mmcm.pdf` Fig. 1
(`b_shape` ~ 84). The data want tighter synchrony than that assumption, not
looser, and the 84 rung is detectably worse than 100, 250, and 400.

#### Bound geometry: not the explanation

Two arms, both **moving** the `[35, 50]` window rather than only widening it.
`cl_move` is `[37.5, 52.5]` -- the same 15 h width centred on the 45.012 h
truth, isolating asymmetry -- and `cl_wide_move` is `[30, 60]`, moved and
widened, asking whether distance from any bound matters.

Cells: paired change in cycle-length bias in hours, arm minus `default` on
the same simulated dataset, **negative means the change reduced the bias**;
`n_pair` is replicates where both converged.

| arm | `n_pair` | mean change | range |
|---|---|---|---|
| `cl_move` | 2 | **+0.059 h** | −0.04 to +0.157 |
| `cl_wide_move` | 1 | **+0.03 h** | -- |
| `wide_bshape` (for scale) | 4 | −0.433 h | −0.531 to −0.228 |

**Bound geometry does not carry the bias.** To account for the ~1.5 h that
survives correcting the nuisance priors, moving the bounds would have to move
it by about −1.5 h. It moves it by +0.06 h, and the sign is *wrong* -- a
centred window is marginally worse, not better. Even the most extreme single
replicate (+0.157 h) is an order of magnitude too small and the wrong way.
Neither asymmetry (`cl_move`) nor distance from the bounds (`cl_wide_move`)
is the mechanism. **Thread 2 is out of cheap suspects.**

**The evidence is thinner than planned: 3 of 6 replicates failed the gate.**
`cl_move-rep2` failed badly (max R-hat 1.32, 513 divergences, min ESS 10),
and `cl_wide_move-rep2` and `-rep3` are near-misses (1.07 and 1.06). So the
contrasts are n=2 and n=1. That is weak for a precise estimate and adequate
for the conclusion drawn, which is only that the effect is nowhere near the
size needed -- the paired contrast's own spread within an arm is ~0.2-0.3 h
across every arm measured, so a −1.5 h effect could not hide in it.

**Moving the bounds degrades sampling**, which is itself consistent with the
package's note that tight bounds are worth real compute and with `max_cl` =
55 having reintroduced a boundary mode. Any future bound experiment should
budget for refits on new seeds.

### The 14 `grp_init` groups do not differ in `b_shape`

`_scripts/bshape-between-group.R`, output `_data/bshape-between-group.log`.
Post-hoc on `_data/wock-prior-panel.rds`, no refitting.

The ladder pinned all 14 groups at a common centre, so it tested **fixed vs
estimated** and said nothing about **one shared value vs 14**. Collapsing
`vector[n_grp_init] b_shape` to a scalar is a separate change and needs its
own evidence.

Cells: over the 14 `grp_init` groups of one fit, `sd_between` is the sd of
the per-group posterior **means** and `sd_within` the mean of the per-group
posterior **sds**; `ratio` is their quotient, so well below 1 means the
groups sit closer together than one group's own uncertainty. `log` rows are
the same on `log(b_shape)`, the scale its prior is written on, with the
within-group sd on that scale taken as `sd/mean` -- a delta-method
approximation, adequate for an order of magnitude and not for an interval.
Dimensionless.

| fit | scale | `sd_between` | `sd_within` | `ratio` | range of means |
|---|---|---|---|---|---|
| `no_pool` | natural | 2.68 | 6.46 | **0.42** | 11.0-19.2 |
| `no_pool` | log | 0.182 | 0.438 | **0.42** | -- |
| `np_wide_bshape` | natural | 14.3 | 50.7 | **0.28** | 44.9-87.8 |
| `np_wide_bshape` | log | 0.223 | 0.799 | **0.28** | -- |
| `np_wide_both` | natural | 15.4 | 49.8 | **0.31** | 36.7-89.7 |
| `np_wide_both` | log | 0.249 | 0.800 | **0.31** | -- |

**The data carry no evidence that the groups differ.** The ratio is 0.28-0.42
in every fit on both scales: the 14 estimates are three to four times closer
together than one group's own posterior sd. The per-group means do span a
factor of ~2 under the widened prior (44.9 to 87.8), but each carries a
posterior sd near 50, so that spread sits well inside noise.

The comparison is conservative in the direction that matters. `b_shape` has
an independent prior per group and no hierarchical pooling, so nothing shrinks
the per-group estimates toward each other; if anything they should be *more*
dispersed than the truth, and they are less dispersed than the within-group
noise.

**So a single shared `b_shape` loses nothing measurable**, which is what
would license moving it to the data block as one user-set value.

### Regression test for `b_shape`-as-data: no change detected, not yet a pass

`_scripts/wockner-bshape-regression.R`, output `_data/bshape-regression.log`.
Fits are SLURM 29526 (tasks 9 and 19 refitted under the rebuilt package, as
`<config>-rebuild`) and 29528 (task 22, `np_bs400_data`), all COMPLETED.
Package branch `fixed-b-shape`, commits `0bf12b0` and `e637d8d`.

The change renames the free parameter to `b_shape_free`, makes it zero-sized
when `b_shape` is data, and rebuilds `b_shape` as a transformed parameter.
With `b_shape = NULL` that is mathematically a rename, so the old path must
be unchanged.

**Structural check, exact.** In all three comparisons the only name present
in the new fit and absent from the old is `b_shape_free`, and nothing is
present in the old and absent from the new. That is the expected signature of
the rename and nothing else moved.

Cells: over the scalar entries both fits share, excluding `lp__`, `log_lik`,
and any entry with zero or missing `se_mean` (a constant has no Monte Carlo
error and no meaningful z). `z` is the difference in posterior means over
`sqrt(se_mean_new^2 + se_mean_old^2)`, so the shift in units of the two runs'
combined Monte Carlo standard error. `n` is entries compared, `frac > 2` the
share beyond 2. Dimensionless.

| | comparison | `n` | max abs `z` | worst entry | `frac > 2` | `frac > 3` |
|---|---|---|---|---|---|---|
| A | `np_wide_total0-rebuild` vs `np_wide_total0` | 206 | **3.52** | `b_shape[2]` | 0.175 | 0.005 |
| B | `np_bs400-rebuild` vs `np_bs400` | 206 | 2.56 | `z_sd_iRBC[21]` | 0.063 | 0 |
| C | `np_bs400_data` vs `np_bs400` | 192 | 2.02 | `z_sd_iRBC[26]` | 0.005 | 0 |

C compares 14 fewer entries because `b_shape` is constant there and drops out
on the zero-`se_mean` rule.

All five fits clear the R-hat gate: 1.0196 and 1.0201 (A), 1.0291 and 1.0242
(B), 1.0468 for `np_bs400_data`, whose worst entry is `lp__`.

**B and C are unremarkable.** No entry moves beyond 3 combined MCSE, and C's
distribution is tighter than either regression pair -- so passing `b_shape`
as data lands where pinning it with a tight prior does, which is what the
`sd_log_b_shape` = 0.05 prior being nearly a point mass predicts.

**A is the one to look at, and this test cannot resolve it.** Its worst entry
is `b_shape[2]` at 3.52, and 17.5% of entries exceed 2 against about 5%
expected if `z` were standard normal. Two explanations are confounded here
and neither is favoured by the data in hand:

- `b_shape` under a free, wide prior is the most heavy-tailed quantity in
  this model -- posterior median 48 against a mean of 65, sd 50 (see
  "What `b_shape` means biologically"). A posterior *mean* is a poor,
  high-variance summary of it, and `se_mean` understates run-to-run spread
  worst exactly there. On that reading A's inflation is an artefact of the
  statistic.
- `b_shape` is also the parameter the code change touches, so it is where a
  real regression would show first.

**The comparison is anti-conservative by construction**, which is why the
inflation cannot be waved away: two runs differ in step-size and mass-matrix
adaptation as well as sampling noise, and MCSE captures none of that. Against
a correctly-scaled null, |z| of 3.5 on the least stable statistic in the
model may be ordinary.

**So: no change detected, and that is not the same as a pass.** By the
standard `wockner-anchor-regression.R` sets, a regression test without a null
run is inconclusive, and none exists for these configs. The null is two fits
of identical code and data differing only in the sampler seed;
`wockner-fit.R` now takes `WOCKFIT_SEED` so one can be produced:

```
sbatch --array=9 --export=ALL,WOCKFIT_SEED=1618033989,WOCKFIT_SUFFIX=-null \
    _scripts/wockner-fit.sh
```

Read A's 3.52 against that run's max |z| before concluding anything. If the
null reaches a similar figure, A is Monte Carlo variation; if it sits near 2,
A is a real shift in `b_shape` and the change is not inert.

### The null run: the change is clear, and something else is not

`_scripts/wockner-bshape-regression.R` (output `_data/bshape-regression.log`)
and `_scripts/bshape-regression-blocks.R` (output
`_data/bshape-regression-blocks.log`). The null is SLURM 29531,
`np_wide_total0-null`: identical code and data to `np_wide_total0-rebuild`,
sampler seed alone differs (`WOCKFIT_SEED` = 1618033989).

**The null is well calibrated**, which is what makes the rest readable. Its
mean |z| is **0.85** against the 0.798 expected for |N(0, 1)|, so combined
Monte Carlo standard error is an accurate yardstick for run-to-run variation
in this model. Max |z| 2.41, 4.1% of entries beyond 2.

Cells: mean over the scalar entries of |difference in posterior means| over
the two runs' combined Monte Carlo standard error. `ratio` is each
comparison's mean |z| over the null's, so 1 means a comparison adds nothing
beyond what the sampler seed alone produces. Dimensionless.

| | comparison | `n` | mean \|z\| | ratio to null | max \|z\| |
|---|---|---|---|---|---|
| N | **null**, seed alone | 220 | 0.85 | 1 | 2.41 |
| A | `np_wide_total0-rebuild` vs `np_wide_total0` | 206 | 1.23 | **1.45** | 3.52 |
| B | `np_bs400-rebuild` vs `np_bs400` | 206 | 0.92 | 1.08 | 2.56 |
| C | `np_bs400_data` vs `np_bs400` | 192 | 0.81 | 0.95 | 2.02 |

**Read mean |z|, not max |z|.** The max is one order statistic over ~200
entries and is far noisier. A's max of 3.52 looked alarming on its own; the
block split below shows the entry it lands on is not where the change is.

**By block**, test over null, for A (`b_shape` free) and B (`b_shape` pinned):

| block | `n` | null mean \|z\| | ratio A | ratio B |
|---|---|---|---|---|
| `mu_logit_R` | 1 | 0.72 | **3.49** | **2.48** |
| `R` / `eta_R` / `logit_R` | 13 each | 0.72 | **2.81** | 1.61–1.63 |
| `mu_logit_cl` | 1 | 0.90 | 2.57 | 0.35 |
| `cycle_length` / `eta_cl` / `logit_cl` | 13 each | 0.67–0.72 | 1.80–1.87 | 1.15–1.21 |
| `log10_total0` | 14 | 0.82 | 1.79 | 1.00 |
| `b_offset` / `b_off_vec` | 14 / 28 | 0.85–1.06 | 1.05–1.08 | 0.68–0.86 |
| `sd_iRBC` / `z_sd_iRBC` | 27 each | 0.89–0.90 | 0.91–0.92 | 1.08–1.09 |
| **`b_shape`** | 14 | 1.00 | **0.90** | **1.12** |
| `sigma_logit_R` / `sigma_logit_cl` | 1 each | 1.31–1.39 | 0.25–0.45 | 0.19–0.47 |

**`b_shape` is the quietest block in the table.** Ratio 0.90 and 1.12 — at
the null in both test pairs. The parameter the code change touches is the one
that does not move. B and C pass outright, and C, the new data path, is the
closest of all three to the null.

**So the `b_shape` change is clear**, as firmly as this design allows: if it
had perturbed anything, `b_shape` is where that would appear first, and
`b_shape` is flat.

**But A's excess is real and is a separate finding.** It sits in the `R`
hyperparameter and the `(log10_total0, R)` ridge — `mu_logit_R` at 3.5x the
null, `R`/`eta_R`/`logit_R` at 2.8x — and it is present in **B as well**,
where `b_shape` is pinned, so `b_shape`'s freedom does not cause it. What A
and B share and the null does not is the **build**: both compare a fit from
the 2026-09 package against one from the 2026-10 package, where the null
compares two fits from the same build.

**What that means for reported numbers.** Build-to-build variation in the
weakly identified directions is **larger than Monte Carlo error**, by roughly
2.5–3.5x in the `R` block. Any `R` or `mu_logit_R` figure therefore carries
irreproducibility across package rebuilds well beyond its stated MCSE. This
is consistent with everything else known about those parameters — `R` is not
identified by these data and sits on a ridge with `log10_total0` — but it has
not been measured before and it is not captured by any interval this project
reports.

**Not separated**: the two builds differ by the `b_shape` change *and* by a
recompile, and the change alters the `parameters` block, so a given seed no
longer yields the same initial values across builds. Isolating pure rebuild
drift needs the pre-change code rebuilt and refitted, which is a ~2.5 h fit
plus a reinstall that would clobber the current one. The block split is what
makes that unnecessary for clearing the change; it is still what would be
needed to characterise the drift itself.

### Thread 1's anchor regression: PASS, on a null that was already on disk

`_scripts/wockner-anchor-regression.R`, output `_data/anchor-regression.log`.
No new fitting: it reads three saved fits.

The thread has sat at "one fit short of a verdict" since 2026-09-24 because
SLURM 28892 was recorded as FAILED. It was — **at the summary stage, after
writing the fit.** `wock-schedsim-fit-default-rep2-seed271828183.rds` has
been on disk the whole time, and the script's own glob finds it.

Cells: `test` is the pre-change fit against the post-change rebuild, same
seed; `null` is two runs of identical post-change code differing only in
sampler seed. `median |log sd ratio|` is over the 208 shared scalar entries,
measuring posterior-width agreement; `|z|` is the difference in posterior
means over its Monte Carlo standard error. Dimensionless.

| | median \|log sd ratio\| | 95% | \|z\| median | \|z\| max | \|z\| > 3 |
|---|---|---|---|---|---|
| test (across the change) | **0.0262** | 0.1248 | **0.73** | 2.67 | 0 of 208 |
| null (seed alone) | 0.0368 | 0.1241 | 1.12 | 2.93 | 0 of 208 |

**The test is tighter than the null on both statistics.** Posteriors agree
better across the anchor change than two runs of identical code differ from
each other, and width inflation at the 95% is 1.01x. Nothing in 208 entries
moves beyond 3 MCSE in either. With the anchor off, the change targets the
same posterior as before it.

`log10_total0`, the parameter the anchor exists to re-prior, is the one to
check specifically: max |z| 1.79, sd ratio 0.964-1.049.

**Thread 1's first half is closed.** What remains of that thread is the
anchored fit itself, which was taken on 2026-09-24 (`np_anchor`) and is
written up under "The `log10_total0` prior was misspecified".

### The posterior median does not help the cycle-length bias

`_scripts/summary-stat-check.R`, output `_data/summary-stat-check.log`.
Post-hoc on four saved fits, no refitting.

Cells: per-trial posterior `cycle_length` averaged over 13 trials, by summary
statistic, in hours; `bias` is that minus the simulated truth of 45.012 h,
**positive means overestimation**; `shift` is median minus mean; `skew` is
the median over trials of each trial's posterior skew, 0 being symmetric.

| fit | bias by mean | bias by median | shift | skew |
|---|---|---|---|---|
| `default-rep1` | +2.58 | **+2.66** | +0.080 | −0.52 |
| `default-rep2` | +1.14 | **+1.17** | +0.031 | −0.18 |
| `default-rep3` | +2.19 | **+2.28** | +0.088 | −0.58 |

**It makes the bias slightly worse, not better**, by 0.03-0.09 h, because the
per-trial `cycle_length` posteriors are **left**-skewed, so the median sits
*above* the mean. That matches the earlier mode check from the other side:
+1.82 h by the mean against +1.87 h by the mode. Mean, median, and mode span
0.05 h on a bias of 1.8 h — a 3% effect on a problem that is 100%.

**The deeper reason is not skew.** Mean, median, and marginal mode are all
**marginal** summaries: each integrates over the other ~200 parameters alike.
If the bias comes from integrating over the `(log10_total0, R)` ridge and two
unidentified nuisances, all three carry it equally and choosing between them
cannot help. Only a **joint** summary -- the MAP -- escapes marginalisation,
which is why `rstan::optimizing` on the same model and priors is the
discriminating test and a change of summary statistic is not.

**Where the median does change a reported number.** Cells: over each
parameter's groups, the mean of the per-group posterior means and medians;
`rel_shift` is (median − mean) / mean; `skew` the median over groups of the
per-group posterior skew. Dimensionless.

| parameter | skew | `rel_shift` |
|---|---|---|
| `b_shape` | **+1.45 to +1.77** | **−8.5% to −10.7%** |
| `sd_iRBC` | +0.54 to +0.60 | −1.1% |
| `R` | +0.12 to +0.55 | −0.4% to −0.7% |
| `log10_total0` | −0.25 to +0.15 | −0.2% to −0.6% |
| `cycle_length` | +0.06 (real), −0.18 to −0.58 (sim) | ~0 |

**`b_shape` is the one parameter where the summary choice matters**, and it
matters by about 10% under the default prior as well as the 65.4-vs-48.3 gap
already recorded under the widened one. Quote its median. Everything else
moves by around 1% or less.

**An observation not previously noted**: the real `cycle_length` posterior is
symmetric (skew +0.06) while the simulated ones are distinctly left-skewed
(−0.18 to −0.58). The simulated estimates sit at 46.2-47.6 h, within 2.4-3.8 h
of `max_cl` = 50, where the real fit at 45.5 h is further from it — so this
looks like proximity to the upper bound shaping the posterior. It is **not**
the bias mechanism: moving the bounds to `[30, 60]`, which puts the estimate
13 h clear of them, changed the paired bias by +0.03 h.

### The MAP check: the joint surface is too rough to have a mode worth reporting

`_scripts/schedsim-map-check.R`, outputs `_data/map-check-rep{1,2,3}.log` and
`_data/wock-map-check-*.rds`. `rstan::optimizing` on the same model, priors,
and simulated data the fits used, from many random starts.

The question was whether the cycle-length bias is a **marginalisation**
effect. Mean, median, and marginal mode cannot distinguish it, because all
three integrate over the other ~200 parameters alike. A joint maximum does
not integrate, so it sits at the other end of a three-point ladder from the
pooled MLE and the posterior mean.

**The answer is that there is no usable joint maximum here.** Cells: `best
lp` is the largest joint log posterior found across random starts; `bias` is
that optimum's `cycle_length` averaged over 13 trials, minus the simulated
truth of 45.012 h; `lp range` and `cl range` are across converged starts;
`found twice` counts starts reaching the best mode within 0.01 lp.

| replicate | starts | best lp | bias | lp range | cl range | found twice |
|---|---|---|---|---|---|---|
| `default-rep1` | 198 | −818.65 | +0.62 | 2097 | 8.2 h | 1 of 198 |
| `default-rep2` | 198 | −854.37 | +2.47 | 2014 | 7.3 h | 1 of 198 |
| `default-rep3` | 199 | −878.23 | **−1.67** | 2009 | 12.6 h | 1 of 199 |
| `default-rep2`, rerun | 120 | **−823.99** | +0.31 | — | — | 1 of 120 |

**The best mode was never found twice, and more searching keeps finding
better ones** -- the 120-start rerun of rep2 beat the 198-start run by 30 lp.
Modes span ~2000 nats and 7-13 h of `cycle_length`. There is no evidence any
of these is the global optimum, so **"report the MAP instead" is not a
remedy**, which was the pre-registered outcome for starts disagreeing.

**What the modes look like.** From the rep2 rerun, top 15 by lp. Two families:
some put all 13 trials at one `cycle_length` (min = max, the hierarchy
collapsed, `sigma_logit_cl` -> 0) and some spread them over 10+ h
(37.6-49.9). `b_shape` ranges **4.3 to 28.2** across modes -- different modes
explain the same data with entirely different synchrony, which is what
non-identification looks like on the joint surface rather than in a marginal.
`norm_boff` is 1 at every optimum, so the `unit_vector` is properly
normalised and radial degeneracy is **not** the cause.

**Suggestive, not conclusive.** The better modes sit nearer the truth than
the posterior mean does -- rep2's best is +0.31 h against a posterior mean of
+1.22 h -- which is the direction the marginalisation hypothesis predicts.
But with no optimum found twice, this cannot carry weight.

**`rstan::optimizing` is not reproducible here at a fixed seed.** Three calls
with identical data, identical `seed`, and `set.seed()` fixed returned
lp −935.5, −955.1, −961.3. It happens at one thread as well as four, so it is
not `reduce_sum`'s summation order. The clean test -- a fixed numeric init --
**segfaults**, because `init = 0` puts the `unit_vector` at the origin, so
the cause is not isolated.

**This does not extend to sampling.** `tests/testthat/test-holdout.R` asserts
that two `archer_fit` runs with the same seed and data return *bit-identical*
draws, and it passes. Whatever this is, it is specific to optimisation, and
no claim about fit reproducibility follows from it.

**What it leaves.** The marginalisation question is not settled; what is
settled is that the joint mode cannot settle it. The remaining separable
component of the MLE-to-posterior gap is the one the notes already name and
have never isolated: **the MLE fixes `sd_iRBC` at truth where the fit
estimates it.** That is a single simulation arm, and it is now the cheapest
untried thing on this question.

### Design A: the hierarchy on `cycle_length` earns its keep

`_scripts/wockner-designA-score.R`, output `_data/designA-score.log`. Fits are
`wockner-fit.R` configs 23-26, SLURM 29538, all COMPLETED. The last third of
every series is masked -- 306 of 1130 observations, every series keeping 3-6
-- all four model variants are fitted on the remainder, and each is scored on
the window none of them saw.

**No importance sampling is involved and none of its diagnostics apply.**
These points were genuinely not fitted, which is the whole reason the design
exists: observation-level `loo` is the weak comparison by
`archer_log_lik()`'s own docs, and trial-level `loo` is broken here, Pareto
k > 0.7 for all 13 units in every model.

Cells: `held_out` is the log pointwise predictive density summed over the 306
observations none of these models were fitted on, computed as
`log(mean_s exp(log_lik[s, i]))` per observation; `in_sample` the same over
the 824 that were fitted. **Higher is better**, both in log units. `vs_ref`
is the **paired** held-out difference against `no_pool` over the same points,
with the standard error of that paired difference; positive favours the row.

| model | `held_out` | `in_sample` | `vs_ref` (se) | z | max R-hat |
|---|---|---|---|---|---|
| `no_pool` (hierarchy kept) | **−332.3** | −604.5 | 0 | — | 1.03 |
| `pooled_cl` (one cycle length) | −335.5 | −606.4 | **−3.27** (0.72) | **−4.53** | 1.02 |
| `pooled_R` | −330.9 | −609.4 | +1.25 (2.70) | 0.46 | **1.08 — VOID** |
| `pooled_both` | −330.4 | −609.9 | +1.80 (2.73) | 0.66 | 1.02 |

**Collapsing `cycle_length` to a single value costs held-out predictive
accuracy, by 3.27 log units at z = −4.5.** Read the direction carefully:
`no_pool` is the model that *keeps* the hierarchy. So the hierarchy earns its
keep, and **this reverses the lean of every weak comparison before it**,
which put the two indistinguishable (observation-level `loo`: −1.38, se 1.03,
z −1.34; and see "Cycle-length hierarchy").

Why this test can see what the others could not is the thing it was designed
for: each series keeps 3-6 points, so its initial conditions and error scale
stay informed, and `no_pool` can adapt `eta_cl` per trial where `pooled_cl`
cannot. A cycle-length error then accumulates as phase drift across the
held-out window instead of being absorbed by a neighbouring observation.

**`pooled_R` is void**: max R-hat 1.08, min n_eff 55.5. Nothing about that row
counts, per the project's own gate.

**`pooled_both` is not evidence against this.** Its +1.80 carries se 2.73 --
four times `pooled_cl`'s -- because it differs from `no_pool` in two
structures at once rather than one, so its predictions are less correlated
with the reference and the paired se is larger. The interval comfortably
contains `pooled_cl`'s −3.27. It neither supports nor contradicts; it is
imprecise.

**Pending**: the quarter-mask sensitivity (configs 27-30, SLURM 29575) is
running. Both fractions were fixed before any fit, so whichever way it reads
it is reported, not consulted and discarded. A verdict holding at one
fraction and not the other is a finding about the design's sensitivity.

### Budget and bound geometry at n=6

`_scripts/wockner-schedule-sim-analyze.R`, output `_data/schedsim-analyze.log`.
SLURM 29543 and 29544, all COMPLETED.

Cells: paired change in cycle-length bias in hours, arm minus `default` on
the same simulated dataset, **negative means the change reduced the bias**;
`n_pair` counts replicates where both converged.

| arm | `n_pair` | mean change | range | was (n=4) |
|---|---|---|---|---|
| `wide_bshape` | 6 | **−0.513 h** | −0.865 to −0.228 | −0.433 |
| `wide_nuis` | 6 | −0.460 h | −1.08 to −0.071 | −0.423 |
| `wide_total0` | 6 | +0.147 h | −0.156 to +0.430 | +0.099 |
| `cl_move` | 4 | **+0.033 h** | −0.111 to +0.157 | +0.059 |
| `cl_wide_move` | 5 | **+0.471 h** | +0.030 to +0.796 | +0.030 |

**The `b_shape` prior's share has gone 0.50 (n=3) → 0.43 (n=4) → 0.513
(n=6)**, with the range widening to 0.64 h. It is not settling: the
replicate spread remains comparable to the effect, which is the honest
headline of this whole decomposition.

**Bound geometry splits in two, and only now is it readable.** `cl_move`
keeps the 15 h width and centres it on the truth, isolating asymmetry: +0.033
h at n=4, nothing. `cl_wide_move` moves *and* widens, and at n=5 it is
**+0.471 h and positive in every replicate** -- widening the window makes the
bias worse. That is the confound named when the arms were submitted: holding
`cl_prior_center` at 48 h fixes the prior's location but not its width in
hours, and a wider window spans more hours at the same `sd_logit_cl`. So
`cl_wide_move` is measuring prior width, not bound distance. **Bound
asymmetry remains ruled out; neither arm is anywhere near the −1.5 h needed.**

**Convergence.** `default-rep4` converged on its **third** seed
(`SCHEDSIM_FIT_SEED=1123581321`: R-hat 1.01, 28 divergences, min ESS 240,
against 1.22/478/13 originally), so the decision rule's "stop after a third
failure" was not triggered. New failures appeared among the new replicates --
`default-rep5`, `default-rep6`, `wide-rep6`, `wide_nuis-rep6`,
`cl_move-rep2` on its new seed -- which is why `n_pair` is 6 and not 7.
Nuisance recovery barely moved at n=6-7: `default` `b_shape` 9.01 (was 8.77),
`R` −22.6%, `log10_total0` +99.5%.

### Design A's verdict is mask-dependent, so it does not settle the hierarchy

`_scripts/wockner-designA-score.R` with `PREFIX=daQ_` (output
`_data/designQ-score.log`) and `_scripts/designA-horizon.R` (output
`_data/designA-horizon.log`). Fits are configs 27-30, SLURM 29575, all
COMPLETED and all clearing R-hat < 1.05 -- including `pooled_R`, which was
void at the third mask.

**The pre-registered sensitivity disagrees with the headline.** Cells: as in
"Design A" -- held-out log pointwise predictive density summed over the
masked observations, and the paired difference against `no_pool` over those
same points with its standard error. Negative means the row predicts worse
than `no_pool`, which is the model that *keeps* the hierarchy.

| mask | held out | `pooled_cl` vs `no_pool` (se) | z |
|---|---|---|---|
| last third | 306 | **−3.27** (0.72) | **−4.53** |
| last quarter | 218 | **−0.04** (0.54) | **−0.08** |

Both fractions were fixed before any fit ran, so neither can be preferred for
reading better. **The disagreement is the result.**

**Where it comes from, and it is not the scoring window.** The masks nest:
both hold out the final `k` observations of each series, `k = floor(n/3)` and
`floor(n/4)`, so the quarter's 218 points are a subset of the third's 306 and
the 88 that differ lie deeper in each series' tail. That makes the third-mask
fits decomposable with no refitting.

Cells: paired difference against `no_pool` in log units, computed from the
**third-mask fits only**, split by subset of their held-out window;
`per obs` divides by the number of points so 218 and 88 are comparable.
Negative means worse than `no_pool`.

| model | shared 218: diff (se) | deeper 88: diff (se) | per obs, shared | per obs, deeper |
|---|---|---|---|---|
| `pooled_cl` | **−2.92** (0.66) | −0.35 (0.28) | −0.0134 | −0.0040 |
| `pooled_R` | −2.17 (2.51) | **+3.42** (0.95) | −0.0100 | +0.0389 |
| `pooled_both` | −1.49 (2.56) | **+3.29** (0.88) | −0.0069 | +0.0374 |

**89% of `pooled_cl`'s deficit sits on the 218 points the quarter mask also
holds out**, and per observation those points are *worse* (−0.0134) than the
88 deeper ones (−0.0040). So the obvious explanation -- that the hierarchy
needs a longer horizon for phase drift to accumulate -- is **wrong**. The
deficit is not concentrated where the extra horizon is.

**The difference is what the models were fitted on.** The same 218
observations score −2.92 under fits that never saw the 88 deepest points, and
−0.04 under fits that did. Withholding those 88 points from training is what
produces the hierarchy's advantage; restoring them removes it.

**So Design A does not establish that the hierarchy earns its keep.** It
gives a strong signal under one pre-registered mask and nothing under the
other, and the difference traces to the training set rather than to the
scoring window. **The previous section's reversal is withdrawn**: the
hierarchy question is not settled, and the honest statement remains the one
the weak comparisons gave -- no robust detectable benefit -- now supported by
a test that is not broken but is mask-sensitive.

**An observation that is not part of this and needs its own test**: `pooled_R`
and `pooled_both` predict the 88 deepest points *better* than `no_pool`, by
+3.42 (se 0.95) and +3.29 (se 0.88). `pooled_R` is void at this mask on the
earlier R-hat failure, but `pooled_both` is not. Pooling `R` helping deep
into a series is a different claim from anything examined here and rests on
one subset of one mask.

### The horizon ladder: no training-set trend, and the one signal is void

`_scripts/wockner-horizon-score.R`, output `_data/horizon-score.log`. Configs
31-36, SLURM 29591, all COMPLETED. Masks hold out the last k = 1, 2, 3
observations of each series, capped at n − 3, so 177, 350 and 479 of 1130.
They **nest**, so the k = 1 points are held out by every rung.

**Every row below is scored on the same 177 observations.** The horizon is
therefore identical across rows and only the training set varies, which is
the quantity Design A's decomposition identified and which neither of its
masks isolates.

Cells: log pointwise predictive density summed over those same 177 held-out
observations, computed per observation as `log(mean_s exp(log_lik[s,i]))`, in
log units, higher better. `withheld` is how many of the 1130 that rung's fits
were denied in training. `diff` is the paired `pooled_cl` minus `no_pool`
over those points with the standard error of the paired difference;
**negative means collapsing `cycle_length` predicts worse**. `max R-hat` is
each rung's worse of the two fits.

| rung | withheld | `no_pool` | `pooled_cl` | diff (se) | z | max R-hat |
|---|---|---|---|---|---|---|
| k = 1 | 177 | −186.1 | −186.4 | −0.29 (0.65) | −0.46 | 1.02 |
| k = 2 | 350 | −241.6 | −243.3 | **−1.73** (0.56) | **−3.08** | **1.15 — VOID** |
| k = 3 | 479 | −298.3 | −297.9 | +0.44 (0.64) | +0.70 | 1.03 |

**There is no training-set trend**: the slope is +0.0019 log units per
observation withheld, which is nothing, and the sign is not even consistent.
At k = 1 and k = 3, both of which converge, `pooled_cl` is indistinguishable
from `no_pool`.

**The only rung with a signal is the only rung that fails the gate.**
`daH2_no_pool` has max R-hat **1.15** and 312 divergences, so by the
project's own rule everything about that row is void. Its −1.73 cannot be
used, and a reseed of that one fit is what would complete the ladder.

**So Design A's third-mask result does not reproduce when the scoring window
is held fixed.** Together with its mask-dependence, the position is that
**no robust evidence has been produced that the hierarchy on `cycle_length`
earns its keep** — and equally none that it costs. That is where the weak
comparisons always were; what is new is that it now rests on tests that are
not broken.

### `fix_sd`: estimating the error scale is not the missing piece

`_scripts/wockner-schedule-sim.R` arm `fix_sd`, SLURM 29592, read by
`_scripts/wockner-schedule-sim-analyze.R`. The arm passes the simulation's own
`sd_iRBC` (27 values, 0.396-0.786) into `archer_stan_data()`, so it differs
from `default` in that alone.

This was the last **named** difference between the maximum-likelihood ladder
and the fitted posterior that had never been isolated: the MLE fixes
`sd_iRBC` at the truth where the fit estimates it.

Cells: paired change in cycle-length bias in hours, `fix_sd` minus `default`
on the same simulated dataset; **negative would mean fixing the error scale
reduced the bias**. `n_pair` is replicates where both converged.

| arm | `n_pair` | mean change | range |
|---|---|---|---|
| `fix_sd` | 2 | **+0.104 h** | +0.027 to +0.181 |
| `wide_bshape` (for scale) | 6 | −0.513 h | −0.865 to −0.228 |

**It does not reduce the bias; it slightly increases it**, and both available
replicates agree in sign. `fix_sd-rep1` failed the gate (max R-hat 1.06), so
this is n = 2 — thin, but the effect is an order of magnitude short of the
~1.5 h at issue and has the wrong sign, which n = 2 is adequate to establish.

**So the MLE-to-posterior gap is not the error scale.** With marginalisation
untestable by the MAP route — the joint surface has no usable mode — **there
is no named, untested component of that gap left.** Thread 2 is now
idea-limited rather than compute-limited, and that is the honest state of it.

### Phase-period compensation: refuted, and the fits mistime the oscillation

`_scripts/schedsim-phase-period.R`, outputs per replicate. Post-hoc on three
converged `default` simulation fits, no refitting.

The hypothesis, suggested by two numbers in the recovery table never read
together: `cycle_length` recovers +1.9 h and `b_offset` +0.115 cycles, and a
longer period loses phase over the window while a larger initial offset gains
it back. If the data pinned the **phase at observation times** and only
weakly separated period from offset, `cycle_length`'s bias would be one
coordinate of a flat ridge rather than an error the data could correct.

Cells: averaged over the 177 series. `err_cl_h` is posterior mean
`cycle_length` minus the truth, in **hours**; `cl_in_cycles` converts it to
the phase it costs over the mean observation span. The remaining columns are
in **cycles**: `err_bo` for the initial offset, `err_first`/`err_last` for the
phase `frac(b_offset + t / cycle_length)` at each series' first and last
observation, `sd_last` the posterior **circular** sd of that phase. `cor` is
the within-draw posterior correlation between a series' `cycle_length` and
its `b_offset`. Phase summaries are circular means.

| replicate | `err_cl_h` | `cl_in_cycles` | `err_bo` | `err_first` | `err_last` | `sd_last` | `cor` |
|---|---|---|---|---|---|---|---|
| rep1 | +2.56 | −0.231 | +0.118 | −0.314 | **−0.408** | 0.174 | +0.16 |
| rep2 | +1.12 | −0.101 | +0.067 | −0.148 | **−0.191** | 0.196 | −0.02 |
| rep3 | +2.13 | −0.192 | +0.134 | −0.258 | **−0.342** | 0.192 | −0.02 |

**The hypothesis is refuted, and backwards.** The phase at the last
observation is recovered about **1.8x worse** than the period that supposedly
produces it (0.408 against 0.231 cycles, 0.191 against 0.101, 0.342 against
0.192 — a strikingly consistent ratio). And `cycle_length` and `b_offset` are
essentially **uncorrelated** in the posterior (|r| <= 0.16), so there is no
ridge between them for the bias to be a coordinate of. `b_offset` does not
cancel the period error; it offsets about 60% of it on the group averages and
leaves a large residual.

**The "growing phase lag" is NOT a separate phenomenon, and an earlier
version of this section wrongly presented it as one.** Since
`phi(t) = b_offset + t / cycle_length`, the phase error decomposes
*algebraically* as `err_bo + t * (1/cl_est - 1/cl_true)` -- checked directly,
`max |predicted - observed| = 0` to machine precision. A period biased long
therefore MUST produce a lag growing linearly in t. Observing one is the
cycle-length bias restated, not evidence of anything further, and there is no
separate thing here to explain or fix.

**What the decomposition does show**, per series and unaveraged, for rep1 at
the last observation. Cells: quartiles over the 177 series of each term, in
cycles. `drift` is the period term `t * (1/cl_est - 1/cl_true)`, `err_bo` the
initial-offset error, `total` their sum.

| term | min | q25 | median | q75 | max |
|---|---|---|---|---|---|
| `drift` | −0.263 | −0.234 | **−0.213** | −0.203 | −0.186 |
| `err_bo` | −0.212 | −0.155 | **−0.016** | +0.365 | +0.495 |
| `total` | −0.418 | −0.389 | −0.249 | +0.151 | +0.252 |

**The period term is tight and systematic; the offset term is wide and
centred near zero.** So the phase error is driven by the period, with
`b_offset` adding spread rather than bias. Note this also qualifies the
nuisance-recovery table's "`b_offset` recovers +45%": that is a MEAN over
groups, and the median per-series error is −0.016 cycles. As with `b_shape`,
the mean is being set by a tail.

**A methodological note worth keeping.** The first version of this script
took an arithmetic mean of wrapped phase differences. On a quantity whose
posterior spread is a quarter of a cycle that is pulled toward zero, and it
reported phase errors roughly a third too small -- small enough to look like
the hypothesis was supported. Circular quantities need circular means; the
numbers above are circular throughout.

### The cycle-length prior's location was never tested, only its width

Arithmetic, no fitting: `qlogis`/`plogis` over the prior
`logit((cl - 35) / 15) ~ normal(logit((center - 35) / 15), sd_logit_cl)`.

Cells: the implied prior on `cycle_length` in **hours**, by Monte Carlo over
2e6 draws. `center` is `cl_prior_center`. The simulated truth is 45.012 h.

| `center` | `sd_logit_cl` | prior mean | prior median | 90% interval |
|---|---|---|---|---|
| 48 | 1 | **47.44** | 48.00 | 43.3-49.6 |
| 48 | 2 | 46.41 | 48.00 | 37.9-49.9 |
| 48 | 3 | 45.57 | 48.00 | 35.7-50.0 |
| 45 | 1 | 44.60 | 45.01 | -- |
| 42 | 1 | 42.09 | 42.00 | -- |

**The prior's mean sits 2.4 h above the truth**, and the observed posterior
means (46.2-47.6 h) lie between the two. **Widening `sd_logit_cl` leaves the
median pinned at 48.00** and drags the mean down only slowly through the
bounds -- so the test that closed this question ("quadrupling the prior
variance moves the estimate 0.15 h") varied the prior's **width** and never
its **location**. The normal-normal weight calculation elsewhere in this file
(+0.271 h) is also computed on the logit scale, where the asymmetry that
produces the 2.4 h gap does not appear.

A crude check from the arms already run: `wide` (sd 2) moves the prior mean
1.03 h and the measured bias by 0.15 h, implying a prior weight near 0.15,
which against a 2.4 h displacement would be ~0.36 h of upward bias -- more
than the 0.15 h on record, less than the ~1.5 h at issue. **That is an
estimate from a linear argument on a non-linear prior and should not be
quoted**; the direct test is an arm with `cl_prior_center` at the truth,
which is diagnostic in simulation exactly as `fix_sd` was, and is the one
cheap thing thread 2 has left. `wockner-fit.R` configs 5 and 6
(`np_center45`, `np_center42`) were written for the real data and never
submitted.

### Literature: the bias is a known property of this estimation problem

Checked before searching outward, per the project's own habit: `mmcm.pdf`
(Greischar & Childs, *Trends in Parasitology* 39(8), 2023) sits at the repo
root, untracked, and is directly about estimation bias in exactly this
system. Supplemented with
an OpenAlex search; the relevant works are recorded in `claude/scripts.md`'s
reference list rather than re-summarised here.

**What `mmcm.pdf` establishes, for PMR rather than cycle length:**

- Estimates "systematically overestimate the multiplication rate in
  **synchronized** infections", and "the estimates for initially
  **asynchronous** infections are close to the true values and do not vary
  depending on the initial median parasite age".
- "the largest errors in PMR estimates when the initial median parasite age
  is **offset from the sampling time by roughly 12 h**", because "some of the
  samples occur when most parasites are sequestered".
- "Even moderate levels of synchrony generate exaggerated PMRs ... estimates
  are especially poor with higher levels of synchrony."

**Why this maps onto thread 2.** The two quantities that paper identifies as
controlling the bias are exactly the two this model parameterises and cannot
identify: synchrony is `b_shape` and initial median parasite age is
`b_offset`. The Wockner sampling grid is **12 h** for 87.8% of intervals,
which is the offset the paper reports as worst. And Design A pinned `b_shape`
at 400 — near-maximal synchrony, which that paper reports as the regime where
estimates are *especially* poor.

**What it does not establish.** It is about PMR (`R` here), not
`cycle_length`, and it uses a different estimator (regression and maximum
observed ratio, not a fitted mechanistic posterior). It makes the
phase-and-synchrony dependence a **hypothesis worth testing here**, not a
result that transfers.

**The testable prediction**: the cycle-length bias should vary systematically
with the true `b_offset` relative to the 12 h sampling grid, and should be
larger at high `b_shape`. Every simulation replicate so far shares one true
`b_offset`, so nothing run to date could have seen this. Varying it is one
new arm dimension.

**Also relevant, and it cuts the other way on `cl_prior_center`**: the IDC is
reported to complete in multiples of 24 h under circadian coordination
(Subudhi et al. 2020, *Nat Commun*), so `cl_prior_center = 48` is a
biologically motivated choice, not an arbitrary one. Whether to move it is a
scientific judgement about which evidence to believe, not a technical fix —
which is thread 7, and separate from measuring the prior's weight (SLURM
29612).

### Phase volume at commensurate periods: refuted

`_scripts/phase-volume-test.R`, output `_data/phase-volume-test.log`, on
`default-rep1`. No MCMC: every parameter held at the simulated truth except a
common phase shift applied to all `b_offset` groups, scanned over
`cycle_length` 42-50 h at 0.25 h by 5 degrees.

The hypothesis was that periods commensurate with the 12 h sampling grid --
87.8% of within-series intervals -- carry extra parameter-space volume,
because fewer distinct phases are observed there and more `b_offset` fits
equally well. A posterior mean integrates over volume where a maximum does
not, which is the split this project keeps seeing.

The coverage asymmetry is real and large. Cells: `eff. phases` is
`1 / sum(p_i^2)` over 36 phase bins of all 1130 observation times at that
period, 36 if uniform; `max gap` the largest gap between sorted phases, in
cycles. Deterministic, from the observation times alone.

| period | eff. phases | max gap |
|---|---|---|
| 45.012 h (truth) | 7.8 | 0.134 |
| 48 h (prior centre) | **3.7** | 0.250 |

**But it does not become likelihood volume.** The profile likelihood (maximum
over phase) and the phase-integrated likelihood peak at **the same place**,
45.75 h, and the integrated-minus-profile difference at 48 h is **−0.10 log
units** relative to the truth -- essentially nothing, and the wrong sign.
Degenerate phase coverage does not mean more of phase space fits; the
likelihood is set by trajectory shape, not by how many distinct phases are
sampled. **Refuted.**

Worth keeping from it: with every nuisance at the truth, the profile
likelihood on this replicate still peaks **+0.74 h** high. Part of the bias
is present in the likelihood itself, before any prior or marginalisation.

### The cycle-length prior's location is worth ~0.44 h, not 0.15 h, and not the answer

SLURM 29612, arms `cl_center_truth` and `cl_center_low`, 6 fits, all
COMPLETED and all converged. Read with
`_scripts/wockner-schedule-sim-analyze.R`.

Cells: paired change in cycle-length bias in hours, arm minus `default` on
the same simulated dataset, **negative means the change reduced the bias**;
`displacement` is that arm's prior MEAN in hours minus the simulated truth of
45.012 h, so `default` sits at +2.43 h by construction and is the reference.
`slope` is the bias change over the displacement change from `default`,
which is the prior's weight.

| arm | displacement | `n_pair` | mean change | range | slope |
|---|---|---|---|---|---|
| `default` | +2.43 h | — | 0 (reference) | — | — |
| `cl_center_truth` | −0.41 h | 3 | **−0.436 h** | −0.51 to −0.29 | 0.154 |
| `cl_center_low` | −2.93 h | 3 | **−0.637 h** | −0.76 to −0.51 | 0.119 |

**The prior's location is worth about three times what the record says.**
Centring it removes 0.44 h where the figure on file was 0.15 h — because
that figure came from varying the prior's **width**, which leaves the median
pinned at 48.00 and barely moves the mean.

**The slope is ~0.12-0.15 and it flattens**: the first 2.84 h of displacement
buys 0.436 h, the next 2.52 h buys only 0.201 h. So the prior behaves like a
weight near 0.15 that saturates, not like a free lever.

**And it is not the explanation.** With the prior centred on the truth the
residual bias is **+2.07, +0.85, +1.68 h**, mean **1.53 h**. Pushing the
prior 2.93 h *below* the truth still leaves ~1.3 h. **A correctly located
prior does not fix this**, which settles a question that has been ambiguous
in this file since the schedule-bias simulation: the prior is worth more than
recorded and is still not the problem.

### Budget, everything measured to date

Cells: paired change in cycle-length bias in hours against `default` on the
same simulated datasets; **negative reduces the bias**. `n` is paired
replicates. These are single interventions against a common baseline and are
**not additive** -- `wide_nuis` already shows that widening both nuisance
priors (−0.46) is not the sum of widening each (−0.51 and +0.15).

| intervention | n | change | status |
|---|---|---|---|
| `b_shape` prior widened | 6 | **−0.513 h** | the largest single lever |
| both nuisance priors widened | 6 | −0.460 h | not additive |
| prior located on the truth | 3 | **−0.436 h** | new; was recorded as 0.15 h |
| prior located 2.9 h below truth | 3 | −0.637 h | saturating |
| hierarchy removed (`no_hier`) | 2 | ~−0.25 h | from the earlier ladder |
| `log10_total0` prior widened | 6 | +0.147 h | wrong sign |
| `sd_iRBC` fixed at truth | 2 | +0.104 h | wrong sign |
| bounds centred, same width | 4 | +0.033 h | nothing |
| bounds moved and widened | 5 | +0.471 h | prior-width confound |

**Against a `default` bias of about +1.97 h, the best single intervention
removes 0.51 h and the two largest together could not remove 1 h.** Every
named mechanism has now been measured. The remainder, ~1.0-1.5 h, has no
candidate attached to it.

### The horizon ladder, completed: flat

`_scripts/wockner-horizon-score.R`. `daH2_no_pool` failed the gate on its
first seed (max R-hat 1.154, 312 divergences, min ESS 28) and was refit as
`daH2_no_pool-seed2` (SLURM 29622: **1.017, 88, 165**). The reseed is used in
its place and the failed fit stays on disk. All six fits now converge.

Cells as before: every row scored on the **same 177 observations** (the k = 1
points, held out by all three rungs), so the horizon is fixed and only the
training set varies. `diff` is paired `pooled_cl` minus `no_pool` over those
points with the standard error of the paired difference; negative means
collapsing `cycle_length` predicts worse.

| rung | withheld | diff (se) | z | was, with the failed fit |
|---|---|---|---|---|
| k = 1 | 177 | −0.29 (0.65) | −0.46 | −0.29 |
| k = 2 | 350 | **−1.26** (0.71) | **−1.78** | −1.73, z −3.08 |
| k = 3 | 479 | +0.44 (0.64) | +0.70 | +0.44 |

**No rung now reaches |z| = 2**, and the slope is +0.0020 log units per
observation withheld — flat, with an inconsistent sign. The one signal in the
earlier version was carried by the fit that failed to converge; on a seed
that converges it falls to z = −1.78.

**So Design A's third-mask result does not reproduce when the scoring window
is held fixed.** With the mask-dependence already on record, the position is
settled as far as this design can settle it: **no robust evidence that the
hierarchy on `cycle_length` earns its keep, and none that it costs.**

### Where the bias sits: the likelihood explains about a fifth of it

`_scripts/phase-volume-test.R` on all three converged `default` replicates.
Every parameter held at the simulated truth except a common phase shift.

Cells: `profile` is where the maximum-over-phase log likelihood peaks and
`integrated` where the phase-integrated likelihood peaks, both in hours, with
the bias against the 45.012 h truth; `volume` is the integrated-minus-profile
difference at 48 h relative to the truth, in log units, positive meaning more
phase-space fits at 48. `posterior` is that replicate's full fitted bias from
the schedule simulation, for scale. Grid 0.25 h by 5 degrees.

| replicate | profile | integrated | volume at 48 | posterior |
|---|---|---|---|---|
| rep1 | +0.74 | +0.74 | −0.10 | +2.58 |
| rep2 | −0.01 | +0.24 | −0.06 | +1.14 |
| rep3 | +0.49 | +0.49 | −0.18 | +2.19 |
| **mean** | **+0.41** | **+0.49** | — | **+1.97** |

**The volume hypothesis is refuted in all three.** Profile and integrated
peak at the same place in two of them and one grid step apart in the third,
and the differential volume at 48 h is negative every time — less phase-space
fits there, not more.

**The likelihood carries about a fifth of the bias.** With every nuisance
known, the likelihood still prefers +0.41 h on average. It is not a constant
offset: +0.74, −0.01, +0.49, so it is noise-driven. And it **tracks the full
bias across replicates** — rep1 highest in both, rep2 lowest in both, rep3
between — so the same noise realisations that make the posterior worse also
make the likelihood worse.

**The other four fifths come from estimating the nuisances.** ~1.56 h appears
only when `b_shape`, `R`, `log10_total0`, `sd_iRBC` and `b_offset` are
estimated rather than known, and the prior-location test accounts for 0.44 h
of that. **This is now the sharpest statement of thread 2**: the bias is
mostly the cost of not knowing the nuisance parameters, not a defect in the
likelihood, the prior, the bounds, the sampling grid, or the error scale.

**The test that follows directly**, and which the package change already
makes possible: fix each nuisance at its simulated truth one at a time and
see which one's estimation carries the 1.56 h. `fix_sd` is done and does not
(+0.104 h). **`b_shape` fixed at the truth is untested**, is the parameter
this project has repeatedly found unidentified, and is one arm.

### `b_shape` is the nuisance whose estimation carries the bias

`_scripts/wockner-schedule-sim-fixbshape.sh` (SLURM 29628, tasks 85-87),
read by `_scripts/wockner-schedule-sim-analyze.R`, saved output
`_data/schedsim-analyze-2026-10-05.txt`.

`fix_bshape` pins `b_shape` at the 14 values the data were simulated from
(range 10.41-21.78) and estimates everything else, so it pairs against
`default` on the identical simulated datasets.

Cells: `delta` is that replicate's mean posterior `cycle_length` minus
`default`'s on the same simulated data, in hours; negative means the arm
reduces the bias. `n_pair` counts replicates where both arms converged.
`default` totals +1.97 h over reps 1-3.

| arm | what it removes | n_pair | mean delta | range |
|---|---|---|---|---|
| **`fix_bshape`** | `b_shape` uncertainty entirely | **1** | **−0.919** | — |
| `wide_bshape` | most of the `b_shape` prior's pull | 6 | **−0.513** | −0.865 to −0.228 |
| `wide_nuis` | all five nuisance priors at once | 6 | −0.46 | −1.08 to −0.071 |
| `cl_center_low` | prior location on `cycle_length` | 3 | −0.637 | −0.757 to −0.508 |
| `cl_center_truth` | prior location, set to truth | 3 | −0.436 | −0.51 to −0.289 |
| `fix_sd` | `sd_iRBC` uncertainty entirely | 2 | +0.104 | +0.027 to +0.181 |
| `wide_total0` | the `log10_total0` prior | 6 | +0.147 | −0.156 to +0.43 |

**Knowing `b_shape` removes about half the bias.** On the one replicate that
converged, −0.919 h of +2.58 h. Across all three unfiltered the paired
differences are −0.919, −1.151, and −0.974, mean −1.015 h: the answer does
not depend on which replicates are admitted, which matters because two of
the three `fix_bshape` fits are excluded as non-converged (below).

**The corroboration that does not rest on those fits is `wide_bshape`**:
6 converged pairs, every one negative, mean −0.513 h. Merely widening the
`b_shape` prior buys half of what pinning it at the truth buys. That is a
dose-response in the right direction, from well-converged fits, and it is
why the conclusion stands despite `fix_bshape`'s n_pair of 1.

**And `b_shape` carries essentially all of the nuisance-prior effect.**
`wide_nuis` widens all five and gains −0.46 h, no more than `wide_bshape`
alone at −0.513 h. `fix_sd` and `wide_total0` are both nil and the wrong
sign. Of the nuisances, only `b_shape` moves `cycle_length`.

**Budget.** Of +1.97 h: `b_shape` estimation ~0.92 h, the likelihood with
every nuisance known +0.41 h, prior location ~0.44 h. Those sum to 1.77 of
1.97. **The arms cannot simply be added** -- `cl_center_truth` and
`fix_bshape` both work by removing uncertainty and plausibly overlap -- so
read this as the pieces being of the right order to close the budget, not
as an exact decomposition.

#### What this does and does not say about thread 5

It says `b_shape` is **weakly identified**, and that paying to estimate it
costs `cycle_length` about an hour. It does **not** say the desynchronisation
rate is wrong: in simulation `n_c` is identical in the generating and fitted
model, so the decay rate is correct by construction, exactly as thread 5
already records. These are two separate claims and only the first is tested
here.

The link is that both concern the same parameter. `b_shape` is the only free
synchrony knob, the chain's own decay dominates the late-window amplitude,
and `b_shape` runs to 65+ when its prior is relaxed. If on **real** data the
fixed decay rate is wrong, `b_shape` is where that error would land -- and
this result shows that errors landing in `b_shape` propagate into
`cycle_length`. That raises thread 5 from a misspecification risk with no
known consequence to one with a measured channel to the headline estimate.

#### Sampler health, and why two of three are being refit

Fixing a parameter at the truth made sampling **worse** in every replicate,
not better:

Cells: divergences out of 4000 post-warmup draws, max R-hat, min bulk ESS.

| rep | `default` | `fix_bshape` |
|---|---|---|
| 1 | 53 div, 1.030, ESS 93 | 67 div, 1.029, ESS 200 |
| 2 | 25 div, 1.015, ESS 204 | **586 div, 1.343, ESS 10** |
| 3 | 35 div, 1.022, ESS 169 | **222 div, 1.066, ESS 53** |

The analyzer's convergence filter drops reps 2 and 3, leaving n_pair = 1.
Reps 2 and 3 are being refit on a second sampler seed
(`_scripts/wockner-schedule-sim-fixbshape-seed2.sh`); the simulated data are
unchanged, so they still pair against the same `default` fits. The reseed is
triggered by the diagnostics alone and would have been run whichever way the
biases pointed -- precedent is the `daH2_no_pool` reseed, where R-hat
1.154 -> 1.017 moved a rung from z −3.08 to z −1.78.

That the geometry gets harder when a parameter is **removed** is itself
worth noting: it is the signature of `b_shape` having been absorbing
something, so that pinning it forces the conflict into the parameters left.

### The `fix_bshape` reseeds: the answer reproduces, the gate still fails

SLURM 29631, `SCHEDSIM_FIT_SEED=2`, reps 2 and 3. Simulated data unchanged, so
they pair against the same `default` fits.

Cells: `bias` is mean posterior `cycle_length` minus the 45.01 h truth, in
hours. `delta` is that minus `default`'s on the same data.

| rep | seed | R-hat | div | min ESS | bias | delta vs `default` |
|---|---|---|---|---|---|---|
| 2 | original | 1.343 | 586 | 10.2 | −0.011 | −1.151 |
| 2 | **seed2** | **1.055** | 108 | 60.1 | **−0.117** | **−1.257** |
| 3 | original | 1.066 | 222 | 52.9 | +1.219 | −0.974 |
| 3 | **seed2** | **1.052** | 255 | 65.9 | **+1.239** | **−0.955** |

**The biases reproduce.** rep2 moves −0.011 → −0.117 and rep3 +1.219 →
+1.239, while R-hat falls 1.343 → 1.055 and ESS rises 10 → 60. So rep2's
near-zero bias was **not** an artefact of the bad geometry, which was the
live alternative. The paired differences across all three replicates are
−0.919, −1.257, −0.955, mean **−1.044 h**, against −1.015 h before.

**But both reseeds still fail the 1.05 convergence gate** (1.055 and 1.052),
so the analyzer still reports `fix_bshape` at **n_pair = 1**. The
pre-registered rule applies and is being followed: report n_pair = 1 with the
caveat, lean on `wide_bshape`'s −0.513 h over 6 clean pairs, and **do not
reseed a third time**. Two seeds that both stall just above the gate on the
same two replicates is itself informative -- pinning `b_shape` makes the
geometry harder, which is consistent with it having been absorbing something.

### `n_c` is a biological assumption, not a numerical setting

SLURM 29633, entries 37-38. `np_wide_both` (`n_c = 96`), `np_nc192`,
`np_nc384` differ in `n_c` and nothing else. Scripts
`_scripts/nc-ladder-read.R` and `_scripts/nc-numerical-check.R`, saved outputs
`_data/nc-ladder-2026-10-06.txt`, `_data/nc-numerical-check-2026-10-06.txt`,
`_data/schedsim-analyze-2026-10-06.txt`.

Cells: `b_shape` and `cycle_length` are posterior means averaged over the 14
and 13 groups; `cycle_length` in hours. `elpd` is observation-level PSIS-LOO.
`stage sd` is the deterministic stage-distribution sd after the 4.80-cycle
Wockner window, `sqrt(k / n_c)` cycles.

| `n_c` | stage sd | `b_shape` | `cycle_length` | elpd | R-hat | div |
|---|---|---|---|---|---|---|
| 96 | 0.224 cyc | 64.45 | 44.25 | −926.2 | 1.019 | 137 |
| 192 | 0.158 cyc | 61.57 | 42.07 | −857.6 | 1.032 | 165 |
| 384 | 0.112 cyc | **25.94** | **41.04** | **−854.4** | 1.055 | 54 |

**`b_shape` falls monotonically, as predicted**: 64.45 → 61.57 → 25.94. The
prediction was recorded in `wockner-fit.R` before the fits ran. `max_shape` is
250 and the largest group mean is 100.6, so this is not bound truncation.

**`cycle_length` moves −3.20 h.** That is larger than the entire +1.97 h
simulated bias this project has spent weeks on, and it is driven by a setting
`check_erlang_window()` picks for numerical accuracy.

**`n_c = 96` is 71.8 elpd worse than 384** (se 13.2, so z = −5.4). 192 against
384 is −3.2 (se 5.0): indistinguishable. **The ladder plateaus at 192.**

#### The effect is structural, not numerical

`n_c` does two jobs -- it sets the Erlang shape, and so the desynchronisation
rate, and it has to be large enough for the series solution to approximate the
matrix exponential. These have opposite implications, so they were separated.

`max_rel_diff` is the model's own cross-check: generated quantities recomputes
the trajectory with `matrix_exp` and reports the largest relative difference
against the series solution the likelihood uses.

Cells: largest and median relative difference over every observation and
every draw of a short probe run at each `n_c`, with `run_check = 1`.

| `n_c` | max `max_rel_diff` | median |
|---|---|---|
| 96 | 8.05e−13 | 5.04e−13 |
| 192 | 1.39e−12 | 1.09e−12 |
| 384 | 2.81e−12 | 2.07e−12 |

**At `n_c = 96` the arithmetic is faithful to one part in 10^12.**
`check_erlang_window()` was doing its job correctly. The `n_c = 96` model is
computed accurately and simply **fits the data worse**, by a wide margin. So
this is a statement about the biology the model assumes, not about its
numerics, and no earlier fit is suspect arithmetic.

**And the error grows mildly WITH `n_c`** -- 8e−13, 1.4e−12, 2.8e−12 -- which
is expected from more compartments accumulating more rounding. This rules the
numerical explanation out rather than merely failing to support it: if the
71.8 elpd gain came from better arithmetic at higher `n_c`, the higher rungs
would have to be more accurate, and they are slightly less.

#### What this does and does not license

It does **not** yet license changing any reported cycle length. These fits run
with `b_shape` **estimated** (`sd_log_b_shape = 1.5`), which is the right
setting for asking whether `b_shape` absorbs the decay rate -- it has to be
free to move -- but it is **not** the configuration the project has settled
on, where `b_shape` is pinned at 400 for +28 to +37 elpd. Entries 39-40
(`np_bs400_nc192`, `np_bs400_nc384`) ask whether the `n_c` effect survives a
pinned `b_shape`, with the three outcomes recorded in advance.

Caveat on the top rung: `np_nc384` has max R-hat 1.055, just above the 1.05
gate used elsewhere in this project. Its elpd is statistically tied with
`np_nc192` (R-hat 1.032), which is the rung the plateau conclusion rests on,
so the headline comparison does not depend on the marginal fit.

#### Sizing lesson

`nc-sizing.R` predicted ~2.9 h at `n_c = 192`; it took **5 h 47**, and 384
took **15 h 59**. The probe measured cost **per leapfrog** (1.81x from 96 to
192) and was right about that, but the leapfrog count **also** rose, 223 → 390
on the production fit. Per-gradient cost and the number of gradients are
separate multipliers and a probe that fixes the iteration count only measures
the first. Multiply them next time, or treat a per-leapfrog probe as a lower
bound.

### The `n_c` 2x2: separating estimator bias from misspecification bias

**IN FLIGHT, SLURM 29637.** Recorded here before the results so the reading
rule is fixed in advance.

The real-data ladder cannot say which end is right, because real data has no
truth to miss, and elpd scores prediction on the observed window -- which an
amplitude/period tradeoff can win while getting the period wrong. These arms
put a truth in. `wockner-schedule-sim.R` now decouples the `n_c` used to
GENERATE (`arm$sim_n_c`) from the one used to FIT (`arm$build$n_c`):

| | fit `n_c` = 96 | fit `n_c` = 192 |
|---|---|---|
| **generate 96** | `default` (done, +1.97 h) | `sim96_fit192` |
| **generate 192** | `sim192_fit96` | `sim192_fit192` |

The diagonal is correctly specified and measures **estimator** bias at each
`n_c`. The off-diagonal measures what **misspecifying** `n_c` manufactures.

**`sim96_fit192` is the cell that bears on the real-data result**, and it is
the one that could overturn it. If fitting ABOVE the true `n_c` drags
`cycle_length` down on data where the truth is 96, then the real-data −3.20 h
is an artefact of over-large `n_c` rather than a correction, and 71.8 elpd
would not settle the period.

Pairing, which the data support rather than assume: `sim96_fit192` generates
at 96 with the same noise seed as `default`, so it shares its simulated data
exactly -- verified, both give y_sim mean log10 2.7590 on rep 1 -- and it is
in `PAIRED_ARMS`. The two gen-192 arms share data with **each other** (y_sim
mean 2.6934), **not** with `default`, so they are deliberately **not** in
`PAIRED_ARMS`; differencing them against it would difference two datasets
rather than two fits. They are read against the true `cycle_length` instead,
by `_scripts/nc-2x2-read.R`.

#### Why this is worth the compute

The real-data move already has an arithmetic coincidence behind it. The
simulation says the `n_c = 96` estimator reads **+1.97 h high**. Applying that
to the real-data fit:

| route | value |
|---|---|
| `np_wide_both` at `n_c = 96` | 44.24 h |
| minus the simulated +1.97 h bias | **42.27 h** |
| `np_nc192`, measured | **42.07 h** |

Two independent routes land within 0.2 h, which would make the +1.97 h bias
and the `n_c` misspecification one phenomenon rather than two. **That is a
hypothesis, not a result**: the +1.97 h was measured generating AND fitting at
96, so it is estimator bias inside a correctly specified model, and the
agreement may be coincidence. The 2x2 is what tells the difference.

#### The mechanism both results share

The model can match the observed oscillation with high amplitude and a short
period, or low amplitude and a long one. `n_c` fixes the decay rate, and at 96
the chain desynchronises fast, killing late-window amplitude; the fit buys it
back by raising `b_shape` and lengthening the period. Raise `n_c` and
amplitude comes free, so both compensations relax -- which is why `b_shape`
and `cycle_length` fell **together** along the ladder. The same trade shows at
fixed `n_c`: pinning `b_shape` at 400 gives 43.65 h against `np_wide_both`'s
44.24 h.

**So `cycle_length` is set by an amplitude-period tradeoff with `n_c` fixed by
assumption on one side of it.** That is the real reason the headline number
was never as well determined as it looked, and it is independent of which end
of the ladder turns out to be right.

### The `n_c` 2x2: misspecifying `n_c` moves `cycle_length` in both directions

SLURM 29637. Scripts `_scripts/nc-2x2-read.R` and
`_scripts/nc-bias-correct.R`; saved outputs `_data/nc-2x2-2026-10-07.txt`,
`_data/nc-bias-correct-2026-10-07.txt`.

Cells: mean posterior `cycle_length` over trials minus the 45.012 h simulated
truth, in hours, averaged over converged replicates (max R-hat < 1.05) among
reps 1-3. The diagonal is correctly specified.

| | fit 96 | fit 192 |
|---|---|---|
| **gen 96** | **+1.970** (n=3) | **+0.928** (n=3) |
| **gen 192** | **+3.000** (n=2) | **+1.547** (n=3) |

**Both misspecification directions are real and roughly symmetric.**
Under-specifying (`gen 192, fit 96`) inflates `cycle_length` by **+1.454 h**
over the correctly specified cell; over-specifying (`gen 96, fit 192`)
deflates it by **−1.042 h**. So an `n_c` wrong by a factor of two shifts the
estimate by about an hour, in the direction you would expect.

**`sim96_fit192` did not come back near zero**, which was one of the three
pre-registered readings. So over-specifying `n_c` **does** drag
`cycle_length` down, and part of the real-data −2.17 h from 96 to 192 could
be that artefact. But only part: the manufactured amount is −1.04 h against
an observed −2.17 h, and the real-data elpd *preferred* 192 by 68.6, which a
genuinely over-specified model should not do.

**Estimator bias is also lower at the higher `n_c` when correctly
specified**: +1.547 at 192/192 against +1.970 at 96/96. Fitting at a higher
`n_c` is mildly better behaved even with no misspecification to fix.

#### Which hypothesis reconciles the real data

Two rungs corrected under the same hypothesis are two measurements of one
quantity and should agree. Real fits are `np_wide_both` 44.24 h and
`np_nc192` 42.07 h.

Cells: real fitted `cycle_length` minus the bias the 2x2 measured for
(true `n_c` = H, fitted `n_c` = that rung).

| true `n_c` | rung 96 | rung 192 | spread |
|---|---|---|---|
| 96 | 42.27 | 41.15 | **1.13 h** |
| **192** | **41.24** | **40.53** | **0.72 h** |

**`true n_c = 192` reconciles the two rungs better**, 0.72 h against 1.13 h,
and it agrees in direction with the elpd ordering. This is suggestive, not
decisive: the gap between 1.13 and 0.72 is well inside the replicate scatter.

**What both hypotheses agree on**: every corrected estimate is **below every
raw one**. The raw ladder runs 44.24 → 42.07 → 41.04; the corrected values
run 40.5 to 42.3. So the real cycle length is probably in the low 40s and
**the +1.97 h bias and the `n_c` sensitivity push the same way**, which is
the first direct support for them being one phenomenon rather than two.

**Caveats that limit this**, both written into the script header: the 2x2
arms use the schedule simulation's priors (`sd_log_b_shape` 0.5) while the
real ladder uses the widened one (1.5), so the biases transfer only
approximately; and the simulated truth of 45.012 h is itself taken from an
`n_c = 96` fit, so if 96 inflates then the biases were measured in the wrong
neighbourhood. `gen 192, fit 96` also has only 2 converged replicates.

### The production `n_c` ladder did not converge, and is not reportable

SLURM 29635, entries 39-40 (`np_bs400_nc192`, `np_bs400_nc384`).

| fit | max R-hat | div | `lp__` by chain |
|---|---|---|---|
| `np_bs400_data` (96) | 1.047 | 191 | — |
| `np_bs400_nc192` | **6.13** | 275 | −929.0, −832.3, −833.2, −836.7 |
| `np_bs400_nc384` | **8.34** | 158 | −996.0, −925.6, −955.3, −957.8 |

**One chain is stuck far below the others** in each, by 96 and 70 log units.
The majority chains agree with each other -- `cycle_length[1]` is 42.18,
42.28, 42.10 with the outlier at 43.09 -- so the information is probably
there, but **a run with a stuck chain is not a posterior**. The fitted values
(41.24 h, 38.73 h) and the `loo_compare` from this job are **not usable**:
`log_lik` itself has R-hat near 6, so the LOO is as unreliable as the
parameters. **Nothing here licenses any statement about whether the `n_c`
effect survives a pinned `b_shape`** -- that question is still open.

Reporting the three agreeing chains would be choosing the chains that give a
tidy answer, so it is not done.

**Why it probably failed.** Pinning `b_shape` is already known to harden the
geometry: `fix_bshape` found exactly that in simulation, where two seeds both
stalled just above the gate. Pinning at 400 -- very tight synchrony -- and
raising `n_c` appears to compound it. The likely mechanism is phase
multimodality: sharply synchronised parasites with a slightly mismatched
period admit more than one local phase alignment, and `b_offset` is an
`array[n] unit_vector[2]`, so those are separate modes rather than a ridge.

**The retry** is entries 41-42, `_scripts/wockner-fit-nc-bs400-retry.sh`:
`adapt_delta` 0.95, `max_treedepth` 12, and a different seed. It changes two
things at once on purpose, because the aim is to get a converged answer
rather than to attribute the failure. **If a chain still sticks at a similar
lp gap, the mode is real**: report it as multimodality, show both modes, and
do not chase it with a third configuration.

### Why `cycle_length` shifts with `n_c`: the observable period is not the parameter

`_scripts/nc-mechanism.R`, saved output `_data/nc-mechanism-2026-10-07.txt`.
Everything here is **deterministic and noiseless** -- no priors, no sampler,
no estimation -- so whatever appears is built into the model's forward map.

Reading the source turned up two channels that `n_c` controls, with opposite
status.

#### Channel A, sequestration-grid discretisation: RULED OUT, wrong sign

`make_log_y_vals` puts the circulating fraction on a logistic in **absolute
developmental age**, centred at `p3 = 18.5802` h, and evaluates it at age
`k * cycle_length / n_c` for compartment `k`. The grid spacing is therefore
`cycle_length / n_c` hours, and the source documents an O(1/n_c) delay in
sequestration onset from the Archer convention `q[1] = 0`. That looked like a
candidate.

Cells: `duty_cycle` is `mean(y_vals)`, the average probability of not being
sequestered over one cycle, at `cycle_length` = 45.012 h.

| `n_c` | step (h) | duty cycle |
|---|---|---|
| 96 | 0.469 | 0.39765 |
| 192 | 0.234 | 0.40019 |
| 384 | 0.117 | 0.40146 |
| 1536 | 0.029 | 0.40242 |

The duty cycle **rises** with `n_c`. Since `d(duty)/d(cycle_length)` is
−0.00895 per hour, offsetting the 96 → 192 change of +0.00254 would require
`cycle_length` to **increase by +0.284 h**. The observed shift is **−2.17 h**.
**Wrong sign and roughly eight times too small**, so this channel is not the
explanation and in fact slightly opposes the effect. It does converge away as
`n_c` grows, as a discretisation artefact should.

#### Channel B, the one-sided observation window: CONFIRMED

`_scripts/nc-period-check.R`, saved output
`_data/nc-period-check-2026-10-07.txt`. Peak times are refined sub-grid by a
parabola through each maximum and its neighbours; without that every interval
lands on the 0.25 h sampling grid, which is coarser than the differences
being measured.

Cells: the third peak-to-peak interval of the detrended log10 trajectory
minus the `cycle_length` that generated it, in hours, so a negative value
means the observable period runs SHORT of the parameter. Amplitude is the sd
of the detrended series across `n_c` 96 → 192 → 384.

| `b_shape` | `n_c` = 96 | 192 | 384 | 96 → 384 | amplitude |
|---|---|---|---|---|---|
| 15 | **−0.77** | −0.28 | −0.05 | **0.72 h** | 0.305 → 0.405 → 0.496 |
| 400 | **−0.74** | −0.14 | +0.24 | **0.98 h** | 0.390 → 0.547 → 0.719 |

**At `n_c = 96` the observable period is about 0.75 h SHORTER than the
`cycle_length` parameter that produced it, and the deficit closes as `n_c`
rises.** So a fit at `n_c = 96` must **inflate** `cycle_length` to reproduce
a given observed period. That is the sign seen on real data.

**The mechanism.** Sequestration at 18.58 h of a ~45 h cycle means only the
first ~40% of the cycle is visible: the observation window is **one-sided**.
The stage distribution broadens at a rate `n_c` fixes -- sd after k cycles is
`cycle_length * sqrt(k / n_c)`. Convolving a **broadening, right-skewed**
age distribution with a **one-sided** visibility window moves the centroid of
the *visible* subpopulation, and keeps moving it as the distribution spreads.
The peak of the circulating signal therefore drifts relative to the true cycle
boundary, and that drift reads as a period change.

It is progressive, as that account requires. At `n_c = 96`, `b_shape` 15 the
successive intervals are 44.617, 44.346, 44.241 rather than constant,
shortening as the spread accumulates -- and still drifting at the third
interval, so the figures above are "by the third cycle", not a converged
steady state. Over a longer window the deficit would grow further.

**This is the Greischar & Childs mechanism applied to period rather than
multiplication rate** -- some developmental ages are easier to sample than
others, and the age distribution changes over the infection (see
`references.md`). Their paper is prior art on the same interaction.

#### How much of the real-data shift this accounts for

The peak-interval deficit above is the **period channel alone**. The least
squares version of the same experiment (`_scripts/nc-mechanism.R`) asks the
fuller question: which `cycle_length` at the new `n_c` best reproduces the
whole trajectory, amplitude included.

Cells: `cycle_length` minus 45.012 h, where that value best reproduces a
trajectory generated at `n_c` = 96, in hours. Window 96 h is closest to the
Wockner span of 72-120 h. `free` re-optimises `R` and `log10_total0` too.

| `b_shape` | fit `n_c` | 48 h | 96 h | 144 h | 96 h, free |
|---|---|---|---|---|---|
| 15 | 192 | −1.15 | **−1.01** | −0.83 | −0.99 |
| 15 | 384 | −1.79 | **−1.64** | −1.33 | −1.58 |
| 400 | 192 | −1.81 | **−1.41** | −1.10 | −1.38 |
| 400 | 384 | −3.02 | **−2.35** | −1.82 | −2.33 |

Against real-data shifts of **−2.17 h** (96 → 192) and **−3.20 h**
(96 → 384). `np_wide_both` fitted `b_shape` at 64.45, between the two rows.
So at the relevant window the forward map accounts for roughly **half to two
thirds** of the real shift, not the quarter the period channel alone
suggested. **Correction to the earlier figure in this file's history.**

Re-optimising the nuisances barely changes anything (−1.01 vs −0.99,
−1.64 vs −1.58), so the shift is not an artefact of `R` or `log10_total0`
absorbing something.

#### The pre-registered reading rule for this script did not work

`nc-mechanism.R` was written to separate the channels by window length:
constant across windows meaning the within-cycle channel, growing meaning
the accumulating one. **The shift SHRINKS with window length** in every cell,
which the rule had no reading for.

The explanation is that window length changes two things at once, so the rule
was badly designed. A longer window does accumulate more desynchronisation,
but it also **pins the period far more tightly**: three peaks have to align
instead of one, so least squares simply cannot move `cycle_length` as far.
The second effect dominates. The window trend therefore measures how
identified the period is, not which channel causes the bias, and the
peak-interval test is the clean instrument for that. Recorded rather than
quietly dropped, because the rule was pre-registered.

#### A prediction for the retry, recorded before it runs

The forward-map effect is **larger at `b_shape` 400 than at 15** in **every
one of the twelve cells** of the least-squares table, by 0.4 to 1.2 h, as
well as in the peak-interval test (0.98 h vs 0.72 h) where it survives
sub-grid refinement.
So the production configuration should be **more** `n_c`-sensitive, not less:
entries 41-42 should show the effect **persisting and somewhat larger**, not
vanishing. If the retry converges and the effect is gone, this account is
wrong.

#### What this means for `n_c` as a modelling choice

`n_c` is doing two jobs that pull in opposite directions. As a numerical
discretisation it should be as large as affordable, and Channel A does
converge away. As the **Erlang shape it is also the desynchronisation rate**,
and that does **not** converge: at `n_c` = infinity the chain is
deterministic and there is no desynchronisation at all, which is
biologically wrong in the other direction. **Raising `n_c` to fix the
numerics simultaneously removes a biological process**, and the model has no
separate parameter for that rate. That is the real finding, and it is why
elpd keeps improving with `n_c` without that settling what the period is.

### The decay law: the test cannot tell, and it says why

SLURM 29684, 14 tasks, one per `grp_init` unit, submitted 2026-10-07 13:20 and
finished 2026-10-08 10:10. All 14 COMPLETED, every `.err` empty. Elapsed 9:07
to 21:04 per task, against a 6-9 h estimate — the extrapolation of model B at
`n_c` = 384 from model A's cost was **2.3x low**. `_scripts/decay-law-test.R`,
read with `_scripts/decay-law-read.R`, output saved to
`_data/decay-law-read-2026-10-08.txt`, per-unit fits in
`_data/decay-law-unit01..14.rds`.

Maximum likelihood, no priors. Model A is the Erlang chain alone, which gives
synchrony decaying as √(cycles), at `n_c` in {96, 192, 384}. Model B keeps the
chain and adds a lognormal mixture over `cycle_length` across 7 Gauss-Hermite
nodes, which gives a component linear in cycles, at `n_c` in {192, 384}. Both
are scored at 6 parameters, A spending one on its choice of `n_c` exactly as B
spends one on sigma. `d_ll = ll_B - ll_A`, positive favours linear.

**Headline: the design cannot distinguish the two laws.** The pre-registered
rule and a correction to it reach that same verdict by different routes.

**The pre-registered reading.** Total `d_ll` is **+1.53** over 14 units, B
ahead in 8 and A ahead in 6 — mixed signs, which the rule called
uninformative. The total is also **not distributed**: one unit (`OZ439|1800`)
supplies +1.449 of it, so the other 13 together give +0.08.

**The correction, which is a property of the models and not of the result.**
B contains A at sigma = 0, and A's maximum is never at `n_c` = 96 in any unit,
so both models maximise over the same {192, 384} and `ll_B >= ll_A` must hold
in exact arithmetic. **A negative `d_ll` is therefore impossible as evidence**
and can only be the Nelder-Mead optimiser stopping short. The reading rule, as
written, counted six optimiser failures as wins for A. The nesting was stated
at the top of `decay-law-read.R` from the start; the rule simply did not use
it.

The violations instead measure the noise floor. Cells below are `ll_B - ll_A`
at **matched** `n_c`, one per unit x `n_c`, in log-likelihood units; every cell
must be non-negative.

| violations | worst | median violation | noise floor |
|---|---|---|---|
| 6 of 28 | −0.521 (`OZ439/DSM265`, `n_c` 384) | −0.0169 | 0.521 |

Against a floor of 0.521, **1 of 14 units clears it**: `OZ439|1800` at +1.449.
The next largest gain is `Piperaquine|1800` at +0.396, below the floor.

**Power, which is the real finding.** The maximised gain is by definition at
least the gain at any particular sigma, so each unit's `d_ll` is an **upper
bound** on what a linear component of the pre-registered size — sigma = 0.033,
the between-parasite CV in cycle duration that reproduces `n_c` = 192's spread
of 0.158 cycles at k = 4.8 — could have bought. Five units' optimisers landed
within 0.01 of that sigma, and their gains are **+1.449, +0.097, +0.069,
+0.069, and +0.053**. Excluding the one unit above the floor, the bound across
all units is **at most +0.396, median +0.0000**.

So a linear component of exactly the size predicted is nearly free in
likelihood terms over 1130 observations. The data do not prefer √; they are
indifferent between the two laws at the amplitude the project has reason to
expect. That is a statement about this design, not about the biology.

**One caveat on the power bound.** It needs B's optimiser to have found its
maximum, which in 6 of 14 units it demonstrably did not. For those units the
bound is not airtight and the true gain could be larger. `fit_unit()` starts B
from two fixed points and never warm-starts it from A's solution; a third
start at A's optimum with a small sigma would make the violation impossible
and cost one extra optimisation. The present result is therefore **biased
toward A in magnitude**, by up to the floor.

**Consequence for the IPM.** The test was built to gate thread 15, and it does
not close it either way. The branch that would have killed a Gaussian-kernel
IPM — B clearly ahead and consistent in sign — did not fire. Neither did the
branch that would have endorsed one, since √ was not shown adequate, only not
beaten. An IPM therefore **cannot be justified on the grounds that the √ law
is established**, and a 1-D IPM with a Gaussian development kernel would still
reproduce √ by construction. The case for an IPM rests entirely on the
structural argument — that `n_c` is simultaneously the mesh and the rate — and
has to be made on biology.

### By-product: the `n_c` result reproduces without priors

Model A on its own is a likelihood profile over `n_c` with no priors anywhere
in it. Cells are the maximised profile log-likelihood of model A, differenced
across `n_c` within a unit, in log-likelihood units, summed over the 14
`grp_init` units that hold all 1130 observations. Positive favours the higher
`n_c`.

| comparison | summed over 14 units | units favouring the higher `n_c` |
|---|---|---|
| 192 over 96 | **+72.4** | **14 of 14** |
| 384 over 96 | +82.8 | 14 of 14 |
| 384 over 192 | +10.4 | 10 of 14 |

The Bayesian comparison on the same 1130 observations put `n_c` = 96 **71.8
elpd** (se 13.2) behind 192, with 192 and 384 effectively tied. The ordering,
the magnitude and the plateau all reproduce here with no priors, a simplified
deterministic forward model and a different inferential machine.

Two things that agreement is not. It is **not independent confirmation of the
number**: an in-sample maximised likelihood and an out-of-sample elpd
difference are different quantities, and matching to 0.6 units is coincidence.
And the per-unit gains are concentrated — `MMV048_PIB|1800` and
`Piperaquine|1800` contribute +10.2 and +10.8 of the +72.4 while
`MMV048_PartB|2800` and `OZ439|1800` contribute +0.30 and +0.33 — so the
result is a consistent sign across all 14 units, not a uniform effect.

What it does rule out is a prior artefact. Nothing in model A has a prior, so
the conclusion that **`n_c` = 96 is structurally wrong** no longer depends on
any part of the Bayesian specification.

### The exact chain is a convolution, and it runs 38x faster

`_scripts/ipm-prototype.R`, output `_data/ipm-prototype-2026-10-08.txt`
(tracked). The forward maps live in `_scripts/convolution-forward-map.R` and
are sourced rather than duplicated. Design reasoning in
`claude/ipm-decision.md`.

Work in **absolute** developmental age rather than age modulo the cycle, and
the model's structure collapses:

- **Transport** is a pure-birth process. A parasite's absolute stage after
  time `t` is its starting stage plus Poisson(`lambda * t`), with
  `lambda = n_c / cycle_length`. No wrap, no boundary condition.
- **Growth** is a weight, not a flux: a parasite at absolute stage `j` has
  divided `d = floor((j-1)/n_c)` times, so it carries `R^d`.
- **Sequestration** is a weight too, because circulating status is a function
  of position within the cycle and **resets at division** — the ODE sends both
  `N[n_c]` and `S[n_c]` into `N[1]`.

So the whole trajectory is a convolution evaluated once per observation time.

Cells below are the maximum over the 7 observation times of
|convolution − `mat_exp_series`| / `mat_exp_series`, at `cycle_length` 45,
`b_shape` 400, `R` 8, `mu` 0; and seconds for one 7-point trajectory, median
of 5 runs.

| `n_c` | max rel err | `mat_exp_series` | convolution | speedup |
|---|---|---|---|---|
| 96 | 4.1e−12 | 0.015 s | 0.006 s | 2x |
| 192 | 1.0e−12 | 0.154 s | 0.009 s | 14x |
| 384 | 8.3e−13 | 1.179 s | 0.028 s | 38x |

The agreement is the same order as the package's own `matrix_exp`
cross-check. This is a **reimplementation from the ODE's structure, not a
refactor**, so it is a second independent check of the forward map. On a whole
unit fit the gain is larger still — **12 s against 1777 s at `n_c` = 384**,
148x, because the matrix build dominates at high `n_c` — and it reproduces
SLURM 29684's maximised log-likelihoods at `n_c` 96/192/384 to **1e−6**.

**Two subtleties, both found by validating rather than by reading the source,
and either of which a fresh implementation would plausibly get wrong in
silence.** The chain applies the sequestration hazard with a **one-stage lag**
— the transition from stage `k−1` to `k` uses `q[k−1]`, the O(1/`n_c`) delay
the Stan source documents — so the circulating fraction is `G[k−1]` and not
`y[k]`; using `y[k]` costs 12% at the first observation. And parasites in
their **first cycle** have not been reset, so their weight depends on where
they started; it is separable, so it costs one extra convolution. My first two
weightings were both wrong and only the third reached 1e−12.

**What this changes.** The speedup is the *same model*, so it needs no rewrite
and invalidates no fit. Cost is therefore no longer a reason to replace the
age structure, which leaves only the decoupling argument.

### A continuous kernel is not a reparameterisation of the chain

Same script. Cells are the maximum over the 7 times of
|log10(IPM) − log10(chain)|, in log10 units — the scale the likelihood works
on — at mesh `M` = `n_c` with the kernel's variance matched to the chain's.
For scale, the residual sd of these ML fits is about **0.48** log10 units.

| `n_c` | gaussian kernel | gamma kernel |
|---|---|---|
| 96 | 0.0135 | 0.0120 |
| 192 | 0.0331 | 0.0310 |
| 384 | 0.0623 | 0.0595 |

The chain's stage at time `t` is exactly **Poisson, a lattice distribution**,
and any continuous kernel differs from it in the tails.

**The gap grows with `n_c`, and skew is not the cause.** A gamma kernel —
`Gamma(shape = n_eff*t/cl, scale = cl/n_eff)`, which is the chain with
*continuous* `n_c` — preserves the right skew exactly and lands within 0.003
of a gaussian. The cause is the **troughs**: sequestration hides ~60% of the
population, so the observable spans four orders of magnitude across the cycle
(trough/peak = 7e−6 at `n_c` = 384), and at a trough the value is set by the
**tail** of the age distribution. At `n_c` = 384 the entire 0.0595 sits at
**one** time, t = 72 h; the other six agree to 0.006 or better.

Troughs are also where **low-end censoring** is an open question, so this is
the one place where the kernel choice and a known data issue coincide. Any IPM
fit must be compared with the chain **at the troughs**, not on an average.

### The mesh and the dispersion rate do separate

Cells are the maximum over the 7 times of |log10 difference| from the chain at
`n_c` = 384, the dispersion the data prefer.

| forward map | difference |
|---|---|
| chain at `n_c` = 96 — mesh and rate locked together | **1.589** |
| gamma IPM, mesh M = 96, `n_eff` = 384 | 0.064 |
| gamma IPM, mesh M = 192, `n_eff` = 384 | 0.061 |

A factor of 25. A coarse mesh with the fine dispersion dialled in is close to
the fine chain; the coarse chain is not. That is the design argument's one
empirical claim, and it holds.

### SLURM 30576: the profiles are not readable, and why

`_scripts/nc-profile-fast.R`, output `_data/nc-profile-fast-2026-10-08.txt`
(tracked). All 14 tasks COMPLETED, 17–38 min each, 392 fits, 393.5
CPU-minutes. **Every unit's built-in regression check passed at 1e−11 or
better** against SLURM 29684's saved log-likelihoods, so the convolution
forward map is sound and this is not a forward-map problem.

**The profiles are rougher than the rule that reads them.** A profile
likelihood in a smooth parameter is smooth, so deviation from a smooth fit is
a lower bound on the optimiser's error. Cells are the maximum over rungs of
|pooled `ll` − loess fit in log(knot)|, in log-likelihood units, against the
**2-unit** currency the pre-registered rule is written in.

| profile | rungs | max deviation | sd |
|---|---|---|---|
| exact chain | 9 | 2.39 | 1.48 |
| gamma IPM, mesh fixed at 192 | 16 | **7.17** | 3.66 |

The reader now measures this first and returns **NO VERDICT** above 2 units.
Both profiles trip it. **No conclusion about the dispersion is drawn from this
run.**

**Three causes, only one of them noise.**

1. **A bound in the screen, not noise.** `b_shape` hit the 5000 cap in 13–14
   of 14 units at every rung with `n_eff` ≤ 157, and in 2 of 14 up to 382.
   Those rungs report a **lower bound** on `ll`, not a maximum. The finding
   behind it is real and worth keeping: **at coarse dispersion the fit wants
   more initial synchrony than the screen allows.** It does not threaten the
   location of the optimum, which sits ~100 units above that arm in the
   bound-free region.
2. **Optimiser noise, about 2.3 units.** An interior dip at `n_eff` 691.7 and
   930.4 (first differences −1.30, −1.01) that recovers by +5.01 at 1251.4,
   with no unit at the bound.
3. **One unit's optimiser failure, worth 11.29 units.** The pooled profile
   dropped 11.10 units at the finest rung, `n_eff` = 4096, which looked like a
   turnover. `_scripts/profile-noise-check.R` refitted that rung with eight
   starts instead of two: **DSM265|1800 alone gained +11.29**, accounting for
   the entire pooled drop.

**A correction, and the reason it is instructive.** The mesh-convergence check
refitted `n_eff` = 4096 at a mesh twice as fine and shifted it by only 0.50,
which reads as confirmation that the drop was real rather than a
discretisation artefact. It is not. **Both meshes ran the same Nelder-Mead
harness from the same two starts**, so both inherited the same failure.
Agreement between two instances of one optimiser is not independent
confirmation of anything. Recorded in `claude/gotchas.md`.

**What the run does establish**, none of it about the dispersion:

- The convolution harness reproduces `mat_exp_series` through a full
  optimisation, not merely on a single trajectory: 1e−11 or better on all 14
  units at `n_c` 96, 192 and 384.
- The screen's `b_shape` cap of 5000 binds at coarse dispersion.
- Two starts are not enough at fine rungs, where the likelihood surface is
  harder. The cost of eight is affordable now and was not before.

**A bug of mine, fixed.** `NE_CHECK` in `nc-profile-fast.R` was not drawn from
`NE_IPM`, so two of three mesh comparisons paired against nothing and the
reader crashed on `NA`. Only the comparison at `n_eff` = 4096 was valid — and
that is the one the correction above shows was misleading anyway.

### The production ladder converged on a reseed: 41.20 h at `n_c` = 192

SLURM 30525, entry 39 (`np_bs400_nc192`), `WOCKFIT_SEED=20261008`,
`WOCKFIT_SUFFIX=-seed2`, 2026-10-08, 24.2 h. **Max R-hat 1.0432, inside the
1.05 gate — PASS.** `lp__` means per chain are −834.7, −828.3, −836.5, −830.1,
a spread of **8.3** against 29635's ~96, so the stuck chain is gone and
reseeding was the right diagnosis. Minimum `n_eff` 118 over 1323 monitored
quantities.

Cells: posterior mean of `cycle_length` in hours for each of the 13 trials
(`grp_cl`), `n_c` = 192, `b_shape` pinned at 400, `log10_total0` prior sd 1.

| statistic over the 13 trials | value |
|---|---|
| **mean** | **41.203 h** |
| range | 40.762 to 41.982 |
| sd across trials | 0.412 |
| worst per-trial R-hat | 1.033 |

**Caveats that travel with this number.** 233 of 4000 post-warmup draws were
divergent (**5.8%**), above the 4.0% that was previously the worst in the
project, so this is a mean over an imperfectly explored posterior. Treedepth
saturated in only 2% of transitions with a maximum of 10 — which contradicts
the 46–61% recorded from 29635 and means the tuned retry
(`wockner-fit-nc-bs400-retry.sh`) is now even less justified than it was.

**What it does not answer.** The question this run was for — does the `n_c`
effect survive a pinned `b_shape`? — needs **both** rungs, and entry 40
(`n_c` = 384) failed in 29635 and has not been reseeded. One converged rung is
not a ladder.

### The eight-start profile: both profiles SATURATE

SLURM 30641 plus 30823 (six tasks re-run after a check of mine failed them for
improving). `_scripts/nc-profile-fast.R`, read with
`_scripts/nc-profile-fast-read.R`, output
`_data/nc-profile-fast-2026-10-09.txt`. 392 fits, 1524.5 CPU-minutes.

**What the extra six starts bought.** Cells are the sum over the 14 units of
(eight-start `ll` − two-start `ll`) at that rung, which can only be ≥ 0.

| rung | pooled gain | worst single unit |
|---|---|---|
| gamma IPM, `n_eff` 4096 | 11.94 | 11.29 |
| chain, `n_c` 768 | 4.84 | 3.69 |
| gamma IPM, `n_eff` 930 | 4.65 | 3.49 |
| gamma IPM, `n_eff` 692 | 2.55 | 1.16 |
| chain, `n_c` 1024 | 1.79 | 1.32 |

**46.2 log-likelihood units in total, 5 of 26 rungs gaining more than 2.**
Concentrated, not uniform.

**The verdicts.** Both use the rule pre-registered in
`nc-dispersion-profile-read.R` before any of this data existed.

| profile | best | 2-unit interval | top-rung gain | unimodality violation | verdict |
|---|---|---|---|---|---|
| exact chain, 9 rungs 64–1024 | 1024 | {512, 1024} | +1.30 | **0.00** | SATURATES |
| gamma IPM, 16 rungs, mesh fixed at 192 | 2264 | {1683, 4096} | +0.69 | 0.69 | SATURATES |

**SATURATES means the data put an upper bound on dispersion and no lower
bound — they do not exclude zero.** An IPM would return `sigma_d` pressed
against zero with an interval touching it: a boundary estimate, not a rate.
That **argues against building one**, which is the opposite of what the
structural argument alone suggests.

The mesh is adequate: refitting three rungs at mesh 384 shifts the pooled `ll`
by at most **0.89** units, inside the 2-unit currency.

**`cycle_length` across the 2-unit interval**, which is the honest uncertainty
from this assumption alone: **40.93 to 41.15 h (spread 0.22)** for the chain
and 41.07 to 41.64 (spread 0.57) for the gamma IPM. Compare the **3.20 h**
span across the whole original ladder: once `n_c` ≥ 384 the estimate is
stable, and the large span came from rungs now known to be wrong.

### Why nothing can pin the dispersion: the total spread is what the data see

This is the result that ties the thread together. Cells are the fitted
stage-distribution spread in **cycles** at the last observation (k = 4.8),
pooled over the 14 units: `sd_init` from the fitted `b_shape`, `sd_acc` =
√(4.8/rung) from the rung, and `sd_tot` their quadrature sum, since the two
are independent. Rungs where `b_shape` hit the 5000 cap are excluded.

| profile | `sd_acc` varies by | `sd_tot` varies by |
|---|---|---|
| exact chain | **2.0×** | **12.7%** of its mean (0.1280–0.1450) |
| gamma IPM | **4.4×** | 23.8% of its mean (0.1199–0.1508) |

**The fit holds the total spread almost fixed and trades initial against
accumulated.** `b_shape` falls monotonically as the rung rises in **14 of 14
units** — Spearman correlation of log(`b_shape`) with log(`n_c`) is
**−1.00**, median over units.

That single fact explains the whole thread:

- **why the decay law was not learnable** (29684): if the spread at the end of
  the window is nearly fixed regardless of the rung, its *growth law* is
  barely constrained;
- **why the `n_c` profile saturates**: any rung fits, because `b_shape`
  compensates;
- **why an IPM's `sigma_d` would sit at a boundary**: the data constrain the
  sum, not the accumulation term;
- **and why `b_shape` is unidentified**: observations begin 1.6 cycles in, so
  `b_shape` and the accumulation rate enter only through their sum.

At the best rungs the split is roughly 0.11 cycles initial against 0.05–0.08
accumulated — **the observed desynchronisation is mostly inherited from t = 0,
not acquired during the observation window.** Note the free `b_shape` here
settles near **9–15**, against the production pin of **400**, whose initial sd
is six times smaller; the ML screen has no priors and no hierarchy, so this is
a flag rather than a contradiction, but it is a flag.

### The two-start `mat_exp_series` path disagreed, and is explained

`_scripts/nc-dispersion-profile-read.R` on the `mat_exp_series` runs (29684,
30524, 30527; seven rungs, two starts) returns **TURNS OVER** at `n_c` = 512,
contradicting the eight-start SATURATES. The disagreement is entirely
optimiser quality. Cells are pooled `ll` over the 14 units.

| `n_c` | `mat_exp_series`, 2 starts | convolution, 8 starts | difference |
|---|---|---|---|
| 96 | −856.130 | −856.111 | +0.018 |
| 192 | −783.707 | −783.478 | +0.230 |
| 384 | −773.351 | −772.368 | +0.983 |
| 512 | −771.824 | −771.145 | +0.679 |
| 768 | −774.639 | −770.535 | **+4.104** |

Every difference is ≥ 0, as it must be when the forward maps agree to 1e−11
and only the start count differs. The two-start path is **4.10 units short at
`n_c` = 768**, which turns its 512 → 768 step from **+0.61 into −2.81** and
manufactures the turnover. **Believe the eight-start result**, not because it
is newer but because more starts can only raise a maximum.

Separately, `_scripts/nc-768-read.R` on the same two-start data returns
PLATEAU (384 → 768 = −1.3 summed, 8/14 units positive), which agrees with
SATURATES. The two readers differ only because one includes `n_c` = 512 and
the other does not.

### Bound asymmetry is ruled out

Migrated 2026-10-07 from the old `claude/CLAUDE.md`, where it was the only
record. Settled 2026-10-03.

Cells: paired change in the simulated cycle-length bias, in hours, against
the `default` arm on the same simulated data; negative means the arm reduces
the bias.

Moving the `[min_cl, max_cl]` window so the truth sits asymmetrically within
it changes the bias by **+0.033 h** (n = 4) — nothing, and the wrong sign for
the ~1.5 h that needed explaining. The `[30, 60]` arm's **+0.471 h** is the
prior-width confound that was named at submission, not bound distance:
widening the window also widens the implied prior on `cycle_length`.

Together with the earlier bound-geometry result (+0.059 h centred, +0.03 h
moved and widened), **the bounds carry none of the cycle-length bias.**

### Standing caveats on the evidence

Migrated 2026-10-07 from the old `claude/CLAUDE.md` ("Known thin spots"),
where they were the only record. These are properties of the design and the
evidence base, not of any one result, and they still hold.

- **The real-data prior comparisons are n = 1 by construction** — one
  dataset, one fit per setting — so a shift gets a Monte Carlo z, not a
  confidence interval. A `mean_z` says a shift beats MCMC noise and nothing
  more.
- **Divergences are non-zero in every panel fit** and reach 4.0% in
  `pl_wide_both`, so every posterior mean in the prior panel is a mean over
  an imperfectly explored posterior.
- **The pooling offset (−0.509 h) is about 0.4 of the within-fit posterior sd**
  of one trial's `cycle_length` (1.28 h under `no_pool` with both priors
  corrected). It is small relative to the uncertainty it sits inside.
- **The schedule simulation rests on three noise realisations**, and the
  posterior-draw variant on three, usable only after a refit on a second
  seed. Several paired arms have n_pair of 1 or 2 after the convergence gate.
- **No cycle-length number is yet biological.** The work so far bounds the
  priors' share of the bias from below, attributes it to `b_shape`, and shows
  `n_c` moves the estimate by more than the bias being chased. It does not
  yet deliver a defensible period.
