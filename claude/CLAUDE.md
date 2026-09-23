# plasmofit-test — scripts and fits

This repo holds the driver scripts and saved fits for the `plasmofit` package;
it is not itself a package. The package lives at `~/GitHub/Cornell/plasmofit`,
and **the notes for package-side work are separate**, in
`~/GitHub/Cornell/plasmofit/claude/`. Read `claude/chat3/CLAUDE.md` there for
the model internals: the gradient-cost optimizations, the Erlang-window guard,
and the `log_lik` / `reduce_sum` details. This file covers only the scripts and
the fits, and the modelling findings needed to make sense of them.

Note the history: these scripts used to live in `plasmofit/_testing/` and were
moved here. Older commits and some stale comments still refer to that path.

## The scripts

Four, after consolidating seven (commit `0538cd2`), plus two added for the
cycle-length hierarchy comparison and four for the schedule-bias simulation
(both below).

| script | runs where | does what |
|---|---|---|
| `wockner-read.R` | local | builds `_data/wockner-cleaned.csv` from the Wockner supplementary `.docx`/`.xlsx` (needs the files in `~/Box/`) |
| `wockner-fit.R` | cluster | fits the Wockner data, one config per SLURM array task |
| `wockner-fit-analyze.R` | local | reads the saved fits and reports diagnostics and `loo` |
| `wockner-fit-kfold.R` | cluster | leave-one-trial-out stability check, one (model, trial) fold per SLURM array task |
| `wockner-kfold-analyze.R` | local | reads the saved kfold summaries and reports fold-to-fold stability |
| `wockner-pooling-offset.R` | local | decomposes the no_pool / pooled_cl cycle-length gap from the saved fits |
| `wockner-cl-clustering.R` | local | tests whether the per-trial cycle_length clustering is real |
| `test-archer-fit.R` | local | simulates from known parameters, fits, checks recovery |
| `wockner-schedule-sim.R` | cluster | simulates all 13 trials from one true `cycle_length` at the real observation times, one (prior arm, replicate) per SLURM array task |
| `wockner-schedule-sim.sh` | cluster | sbatch wrapper for the above, `--array=1-6` |
| `wockner-schedule-sim-analyze.R` | cluster | pools the simulation replicates and compares them against the real fit |
| `wockner-schedule-sim-check.R` | cluster | `run_check = 1` cross-check that the simulator and the likelihood share a forward model |
| `wockner-schedule-bias-profile.R` | cluster | maximum-likelihood `cycle_length` at the real schedules, per trial and pooled -- no prior, no hierarchy, no MCMC |
| `wockner-schedule-sim-mode.R` | cluster | posterior mean vs mode of the population `cycle_length` in the saved simulation fits |

The cluster workflow is in the header comment of `wockner-fit.R` (and
`wockner-fit-kfold.R`, which follows the same pattern): `scp` the script and
CSV up, `sbatch`, `scp` the `.rds` files back into `_data/`. Host alias
`biohpc` (`cbsugreischar.biohpc.cornell.edu`, user `lan68`).

`wockner-fit.R` is driven by a `CONFIGS` list; each entry gives a `model`
(passed to `archer_fit()`) and `data` overrides (passed to
`archer_stan_data()`), and the array range in the sbatch script has to match
its length. It currently fits `no_pool` vs `pooled_cl` to test the
cycle-length hierarchy — see "Cycle-length hierarchy: comparison in
progress" below. Earlier runs used it to establish `max_cl = 50` /
`sd_bs_cl = 0.5` as the settled data overrides (still the default when `data`
is empty); those configs (`wide_max_cl`, `tight_sd_bs`) were removed once
the evidence was in, and are in git history if needed.

### What was removed and why

`wock-fit-2.R`, `wockner-fit-cnc.R` and `cluster-test-archer-fit.R` all
hand-built the Stan data list, which is now `archer_stan_data()`'s job. They
are in git history if needed. `test-archer-fit.R` was kept but largely
rewritten: its simulation front half was sound, while the back half read a fit
from the pre-move `_testing/` path, referenced variables that no longer
existed, and plotted `y_hat_full`, which had been removed from the model
entirely. That was dead code, not working analysis.

## Background needed to read the configs

The defaults in `plasmofit` (`max_cl = 50`, `sd_bs_cl = 0.5`) were not
arbitrary; they were settled by the runs saved in `_data/`.

- **The posterior is multimodal.** With `max_cl = 55` there is a second mode at
  the upper bound (`cycle_length` ~54 h, with `R` correspondingly larger) which
  fits about 16 log units *worse* than the ~45 h mode. Chains split between
  them and R-hat blows up to 1.5–1.7. `max_cl = 50` removes it, and the
  posterior then sits at ~44.8–45.3 h, well clear of the bound, so the bound is
  not doing the work.
- **Because of that, always set a `seed`.** Without one, two runs differ by
  which mode each chain initialized into, not just by Monte Carlo error. The
  runs before `wock-fit-none.rds` did not set one, which made a "did we break
  something" scare much harder to diagnose than it needed to be. Current runs
  use `538065874`.
- **There is a funnel in the cycle-length hierarchy**, and the tight
  `sd_bs_cl = 0.1` prior was partly causing it by forcing `sigma_logit_cl`
  toward zero, into the neck. Widening to 0.5 pulled it out: divergences
  dropped and the `div_by_sigma` gap narrowed from ~11x to ~1.6x.
- **`center_cl = 0` (non-centred) made things worse**, not better (R-hat 1.74
  vs 1.55). The funnel diagnosis was right, but the dominant problem is
  separated modes in the *global* level `mu_logit_cl`, which non-centring does
  nothing about. Keep `center_cl = 1`.
- **Whether the cycle-length hierarchy earns its keep: tested, still not
  settled, but leaning "no".** See "Cycle-length hierarchy: no_pool vs
  pooled_cl" below for the full comparison. Short version: observation-level
  `loo` can't distinguish `no_pool` from `pooled_cl` (`elpd_diff` -1.4 ± 1.0),
  and `sigma_logit_cl` sits at essentially its prior scale (~0.49) across
  every leave-one-trial-out refit, not just the full-data fit -- so, as
  before, the data cannot rule out zero between-trial variation, and now that
  looks robust rather than a one-off. The one comparison sharp enough to
  settle it (trial-level `loo`) is uninterpretable: Pareto k > 0.7 for all 13
  trial-units, in both models.
- **Threading buys ~nothing here.** Measured ~1.0x: there are only ~14 distinct
  trajectory combinations to split over. `wockner-fit.R` sets
  `threads_per_chain = 1` deliberately; cores are better spent on chains. The
  old `n_threads %/% 4` setting was spending 3 cores per chain for nothing.

## `_data/` file naming

Three naming generations are mixed in there, which is worth knowing before
comparing anything.

| pattern | what it is |
|---|---|
| `fit3.rds` … `fit8.rds` | oldest, provenance unclear |
| `wock-fit-400.rds`, `wock-fit-700.rds` | pre-optimization runs, no leading zero |
| `wock-fit-0400/0700/1000.rds` | warmup-length comparison, `%04i` naming, no fixed seed |
| `wock-fit-<config>.rds` | current: config-named, seed `538065874` |
| `wock-fit-RES-<config>.rds` | `summarize_fit()` output (diagnostics, needs a live DSO to produce) |
| `wock-fit-LOO-<config>.rds` | `loo` objects, only from runs with `calc_log_lik = 1` |
| `wock-data-<config>.rds` | the Stan data list, written alongside the fit |
| `wock-schedsim-fit-<arm>-rep<n>.rds` | schedule-bias simulation fits |
| `wock-schedsim-RES-<arm>-rep<n>.rds` | their summaries, read by `wockner-schedule-sim-analyze.R` |

**None of the saved fits contain `log_lik`** — all of them predate it. Any
`loo` work needs fresh fits with `calc_log_lik = 1L`.

Note also that `lp__` is not comparable across configs that change a prior or a
bound (`sd_bs_cl`, `max_cl`, `center_cl` all do), because the prior density term
differs. Compare on the constrained scale instead.

## Gotchas that have bitten more than once

- **`.libPaths()` at the top of `wockner-fit.R` points at the cluster library.**
  Locally that path does not exist and the call drops the user library, so
  `library(rstan)` fails. The script is cluster-only as written; neutralize that
  line to source it on a laptop.
- **DSO liveness.** `unconstrain_pars`, `log_prob`, `grad_log_prob` and
  `get_num_upars` need a live compiled model. A `readRDS`'d fit gives "the model
  object is not created or not valid". Workaround is a throwaway fit in the
  current session used as the handle. This is why `summarize_fit()` computes
  `sec_per_grad` and `fit_cond` on the cluster rather than later.
- **`packageVersion("plasmofit")` can lag the source tree.** Do not trust it to
  tell you which Stan code is compiled in.
- **Single-chain runs break array indexing.** `rstan::extract(fit,
  permuted = FALSE)[, , "lp__"]` collapses to a vector when there is one chain,
  so `ncol()` returns `NULL`. Fixed in `summarize_fit()`, but the same trap is
  easy to reintroduce; use `drop = FALSE` or `apply(..., 2, ...)`.
- **Do not lean on `fit_cond`.** It is reported because it is cheap, but it has
  already failed to flag real trouble in this model. R-hat, divergence counts,
  `div_by_sigma` and per-chain parameter means are what actually caught the
  bimodality.
- **Leave-one-trial-out strains `loo`.** Dropping a whole trial perturbs the
  posterior far more than dropping one observation, so Pareto k goes bad (8 of
  13 trials on a trial run). High k there means the approximation failed, not
  that the model is bad; a trustworthy answer needs K-fold refitting, 13 fits.

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

## Running this on the cluster directly

The workflow above assumes editing locally and `scp`-ing up. If instead you
are working *on* `biohpc`, three things change:

- **The fits are already there**, and they are the large artifacts git does
  not carry (`_data/.gitignore` excludes `*.rds`). Baseline and kfold output
  live in `/home2/lan68/plasmofit/wock-fit/` and
  `/home2/lan68/plasmofit/wock-fit-kfold/`. The analysis scripts all read
  `_data/`, so symlink or copy them in rather than editing paths:
  `ln -s /home2/lan68/plasmofit/wock-fit/*.rds _data/`. `wockner-cleaned.csv`
  *is* tracked, so a clone has it.
- **`.libPaths()` at the top of `wockner-fit.R` is correct there**, and the
  "neutralize it to source locally" gotcha below stops applying. It is only a
  problem in the other direction.
- **Do the analysis in an interactive job, not on the login node.** Each
  saved fit is ~54 MB compressed and `rstan::extract` on it is memory-hungry;
  `wockner-pooling-offset.R` and `wockner-cl-clustering.R` both load a full
  fit. `srun -N 1 -n 1 -c 4 --mem=8G --pty R --vanilla` (or `Rscript`).

## Open threads

Roughly in priority order.

1. **Find the hierarchy-induced bias.** The largest open question, and it is
   about the estimate itself rather than about the hierarchy earning its
   keep. On data generated from a known `cycle_length`, `no_pool`'s
   population estimate runs +1.6 to +1.7 h high, and elimination puts only
   ~0.15 h of that on the prior and at most ~0.5 h on the likelihood (see
   "Where the bias is, by elimination"). The remaining ~1.1 h is the
   hierarchical structure and has no mechanism attached to it yet. Three
   candidates, in the order they are cheap to test:
   - ~~`sigma_logit_cl` estimates ~0.41 when the truth is 0.~~ **Largely
     answered**: the `no_hier` arm (`pooled_cl`, hierarchy removed) still
     biases +0.93 and +1.97 h on replicates 2 and 3, so the hierarchy
     accounts for only ~0.25 h. The `tight_sigma` cross-check was cancelled
     as redundant and pathologically slow; a reading of the two Stan
     programs replaces it (see the section above). Original plan:
     `wockner-schedule-sim.R` arms `tight_sigma` (`sd_bs_cl = 0.001`, same
     model, hierarchy width removed) and `no_hier` (`pooled_cl`, hierarchy
     removed outright), replicates 2-3 of each, tasks 8-9 and 11-12. Two
     directions because `tight_sigma` pins a centred parameterization at a
     near-zero scale and may sample badly, while `no_hier` is well
     conditioned but changes the model; agreement between them is the point.
   - ~~Posterior mean versus maximum on a bounded parameter.~~ **Done and
     eliminated** by `wockner-schedule-sim-mode.R`: mean and mode differ by
     0.044 h and the posterior is symmetric (skew -0.016).
   - Bound geometry: `[35, 50]` leaves 5 h above a truth of 45 and 10 h
     below. Re-simulate with wider bounds, e.g. `[30, 60]`, and see whether
     the bias tracks the asymmetry. Note `max_cl = 55` reintroduces the
     boundary mode, so this needs the bounds moved, not just widened.
   Until one of these lands, no cycle-length number should be reported as an
   estimate of anything biological.
2. **`hold_out` masking in `plasmofit`, then Design A.** The enabling change
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
3. **Design B, only if A is ambiguous.** True leave-one-trial-out K-fold,
   13 folds x 2-4 models. Note it is structurally near-rigged against the
   hierarchy: for a never-seen trial, `no_pool`'s point prediction collapses
   to the population mean, the same location `pooled_cl` gives, so it can
   only win on calibration. That likely explains why trial-level `loo` put
   `pooled_cl` marginally ahead.
4. **`cl_prior_center` decision.** Decide whether the default should move
   off 48 h. Demoted: the simulation showed the prior carries less of the
   error than thought (weight 0.113), so this mostly does not fix anything.
   Matters for reporting a cycle-length number; mostly cancels for model
   comparison.
5. `test-archer-fit.R` has not been run to completion with a full-length fit;
   it has only been smoke-tested with a short one, where recovery was good
   (23/23 parameters inside their 95% intervals, max |z| 0.76). Worth
   revisiting in light of the schedule-bias finding: that smoke test used a
   denser design (8 observations per series, against 4-8 in Wockner).

### Lower priority

6. **More schedule-simulation replicates.** The correlation question is
   limited by having only three noise realizations, not by compute. Another
   6-9 replicates in the default arm would say whether the schedule-induced
   correlation routinely reaches the real -0.93 or only occasionally brushes
   it. Add seeds to `REP_SEEDS` in `wockner-schedule-sim.R` and widen the
   array; ~2 h wall clock.
7. **Re-run `wide-rep1` with a different fit seed.** It failed at R-hat 1.23,
   11.5% divergences, ESS 13, chains 22 `lp__` units apart, and is excluded
   from all the numbers above. It shares a noise seed with `default-rep1`,
   which was also the shakiest default replicate (R-hat 1.030), so it is
   worth knowing whether that simulated dataset is hard or the wide prior is.
8. **Submit the cycle-length prior sensitivity runs.** `wockner-fit.R`
   entries 3-7, `--array=3-7`. Written but never submitted. Demoted: these
   were to discriminate explanation 1 from 2-3 for the sampling-density
   correlation. The simulation has since killed 1 and found against 2, so
   these would now be confirming on real data that the prior is not the
   story -- worth something, no longer decisive.
