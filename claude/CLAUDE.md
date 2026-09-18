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
cycle-length hierarchy comparison (see below).

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
Near-deterministic, n = 13. The sparsest-sampled trial (OZ439, 4.8 obs/series)
has the highest estimate (46.6 h); the densest (MMV048_PartB, 7.7) has the
lowest (44.1 h). Related measures are far weaker -- `span_h` and `n_cycles`
-0.40, `n_obs` -0.33, `n_series` -0.03 -- so it is specifically **how densely
each series is sampled**, not how much data the trial has overall. That fits
the mechanism: cycle length is identified from the phase of the oscillation
*within* a series, and adding more sparsely-sampled series doesn't add phase
resolution the way adding observations within a series does.

Three candidate explanations, not yet separated:

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

The prior-sensitivity configs (entries 3-7 in `wockner-fit.R`) now do double
duty: if `np_wide_prior` still shows corr ~-0.93, explanation 1 is dead and
it is 2 or 3. The decisive test for 2 is simulation: generate all 13 trials
from a *single* true cycle length using each trial's real observation times,
fit `no_pool`, and see whether the -0.93 correlation reappears. If it does,
the pattern is manufactured by the schedules. `test-archer-fit.R` already has
the simulation machinery for this.

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

1. **Submit the cycle-length prior sensitivity runs.** `wockner-fit.R`
   entries 3-7, `--array=3-7` (entries 1-2 are already fit). Written but
   never submitted. Answers whether the 48 h prior is pulling the estimates,
   and *also* discriminates explanation 1 from 2-3 for the sampling-density
   correlation: if `np_wide_prior` still shows `corr(cl, obs_per_series)`
   ~ -0.93, prior pull is dead. Re-run `wockner-cl-clustering.R` against
   that fit to check.
2. **The schedule-bias simulation.** The decisive test for whether the
   sampling-density correlation is manufactured by the observation schedules:
   simulate all 13 trials from a *single* true `cycle_length`, using each
   trial's real observation times, fit `no_pool`, and see whether the -0.93
   correlation reappears. If it does, between-trial "variation" is partly a
   design artifact and the earn-its-keep question needs reframing.
   `test-archer-fit.R` has the simulation machinery; `plasmofit:::
   generate_starts()` / `plasmofit:::full_mat_exp_series()` are the
   generators (see `plasmofit`'s `claude/chat4/CLAUDE.md` for the recipe).
   Local/interactive, no array job.
3. **`hold_out` masking in `plasmofit`, then Design A.** The enabling change
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
4. **Design B, only if A is ambiguous.** True leave-one-trial-out K-fold,
   13 folds x 2-4 models. Note it is structurally near-rigged against the
   hierarchy: for a never-seen trial, `no_pool`'s point prediction collapses
   to the population mean, the same location `pooled_cl` gives, so it can
   only win on calibration. That likely explains why trial-level `loo` put
   `pooled_cl` marginally ahead.
5. **`cl_prior_center` decision.** Once 1 lands, decide whether the default
   should move off 48 h. Matters for reporting a cycle-length number; mostly
   cancels for model comparison.
6. `test-archer-fit.R` has not been run to completion with a full-length fit;
   it has only been smoke-tested with a short one, where recovery was good
   (23/23 parameters inside their 95% intervals, max |z| 0.76).
