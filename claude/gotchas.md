# Gotchas that have bitten more than once

*Split out of `CLAUDE.md`, which is now an index. See `CLAUDE.md` for
orientation and the current state of play.*

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
- **`~/.Renviron` on the cluster overrides `R_LIBS` passed on the command
  line.** It hard-sets `R_LIBS` and `R_LIBS_USER` to `/home/lan68/Rlibrary`,
  while the Stan toolchain (`rstan`, `rstantools`, `StanHeaders`) lives in
  `/home/lan68/R/x86_64-pc-linux-gnu-library/4.6`. `R CMD INSTALL` of
  `plasmofit` therefore fails to find `rstantools` no matter what `R_LIBS` is
  exported. Prefix the command with `R_ENVIRON_USER=/dev/null`; do not edit
  `.Renviron`, other work depends on it.
- **Do not reinstall `plasmofit` while a fit job is running.** `R CMD INSTALL`
  moves the live package directory into `00LOCK-plasmofit` before writing the
  new one, so for a window the path does not exist; a running R process that
  lazy-loads from `R/plasmofit.rdb` in that window dies. Wait for `squeue` to
  clear, or sequence the install behind the job. There is a second reason on
  this project specifically: a reinstall recompiles the Stan models, and the
  anchor regression test needs its null run to differ from the rebuilt fit by
  the sampler seed *alone*. Installing mid-run would confound seed with
  recompilation and invalidate the null.
- **Never interrupt `R CMD INSTALL`.** It moves the live package directory
  aside into `00LOCK-<pkg>` and restores it at the end. Killing it mid-way
  leaves an empty live directory and the only good copy inside the lock.
  Recovery is `rmdir` the empty directory, `mv` the lock's copy back, then
  remove the lock.
- **`tapply()` returns a 1-d array, and `unname()` does not strip `dim`.**
  Anything built that way and handed to Stan arrives as an array rather than
  a vector. Use `as.numeric()`. This bit `anchor_log10_total0`.
- **A `pooled_cl` fit has no between-trial `cycle_length` variation**, so
  `cor(per_trial$cl_mean, obs_per_series)` is 0/0 and returns `NA`, and
  `quantile()` on the result errors. `wockner-schedule-sim.R` guards this
  with an `sd(...) > 1e-8` check; the analyze script carries such replicates
  through the recovery sections and drops them from the correlation ones.
- **Refitting to rebuild a summary is unnecessary.** `SCHEDSIM_REBUILD=1`
  makes `wockner-schedule-sim.R` read the saved fit instead of sampling,
  which turns a 4 h job into a minute.
- **Do not compare a point-estimate correlation against a distribution of
  per-draw correlations.** The first conditions on the per-trial estimates as
  if known; the second is attenuated by their uncertainty. Mixing them
  produced a `P = 0.0000` that meant nothing. The correct comparisons are
  0.16-0.19. The analyze script reports the two separately and says why.
- **Do not index one list by a gate computed from a re-sorted copy of it.**
  `wockner-schedule-sim-analyze.R` built its health table, `arrange()`d it,
  then used `res[ok]` on the unsorted `res`. `list.files()` orders by locale
  collation, where `wide_bshape-rep1` precedes `wide-rep1`, while
  `arrange(arm, rep)` puts it after; 13 of 22 rows were misaligned, so the
  convergence gate **kept** a replicate with R-hat 1.23 and dropped one that
  had converged. Latent until arms with underscores existed, so results from
  before that are unaffected. The script now keeps the health table in `res`
  order, sorts only for printing, and `stopifnot`s the alignment on both
  sides of the filter.
- **Never edit a script while a job is reading it.** `Rscript` reads a source
  file incrementally rather than parsing it all up front, so editing the file
  shifts byte offsets under a running process and it parses garbage. This ate
  two 1.5-hour fits: jobs 28896 and 28897 started at 10:11:01,
  `wockner-schedule-sim.R` was edited at 10:11:49, and both died at the
  summary stage an hour and a half later with `Error: unexpected ')' in
  "... span_h = max(time) - min(time)) |>P)"` -- text that appears nowhere in
  the file. The same hazard applies to bash scripts. Check `squeue` before
  editing anything a running job sources, or copy the script and submit the
  copy. **Recovery is cheap**: the fit is written before the summary, so
  `SCHEDSIM_REBUILD=1` with the same environment rebuilds the summary in a
  minute instead of refitting.
- **`unconstrain_pars()` needs an entry for every DECLARED parameter, even a
  zero-sized one.** Since the inoculum anchor was added, all four Stan
  programs always declare `delta_total0` and `sigma_total0`; with the anchor
  off they are `array[0]`, contribute no draws, and `rstan::extract` returns
  nothing for them -- but omitting them from the list still fails with
  "variable does not exist". `wockner-fit.R`'s `draw_as_list()` adds them at
  length 0 when `fit@model_pars` says the model has them, which also keeps
  fits compiled before the change working. This broke `grad_time()` and
  `fit_cond()` for **every** config, anchored or not, and the posterior
  regression test could not have caught it because that only compares draws.
- **A length-one container comes back from `extract()` without a `dim`.**
  `x[i, ]` on an (iterations x 1) matrix is a bare scalar in R, so Stan reads
  "dims declared=(1); dims found=()". Use `array(x[i, ], d[2])`, which is
  correct for every length and not a special case.
- **Data-list overrides bypass `archer_stan_data()`'s recycling.**
  `wockner-schedule-sim.R` writes each arm's overrides into the already-built
  list. Some entries are per-`grp_init` **vectors** in the Stan data block
  (`sd_log_b_shape`, `sd_log10_total0`, `mean_log_b_shape`,
  `mean_log10_total0`) and some are true scalars (`sd_logit_cl`, `sd_bs_cl`).
  A scalar written over a vector entry gets as far as data initialization and
  dies with "mismatch in number dimensions declared and found in context".
  The script now recycles to the existing entry's length and errors on an
  override that is not in the list at all.
- **Two fits of the same model are never bit-identical across a recompile.**
  Rebuilding the DSO can reorder floating-point operations and HMC turns a
  last-bit difference into a different trajectory within a few leapfrog
  steps. Any regression test on a Stan change has to compare *posteriors*,
  scaled by Monte Carlo error, against a null of two runs of identical code
  with different seeds -- not draws, and not a flat band on posterior sds.
  See `findings.md`, "Regression test for the anchor".
- **Leave-one-trial-out strains `loo`.** Dropping a whole trial perturbs the
  posterior far more than dropping one observation, so Pareto k goes bad (8 of
  13 trials on a trial run). High k there means the approximation failed, not
  that the model is bad; a trustworthy answer needs K-fold refitting, 13 fits.
