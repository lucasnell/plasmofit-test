# Gotchas that have bitten more than once

*Detail file. Orientation and settled decisions are in the repo-root
`CLAUDE.md`; current status in `PROJECT_INDEX.md`; live work in `TODO.md`;
the last session in `handoff.md`. (These four replaced `claude/CLAUDE.md`
as the top layer on 2026-10-07.)*

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
- **The saved `wock-fit-LOO-<cfg>.rds` is NOT a bare `loo` object.** It is a
  two-element list, `observation` and `trial`, and handing it straight to
  `loo::loo_compare()` fails with "All inputs should have class 'loo'".
  `wockner-prior-panel.R` has a `pick()` helper for exactly this; copy it
  rather than rediscovering the error. Take `observation` -- trial-level
  `loo` is broken on these data anyway, Pareto k > 0.7 for all 13 units in
  every model.
- **`wockner-fit.R` does not pin the nuisance priors.** It inherits
  `archer_stan_data()`'s defaults, so a change to
  `mean_log_b_shape`/`sd_log_b_shape`/`mean_log10_total0`/`sd_log10_total0`
  in the package would silently redefine what `no_pool` and `pooled_cl` mean
  and break the eight-fit panel's baseline without any error.
  `wockner-schedule-sim.R` and `wockner-ppc-trend.R` *are* pinned; this one
  is not, deliberately, because `np_anchor` passes `inoc_size`, which
  replaces the `log10_total0` prior, and whether `archer_stan_data()` accepts
  both together has never been tested. **Test that before pinning it.**
  Related: the package default `sd_log10_total0 = 0.25` is known to be
  indefensible -- it implies a 647% establishment fraction and costs 104.7
  elpd -- and was left in place by decision, not oversight.
- **Fits from different package builds are not bit-comparable, and this has
  already happened mid-project.** Schedule-simulation replicates 1-3 predate
  the inoculum anchor; replicates 4-5 postdate it, and the anchor added
  `delta_total0`/`sigma_total0` declarations to all four Stan programs. With
  the anchor off those are `array[0]` and contribute nothing, so the model is
  mathematically identical and pooling across noise replicates is fine. Do
  not use replicates from both sides of that change in any bit-exact
  regression test.
- **A new arm in `wockner-schedule-sim.R` does not automatically reach the
  paired tables.** `wockner-schedule-sim-analyze.R` filtered on a literal
  `arm %in% c("default", "wide_bshape", "wide_total0", "wide_nuis")`, inlined
  in two places, so `cl_move` and `cl_wide_move` ran, converged, and appeared
  in the health and recovery sections while being silently absent from the
  thread 2 paired contrast -- no error, just missing rows. It is now the
  named constant `PAIRED_ARMS`, defined once. **Add a new arm there when you
  add it to `ARMS`.**
- **Bound overrides must go through `archer_stan_data()`, not the arm's
  `data` list.** `min_cl`/`max_cl` feed
  `mean_logit_cl = logit((cl_prior_center - min_cl) / (max_cl - min_cl))` and
  the Erlang-window check, so writing them over the already-built data list
  leaves every derived field describing the old bounds and produces a wrong
  fit with no error. `wockner-schedule-sim.R` arms take a `build` list for
  this; `data` is still applied afterwards for true pass-through entries.
- **Moving the cycle-length bounds costs convergence.** Three of six
  bound-geometry replicates failed the R-hat gate, one badly (1.32, 513
  divergences). Tight bounds are worth real compute per the package docs, and
  `max_cl = 55` reintroduced a boundary mode once. Budget for refits on new
  seeds in any bound experiment.
- **`R CMD INSTALL` reuses stale `src/*.o`, so a Stan change can install a
  binary that ignores it.** After editing a `.stan` file, the regenerated
  `src/stanExports_*.h` was correct and `stanmodels[[...]]@model_code`
  contained the new code, but the compiled object files were reused from an
  earlier build and the fitted model silently ignored the new data entries --
  Stan discards unknown data fields without complaint, so there was no error
  anywhere. A fit with `b_shape = 400` returned `b_shape` = 114, not 400.
  **Use `R CMD INSTALL --preclean`, or `rm src/*.o src/*.so` first**, and
  **verify against the installed binary, never the source**: check that the
  new parameter appears in `fit@model_pars` and that a fit actually honours
  the new argument. Checking `@model_code` is not enough -- it was right while
  the binary was wrong.
- **A `00LOCK-plasmofit` directory during a running install is normal.**
  `R CMD INSTALL` creates it at the start and removes it on a clean finish.
  Do not delete it to "clean up" while an install is in flight. It is only a
  problem if it outlives the install, and then the recovery is the one in the
  "Never interrupt `R CMD INSTALL`" entry above: `rmdir` the empty live
  directory, move the lock's copy back, then remove the lock.
- **A FAILED SLURM state does not mean no output.** Thread 1 sat blocked for
  over a week on a null run recorded as FAILED; the job had died at the
  summary stage *after* writing its fit, and the file had been on disk the
  whole time. The fit is always written before the summary, and
  `SCHEDSIM_REBUILD=1` regenerates a missing summary in a minute. **Check
  `_data/` before re-running anything a failed job was supposed to produce.**
- **Editing the data-building step of `wockner-fit.R` changes the group level
  ORDER, which silently invalidates every entry-wise comparison across that
  change.** Adding `group_by(id) |> arrange(time, .by_group = TRUE)` for the
  Design A hold-out mask reordered the rows, and although
  `archer_stan_data()` sorts internally, the factor levels it derives came
  out in a different order. The level SETS are identical, so every fit is
  internally correct and anything averaged over groups is unaffected --
  `cycle_length` agreed to 0.06 h. But `sd_iRBC[6]` then names a different
  (trial, cohort) group in the two fits, and comparing them index by index
  compares unrelated quantities: thread 6's drift read **mean |z| 25.5, max
  272** against a null of 0.85. Realigned by level name it is **0.79**.
  - Check `identical(attr(d, "levels"), attr(d2, "levels"))` before any
    entry-wise comparison of two fits. `_scripts/drift-realign.R` does the
    remapping, and note it must key on the FIRST index only -- `b_off_vec` is
    `array[n] unit_vector[2]`, so keying on block length skips it and leaves
    max |z| at 53.
  - **This happened because the script was edited while a job was reading
    it**, which is the hazard recorded three entries above. The job did not
    die, as that entry warns it might; it completed and produced a fit that
    was correct but not comparable. **A clean exit is not evidence the edit
    was safe.**

## A per-leapfrog cost probe is a lower bound, not an estimate

`_scripts/nc-sizing.R` timed short fits at `n_c = 96` and 192 and reported the
per-leapfrog cost ratio, 1.81x, which it got right. The production fit still
took **2x longer than predicted** (5 h 47 against ~2.9 h) because the leapfrog
count per iteration ALSO rose, 223 -> 390: a harder posterior geometry needs
more gradient evaluations, not just costlier ones.

**Cost = (cost per gradient) x (gradients per iteration) x iterations.** A
probe with a fixed short iteration count measures only the first factor, and
the second moves in the same direction whenever the change makes sampling
harder. Treat such a probe as a lower bound, or run it long enough for the
step size to adapt and read the leapfrog count too.

## `gqs()` needs the parameters block, which the saved fits do not keep

Re-running generated quantities on an existing fit (to get `max_rel_diff`
without refitting) fails with `subscript out of bounds`. `as.matrix(fit)`
returns parameters, transformed parameters, and generated quantities, but
`rstan::gqs()` wants exactly the `parameters` block -- and `archer_fit()` does
not monitor the raw ones (`b_shape_free` is absent; only the transformed
`b_shape` is saved), so they cannot be reconstructed from a saved fit.

If a generated quantity has to be evaluated at the posterior, monitor the raw
parameters at fit time. Otherwise fall back to a short fresh run with the flag
on, and say in the write-up that it is a probe over the sampled region rather
than the posterior.

## A high R-hat can be one stuck chain, and the majority is not the answer

SLURM 29635 returned max R-hat 6.13 and 8.34. The cause was not general
mixing failure: three of four chains agreed closely (`cycle_length[1]` =
42.18, 42.28, 42.10) while one sat 96 log-posterior units below them at
43.09. `lp__` **by chain** showed this immediately where the summary R-hat
did not.

Two things follow. First, **check `lp__` per chain before diagnosing**: a
single number cannot distinguish "all chains lost" from "one chain stuck",
and the fixes differ -- the first wants reparameterisation, the second wants
inits, a seed, or more warmup. Second, **the agreeing majority is not a
result**. Dropping the outlier chain and reporting the other three is
choosing the chains that give a tidy answer, which is exactly the thing this
project does not do. Refit.

Also: when `calc_log_lik` is on, a stuck chain breaks the LOO as badly as the
parameters -- `log_lik` entries here had R-hat near 6, so `loo_compare` on
that fit was meaningless even though it printed cleanly with only a quiet
`k_psis > 0.7` flag.

## `mat_exp_series` needs strictly increasing times, and fails silently inside an optimiser

A unit's rows are many series sharing a handful of sampling times, so passing
the observation times straight in gives a vector that is neither sorted nor
unique. `mat_exp_series` rejects it:

> Error: The array 'ts' is not strictly increasing. Element 7 is 120 while
> element 6 is 192

Inside a `tryCatch` that returns a large sentinel on error, which is the
normal way to write an objective function, **this does not look like an
error**. Every parameter set returns the sentinel, the optimiser reports
convergence, and the fit comes back with an identical implausible
log-likelihood for every model and every `n_c`. The first run of
`_scripts/decay-law-test.R` produced exactly that: five rows of −1e10.

Pass `sort(unique(times))` and index back with `match()`. And when an
objective returns the same value for every configuration, suspect the
objective before the model.

## Trajectory cost is roughly CUBIC in `n_c`

Measured, one `mat_exp_series` call: **0.0148 s at `n_c` = 96, 0.153 at 192,
1.18 at 384, 9.21 at 768** — about 8x per doubling, because the matrix
exponential is over a `2*n_c` square.

This is the budget that governs every `n_c` experiment. It is why the
production fit at `n_c` = 384 took 24.5 h against 1.6 h at 96, why the
decay-law test is a 14-task array rather than one job, and why `n_c` = 768 is
not in any grid. Anything that multiplies the trajectory count — quadrature
nodes, replicates, a profile scan — multiplies that cost on top.

Related: `claude/gotchas.md`, "A per-leapfrog cost probe is a lower bound".
The two compound, since a harder geometry also needs more gradients.

## A killed run can still leave a results file that looks valid

`_scripts/decay-law-test.R` writes one RDS per unit, with five rows, at the
very end. A smoke test run with a reduced `n_c` grid completed its three fits
and wrote a **3-row** file under the same name before being killed. The real
SLURM task for that unit was still running and had not written yet, so for
about an hour a stale file sat in `_data/` that the reader would have treated
as a finished unit -- silently dropping both `n_c` = 384 rungs from that
unit's comparison and changing which model won there.

Two habits, both now built in:

- **Name scratch runs differently from production runs**, or run them in a
  separate directory. A smoke test that writes to the production path is a
  trap even when it is deleted afterwards, because the window between writing
  and deleting is live.
- **Make the reader assert the expected shape**, not just the presence of a
  file. `decay-law-read.R` now stops if any unit does not have exactly five
  rows. A file existing is not evidence that the run that produced it
  finished, or that it was configured the way the current script is.

## A nested model that scores worse is an optimiser failure, not evidence

In the decay-law test (SLURM 29684) model B contains model A at sigma = 0, and
both maximise over the same `n_c` grid in practice, so `ll_B >= ll_A` must hold
in exact arithmetic. Six of 14 units came back with `ll_B < ll_A`, the worst by
0.521 log-likelihood units. The pre-registered reading rule counted those as
wins for A, which they cannot be.

Two lessons, the second more useful than the first.

- **When one model nests another, state the inequality in the reader and
  assert it.** `_scripts/decay-law-read.R` had the nesting written in its own
  header comment from the start and still did not use it. The check is one
  line and it converts a silent misreading into a visible diagnostic.
- **The violations are the measurement.** Their magnitude is the optimiser's
  noise floor on that surface, and nothing smaller than the floor counts as
  signal. That turned an ambiguous 8-to-6 sign split into a clean statement:
  1 of 14 units clears the floor.

The cause here is warm starts, or the lack of them: `fit_unit()` starts the
larger model from two fixed points rather than from the smaller model's
solution. **Start a nested model's optimiser at its special case**, with the
extra parameter at the boundary, and the violation becomes impossible.

## Every cost extrapolation in this project has come in low

The decay-law array was estimated at 6-9 h per task from model A's measured
cost times a factor for model B's quadrature nodes. It ran **9-21 h**, 2.3x
low on the slowest task. The quadrature multiplies the cost per likelihood
evaluation, but adding a parameter also changes the optimiser's path and so
the *number* of evaluations, which the factor does not capture. This is the
same failure mode as the per-leapfrog probe recorded above: a cost model that
holds the iteration count fixed measures only one of the two terms.

**Three misses, all in the same direction**, which is why this is a rule and
not an anecdote.

| what was extrapolated | predicted | actual | factor |
|---|---|---|---|
| production fit, from a per-leapfrog probe | — | — | 2x |
| model B at `n_c` = 384, from model A's cost | 6-9 h | 9-21 h | 2.3x |
| `n_c` = 1024, from the cost at 128-384 | ~25 s | 237 s | ~9x |

The third is the clearest: it is the **same model, same machine, same code**,
extrapolated only across `n_c`, and it still missed by an order of magnitude.
Per-trajectory cost scaled as expected; the **optimiser's iteration count**
did not, because the likelihood surface gets harder at fine rungs. Every one
of these misses has the same shape: a cost model that holds the number of
evaluations fixed measures only one of the two terms.

**Measure the expensive configuration directly on one unit before sizing an
array around it**, and when that is not practical, set the walltime from the
measurement times a factor of three rather than from the extrapolation.

## `wockner-fit.R` output names carry no seed, so a reseed overwrites

The output paths are built from `cfg_name` alone — `wock-fit-<cfg>.rds`,
`wock-data-`, `wock-fit-LOO-`, `wock-fit-RES-`. `WOCKFIT_SEED` changes the
sampler and **nothing in the filename**. So rerunning a config with a new seed
and no suffix silently overwrites the earlier fit, which in a reseed is
precisely the non-converged run that justifies the rerun.

**Always set `WOCKFIT_SUFFIX` when reseeding.** It appends to `cfg_name`
before any path is built, so all four outputs land beside the originals.

A related trap caught at the same time. "Reseed-only run of entry 41" is a
contradiction: entries 39-40 are `np_bs400_nc192`/`_nc384` plain, and **41-42
are the same configs with `adapt_delta` 0.95 and `max_treedepth` 12 already
baked in**. Reseeding without other changes means rerunning **39**, not 41.
Earlier notes said 41 and were wrong. **Check a config's index against its
contents before submitting** — `wockner-fit.R` is 42 entries and the index is
the only thing the array range knows about:

```r
p <- parse("_scripts/wockner-fit.R")
k <- which(vapply(p, function(e) is.call(e) && identical(e[[1]], as.name("<-")) &&
                  identical(e[[2]], as.name("CONFIGS")), logical(1)))
str(eval(p[[k]][[3]])[[39]])
```

## Two runs of the same optimiser agreeing is not independent confirmation

SLURM 30576's gamma-IPM profile dropped 11.10 log-likelihood units at its
finest rung, `n_eff` = 4096. The mesh-convergence check refitted that rung at
a mesh twice as fine and shifted it by only 0.50, which looked like
confirmation that the drop was a real feature of the likelihood rather than a
discretisation artefact. It was not. **Both meshes ran the same Nelder-Mead
harness from the same two starts**, so both inherited the same failure.
`_scripts/profile-noise-check.R` refitted that rung with eight starts and one
unit, DSM265|1800, gained **+11.29** — the whole pooled drop, from one unit.

The mesh check was testing the right thing and could never have found this,
because the thing it varies is not the thing that failed. **A resolution
check only detects discretisation error; it cannot detect optimiser error, and
agreement between two instances of the same optimiser says nothing about
either.** Vary the optimiser as well: more starts, a different method, or a
warm start from a neighbouring rung.

## Read a profile likelihood for SMOOTHNESS before reading its shape

A profile likelihood in a smooth parameter is smooth, so any wiggle is a lower
bound on the optimiser's error. Measure it before applying any rule that works
in log-likelihood units, because under-optimisation can only push `ll` DOWN
and so can both **invent** a turnover (by depressing a neighbouring rung) and
**erase** one (by depressing the peak).

`_scripts/nc-profile-fast-read.R` now computes the maximum deviation from a
loess fit in log(knot) and returns **NO VERDICT** when that exceeds the
2-log-likelihood currency its rule is written in. On 30576 it was 2.39 for the
chain and 7.17 for the gamma IPM, so neither profile was readable.

Three distinct causes were separated, and only one was noise:

1. **A bound, not noise.** `b_shape` hit the screen's 5000 cap in 13-14 of 14
   units at every rung with `n_eff` <= 157. Those rungs report a LOWER BOUND
   on `ll`, not a maximum. At coarse dispersion the fit wants more initial
   synchrony than the screen allows.
2. **Optimiser noise**, about 2.3 units: an interior dip that recovers, with
   no unit at the bound.
3. **One unit's optimiser failure** worth 11.29 units, which the pooled sum
   presented as a feature of the curve.

Note that loess residuals at the two END rungs are unreliable, so judge the
interior by whether the first differences reverse sign where the profile
should be monotone.

## Validate AFTER saving, and make a regression check one-sided when it should be

`_scripts/nc-profile-fast.R` carried a built-in check that its chain rungs
reproduced SLURM 29684's saved log-likelihoods, and stopped if any differed by
more than 0.01. Two things were wrong with it, and together they cost six
tasks of SLURM 30641 about 2.5 hours each.

**It was two-sided when the situation had become one-sided.** The check was
written for a run with the same two optimiser starts as 29684, where the
values should match exactly. The re-run used eight starts, and **more starts
can only raise a maximised log-likelihood**, so a higher value is an
improvement and not drift. Six tasks improved by up to **+0.757** ll units at
a single rung and were failed for it. Only a DROP below the stored value means
the harness has drifted. **Ask which direction a difference can legitimately
go before writing the comparison**, and re-ask it whenever the thing being
compared changes.

**It ran before `saveRDS`.** Every fit had completed; the stop threw all of it
away. **A check that destroys the work it validates is worse than no check.**
Write the results first, then judge them — a bad result on disk can be ignored
or deleted, while a result never written has to be recomputed.

The improvements are themselves a measurement and are now reported rather than
treated as a fault: with eight starts the chain rungs beat the two-start run
by up to 0.757 ll units on one unit, which is a third of the 2.39 roughness
the chain profile showed.

## A loess residual measures curvature, not optimiser noise

The noise gate added to `_scripts/nc-profile-fast-read.R` on 2026-10-08 used
the maximum deviation from a loess fit in log(knot). It was wrong, and
demonstrably so without reference to any result: the eight-start chain profile
is **strictly monotone** — zero sign reversals in its first differences — and
the loess metric still scored it **2.51**, because a span-0.75 quadratic
cannot follow a curve that falls 150 log-likelihood units over nine rungs.

**Use a property the quantity actually has.** A profile likelihood is
**unimodal**, so the diagnostic is how far the sequence departs from
unimodality: how much it ever falls before its maximum, or rises after it. A
clean profile scores zero whatever its curvature. On the same data that scored
2.51 and 7.46 under loess, the unimodality violations are **0.00 and 0.69**,
both inside the 2-unit currency, and the profiles are readable.

The loess number is still printed as a curvature summary, but it gates
nothing. Note this gate was added mid-session and was never pre-registered;
the reading rule it guards — turns over / saturates / still climbing / flat —
is unchanged.

## More starts can only raise a maximum, so use that to arbitrate

When two runs of the same model disagree and differ only in optimiser effort,
there is nothing to adjudicate: the higher log-likelihood is correct. This
settled a real conflict. The two-start `mat_exp_series` profile said the `n_c`
profile TURNS OVER at 512; the eight-start convolution profile said it
SATURATES. The forward maps agree to 1e−11, every rung-wise difference was
≥ 0, and the two-start path was **4.10 units short at `n_c` = 768** — which
alone turned its 512 → 768 step from +0.61 into −2.81 and invented the
turnover.

**Check the sign of every difference before concluding**, because all-positive
is what this explanation predicts and a single negative would refute it.

## `pgrep -f` and `pkill -f` match their own command line

A watcher built on `pgrep -f "R CMD INSTALL"` reported an install as running
for **two hours after it had failed**, because the pattern appears in the
watcher's own arguments, so `pgrep` always found itself. Two `pkill -f`
attempts on the same pattern then killed the calling shell rather than the
target.

**Poll the artefact, not the process list.** For a build, wait on a terminal
line in the log:

```bash
until grep -qE "^\* DONE|^ERROR:" build.log; do sleep 20; done
```

That cannot self-match and it reports the outcome as well as the completion.
If a process really must be matched, break the literal — `"[R] CMD INSTALL"` —
or kill by PID from `ps`.

## `Rscript --vanilla` hides a broken `~/.Renviron`

Four package installs failed in a row with missing-dependency errors while
every probe command succeeded, because `~/.Renviron` set `R_LIBS` to an **empty
directory** and `--vanilla` implies `--no-environ`, so the probes skipped the
file and the installs did not. `.Renviron` is also read *after* the process
environment, so setting `R_LIBS` on the command line did not override it.

**Diagnose the environment with the same flags the failing command uses.** To
override without editing the user's file, point `R_ENVIRON_USER` at a
replacement:

```bash
R_ENVIRON_USER=/path/to/dev.Renviron R CMD INSTALL --preclean .
```

Note `R CMD INSTALL -l DIR` is not a substitute: it drops the user library from
the search path, so dependencies vanish.

## Never reinstall the package while a production job is sampling

A running `rstan` job holds `plasmofit.so` loaded. Replacing it is *usually*
survivable on Linux, since the old inode persists for the running process, but
"usually" is not a trade worth making against a 24-hour fit at the project's
critical path. **Install to a scratch library and validate there**, then
install to the live one once `squeue` is clear:

```bash
mkdir -p /home2/lan68/plasmofit/.Rlib-dev
# symlink every package except plasmofit into it, then:
R_ENVIRON_USER=<file pointing R_LIBS at the scratch lib> R CMD INSTALL --preclean .
```

`_scripts/conv-series-validate.R` takes `PLASMOFIT_LIB` for exactly this.

## Benchmark against what the code actually calls

`conv_series()` was measured at **394x faster than `mat_exp_series`** and that
number was used to argue it would make `n_c` = 768 production fits affordable.
The likelihood does not call `mat_exp_series`. It calls `ew_poly_eval`, the
Erlang-window polynomial series; `mat_exp_series` appears only in generated
quantities as a verification path. Measured in the fitted model, `conv_series`
is **5x SLOWER** at `n_c` = 96 and 2.8x slower at 768.

The same mistake produced a second wrong number: an `n_c` = 768 production fit
was said to cost ~540 h, extrapolated from the matrix exponential's cubic
scaling. Real fits scale about 2.9x per doubling of `n_c` (29635: 8.3 h at
192, 24.6 h at 384), so 768 is on the order of 3 days.

Both were avoidable by reading the model's own likelihood before quoting a
speedup. **Before claiming an optimisation helps, grep the hot path for the
function being replaced**, and benchmark through the interface the production
code uses -- here `grad_log_prob` on the real data, not the exposed function in
isolation. A function-level benchmark also amortises per-call setup over a
whole trajectory, which the fitted model never does: it calls the forward map
once per trajectory combo with a handful of times each.
