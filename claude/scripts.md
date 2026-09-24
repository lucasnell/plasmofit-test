# Scripts, configs, file naming, and running on the cluster

*Split out of `CLAUDE.md`, which is now an index. See `CLAUDE.md` for
orientation and the current state of play.*

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
| `wockner-fit.sh` | cluster | sbatch wrapper for `wockner-fit.R`; runs with the working directory set to `_data/`, so outputs land beside the other saved fits |
| `wockner-schedule-sim-analyze.R` | cluster | pools the simulation replicates and compares them against the real fit |
| `wockner-schedule-sim-check.R` | cluster | `run_check = 1` cross-check that the simulator and the likelihood share a forward model |
| `wockner-schedule-bias-profile.R` | cluster | maximum-likelihood `cycle_length` at the real schedules, per trial and pooled -- no prior, no hierarchy, no MCMC |
| `wockner-schedule-sim-mode.R` | cluster | posterior mean vs mode of the population `cycle_length` in the saved simulation fits |
| `wockner-ridge.R` | cluster | posterior correlations and prior-vs-posterior in the real fit; tests whether the simulation's weak identification is real |
| `wockner-inoc-prior.R` | cluster | compares fitted `log10_total0` against `log10(inoculum / 5000 mL)`; sizes the offset the anchored prior has to absorb |
| `wockner-anchor-regression.R` | cluster | checks that `total0_anchor = 0` targets the same posterior as the pre-anchor code |
| `wockner-sim-information.R` | cluster | compares the information the simulated and real designs carry about the oscillation, post-hoc |
| `wockner-anchor-check.R` | cluster | reads the anchored fit against its predictions, against the unanchored fit, and against the widened-prior control |
| `wockner-prior-panel.R` | cluster | reads threads 2 and 3 off the eight-fit real-data panel: sampler health, posterior means, `loo`, and the paired prior and pooling contrasts |

The cluster workflow is in the header comment of `wockner-fit.R` (and
`wockner-fit-kfold.R`, which follows the same pattern): `scp` the script and
CSV up, `sbatch`, `scp` the `.rds` files back into `_data/`. Host alias
`biohpc` (`cbsugreischar.biohpc.cornell.edu`, user `lan68`).

`wockner-fit.R` is driven by a `CONFIGS` list; each entry gives a `model`
(passed to `archer_fit()`) and `data` overrides (passed to
`archer_stan_data()`), and the array range in the sbatch script has to match
its length. It currently fits `no_pool` vs `pooled_cl` to test the
cycle-length hierarchy — see `findings.md`, "Cycle-length hierarchy:
no_pool vs pooled_cl". Earlier runs used it to establish `max_cl = 50` /
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
  settled, but leaning "no".** See `findings.md`, "Cycle-length hierarchy:
  no_pool vs pooled_cl" for the full comparison. Short version:
  observation-level `loo` can't distinguish `no_pool` from `pooled_cl`
  (`elpd_diff` -1.4 ± 1.0), and `sigma_logit_cl` sits at essentially its
  prior scale (~0.49) across every leave-one-trial-out refit, not just the
  full-data fit -- so, as before, the data cannot rule out zero
  between-trial variation, and now that looks robust rather than a one-off.
  The one comparison sharp enough to settle it (trial-level `loo`) is
  uninterpretable: Pareto k > 0.7 for all 13 trial-units, in both models.
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
| `wock-schedsim-{fit,RES}-default-rep2-draw<i>.rds` | same, but simulated from posterior **draw** `i` of `wock-fit-pooled_cl.rds` rather than the posterior mean vector. These carry `truth_source` and `truth_nuisance` in the summary; mean-vector replicates do not, and the analyze script falls back to the `pooled_cl` means for those |
| `regress-ref-{fit,RES}-default-rep2.rds` | frozen copy of `default-rep2` from before the inoculum-anchor package change, kept as the regression reference |
| `wock-ridge.rds` | `wockner-ridge.R` output |
| `wock-inoc-prior.rds` | `wockner-inoc-prior.R` output |

**None of the saved fits contain `log_lik`** — all of them predate it. Any
`loo` work needs fresh fits with `calc_log_lik = 1L`.

Note also that `lp__` is not comparable across configs that change a prior or a
bound (`sd_bs_cl`, `max_cl`, `center_cl` all do), because the prior density term
differs. Compare on the constrained scale instead.

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
