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

Four, after consolidating seven (commit `0538cd2`).

| script | runs where | does what |
|---|---|---|
| `wockner-read.R` | local | builds `_data/wockner-cleaned.csv` from the Wockner supplementary `.docx`/`.xlsx` (needs the files in `~/Box/`) |
| `wockner-fit.R` | cluster | fits the Wockner data, one config per SLURM array task |
| `wockner-fit-analyze.R` | local | reads the saved fits and reports diagnostics and `loo` |
| `test-archer-fit.R` | local | simulates from known parameters, fits, checks recovery |

The cluster workflow is in the header comment of `wockner-fit.R`: `scp` the
script and CSV up, `sbatch`, `scp` the `.rds` files back into `_data/`. Host
alias `biohpc` (`cbsugreischar.biohpc.cornell.edu`, user `lan68`).

`wockner-fit.R` is driven by a `CONFIGS` list; each entry is a set of overrides
passed to `archer_stan_data()`, and the array range in the sbatch script has to
match its length. The `wide_max_cl` and `tight_sd_bs` entries exist to keep the
evidence for the current defaults reproducible, not because those settings are
wanted — see below.

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
- **Whether the cycle-length hierarchy earns its keep is still open.** With the
  boundary mode gone and `sd_bs_cl` given room, `sigma_logit_cl`'s posterior
  lower tail is essentially its prior's, so the data cannot rule out zero
  between-trial variation. But they cannot establish it either, and the
  hierarchical version samples best of everything tried, so collapsing to a
  single `cycle_length` would buy a falsely precise estimate. Settling it
  properly is a model comparison, which is what `log_lik` was added for.
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

## Open threads

- Run a real fit with `calc_log_lik = 1L` — nothing saved so far has it.
- Decide the cycle-length hierarchy question, probably via K-fold rather than
  PSIS at the trial grouping.
- `test-archer-fit.R` has not been run to completion with a full-length fit;
  it has only been smoke-tested with a short one, where recovery was good
  (23/23 parameters inside their 95% intervals, max |z| 0.76).
