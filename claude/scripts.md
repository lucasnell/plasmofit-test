# Scripts, configs, file naming, and running on the cluster

*Detail file. Orientation and settled decisions are in the repo-root
`CLAUDE.md`; current status in `PROJECT_INDEX.md`; live work in `TODO.md`;
the last session in `handoff.md`. (These four replaced `claude/CLAUDE.md`
as the top layer on 2026-10-07.)*

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
| `wockner-ppc-trend.R` (`.sh`) | cluster | option 1: posterior predictive check on the within-series trend statistic, over the four fits spanning the `b_shape` contrast. Post-hoc on saved fits, no refitting, but slow -- 200 draws x 177 series x 4 fits |
| `wockner-bshape-ladder.R` | cluster | reads the pinned-`b_shape` rungs against the free-`b_shape` fit: sampler health, posterior means with the implied starting-stage age range, and paired `loo` against two references |
| `wockner-schedule-sim-bounds.sh` | cluster | thread 2's bound-geometry arms, `cl_move` and `cl_wide_move`, tasks 36-38 and 41-43 |
| `bshape-upper-range.R` | cluster | converts `b_shape` to a starting-stage age range and checks the real posterior against the `max_shape` bound |
| `schedsim-truthfit-check.R` | cluster | reads the `SCHEDSIM_TRUTH_FIT` replicates one at a time against the same arm's `pooled_cl`-truth replicates, so `b_shape` recovery at two truths can be read side by side under identical priors |
| `wockner-schedule-sim-truth.sh` | cluster | option 2: submits `wockner-schedule-sim.R` with `SCHEDSIM_TRUTH_FIT=pl_wide_both`, so the truth has `b_shape` ~ 65 and `cycle_length` 43.7361 h. Tasks 31-33 |
| `panel-bshape-spread.R` | anywhere | reads `_data/wock-prior-panel.rds` only (87 KB, no fits): posterior sd and CV of `b_shape` against its prior's CV, and the within-fit posterior sd of `cycle_length`; seconds to run |
| `_data/panel-watch.sh` | *deleted* | scaffolding, in history at `071cd34`: waits for a SLURM array to drain, logs exit states and which fits landed, then runs the analysis. Run under `setsid` it survives a cleared session or a dropped SSH -- worth copying when arming a long run |

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

## Reference library

Papers bearing on thread 2, kept here so a search is not repeated.

| work | why it matters |
|---|---|
| Greischar & Childs (2023) *Trends Parasitol* 39(8), doi 10.1016/j.pt.2023.05.006 — PDF at the repo root as `mmcm.pdf`, supplement `mmc1.pdf`, both **untracked**, so a fresh clone will not have them | PMR estimates are biased in SYNCHRONOUS infections, worst when initial median parasite age is offset from sampling by ~12 h, because samples land where parasites are sequestered. Asynchronous infections estimate fine. The two controlling quantities are this model's `b_shape` and `b_offset`. |
| Greischar, Reece, Savill, Mideo (2019) *Trends Parasitol*, "The Challenge of Quantifying Synchrony in Malaria Parasites" | Synchrony is hard to quantify from data of this kind — corroborates `b_shape` being unidentified here rather than badly fitted. |
| Subudhi et al. (2020) *Nat Commun* 11, "Malaria parasites regulate intra-erythrocytic development duration via serpentine receptor 10" | The IDC completes in multiples of 24 h under circadian coordination. Makes `cl_prior_center = 48` a biologically motivated choice rather than an arbitrary one. |
| O'Donnell, Greischar, Reece (2021) *Parasite Immunol* 43 | IDC duration is plastic and shortens when mistimed relative to host rhythms — relevant to whether a single `cycle_length` per trial is the right structure at all. |

## `n_c` work, 2026-10-05/06

- `_scripts/nc-sizing.R` — times probe fits at each `n_c` and reports the
  per-leapfrog cost ratio. **A lower bound, not an estimate**: the leapfrog
  count also rises (see `gotchas.md`).
- `_scripts/wockner-fit-nc.sh` — entries 37-38, the real-data `n_c` ladder
  with `b_shape` estimated. SLURM 29633, done.
- `_scripts/nc-ladder-read.R` — reads that ladder: `b_shape`, `cycle_length`,
  sampler health across `n_c`.
- `_scripts/nc-numerical-check.R` — short probe fits with `run_check = 1`, to
  separate a numerical explanation from a structural one via `max_rel_diff`.
  Answered: structural.
- `_scripts/wockner-fit-nc-bs400.sh` — entries 39-40, the same ladder with
  `b_shape` pinned at 400. SLURM 29635.
- `_scripts/wockner-schedule-sim-nc2x2.sh` — tasks 92-94, 99-101, 106-108,
  the simulation 2x2. SLURM 29637. Arms use `sim_n_c` (generate) and
  `build$n_c` (fit), which `wockner-schedule-sim.R` now decouples.
- `_scripts/nc-2x2-read.R` — reads the 2x2 and prints the two contrasts with
  the rule for reading each.
- `_scripts/total0-bloodvol-scale.R` — sizes the between-individual spread in
  initial density implied by blood-volume variation against the observation
  error it would have to be detected over. Backs thread 14. Reads a saved fit
  and the cleaned csv; no fitting.
- `_scripts/nc-2x2-read.R` — reads the `n_c` 2x2 and prints both
  misspecification contrasts with the rule for reading each.
- `_scripts/nc-bias-correct.R` — combines the 2x2 biases with the real-data
  ladder to ask which hypothesis about the true `n_c` reconciles the rungs.
  Restricted to replicates 1-3, which all four cells share.
- `_scripts/wockner-fit-nc-bs400-retry.sh` — entries 41-42, the retry of the
  non-converged production ladder.
- `_scripts/nc-period-check.R` — deterministic, noiseless: measures the
  observable peak-to-peak period against the `cycle_length` that generated
  it, across `n_c` and `b_shape`. Peaks refined sub-grid. This is the script
  behind the mechanism for thread 5.
- `_scripts/nc-mechanism.R` — the same question by least squares: which
  `cycle_length` at one `n_c` best reproduces a trajectory generated at
  another. Also computes the sequestration duty cycle that rules out the
  discretisation channel. Slow (24 optimisations).
- `_scripts/decay-law-test.R` / `.sh` — thread 15's decay-law screen, one
  SLURM task per `grp_init` unit. Maximum likelihood with the observation sd
  profiled out; Gauss-Hermite nodes computed in-script by Golub-Welsch, so
  no extra package. Note `mat_exp_series` needs STRICTLY INCREASING times,
  so the trajectory is evaluated at sorted unique times and indexed back.
- `_scripts/decay-law-read.R` — aggregates the 14 units and reports the
  matched-parameter comparison, a nesting check (B contains A at sigma = 0,
  so a negative `d_ll` is optimiser failure and its magnitude is the noise
  floor), a power bound on what a linear component at the predicted sigma
  could buy, and model A read as a priors-free likelihood profile over
  `n_c`. Output saved as `_data/decay-law-read-<date>.txt`.
- `_scripts/decay-law-768.sh` — extends model A alone to `n_c` = 768 by
  setting `DECAY_MODE=a768`, which switches `decay-law-test.R` to model A at
  768 only and writes `_data/decay-law-a768-unit*.rds`. The finished screen
  is untouched. The pre-registered reading rule is in the sbatch header.
- `_scripts/nc-768-read.R` — merges the screen and the extension and applies
  that rule, reporting the gain at each rung of the `n_c` profile. Guards on
  exactly 4 model-A rows per unit.
- `_scripts/decay-law-profile.sh` — `DECAY_MODE=profile`, model A at `n_c` =
  128, 256, 512, writing `_data/decay-law-prof-unit*.rds`. Fills the gaps so
  the `n_c` profile shows curvature rather than just an ordering.
- `_scripts/nc-dispersion-profile-read.R` — pools model A across every rung
  on disk and reads it as a profile likelihood in the dispersion magnitude,
  since transit CV is `1/√n_c` and `b_shape` is free at every rung. Reports
  the `b_shape`/dispersion trade-off and applies a pre-registered rule with
  four branches: turns over, saturates, still climbing, flat.
- `_scripts/wockner-fit-nc-bs400-reseed.sh` — reseed-only retry of **entry
  39** (`np_bs400_nc192`), the plain config. Entries 41-42 already carry
  `adapt_delta` 0.95 and `max_treedepth` 12, so they are not reseed-only.
  Sets `WOCKFIT_SUFFIX=-seed2`, without which the rerun would overwrite the
  non-converged fit — see `claude/gotchas.md`.
- `_scripts/decay-law-validate.R` — checks the `cycle_length` mixture really
  produces a linear-in-cycles spread before the test is believed.
- `_scripts/ipm-prototype.R` — prototypes and validates a
  transport-with-dispersion forward map against the Erlang chain. Its main
  by-product is that **the exact chain is a convolution**: in absolute
  developmental age, transport is a pure-birth process, growth is the weight
  `R^divisions`, and sequestration is a weight too because circulating status
  resets at division. Reproduces `mat_exp_series` to 1e-12 and runs 38x faster
  at `n_c` = 384. Also measures where a continuous kernel departs from the
  chain's lattice Poisson, which is the troughs. Output
  `_data/ipm-prototype-<date>.txt`, tracked.
- `_scripts/convolution-forward-map.R` — the chain and its continuous
  relaxations as convolutions in absolute developmental age. **Sourced**, not
  duplicated, by `ipm-prototype.R` and `nc-profile-fast.R`, so there is one
  implementation. `cfm_chain` is the exact chain (integer `n_c`, equals
  `mat_exp_series` to 1e-12); `cfm_gamma` is the chain with **continuous**
  `n_eff` at an independent mesh; `cfm_gauss` drops the skew.
- `_scripts/nc-profile-fast.R` / `.sh` — the dispersion profile, densely: the
  exact chain at nine integer rungs 64–1024, the gamma IPM at a fixed mesh of
  192 with sixteen continuous `n_eff` from 48 to 4096, and a mesh-convergence
  check at M = 384. **Stops** if the chain rungs disagree with SLURM 29684's
  saved log-likelihoods by more than 0.01, so the convolution harness cannot
  drift from `mat_exp_series` unnoticed.
- `_scripts/nc-profile-fast-read.R` — applies the rule pre-registered in
  `nc-dispersion-profile-read.R` to both profiles, reports the mesh check, and
  gives the `cycle_length` spread across the 2-log-likelihood interval, which
  is the honest uncertainty from this assumption alone.
