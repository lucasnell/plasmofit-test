# Are the per-trial cycle_length differences manufactured by the observation
# schedules?
#
# The no_pool fit's per-trial cycle_length posterior means correlate -0.93 with
# how densely each trial's series are sampled (observations per series). Three
# explanations were on the table (see claude/CLAUDE.md, "Per-trial cycle_length
# is almost entirely predicted by sampling density"):
#
#   1. prior pull    -- weakly identified trials drift toward the 48 h prior,
#                       which sits above every trial
#   2. schedule bias -- sparse sampling aliases the oscillation, biasing the
#                       recovered phase and hence cycle length
#   3. confounding   -- sampling density proxies something real about the
#                       studies run that way
#
# This script is the decisive test for 1 and 2 together, and separates them.
# Simulate every trial from a SINGLE true cycle_length, using each series'
# real observation times and the real group structure, then fit no_pool and
# ask whether the -0.93 correlation reappears. There is no between-trial
# variation in the truth, so any correlation that shows up is manufactured by
# the design plus the prior, not by biology.
#
# Two prior arms separate 1 from 2:
#   default (sd_logit_cl = 1) -- as fitted; 1 and 2 both active
#   wide    (sd_logit_cl = 2) -- prior pull largely removed; a correlation
#                               surviving here is 2, not 1
#
# If the correlation vanishes in both arms, the real -0.93 is not a schedule
# artifact and explanation 3 is what is left.
#
# The generative parameters are the pooled_cl fit's posterior means, all of
# them, giving a true cycle_length of 45.012 h shared by every trial. That fit
# is the single-cycle-length description of these data, so its phases, growth
# rates and error scales are mutually consistent with that one period. The
# choice was fixed before any simulation was run and is not tuned to the
# outcome.
#
#
# CLUSTER WORKFLOW
#
# Unlike wockner-fit.R this reads the saved fits, so it runs where they are.
# Working on biohpc directly, from the repo root with _data/ populated (see
# claude/CLAUDE.md, "Running this on the cluster directly"):
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   sbatch _scripts/wockner-schedule-sim.sh
#
# Then: Rscript --vanilla _scripts/wockner-schedule-sim-analyze.R
#
# Each task is one (prior arm, replicate) pair and takes ~2 h, the same cost as
# a real fit. Three replicates per arm: a single simulated dataset gives a
# single correlation, and with 13 trials that number is noisy enough that one
# draw of it could not be told from zero either way.


## Cluster library path. Locally this path does not exist and the call drops
## the user library, so neutralize this line to source the script on a laptop.
.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")

suppressPackageStartupMessages({
    library(plasmofit)
    library(tidyverse)
    library(rstan)
    library(posterior)
})


# ------------------------------------------------------------------------ #
# Configuration
# ------------------------------------------------------------------------ #

N_CHAINS <- 4L
ITER <- 1700L
WARMUP <- 700L

## Fit seed held fixed across replicates, so replicate-to-replicate spread in
## the answer is the data realization and not chain initialization. Same seed
## as every other current fit.
SEED_FIT <- 538065874L

## One seed per replicate, drawn once and written down rather than derived
## from the task id, so a replicate is reproducible on its own.
REP_SEEDS <- c(742160183L, 1908445027L, 355721694L)

## Each arm is a model plus data overrides, following wockner-fit.R's CONFIGS.
## default/wide answer where the bias is not (prior); tight_sigma/no_hier ask
## whether it is the hierarchy, from two directions: tight_sigma keeps the
## model identical and only removes the hierarchy's width, no_hier removes the
## hierarchy itself. tight_sigma is the controlled comparison but pins a
## centred parameterization at a near-zero scale, which may sample badly;
## no_hier is well conditioned but changes the model. Agreement between them
## is the point.
ARMS <- list(
    default     = list(model = "no_pool",   data = list(sd_logit_cl = 1)),
    wide        = list(model = "no_pool",   data = list(sd_logit_cl = 2)),
    tight_sigma = list(model = "no_pool",   data = list(sd_logit_cl = 1,
                                                        sd_bs_cl = 0.001)),
    no_hier     = list(model = "pooled_cl", data = list(sd_logit_cl = 1))
)

CONFIGS <- expand_grid(arm = names(ARMS), rep = seq_along(REP_SEEDS)) |>
    mutate(name = sprintf("%s-rep%d", arm, rep))

options(mc.cores = N_CHAINS)

task <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
if (is.na(task) || task < 1 || task > nrow(CONFIGS)) {
    stop("SLURM_ARRAY_TASK_ID must be 1-", nrow(CONFIGS), ", got ", task)
}
cfg <- CONFIGS[task, ]
cfg_name <- cfg$name
seed_noise <- REP_SEEDS[cfg$rep]

arm <- ARMS[[cfg$arm]]

cat("=== schedule-bias simulation:", cfg_name, "===\n")
cat("  arm:", cfg$arm, "| model", arm$model, "|",
    paste(names(arm$data), unlist(arm$data), sep = " = ", collapse = ", "), "\n")
cat("  replicate:", cfg$rep, "| noise seed", seed_noise, "\n\n")


# ------------------------------------------------------------------------ #
# The real design: observation times, series, and group structure
#
# Built exactly as wockner-fit.R and wockner-cl-clustering.R build it, so the
# integer group codes and the series ordering match the saved fits and
# cycle_length[j] refers to the same trial throughout.
# ------------------------------------------------------------------------ #

paras_df <- read_csv("_data/wockner-cleaned.csv", col_types = "cccdcdd") |>
    mutate(inoc = interaction(trial, inoc_size, drop = TRUE) |> paste()) |>
    mutate(obs_error = interaction(trial, cohort, drop = TRUE) |> paste())

d <- archer_stan_data(paras_df,
                      series = "id", time = "time", abundance = "para",
                      grp_init = "inoc", grp_R = "trial", grp_cl = "trial",
                      grp_sd = "obs_error", calc_log_lik = 0L)

trial_names <- attr(d, "levels")[[attr(d, "grp_columns")$grp_cl]]

cat("design:", d$n_ts, "series,", d$n_total_obs, "observations,",
    d$n_grp_cl, "trials\n")


# ------------------------------------------------------------------------ #
# Truth: the pooled_cl fit's posterior means, in full.
#
# Every parameter comes from the same fit, and it has to. The obvious
# alternative -- no_pool's parameters with its cycle_lengths replaced by one
# value -- is not a coherent parameter set: forcing a single cycle length is
# exactly what pushes per-trial timing differences into b_offset (see
# claude/CLAUDE.md, "The no_pool / pooled_cl cycle-length offset"; b_offset is
# the only other parameter that moves between the two fits). Pairing no_pool's
# phases with a pooled period would simulate series whose peaks sit in the
# wrong place, and the fit would then be recovering cycle length from data no
# model generated.
#
# pooled_cl is the single-cycle-length description of these data, fitted as
# such. Simulating from it is the null this script is built to test.
# ------------------------------------------------------------------------ #

## SCHEDSIM_TRUTH_DRAW=<i> takes the truth from posterior draw i instead of
## the posterior mean. The mean vector is NOT a coherent parameter set when
## parameters are correlated -- log10_total0 and R correlate -0.89 in every
## group (see claude/CLAUDE.md, "The ridge is real, but the simulation's weak
## identification is not") -- so a mean-vector truth can sit where no draw
## sits and generate data less informative than the real data. A single draw
## keeps the correlations intact. Output names carry the draw index.
draw_i <- Sys.getenv("SCHEDSIM_TRUTH_DRAW", "")
use_draw <- nzchar(draw_i)
if (use_draw) {
    draw_i <- as.integer(draw_i)
    if (is.na(draw_i) || draw_i < 1) stop("SCHEDSIM_TRUTH_DRAW must be a positive integer")
    cfg_name <- sprintf("%s-draw%d", cfg_name, draw_i)
}

## SCHEDSIM_FIT_SEED=<n> refits the SAME simulated data with a different
## sampler seed. Nothing about the data or the truth changes, so two runs
## that differ only in this seed are two independent samples from one
## posterior -- the null needed to judge whether a package change moved the
## target. Output names carry the seed so the runs sit side by side.
seed_env <- Sys.getenv("SCHEDSIM_FIT_SEED", "")
seed_fit <- SEED_FIT
if (nzchar(seed_env)) {
    seed_fit <- as.integer(seed_env)
    if (is.na(seed_fit)) stop("SCHEDSIM_FIT_SEED must be an integer")
    cfg_name <- sprintf("%s-seed%d", cfg_name, seed_fit)
}

post_mean <- function(f, par) unname(colMeans(as.matrix(f, pars = par)))
post_draw <- function(f, par, i) unname(as.matrix(f, pars = par)[i, ])

f_pl <- read_rds("_data/wock-fit-pooled_cl.rds")
take <- if (use_draw) {
    n_draw <- nrow(as.matrix(f_pl, pars = "cycle_length"))
    if (draw_i > n_draw) stop("SCHEDSIM_TRUTH_DRAW exceeds ", n_draw, " draws")
    \(par) post_draw(f_pl, par, draw_i)
} else {
    \(par) post_mean(f_pl, par)
}

truth <- list(
    b_shape      = take("b_shape"),
    b_offset     = take("b_offset"),
    log10_total0 = take("log10_total0"),
    R            = take("R"),
    sd_iRBC      = take("sd_iRBC")
)
## pooled_cl reports cycle_length as n_grp_cl identical copies of one value
cl_pooled <- take("cycle_length")
stopifnot(length(unique(round(cl_pooled, 10))) == 1L)
truth$cycle_length <- cl_pooled[1]
rm(f_pl); invisible(gc())

cat("truth source:",
    if (use_draw) sprintf("posterior draw %d", draw_i) else "posterior mean",
    "\n")

## A silent length mismatch here would look exactly like a simulation failure
stopifnot(length(truth$b_shape) == d$n_grp_init,
          length(truth$b_offset) == d$n_grp_init,
          length(truth$log10_total0) == d$n_grp_init,
          length(truth$R) == d$n_grp_R,
          length(truth$sd_iRBC) == d$n_grp_sd,
          truth$cycle_length > d$min_cl, truth$cycle_length < d$max_cl)

cat(sprintf("true cycle_length: %.4f h, identical for all %d trials\n\n",
            truth$cycle_length, d$n_grp_cl))


# ------------------------------------------------------------------------ #
# Simulate
#
# Uses the matrix_exp path rather than the series solution the fit inverts, so
# the simulation is not generated by the same code being tested. Observation
# error matches the likelihood: normal on log10(y + 1) with sd_iRBC[grp_sd].
# ------------------------------------------------------------------------ #

ends <- cumsum(d$n_obs)
starts <- ends - d$n_obs + 1L

y_hat <- numeric(d$n_total_obs)
for (i in seq_len(d$n_ts)) {
    ix <- starts[i]:ends[i]
    y0 <- plasmofit:::generate_starts(truth$cycle_length, d$n_c,
                                      truth$b_shape[d$grp_init[i]],
                                      truth$b_offset[d$grp_init[i]],
                                      truth$log10_total0[d$grp_init[i]])
    y_hat[ix] <- plasmofit:::mat_exp_series(y0, truth$cycle_length, d$n_c,
                                            truth$R[d$grp_R[i]], d$mu,
                                            d$ts[ix], d$dt_full)
}
stopifnot(all(is.finite(y_hat)), all(y_hat >= 0))

set.seed(seed_noise)
sd_obs <- truth$sd_iRBC[rep(d$grp_sd, d$n_obs)]
y_sim <- pmax(0, 10^(log10(y_hat + 1) + rnorm(d$n_total_obs, 0, sd_obs)) - 1)

## Does this parameter set actually describe the real data? If the noiseless
## trajectories did not track the real observations, the simulation would be
## testing a design that never produced the -0.93 in the first place.
cat(sprintf("noiseless trajectories vs real observations, log10(y + 1): r = %.3f\n",
            cor(log10(y_hat + 1), log10(d$y + 1))))

cat("simulated vs real, log10(y + 1):\n")
cat(sprintf("  real: mean %.3f sd %.3f range %.2f-%.2f\n",
            mean(log10(d$y + 1)), sd(log10(d$y + 1)),
            min(log10(d$y + 1)), max(log10(d$y + 1))))
cat(sprintf("  sim : mean %.3f sd %.3f range %.2f-%.2f\n\n",
            mean(log10(y_sim + 1)), sd(log10(y_sim + 1)),
            min(log10(y_sim + 1)), max(log10(y_sim + 1))))

d_sim <- d
d_sim$y <- y_sim
for (nm in names(arm$data)) d_sim[[nm]] <- arm$data[[nm]]


# ------------------------------------------------------------------------ #
# Fit
# ------------------------------------------------------------------------ #

## SCHEDSIM_REBUILD=1 regenerates the summary from an already-saved fit,
## for when sampling succeeded but the summary below did not. The data and
## truth above are rebuilt deterministically from the same seeds, so the fit
## being read is the fit those inputs produced.
fit_path <- sprintf("_data/wock-schedsim-fit-%s.rds", cfg_name)

if (identical(Sys.getenv("SCHEDSIM_REBUILD"), "1")) {
    if (!file.exists(fit_path)) stop("no saved fit at ", fit_path)
    cat("rebuilding summary from", fit_path, "-- not refitting\n")
    f <- read_rds(fit_path)
} else {
    f <- archer_fit(d_sim, model = arm$model, chains = N_CHAINS, iter = ITER,
                    warmup = WARMUP, seed = seed_fit, threads_per_chain = 1L)
    write_rds(f, fit_path)
}


# ------------------------------------------------------------------------ #
# Sampler health first: the correlation means nothing if it did not converge
# ------------------------------------------------------------------------ #

has_hier <- arm$model != "pooled_cl"    # pooled_cl has no cycle-length hierarchy
diag_pars <- c("cycle_length", "R", "sd_iRBC")
if (has_hier) diag_pars <- c(diag_pars, "mu_logit_cl", "sigma_logit_cl")
diags <- map(diag_pars, \(p) {
    sims <- rstan::extract(f, p, permuted = FALSE)
    tibble(par = dimnames(sims)$parameters,
           rhat = apply(sims, 3, posterior::rhat),
           ess_bulk = apply(sims, 3, posterior::ess_bulk))
}) |> list_rbind()

lp <- rstan::extract(f, "lp__", permuted = FALSE)
lp_by_chain <- apply(lp, 2, mean)

n_div <- sum(rstan::get_divergent_iterations(f))

cat("=== sampler ===\n")
cat("  divergences:", n_div, sprintf("(%.1f%%)\n",
                                     100 * n_div / (N_CHAINS * (ITER - WARMUP))))
cat("  max R-hat:  ", round(max(diags$rhat, na.rm = TRUE), 4), "\n")
cat("  min ESS:    ", round(min(diags$ess_bulk, na.rm = TRUE)), "\n")
cat("  lp by chain:", paste(round(lp_by_chain, 1), collapse = "  "),
    sprintf("(spread %.2f)\n", diff(range(lp_by_chain))))


# ------------------------------------------------------------------------ #
# The question: does the sampling-density correlation reappear?
# ------------------------------------------------------------------------ #

meta <- paras_df |>
    summarise(.by = trial,
              n_series = n_distinct(id),
              n_obs = n(),
              span_h = max(time) - min(time)) |>
    mutate(obs_per_series = n_obs / n_series)

cl_draws <- rstan::extract(f, "cycle_length")[[1]]          # draws x n_trial
stopifnot(ncol(cl_draws) == length(trial_names))

per_trial <- tibble(trial = trial_names,
                    cl_mean = colMeans(cl_draws),
                    cl_sd = apply(cl_draws, 2, sd)) |>
    left_join(meta, by = "trial")

## Undefined when the model has no between-trial variation to correlate:
## pooled_cl reports n_grp_cl identical copies of one value, so the per-trial
## spread is exactly zero and cor() is 0/0. That is a property of the model,
## not a failure, so the correlation is reported as NA and the recovery
## numbers above carry the arm on their own.
has_spread <- sd(per_trial$cl_mean) > 1e-8

## Headline statistic, computed exactly as wockner-cl-clustering.R computes it
## on the real fit: Pearson correlation of per-trial posterior means against
## observations per series.
r_means <- if (has_spread) {
    cor(per_trial$cl_mean, per_trial$obs_per_series)
} else NA_real_

## The same correlation within each posterior draw. The point-estimate version
## above hides how much of its own sampling uncertainty it has; this does not.
r_draws <- apply(cl_draws, 1, \(x) if (sd(x) > 0)
                     cor(x, per_trial$obs_per_series) else NA_real_)

sig <- if (has_hier) as.numeric(rstan::extract(f, "sigma_logit_cl")[[1]]) else NA_real_

cat("\n=== correlation of per-trial cycle_length with obs_per_series ===\n")
if (has_spread) {
    cat(sprintf("  posterior means:  %+.3f   (real data: -0.93)\n", r_means))
    cat(sprintf("  per draw:         %+.3f  (95%% %+.3f to %+.3f), P(r < 0) = %.3f\n",
                mean(r_draws), quantile(r_draws, 0.025), quantile(r_draws, 0.975),
                mean(r_draws < 0)))
    cat(sprintf("  slope:            %+.3f h per obs/series\n",
                coef(lm(cl_mean ~ obs_per_series, data = per_trial))[2]))
} else {
    cat("  n/a: this model gives every trial the same cycle_length, so there\n")
    cat("  is no between-trial variation to correlate against sampling density.\n")
}
cat(sprintf("  per-trial spread: %.2f-%.2f h (range %.2f), true value %.2f\n",
            min(per_trial$cl_mean), max(per_trial$cl_mean),
            diff(range(per_trial$cl_mean)), truth$cycle_length))
if (has_hier) {
    cat(sprintf("  sigma_logit_cl:   %.3f (95%% %.3f-%.3f), prior scale %g\n",
                mean(sig), quantile(sig, 0.025), quantile(sig, 0.975),
                d_sim$sd_bs_cl))
} else {
    cat("  sigma_logit_cl:   n/a (pooled_cl has no cycle-length hierarchy)\n")
}

cat("\n=== per trial ===\n")
per_trial |>
    mutate(bias = cl_mean - truth$cycle_length) |>
    arrange(obs_per_series) |>
    mutate(across(where(is.numeric), \(x) round(x, 3))) |>
    print(n = Inf)

summ <- list(
    config = cfg_name, arm = cfg$arm, rep = cfg$rep,
    model = arm$model, overrides = arm$data,
    truth_source = if (use_draw) sprintf("draw%d", draw_i) else "mean",
    truth_nuisance = truth[c("b_shape", "b_offset", "log10_total0", "R",
                             "sd_iRBC")],
    sd_logit_cl = arm$data$sd_logit_cl,
    seed_noise = seed_noise, seed_fit = seed_fit,
    true_cl = truth$cycle_length,
    per_trial = per_trial,
    r_means = r_means,
    r_draws = r_draws,
    sigma_logit_cl = if (has_hier)
        c(mean = mean(sig), q025 = unname(quantile(sig, 0.025)),
          q975 = unname(quantile(sig, 0.975)))
    else c(mean = NA_real_, q025 = NA_real_, q975 = NA_real_),
    n_div = n_div,
    max_rhat = max(diags$rhat, na.rm = TRUE),
    min_ess = min(diags$ess_bulk, na.rm = TRUE),
    lp_by_chain = lp_by_chain
)
write_rds(summ, sprintf("_data/wock-schedsim-RES-%s.rds", cfg_name))

cat("\nFinished", cfg_name, "\n")
