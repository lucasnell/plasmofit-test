# Leave-one-trial-out stability check for the cycle-length (and R) hierarchy,
# follow-up to wockner-fit.R's no_pool vs pooled_cl comparison.
#
# This is NOT full predictive K-fold cross-validation. Scoring a held-out
# trial's data under a posterior that never saw it needs a way to predict a
# brand-new hierarchy group, and archer_fit.stan has no machinery for that
# (see claude/CLAUDE.md and the discussion that produced this script). What
# this does instead: refit no_pool and pooled_cl on each of the 13
# leave-one-trial-out subsets of the Wockner data, and see how much the
# population-level cycle-length (and R) estimates move around depending on
# which trial is excluded. If they're stable across all 13 refits, that's
# evidence the full-data finding (sigma_logit_cl's posterior lower tail is
# essentially its prior's -- see CLAUDE.md) isn't an artifact of any single
# influential trial. It is a robustness/sensitivity check, not a predictive
# comparison -- do not read an elpd or a "winner" out of this script.
#
# One array task per (model, held-out trial) combination: 2 models x 13
# trials = 26 tasks. Each fit uses the same chains/warmup/iter as
# wockner-fit.R by default (see KFOLD_WARMUP/KFOLD_ITER below) -- lower them
# if 26 fits is too slow for your queue; the quantities read off here are a
# handful of population-level scalars, which likely converge with less
# warmup than the full per-series fit needs, but that hasn't been checked.
#
# Only a lightweight summary is saved per task (not the full stanfit): with
# 26 tasks, saving every raw fit would multiply the storage cost of the main
# run for information this script only reads a few scalars out of.
#
#
# First on local machine:
#
# cd ~/GitHub/Cornell/plasmofit-test
# scp ./_scripts/wockner-fit-kfold.R ./_data/wockner-cleaned.csv lan68@cbsugreischar.biohpc.cornell.edu:/home2/lan68/plasmofit/wock-fit-kfold/
#
# Now on cluster:
#
# cd /home2/lan68/plasmofit/wock-fit-kfold
#
# cat << EOF > wockner-fit-kfold.sh
# #!/bin/bash -l
#
# #SBATCH --nodes=1
# #SBATCH --ntasks=1
# #SBATCH --array=1-26
# #SBATCH --cpus-per-task=4
# #SBATCH --mem=8G
# #SBATCH --time=2-00:00:00
# #SBATCH --job-name=wock-kfold
# #SBATCH --output=wock-kfold-%a.out
# #SBATCH --error=wock-kfold-%a.err
# #SBATCH --mail-user=lan68@cornell.edu
# #SBATCH --mail-type=END,FAIL
#
# Rscript --vanilla wockner-fit-kfold.R
#
# EOF
#
# sbatch wockner-fit-kfold.sh
#
# Then when back on local computer:
# cd ~/GitHub/Cornell/plasmofit-test
# scp lan68@cbsugreischar.biohpc.cornell.edu:/home2/lan68/plasmofit/wock-fit-kfold/wock-fit-kfold-*.rds ./_data/
#
# Then analyze with wockner-kfold-analyze.R.
#


.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")

suppressPackageStartupMessages({
    library(tidyverse)
    library(rstan)
    library(posterior)
    library(plasmofit)
})


# ------------------------------------------------------------------------ #
# Fold / model grid
# ------------------------------------------------------------------------ #

KFOLD_MODELS <- c("no_pool", "pooled_cl")

# Same tuning as wockner-fit.R's `no_pool` config (default adapt_delta/
# max_treedepth for both models; see claude/CLAUDE.md on why pooled_cl isn't
# preemptively given different tuning here).
SEED <- 538065874
N_CHAINS <- 4L
KFOLD_WARMUP <- 700L
KFOLD_ITER <- 1000L + KFOLD_WARMUP

n_threads <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
if (is.na(n_threads)) stop("n_threads cannot be NA")
threads_per_chain <- 1L
options(mc.cores = min(N_CHAINS, max(1L, n_threads %/% threads_per_chain)))


# ------------------------------------------------------------------------ #
# Data
# ------------------------------------------------------------------------ #

paras_df_full <- read_csv("wockner-cleaned.csv", col_types = "cccdcdd")
trials <- sort(unique(paras_df_full$trial))
if (length(trials) < 3L) stop("need >= 3 trials for leave-one-trial-out (grp_R/grp_cl need >= 2 groups after holding one out)")

## deterministic (model, held-out trial) grid, one row per array task
FOLDS <- expand.grid(trial = trials, model = KFOLD_MODELS,
                     stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)

curr_idx <- suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")))
if (is.na(curr_idx)) stop("SLURM_ARRAY_TASK_ID must be an integer")
if (curr_idx > nrow(FOLDS)) {
    stop("SLURM_ARRAY_TASK_ID ", curr_idx, " exceeds the ", nrow(FOLDS),
         " (model, trial) folds")
}

held_trial <- FOLDS$trial[curr_idx]
model <- FOLDS$model[curr_idx]
## "/" appears in a real trial name (OZ439/DSM265) and can't go in a filename
held_trial_safe <- gsub("[/\\\\]", "_", held_trial)

cat("fold", curr_idx, "of", nrow(FOLDS), "| model:", model,
    "| held out:", held_trial, "\n")

paras_df <- paras_df_full |>
    filter(trial != held_trial) |>
    mutate(inoc = interaction(trial, inoc_size, drop = TRUE) |> paste()) |>
    mutate(obs_error = interaction(trial, cohort, drop = TRUE) |> paste())

d <- archer_stan_data(paras_df,
                      series = "id", time = "time", abundance = "para",
                      grp_init = "inoc", grp_R = "trial", grp_cl = "trial",
                      grp_sd = "obs_error",
                      calc_log_lik = 0L) # not used here; see header

cat(sprintf("  %d series, %d observations, %d trials in training\n",
            d$n_ts, d$n_total_obs, d$n_grp_cl))


# ------------------------------------------------------------------------ #
# Fit
# ------------------------------------------------------------------------ #

f <- archer_fit(d, model = model, chains = N_CHAINS, iter = KFOLD_ITER,
                warmup = KFOLD_WARMUP, seed = SEED,
                threads_per_chain = threads_per_chain)

cat("leapfrog", round(mean(rstan::get_num_leapfrog_per_iteration(f)), 1),
    "| div", sum(rstan::get_divergent_iterations(f)), "\n")


# ------------------------------------------------------------------------ #
# Summary: population-level cycle_length (and R) implied by this fold, plus
# enough diagnostics to tell whether to trust it. See archer_fit_pooled_cl
# .stan / archer_fit.stan for why the parameter names differ by model:
# no_pool has mu_logit_cl/sigma_logit_cl (hierarchy); pooled_cl has a single
# logit_cl. R stays hierarchical (mu_logit_R/sigma_logit_R) in both, since
# only cycle_length is pooled here.
# ------------------------------------------------------------------------ #

pop_cl_par <- if (model == "pooled_cl") "logit_cl" else "mu_logit_cl"
cl_draws <- as.numeric(rstan::extract(f, pop_cl_par)[[1]])
cl_h <- d$min_cl + (d$max_cl - d$min_cl) * inv_logit(cl_draws)

R_draws <- as.numeric(rstan::extract(f, "mu_logit_R")[[1]])
R_pop <- d$max_R * inv_logit(R_draws)

sigma_cl_draws <- if (model == "no_pool") {
    as.numeric(rstan::extract(f, "sigma_logit_cl")[[1]])
} else NULL
sigma_R_draws <- as.numeric(rstan::extract(f, "sigma_logit_R")[[1]])

rhat_pars <- c(pop_cl_par, "mu_logit_R", "sigma_logit_R",
               if (model == "no_pool") "sigma_logit_cl")
rhats <- sapply(rhat_pars, function(p) {
    posterior::rhat(rstan::extract(f, p, permuted = FALSE)[, , 1])
})

res <- list(
    model         = model,
    held_trial    = held_trial,
    n_grp_cl      = d$n_grp_cl,
    n_divergent   = sum(rstan::get_divergent_iterations(f)),
    max_rhat      = max(rhats, na.rm = TRUE),
    cl_pop_mean_h = mean(cl_h),
    cl_pop_q      = quantile(cl_h, c(0.025, 0.975)),
    R_pop_mean    = mean(R_pop),
    R_pop_q       = quantile(R_pop, c(0.025, 0.975)),
    sigma_cl_mean = if (!is.null(sigma_cl_draws)) mean(sigma_cl_draws) else NA_real_,
    sigma_cl_q    = if (!is.null(sigma_cl_draws)) quantile(sigma_cl_draws, c(0.025, 0.975)) else c(NA, NA),
    sigma_R_mean  = mean(sigma_R_draws),
    sigma_R_q     = quantile(sigma_R_draws, c(0.025, 0.975))
)

write_rds(res, sprintf("wock-fit-kfold-%s-%s.rds", model, held_trial_safe))

cat("\ncycle_length (population): ", round(res$cl_pop_mean_h, 2), " h [",
    round(res$cl_pop_q[1], 2), ", ", round(res$cl_pop_q[2], 2), "]\n", sep = "")
if (!is.na(res$sigma_cl_mean)) {
    cat("sigma_logit_cl: ", round(res$sigma_cl_mean, 3), " [",
        round(res$sigma_cl_q[1], 3), ", ", round(res$sigma_cl_q[2], 3), "]\n", sep = "")
}
cat("Finished fold", curr_idx, "(", model, "/", held_trial, ")\n")
