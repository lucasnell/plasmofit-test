# Is the apparent clustering of per-trial cycle_length real?
#
# The no_pool fit's per-trial posterior means sort into what looks like two
# groups -- seven trials at 44.1-45.0 h, a ~1.1 h gap, then six at 46.1-46.6 h.
# With 13 points that could easily be noise. This script asks three questions,
# all post-hoc on the saved full-data fit:
#
#  1. Is the gap larger than the fitted hierarchy itself would produce? The
#     model says eta_cl ~ normal(mu_logit_cl, sigma_logit_cl), which is
#     unimodal. So within each posterior draw, compare the spacing structure
#     of the 13 fitted cycle_length values against 13 fresh values simulated
#     from that same draw's population distribution. Comparing within a draw
#     keeps it apples-to-apples: both are 13 values on the hours scale, and
#     neither is shrunk relative to the other the way posterior means are.
#
#  2. Do the groups survive the per-trial uncertainty, or do the posteriors
#     overlap enough that the split is an artifact of reading point estimates?
#
#  3. Does group membership line up with anything about the trials? The
#     worry: cycle_length is identified from the phase of the oscillation, so
#     a trial observed briefly or sparsely is weakly identified and gets
#     pulled toward the prior -- which is centred at 48 h, i.e. UPWARD, and so
#     would manufacture exactly this pattern. If the high group is the poorly
#     identified trials, the clustering is an artifact of the prior rather
#     than a property of the parasites.


suppressPackageStartupMessages({
    library(tidyverse)
    library(rstan)
    library(plasmofit)
})

f <- read_rds("_data/wock-fit-no_pool.rds")

paras_df <- read_csv("_data/wockner-cleaned.csv", col_types = "cccdcdd") |>
    mutate(inoc = interaction(trial, inoc_size, drop = TRUE) |> paste()) |>
    mutate(obs_error = interaction(trial, cohort, drop = TRUE) |> paste())

d <- archer_stan_data(paras_df,
                      series = "id", time = "time", abundance = "para",
                      grp_init = "inoc", grp_R = "trial", grp_cl = "trial",
                      grp_sd = "obs_error", calc_log_lik = 1L)

trial_names <- attr(d, "levels")[[attr(d, "grp_columns")$grp_cl]]

cl_of <- function(x) d$min_cl + (d$max_cl - d$min_cl) * inv_logit(x)

cl_draws <- rstan::extract(f, "cycle_length")[[1]]        # draws x 13
mu_draws <- as.numeric(rstan::extract(f, "mu_logit_cl")[[1]])
sg_draws <- as.numeric(rstan::extract(f, "sigma_logit_cl")[[1]])

n_draw <- nrow(cl_draws)
n_trial <- ncol(cl_draws)


# ------------------------------------------------------------------------ #
# 1. Gap structure vs what the fitted hierarchy produces
#
# Statistic: largest spacing between consecutive order statistics, divided by
# the sd of the 13 values. Scale-free, and large when the values fall into
# separated clumps rather than spreading out evenly.
# ------------------------------------------------------------------------ #

gap_stat <- function(x) max(diff(sort(x))) / sd(x)

set.seed(20260918)
obs_stat <- apply(cl_draws, 1, gap_stat)
sim_stat <- vapply(seq_len(n_draw), function(s) {
    gap_stat(cl_of(rnorm(n_trial, mu_draws[s], sg_draws[s])))
}, numeric(1))

cat("=== 1. gap statistic: fitted per-trial values vs the fitted hierarchy ===\n")
cat(sprintf("  fitted values,   mean %.3f (95%% %.3f-%.3f)\n",
            mean(obs_stat), quantile(obs_stat, 0.025), quantile(obs_stat, 0.975)))
cat(sprintf("  simulated draws, mean %.3f (95%% %.3f-%.3f)\n",
            mean(sim_stat), quantile(sim_stat, 0.025), quantile(sim_stat, 0.975)))
cat(sprintf("  P(fitted > simulated), paired within draw: %.3f\n", mean(obs_stat > sim_stat)))
cat("  (~0.5 means the spacing looks like an ordinary unimodal hierarchy;\n",
    "   near 1 means the fitted values clump more than the model expects)\n\n", sep = "")


# ------------------------------------------------------------------------ #
# 2. Does the split survive per-trial uncertainty?
# ------------------------------------------------------------------------ #

split_at <- 45.5
above <- rowSums(cl_draws > split_at)

per_trial <- tibble(
    trial = trial_names,
    cl_mean = colMeans(cl_draws),
    cl_sd = apply(cl_draws, 2, sd),
    p_above = colMeans(cl_draws > split_at)
) |> arrange(cl_mean)

cat("=== 2. per-trial posteriors (split at", split_at, "h) ===\n")
per_trial |>
    mutate(across(c(cl_mean, cl_sd, p_above), \(x) round(x, 3))) |>
    print(n = Inf)

cat(sprintf("\n  typical per-trial posterior sd: %.3f h; gap between groups: %.3f h\n",
            mean(per_trial$cl_sd),
            min(per_trial$cl_mean[per_trial$cl_mean > split_at]) -
                max(per_trial$cl_mean[per_trial$cl_mean < split_at])))
cat("  posterior distribution of how many trials exceed the split:\n")
print(round(prop.table(table(above)), 3))
cat("  (a sharp spike on one value means the split is a stable feature;\n",
    "   a broad spread means it is an artifact of reading point estimates)\n\n", sep = "")


# ------------------------------------------------------------------------ #
# 3. Does group membership track identifiability, or anything else?
#
# The prior sits at 48 h, above every trial, so anything that weakens a
# trial's identification drags it upward. n_cycles is the direct measure of
# how much phase information a trial carries.
# ------------------------------------------------------------------------ #

meta <- paras_df |>
    summarise(.by = trial,
              n_series = n_distinct(id),
              n_obs = n(),
              t_min = min(time),
              t_max = max(time),
              n_times = n_distinct(time),
              n_cohort = n_distinct(cohort),
              inoc_sizes = n_distinct(inoc_size)) |>
    mutate(span_h = t_max - t_min,
           n_cycles = span_h / 45,
           obs_per_series = n_obs / n_series)

tab <- per_trial |>
    left_join(meta, by = "trial") |>
    mutate(group = if_else(cl_mean > split_at, "high", "low"))

cat("=== 3. group membership vs trial characteristics ===\n")
tab |>
    select(trial, group, cl_mean, cl_sd, n_series, n_obs, span_h, n_cycles,
           obs_per_series, n_cohort) |>
    mutate(across(where(is.numeric), \(x) round(x, 2))) |>
    print(n = Inf)

cat("\n  group means:\n")
tab |>
    summarise(.by = group,
              n_trials = n(),
              cl = round(mean(cl_mean), 2),
              cl_sd = round(mean(cl_sd), 3),
              n_obs = round(mean(n_obs), 1),
              span_h = round(mean(span_h), 1),
              n_cycles = round(mean(n_cycles), 2),
              obs_per_series = round(mean(obs_per_series), 2)) |>
    print(width = Inf)

cat("\n  correlation of per-trial cycle_length with:\n")
for (v in c("cl_sd", "n_obs", "n_series", "span_h", "n_cycles", "obs_per_series")) {
    cat(sprintf("    %-16s %+.3f\n", v, cor(tab$cl_mean, tab[[v]])))
}

cat("\n  A strong positive corr with cl_sd, or negative with n_cycles/span_h,\n",
    "  means the high group is the weakly identified trials drifting toward the\n",
    "  48 h prior -- an artifact. No such pattern leaves the clustering as a\n",
    "  real feature of the data, worth explaining on its own terms.\n", sep = "")
