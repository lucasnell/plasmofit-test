# Why does no_pool report a higher cycle_length than pooled_cl?
#
# wockner-fit-kfold.R found no_pool's population cycle length ~0.3-0.7 h above
# pooled_cl's single value, in all 13 leave-one-trial-out folds -- small, but
# far too consistent to be Monte Carlo noise. This script decomposes that gap.
# Everything here is post-hoc on the saved full-data fits; nothing refits.
#
# Expects _data/wock-fit-no_pool.rds and _data/wock-fit-pooled_cl.rds.
#
# Four candidate mechanisms, in the order tested below:
#
#  A. The statistic, not the model. For no_pool the reported number was
#     transform(mu_logit_cl) -- a logit-scale hyper-mean pushed through a
#     nonlinear transform. The quantity actually comparable to pooled_cl's
#     single value is the hours-scale average of the per-trial cycle_length.
#     inv_logit is concave where we sit (p ~ 0.67), so by Jensen the trial
#     average is LOWER than transform(mu_logit_cl); part of the "offset" is
#     therefore an artifact of which summary got reported.
#  B. Differential prior leverage. Both models put the same prior --
#     normal(mean_logit_cl, sd_logit_cl), centred at cl_prior_center = 48 h --
#     on their population parameter. But pooled_cl's logit_cl sits directly in
#     the likelihood for every observation, while no_pool's mu_logit_cl never
#     enters the likelihood at all (it reaches the data only through
#     eta_cl ~ normal(mu_logit_cl, sigma_logit_cl)). The same prior therefore
#     has more pull on no_pool, and it pulls UP, toward 48 h.
#  C. Trial weighting vs observation weighting. The hierarchical mean weights
#     trials roughly equally; a single pooled value is an observation-weighted
#     compromise. If larger trials have shorter cycles, pooled_cl lands lower.
#  D. Compensation. Whatever cycle_length gives up under forced pooling has to
#     go somewhere: R, the b_offset phase, or the error scale.


suppressPackageStartupMessages({
    library(tidyverse)
    library(rstan)
    library(plasmofit)
})

f_np <- read_rds("_data/wock-fit-no_pool.rds")
f_pl <- read_rds("_data/wock-fit-pooled_cl.rds")

## The Stan data list was written on the cluster but not copied back; rebuild
## it here. archer_stan_data() is pure R and deterministic, and the inputs are
## the same ones wockner-fit.R used, so this reproduces the list that was fit.
paras_df <- read_csv("_data/wockner-cleaned.csv", col_types = "cccdcdd") |>
    mutate(inoc = interaction(trial, inoc_size, drop = TRUE) |> paste()) |>
    mutate(obs_error = interaction(trial, cohort, drop = TRUE) |> paste())

d <- archer_stan_data(paras_df,
                      series = "id", time = "time", abundance = "para",
                      grp_init = "inoc", grp_R = "trial", grp_cl = "trial",
                      grp_sd = "obs_error", calc_log_lik = 1L)

trial_names <- attr(d, "levels")[[attr(d, "grp_columns")$grp_cl]]

## cycle_length(x) for x on the logit scale, and its first two derivatives,
## which are what the Jensen correction below is built from
cl_of <- function(x) d$min_cl + (d$max_cl - d$min_cl) * inv_logit(x)
cl_d1 <- function(x) { p <- inv_logit(x); (d$max_cl - d$min_cl) * p * (1 - p) }
cl_d2 <- function(x) { p <- inv_logit(x); (d$max_cl - d$min_cl) * p * (1 - p) * (1 - 2 * p) }

mu_cl_np    <- as.numeric(rstan::extract(f_np, "mu_logit_cl")[[1]])
sigma_cl_np <- as.numeric(rstan::extract(f_np, "sigma_logit_cl")[[1]])
cl_np       <- rstan::extract(f_np, "cycle_length")[[1]]   # draws x n_grp_cl
logit_cl_pl <- as.numeric(rstan::extract(f_pl, "logit_cl")[[1]])
cl_pl       <- rstan::extract(f_pl, "cycle_length")[[1]][, 1]

cat("prior on the population parameter (identical in both models):\n")
cat(sprintf("  normal(%.4f, %.2f) on the logit scale = centred at %.2f h\n\n",
            d$mean_logit_cl, d$sd_logit_cl, cl_of(d$mean_logit_cl)))


# ------------------------------------------------------------------------ #
# A. How much of the gap is the choice of summary statistic?
# ------------------------------------------------------------------------ #

from_mu   <- mean(cl_of(mu_cl_np))          # what the kfold script reported
trial_avg <- mean(rowMeans(cl_np))          # hours-scale average across trials
pooled    <- mean(cl_pl)

jensen_pred <- 0.5 * cl_d2(mean(mu_cl_np)) * mean(sigma_cl_np)^2

cat("=== A. summary statistic ===\n")
cat(sprintf("  no_pool, transform(mu_logit_cl):      %.3f h  <- what was reported\n", from_mu))
cat(sprintf("  no_pool, mean over trials of cl[j]:   %.3f h  <- comparable to pooled\n", trial_avg))
cat(sprintf("  pooled_cl, single cycle_length:       %.3f h\n", pooled))
cat(sprintf("  Jensen gap (observed):                %+.3f h\n", trial_avg - from_mu))
cat(sprintf("  Jensen gap (2nd-order prediction):    %+.3f h\n", jensen_pred))
cat(sprintf("  => of the %.3f h reported gap, %.3f h is the statistic, %.3f h is real\n\n",
            from_mu - pooled, from_mu - trial_avg, trial_avg - pooled))


# ------------------------------------------------------------------------ #
# B. Differential prior leverage
#
# Under a normal-normal approximation the posterior mean is a precision
# weighted blend of prior and data, and the weight the prior carries is
# (posterior sd / prior sd)^2. Applying that to each model's population
# parameter says how far each one was pulled toward the 48 h prior centre.
# ------------------------------------------------------------------------ #

prior_pull <- function(draws, label) {
    post_m <- mean(draws); post_s <- sd(draws)
    w <- (post_s / d$sd_logit_cl)^2
    ## post_m = (1-w)*data_loc + w*prior_centre  =>  shift = post_m - data_loc
    data_loc <- (post_m - w * d$mean_logit_cl) / (1 - w)
    tibble(parameter = label,
           post_sd_logit = round(post_s, 4),
           prior_weight = round(w, 4),
           shift_h = round(cl_of(post_m) - cl_of(data_loc), 3))
}

cat("=== B. prior leverage (pull toward the 48 h prior centre) ===\n")
bind_rows(prior_pull(mu_cl_np, "no_pool mu_logit_cl"),
          prior_pull(logit_cl_pl, "pooled_cl logit_cl")) |>
    print(width = Inf)
cat("\n")


# ------------------------------------------------------------------------ #
# C. Trial weighting vs observation weighting
# ------------------------------------------------------------------------ #

per_trial <- tibble(
    trial = trial_names,
    cl = colMeans(cl_np),
    n_series = as.integer(table(factor(d$grp_cl, levels = seq_along(trial_names)))),
    n_obs = as.integer(tapply(d$n_obs, factor(d$grp_cl, levels = seq_along(trial_names)), sum))
) |> arrange(cl)

cat("=== C. per-trial cycle_length (no_pool) vs trial size ===\n")
print(per_trial, n = Inf)

w_none   <- mean(per_trial$cl)
w_series <- weighted.mean(per_trial$cl, per_trial$n_series)
w_obs    <- weighted.mean(per_trial$cl, per_trial$n_obs)

cat(sprintf("\n  unweighted mean of cl[j]:      %.3f h  (diff from pooled %+.3f)\n",
            w_none, w_none - pooled))
cat(sprintf("  series-weighted mean:          %.3f h  (diff from pooled %+.3f)\n",
            w_series, w_series - pooled))
cat(sprintf("  observation-weighted mean:     %.3f h  (diff from pooled %+.3f)\n",
            w_obs, w_obs - pooled))
cat(sprintf("  corr(cl, n_obs) = %+.3f ; corr(cl, n_series) = %+.3f\n\n",
            cor(per_trial$cl, per_trial$n_obs),
            cor(per_trial$cl, per_trial$n_series)))


# ------------------------------------------------------------------------ #
# D. What absorbed the difference?
# ------------------------------------------------------------------------ #

cmp <- function(par, n_expected = NULL) {
    a <- rstan::extract(f_np, par)[[1]]
    b <- rstan::extract(f_pl, par)[[1]]
    am <- if (is.null(dim(a)) || length(dim(a)) == 1L) mean(a) else colMeans(a)
    bm <- if (is.null(dim(b)) || length(dim(b)) == 1L) mean(b) else colMeans(b)
    tibble(parameter = par,
           n = length(am),
           mean_no_pool = round(mean(am), 4),
           mean_pooled_cl = round(mean(bm), 4),
           mean_diff = round(mean(am) - mean(bm), 4),
           max_abs_diff = round(max(abs(am - bm)), 4))
}

cat("=== D. where the compromise went (posterior means, no_pool - pooled_cl) ===\n")
map(c("R", "b_shape", "b_offset", "log10_total0", "sd_iRBC",
      "mu_logit_R", "sigma_logit_R"), cmp) |>
    list_rbind() |>
    print(width = Inf)

cat("\nA large mean_diff in R or b_offset means forced pooling of cycle_length\n",
    "was partly compensated elsewhere rather than simply splitting the\n",
    "difference between trials.\n", sep = "")
