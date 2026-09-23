# Pools the schedule-bias simulation replicates and compares them against the
# real data. Reads _data/wock-schedsim-RES-*.rds (written by
# wockner-schedule-sim.R) and _data/wock-fit-no_pool.rds.
#
# The question, restated: on the real data, per-trial cycle_length posterior
# means correlate -0.93 with observations per series. Each simulated replicate
# has ONE true cycle_length shared by all 13 trials, fitted at the real
# observation times, so any correlation it shows is manufactured by the
# observation schedules and the prior rather than by biology.
#
# Reading the output:
#   correlation near -0.93 in both arms  -> the real pattern is a schedule
#                                           artifact (explanation 2)
#   near -0.93 at sd_logit_cl = 1 only   -> prior pull (explanation 1)
#   near zero in both arms               -> not manufactured; the real -0.93
#                                           needs a different explanation
#                                           (explanation 3)
#
# Runs where the fits are; loads one full fit, so use an interactive job:
#   srun -N 1 -n 1 -c 4 --mem=16G --pty Rscript --vanilla \
#       _scripts/wockner-schedule-sim-analyze.R

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")

suppressPackageStartupMessages({
    library(plasmofit)
    library(tidyverse)
    library(rstan)
})

res_files <- list.files("_data", "^wock-schedsim-RES-.*[.]rds$", full.names = TRUE)
if (length(res_files) == 0) stop("no wock-schedsim-RES-*.rds found in _data/")

res <- map(res_files, read_rds)
cat("read", length(res), "replicate summaries\n\n")


# ------------------------------------------------------------------------ #
# Sampler health first. A replicate that did not converge does not get to
# contribute a correlation.
# ------------------------------------------------------------------------ #

health <- map(res, \(r) tibble(
    config = r$config, arm = r$arm, rep = r$rep,
    div = r$n_div, max_rhat = r$max_rhat, min_ess = r$min_ess,
    lp_spread = diff(range(r$lp_by_chain))
)) |> list_rbind() |> arrange(arm, rep)

cat("=== sampler health ===\n")
health |> mutate(across(where(is.numeric), \(x) round(x, 3))) |> print(n = Inf)

ok <- health$max_rhat < 1.05
if (!all(ok)) {
    cat("\n  NOT CONVERGED:", paste(health$config[!ok], collapse = ", "),
        "-- excluded below\n")
}
cat("\n")

res <- res[ok]
if (length(res) == 0) stop("no converged replicates to analyze")


# ------------------------------------------------------------------------ #
# The real data, computed the same way for comparison
# ------------------------------------------------------------------------ #

paras_df <- read_csv("_data/wockner-cleaned.csv", col_types = "cccdcdd") |>
    mutate(inoc = interaction(trial, inoc_size, drop = TRUE) |> paste()) |>
    mutate(obs_error = interaction(trial, cohort, drop = TRUE) |> paste())

d <- archer_stan_data(paras_df,
                      series = "id", time = "time", abundance = "para",
                      grp_init = "inoc", grp_R = "trial", grp_cl = "trial",
                      grp_sd = "obs_error", calc_log_lik = 0L)
trial_names <- attr(d, "levels")[[attr(d, "grp_columns")$grp_cl]]

meta <- paras_df |>
    summarise(.by = trial, n_series = n_distinct(id), n_obs = n()) |>
    mutate(obs_per_series = n_obs / n_series)
ops <- meta$obs_per_series[match(trial_names, meta$trial)]

f_real <- read_rds("_data/wock-fit-no_pool.rds")
cl_real <- rstan::extract(f_real, "cycle_length")[[1]]
sig_real <- as.numeric(rstan::extract(f_real, "sigma_logit_cl")[[1]])
rm(f_real); invisible(gc())

real <- list(
    r_means = cor(colMeans(cl_real), ops),
    r_draws = apply(cl_real, 1, \(x) cor(x, ops)),
    slope = unname(coef(lm(colMeans(cl_real) ~ ops))[2]),
    spread = diff(range(colMeans(cl_real))),
    sigma = mean(sig_real)
)


# ------------------------------------------------------------------------ #
# Headline comparison
# ------------------------------------------------------------------------ #

tab <- map(res, \(r) {
    pt <- r$per_trial
    tibble(config = r$config, arm = r$arm, rep = r$rep,
           sd_logit_cl = r$sd_logit_cl,
           r_means = r$r_means,
           r_draw_mean = mean(r$r_draws),
           r_draw_lo = unname(quantile(r$r_draws, 0.025)),
           r_draw_hi = unname(quantile(r$r_draws, 0.975)),
           slope = unname(coef(lm(pt$cl_mean ~ pt$obs_per_series))[2]),
           spread = diff(range(pt$cl_mean)),
           bias = mean(pt$cl_mean) - r$true_cl,
           sigma = r$sigma_logit_cl[["mean"]])
}) |> list_rbind() |> arrange(arm, rep)

cat("=== per replicate (true cycle_length identical across trials) ===\n")
tab |> mutate(across(where(is.numeric), \(x) round(x, 3))) |> print(n = Inf)

cat("\n=== by prior arm, vs the real data ===\n")
by_arm <- tab |>
    summarise(.by = c(arm, sd_logit_cl),
              n = n(),
              r_mean = mean(r_means),
              r_min = min(r_means), r_max = max(r_means),
              slope = mean(slope),
              spread = mean(spread),
              sigma = mean(sigma))
by_arm |> mutate(across(where(is.numeric), \(x) round(x, 3))) |> print(width = Inf)

cat(sprintf("\n  real data:  r = %+.3f (per-draw mean %+.3f, 95%% %+.3f to %+.3f)\n",
            real$r_means, mean(real$r_draws),
            quantile(real$r_draws, 0.025), quantile(real$r_draws, 0.975)))
cat(sprintf("              slope %+.3f h per obs/series, per-trial spread %.2f h,\n",
            real$slope, real$spread))
cat(sprintf("              sigma_logit_cl %.3f\n", real$sigma))


# ------------------------------------------------------------------------ #
# Recovery of the known cycle_length, and what the two prior arms say about
# why it misses.
#
# This comes before the correlation because it is the stronger result: the
# truth is known here, so a systematic miss is measurable rather than
# inferred. Replicates 2 and 3 ran in both arms on the same simulated data
# (noise seeds are shared across arms), so the arm difference is the prior and
# nothing else.
# ------------------------------------------------------------------------ #

lg <- function(cl, lo, hi) log((cl - lo) / (hi - cl))
inv_lg <- function(x, lo, hi) lo + (hi - lo) / (1 + exp(-x))

cat("\n=== recovery of the true cycle_length ===\n")
rec <- tab |>
    select(config, arm, rep, sd_logit_cl, bias) |>
    left_join(map(res, \(r) tibble(config = r$config,
                                   cl_mean = mean(r$per_trial$cl_mean),
                                   true_cl = r$true_cl)) |> list_rbind(),
              by = "config") |>
    arrange(rep, arm)
rec |> mutate(across(where(is.numeric), \(x) round(x, 3))) |> print(n = Inf)

paired <- rec |> filter(rep %in% rec$rep[duplicated(rec$rep)])

if (n_distinct(paired$arm) == 2L) {
    m_def <- mean(paired$cl_mean[paired$arm == "default"])
    m_wid <- mean(paired$cl_mean[paired$arm == "wide"])
    tr <- paired$true_cl[1]
    sd_def <- paired$sd_logit_cl[paired$arm == "default"][1]
    sd_wid <- paired$sd_logit_cl[paired$arm == "wide"][1]
    lo <- d$min_cl; hi <- d$max_cl
    prior_ctr <- inv_lg(d$mean_logit_cl, lo, hi)

    cat(sprintf("\n  paired on replicates %s\n",
                paste(sort(unique(paired$rep)), collapse = ", ")))
    cat(sprintf("  truth           %.3f h (logit %.4f)\n", tr, lg(tr, lo, hi)))
    cat(sprintf("  prior centre    %.3f h (logit %.4f)\n",
                prior_ctr, d$mean_logit_cl))
    cat(sprintf("  sd_logit_cl %-4g %.3f h (logit %.4f), bias %+.3f h\n",
                sd_def, m_def, lg(m_def, lo, hi), m_def - tr))
    cat(sprintf("  sd_logit_cl %-4g %.3f h (logit %.4f), bias %+.3f h\n",
                sd_wid, m_wid, lg(m_wid, lo, hi), m_wid - tr))

    ## Normal-normal on the logit scale: posterior = w * prior + (1 - w) * MLE
    ## with w = v / (v + sd_logit_cl^2) and the same likelihood variance v in
    ## both arms. Two arms, two unknowns. This is an approximation on two
    ## replicates -- read the sign and rough size, not the third digit.
    a <- lg(m_def, lo, hi); b <- lg(m_wid, lo, hi); pr <- d$mean_logit_cl
    root <- uniroot(\(v) (a * (v + sd_def^2) - v * pr) / sd_def^2 -
                         (b * (v + sd_wid^2) - v * pr) / sd_wid^2,
                    c(1e-8, 1e4))
    v <- root$root
    mle <- a * (v + sd_def^2) - v * pr
    mle <- mle / 1  ## logit-scale likelihood-only location

    cat(sprintf("\n  implied prior weight: %.3f at sd %g, %.3f at sd %g\n",
                v / (v + sd_def^2), sd_def, v / (v + sd_wid^2), sd_wid))
    cat(sprintf("  implied likelihood-only estimate: %.3f h (bias %+.3f h)\n",
                inv_lg(mle, lo, hi), inv_lg(mle, lo, hi) - tr))
    cat("\n  Quadrupling the prior variance barely moves the estimate, so most\n",
        "  of the miss is the likelihood's own, not the prior pulling toward\n",
        "  the prior centre. Explanation 1 cannot carry this.\n", sep = "")
}


# ------------------------------------------------------------------------ #
# Is the simulated correlation as extreme as the real one?
#
# Two statistics, kept apart, because mixing them is misleading. The
# correlation of posterior MEANS conditions on the per-trial estimates as if
# known. The correlation WITHIN a draw carries each trial's uncertainty and is
# attenuated by it. A simulated value is only comparable to the real value
# computed the same way; comparing the real point estimate against a simulated
# per-draw distribution would guarantee an extreme-looking answer and mean
# nothing.
# ------------------------------------------------------------------------ #

cat("\n=== correlation with obs_per_series: simulated vs real ===\n")

cat(sprintf("\n  posterior means -- real %+.3f\n", real$r_means))
for (a in unique(tab$arm)) {
    rr <- sort(tab$r_means[tab$arm == a])
    cat(sprintf("    %-8s n=%d: %s | most negative %+.3f\n", a, length(rr),
                paste(sprintf("%+.3f", rr), collapse = " "), min(rr)))
}
n_le <- sum(tab$r_means <= real$r_means)
cat(sprintf("    simulated replicates at least as negative as real: %d of %d\n",
            n_le, nrow(tab)))

cat(sprintf("\n  within draw -- real mean %+.3f (95%% %+.3f to %+.3f)\n",
            mean(real$r_draws), quantile(real$r_draws, 0.025),
            quantile(real$r_draws, 0.975)))
for (a in unique(tab$arm)) {
    pooled <- unlist(map(res[map_chr(res, "arm") == a], "r_draws"))
    cat(sprintf("    %-8s mean %+.3f (95%% %+.3f to %+.3f)\n", a,
                mean(pooled), quantile(pooled, 0.025), quantile(pooled, 0.975)))
    cat(sprintf("             P(simulated draw <= real draw, both sampled) = %.3f\n",
                mean(sample(pooled, 2e5, replace = TRUE) <=
                     sample(real$r_draws, 2e5, replace = TRUE))))
}

cat(sprintf(paste0("\n  With only %d noise realizations the spread across replicates is\n",
                   "  itself poorly pinned down. Read whether the simulated correlations\n",
                   "  reach the real one, not a p-value.\n"),
            n_distinct(tab$rep)))


# ------------------------------------------------------------------------ #
# Per-trial detail: is the bias structured by sampling density?
# ------------------------------------------------------------------------ #

per_trial_all <- map(res, \(r) r$per_trial |>
                         mutate(config = r$config, arm = r$arm,
                                bias = cl_mean - r$true_cl)) |>
    list_rbind()

cat("\n=== per-trial bias, averaged over replicates within arm ===\n")
per_trial_all |>
    summarise(.by = c(arm, trial, obs_per_series),
              bias = mean(bias), cl_sd = mean(cl_sd)) |>
    arrange(arm, obs_per_series) |>
    mutate(across(where(is.numeric), \(x) round(x, 3))) |>
    print(n = Inf)

cat("\nreal per-trial means, for reference:\n")
tibble(trial = trial_names, cl_mean = colMeans(cl_real),
       obs_per_series = ops) |>
    arrange(obs_per_series) |>
    mutate(across(where(is.numeric), \(x) round(x, 3))) |>
    print(n = Inf)
