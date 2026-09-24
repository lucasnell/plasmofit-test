# Can the inoculation sizes supply an informative prior for log10_total0?
#
# The Wockner data records inoc_size, viable parasites inoculated per subject.
# The model's log10_total0 is log10 of total parasite density at t = 0 in the
# units of y (generate_starts: total0 = pow(10, log10_total0), and
# beta_starts returns total0 * softmax(...), so the compartments sum to it).
# The data are iRBC per mL. So an inoculum of N parasites distributed through
# a 5 L blood volume predicts
#
#     log10_total0 ~ log10(N / 5000 mL)
#
# per grp_init group, replacing the current generic normal(1, 0.25) that is
# the same for every group regardless of how much was injected.
#
# This checks the idea against the saved real fit before anything is changed:
# does the fitted log10_total0 actually track the inoculum across groups, and
# with what offset? A consistent offset is interpretable (the fraction of the
# inoculum that establishes, plus anything t = 0 does not capture). No
# relationship at all would mean the two quantities are not measuring the same
# thing and the prior would fight the likelihood.
#
# Assumes time is measured from inoculation, so that t = 0 is the injection.
# Observations start at 72 h, so nothing in the data pins t = 0 directly.
#
#   srun -N 1 -n 1 -c 4 --mem=16G --time=00:30:00 \
#       Rscript --vanilla _scripts/wockner-inoc-prior.R

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")

suppressPackageStartupMessages({
    library(plasmofit)
    library(tidyverse)
    library(rstan)
})

BLOOD_ML <- 5000        # 5 L, adult; varies ~10% with body size, i.e. ~0.04 log10

paras_df <- read_csv("_data/wockner-cleaned.csv", col_types = "cccdcdd") |>
    mutate(inoc = interaction(trial, inoc_size, drop = TRUE) |> paste()) |>
    mutate(obs_error = interaction(trial, cohort, drop = TRUE) |> paste())

d <- archer_stan_data(paras_df,
                      series = "id", time = "time", abundance = "para",
                      grp_init = "inoc", grp_R = "trial", grp_cl = "trial",
                      grp_sd = "obs_error", calc_log_lik = 0L)

## grp_init codes follow order of first appearance; the levels attribute maps
## them back to the trial.inoc_size labels
inoc_levels <- attr(d, "levels")[[attr(d, "grp_columns")$grp_init]]

meta <- paras_df |>
    summarise(.by = inoc,
              trial = first(trial),
              inoc_size = first(inoc_size),
              n_series = n_distinct(id),
              n_subj = n_distinct(subject)) |>
    mutate(grp_init = match(inoc, inoc_levels)) |>
    arrange(grp_init)
stopifnot(!anyNA(meta$grp_init), nrow(meta) == d$n_grp_init)

f <- read_rds("_data/wock-fit-no_pool.rds")
lt <- as.matrix(f, pars = "log10_total0")
rm(f); invisible(gc())

tab <- meta |>
    mutate(pred = log10(inoc_size / BLOOD_ML),
           fitted = colMeans(lt),
           fitted_sd = apply(lt, 2, sd),
           diff = fitted - pred,
           fold = 10^(fitted - pred))

cat("=== fitted log10_total0 vs the inoculum prediction ===\n")
tab |>
    select(grp_init, trial, inoc_size, n_series, pred, fitted, fitted_sd,
           diff, fold) |>
    mutate(across(where(is.numeric), \(x) round(x, 3))) |>
    print(n = Inf, width = Inf)

cat(sprintf("\n  distinct inoculum sizes: %s\n",
            paste(sort(unique(tab$inoc_size)), collapse = ", ")))
cat(sprintf("  inoculum spans %.2f log10 units; fitted spans %.2f\n",
            diff(range(tab$pred)), diff(range(tab$fitted))))

if (n_distinct(tab$inoc_size) > 1) {
    cat(sprintf("  corr(fitted, predicted) across %d groups: %+.3f\n",
                nrow(tab), cor(tab$fitted, tab$pred)))
    fit_lm <- lm(fitted ~ pred, data = tab)
    cat(sprintf("  slope %.3f (se %.3f), intercept %.3f\n",
                coef(fit_lm)[2], summary(fit_lm)$coefficients[2, 2],
                coef(fit_lm)[1]))
    cat("  (slope 1 would mean the fit tracks the inoculum one-for-one)\n")
} else {
    cat("  only one distinct inoculum size: the prediction is a constant,\n")
    cat("  so it cannot be checked by correlation across groups. The offset\n")
    cat("  below is still informative about its level.\n")
}

cat(sprintf("\n  offset (fitted - predicted): mean %+.3f log10, sd %.3f\n",
            mean(tab$diff), sd(tab$diff)))
cat(sprintf("  implied multiplier on the inoculum: %.2fx (range %.2f-%.2f)\n",
            10^mean(tab$diff), min(tab$fold), max(tab$fold)))
cat(sprintf("  current prior: normal(%g, %g), the same for all %d groups\n",
            d$mean_log10_total0[1], d$sd_log10_total0[1], d$n_grp_init))

cat("\n  A multiplier below 1 is the fraction of the inoculum that\n",
    "  establishes, which is the expected direction. Above 1 means the fit\n",
    "  wants more parasites at t = 0 than were injected, which the inoculum\n",
    "  cannot explain and would need a different account -- growth before\n",
    "  t = 0 is not one, since t = 0 is the injection.\n", sep = "")



# ------------------------------------------------------------------------ #
# What would the anchor do to R?
#
# log10_total0 and R correlate about -0.89 in every grp_init group (see
# wockner-ridge.R), so what the data pin is roughly a combination of the two,
# not either alone. Moving log10_total0 to the inoculum anchor therefore
# moves R. Read the trade-off off the posterior itself -- regress log10(R) on
# log10_total0 across draws -- rather than from cycle arithmetic.
#
# This is an extrapolation well outside the posterior's support: the anchor
# sits ~0.8 log10 below the posterior mean, which is several posterior sds.
# It gives the direction and rough size, not a prediction.
#
# Burst size bounds it. Greischar & Childs 2023 (mmcm.pdf): median burst size
# 15-18, maximum 32 for P. falciparum, in vitro PMRs typically below 15, and
# ~8 for 3D7 specifically. Wockner et al.'s own in vivo estimates of 28.7-35.4
# are what that paper argues are sequestration artifacts.
# ------------------------------------------------------------------------ #

f2 <- read_rds("_data/wock-fit-no_pool.rds")
R_draws <- as.matrix(f2, pars = "R")
lt_draws <- as.matrix(f2, pars = "log10_total0")
rm(f2); invisible(gc())

grpR_of <- map_int(seq_len(d$n_grp_init),
                   \(g) d$grp_R[which(d$grp_init == g)[1]])

implied <- map(seq_len(d$n_grp_init), \(g) {
    x <- lt_draws[, g]; y <- log10(R_draws[, grpR_of[g]])
    b <- coef(lm(y ~ x))
    tibble(grp_init = g, trial = tab$trial[tab$grp_init == g],
           slope = b[2],
           R_now = mean(R_draws[, grpR_of[g]]),
           R_at_anchor = 10^(b[1] + b[2] * tab$pred[tab$grp_init == g]))
}) |> list_rbind()

cat("\n=== what the inoculum anchor implies for R ===\n")
implied |> mutate(across(where(is.numeric), \(x) round(x, 2))) |> print(n = Inf)

cat(sprintf("\n  R now:        mean %.2f (range %.2f-%.2f)\n",
            mean(implied$R_now), min(implied$R_now), max(implied$R_now)))
cat(sprintf("  R at anchor:  mean %.2f (range %.2f-%.2f)\n",
            mean(implied$R_at_anchor), min(implied$R_at_anchor),
            max(implied$R_at_anchor)))
cat(sprintf("  groups with R > 16 at the anchor: %d of %d\n",
            sum(implied$R_at_anchor > 16), nrow(implied)))
cat(sprintf("  groups with R > 32 at the anchor: %d of %d\n",
            sum(implied$R_at_anchor > 32), nrow(implied)))
cat(sprintf("\n  current max_R = %g, above the burst-size maximum of 32\n",
            d$max_R))
cat("  reference points: in vitro 3D7 ~8, median burst size 15-18,\n")
cat("  maximum burst size 32, Wockner et al.'s own in vivo 28.7-35.4\n")

write_rds(list(tab = tab, implied = implied), "_data/wock-inoc-prior.rds")
cat("\nwrote _data/wock-inoc-prior.rds\n")
