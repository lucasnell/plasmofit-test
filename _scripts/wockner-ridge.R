# Is the (log10_total0, R, b_shape) ridge a property of the real fit, or only
# of the simulation?
#
# The schedule-bias simulation does not recover its own nuisance parameters:
# b_shape lands on its prior median instead of the truth, log10_total0 moves
# halfway to its prior mean, and the two trade off against R (see
# claude/CLAUDE.md, "The simulation does not recover its own nuisance
# parameters"). Two readings, with different consequences:
#
#   a) the simulation was built badly -- generating from the pooled_cl fit's
#      posterior MEAN vector broke the parameter correlations, so the
#      simulated data are less informative than the real data and the ridge
#      is an artifact of that
#   b) the model is weakly identified along this ridge on the real data too,
#      in which case b_shape, log10_total0 and R in the saved real fits are
#      substantially prior-driven and should not be reported as estimates
#
# The discriminator is cheap and needs no refitting: if the ridge is real, the
# posterior draws of these parameters are strongly correlated *within* the
# real fit. A posterior that pins each parameter independently cannot produce
# the trade-off that (a) and (b) both describe.
#
# Also compares each parameter's posterior against its own prior. A posterior
# sitting on its prior with the prior's width has been told little by the
# data, whichever reading is right.
#
#   srun -N 1 -n 1 -c 4 --mem=16G --time=00:30:00 \
#       Rscript --vanilla _scripts/wockner-ridge.R

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")

suppressPackageStartupMessages({
    library(plasmofit)
    library(tidyverse)
    library(rstan)
})

paras_df <- read_csv("_data/wockner-cleaned.csv", col_types = "cccdcdd") |>
    mutate(inoc = interaction(trial, inoc_size, drop = TRUE) |> paste()) |>
    mutate(obs_error = interaction(trial, cohort, drop = TRUE) |> paste())

d <- archer_stan_data(paras_df,
                      series = "id", time = "time", abundance = "para",
                      grp_init = "inoc", grp_R = "trial", grp_cl = "trial",
                      grp_sd = "obs_error", calc_log_lik = 0L)

f <- read_rds("_data/wock-fit-no_pool.rds")

## grp_init group 1, the one b_shape/b_offset/log10_total0 are indexed by;
## R is per trial, so pair with the trial that group's series belong to.
gi <- 1L
series_gi <- which(d$grp_init == gi)
gR <- d$grp_R[series_gi[1]]
cat(sprintf("using grp_init %d (%d series) and its grp_R %d\n\n",
            gi, length(series_gi), gR))

draws <- tibble(
    b_shape      = as.matrix(f, pars = "b_shape")[, gi],
    b_offset     = as.matrix(f, pars = "b_offset")[, gi],
    log10_total0 = as.matrix(f, pars = "log10_total0")[, gi],
    R            = as.matrix(f, pars = "R")[, gR],
    cycle_length = as.matrix(f, pars = "cycle_length")[, gR]
)


# ------------------------------------------------------------------------ #
# Is there a ridge?
# ------------------------------------------------------------------------ #

cat("=== posterior correlations among nuisance parameters, real no_pool fit ===\n")
cm <- cor(draws)
print(round(cm, 3))

cat("\n  |r| > 0.7 pairs:\n")
idx <- which(abs(cm) > 0.7 & upper.tri(cm), arr.ind = TRUE)
if (nrow(idx) == 0) {
    cat("    none -- no strong ridge in the real posterior\n")
} else {
    for (k in seq_len(nrow(idx))) {
        cat(sprintf("    %-13s %-13s %+.3f\n", rownames(cm)[idx[k, 1]],
                    colnames(cm)[idx[k, 2]], cm[idx[k, 1], idx[k, 2]]))
    }
}


# ------------------------------------------------------------------------ #
# How much did the data say? Posterior against prior, on each parameter's
# own scale.
# ------------------------------------------------------------------------ #

cat("\n=== posterior vs prior ===\n")

## b_shape ~ lognormal(mean_log_b_shape, sd_log_b_shape)
## log10_total0 ~ normal(mean_log10_total0, sd_log10_total0)
## R = max_R * inv_logit(eta_R), eta_R ~ normal(mu_logit_R, sigma_logit_R),
##   mu_logit_R ~ normal(mean_logit_R, sd_logit_R)
prior <- tibble(
    par = c("b_shape", "log10_total0", "R"),
    prior_loc = c(exp(d$mean_log_b_shape[gi]),
                  d$mean_log10_total0[gi],
                  d$max_R / (1 + exp(-d$mean_logit_R))),
    prior_sd_desc = c(sprintf("lognormal(%g, %g), median %.2f",
                              d$mean_log_b_shape[gi], d$sd_log_b_shape[gi],
                              exp(d$mean_log_b_shape[gi])),
                      sprintf("normal(%g, %g)", d$mean_log10_total0[gi],
                              d$sd_log10_total0[gi]),
                      sprintf("max_R * inv_logit(normal(%g, %g)), median %.2f",
                              d$mean_logit_R, d$sd_logit_R,
                              d$max_R / (1 + exp(-d$mean_logit_R))))
)

post <- map(prior$par, \(nm) {
    x <- draws[[nm]]
    tibble(par = nm, post_mean = mean(x), post_sd = sd(x),
           q025 = quantile(x, 0.025), q975 = quantile(x, 0.975))
}) |> list_rbind()

left_join(post, prior, by = "par") |>
    mutate(post_minus_prior = post_mean - prior_loc) |>
    select(par, prior_loc, post_mean, post_minus_prior, post_sd, q025, q975,
           prior_sd_desc) |>
    mutate(across(where(is.numeric), \(x) round(x, 3))) |>
    print(width = Inf)

## Prior-predictive draws on the same scale, to say whether the posterior is
## narrower than the prior at all. Same comparison, without assuming
## normality on any transformed scale.
set.seed(20260924)
n <- 20000
pp <- tibble(
    b_shape = rlnorm(n, d$mean_log_b_shape[gi], d$sd_log_b_shape[gi]),
    log10_total0 = rnorm(n, d$mean_log10_total0[gi], d$sd_log10_total0[gi]),
    R = d$max_R / (1 + exp(-rnorm(n, d$mean_logit_R, d$sd_logit_R)))
)

## Compared on the scale each prior is written on. For b_shape that is the
## log scale: a lognormal's sd on the natural scale grows with its location,
## so a posterior centred well above the prior would look "wider" than the
## prior even when the data have sharply constrained it.
sd_ratio <- c(
    b_shape      = sd(log(draws$b_shape)) / d$sd_log_b_shape[gi],
    log10_total0 = sd(draws$log10_total0) / d$sd_log10_total0[gi],
    R            = sd(draws$R) / sd(pp$R)
)

cat("\n  posterior sd as a fraction of prior sd (1 = data said nothing):\n")
cat(sprintf("    %-13s %.3f   (log scale: post %.3f, prior %.3f)\n", "b_shape",
            sd_ratio[["b_shape"]], sd(log(draws$b_shape)), d$sd_log_b_shape[gi]))
cat(sprintf("    %-13s %.3f   (post %.3f, prior %.3f)\n", "log10_total0",
            sd_ratio[["log10_total0"]], sd(draws$log10_total0),
            d$sd_log10_total0[gi]))
cat(sprintf("    %-13s %.3f   (post %.3f, prior %.3f)\n", "R",
            sd_ratio[["R"]], sd(draws$R), sd(pp$R)))

cat("\n=== does cycle_length ride on them? ===\n")
for (nm in c("b_shape", "b_offset", "log10_total0", "R")) {
    cat(sprintf("    corr(cycle_length, %-13s) %+.3f\n", nm,
                cor(draws$cycle_length, draws[[nm]])))
}


# ------------------------------------------------------------------------ #
# Across every grp_init group, not just one
#
# One group could be unrepresentative, and the pairs above are the whole
# finding, so check them everywhere before concluding anything.
# ------------------------------------------------------------------------ #

bs_all <- as.matrix(f, pars = "b_shape")
bo_all <- as.matrix(f, pars = "b_offset")
lt_all <- as.matrix(f, pars = "log10_total0")
R_all  <- as.matrix(f, pars = "R")
cl_all <- as.matrix(f, pars = "cycle_length")

all_grp <- map(seq_len(d$n_grp_init), \(g) {
    sers <- which(d$grp_init == g)
    r_idx <- d$grp_R[sers[1]]
    tibble(grp_init = g, n_series = length(sers), grp_R = r_idx,
           lt_R = cor(lt_all[, g], R_all[, r_idx]),
           bo_cl = cor(bo_all[, g], cl_all[, r_idx]),
           bs_lt = cor(bs_all[, g], lt_all[, g]),
           bs_R = cor(bs_all[, g], R_all[, r_idx]),
           bs_cl = cor(bs_all[, g], cl_all[, r_idx]))
}) |> list_rbind()

cat("\n=== the same correlations in every grp_init group ===\n")
all_grp |> mutate(across(where(is.numeric), \(x) round(x, 3))) |> print(n = Inf)

cat("\n  summary over the", d$n_grp_init, "groups:\n")
for (nm in c("lt_R", "bo_cl", "bs_lt", "bs_R", "bs_cl")) {
    v <- all_grp[[nm]]
    cat(sprintf("    %-6s median %+.3f  range %+.3f to %+.3f\n",
                nm, median(v), min(v), max(v)))
}
cat("\n  lt_R  = log10_total0 vs R      bo_cl = b_offset vs cycle_length\n")
cat("  bs_*  = b_shape vs each of log10_total0, R, cycle_length\n")

write_rds(list(cor = cm, post = post, prior = prior, sd_ratio = sd_ratio,
               all_grp = all_grp),
          "_data/wock-ridge.rds")
cat("\nwrote _data/wock-ridge.rds\n")
