# Reads the first inoculum-anchored fit against the predictions made from the
# unanchored one, and against the unanchored fit itself.
#
# The predictions, from wockner-inoc-prior.R on _data/wock-fit-no_pool.rds:
#   delta_total0 ~ +0.811  (the fit wants 6.47x more parasites at t = 0 than
#                           were inoculated; blood volume cannot absorb it,
#                           +-10% moves the anchor only 0.041)
#   R            ~ 9.96    (implied at the anchor by regressing log10(R) on
#                           log10_total0 across draws, slope -0.25 along the
#                           log10_total0/R ridge)
#
# The second question is thread 2's: log10_total0's prior is a suspect for
# dragging cycle_length along that ridge. Replacing it with a data-derived
# anchor is an independent route at the same target, so cycle_length is
# compared too.
#
#   srun -N 1 -n 1 -c 4 --mem=32G Rscript --vanilla \
#       _scripts/wockner-anchor-check.R

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")

suppressPackageStartupMessages({
    library(plasmofit)
    library(tidyverse)
    library(rstan)
})

PRED_DELTA <- 0.811
PRED_R <- 9.96

d <- read_rds("_data/wock-data-np_anchor.rds")
f_a <- read_rds("_data/wock-fit-np_anchor.rds")

smry <- function(f, par) {
    m <- as.matrix(f, pars = par)
    tibble(par = colnames(m), mean = colMeans(m),
           sd = apply(m, 2, sd),
           q025 = apply(m, 2, quantile, 0.025),
           q975 = apply(m, 2, quantile, 0.975))
}

cat("=== sampler ===\n")
sp <- get_sampler_params(f_a, inc_warmup = FALSE)
rh <- rstan::summary(f_a)$summary[, "Rhat"]
cat(sprintf("  divergences %d (%.1f%%), max R-hat %.4f, min n_eff %.0f\n",
            sum(sapply(sp, \(x) sum(x[, "divergent__"]))),
            100 * sum(sapply(sp, \(x) sum(x[, "divergent__"]))) / (4 * 1000),
            max(rh, na.rm = TRUE),
            min(rstan::summary(f_a)$summary[, "n_eff"], na.rm = TRUE)))


# ------------------------------------------------------------------------ #
# The offset: does it land on +0.811?
# ------------------------------------------------------------------------ #

dl <- smry(f_a, "delta_total0")
sg <- smry(f_a, "sigma_total0")

cat("\n=== the shared offset ===\n")
cat(sprintf("  delta_total0  %+.3f (95%% %+.3f to %+.3f), sd %.3f\n",
            dl$mean, dl$q025, dl$q975, dl$sd))
cat(sprintf("  predicted     %+.3f  -> %s\n", PRED_DELTA,
            if (PRED_DELTA >= dl$q025 && PRED_DELTA <= dl$q975)
                "INSIDE the 95% interval" else "OUTSIDE the 95% interval"))
cat(sprintf("  as a factor   %.2fx more parasites at t = 0 than inoculated\n",
            10^dl$mean))
cat(sprintf("  sigma_total0  %.3f (95%% %.3f to %.3f)\n",
            sg$mean, sg$q025, sg$q975))
cat(sprintf("  prior on sigma_total0: half-normal(0, %g)\n", d$sd_bs_total0))


# ------------------------------------------------------------------------ #
# Per-group: how far does each group sit from its own anchor?
# ------------------------------------------------------------------------ #

t0_a <- smry(f_a, "log10_total0")
cat("\n=== log10_total0 per grp_init, against its anchor ===\n")
tibble(group = seq_len(d$n_grp_init),
       anchor = d$anchor_log10_total0,
       anchored_fit = t0_a$mean,
       offset = t0_a$mean - d$anchor_log10_total0) |>
    mutate(across(where(is.numeric), \(x) round(x, 3))) |> print(n = Inf)
cat(sprintf("\n  sd of the per-group offsets: %.3f (sigma_total0 posterior %.3f)\n",
            sd(t0_a$mean - d$anchor_log10_total0), sg$mean))


# ------------------------------------------------------------------------ #
# Against the unanchored fit: did the anchor move anything else?
# ------------------------------------------------------------------------ #

f_u <- read_rds("_data/wock-fit-no_pool.rds")

cmp <- function(par) {
    a <- smry(f_a, par); u <- smry(f_u, par)
    tibble(par = par,
           unanchored = mean(u$mean), anchored = mean(a$mean),
           diff = mean(a$mean) - mean(u$mean),
           sd_un = mean(u$sd), sd_an = mean(a$sd),
           sd_ratio = mean(a$sd) / mean(u$sd))
}

cat("\n=== anchored vs unanchored, averaged over groups ===\n")
map(c("log10_total0", "R", "cycle_length", "b_shape", "b_offset", "sd_iRBC"),
    cmp) |> list_rbind() |>
    mutate(across(where(is.numeric), \(x) round(x, 4))) |> print(width = Inf)

## Three-way, including the control.
f_w <- read_rds("_data/wock-fit-np_wide_total0.rds")
cat("\n=== three-way: does the inoculum add anything over just widening? ===\n")
map(c("log10_total0", "R", "cycle_length", "b_shape", "sd_iRBC"), \(par) {
    tibble(par = par,
           unanchored = mean(smry(f_u, par)$mean),
           wide_total0 = mean(smry(f_w, par)$mean),
           anchored = mean(smry(f_a, par)$mean))
}) |> list_rbind() |>
    mutate(across(where(is.numeric), \(x) round(x, 3))) |> print(width = Inf)
sp_w <- get_sampler_params(f_w, inc_warmup = FALSE)
cat(sprintf("\n  np_wide_total0 sampler: divergences %d, max R-hat %.4f\n",
            sum(sapply(sp_w, \(x) sum(x[, "divergent__"]))),
            max(rstan::summary(f_w)$summary[, "Rhat"], na.rm = TRUE)))
rm(f_w); invisible(gc())

r_a <- smry(f_a, "R"); r_u <- smry(f_u, "R")
cat(sprintf("\n  R: unanchored %.2f, anchored %.2f, predicted at the anchor %.2f\n",
            mean(r_u$mean), mean(r_a$mean), PRED_R))
cat(sprintf("  R range across trials: unanchored %.2f-%.2f, anchored %.2f-%.2f\n",
            min(r_u$mean), max(r_u$mean), min(r_a$mean), max(r_a$mean)))

cl_a <- smry(f_a, "cycle_length"); cl_u <- smry(f_u, "cycle_length")
cat(sprintf("\n  cycle_length: unanchored %.3f h, anchored %.3f h, shift %+.3f h\n",
            mean(cl_u$mean), mean(cl_a$mean),
            mean(cl_a$mean) - mean(cl_u$mean)))
cat("  Thread 2 reads this shift: if the log10_total0 prior were dragging\n")
cat("  cycle_length along the ridge, replacing it with a data-derived anchor\n")
cat("  should move cycle_length. Compare against the ~0.25 h the hierarchy\n")
cat("  contributes and the ~0.15 h the cycle-length prior contributes.\n")

# ------------------------------------------------------------------------ #
# Does the ridge explain where R went?
#
# wockner-inoc-prior.R predicted R ~ 9.96 by assuming log10_total0 would land
# ON the anchor and sliding along the log10_total0/R ridge from there. If
# log10_total0 instead lands somewhere else, the same slope applied to the
# distance actually travelled should still predict R. That separates "the
# ridge is wrong" from "the predicted landing point was wrong".
# ------------------------------------------------------------------------ #

t0_u <- smry(f_u, "log10_total0")
dr_u <- as.matrix(f_u, pars = c("log10_total0", "R"))
i_t0 <- grep("^log10_total0", colnames(dr_u))
i_R <- grep("^R\\[", colnames(dr_u))
slope <- unname(coef(lm(rowMeans(log10(dr_u[, i_R])) ~
                        rowMeans(dr_u[, i_t0])))[2])

moved <- mean(t0_a$mean) - mean(t0_u$mean)
to_anchor <- mean(d$anchor_log10_total0) - mean(t0_u$mean)
cat("\n=== does the ridge account for R? ===\n")
cat(sprintf("  ridge slope, d log10(R) / d log10_total0: %+.3f\n", slope))
cat(sprintf("  if log10_total0 moved only to the anchor (%+.3f): R = %.2f\n",
            to_anchor, mean(r_u$mean) * 10^(slope * to_anchor)))
cat(sprintf("  it actually moved %+.3f: predicted R = %.2f, observed %.2f\n",
            moved, mean(r_u$mean) * 10^(slope * moved), mean(r_a$mean)))
cat("  The ridge is not what failed; the predicted landing point was.\n")


# ------------------------------------------------------------------------ #
# Predictive comparison. Both fits have calc_log_lik = 1.
# ------------------------------------------------------------------------ #

## The control np_anchor needs: np_wide_total0 keeps the normal(1, .) centre
## and only widens the sd to 1, so it relaxes the tight prior WITHOUT using
## the inoculum. If it matches np_anchor, the inoculum adds nothing and the
## finding is that normal(1, 0.25) was misspecified.
cat("\n=== predictive comparison (loo, per observation) ===\n")
pick <- function(x) if (inherits(x, "loo")) x else x[[1]]
loos <- list(
    unanchored = pick(read_rds("_data/wock-fit-LOO-no_pool.rds")),
    wide_total0 = pick(read_rds("_data/wock-fit-LOO-np_wide_total0.rds")),
    anchored = pick(read_rds("_data/wock-fit-LOO-np_anchor.rds")))
for (nm in names(loos))
    cat(sprintf("  %-12s elpd %.1f (se %.1f)\n", nm,
                loos[[nm]]$estimates["elpd_loo", 1],
                loos[[nm]]$estimates["elpd_loo", 2]))
## loo_compare() can return a data frame carrying non-numeric columns, so
## round only what is numeric rather than the whole object.
cmpl <- loo::loo_compare(loos)
cmp_print <- as.data.frame(cmpl)
num <- vapply(cmp_print, is.numeric, logical(1))
cmp_print[num] <- lapply(cmp_print[num], \(x) round(x, 2))
print(cmp_print[, c("model", "elpd_diff", "se_diff")[
    c("model", "elpd_diff", "se_diff") %in% names(cmp_print)], drop = FALSE])
cat("\n  A difference of a few se is a real predictive gain. Read the sign:\n")
cat("  the row at the top is preferred.\n")

out <- list(delta_total0 = dl, sigma_total0 = sg,
            ridge_slope = slope, loo_compare = cmpl,
            anchor = d$anchor_log10_total0,
            log10_total0 = t0_a, R = r_a, cycle_length = cl_a,
            pred_delta = PRED_DELTA, pred_R = PRED_R)
write_rds(out, "_data/wock-anchor-check.rds")
cat("\nwrote _data/wock-anchor-check.rds\n")
