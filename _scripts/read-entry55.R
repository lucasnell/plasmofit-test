## Read wockner-fit.R entry 55, np_start13_nc384_t (SLURM 35287, 2026-11-18),
## against entry 54, np_wide_total0_nc384_t (SLURM 30842), which it differs
## from ONLY in the starting age distribution:
##   54: 14 per-unit b_shape (lognormal(2, 0.5)) and b_offset, a wrapped
##       symmetric Beta over the whole cycle, split circulating/sequestered
##       by age at t = 0;
##   55: one Beta(start_alpha, start_beta) on [0, 13] h since invasion,
##       shared by all 14 units, all circulating at t = 0, flat priors on the
##       shapes over [1, 3] (claude/start-dist-design.md).
## Everything else is the same: no_pool, n_c 384, sd_log10_total0 1,
## adapt_delta 0.95, max_treedepth 12, seed 538065874, 4 x 1000 draws.
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   Rscript --vanilla _scripts/read-entry55.R > _data/read-entry55-<date>.txt

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
suppressPackageStartupMessages({library(rstan); library(posterior); library(loo)})
options(width = 140)
D <- "_data/"
f55 <- readRDS(paste0(D, "wock-fit-np_start13_nc384_t.rds"))
d55 <- readRDS(paste0(D, "wock-data-np_start13_nc384_t.rds"))
f54 <- readRDS(paste0(D, "wock-fit-np_wide_total0_nc384_t.rds"))
d54 <- readRDS(paste0(D, "wock-data-np_wide_total0_nc384_t.rds"))
lv <- attr(d55, "levels")
stopifnot(identical(lv, attr(d54, "levels")))
W <- d55$start_window; CAP <- d55$start_shape_cap

## ---- 1. Sampler health ----------------------------------------------------
health <- function(f, label) {
    s <- summary(f)$summary
    rh <- s[, "Rhat"]; rh <- rh[is.finite(rh)]
    sp <- get_sampler_params(f, inc_warmup = FALSE)
    lp <- apply(extract(f, "lp__", permuted = FALSE)[, , 1, drop = FALSE], 2, mean)
    el <- rowSums(get_elapsed_time(f)) / 3600
    cat(sprintf("%s: max R-hat %.3f (%s) -> %s; divergent %d of %d post-warmup iterations;\n",
                label, max(rh), names(which.max(rh)),
                if (max(rh) < 1.05) "PASSES the gate" else "FAILS the gate",
                sum(unlist(lapply(sp, \(x) x[, "divergent__"]))),
                sum(vapply(sp, nrow, integer(1)))))
    cat(sprintf("   at treedepth cap 12: %.1f%%; leapfrog steps per iteration %.1f; chain hours %s\n",
                100 * mean(unlist(lapply(sp, \(x) x[, "treedepth__"] >= 12))),
                mean(unlist(lapply(sp, \(x) x[, "n_leapfrog__"]))),
                paste(sprintf("%.2f", el), collapse = ", ")))
    cat(sprintf("   per-chain mean lp__: %s (spread %.2f)\n",
                paste(sprintf("%.1f", lp), collapse = ", "), diff(range(lp))))
    cat(sprintf("   bulk ESS min %.0f, tail ESS min %.0f\n\n",
                min(s[, "n_eff"], na.rm = TRUE),
                min(apply(as_draws_array(f), 3, ess_tail), na.rm = TRUE)))
}
cat("1. SAMPLER HEALTH (4 chains x 1000 post-warmup draws each)\n\n")
health(f55, "entry 55 (windowed, shared start)")
health(f54, "entry 54 (per-unit b_shape, b_offset)")

## ---- 2. cycle_length ---------------------------------------------------------
cat("2. CYCLE LENGTH. Per trial: posterior mean and sd of cycle_length (hours).\n",
    "   'mean over trials' = posterior mean (sd) of the per-draw mean over the 13\n",
    "   trials. 55 - 54 = difference of posterior means; unpaired (different fits\n",
    "   of the same data).\n\n", sep = "")
c55 <- as.matrix(f55, pars = "cycle_length"); c54 <- as.matrix(f54, pars = "cycle_length")
tab <- data.frame(trial = lv$trial,
                  e55_mean = colMeans(c55), e55_sd = apply(c55, 2, sd),
                  e54_mean = colMeans(c54), e54_sd = apply(c54, 2, sd))
tab$diff <- tab$e55_mean - tab$e54_mean
tab[-1] <- round(tab[-1], 2)
print(tab, row.names = FALSE)
m55 <- rowMeans(c55); m54 <- rowMeans(c54)
cat(sprintf("\n   mean over trials: entry 55 %.2f h (sd %.2f), entry 54 %.2f h (sd %.2f); 55 - 54 = %+.2f h\n",
            mean(m55), sd(m55), mean(m54), sd(m54), mean(m55) - mean(m54)))
cat(sprintf("   trials with 55 shorter than 54: %d of 13; range of 55 - 54 %+.2f to %+.2f h\n\n",
            sum(tab$diff < 0), min(tab$diff), max(tab$diff)))

## ---- 3. The starting distribution -------------------------------------------
cat("3. THE STARTING AGE DISTRIBUTION (entry 55, one for all units).\n",
    "   start_alpha, start_beta: posterior mean and 90% interval; flat prior on\n",
    "   [1, ", CAP, "]. Implied start, hours since invasion on [0, ", W, "], per draw:\n",
    "   mean age W a/(a+b), sd W sqrt(ab/((a+b)^2(a+b+1))), mode W (a-1)/(a+b-2).\n\n", sep = "")
a <- as.vector(as.matrix(f55, pars = "start_alpha")); b <- as.vector(as.matrix(f55, pars = "start_beta"))
q <- function(v) sprintf("%.2f [%.2f, %.2f]", mean(v), quantile(v, 0.05), quantile(v, 0.95))
cat("   start_alpha ", q(a), "\n   start_beta  ", q(b), "\n", sep = "")
cat(sprintf("   share of draws within 0.1 of a bound: alpha %.2f (at 1) %.2f (at %g); beta %.2f (at 1) %.2f (at %g)\n",
            mean(a < 1.1), mean(a > CAP - 0.1), CAP, mean(b < 1.1), mean(b > CAP - 0.1), CAP))
mean_age <- W * a / (a + b); sd_age <- W * sqrt(a * b / ((a + b)^2 * (a + b + 1)))
mode_age <- ifelse(a > 1 & b > 1, W * (a - 1) / (a + b - 2), NA)
cat("   mean age (h) ", q(mean_age), "\n   sd of ages (h) ", q(sd_age), "\n",
    "   mode (h)     ", q(mode_age[is.finite(mode_age)]), "\n", sep = "")
cat(sprintf("   (a uniform start on [0, %g] has mean %.2f h and sd %.2f h)\n", W, W / 2, W / sqrt(12)))
cat("   For comparison, entry 54's 14 per-unit starts (posterior means; from\n",
    "   _data/inoculum-age-check-2026-11-09.txt): modal age 1.1-11.7 h, sd 3.3-6.1 h,\n",
    "   with 7-42% of the starting mass in sequestered compartments.\n\n", sep = "")

## ---- 4. The other parameters --------------------------------------------------
cat("4. R (per trial) AND log10_total0 (per unit): posterior means, 55 against 54.\n\n")
R55 <- colMeans(as.matrix(f55, pars = "R")); R54 <- colMeans(as.matrix(f54, pars = "R"))
cat(sprintf("   R: entry 55 %.2f to %.2f (median %.2f); entry 54 %.2f to %.2f (median %.2f); 55 - 54 median %+.2f\n",
            min(R55), max(R55), median(R55), min(R54), max(R54), median(R54), median(R55 - R54)))
t55 <- colMeans(as.matrix(f55, pars = "log10_total0")); t54 <- colMeans(as.matrix(f54, pars = "log10_total0"))
cat(sprintf("   log10_total0: entry 55 %.2f to %.2f (median %.2f); entry 54 %.2f to %.2f (median %.2f); 55 - 54 median %+.2f\n",
            min(t55), max(t55), median(t55), min(t54), max(t54), median(t54), median(t55 - t54)))
s55 <- colMeans(as.matrix(f55, pars = "sd_iRBC")); s54 <- colMeans(as.matrix(f54, pars = "sd_iRBC"))
cat(sprintf("   sd_iRBC: median over 27 cohorts entry 55 %.3f, entry 54 %.3f; 55 - 54 median %+.3f\n\n",
            median(s55), median(s54), median(s55 - s54)))

## ---- 5. Predictive comparison --------------------------------------------------
cat("5. PSIS-LOO over the same 1130 observations (the 'observation' element of\n",
    "   each LOO file). elpd_diff is relative to the better model; positive\n",
    "   elpd is better.\n\n", sep = "")
L55 <- readRDS(paste0(D, "wock-fit-LOO-np_start13_nc384_t.rds"))$observation
L54 <- readRDS(paste0(D, "wock-fit-LOO-np_wide_total0_nc384_t.rds"))$observation
print(loo_compare(list(entry55_windowed_start = L55, entry54_per_unit_start = L54)))
cat(sprintf("\n   elpd_loo: entry 55 %.1f (se %.1f), entry 54 %.1f (se %.1f)\n",
            L55$estimates["elpd_loo", 1], L55$estimates["elpd_loo", 2],
            L54$estimates["elpd_loo", 1], L54$estimates["elpd_loo", 2]))
cat(sprintf("   p_loo: entry 55 %.1f, entry 54 %.1f; Pareto k > 0.7: entry 55 %d, entry 54 %d\n",
            L55$estimates["p_loo", 1], L54$estimates["p_loo", 1],
            sum(L55$diagnostics$pareto_k > 0.7), sum(L54$diagnostics$pareto_k > 0.7)))

## ---- 6. Per chain (the gate failed: which chains are where) -----------------
## Added after the first read showed max R-hat 18.3. Per chain: mean lp__,
## mean of the SUMMED pointwise log-likelihood over the 1130 observations
## (log_lik, priors excluded, so comparable across the two models where lp__
## is not), the start shapes, the mean cycle_length over 13 trials, and the
## trials whose entry-55 posterior sd exceeded 1 h.
cat("\n6. PER CHAIN. lp__ and summed log_lik: mean over the chain's 1000 draws.\n",
    "   cycle_length columns: chain mean, hours.\n\n", sep = "")
per_chain <- function(f, start) {
    arr <- function(p) extract(f, p, permuted = FALSE)
    lp <- apply(arr("lp__")[, , 1, drop = FALSE], 2, mean)
    ll <- apply(arr("log_lik"), 2, \(m) mean(rowSums(m)))
    cl <- arr("cycle_length")
    out <- data.frame(chain = seq_along(lp), lp = round(lp, 1), sum_log_lik = round(ll, 1),
                      cl_mean13 = round(apply(cl, 2, \(m) mean(rowMeans(m))), 2))
    if (start) {
        out$start_alpha <- round(apply(arr("start_alpha")[, , 1, drop = FALSE], 2, mean), 2)
        out$start_beta <- round(apply(arr("start_beta")[, , 1, drop = FALSE], 2, mean), 2)
    }
    wide <- which(apply(as.matrix(f55, pars = "cycle_length"), 2, sd) > 1)
    for (k in wide) out[[lv$trial[k]]] <- round(apply(cl[, , k], 2, mean), 2)
    out
}
cat("   entry 55\n"); print(per_chain(f55, TRUE), row.names = FALSE)
cat("\n   entry 54 (converged)\n"); print(per_chain(f54, FALSE), row.names = FALSE)
