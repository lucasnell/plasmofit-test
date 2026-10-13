## Validate the pinned one-chain-per-L3-group driver (_scripts/wockner-fit-pinned.sh)
## against the default layout, on entry 43 (np_bs15_nc384):
##   original  SLURM 30830_43, default layout, WOCKFIT_SEED 20261010
##   pinned    SLURM 30856,    pinned layout,  WOCKFIT_SEED 20261010
## Same model, data, seed, and sampler settings. Draws cannot match bit for
## bit (threads change floating-point summation order, and HMC amplifies a
## last-bit difference within a few steps), so the posteriors are compared
## statistically, not draw by draw.
##
## Pass criteria, fixed before the pinned run finished:
##   1. the pinned fit passes the gate (max R-hat < 1.05);
##   2. per-trial posterior means of cycle_length agree: |z| < 3 for all 13,
##      where z = difference / sqrt(mcse_1^2 + mcse_2^2). MCSE understates
##      run-to-run variation (claude/gotchas.md), so this is lenient on purpose
##      and a failure is a real signal;
##   3. observation-level elpd differs by less than 2 x its paired se.
## Speed is reported, not gated.
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   Rscript --vanilla _scripts/pinned-validate.R > _data/pinned-validate-<date>.txt

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
suppressPackageStartupMessages({library(rstan); library(posterior); library(loo)})

NM <- c(original = "np_bs15_nc384", pinned = "np_bs15_nc384-pinned")
fits <- lapply(NM, function(nm) readRDS(sprintf("_data/wock-fit-%s.rds", nm)))
dats <- lapply(NM, function(nm) readRDS(sprintf("_data/wock-data-%s.rds", nm)))
stopifnot(identical(dats$original, dats$pinned))
trials <- attr(dats$original, "levels")[["trial"]]

cat("=== cost ===\n")
cat("Cells: per fit, over 4 chains. hours = longest chain (warmup + sampling);\n",
    "leapfrog_M = total gradient evaluations in millions; ms_per_grad = total\n",
    "chain time / total gradients, milliseconds.\n\n", sep = "")
cost <- do.call(rbind, lapply(names(fits), function(k) {
    f <- fits[[k]]; sp <- get_sampler_params(f, inc_warmup = TRUE); et <- get_elapsed_time(f)
    lf <- sum(sapply(sp, function(x) sum(x[, "n_leapfrog__"])))
    data.frame(fit = k, hours = round(max(rowSums(et)) / 3600, 2),
               leapfrog_M = round(lf / 1e6, 3),
               ms_per_grad = round(1000 * sum(et) / lf, 1))
}))
print(cost, row.names = FALSE)
cat(sprintf("\nspeed-up in wall time: %.1fx; in time per gradient: %.1fx\n",
            cost$hours[1] / cost$hours[2], cost$ms_per_grad[1] / cost$ms_per_grad[2]))

cat("\n=== health ===\n")
for (k in names(fits)) {
    s <- summary(fits[[k]])$summary
    lp <- apply(extract(fits[[k]], "lp__", permuted = FALSE)[, , 1, drop = FALSE], 2, mean)
    cat(sprintf("%-8s max R-hat %.4f | lp__ by chain %s | divergent %d\n", k,
                max(s[, "Rhat"], na.rm = TRUE), paste(sprintf("%.1f", lp), collapse = ", "),
                sum(get_divergent_iterations(fits[[k]]))))
}
gate <- max(summary(fits$pinned)$summary[, "Rhat"], na.rm = TRUE) < 1.05

cat("\n=== cycle_length by trial ===\n")
cat("Cells: posterior mean (MCSE of the mean) in hours, per trial (grp_cl);\n",
    "z = (pinned - original) / sqrt(mcse_o^2 + mcse_p^2).\n\n", sep = "")
dr <- lapply(fits, function(f) as_draws_array(extract(f, "cycle_length", permuted = FALSE)))
sm <- lapply(dr, function(x) summarise_draws(x, mean, mcse_mean))
z <- (sm$pinned$mean - sm$original$mean) / sqrt(sm$original$mcse_mean^2 + sm$pinned$mcse_mean^2)
tab <- data.frame(trial = trials,
                  original = sprintf("%.3f (%.3f)", sm$original$mean, sm$original$mcse_mean),
                  pinned = sprintf("%.3f (%.3f)", sm$pinned$mean, sm$pinned$mcse_mean),
                  z = round(z, 2))
print(tab, row.names = FALSE)
cat(sprintf("\nmax |z| %.2f; mean over trials: original %.3f, pinned %.3f h\n",
            max(abs(z)), mean(sm$original$mean), mean(sm$pinned$mean)))

cat("\n=== elpd ===\n")
L <- lapply(NM, function(nm) readRDS(sprintf("_data/wock-fit-LOO-%s.rds", nm))$observation)
dd <- L$pinned$pointwise[, "elpd_loo"] - L$original$pointwise[, "elpd_loo"]
se <- sqrt(length(dd)) * sd(dd)
cat(sprintf("Cells: summed over 1130 observations. elpd original %.2f, pinned %.2f;\n",
            L$original$estimates["elpd_loo", 1], L$pinned$estimates["elpd_loo", 1]))
cat(sprintf("pinned - original %+.2f (paired se %.2f)\n", sum(dd), se))

cat("\n=== verdict ===\n")
cat(sprintf("1. gate: %s\n2. cycle_length max |z| < 3: %s\n3. |elpd diff| < 2 se: %s\n",
            if (gate) "PASS" else "FAIL", if (max(abs(z)) < 3) "PASS" else "FAIL",
            if (abs(sum(dd)) < 2 * se) "PASS" else "FAIL"))
