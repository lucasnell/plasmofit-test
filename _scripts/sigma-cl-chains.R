## Entries 54-59 (2026-11-20): the cycle_length hierarchy's settings, and per
## chain the posterior mean of sigma_logit_cl (between-trial sd of
## logit-scale cycle_length) and mu_logit_cl, to see whether chains with
## trials stuck at the 35/50 h bounds also carry an inflated sigma.
##   Rscript --vanilla _scripts/sigma-cl-chains.R > _data/sigma-cl-chains-<date>.txt
.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
suppressPackageStartupMessages(library(rstan))
CFG <- c(`54` = "np_wide_total0_nc384_t", `55` = "np_start13_nc384_t",
         `56` = "np_start13mode_nc384_t", `57` = "np_start13unit_nc384_t",
         `58` = "np_start13below6unit_nc384_t", `59` = "np_start13delay12_nc384_t")
d <- readRDS("_data/wock-data-np_start13delay12_nc384_t.rds")
cat(sprintf("Hierarchy settings (same data list fields in all six entries): center_cl %d, n_grp_cl %d,\n",
            d$center_cl, d$n_grp_cl))
cat(sprintf("  mu_logit_cl ~ normal(%.3f, %.3f), sigma_logit_cl ~ half-normal(0, %.3f); cycle_length = 35 + 15 inv_logit(logit_cl)\n",
            d$mean_logit_cl, d$sd_logit_cl, d$sd_bs_cl))
cat(sprintf("  prior centre: %.2f h\n\n", 35 + 15 * plogis(d$mean_logit_cl)))
for (e in names(CFG)) {
    dd <- readRDS(sprintf("_data/wock-data-%s.rds", CFG[[e]]))
    stopifnot(dd$center_cl == d$center_cl, dd$mean_logit_cl == d$mean_logit_cl,
              dd$sd_logit_cl == d$sd_logit_cl, dd$sd_bs_cl == d$sd_bs_cl)
}
cat("Per chain (4 x 1000 draws): posterior mean of sigma_logit_cl and of mu_logit_cl\n",
    "converted to hours (35 + 15 inv_logit), and the number of trials whose chain-mean\n",
    "cycle_length is within 0.5 h of 35 or 50 h.\n\n", sep = "")
for (e in names(CFG)) {
    f <- readRDS(sprintf("_data/wock-fit-%s.rds", CFG[[e]]))
    sg <- apply(extract(f, "sigma_logit_cl", permuted = FALSE)[, , 1], 2, mean)
    mu <- apply(extract(f, "mu_logit_cl", permuted = FALSE)[, , 1], 2, \(v) mean(35 + 15 * plogis(v)))
    cm <- apply(extract(f, "cycle_length", permuted = FALSE), c(2, 3), mean)
    nb <- rowSums(cm < 35.5 | cm > 49.5)
    cat(sprintf("   %s: sigma %s | mu (h) %s | trials at a bound %s\n", e,
                paste(sprintf("%5.2f", sg), collapse = " "), paste(sprintf("%5.2f", mu), collapse = " "),
                paste(nb, collapse = " ")))
}
