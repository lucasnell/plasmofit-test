## Follow-up to read-entry5859.R (2026-11-20): the chain-MEAN summed
## log-likelihood penalises models with more parameters (roughly p/2 below the
## best-fitting point when the posterior is near-normal), so across entries
## with different parameter counts compare the MAXIMUM over each chain's 1000
## draws as well. Entries 54-59, same 1130 observations, priors excluded.
## Nesting: 59 contains 55 (delay -> 0, checked in start-dist-validate.R 7e);
## 57 contains 55 (any Beta with shapes in [1, 3] has mode in [0, 1] and
## a + b in [2, 6]); 58 contains neither (its Beta always spans [-6, 13] h).
##   Rscript --vanilla _scripts/read-entry5859-maxll.R > _data/read-entry5859-maxll-<date>.txt
.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
suppressPackageStartupMessages(library(rstan))
CFG <- c(`54` = "np_wide_total0_nc384_t", `55` = "np_start13_nc384_t",
         `56` = "np_start13mode_nc384_t", `57` = "np_start13unit_nc384_t",
         `58` = "np_start13below6unit_nc384_t", `59` = "np_start13delay12_nc384_t")
cat("Per chain: mean and max over 1000 draws of the summed pointwise log-likelihood\n",
    "(1130 observations, priors excluded), chains 1-4 in order.\n\n", sep = "")
for (e in names(CFG)) {
    f <- readRDS(sprintf("_data/wock-fit-%s.rds", CFG[[e]]))
    s <- apply(extract(f, "log_lik", permuted = FALSE), c(1, 2), sum)   # iter x chain
    cat(sprintf("   %s: mean %s | max %s\n", e,
                paste(sprintf("%7.1f", colMeans(s)), collapse = " "),
                paste(sprintf("%7.1f", apply(s, 2, max)), collapse = " ")))
}
