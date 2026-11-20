## Follow-up to start-ml.R (2026-11-20): cycle_length at each model's
## per-trial maximum (best start; 55 and 59 at their best (alpha, beta) grid
## point) against entry 54's posterior mean, which includes the priors and the
## no_pool hierarchy on cycle_length that the maxima ignore.
##   Rscript --vanilla _scripts/start-ml-cl.R > _data/start-ml-cl-<date>.txt
.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
suppressPackageStartupMessages(library(rstan))
d <- readRDS("_data/wock-data-np_start13unit_nc384_t.rds"); TR <- attr(d, "levels")$trial
best_key <- function(m) if (m %in% c("54", "57", "58")) "none" else {
    R <- lapply(1:13, \(c) readRDS(sprintf("_data/start-ml/task-%s-%02d.rds", m, c)))
    ks <- names(R[[1]]$runs)
    ks[which.max(vapply(ks, \(k) sum(vapply(R, \(x) max(vapply(x$runs[[k]], `[[`, 0, "ll")), 0)), 0))]
}
tab <- sapply(c("54", "55", "57", "58", "59"), \(m) { k <- best_key(m)
    vapply(1:13, \(c) { r <- readRDS(sprintf("_data/start-ml/task-%s-%02d.rds", m, c))$runs[[k]]
        r[[which.max(vapply(r, `[[`, 0, "ll"))]]$cl }, 0) })
c54 <- colMeans(as.matrix(readRDS("_data/wock-fit-np_wide_total0_nc384_t.rds"), pars = "cycle_length"))
cat("cycle_length (hours) at each model's per-trial MAXIMUM likelihood (priors\n",
    "ignored; bounds 35 and 50 h), and entry 54's POSTERIOR MEAN (priors and\n",
    "cycle_length hierarchy included). Last row: mean over the 13 trials.\n\n", sep = "")
out <- cbind(round(tab, 2), e54_posterior = round(c54, 2)); rownames(out) <- TR
print(rbind(out, mean13 = round(colMeans(out), 2)))
