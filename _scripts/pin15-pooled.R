## Pool the five pin-15 runs at n_c = 384 (entry 43) into one 20-chain set and
## recompute the diagnostics. The runs are independent chains on the same
## target (same model, data, and sampler settings; different seeds or
## layouts), so pooling is legitimate. It cannot fix a failure to move between
## regions: a region's share of pooled draws reflects how many chains happened
## to find it, not its posterior mass.
##
## Runs: 30830_43 (seed 20261010, default layout), 30856 (seed 20261010,
## pinned), 30858 tasks 1-3 (seeds 20261013-15, pinned).
##
## Cells:
##   - max R-hat and min bulk ESS: over every monitored quantity, 20 chains
##     (rank-normalised split R-hat, posterior::rhat).
##   - per trial: pooled posterior mean and sd of cycle_length (h), R-hat and
##     bulk ESS over the 20 chains, the range of the five per-run means, and
##     the range of the 20 per-chain means.
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   Rscript --vanilla _scripts/pin15-pooled.R > _data/pin15-pooled-<date>.txt

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
suppressPackageStartupMessages({library(rstan); library(posterior)})

RUNS <- c("np_bs15_nc384", "np_bs15_nc384-pinned",
          sprintf("np_bs15_nc384-pinned-s%d", 20261013:20261015))
dat <- lapply(RUNS, function(nm) readRDS(sprintf("_data/wock-data-%s.rds", nm)))
for (x in dat[-1]) stopifnot(identical(x, dat[[1]]))
trials <- attr(dat[[1]], "levels")[["trial"]]

arr <- lapply(RUNS, function(nm) {
    f <- readRDS(sprintf("_data/wock-fit-%s.rds", nm))
    as_draws_array(extract(f, permuted = FALSE))
})
pooled <- do.call(bind_draws, c(arr, along = "chain"))
cat(sprintf("pooled: %d chains x %d draws, %d quantities\n\n",
            nchains(pooled), niterations(pooled), nvariables(pooled)))

## b_shape is data here, so its draws are constant: drop zero-variance
## quantities before computing R-hat.
v <- variables(pooled)
keep <- v[vapply(v, function(nm) sd(as.vector(pooled[, , nm])) > 0, logical(1))]
dg <- summarise_draws(subset_draws(pooled, variable = keep), rhat, ess_bulk)
cat(sprintf("max R-hat %.3f over %d quantities; %d above 1.05; min bulk ESS %.0f\n",
            max(dg$rhat, na.rm = TRUE), nrow(dg), sum(dg$rhat > 1.05, na.rm = TRUE),
            min(dg$ess_bulk, na.rm = TRUE)))
w <- dg[order(-dg$rhat), ][1:10, ]
cat("ten worst:", paste(sprintf("%s %.3f", w$variable, w$rhat), collapse = "; "), "\n\n")

cl <- subset_draws(pooled, variable = "cycle_length")
s <- summarise_draws(cl, mean, sd, rhat, ess_bulk)
run_means <- sapply(arr, function(a) colMeans(apply(as.array(subset_draws(a, variable = "cycle_length")), 3, c)))
chain_means <- apply(as.array(cl), c(2, 3), mean)            # chain x 13
tab <- data.frame(trial = trials, mean = round(s$mean, 2), sd = round(s$sd, 2),
                  rhat = round(s$rhat, 3), ess_bulk = round(s$ess_bulk),
                  run_means = apply(run_means, 1, function(x) sprintf("%.2f-%.2f", min(x), max(x))),
                  chain_means = apply(chain_means, 2, function(x) sprintf("%.2f-%.2f", min(x), max(x))))
print(tab, row.names = FALSE)
tm <- rowMeans(apply(as.array(cl), 3, c))
cat(sprintf("\nmean over 13 trials: pooled posterior mean %.3f h (sd %.3f); per-run means %s\n",
            mean(tm), sd(tm), paste(sprintf("%.3f", colMeans(run_means)), collapse = ", ")))
