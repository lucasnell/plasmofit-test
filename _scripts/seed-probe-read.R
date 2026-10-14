## Read the seed probe (2026-10-13): three new seeds each of entry 43 (pin 15,
## default sampler settings) and entry 54 (free b_shape, tuned settings), run
## with the pinned driver (_scripts/wockner-fit-pinned.sh, PIN_PLAN). The
## question: how often do chains reach the low region of MMV048_PIB's
## cycle_length that one chain of the pinned pin-15 refit (30856) found, and
## does anything similar happen in other trials?
##
## Definitions, FIXED BEFORE THE RUNS (from the two existing pin-15 fits,
## 30830_43 and 30856, whose main-region 5% quantile for MMV048_PIB was
## 41.7-41.9 h in all 8 chains and whose excursion averaged 39.9 h):
##   - LOW = cycle_length[7] < 41.0 h.
##   - A chain VISITS the low region if >= 5% of its post-warmup draws are LOW.
##   - Other trials: a chain is an OUTLIER for a trial if its posterior mean
##     differs by > 1 h from the median of the other chains' means in that run.
## Reported, not gated: the gate per run; visits per chain; pooled share of
## LOW draws per config (a visit frequency, NOT a posterior mass, since
## chains that do not move between regions cannot weigh them); per-trial
## cycle_length across runs that pass the gate.
## The two earlier runs (30830_43 original, 30856 pinned) are read alongside.
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   Rscript --vanilla _scripts/seed-probe-read.R > _data/seed-probe-read-<date>.txt

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
suppressPackageStartupMessages(library(rstan))
LOW <- 41.0; VISIT <- 0.05; OUTLIER_H <- 1
SEEDS <- c(20261013, 20261014, 20261015)
RUNS <- rbind(
    data.frame(config = "pin15", name = c("np_bs15_nc384", "np_bs15_nc384-pinned",
               sprintf("np_bs15_nc384-pinned-s%d", SEEDS))),
    data.frame(config = "free", name = c("np_wide_total0_nc384_t",
               sprintf("np_wide_total0_nc384_t-pinned-s%d", SEEDS))))
RUNS <- RUNS[file.exists(sprintf("_data/wock-fit-%s.rds", RUNS$name)), ]
cat("runs present:", nrow(RUNS), "\n\n")

trials <- attr(readRDS(sprintf("_data/wock-data-%s.rds", RUNS$name[1])), "levels")[["trial"]]
stopifnot(trials[7] == "MMV048_PIB")

out <- list(); pass <- c(); clm <- list(); outl <- list()
for (r in seq_len(nrow(RUNS))) {
    nm <- RUNS$name[r]; f <- readRDS(sprintf("_data/wock-fit-%s.rds", nm))
    s <- summary(f)$summary; mx <- max(s[, "Rhat"], na.rm = TRUE); pass[nm] <- mx < 1.05
    lp <- apply(extract(f, "lp__", permuted = FALSE)[, , 1, drop = FALSE], 2, mean)
    cl <- extract(f, "cycle_length", permuted = FALSE)        # iter x chain x 13
    x7 <- cl[, , 7]
    out[[nm]] <- data.frame(config = RUNS$config[r], run = nm,
        max_rhat = round(mx, 3), gate = if (mx < 1.05) "PASS" else "FAIL",
        chain = seq_len(ncol(x7)), lp = round(lp, 1),
        mean_cl7 = round(colMeans(x7), 2), low_pct = round(100 * colMeans(x7 < LOW), 1),
        visits = colMeans(x7 < LOW) >= VISIT)
    cm <- apply(cl, c(2, 3), mean)                            # chain x 13
    for (t in seq_along(trials)) for (ch in seq_len(nrow(cm))) {
        dev <- cm[ch, t] - median(cm[-ch, t])
        if (abs(dev) > OUTLIER_H) outl[[length(outl) + 1]] <- data.frame(
            run = nm, trial = trials[t], chain = ch, chain_mean = round(cm[ch, t], 2),
            others_median = round(median(cm[-ch, t]), 2))
    }
    clm[[nm]] <- colMeans(apply(cl, 3, c))
}
O <- do.call(rbind, out); rownames(O) <- NULL

cat("=== MMV048_PIB by chain ===\n")
cat("Cells: per chain. lp = mean lp__; mean_cl7 = posterior mean of MMV048_PIB's\n",
    "cycle_length (h); low_pct = % of 1000 post-warmup draws below ", LOW, " h;\n",
    "visits = low_pct >= ", 100 * VISIT, "%. max_rhat and gate are per run.\n\n", sep = "")
print(O, row.names = FALSE)

cat("\n=== visit frequency by config ===\n")
cat("Cells: chains that visit / chains; pooled % of LOW draws over all chains.\n",
    "A frequency of visits, not a posterior mass.\n\n", sep = "")
for (cf in unique(O$config)) {
    x <- O[O$config == cf, ]
    cat(sprintf("%-6s %d of %d chains visit, in %d of %d runs; pooled LOW %.1f%%\n", cf,
                sum(x$visits), nrow(x), length(unique(x$run[x$visits])),
                length(unique(x$run)), mean(x$low_pct)))
}

cat("\n=== chains more than ", OUTLIER_H, " h from the others, any trial ===\n", sep = "")
if (length(outl)) print(do.call(rbind, outl), row.names = FALSE) else cat("none\n")

cat("\n=== cycle_length by trial, runs that pass the gate ===\n")
cat("Cells: posterior mean in hours per trial; last row the mean over 13 trials.\n\n")
ok <- names(pass)[pass]
if (length(ok)) {
    M <- sapply(ok, function(nm) clm[[nm]]); rownames(M) <- trials
    M <- rbind(M, mean = colMeans(M)); print(round(M, 2))
} else cat("no run passes\n")

## Added 2026-10-14 after the first partial read, descriptive only (no
## definition above changed): which quantities fail the gate in each run, by
## the parameter's index, and per-chain means of the two MMV048 trials.
cat("\n=== worst R-hats per run (added after the first read) ===\n")
cat("Cells: the 6 quantities with the highest R-hat in each run, and how many\n",
    "of the run's quantities exceed 1.05.\n\n", sep = "")
for (nm in RUNS$name) {
    f <- readRDS(sprintf("_data/wock-fit-%s.rds", nm))
    r <- sort(summary(f)$summary[, "Rhat"], decreasing = TRUE)
    cat(sprintf("%s (%d > 1.05): %s\n", nm, sum(r > 1.05, na.rm = TRUE),
                paste(sprintf("%s %.3f", names(r)[1:6], r[1:6]), collapse = "; ")))
    cl <- extract(f, "cycle_length", permuted = FALSE)
    cat(sprintf("   chain means, MMV048_PIB: %s | MMV048_PartB: %s\n",
                paste(sprintf("%.2f", colMeans(cl[, , 7])), collapse = " "),
                paste(sprintf("%.2f", colMeans(cl[, , 13])), collapse = " ")))
}
