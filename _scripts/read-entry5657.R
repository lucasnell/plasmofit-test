## Read wockner-fit.R entries 56 and 57 (SLURM 35288, 2026-11-19) beside 54
## and 55. All four share data, n_c 384, sd_log10_total0 1, adapt_delta 0.95,
## max_treedepth 12, seed and draw counts; they differ only in the start:
##   54  per-unit b_shape and b_offset, wrapped over the whole cycle (passed)
##   55  start_dist 1: one Beta on [0, 13] h shared by all units (failed)
##   56  start_dist 2: a mode per unit, one concentration shared
##   57  start_dist 3: mode and concentration per unit
## 55-57 load the start_dist builds (55: live 13db2bd; 56-57: scratch fc7db77)
## and keep the whole inoculum circulating at t = 0, window 13 h, cap 3.
## Question (claude/start-dist-design.md, "start_dist 2 and 3"): did 55 fail
## because sharing removed the per-unit PHASE?
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   Rscript --vanilla _scripts/read-entry5657.R > _data/read-entry5657-<date>.txt

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
suppressPackageStartupMessages({library(rstan); library(posterior); library(loo)})
options(width = 150)
D <- "_data/"
CFG <- c(`54` = "np_wide_total0_nc384_t", `55` = "np_start13_nc384_t",
         `56` = "np_start13mode_nc384_t", `57` = "np_start13unit_nc384_t")
F <- lapply(CFG, \(c) readRDS(sprintf("%swock-fit-%s.rds", D, c)))
d <- readRDS(sprintf("%swock-data-%s.rds", D, CFG[["56"]]))
lv <- attr(d, "levels"); W <- d$start_window

## ---- 1. Gate and sampler health -----------------------------------------------
cat("1. GATE. max R-hat over all monitored quantities (pass < 1.05), bulk ESS\n",
    "   minimum, divergent post-warmup iterations of 4000, leapfrog steps per\n",
    "   iteration, longest chain in hours, and per-chain mean lp__ (4 chains x\n",
    "   1000 draws). lp__ is comparable within an entry, not across entries.\n\n", sep = "")
for (e in names(F)) {
    f <- F[[e]]; s <- summary(f)$summary
    rh <- s[, "Rhat"]; rh <- rh[is.finite(rh)]
    sp <- get_sampler_params(f, inc_warmup = FALSE)
    lp <- apply(extract(f, "lp__", permuted = FALSE)[, , 1, drop = FALSE], 2, mean)
    cat(sprintf("   %s: max R-hat %6.3f (%-18s) %s | ESS min %4.0f | div %3d | leapfrog %6.1f | longest chain %5.2f h | lp__ %s\n",
                e, max(rh), names(which.max(rh)), if (max(rh) < 1.05) "PASS" else "FAIL",
                min(s[, "n_eff"], na.rm = TRUE),
                sum(unlist(lapply(sp, \(x) x[, "divergent__"]))),
                mean(unlist(lapply(sp, \(x) x[, "n_leapfrog__"]))),
                max(rowSums(get_elapsed_time(f))) / 3600,
                paste(sprintf("%.1f", lp), collapse = ", ")))
}

## ---- 2. Per chain: where each chain sits --------------------------------------
cat("\n2. PER CHAIN. sum_log_lik = mean over the chain's draws of the pointwise\n",
    "   log-likelihood summed over 1130 observations (priors excluded; comparable\n",
    "   across entries; entry 54 is -807.2 in every chain). at_bound = trials whose\n",
    "   chain-mean cycle_length is within 0.5 h of 35 or 50 h. cl_mean13 = chain\n",
    "   mean of the mean cycle_length over 13 trials (hours).\n\n", sep = "")
for (e in names(F)) {
    f <- F[[e]]
    ll <- apply(extract(f, "log_lik", permuted = FALSE), 2, \(m) mean(rowSums(m)))
    cl <- extract(f, "cycle_length", permuted = FALSE)
    cm <- apply(cl, c(2, 3), mean)                       # chain x trial
    for (k in seq_along(ll)) {
        b <- which(cm[k, ] < 35.5 | cm[k, ] > 49.5)
        cat(sprintf("   %s chain %d: sum_log_lik %7.1f | cl_mean13 %5.2f | at_bound: %s\n",
                    e, k, ll[k], mean(cm[k, ]),
                    if (length(b)) paste(sprintf("%s %.1f", lv$trial[b], cm[k, b]), collapse = ", ") else "none"))
    }
}

## ---- 3. Which trials split the chains -----------------------------------------
cat("\n3. TRIALS WHOSE cycle_length SPLITS THE CHAINS (between-chain range of\n",
    "   chain means > 2 h), chain means in hours.\n\n", sep = "")
for (e in c("55", "56", "57")) {
    cm <- apply(extract(F[[e]], "cycle_length", permuted = FALSE), c(2, 3), mean)
    split <- which(apply(cm, 2, \(v) diff(range(v))) > 2)
    cat(sprintf("   %s: %d of 13 trials split", e, length(split)))
    if (length(split)) cat(":", paste(sprintf("%s [%s]", lv$trial[split],
        apply(cm[, split, drop = FALSE], 2, \(v) paste(sprintf("%.1f", v), collapse = " "))), collapse = "; "))
    cat("\n")
}

## ---- 4. The starts -------------------------------------------------------------
cat("\n4. STARTS. Per unit, posterior mean of the mode in hours (start_mode x 13),\n",
    "   pooled over all chains; ALSO the between-chain range of chain means, since\n",
    "   pooling chains that disagree hides the disagreement. Entry 54's per-unit\n",
    "   modal ages (posterior means, _data/inoculum-age-check-2026-11-09.txt) for\n",
    "   comparison.\n\n", sep = "")
m54 <- c(3.651, 10.922, 1.129, 2.278, 2.829, 10.957, 11.739, 5.322, 5.464, 6.503,
         10.891, 3.109, 1.324, 4.620)
for (e in c("56", "57")) {
    mo <- extract(F[[e]], "start_mode", permuted = FALSE) * W      # iter x chain x unit
    pooled <- apply(mo, 3, mean); rng <- apply(apply(mo, c(2, 3), mean), 2, \(v) diff(range(v)))
    cat(sprintf("   entry %s\n", e))
    tab <- data.frame(unit = lv$inoc, mode_h = round(pooled, 2), chain_range_h = round(rng, 2),
                      e54_modal_h = m54)
    if (e == "57") {
        kc <- extract(F[[e]], "start_conc", permuted = FALSE)
        tab$conc <- round(apply(kc, 3, mean), 2)
        tab$conc_chain_range <- round(apply(apply(kc, c(2, 3), mean), 2, \(v) diff(range(v))), 2)
    }
    print(tab, row.names = FALSE)
    if (e == "56") {
        kc <- as.vector(extract(F[[e]], "start_conc", permuted = FALSE)[, , 1])
        kch <- apply(extract(F[[e]], "start_conc", permuted = FALSE)[, , 1, drop = FALSE], 2, mean)
        cat(sprintf("   shared start_conc: pooled mean %.2f; chain means %s (prior flat on [2, 6])\n",
                    mean(kc), paste(sprintf("%.2f", kch), collapse = ", ")))
    }
    cat(sprintf("   Spearman, pooled mode against entry 54's modal age over 14 units: %.2f\n\n",
                cor(pooled, m54, method = "spearman")))
}

## ---- 5. LOO, for the record ------------------------------------------------------
cat("5. PSIS-LOO elpd over the same 1130 observations. NOT a model comparison for\n",
    "   any entry that fails the gate: there it averages over separate modes.\n\n", sep = "")
for (e in names(F)) {
    L <- readRDS(sprintf("%swock-fit-LOO-%s.rds", D, CFG[[e]]))$observation
    cat(sprintf("   %s: elpd_loo %7.1f (se %4.1f), p_loo %5.1f, Pareto k > 0.7: %d\n",
                e, L$estimates["elpd_loo", 1], L$estimates["elpd_loo", 2],
                L$estimates["p_loo", 1], sum(L$diagnostics$pareto_k > 0.7)))
}
