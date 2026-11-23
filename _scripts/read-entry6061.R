## Read wockner-fit.R entries 60 and 61 beside 54 and 59. 2026-11-23.
## All four are no_pool, n_c 384, sd_log10_total0 1, adapt_delta 0.95,
## max_treedepth 12, seed 538065874, 4 x 1000 draws. 59, 60 and 61 share the
## start: one Beta on [0, 13] h for every unit plus a per-unit delay flat on
## [0, 12] h, whole inoculum circulating at t = 0.
##   54  per-unit wrapped Beta over the whole cycle (the only entry that passes)
##   59  the delay model, centred cycle_length hierarchy, default inits
##   60  59 with the NON-CENTRED cycle_length hierarchy
##   61  59 with chains STARTED AT its maximum-likelihood point
##       (-766.28 summed over 1130 observations, profiled sds;
##       _scripts/start-ml-init.R). DIAGNOSTIC: chains started in one place
##       can agree without exploring, so a passing gate here does NOT license
##       reporting this posterior.
## Question: is entry 59's ~87-unit gap below its own maximum the SAMPLER or
## the MODEL? 60 asks whether the hierarchy's parameterisation traps chains;
## 61 whether the good region is a posterior mode at all.
##   Rscript --vanilla _scripts/read-entry6061.R > _data/read-entry6061-<date>.txt

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
suppressPackageStartupMessages({library(rstan); library(posterior)})
options(width = 150)
D <- "_data/"
ML59 <- -766.28
CFG <- c(`54` = "np_wide_total0_nc384_t", `59` = "np_start13delay12_nc384_t",
         `60` = "np_start13delay12_ncp_nc384_t", `61` = "np_start13delay12_mlinit_nc384_t")
have <- vapply(CFG, \(c) file.exists(sprintf("%swock-fit-%s.rds", D, c)), TRUE)
if (any(!have)) cat("NOT YET AVAILABLE:", paste(names(CFG)[!have], collapse = ", "), "\n\n")
CFG <- CFG[have]
F <- lapply(CFG, \(c) readRDS(sprintf("%swock-fit-%s.rds", D, c)))
lv <- attr(readRDS(sprintf("%swock-data-%s.rds", D, CFG[[1]])), "levels")

cat("1. GATE. max R-hat over all monitored quantities (pass < 1.05), bulk ESS\n",
    "   minimum, divergent post-warmup iterations of 4000, longest chain in hours.\n\n", sep = "")
for (e in names(F)) {
    f <- F[[e]]; s <- summary(f)$summary; rh <- s[, "Rhat"]; rh <- rh[is.finite(rh)]
    sp <- get_sampler_params(f, inc_warmup = FALSE)
    cat(sprintf("   %s: max R-hat %7.3f (%-18s) %s | ESS min %4.0f | div %3d | longest chain %5.2f h\n",
                e, max(rh), names(which.max(rh)), if (max(rh) < 1.05) "PASS" else "FAIL",
                min(s[, "n_eff"], na.rm = TRUE),
                sum(unlist(lapply(sp, \(x) x[, "divergent__"]))),
                max(rowSums(get_elapsed_time(f))) / 3600))
}

cat("\n2. PER CHAIN. sum_log_lik = the pointwise log-likelihood summed over 1130\n",
    "   observations (priors excluded), as a chain mean and as the best single\n",
    "   draw. For reference: entry 54 means -807, entry 59's best draw -860.7,\n",
    "   and the multistart optimiser's maximum for this model ", sprintf("%.2f", ML59), "\n",
    "   (profiled sds, so the draws cannot exceed it). sigma = chain mean of\n",
    "   sigma_logit_cl, the between-trial spread of cycle_length on the logit\n",
    "   scale. at_bound = trials whose chain-mean cycle_length is within 0.5 h\n",
    "   of 35 or 50 h.\n\n", sep = "")
for (e in names(F)) {
    f <- F[[e]]
    ll <- apply(extract(f, "log_lik", permuted = FALSE), c(1, 2), sum)
    sg <- apply(extract(f, "sigma_logit_cl", permuted = FALSE)[, , 1], 2, mean)
    cm <- apply(extract(f, "cycle_length", permuted = FALSE), c(2, 3), mean)
    for (k in seq_len(ncol(ll))) {
        b <- which(cm[k, ] < 35.5 | cm[k, ] > 49.5)
        cat(sprintf("   %s chain %d: mean %8.1f | best %8.1f | gap to max %6.1f | sigma %5.2f | cl_mean13 %5.2f | at_bound: %s\n",
                    e, k, mean(ll[, k]), max(ll[, k]), ML59 - max(ll[, k]), sg[k], mean(cm[k, ]),
                    if (length(b)) paste(sprintf("%s %.1f", lv$trial[b], cm[k, b]), collapse = ", ") else "none"))
    }
}

cat("\n3. TRIALS WHOSE cycle_length SPLITS THE CHAINS (between-chain range of chain\n",
    "   means > 2 h), chain means in hours.\n\n", sep = "")
for (e in names(F)) {
    cm <- apply(extract(F[[e]], "cycle_length", permuted = FALSE), c(2, 3), mean)
    sp <- which(apply(cm, 2, \(v) diff(range(v))) > 2)
    cat(sprintf("   %s: %d of 13 trials split", e, length(sp)))
    if (length(sp)) cat(":", paste(sprintf("%s [%s]", lv$trial[sp],
        apply(cm[, sp, drop = FALSE], 2, \(v) paste(sprintf("%.1f", v), collapse = " "))), collapse = "; "))
    cat("\n")
}

cat("\n4. DELAYS, for entries that have them. Per unit, pooled over chains:\n",
    "   posterior mean (hours, prior flat on [0, 12]) and the between-chain range\n",
    "   of chain means. ml = the maximum-likelihood delay (_data/start-ml-init-59.rds).\n\n", sep = "")
mlv <- readRDS(paste0(D, "start-ml-init-59.rds"))
tab <- data.frame(unit = lv$inoc, ml = round(mlv$delay, 2))
for (e in intersect(names(F), c("59", "60", "61"))) {
    dl <- extract(F[[e]], "delay", permuted = FALSE)
    tab[[paste0("e", e)]] <- round(apply(dl, 3, mean), 2)
    tab[[paste0("e", e, "_rng")]] <- round(apply(apply(dl, c(2, 3), mean), 2, \(v) diff(range(v))), 2)
}
print(tab, row.names = FALSE)

cat("\n5. CYCLE LENGTH, only for entries that PASS the gate. Posterior mean per\n",
    "   trial (hours); e54 = entry 54's posterior mean; ml = the maximum-likelihood\n",
    "   value. mean13 = posterior mean of the per-draw mean over 13 trials.\n\n", sep = "")
pass <- vapply(F, \(f) { rh <- summary(f)$summary[, "Rhat"]; max(rh[is.finite(rh)]) < 1.05 }, TRUE)
ct <- data.frame(trial = lv$trial, ml = round(mlv$cycle_length, 2))
if ("54" %in% names(F)) ct$e54 <- round(colMeans(as.matrix(F[["54"]], pars = "cycle_length")), 2)
for (e in setdiff(names(F)[pass], "54")) ct[[paste0("e", e)]] <- round(colMeans(as.matrix(F[[e]], pars = "cycle_length")), 2)
print(ct, row.names = FALSE)
for (e in names(F)[pass]) {
    m <- rowMeans(as.matrix(F[[e]], pars = "cycle_length"))
    cat(sprintf("   %s mean13 %.2f h (sd %.2f)\n", e, mean(m), sd(m)))
}
if (!all(pass)) cat("   not reported (fails the gate):", paste(names(F)[!pass], collapse = ", "), "\n")
