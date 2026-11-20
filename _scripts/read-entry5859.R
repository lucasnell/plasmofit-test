## Read wockner-fit.R entries 58 and 59 beside 54, 55 and 57. All share data,
## no_pool, n_c 384, sd_log10_total0 1, adapt_delta 0.95, max_treedepth 12,
## seed and draw counts (4 x 1000); they differ only in the start:
##   54  per-unit b_shape and b_offset, wrapped over the whole cycle (passed)
##   55  start_dist 1: one Beta on [0, 13] h shared by all units (failed)
##   57  start_dist 3: mode and concentration per unit on [0, 13] h (failed)
##   58  57 with the window reaching below age 0: [-6, 13] h, ages below 0
##       being late schizonts that burst soon after injection
##   59  55 plus a per-unit post-injection delay, flat on [0, 12] h:
##       y_obs(t) = y_model(t - delay)
## 55-59 keep the whole inoculum circulating at t = 0, cap 3. 58-59 need the
## plasmofit build with start_below/start_delay (claude/start-dist-design.md,
## "Below age 0 and a delay").
## Question: do 55-57 fail because the 13 h window cannot reach the phase the
## data want (entry 54's extra start mass sits in the cycle's last 6 h)?
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   Rscript --vanilla _scripts/read-entry5859.R > _data/read-entry5859-<date>.txt

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
suppressPackageStartupMessages({library(rstan); library(posterior); library(loo)})
options(width = 150)
D <- "_data/"
CFG <- c(`54` = "np_wide_total0_nc384_t", `55` = "np_start13_nc384_t",
         `57` = "np_start13unit_nc384_t", `58` = "np_start13below6unit_nc384_t",
         `59` = "np_start13delay12_nc384_t")
F <- lapply(CFG, \(c) readRDS(sprintf("%swock-fit-%s.rds", D, c)))
d58 <- readRDS(sprintf("%swock-data-%s.rds", D, CFG[["58"]]))
d59 <- readRDS(sprintf("%swock-data-%s.rds", D, CFG[["59"]]))
lv <- attr(d58, "levels"); stopifnot(identical(lv, attr(d59, "levels")))
W <- d58$start_window; B <- d58$start_below; DMAX <- d59$delay_max
pass <- c()

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
    pass[e] <- max(rh) < 1.05
    cat(sprintf("   %s: max R-hat %6.3f (%-18s) %s | ESS min %4.0f | div %3d | leapfrog %6.1f | longest chain %5.2f h | lp__ %s\n",
                e, max(rh), names(which.max(rh)), if (pass[e]) "PASS" else "FAIL",
                min(s[, "n_eff"], na.rm = TRUE),
                sum(unlist(lapply(sp, \(x) x[, "divergent__"]))),
                mean(unlist(lapply(sp, \(x) x[, "n_leapfrog__"]))),
                max(rowSums(get_elapsed_time(f))) / 3600,
                paste(sprintf("%.1f", lp), collapse = ", ")))
}

## ---- 2. Per chain ------------------------------------------------------------
cat("\n2. PER CHAIN. sum_log_lik = mean over the chain's draws of the pointwise\n",
    "   log-likelihood summed over 1130 observations (priors excluded; comparable\n",
    "   across entries; entry 54 is -807.2 in every chain; 57's best chains -825).\n",
    "   at_bound = trials whose chain-mean cycle_length is within 0.5 h of 35 or\n",
    "   50 h. cl_mean13 = chain mean of the mean cycle_length over 13 trials (h).\n\n", sep = "")
for (e in c("58", "59")) {
    f <- F[[e]]
    ll <- apply(extract(f, "log_lik", permuted = FALSE), 2, \(m) mean(rowSums(m)))
    cm <- apply(extract(f, "cycle_length", permuted = FALSE), c(2, 3), mean)
    for (k in seq_along(ll)) {
        b <- which(cm[k, ] < 35.5 | cm[k, ] > 49.5)
        cat(sprintf("   %s chain %d: sum_log_lik %7.1f | cl_mean13 %5.2f | at_bound: %s\n",
                    e, k, ll[k], mean(cm[k, ]),
                    if (length(b)) paste(sprintf("%s %.1f", lv$trial[b], cm[k, b]), collapse = ", ") else "none"))
    }
}

## ---- 3. Trials that split the chains --------------------------------------------
cat("\n3. TRIALS WHOSE cycle_length SPLITS THE CHAINS (between-chain range of\n",
    "   chain means > 2 h), chain means in hours.\n\n", sep = "")
for (e in c("58", "59")) {
    cm <- apply(extract(F[[e]], "cycle_length", permuted = FALSE), c(2, 3), mean)
    split <- which(apply(cm, 2, \(v) diff(range(v))) > 2)
    cat(sprintf("   %s: %d of 13 trials split", e, length(split)))
    if (length(split)) cat(":", paste(sprintf("%s [%s]", lv$trial[split],
        apply(cm[, split, drop = FALSE], 2, \(v) paste(sprintf("%.1f", v), collapse = " "))), collapse = "; "))
    cat("\n")
}

## ---- 4. Starts and delays ---------------------------------------------------------
m54 <- c(3.651, 10.922, 1.129, 2.278, 2.829, 10.957, 11.739, 5.322, 5.464, 6.503,
         10.891, 3.109, 1.324, 4.620)
cat("\n4a. ENTRY 58 STARTS. Per unit, pooled over all chains: posterior mean of the\n",
    "   modal age in hours (start_mode x ", W + B, " - ", B, "; negative = before bursting),\n",
    "   the between-chain range of chain means of that mode (h), the posterior\n",
    "   mean share of the start below age 0 (pbeta(", B, "/", W + B, ", a, b) per draw), and\n",
    "   start_conc (prior flat on [2, ", 2 * d58$start_shape_cap, "]). e54_modal_h: entry 54's\n",
    "   per-unit modal ages (posterior means, _data/inoculum-age-check-2026-11-09.txt).\n\n", sep = "")
mo <- extract(F[["58"]], "start_mode", permuted = FALSE) * (W + B) - B
a <- extract(F[["58"]], "start_a", permuted = FALSE); b <- extract(F[["58"]], "start_b", permuted = FALSE)
below <- array(pbeta(B / (W + B), a, b), dim(a))
kc <- extract(F[["58"]], "start_conc", permuted = FALSE)
print(data.frame(unit = lv$inoc, mode_h = round(apply(mo, 3, mean), 2),
                 mode_chain_range_h = round(apply(apply(mo, c(2, 3), mean), 2, \(v) diff(range(v))), 2),
                 share_below0 = round(apply(below, 3, mean), 3),
                 conc = round(apply(kc, 3, mean), 2), e54_modal_h = m54), row.names = FALSE)

cat("\n4b. ENTRY 59 DELAYS AND SHARED START. Per unit, pooled over chains: posterior\n",
    "   mean and 90% interval of the delay (hours, prior flat on [0, ", DMAX, "]), the\n",
    "   between-chain range of chain means (h), and the share of draws within\n",
    "   0.25 h of either bound.\n\n", sep = "")
dl <- extract(F[["59"]], "delay", permuted = FALSE)
q <- apply(dl, 3, \(v) quantile(v, c(0.05, 0.95)))
print(data.frame(unit = lv$inoc, delay_h = round(apply(dl, 3, mean), 2),
                 q05 = round(q[1, ], 2), q95 = round(q[2, ], 2),
                 chain_range_h = round(apply(apply(dl, c(2, 3), mean), 2, \(v) diff(range(v))), 2),
                 at_0 = round(apply(dl, 3, \(v) mean(v < 0.25)), 2),
                 at_max = round(apply(dl, 3, \(v) mean(v > DMAX - 0.25)), 2)), row.names = FALSE)
al <- as.vector(extract(F[["59"]], "start_alpha", permuted = FALSE))
be <- as.vector(extract(F[["59"]], "start_beta", permuted = FALSE))
cat(sprintf("\n   shared start: start_alpha %.2f [%.2f, %.2f], start_beta %.2f [%.2f, %.2f] (posterior mean [90%%]; flat on [1, %g])\n",
            mean(al), quantile(al, 0.05), quantile(al, 0.95),
            mean(be), quantile(be, 0.05), quantile(be, 0.95), d59$start_shape_cap))
cat(sprintf("   Spearman, unit delay against entry 54's modal age over 14 units: %.2f\n",
            cor(apply(dl, 3, mean), m54, method = "spearman")))

## ---- 5. cycle_length, only for entries that pass ------------------------------------
cat("\n5. CYCLE LENGTH, entries that pass the gate only. Per trial: posterior mean\n",
    "   (hours); 'mean13' = posterior mean of the per-draw mean over 13 trials;\n",
    "   diff = this entry minus entry 54, difference of posterior means, unpaired.\n\n", sep = "")
c54 <- as.matrix(F[["54"]], pars = "cycle_length")
for (e in c("58", "59")) {
    if (!pass[e]) { cat(sprintf("   %s: fails the gate; not reported\n", e)); next }
    ce <- as.matrix(F[[e]], pars = "cycle_length")
    print(data.frame(trial = lv$trial, mean = round(colMeans(ce), 2), sd = round(apply(ce, 2, sd), 2),
                     e54 = round(colMeans(c54), 2), diff = round(colMeans(ce) - colMeans(c54), 2)),
          row.names = FALSE)
    cat(sprintf("   %s mean13 %.2f h (sd %.2f); entry 54 %.2f h (sd %.2f)\n\n", e,
                mean(rowMeans(ce)), sd(rowMeans(ce)), mean(rowMeans(c54)), sd(rowMeans(c54))))
}

## ---- 6. LOO ----------------------------------------------------------------------------
cat("6. PSIS-LOO elpd over the same 1130 observations. NOT a model comparison for\n",
    "   any entry that fails the gate: there it averages over separate modes.\n\n", sep = "")
for (e in names(F)) {
    L <- readRDS(sprintf("%swock-fit-LOO-%s.rds", D, CFG[[e]]))$observation
    cat(sprintf("   %s: elpd_loo %7.1f (se %4.1f), p_loo %5.1f, Pareto k > 0.7: %d%s\n",
                e, L$estimates["elpd_loo", 1], L$estimates["elpd_loo", 2],
                L$estimates["p_loo", 1], sum(L$diagnostics$pareto_k > 0.7),
                if (pass[e]) "" else "  (fails the gate)"))
}
