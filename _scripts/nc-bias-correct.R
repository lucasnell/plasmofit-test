## Combine the simulation 2x2 with the real-data n_c ladder.
##
## The 2x2 measures how much cycle_length bias each (generating, fitting) n_c
## pair produces against a known truth. The real-data ladder gives fitted
## cycle_lengths with no truth. Putting them together asks: under each
## hypothesis about the TRUE n_c, what is the bias-corrected real-data
## cycle_length, and do the rungs then agree with each other?
##
## Internal consistency is the test. If the true n_c is H, then correcting the
## n_c = 96 fit by the (H, 96) bias and the n_c = 192 fit by the (H, 192) bias
## should land on the SAME number, because they are two measurements of one
## quantity. The hypothesis that makes them agree is the better supported one.
##
## CAVEAT, stated because it limits what this can claim: the 2x2 arms use the
## schedule simulation's prior configuration (mean_log_b_shape 2,
## sd_log_b_shape 0.5, sd_logit_cl 1) while the real-data ladder uses the
## widened b_shape prior (sd_log_b_shape 1.5). The bias estimates are
## therefore approximate when transferred. The simulated truth, 45.012 h, is
## itself taken from an n_c = 96 fit, so if 96 inflates cycle_length the truth
## sits high too and the biases are measured in the wrong neighbourhood.
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   Rscript --vanilla _scripts/nc-bias-correct.R

suppressPackageStartupMessages({library(rstan); library(dplyr); library(tibble)})

## --- real-data ladder, b_shape estimated -------------------------------- ##
REAL <- c(np_wide_both = 96L, np_nc192 = 192L, np_nc384 = 384L)
real <- bind_rows(lapply(names(REAL), function(nm) {
    s <- summary(readRDS(sprintf("_data/wock-fit-%s.rds", nm)))$summary
    i <- grep("^cycle_length\\[", rownames(s))
    tibble(fit_n_c = REAL[[nm]], fitted_cl = mean(s[i, "mean"]))
}))

## --- 2x2 biases, converged replicates only ------------------------------ ##
f <- list.files("_data", "^wock-schedsim-RES-.*[.]rds$", full.names = TRUE)
CELLS <- c("default", "sim96_fit192", "sim192_fit96", "sim192_fit192")
sim <- bind_rows(lapply(f, function(p) {
    x <- readRDS(p)
    if (!x$arm %in% CELLS) return(NULL)
    if (x$config != sprintf("%s-rep%d", x$arm, x$rep)) return(NULL)
    ## `default` has replicates 1-7 while the new arms have only 1-3. Using
    ## all of default's would compare cells built on different replicate
    ## sets, so the contrast is restricted to the shared 1-3.
    if (x$rep > 3) return(NULL)
    if (x$max_rhat >= 1.05) return(NULL)
    tibble(gen_n_c = as.integer(x$sim_n_c %||% 96L),
           fit_n_c = as.integer(x$fit_n_c %||% 96L),
           bias = mean(x$per_trial$cl_mean) - x$true_cl)
})) |> group_by(gen_n_c, fit_n_c) |>
    summarise(n = n(), bias = mean(bias), .groups = "drop")

cat("=== measured bias by cell (converged only) ===\n")
cat("Cells: mean posterior cycle_length minus the 45.012 h simulated truth.\n\n")
print(as.data.frame(sim), digits = 4)

cat("\n=== bias-corrected real-data cycle_length, by hypothesis ===\n")
cat("Cells: the real fit at that n_c, minus the bias the 2x2 measured for\n",
    "(true n_c = H, fitted n_c = that rung). Two rungs corrected under the\n",
    "same H are two measurements of one quantity and should agree.\n\n", sep = "")
out <- bind_rows(lapply(c(96L, 192L), function(H) {
    bind_rows(lapply(c(96L, 192L), function(FN) {
        b <- sim$bias[sim$gen_n_c == H & sim$fit_n_c == FN]
        if (!length(b)) return(NULL)
        tibble(true_n_c = H, rung = FN,
               fitted = real$fitted_cl[real$fit_n_c == FN],
               bias = b, corrected = real$fitted_cl[real$fit_n_c == FN] - b)
    }))
}))
print(as.data.frame(out), digits = 4)

cat("\n=== internal consistency ===\n")
for (H in unique(out$true_n_c)) {
    v <- out$corrected[out$true_n_c == H]
    cat(sprintf("true n_c = %3d: corrected estimates %s | spread %.2f h\n",
                H, paste(sprintf("%.2f", v), collapse = ", "),
                max(v) - min(v)))
}
cat("\nThe hypothesis with the SMALLER spread reconciles the two rungs better.\n",
    "Read it together with the elpd ordering, not instead of it, and against\n",
    "the caveats at the top of this file.\n", sep = "")
