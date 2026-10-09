## Does the forward map lose precision at high n_c?
##
## _scripts/forward-map-scaling.R found that log_prob under the two forward
## maps, which agree to 5e-15 at n_c = 96, drift apart as n_c rises: 1.0e-11 at
## 1536 and 8.3e-10 at 6144. Something is losing accuracy. This matters beyond
## the switch, because the ML dispersion profile prefers n_c >= 384 and still
## gains slightly to 1024, so high-n_c fits are a live option.
##
## mat_exp_series is the independent arbiter: it shares no code with either the
## Erlang-window series or the convolution. If conv_series still matches it at
## high n_c, the convolution is fine and the production series is the one
## drifting; if both drift, the comparison is uninformative and more work is
## needed. mat_exp_series is cubic, so this is one trajectory per rung and the
## top rung costs a minute or so.
##
##
## NOTE, 2026-10-09: conv_series() has been REMOVED from the package (revert
## commit bba580f) because it is slower and less accurate than the
## Erlang-window series at every n_c this project uses. This script therefore
## needs a plasmofit built from commit 05c9c1f to run. It is kept because the
## numbers in claude/findings.md come from it and a claim should name the
## script that produced it; the saved output beside it is the record.
##   cd /home2/lan68/plasmofit/plasmofit-test
##   PLASMOFIT_LIB=/home2/lan68/plasmofit/.Rlib-dev \
##     Rscript --vanilla _scripts/high-nc-accuracy.R \
##     | tee _data/high-nc-accuracy-$(date +%F).txt

LIB <- Sys.getenv("PLASMOFIT_LIB", "/home2/lan68/plasmofit/.Rlib-dev")
if (nzchar(LIB)) .libPaths(c(LIB, .libPaths()))
suppressPackageStartupMessages(library(plasmofit))
cat("plasmofit from:", dirname(find.package("plasmofit")), "\n")

TS <- c(72, 96, 120, 144, 168, 192, 216)
CL <- 45; BS <- 400; BO <- 0.3; LT0 <- 5.5; RR <- 8; MU <- 0

cat("\nCells: max over the 7 observation times of |conv - mat_exp| / |mat_exp|,\n",
    "at cycle_length 45, b_shape 400, R 8, mu 0. mat_exp_series shares no code\n",
    "with conv_series, so this isolates the convolution. secs are for one\n",
    "trajectory each; mat_exp is cubic in n_c, which is why the grid stops.\n\n",
    sep = "")
res <- NULL
for (n_c in c(192L, 384L, 768L, 1536L)) {
    M <- plasmofit:::conv_min_M(n_c, max(TS), 35)
    y0 <- plasmofit:::generate_starts(CL, n_c, BS, BO, LT0)
    t0 <- proc.time()[["elapsed"]]
    a <- plasmofit:::mat_exp_series(y0, CL, n_c, RR, MU, TS, 12)
    t1 <- proc.time()[["elapsed"]]
    b <- plasmofit:::conv_series(y0, CL, n_c, RR, MU, TS, M)
    t2 <- proc.time()[["elapsed"]]
    res <- rbind(res, data.frame(n_c = n_c, fft_M = M,
                                 max_rel_diff = max(abs(b - a) / abs(a)),
                                 mat_exp_s = t1 - t0, conv_s = t2 - t1))
    print(res[nrow(res), ], digits = 4, row.names = FALSE); flush.console()
}
cat("\n")
cat(sprintf("conv_series and mat_exp_series agree to %.1e at n_c <= 384 and drift\n",
            max(res$max_rel_diff[res$n_c <= 384])),
    sprintf("to %.1e by n_c = %d. That alone does NOT say which is wrong, so the\n",
            max(res$max_rel_diff), max(res$n_c)),
    "next section tests mat_exp_series against itself.\n", sep = "")
cat("\nFor scale: a log_prob relative difference of 1e-9 on a density of ~1e5 is\n",
    "1e-4 in absolute log-density units, far below anything that affects a fit.\n",
    "The question here is which component is losing digits, not whether the\n",
    "fits are wrong.\n", sep = "")

## ---- is mat_exp_series itself the one drifting? -------------------------
## conv_series and mat_exp_series disagree by ~1e-9 at high n_c, but that does
## not say which is wrong. mat_exp_series exponentiates a 2*n_c square matrix
## ONCE and then multiplies by it max(ts)/dt times, so its error has two
## sources that both grow with n_c: the matrix exponential's own accuracy at
## that dimension, and the accumulated matrix-vector products.
##
## The step count is the testable part. The exact answer cannot depend on dt,
## so if mat_exp_series moves when dt changes, the movement is its own error.
## ts = 72..216 is divisible by both 12 and 24, which halves the step count.
cat("\n=== is mat_exp_series self-consistent across step sizes? ===\n")
cat("Cells: max over the 7 times of |dt=12 - dt=24| / |dt=12| for",
    " mat_exp_series,\nand for conv_series the difference from each. The exact",
    " answer is dt-free,\nso any dt dependence is mat_exp_series its own",
    " accumulated error.\n\n", sep = "")
dtdep <- numeric(0); convdiff <- numeric(0)
for (n_c in c(192L, 384L, 768L)) {
    M <- plasmofit:::conv_min_M(n_c, max(TS), 35)
    y0 <- plasmofit:::generate_starts(CL, n_c, BS, BO, LT0)
    a12 <- plasmofit:::mat_exp_series(y0, CL, n_c, RR, MU, TS, 12)
    a24 <- plasmofit:::mat_exp_series(y0, CL, n_c, RR, MU, TS, 24)
    b   <- plasmofit:::conv_series(y0, CL, n_c, RR, MU, TS, M)
    cat(sprintf("  n_c = %4d : mat_exp dt12 vs dt24 = %.2e | conv vs dt12 = %.2e | conv vs dt24 = %.2e\n",
                n_c, max(abs(a12 - a24) / abs(a12)),
                max(abs(b - a12) / abs(a12)), max(abs(b - a24) / abs(a24))))
    dtdep <- c(dtdep, max(abs(a12 - a24) / abs(a12)))
    convdiff <- c(convdiff, max(abs(b - a12) / abs(a12)))
    flush.console()
}
cat("\n")
if (max(dtdep) < 1e-14 && max(convdiff) > 100 * max(dtdep)) {
    cat("VERDICT: mat_exp_series is stable -- halving its step count, which also\n",
        "changes the matrix it exponentiates, moves it by ~1e-16 and does NOT\n",
        "grow with n_c. So conv_series is what drifts, and the project reading of\n",
        "max_rel_diff (findings.md: 8e-13, 1.4e-12, 2.8e-12 at n_c 96/192/384) as\n",
        "the SERIES degrading rather than the reference stands.\n\n",
        "Why: an FFT convolution has an error floor of ~eps * max|u| in ABSOLUTE\n",
        "terms, while the observable at a trough is a small sum of small positive\n",
        "terms. As n_c rises the age distribution narrows and the troughs deepen,\n",
        "so that fixed absolute floor becomes a larger relative error exactly\n",
        "where the observable is smallest. mat_exp_series propagates the state\n",
        "vector directly, with no global transform, so small entries keep their\n",
        "relative accuracy.\n", sep = "")
} else {
    cat("mat_exp_series shows its own dt dependence, so it is not a clean arbiter\n",
        "here and neither map can be blamed without a higher-precision reference.\n",
        sep = "")
}
