## Is the OBSERVABLE period equal to the `cycle_length` parameter?
##
## Deterministic and noiseless: generate a trajectory at a known cycle_length
## and measure the interval between successive peaks of the detrended log10
## signal. If the observable period differs from the parameter, and the
## difference depends on n_c, then a fit must move cycle_length to compensate
## when n_c changes -- which is the real-data effect, with no estimation
## involved.
##
## Why it could differ: sequestration begins at a FIXED 18.58 h of a ~45 h
## cycle, so only the first ~40% of the cycle is visible and the observation
## window is ONE-SIDED. The stage distribution broadens at a rate n_c sets
## (sd after k cycles = cycle_length * sqrt(k / n_c)) and is right-skewed.
## Convolving a broadening, skewed age distribution with a one-sided window
## moves the centroid of the VISIBLE subpopulation, progressively.
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   Rscript --vanilla _scripts/nc-period-check.R

suppressPackageStartupMessages(library(plasmofit))

CL <- 45.012; R <- 6.5; MU <- 0; DT <- 0.25
ts <- seq(DT, 200, by = DT)

peaks <- function(nc, bs) {
    y0 <- plasmofit:::generate_starts(CL, nc, bs, 0.25, 1)
    y  <- plasmofit:::mat_exp_series(y0, CL, nc, R, MU, ts, DT)
    l  <- log10(y + 1)
    ## remove the exponential growth so the oscillation is what is measured
    r  <- residuals(lm(l ~ ts))
    i  <- which(diff(sign(diff(r))) == -2) + 1
    ## Peak times refined by fitting a parabola through the maximum and its
    ## two neighbours. Without this the peaks land on the DT grid and every
    ## interval is a multiple of DT, which is coarser than the differences
    ## being measured.
    refine <- function(k) {
        if (k < 2 || k > length(r) - 1) return(ts[k])
        y1 <- r[k - 1]; y2 <- r[k]; y3 <- r[k + 1]
        denom <- y1 - 2 * y2 + y3
        if (!is.finite(denom) || denom == 0) return(ts[k])
        ts[k] + DT * 0.5 * (y1 - y3) / denom
    }
    list(times = vapply(i, refine, numeric(1)), resid_sd = sd(r))
}

cat("=== observable period vs the cycle_length parameter ===\n")
cat("Cells: peak times and successive intervals of the detrended log10\n",
    "trajectory, in hours, generated noiselessly at cycle_length = ", CL,
    " h.\n`amplitude` is the sd of the detrended series.\n\n", sep = "")

for (bs in c(15, 400)) {
    cat("--- b_shape =", bs, "---\n")
    for (nc in c(96L, 192L, 384L)) {
        p <- peaks(nc, bs)
        iv <- diff(head(p$times, 5))
        cat(sprintf("n_c %3d | peaks %s\n", nc,
                    paste(sprintf("%.2f", head(p$times, 5)), collapse = ", ")))
        cat(sprintf("        | intervals %s | steady-state %+.2f h vs parameter | amplitude %.4f\n",
                    paste(sprintf("%.3f", iv), collapse = ", "),
                    iv[length(iv)] - CL, p$resid_sd))
    }
}

cat("\nA NEGATIVE steady-state difference means the observable period runs\n",
    "SHORT of the parameter, so a fit at that n_c must INFLATE cycle_length\n",
    "to reproduce a given observed period. If that deficit shrinks as n_c\n",
    "rises, the fitted cycle_length falls with n_c -- the real-data effect.\n",
    sep = "")
