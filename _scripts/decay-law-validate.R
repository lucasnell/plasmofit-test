## Does the cycle_length mixture actually implement a LINEAR decay law?
##
## The claim behind model B in decay-law-test.R: if cycle duration varies
## between parasites with coefficient of variation sigma, a parasite that is
## 1% fast is 1% further ahead every cycle, so the phase spread after k
## cycles is sigma * k CYCLES -- linear in k, against the chain's
## sqrt(k / n_c). This checks the implementation does that, rather than
## assuming it.
##
## Method: generate noiselessly at a high n_c, so the chain's own sqrt
## contribution is small, and measure the oscillation amplitude by cycle.
## More spread means more damping, so amplitude is the observable proxy.
##
##   Rscript --vanilla _scripts/decay-law-validate.R

suppressPackageStartupMessages(library(plasmofit))

gauss_hermite <- function(n) {
    k <- seq_len(n - 1); J <- diag(0, n)
    J[cbind(k, k + 1)] <- sqrt(k / 2); J[cbind(k + 1, k)] <- sqrt(k / 2)
    e <- eigen(J, symmetric = TRUE); o <- order(e$values)
    list(nodes = e$values[o], weights = sqrt(pi) * (e$vectors[1, o])^2)
}
gh <- gauss_hermite(7L)

CL <- 45.012; NC <- 192L; R <- 6.5; DT <- 0.5
ts <- seq(DT, 225, by = DT)

mix <- function(sigma) {
    if (sigma <= 0)
        return(plasmofit:::mat_exp_series(
            plasmofit:::generate_starts(CL, NC, 15, 0.25, 1),
            CL, NC, R, 0, ts, DT))
    nodes <- CL * exp(sqrt(2) * sigma * gh$nodes)
    w <- gh$weights / sqrt(pi)
    out <- numeric(length(ts))
    for (q in seq_along(nodes))
        out <- out + w[q] * plasmofit:::mat_exp_series(
            plasmofit:::generate_starts(nodes[q], NC, 15, 0.25, 0),
            nodes[q], NC, R, 0, ts, DT)
    out * 10
}

## amplitude per cycle: sd of the detrended log10 signal within each cycle
amp_by_cycle <- function(y) {
    l <- log10(y + 1); r <- residuals(lm(l ~ ts))
    cyc <- floor(ts / CL)
    tapply(r, cyc, function(z) diff(range(z)) / 2)
}

cat("=== amplitude by cycle, n_c =", NC, "===\n")
cat("Cells: half-range of the detrended log10 trajectory within each cycle,\n",
    "a damping proxy. More between-parasite spread means faster damping.\n\n",
    sep = "")
res <- sapply(c(0, 0.02, 0.05), function(s) amp_by_cycle(mix(s)))
colnames(res) <- paste0("sigma=", c(0, 0.02, 0.05))
print(round(res, 4))

cat("\n=== implied extra spread, relative to sigma = 0 ===\n")
cat("Cells: log of the amplitude ratio to the sigma = 0 column. For a\n",
    "Gaussian phase spread s, amplitude is damped by exp(-2*pi^2*s^2), so\n",
    "-log(ratio) is proportional to s^2. A LINEAR law (s = sigma*k) makes\n",
    "that quantity grow as k^2; a sqrt law makes it grow as k.\n\n", sep = "")
for (j in 2:3) {
    lr <- -log(res[, j] / res[, 1])
    k  <- as.numeric(rownames(res)) + 0.5
    ok <- is.finite(lr) & lr > 0
    cat(sprintf("%s: -log ratio by cycle %s\n", colnames(res)[j],
                paste(sprintf("%.4f", lr), collapse = ", ")))
    if (sum(ok) >= 3)
        cat(sprintf("   slope of log(-log ratio) on log(k) = %.2f (2 = linear law, 1 = sqrt law)\n",
                    coef(lm(log(lr[ok]) ~ log(k[ok])))[2]))
}
