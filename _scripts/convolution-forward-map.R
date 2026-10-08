## The chain, and its continuous relaxations, as convolutions in ABSOLUTE
## developmental age. Thread 15. Sourced by _scripts/ipm-prototype.R and
## _scripts/nc-profile-fast.R so there is ONE implementation; the derivation,
## the validation against mat_exp_series, and the caveats are in
## claude/ipm-decision.md.
##
## Working in absolute age rather than age modulo the cycle removes the wrap:
## transport is then a pure-birth process, growth is the weight R^divisions,
## and sequestration is a weight too, because circulating status is a function
## of position within the cycle and RESETS at division.
##
## Two things that are easy to get wrong and are both load-bearing:
##   - the chain applies the sequestration hazard with a ONE-STAGE LAG, so the
##     circulating fraction at stage k is G[k-1], not y[k];
##   - parasites in their FIRST cycle have not been reset, so their weight
##     depends on where they started. It is separable, hence w0 below.
## Using y[k], or skipping w0, costs ~12% and ~2.5% at the first observation.

## Per-mesh weights. All of these depend on cycle_length, so they must be
## rebuilt inside a likelihood, not cached across parameter values.
cfm_weights <- function(M, cl, bs, bo, lt0) {
    y <- plasmofit:::make_y_vals(cl, M)
    w <- plasmofit:::beta_starts(M, bs, bo, 10^lt0)
    G <- c(1, y / y[1])                 # G[k+1] = prod_{j<=k}(1-q_j), G[1] = 1
    list(y = y, w = w, G = G, w0 = w * y / G[seq_len(M)])
}

## Fold a kernel over non-negative cell offsets into the observable. Shared by
## every forward map here, so the ONLY difference between them is the kernel.
cfm_fold <- function(k, M, wts, R, mu, t) {
    u0 <- convolve(wts$w0, rev(k), type = "open")
    u  <- convolve(wts$w,  rev(k), type = "open")
    j <- seq_along(u); d <- (j - 1) %/% M; p <- (j - 1) %% M + 1
    g0 <- numeric(length(u)); g0[j <= M] <- wts$G[j[j <= M]]
    exp(-mu * t) * (sum(u0 * g0) + sum((u * R^d * wts$G[p])[d >= 1]))
}

## The chain, exactly. Equals mat_exp_series to ~1e-12; n_c must be an integer
## because it is simultaneously the lattice.
cfm_chain <- function(ts, cl, n_c, R, mu, bs, bo, lt0) {
    wts <- cfm_weights(n_c, cl, bs, bo, lt0)
    lam <- n_c / cl
    vapply(ts, function(t) {
        m <- lam * t; hi <- ceiling(m + 12 * sqrt(m) + 12)
        cfm_fold(dpois(0:hi, m), n_c, wts, R, mu, t)
    }, 0)
}

## The IPM with a gamma kernel: Gamma(shape = n_eff*t/cl, scale = cl/n_eff) has
## the chain's mean AND variance AND right skew, but n_eff is CONTINUOUS and
## independent of the mesh M. This is the forward map an IPM would use.
cfm_gamma <- function(ts, cl, M, n_eff, R, mu, bs, bo, lt0) {
    wts <- cfm_weights(M, cl, bs, bo, lt0)
    h <- cl / M
    vapply(ts, function(t) {
        hi <- ceiling((t + 12 * sqrt(cl * t / n_eff)) / h) + 12
        e <- (seq.int(0, hi + 1) - 0.5) * h; e[1] <- 0
        k <- diff(pgamma(e, shape = n_eff * t / cl, scale = cl / n_eff))
        cfm_fold(k / sum(k), M, wts, R, mu, t)
    }, 0)
}

## The IPM with a gaussian kernel: same mean and variance, no skew. Kept
## because the difference between this and cfm_gamma is how much the skew is
## worth, which turns out to be under 0.003 log10 units.
cfm_gauss <- function(ts, cl, M, sigma_c, R, mu, bs, bo, lt0) {
    wts <- cfm_weights(M, cl, bs, bo, lt0)
    h <- cl / M
    vapply(ts, function(t) {
        sd_a <- sigma_c * sqrt(t / cl)
        hi <- ceiling((t + 12 * sd_a) / h) + 12
        k <- dnorm(seq.int(0, hi) * h, mean = t, sd = sd_a) * h
        cfm_fold(k / sum(k), M, wts, R, mu, t)
    }, 0)
}
