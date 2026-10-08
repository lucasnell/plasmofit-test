## Which law does synchrony decay by? A screening test before any rewrite.
##
## Thread 15. The Erlang chain makes the stage distribution spread as
## sd = cycle_length * sqrt(k / n_c) after k cycles -- DIFFUSIVE, sd ~ sqrt(k).
## A 1-D IPM with a Gaussian development kernel gives the same sqrt law, so
## switching to one would decouple the rate from the mesh without testing the
## FORM. The competing form is between-parasite heterogeneity in cycle
## duration: a parasite 1% fast is 1% further ahead every cycle, so the spread
## grows LINEARLY in k. These are different models of the same phenomenon and
## the choice should be made on evidence, not on machinery.
##
## This is a maximum-likelihood screen, not a replacement for the Stan fit.
## It exploits a structural fact: b_shape, b_offset and log10_total0 are per
## grp_init (trial x inoc_size) and R and cycle_length are per trial, so
## within a grp_init unit EVERY series shares EVERY parameter and the unit is
## a single trajectory with its observations as replicates. 14 units, ~6
## parameters each.
##
## CAVEATS, which bound what this can claim:
##  - No hierarchy. Each unit is fitted independently, so SJ733IBSMCS's two
##    cohorts get their own cycle_length and R although the real model shares
##    them. This is a screen for a LAW, not an estimate of anything.
##  - Maximum likelihood with sd profiled out, so no priors and no elpd. A
##    clear winner here justifies the Stan work; a tie says do not bother.
##  - Model B still contains the chain, so it NESTS model A at sigma = 0.
##    Run at a high n_c so the chain's own sqrt contribution is small and
##    sigma carries the spread.
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   Rscript --vanilla _scripts/decay-law-test.R

suppressPackageStartupMessages({
    library(readr); library(dplyr); library(tibble)
})

## Gauss-Hermite nodes and weights by Golub-Welsch, so no extra package is
## needed: the Jacobi matrix for the physicists' Hermite weight exp(-x^2) is
## symmetric tridiagonal with zero diagonal and off-diagonal sqrt(k/2).
## Returns nodes and weights for integral f(x) exp(-x^2) dx.
gauss_hermite <- function(n) {
    k <- seq_len(n - 1)
    J <- diag(0, n)
    J[cbind(k, k + 1)] <- sqrt(k / 2)
    J[cbind(k + 1, k)] <- sqrt(k / 2)
    e <- eigen(J, symmetric = TRUE)
    o <- order(e$values)
    list(nodes = e$values[o], weights = sqrt(pi) * (e$vectors[1, o])^2)
}

## Cost is roughly CUBIC in n_c -- one trajectory takes 0.0148 s at 96,
## 0.153 at 192, 1.18 at 384, 9.21 at 768 (measured) -- because the matrix
## exponential is over a 2*n_c square. Model B multiplies that by QNODE, so
## the grid is chosen to keep one unit inside a few hours and the work is
## split one SLURM task per grp_init unit.
##
## MODE, added 2026-10-08. The default reruns the whole screen. Mode "a768"
## runs model A alone at n_c = 768 and writes to a separate file, to extend
## the n_c profile by one rung without touching the finished results. It is
## model A only because the question it answers is whether the Erlang family
## is still gaining as the mesh refines -- model B is irrelevant to that and
## costs 7x.
MODE <- Sys.getenv("DECAY_MODE", "full")
if (!MODE %in% c("full", "a768", "profile")) stop("DECAY_MODE must be full, a768 or profile")

if (MODE == "full") {
    NC_A <- c(96L, 192L, 384L)  # model A: the sqrt law, slope set by n_c
    NC_B <- c(192L, 384L)       # model B: sigma adds a linear component
    OUT  <- "_data/decay-law-unit%02d.rds"
} else if (MODE == "a768") {
    NC_A <- 768L
    NC_B <- integer(0)
    OUT  <- "_data/decay-law-a768-unit%02d.rds"
} else if (MODE == "profile") {
    ## Intermediate rungs, so the n_c profile has enough points to show
    ## CURVATURE rather than just an ordering. n_c maps exactly to dispersion
    ## -- transit CV is 1/sqrt(n_c) -- so profiling over n_c with b_shape free
    ## at every rung IS a profile likelihood in the dispersion magnitude with
    ## the initial-spread nuisance concentrated out. Staying inside the sqrt
    ## family costs nothing, because 29684 showed the form is not learnable.
    NC_A <- c(128L, 256L, 512L)
    NC_B <- integer(0)
    OUT  <- "_data/decay-law-prof-unit%02d.rds"
}
QNODE  <- 7L
MU     <- 0
DTF    <- 12                   # gcd of the observation times

d <- read_csv("_data/wockner-cleaned.csv", col_types = "cccdcdd") |>
    mutate(unit = paste(trial, inoc_size, sep = "|"))

gh <- gauss_hermite(QNODE)
## Check the rule before trusting it: E[exp(X)] for X ~ N(0, s^2) is
## exp(s^2/2), which a correct rule reproduces to many digits.
local({
    s_ <- 0.3
    got <- sum(gh$weights / sqrt(pi) * exp(sqrt(2) * s_ * gh$nodes))
    want <- exp(s_^2 / 2)
    cat(sprintf("Gauss-Hermite check: %.10f vs %.10f (rel err %.2e)\n",
                got, want, abs(got - want) / want))
    stopifnot(abs(got - want) / want < 1e-8)
})

## trajectory for one unit at a set of times, total0 = 10^lt0
traj1 <- function(ts, cl, n_c, bs, bo, lt0, R) {
    y0 <- plasmofit:::generate_starts(cl, n_c, bs, bo, lt0)
    plasmofit:::mat_exp_series(y0, cl, n_c, R, MU, ts, DTF)
}

## mixture over cycle_length: lognormal with log-sd sigma, Gauss-Hermite
traj_mix <- function(ts, cl, n_c, bs, bo, lt0, R, sigma) {
    if (sigma <= 1e-8) return(traj1(ts, cl, n_c, bs, bo, lt0, R))
    nodes <- cl * exp(sqrt(2) * sigma * gh$nodes)
    w     <- gh$weights / sqrt(pi)
    out   <- numeric(length(ts))
    for (q in seq_len(QNODE)) {
        if (nodes[q] <= 30 || nodes[q] >= 60) return(rep(NA_real_, length(ts)))
        out <- out + w[q] * traj1(ts, nodes[q], n_c, bs, bo, 0, R)
    }
    out * 10^lt0
}

## mat_exp_series requires STRICTLY INCREASING times, and a unit's rows are
## many series sharing a handful of sampling times (10 distinct for
## Mefloquine's 132 observations). So the trajectory is evaluated once at the
## sorted unique times and indexed back out to the observations.

## profile log-likelihood: normal on log10(y+1), sd profiled out at its MLE
nll <- function(par, ut, idx, yobs, n_c, with_sigma) {
    cl  <- 35 + 15 / (1 + exp(-par[1]))
    bs  <- exp(par[2]); bo <- par[3] %% 1
    lt0 <- par[4];      R  <- 1 + exp(par[5])
    sg  <- if (with_sigma) exp(par[6]) else 0
    if (!is.finite(bs) || bs > 5000 || !is.finite(R) || R > 200) return(1e10)
    yu <- tryCatch(traj_mix(ut, cl, n_c, bs, bo, lt0, R, sg),
                   error = function(e) NULL)
    if (is.null(yu) || any(!is.finite(yu)) || any(yu < 0)) return(1e10)
    yh <- yu[idx]
    r <- log10(yh + 1) - log10(yobs + 1)
    n <- length(r); s2 <- sum(r^2) / n
    if (!is.finite(s2) || s2 <= 0) return(1e10)
    0.5 * n * (log(2 * pi * s2) + 1)
}

fit_unit <- function(ut, idx, yobs, n_c, with_sigma) {
    ## Two starts, not more: cost is cubic in n_c and model B multiplies it
    ## by QNODE. Starts bracket the b_shape range this project has seen
    ## (posteriors of 15 and of 65+ depending on the prior).
    starts <- list(c(0, log(15), 0.25, 1, log(5.5), log(0.02)),
                   c(0.5, log(60), 0.6, 0.5, log(8), log(0.05)))
    best <- NULL
    for (s0 in starts) {
        p0 <- if (with_sigma) s0 else s0[1:5]
        o <- tryCatch(optim(p0, nll, ut = ut, idx = idx, yobs = yobs, n_c = n_c,
                            with_sigma = with_sigma, method = "Nelder-Mead",
                            control = list(maxit = 2000, reltol = 1e-9)),
                      error = function(e) NULL)
        if (!is.null(o) && (is.null(best) || o$value < best$value)) best <- o
    }
    best
}

units <- sort(unique(d$unit))
task <- suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", NA)))
if (is.na(task)) stop("set SLURM_ARRAY_TASK_ID to 1-", length(units))
if (task < 1 || task > length(units))
    stop("SLURM_ARRAY_TASK_ID must be 1-", length(units), ", got ", task)

u  <- units[task]
du <- d |> filter(unit == u)
yo  <- du$para
ut  <- sort(unique(du$time))
idx <- match(du$time, ut)
stopifnot(!is.unsorted(ut, strictly = TRUE), !anyNA(idx))
cat("unit", task, "of", length(units), ":", u,
    "|", nrow(du), "observations,", n_distinct(du$id), "series,",
    length(ut), "distinct times\n")

res <- list()
for (nc in NC_A) {
    t0 <- proc.time()[["elapsed"]]
    o <- fit_unit(ut, idx, yo, nc, FALSE)
    res[[length(res) + 1]] <- tibble(unit = u, model = "A_sqrt",
        n_c = nc, sigma = 0, npar = 6, ll = -o$value, n = length(yo),
        cl_hat = 35 + 15 / (1 + exp(-o$par[1])), bs_hat = exp(o$par[2]))
    cat(sprintf("  A n_c=%3d ll %+.3f (%.0f s)\n", nc, -o$value,
                proc.time()[["elapsed"]] - t0)); flush.console()
}
for (nc in NC_B) {
    t0 <- proc.time()[["elapsed"]]
    o <- fit_unit(ut, idx, yo, nc, TRUE)
    res[[length(res) + 1]] <- tibble(unit = u, model = "B_linear",
        n_c = nc, sigma = exp(o$par[6]), npar = 6, ll = -o$value,
        n = length(yo),
        cl_hat = 35 + 15 / (1 + exp(-o$par[1])), bs_hat = exp(o$par[2]))
    cat(sprintf("  B n_c=%3d ll %+.3f sigma %.4f (%.0f s)\n", nc, -o$value,
                exp(o$par[6]), proc.time()[["elapsed"]] - t0)); flush.console()
}
r <- bind_rows(res)
saveRDS(r, sprintf(OUT, task))
cat("\nwrote ", sprintf(OUT, task), "\n", sep = "")
print(as.data.frame(r), digits = 5)

## Aggregation lives in _scripts/decay-law-read.R, because each task here
## handles one unit and the comparison is across all of them.
