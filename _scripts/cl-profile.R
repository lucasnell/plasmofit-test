## Profile likelihood of cycle_length for one grp_init unit, at n_c = 384,
## prior-free: fix cycle_length on a grid, maximise over everything else.
## Asks whether the two regions ~3 h apart that MCMC chains found for the two
## MMV048 trials (claude/findings.md, "Seed probe, first read") are two peaks
## of the likelihood itself, and how their heights compare. A profile does not
## depend on chains moving between regions, which is exactly what failed.
##
## Same forward map (the exact chain as a convolution, equal to the Stan
## model's to ~1e-12; _scripts/convolution-forward-map.R) and the same
## likelihood as _scripts/nc-profile-fast.R: normal on log10(y + 1) with ONE
## residual variance per unit, profiled out. This differs from the Stan model
## (per-cohort sd, priors, hierarchy), so heights compare regions WITHIN a
## unit; they are not posterior masses.
##
## Design, fixed before any profile existed (2026-10-14):
##   units   MMV048_PIB|1800, MMV048_PartB|2800, and KAE609|1800 as a control
##           (the trial whose 20 pin-15 chains agreed best)
##   b_shape fixed at 15 (as entry 43) or free
##   grid    cycle_length 36.0-49.0 h by 0.5 h
##   optimiser at each grid point 8 cold starts (nc-profile-fast.R's start set,
##           minus the cycle_length coordinate), then warm-start sweeps from
##           neighbouring grid points, up and down, repeated until no point
##           improves by more than 1e-3 (at most 4 sweeps). More starts can
##           only raise a maximum, so the profile is the best found.
##   reading TWO PEAKS if two local maxima are separated by a dip of at least
##           2 log-likelihood units below the lower of them; otherwise ONE.
##
## Environment: CLP_TASK = 1-6 (unit x mode, see PLAN below).
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   sbatch _scripts/cl-profile.sh

suppressPackageStartupMessages({library(readr); library(dplyr)})
.libPaths(c("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6", .libPaths()))
source("_scripts/convolution-forward-map.R")

N_C <- 384L; MU <- 0
GRID <- seq(36, 49, by = 0.5)
PLAN <- expand.grid(mode = c("bs15", "free"),
                    unit = c("MMV048_PIB|1800", "MMV048_PartB|2800", "KAE609|1800"),
                    stringsAsFactors = FALSE)
task <- as.integer(Sys.getenv("CLP_TASK", Sys.getenv("SLURM_ARRAY_TASK_ID", NA)))
stopifnot(!is.na(task), task >= 1, task <= nrow(PLAN))
u <- PLAN$unit[task]; mode <- PLAN$mode[task]
OUT <- sprintf("_data/cl-profile-%s-%s.rds", gsub("[|]", "_", u), mode)
if (file.exists(OUT)) stop(OUT, " exists; refusing to overwrite")

d <- read_csv("_data/wockner-cleaned.csv", col_types = "cccdcdd") |>
    mutate(unit = paste(trial, inoc_size, sep = "|"))
du <- filter(d, unit == u); stopifnot(nrow(du) > 0)
yo <- du$para; ut <- sort(unique(du$time)); idx <- match(du$time, ut)
cat(sprintf("task %d: %s, b_shape %s, %d observations, %d times, grid %d points\n",
            task, u, mode, length(yo), length(ut), length(GRID))); flush.console()

## Free parameters: [log bs], bo, lt0, log(R - 1). Same transforms and
## guards as nc-profile-fast.R's nll.
nll <- function(par, cl) {
    if (mode == "bs15") { bs <- 15; p <- par } else { bs <- exp(par[1]); p <- par[-1] }
    bo <- p[1] %% 1; lt0 <- p[2]; R <- 1 + exp(p[3])
    if (!is.finite(bs) || bs > 5000 || bs <= 2 || !is.finite(R) || R > 200) return(1e10)
    yu <- tryCatch(cfm_chain(ut, cl, N_C, R, MU, bs, bo, lt0), error = function(e) NULL)
    if (is.null(yu) || any(!is.finite(yu)) || any(yu < 0)) return(1e10)
    r <- log10(yu[idx] + 1) - log10(yo + 1)
    n <- length(r); s2 <- sum(r^2) / n
    if (!is.finite(s2) || s2 <= 0) return(1e10)
    0.5 * n * (log(2 * pi * s2) + 1)
}
STARTS_ALL <- list(                     # nc-profile-fast.R, cl coordinate dropped
    c(log(15),   0.25, 1.0, log(5.5)), c(log(60),   0.60, 0.5, log(8)),
    c(log(8),    0.05, 2.0, log(3)),   c(log(300),  0.85, 0.0, log(15)),
    c(log(1500), 0.40, 1.5, log(6)),   c(log(40),   0.70, 0.8, log(10)),
    c(log(5),    0.15, 1.2, log(20)),  c(log(900),  0.95, 0.3, log(4)))
STARTS <- if (mode == "bs15") lapply(STARTS_ALL, `[`, -1) else STARTS_ALL
opt <- function(s0, cl) tryCatch(optim(s0, nll, cl = cl, method = "Nelder-Mead",
                                       control = list(maxit = 3000, reltol = 1e-10)),
                                 error = function(e) NULL)

best <- vector("list", length(GRID))
t0 <- proc.time()[["elapsed"]]
for (i in seq_along(GRID)) {
    for (s0 in STARTS) {
        o <- opt(s0, GRID[i])
        if (!is.null(o) && (is.null(best[[i]]) || o$value < best[[i]]$value)) best[[i]] <- o
    }
    cat(sprintf("  cold  cl %.1f  ll %+.4f  (%.0f s)\n", GRID[i], -best[[i]]$value,
                proc.time()[["elapsed"]] - t0)); flush.console()
}
cold_ll <- -vapply(best, `[[`, 0, "value")

for (sweep in 1:4) {
    gain <- 0
    for (ord in list(seq_along(GRID)[-1], rev(seq_along(GRID))[-1])) {
        for (i in ord) {
            j <- if (identical(ord, seq_along(GRID)[-1])) i - 1 else i + 1
            o <- opt(best[[j]]$par, GRID[i])
            if (!is.null(o) && o$value < best[[i]]$value - 1e-9) {
                gain <- max(gain, best[[i]]$value - o$value); best[[i]] <- o
            }
        }
    }
    cat(sprintf("  sweep %d: largest gain %.4f  (%.0f s)\n", sweep, gain,
                proc.time()[["elapsed"]] - t0)); flush.console()
    if (gain < 1e-3) break
}

res <- data.frame(unit = u, mode = mode, cl = GRID,
                  ll = -vapply(best, `[[`, 0, "value"), cold_ll = cold_ll,
                  bs = if (mode == "bs15") 15 else exp(vapply(best, function(o) o$par[1], 0)),
                  bo = vapply(best, function(o) tail(o$par, 3)[1] %% 1, 0),
                  lt0 = vapply(best, function(o) tail(o$par, 3)[2], 0),
                  R = 1 + exp(vapply(best, function(o) tail(o$par, 3)[3], 0)),
                  n = length(yo))
saveRDS(res, OUT)
cat("wrote", OUT, "\n")
