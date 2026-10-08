## Is the roughness in the nc-profile-fast profiles optimiser noise, and do
## more starts remove it? Thread 15.
##
## A profile likelihood in a smooth parameter should be smooth. Both profiles
## from SLURM 30576 deviate from a loess fit in log(knot) by up to 2.39
## (chain) and 7.17 (gamma IPM) log-likelihood units, against the 2-unit
## currency the pre-registered reading rule uses. Under-optimisation can only
## push ll DOWN, so the rungs with large NEGATIVE residuals are the suspects.
##
## This refits a few of them with more starts and reports what that buys. If
## the gain is near zero the roughness is real structure; if it is several
## units the harness is under-optimised and the whole profile needs redoing.
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   Rscript --vanilla _scripts/profile-noise-check.R | tee _data/profile-noise-check-$(date +%F).txt

suppressPackageStartupMessages({library(readr); library(dplyr)})
source("_scripts/convolution-forward-map.R")

UNITS  <- c(1L, 2L)                       # two units, enough to see a pattern
NE     <- c(930.37, 2264.01, 4096.00)     # rungs with large negative residuals
M      <- 192L
MU     <- 0

d <- read_csv("_data/wockner-cleaned.csv", col_types = "cccdcdd") |>
    mutate(unit = paste(trial, inoc_size, sep = "|"))
units <- sort(unique(d$unit))

## The two starts the production script uses, then six more spread over the
## same box. Not random: a fixed list, so this is reproducible.
STARTS2 <- list(c(0, log(15), 0.25, 1, log(5.5)),
                c(0.5, log(60), 0.6, 0.5, log(8)))
STARTS8 <- c(STARTS2, list(
    c(-0.5, log(8),   0.05, 2.0, log(3)),
    c( 0.9, log(300), 0.85, 0.0, log(15)),
    c( 0.2, log(1500),0.40, 1.5, log(6)),
    c(-0.2, log(40),  0.70, 0.8, log(10)),
    c( 0.7, log(5),   0.15, 1.2, log(20)),
    c( 0.0, log(900), 0.95, 0.3, log(4))))

cat("Cells: ll is the maximised profile log-likelihood for ONE unit at one\n",
    "n_eff, gamma kernel, mesh 192. ll2 uses the two starts the production\n",
    "script uses; ll8 adds six more. gain = ll8 - ll2, and can only be >= 0.\n",
    "Any gain near or above the 2-unit currency means the production profile\n",
    "is under-optimised at that rung.\n\n", sep = "")

out <- list()
for (ti in UNITS) {
    u <- units[ti]; du <- filter(d, unit == u)
    yo <- du$para; ut <- sort(unique(du$time)); idx <- match(du$time, ut)
    nll <- function(par, ne) {
        cl <- 35 + 15/(1 + exp(-par[1])); bs <- exp(par[2]); bo <- par[3] %% 1
        lt0 <- par[4]; R <- 1 + exp(par[5])
        if (!is.finite(bs) || bs > 5000 || bs <= 2 || !is.finite(R) || R > 200)
            return(1e10)
        yu <- tryCatch(cfm_gamma(ut, cl, M, ne, R, MU, bs, bo, lt0),
                       error = function(e) NULL)
        if (is.null(yu) || any(!is.finite(yu)) || any(yu < 0)) return(1e10)
        r <- log10(yu[idx] + 1) - log10(yo + 1)
        n <- length(r); s2 <- sum(r^2)/n
        if (!is.finite(s2) || s2 <= 0) return(1e10)
        0.5 * n * (log(2*pi*s2) + 1)
    }
    best <- function(st, ne) {
        b <- Inf
        for (s0 in st) {
            o <- tryCatch(optim(s0, nll, ne = ne, method = "Nelder-Mead",
                                control = list(maxit = 2000, reltol = 1e-9)),
                          error = function(e) NULL)
            if (!is.null(o) && o$value < b) b <- o$value
        }
        -b
    }
    for (ne in NE) {
        t0 <- proc.time()[["elapsed"]]
        a <- best(STARTS2, ne); b <- best(STARTS8, ne)
        out[[length(out)+1]] <- tibble(unit = u, n_eff = ne, ll2 = a, ll8 = b,
                                       gain = b - a,
                                       secs = proc.time()[["elapsed"]] - t0)
        cat(sprintf("  %-20s n_eff %7.1f  ll2 %+9.4f  ll8 %+9.4f  gain %+7.4f  (%.0f s)\n",
                    u, ne, a, b, b - a, proc.time()[["elapsed"]] - t0))
        flush.console()
    }
}
o <- bind_rows(out)
cat(sprintf("\nmax gain %.3f | mean gain %.3f | rungs with gain > 1 ll unit: %d of %d\n",
            max(o$gain), mean(o$gain), sum(o$gain > 1), nrow(o)))
if (max(o$gain) > 1) {
    cat("UNDER-OPTIMISED. The production profile understates ll at these rungs,\n",
        "so its shape is partly an artefact of where the optimiser stopped.\n",
        "Re-run nc-profile-fast.R with the larger start list before reading any\n",
        "verdict from it.\n", sep = "")
} else {
    cat("NOT under-optimisation. The roughness is structure, or noise of a\n",
        "different kind; more starts will not fix it.\n", sep = "")
}
saveRDS(o, "_data/profile-noise-check.rds")
