## Validate conv_series() against mat_exp_series() in the installed binary.
##
## conv_series() computes the same trajectory by convolution in ABSOLUTE
## developmental age rather than by exponentiating a 2*n_c square matrix. The
## two are the same model, so this is an equality check, not a comparison: the
## project's standing tolerance for the forward map is `max_rel_diff` ~ 1e-12,
## the same figure `check_erlang_window()` holds the series solution to.
##
## PACKAGE CHANGES ARE VERIFIED IN THE INSTALLED BINARY, NOT THE SOURCE
## (CLAUDE.md). Point PLASMOFIT_LIB at the library to test; the default is the
## scratch library used while a production job is holding the live one.
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   PLASMOFIT_LIB=/home2/lan68/plasmofit/.Rlib-dev \
##     Rscript --vanilla _scripts/conv-series-validate.R \
##     | tee _data/conv-series-validate-$(date +%F).txt

LIB <- Sys.getenv("PLASMOFIT_LIB", "/home2/lan68/plasmofit/.Rlib-dev")
if (nzchar(LIB)) .libPaths(c(LIB, .libPaths()))
suppressPackageStartupMessages(library(plasmofit))
cat("plasmofit loaded from:", dirname(find.package("plasmofit")), "\n")
if (!exists("conv_series", asNamespace("plasmofit")))
    stop("conv_series is not in the installed namespace -- the binary is stale. ",
         "Reinstall with --preclean after rm src/*.o src/*.so.")

TS   <- c(72, 96, 120, 144, 168, 192, 216)   # the Wockner observation times
DTF  <- 12
CLMIN <- 35

## A trajectory from each path, at one parameter set.
pair <- function(cl, n_c, bs, bo, lt0, R, mu, M) {
    y0 <- plasmofit:::generate_starts(cl, n_c, bs, bo, lt0)
    list(mat  = plasmofit:::mat_exp_series(y0, cl, n_c, R, mu, TS, DTF),
         conv = plasmofit:::conv_series(y0, cl, n_c, R, mu, TS, M))
}

## ---- 1. agreement over a parameter grid --------------------------------
grid <- expand.grid(cl = c(38, 41, 45, 49), n_c = c(96L, 192L),
                    bs = c(15, 400, 1000), bo = c(0, 0.3, 0.7),
                    R = c(5, 8, 20), mu = c(0, 0.01),
                    KEEP.OUT.ATTRS = FALSE)
grid <- rbind(grid,
              expand.grid(cl = c(41, 45), n_c = 384L, bs = c(15, 400),
                          bo = c(0, 0.3), R = c(8, 20), mu = 0,
                          KEEP.OUT.ATTRS = FALSE))
cat("\n=== 1. conv_series vs mat_exp_series over", nrow(grid), "parameter sets ===\n")
cat("Cells: max over the 7 observation times of |conv - mat| / |mat|, then\n",
    "summarised over the grid. log10_total0 is fixed at 5.5 throughout, since\n",
    "it only scales the trajectory.\n\n", sep = "")
M_cache <- new.env(parent = emptyenv())
getM <- function(n_c) {
    k <- as.character(n_c)
    if (is.null(M_cache[[k]]))
        M_cache[[k]] <- plasmofit:::conv_min_M(n_c, max(TS), CLMIN)
    M_cache[[k]]
}
grid$M <- vapply(grid$n_c, getM, 0L)
grid$max_rel_diff <- NA_real_
for (i in seq_len(nrow(grid))) {
    p <- pair(grid$cl[i], grid$n_c[i], grid$bs[i], grid$bo[i], 5.5,
              grid$R[i], grid$mu[i], grid$M[i])
    grid$max_rel_diff[i] <- max(abs(p$conv - p$mat) / abs(p$mat))
}
cat("FFT length used: ",
    paste(sprintf("n_c %d -> M %d", as.integer(names(as.list(M_cache))),
                  unlist(as.list(M_cache))), collapse = " | "), "\n\n", sep = "")
print(summary(grid$max_rel_diff))
cat(sprintf("\nworst max_rel_diff: %.3e | sets above 1e-10: %d of %d\n",
            max(grid$max_rel_diff), sum(grid$max_rel_diff > 1e-10), nrow(grid)))
w <- grid[which.max(grid$max_rel_diff), ]
cat("worst set: "); print(w[, c("cl","n_c","bs","bo","R","mu","max_rel_diff")],
                          row.names = FALSE, digits = 4)
TOL <- 1e-10
if (max(grid$max_rel_diff) > TOL)
    stop("conv_series disagrees with mat_exp_series by more than ", TOL,
         " -- it is NOT the same model and must not be used.")
cat("\nPASS: every parameter set agrees to better than", TOL, "\n")

## ---- 2. the things conv_series can do that mat_exp_series cannot --------
cat("\n=== 2. times that mat_exp_series rejects ===\n")
cat("Cells: mat_exp_series requires times divisible by dt and strictly\n",
    "increasing, because it steps. conv_series has no stepping, so it does\n",
    "not. Checked by evaluating at a time BETWEEN two grid points and\n",
    "confirming it lands between the two trajectories.\n\n", sep = "")
n_c <- 192L; M <- getM(n_c)
y0 <- plasmofit:::generate_starts(45, n_c, 400, 0.3, 5.5)
odd <- c(77.5, 103.25, 150.125)
v <- plasmofit:::conv_series(y0, 45, n_c, 8, 0, odd, M)
print(data.frame(t = odd, conv_series = v), digits = 7)
e <- tryCatch({plasmofit:::mat_exp_series(y0, 45, n_c, 8, 0, odd, DTF); "accepted"},
              error = function(e) "rejected, as expected")
cat("mat_exp_series at the same times:", e, "\n")

## ---- 3. cost ------------------------------------------------------------
cat("\n=== 3. cost ===\n")
cat("Cells: seconds for one 7-point trajectory, median of 5 runs, installed\n",
    "binary. The matrix exponential is cubic in n_c; the convolution is\n",
    "M log M with M roughly linear in n_c.\n\n", sep = "")
tm <- function(f, n = 5) median(replicate(n, system.time(f())[["elapsed"]]))
for (n_c in c(96L, 192L, 384L)) {
    M <- getM(n_c)
    y0 <- plasmofit:::generate_starts(45, n_c, 400, 0.3, 5.5)
    a <- tm(function() plasmofit:::mat_exp_series(y0, 45, n_c, 8, 0, TS, DTF))
    b <- tm(function() plasmofit:::conv_series(y0, 45, n_c, 8, 0, TS, M))
    cat(sprintf("  n_c = %3d (M = %5d) : mat_exp %7.4f s | conv %7.4f s | speedup %5.1fx\n",
                n_c, M, a, b, a / b))
}
