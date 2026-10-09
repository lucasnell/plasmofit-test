## Does conv_series() differentiate correctly inside Stan?
##
## The numerical agreement in _scripts/conv-series-validate.R is on the EXPOSED
## functions, which are evaluated in double precision with no autodiff. A
## production fit needs gradients, and those go through a complex FFT whose
## adjoint is newer and less exercised than matrix_exp. So this compiles both
## paths behind a data switch and compares, at identical parameter values:
##
##   1. log_prob      -- the two should agree to ~1e-10, which also confirms the
##                       functions behave the same inside Stan's own evaluation
##                       rather than only in the R bindings;
##   2. grad_log_prob -- conv_series against mat_exp_series, which is the real
##                       test of the FFT adjoint, and against central finite
##                       differences, which catches an error common to both.
##
## Nothing here touches the installed package: stan_model() compiles to a
## temporary directory, and the plasmofit namespace is only read.
##
##
## NOTE, 2026-10-09: conv_series() has been REMOVED from the package (revert
## commit bba580f) because it is slower and less accurate than the
## Erlang-window series at every n_c this project uses. This script therefore
## needs a plasmofit built from commit 05c9c1f to run. It is kept because the
## numbers in claude/findings.md come from it and a claim should name the
## script that produced it; the saved output beside it is the record.
##   cd /home2/lan68/plasmofit/plasmofit-test
##   Rscript --vanilla _scripts/conv-series-grad.R \
##     | tee _data/conv-series-grad-$(date +%F).txt

LIB <- Sys.getenv("PLASMOFIT_LIB", "/home2/lan68/plasmofit/.Rlib-dev")
if (nzchar(LIB)) .libPaths(c(LIB, .libPaths()))
suppressPackageStartupMessages({library(rstan); library(plasmofit)})
cat("plasmofit from:", dirname(find.package("plasmofit")), "\n")
PKG <- "/home2/lan68/plasmofit/plasmofit/inst/stan"
rstan_options(auto_write = FALSE)

TS <- c(72, 96, 120, 144, 168, 192, 216)
N_C <- 96L                      # small, because mat_exp_series is cubic
M <- plasmofit:::conv_min_M(N_C, max(TS), 35)
TRUTH <- list(cl = 44, bs = 60, bo = 0.3, lt0 = 5.5, R = 8)
y0 <- plasmofit:::generate_starts(TRUTH$cl, N_C, TRUTH$bs, TRUTH$bo, TRUTH$lt0)
set.seed(20261009)
y_obs <- plasmofit:::mat_exp_series(y0, TRUTH$cl, N_C, TRUTH$R, 0, TS, 12)
y_obs <- pmax(0, y_obs * 10^rnorm(length(TS), 0, 0.2))

cat(sprintf("n_c = %d, FFT length M = %d, %d observation times\n", N_C, M, length(TS)))
cat("compiling both paths (a few minutes)...\n"); flush.console()
mod <- stan_model(file = "_scripts/conv-series-grad.stan", isystem = PKG,
                  model_name = "conv_grad")

mk <- function(use_conv) {
    d <- list(n_ts = length(TS), ts = TS, y_obs = y_obs, n_c = N_C, M = M,
              dt = 12, use_conv = use_conv)
    suppressMessages(sampling(mod, data = d, chains = 0))
}
f_conv <- mk(1L); f_mat <- mk(0L)

## unconstrained parameter vectors to test at: the truth plus random draws
up <- function(p) unconstrain_pars(f_conv, p)
pars <- list(list(cycle_length = 44, b_shape = 60, b_offset = 0.3,
                  log10_total0 = 5.5, R = 8, sigma = 0.2))
set.seed(11)
for (i in 1:4) pars[[length(pars) + 1]] <-
    list(cycle_length = runif(1, 38, 49), b_shape = exp(runif(1, log(3), log(2000))),
         b_offset = runif(1), log10_total0 = runif(1, 4, 7),
         R = runif(1, 2, 40), sigma = runif(1, 0.1, 0.6))
U <- lapply(pars, up)

cat("\n=== 1. log_prob: conv_series vs mat_exp_series ===\n")
cat("Cells: log posterior density at identical unconstrained parameter values.\n",
    "The two are the same model, so this is an equality check.\n\n", sep = "")
lp <- t(vapply(U, function(u) c(conv = log_prob(f_conv, u), mat = log_prob(f_mat, u)), c(0, 0)))
lp <- cbind(lp, abs_diff = abs(lp[, "conv"] - lp[, "mat"]))
print(as.data.frame(lp), digits = 10)
cat(sprintf("\nworst |difference| = %.3e\n", max(lp[, "abs_diff"])))

cat("\n=== 2. grad_log_prob: conv_series vs mat_exp_series ===\n")
cat("Cells: max over the 6 unconstrained parameters of |gradient difference|,\n",
    "and the same relative to the gradient's own scale. mat_exp_series is the\n",
    "reference because its gradients come from Stan's matrix_exp.\n\n", sep = "")
g <- t(vapply(U, function(u) {
    a <- grad_log_prob(f_conv, u); b <- grad_log_prob(f_mat, u)
    c(max_abs = max(abs(a - b)), rel = max(abs(a - b)) / max(abs(b)),
      grad_scale = max(abs(b)))
}, c(0, 0, 0)))
print(as.data.frame(g), digits = 6)
cat(sprintf("\nworst relative gradient difference = %.3e\n", max(g[, "rel"])))

cat("\n=== 3. grad_log_prob vs central finite differences ===\n")
cat("Cells: max over parameters of |analytic - finite difference| relative to\n",
    "the gradient scale, for conv_series. This catches an error the two paths\n",
    "could share. h is 1e-4 on the unconstrained scale.\n\n", sep = "")
fd <- function(fit, u, h = 1e-4) {
    vapply(seq_along(u), function(j) {
        up_ <- u; dn <- u; up_[j] <- u[j] + h; dn[j] <- u[j] - h
        (log_prob(fit, up_) - log_prob(fit, dn)) / (2 * h)
    }, 0)
}
d <- t(vapply(U, function(u) {
    a <- grad_log_prob(f_conv, u); b <- fd(f_conv, u)
    c(max_abs = max(abs(a - b)), rel = max(abs(a - b)) / max(abs(a)))
}, c(0, 0)))
print(as.data.frame(d), digits = 6)
cat(sprintf("\nworst relative difference from finite differences = %.3e\n", max(d[, "rel"])))

cat("\n=== verdict ===\n")
ok1 <- max(lp[, "abs_diff"]) < 1e-6
ok2 <- max(g[, "rel"]) < 1e-6
ok3 <- max(d[, "rel"]) < 1e-4            # finite differences are the loose one
cat(sprintf("log_prob agreement      : %s\n", if (ok1) "PASS" else "FAIL"))
cat(sprintf("gradient vs mat_exp     : %s\n", if (ok2) "PASS" else "FAIL"))
cat(sprintf("gradient vs finite diff : %s\n", if (ok3) "PASS" else "FAIL"))
if (ok1 && ok2 && ok3) {
    cat("\nconv_series differentiates correctly inside Stan. The FFT adjoint is\n",
        "sound and the function is safe to put in a fitted model.\n", sep = "")
} else {
    cat("\nDO NOT put conv_series in a fitted model until this passes.\n")
}
