# Regression test for the inoculum-anchored log10_total0 prior.
#
# The anchor (plasmofit 0.0.0.9008) adds five data fields and two parameters
# to all four Stan programs. With total0_anchor = 0 the new parameters are
# declared array[0] -- zero-sized -- and the prior statement falls through to
# the original
#
#     log10_total0 ~ normal(mean_log10_total0, sd_log10_total0);
#
# so a fit with the anchor off must describe the same posterior as the same
# fit before the change.
#
# NOT bit-identical draws. The reference predates the change and was sampled
# by a DIFFERENT COMPILED DSO; recompiling can reorder floating-point
# operations, and HMC amplifies a last-bit difference into an entirely
# different trajectory within a few leapfrog steps. Identical draws are
# unachievable here and demanding them only ever prints FAIL. The test that
# means something is posterior equivalence within Monte Carlo error.
#
# The structural checks still have to be exact: the parameter name sets must
# match, because a zero-sized array contributes no columns, and any extra or
# missing column means the anchor is not really switched off.
#
# Three fits, when all are available:
#   reference   pre-change code, seed A   _data/regress-ref-fit-default-rep2.rds
#   rebuilt     post-change code, seed A  _data/wock-schedsim-fit-default-rep2.rds
#   null run    post-change code, seed B  _data/wock-schedsim-fit-default-rep2-seed*.rds
#
# reference-vs-rebuilt is the TEST: it differs by the code change, and by the
# sampler trajectory, which a recompile alone changes.
# rebuilt-vs-null-run is the NULL: identical code and data, different seed,
# so every difference it shows is run-to-run Monte Carlo variation including
# the separate step-size and mass-matrix adaptation each run performs.
#
# Without the null run there is nothing to judge the test against: splitting
# one fit's chains holds the adaptation constant and so understates run-to-run
# spread. In that case this script reports INCONCLUSIVE rather than guessing.
# Produce the null run with
#
#   sbatch --array=2 --export=ALL,SCHEDSIM_FIT_SEED=271828183 \
#          _scripts/wockner-schedule-sim.sh
#
# Inputs, both written by wockner-schedule-sim.R:
#   _data/regress-ref-fit-default-rep2.rds  frozen pre-change fit
#   _data/wock-schedsim-fit-default-rep2.rds  rebuilt under 0.0.0.9008
#
# Run where the fits are; loads two full fits.
#   srun -N 1 -n 1 -c 4 --mem=24G Rscript --vanilla \
#       _scripts/wockner-anchor-regression.R

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")

suppressPackageStartupMessages({
    library(plasmofit)
    library(rstan)
})

ref_path <- "_data/regress-ref-fit-default-rep2.rds"
new_path <- "_data/wock-schedsim-fit-default-rep2.rds"
for (p in c(ref_path, new_path)) if (!file.exists(p)) stop("missing: ", p)

cat("plasmofit", as.character(packageVersion("plasmofit")), "\n")
cat("reference:", ref_path, format(file.mtime(ref_path)), "\n")
cat("rebuilt:  ", new_path, format(file.mtime(new_path)), "\n\n")

ref <- readRDS(ref_path)
new <- readRDS(new_path)

m_ref <- as.matrix(ref)
m_new <- as.matrix(new)

cat("=== dimensions ===\n")
cat(sprintf("  reference: %d draws x %d parameters\n", nrow(m_ref), ncol(m_ref)))
cat(sprintf("  rebuilt:   %d draws x %d parameters\n", nrow(m_new), ncol(m_new)))

## Zero-sized parameters contribute no columns, so the name sets should match
## exactly. Report any difference rather than quietly intersecting.
only_ref <- setdiff(colnames(m_ref), colnames(m_new))
only_new <- setdiff(colnames(m_new), colnames(m_ref))
cat("\n=== parameter names ===\n")
if (length(only_ref) == 0 && length(only_new) == 0) {
    cat("  identical\n")
} else {
    cat("  only in reference:", paste(only_ref, collapse = ", "), "\n")
    cat("  only in rebuilt:  ", paste(only_new, collapse = ", "), "\n")
}

shared <- intersect(colnames(m_ref), colnames(m_new))
a <- m_ref[, shared, drop = FALSE]
b <- m_new[, shared, drop = FALSE]

## Bit-identity is reported for the record, not used as the verdict: see the
## header. Expect FALSE whenever the DSO was rebuilt.
cat("\n=== draw matrices (informational) ===\n")
cat("  identical():", identical(a, b), "\n")

## ---------------------------------------------------------------------- #
## Posterior equivalence.
##
## Two independent samplers targeting the SAME posterior give posterior means
## that differ only by Monte Carlo error. For parameter p,
##
##     z_p = (mean_ref - mean_new) / sqrt(mcse_ref^2 + mcse_new^2)
##
## where mcse = posterior sd / sqrt(effective sample size). If the target is
## unchanged the z_p behave like standard normals, so |z| up to ~3 is ordinary
## and a handful above that across 208 parameters is expected. A changed
## target shows up as z values in the tens, and as a systematic shift rather
## than a scatter around zero.
##
## Posterior sds are compared as a ratio for the same reason: a changed prior
## on log10_total0 would widen or narrow it, which a mean-only test could miss
## if the shift happened to be small.
## ---------------------------------------------------------------------- #

## ESS and MCSE per parameter, from rstan's own summary so the ESS definition
## matches the diagnostics used everywhere else in this project.
summ <- function(f) {
    x <- rstan::summary(f)$summary
    data.frame(par = rownames(x), mean = x[, "mean"], sd = x[, "sd"],
               ess = x[, "n_eff"], rhat = x[, "Rhat"],
               row.names = NULL, stringsAsFactors = FALSE)
}
s_ref <- summ(ref)
s_new <- summ(new)
s_ref <- s_ref[match(shared, s_ref$par), ]
s_new <- s_new[match(shared, s_new$par), ]

mcse_ref <- s_ref$sd / sqrt(s_ref$ess)
mcse_new <- s_new$sd / sqrt(s_new$ess)
z <- (s_ref$mean - s_new$mean) / sqrt(mcse_ref^2 + mcse_new^2)
sd_ratio <- s_new$sd / s_ref$sd

cat("\n=== posterior equivalence ===\n")
cat(sprintf("  parameters compared: %d\n", length(shared)))
cat(sprintf("  max R-hat: reference %.4f, rebuilt %.4f\n",
            max(s_ref$rhat, na.rm = TRUE), max(s_new$rhat, na.rm = TRUE)))
cat(sprintf("  |z| on posterior means: median %.2f, 90%% %.2f, max %.2f\n",
            median(abs(z), na.rm = TRUE),
            quantile(abs(z), 0.9, na.rm = TRUE), max(abs(z), na.rm = TRUE)))
cat(sprintf("  z is signed: mean %+.3f (a changed target shifts this off 0)\n",
            mean(z, na.rm = TRUE)))
cat(sprintf("  |z| > 3: %d of %d (%.1f%%)\n", sum(abs(z) > 3, na.rm = TRUE),
            length(z), 100 * mean(abs(z) > 3, na.rm = TRUE)))
cat(sprintf("  posterior sd ratio (new/ref): median %.3f, range %.3f-%.3f\n",
            median(sd_ratio, na.rm = TRUE), min(sd_ratio, na.rm = TRUE),
            max(sd_ratio, na.rm = TRUE)))

ord <- order(abs(z), decreasing = TRUE)[1:10]
cat("\n  largest |z|:\n")
cat(sprintf("    %-22s %9s %9s %7s %7s %7s\n",
            "parameter", "ref mean", "new mean", "z", "sd rat", "min ESS"))
for (k in ord) {
    cat(sprintf("    %-22s %9.4g %9.4g %7.2f %7.3f %7.0f\n",
                s_ref$par[k], s_ref$mean[k], s_new$mean[k], z[k],
                sd_ratio[k], min(s_ref$ess[k], s_new$ess[k])))
}

## log10_total0 gets its own look: it is the parameter whose prior statement
## the change touches, so if anything moved it should move here first.
i_t0 <- grep("^log10_total0\\[", shared)
cat(sprintf("\n  log10_total0 specifically: max |z| %.2f, sd ratio %.3f-%.3f\n",
            max(abs(z[i_t0])), min(sd_ratio[i_t0]), max(sd_ratio[i_t0])))

## ---------------------------------------------------------------------- #
## Are the posterior WIDTHS the same?
##
## A flat band on sd_new/sd_ref is the wrong test: a posterior sd estimated
## from an autocorrelated chain carries Monte Carlo error that grows as
## effective sample size falls, so the same band is lax for a well-mixed
## parameter and impossibly tight for a poorly-mixed one.
##
## Instead measure |log(sd_1 / sd_2)| over parameters for the test pair and
## for the null pair, and ask whether the test is inflated. Both pairs use
## whole fits, so the two statistics are computed the same way and are
## directly comparable.
## ---------------------------------------------------------------------- #

null_path <- Sys.glob("_data/wock-schedsim-fit-default-rep2-seed*.rds")
have_null <- length(null_path) > 0

lr <- function(x, y) log(x / y)
test_lr <- lr(s_new$sd, s_ref$sd)

cat("\n=== posterior widths ===\n")
cat(sprintf("  test  (pre-change vs post-change): median |log sd ratio| %.4f, 95%% %.4f\n",
            median(abs(test_lr), na.rm = TRUE),
            quantile(abs(test_lr), 0.95, na.rm = TRUE)))

if (have_null) {
    null_path <- null_path[1]
    cat("  null run:", null_path, "\n")
    nullfit <- readRDS(null_path)
    s_nul <- summ(nullfit)
    s_nul <- s_nul[match(shared, s_nul$par), ]
    if (anyNA(s_nul$sd)) stop("null run does not carry the same parameters")

    null_lr <- lr(s_nul$sd, s_new$sd)
    z_null <- (s_new$mean - s_nul$mean) /
        sqrt((s_new$sd / sqrt(s_new$ess))^2 + (s_nul$sd / sqrt(s_nul$ess))^2)

    cat(sprintf("  null  (same code, seed %s):       median |log sd ratio| %.4f, 95%% %.4f\n",
                sub(".*seed([0-9]+).*", "\\1", null_path),
                median(abs(null_lr), na.rm = TRUE),
                quantile(abs(null_lr), 0.95, na.rm = TRUE)))
    width_infl <- quantile(abs(test_lr), 0.95, na.rm = TRUE) /
                  quantile(abs(null_lr), 0.95, na.rm = TRUE)
    cat(sprintf("  inflation (test / null) at the 95%%: %.2fx\n", width_infl))

    cat(sprintf("\n  means, null pair: |z| median %.2f, max %.2f, |z| > 3: %d of %d\n",
                median(abs(z_null), na.rm = TRUE), max(abs(z_null), na.rm = TRUE),
                sum(abs(z_null) > 3, na.rm = TRUE), length(z_null)))
    cat(sprintf("  means, test pair: |z| median %.2f, max %.2f, |z| > 3: %d of %d\n",
                median(abs(z), na.rm = TRUE), max(abs(z), na.rm = TRUE),
                sum(abs(z) > 3, na.rm = TRUE), length(z)))

    ## Thresholds fixed before the null run existed: the test pair must not be
    ## more than half again as spread as two runs of identical code.
    sd_ok <- is.finite(width_infl) && width_infl < 1.5
    rm(nullfit); invisible(gc())
} else {
    cat("  null run NOT FOUND -- nothing to judge the test against.\n")
    cat("  Splitting one fit's chains is not a substitute: it holds the\n")
    cat("  step-size and mass-matrix adaptation constant, and those differ\n")
    cat("  between runs, so it understates run-to-run spread.\n")
    width_infl <- NA_real_
    sd_ok <- NA
}

## Sampler behaviour. Differs freely between two runs of the same target;
## reported so a wild difference is visible, not used as the verdict.
sp <- function(f) {
    s <- get_sampler_params(f, inc_warmup = FALSE)
    c(div = sum(sapply(s, \(x) sum(x[, "divergent__"]))),
      leapfrog = sum(sapply(s, \(x) sum(x[, "n_leapfrog__"]))),
      stepsize = mean(sapply(s, \(x) mean(x[, "stepsize__"]))))
}
cat("\n=== sampler (informational) ===\n")
print(rbind(reference = sp(ref), rebuilt = sp(new)))

cat("\n=== verdict ===\n")
names_ok <- length(only_ref) == 0 && length(only_new) == 0
## Fixed before seeing any result: a shared target puts almost every |z| under
## 3 and leaves the signed mean near zero.
z_ok <- mean(abs(z) > 3, na.rm = TRUE) < 0.05 && abs(mean(z, na.rm = TRUE)) < 0.5

if (!names_ok) {
    cat("  FAIL -- parameter sets differ, so the anchor is not really off.\n")
} else if (!z_ok) {
    cat("  FAIL -- posterior means moved beyond Monte Carlo error.\n")
    cat("  Do not interpret any anchored fit until this is understood.\n")
} else if (is.na(sd_ok)) {
    cat("  INCONCLUSIVE -- parameter sets match and posterior means agree\n")
    cat("  within Monte Carlo error, which is the substantive check. Posterior\n")
    cat("  widths cannot be judged without the null run; produce it as shown in\n")
    cat("  the header and re-run. Draws are not bit-identical, which is expected\n")
    cat("  across a rebuilt DSO and is not evidence of anything.\n")
} else if (!sd_ok) {
    cat(sprintf("  FAIL -- posterior widths are %.2fx as spread as two runs of\n",
                width_infl))
    cat("  identical code, so the change moved the target.\n")
} else {
    cat("  PASS -- anchor off targets the same posterior as before the change.\n")
    cat("  Means agree within Monte Carlo error and widths are no more spread\n")
    cat("  than two runs of identical code. Draws are not bit-identical, which\n")
    cat("  is expected across a rebuilt DSO.\n")
}
