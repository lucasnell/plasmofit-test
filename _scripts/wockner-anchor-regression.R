# Regression test for the inoculum-anchored log10_total0 prior.
#
# The anchor (plasmofit 0.0.0.9008) adds five data fields and two parameters
# to all four Stan programs. With total0_anchor = 0 the new parameters are
# declared array[0] -- zero-sized -- and the prior statement falls through to
# the original
#
#     log10_total0 ~ normal(mean_log10_total0, sd_log10_total0);
#
# so a fit with the anchor off must be bit-identical to the same fit before
# the change, not merely similar. Anything less means a zero-sized parameter
# is perturbing the sampler, and no anchored result can be trusted.
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

cat("\n=== draw matrices ===\n")
if (nrow(a) != nrow(b)) {
    cat("  DIFFERENT NUMBER OF DRAWS -- cannot compare element by element\n")
    bit_identical <- FALSE
} else {
    bit_identical <- identical(a, b)
    cat("  identical():", bit_identical, "\n")
    d <- abs(a - b)
    cat(sprintf("  max abs difference: %.3g\n", max(d)))
    if (!bit_identical) {
        worst <- sort(apply(d, 2, max), decreasing = TRUE)[1:10]
        cat("  worst parameters:\n")
        for (nm in names(worst)) cat(sprintf("    %-24s %.3g\n", nm, worst[[nm]]))
    }
}

## Sampler behaviour, which would move before the draws do if something
## subtle changed in the gradient.
sp <- function(f) {
    s <- get_sampler_params(f, inc_warmup = FALSE)
    c(div = sum(sapply(s, \(x) sum(x[, "divergent__"]))),
      leapfrog = sum(sapply(s, \(x) sum(x[, "n_leapfrog__"]))),
      stepsize = mean(sapply(s, \(x) mean(x[, "stepsize__"]))))
}
cat("\n=== sampler ===\n")
print(rbind(reference = sp(ref), rebuilt = sp(new)))

cat("\n=== verdict ===\n")
if (bit_identical && length(only_ref) == 0 && length(only_new) == 0) {
    cat("  PASS -- anchor off reproduces the pre-change fit exactly.\n")
} else {
    cat("  FAIL -- the anchor change perturbs a fit it should not touch.\n",
        "  Do not interpret any anchored fit until this is understood.\n", sep = "")
}
