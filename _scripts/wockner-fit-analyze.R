# Local analysis of the fits produced by wockner-fit.R on the cluster.
#
# Expects wock-fit-RES-*.rds (and optionally wock-fit-LOO-*.rds and
# wock-fit-*.rds) to have been copied into _data/.
#
# Note on which diagnostics to trust here. fit_cond (condition number) is
# reported because it is cheap, but it has already failed to flag real trouble
# in this model, so do not read much into it on its own. R-hat, divergences and
# the per-chain parameter means are what actually caught the bimodality.


suppressPackageStartupMessages({
    library(tidyverse)
    library(rstan)
    library(plasmofit)
})

par_names <- c("b_shape", "b_offset", "R", "log10_total0", "cycle_length",
               "sd_iRBC")

res_files <- list.files("_data", "^wock-fit-RES-.*[.]rds$", full.names = TRUE)
if (length(res_files) == 0) stop("no wock-fit-RES-*.rds found in _data/")
res <- set_names(map(res_files, read_rds),
                 str_remove_all(basename(res_files), "^wock-fit-RES-|[.]rds$"))


# ------------------------------------------------------------------------ #
# Headline comparison
# ------------------------------------------------------------------------ #

summary_tbl <- imap(res, \(r, nm) {
    tibble(config = nm,
           max_cl = r$max_cl %||% NA_real_,
           sd_bs_cl = r$sd_bs_cl %||% NA_real_,
           center_cl = r$center_cl,
           divergences = r$n_divergent,
           leapfrog = round(r$leapfrog, 1),
           cond = round(r$cond, 2),
           max_rhat = round(max(r$diagnostics$rhat, na.rm = TRUE), 4),
           min_ess = round(min(r$diagnostics$ess_bulk, na.rm = TRUE)),
           sec_per_grad = round(r$sec_per_grad, 5),
           hours = round(max(rowSums(r$elapsed)) / 3600, 2))
}) |> list_rbind()

cat("=== configurations ===\n")
print(summary_tbl, width = Inf)


# ------------------------------------------------------------------------ #
# Chain agreement. Chains splitting into groups by lp__ is the signature of
# multimodality, which is what R-hat is reacting to when it blows up.
# ------------------------------------------------------------------------ #

cat("\n=== lp__ by chain (clusters here mean separated modes) ===\n")
iwalk(res, \(r, nm) cat(sprintf("  %-14s %s\n", nm,
                                paste(sprintf("%9.1f", r$lp_by_chain), collapse = " "))))


# ------------------------------------------------------------------------ #
# Funnel check: divergences concentrating where a scale parameter is small
# means the geometry, not the step size, is the problem.
# ------------------------------------------------------------------------ #

cat("\n=== divergences vs sigma (divergent / ok means) ===\n")
iwalk(res, \(r, nm) {
    cat("  ", nm, "\n", sep = "")
    if (is.null(r$div_sigma)) {
        cat("    (no divergences)\n")
    } else {
        print(round(r$div_sigma, 4))
    }
})


# ------------------------------------------------------------------------ #
# Per-chain cycle_length: the direct read on whether chains agree about the
# mode, and whether any trial is pinned against a bound.
# ------------------------------------------------------------------------ #

fit_files <- list.files("_data", "^wock-fit-[^R].*[.]rds$", full.names = TRUE)
fit_files <- fit_files[!str_detect(basename(fit_files), "^wock-fit-(RES|LOO)-")]

if (length(fit_files) > 0) {
    cat("\n=== cycle_length by chain, and per-trial spread ===\n")
    walk(fit_files, \(fl) {
        nm <- str_remove_all(basename(fl), "^wock-fit-|[.]rds$")
        f <- read_rds(fl)
        cl <- rstan::extract(f, "cycle_length", permuted = FALSE)
        ch_means <- map_dbl(seq_len(dim(cl)[2]), \(i) mean(colMeans(cl[, i, ])))
        per_trial <- colMeans(as.matrix(f, pars = "cycle_length"))
        cat(sprintf("  %-14s chains: %s | trials %.2f-%.2f h (spread %.2f)\n",
                    nm, paste(sprintf("%6.2f", ch_means), collapse = " "),
                    min(per_trial), max(per_trial),
                    diff(range(per_trial))))
    })
}


# ------------------------------------------------------------------------ #
# Model comparison, if log_lik was written
# ------------------------------------------------------------------------ #

loo_files <- list.files("_data", "^wock-fit-LOO-.*[.]rds$", full.names = TRUE)

if (length(loo_files) > 0 && requireNamespace("loo", quietly = TRUE)) {
    loos <- set_names(map(loo_files, read_rds),
                      str_remove_all(basename(loo_files),
                                     "^wock-fit-LOO-|[.]rds$"))
    cat("\n=== loo ===\n")
    imap(loos, \(l, nm) {
        map(names(l), \(by) {
            k <- loo::pareto_k_values(l[[by]])
            tibble(config = nm, by = by,
                   elpd = round(l[[by]]$estimates["elpd_loo", "Estimate"], 1),
                   se = round(l[[by]]$estimates["elpd_loo", "SE"], 1),
                   units = length(k),
                   bad_k = sum(k > 0.7))
        }) |> list_rbind()
    }) |> list_rbind() |> print()

    # Comparisons are only meaningful within a grouping, and only across models
    # fitted to the same observations.
    for (by in c("observation", "trial")) {
        ll <- map(loos, by) |> compact()
        if (length(ll) > 1) {
            cat("\n  loo_compare, by =", by, "\n")
            print(loo::loo_compare(ll))
        }
    }
    cat("\n  A large bad_k count at by = 'trial' means the importance sampling\n",
        "  broke down, not that the model is bad. That comparison then needs\n",
        "  K-fold refitting rather than this approximation.\n", sep = "")
}
