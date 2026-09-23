# Is the hierarchical cycle_length bias partly a posterior-mean artifact?
#
# The bias reported in claude/CLAUDE.md ("Where the bias is, by elimination")
# is measured on posterior MEANS of cycle_length. cycle_length is a bounded,
# nonlinearly transformed parameter -- min_cl + (max_cl - min_cl) *
# inv_logit(eta_cl) on [35, 50] -- so its posterior can be skewed, and the
# mean of a skewed posterior sits away from its maximum. If the population
# posterior is right-skewed, part of the ~1.1 h attributed to the hierarchy is
# just the choice of summary.
#
# Post-hoc on the saved simulation fits, no refitting. Reports mean, median
# and mode of the population cycle_length on the hours scale, plus the same on
# the logit scale pushed through the transform, which separates skew created
# by the transform from skew already present in eta_cl.
#
#   srun -N 1 -n 1 -c 4 --mem=16G --time=00:30:00 \
#       Rscript --vanilla _scripts/wockner-schedule-sim-mode.R

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")

suppressPackageStartupMessages({
    library(plasmofit)
    library(tidyverse)
    library(rstan)
})

res_files <- list.files("_data", "^wock-schedsim-RES-.*[.]rds$", full.names = TRUE)
if (length(res_files) == 0) stop("no wock-schedsim-RES-*.rds found in _data/")
res <- map(res_files, read_rds)

ok <- map_dbl(res, "max_rhat") < 1.05
if (!all(ok)) {
    cat("excluding non-converged:",
        paste(map_chr(res[!ok], "config"), collapse = ", "), "\n\n")
}
res <- res[ok]

## kernel-density mode; Sheather-Jones bandwidth, falling back to nrd0 when
## the draws are too few or too tied for it
post_mode <- function(x) {
    bw <- tryCatch(bw.SJ(x), error = function(e) bw.nrd0(x))
    dd <- density(x, bw = bw, n = 4096)
    dd$x[which.max(dd$y)]
}

out <- map(res, \(r) {
    f <- read_rds(sprintf("_data/wock-schedsim-fit-%s.rds", r$config))

    ## Population location on the hours scale, as the bias numbers define it:
    ## the trial-averaged cycle_length within each draw. (Not
    ## transform(mu_logit_cl) -- see "The no_pool / pooled_cl cycle-length
    ## offset" for why that differs.)
    cl_draws <- rstan::extract(f, "cycle_length")[[1]]
    cl_pop <- rowMeans(cl_draws)

    ## The same location on the logit scale, transformed after summarizing,
    ## so transform-induced skew can be told from skew in eta_cl itself.
    mu <- as.numeric(rstan::extract(f, "mu_logit_cl")[[1]])
    tr <- function(x) 35 + (50 - 35) / (1 + exp(-x))
    rm(f); invisible(gc())

    tibble(config = r$config, arm = r$arm, rep = r$rep, true_cl = r$true_cl,
           mean_h = mean(cl_pop), median_h = median(cl_pop),
           mode_h = post_mode(cl_pop),
           skew = mean((cl_pop - mean(cl_pop))^3) / sd(cl_pop)^3,
           mu_mean_tr = tr(mean(mu)), mu_mode_tr = tr(post_mode(mu)))
}) |> list_rbind()

cat("=== population cycle_length: mean vs median vs mode ===\n")
out |>
    mutate(bias_mean = mean_h - true_cl,
           bias_mode = mode_h - true_cl,
           mean_minus_mode = mean_h - mode_h) |>
    select(config, true_cl, mean_h, median_h, mode_h, mean_minus_mode,
           bias_mean, bias_mode, skew) |>
    mutate(across(where(is.numeric), \(x) round(x, 3))) |>
    print(n = Inf, width = Inf)

cat(sprintf("\n  mean(mean - mode) across fits: %+.3f h\n",
            mean(out$mean_h - out$mode_h)))
cat(sprintf("  mean bias using posterior means: %+.3f h\n",
            mean(out$mean_h - out$true_cl)))
cat(sprintf("  mean bias using posterior modes: %+.3f h\n",
            mean(out$mode_h - out$true_cl)))
cat(sprintf("  mean posterior skew: %+.3f\n", mean(out$skew)))

cat("\n=== transform applied after summarizing mu_logit_cl ===\n")
out |>
    mutate(mu_mean_bias = mu_mean_tr - true_cl,
           mu_mode_bias = mu_mode_tr - true_cl) |>
    select(config, mu_mean_tr, mu_mode_tr, mu_mean_bias, mu_mode_bias) |>
    mutate(across(where(is.numeric), \(x) round(x, 3))) |>
    print(n = Inf)

cat("\n  If mode and mean agree, the bias is not a summary artifact and the\n",
    "  ~1.1 h attributed to the hierarchy stands. If the mode sits at the\n",
    "  truth and only the mean is high, the posterior is skewed and the bias\n",
    "  is partly the choice of summary.\n", sep = "")

write_rds(out, "_data/wock-schedsim-mode.rds")
cat("\nwrote _data/wock-schedsim-mode.rds\n")
