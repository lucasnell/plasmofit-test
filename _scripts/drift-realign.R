# Thread 6's drift measurement, with the group indices realigned by NAME.
#
# WHY THIS IS NEEDED. The drift fits were produced after wockner-fit.R's
# data-building step was edited -- a `group_by(id) |> arrange(time)` added to
# build the Design A hold-out mask -- which changed the ORDER of the group
# factor levels. The level SETS are identical, so every fit is internally
# correct and every quantity averaged over groups is unaffected (cycle_length
# agrees to 0.06 h). But `sd_iRBC[6]` means a different (trial, cohort) group
# in the two fits, so comparing them entry by entry compares unrelated
# quantities: it reported mean |z| of 25 and a max of 272, against a null of
# 0.85, which is not drift but a relabelling.
#
# The saved data lists carry attr(d, "levels"), so the mapping is recoverable
# exactly and no refit is needed. Each group-indexed parameter is permuted
# back into the reference fit's order before the comparison.
#
# The ordering change is cosmetic for the science and fatal for entry-wise
# comparison. See claude/gotchas.md.
#
#     srun -N 1 -n 1 -c 4 --mem=48G Rscript --vanilla \
#         _scripts/drift-realign.R 2>&1 | tee _data/drift-realign.log

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
suppressPackageStartupMessages({ library(tidyverse); library(rstan) })

## which group each parameter block is indexed by
## b_off_vec is array[n_grp_init] unit_vector[2], so it is indexed [i, j] and
## only the FIRST index is the group -- which is why the first version of this
## script, keying on block length, silently skipped it and left max |z| at 53.
PAR_GRP <- c(log10_total0 = "inoc", b_shape = "inoc", b_offset = "inoc",
             b_off_vec = "inoc",
             R = "trial", logit_R = "trial", eta_R = "trial",
             cycle_length = "trial", logit_cl = "trial", eta_cl = "trial",
             sd_iRBC = "obs_error", z_sd_iRBC = "obs_error")

PAIRS <- tribble(
    ~cmp, ~new,                    ~old,              ~what,
    "D1", "np_wide_total0-drift",  "np_wide_total0",  "rebuild only, b_shape free",
    "D2", "np_bs400-drift",        "np_bs400",        "rebuild only, b_shape pinned")

summ <- function(cfg) {
    d <- read_rds(sprintf("_data/wock-data-%s.rds", cfg))
    f <- read_rds(sprintf("_data/wock-fit-%s.rds", cfg))
    ss <- rstan::summary(f)$summary
    rm(f); invisible(gc())
    list(tab = tibble(par = rownames(ss), mean = ss[, "mean"],
                      se = ss[, "se_mean"]),
         levels = attr(d, "levels"))
}

## realign `new` onto `old`'s group ordering, by level name
realign <- function(new, old) {
    out <- new$tab
    for (p in names(PAR_GRP)) {
        g <- PAR_GRP[[p]]
        lo <- old$levels[[g]]; ln <- new$levels[[g]]
        if (is.null(lo) || is.null(ln)) next
        stopifnot(setequal(lo, ln))
        perm <- match(lo, ln)             # old group i is new group perm[i]
        rows <- which(grepl(sprintf("^%s\\[", p), out$par))
        if (!length(rows)) next
        ## Parse the leading index and keep any trailing one, so a matrix
        ## parameter is permuted on its group index alone. Keying on block
        ## length instead would skip b_off_vec, which has 28 entries for 14
        ## groups.
        nmr  <- out$par[rows]
        i1   <- as.integer(sub("^[^\\[]+\\[([0-9]+).*$", "\\1", nmr))
        rest <- sub("^[^\\[]+\\[[0-9]+", "", nmr)
        stopifnot(max(i1) == length(perm))
        ## for each stored row, which row should supply its value
        want <- match(paste0(perm[i1], rest), paste0(i1, rest))
        stopifnot(!anyNA(want))
        out$mean[rows] <- out$mean[rows][want]
        out$se[rows]   <- out$se[rows][want]
    }
    out
}

cat("=== reading ===\n")
need <- unique(c(PAIRS$new, PAIRS$old))
S <- map(set_names(need), \(c) { cat("  ", c, "\n"); summ(c) })

cat("\nCells: over the scalar entries both fits share, excluding lp__ and\n")
cat("log_lik. `z` is the difference in posterior means over the two runs'\n")
cat("combined Monte Carlo standard error. `raw` is as the fits are stored;\n")
cat("`realigned` maps each group-indexed parameter back onto the reference\n")
cat("fit's level ORDER by name first. The null (identical code and data,\n")
cat("sampler seed alone differing) is mean |z| 0.85 -- a drift at or near\n")
cat("that is no more than one run differs from another.\n\n")

out <- map(seq_len(nrow(PAIRS)), \(i) {
    p <- PAIRS[i, ]
    base <- S[[p$old]]$tab
    one <- function(tb, lab) {
        j <- inner_join(tb, base, by = "par", suffix = c("_n", "_o")) |>
            filter(par != "lp__", !startsWith(par, "log_lik"),
                   is.finite(se_n), is.finite(se_o), se_n > 0, se_o > 0) |>
            mutate(z = (mean_n - mean_o) / sqrt(se_n^2 + se_o^2))
        tibble(cmp = p$cmp, what = p$what, align = lab, n = nrow(j),
               mean_abs_z = mean(abs(j$z)), max_abs_z = max(abs(j$z)),
               frac_gt2 = mean(abs(j$z) > 2))
    }
    bind_rows(one(S[[p$new]]$tab, "raw"),
              one(realign(S[[p$new]], S[[p$old]]), "realigned"))
}) |> list_rbind()

out |> mutate(across(where(is.numeric), \(x) round(x, 3))) |> print(width = Inf)

cat("\n=== verdict ===\n")
for (c in unique(out$cmp)) {
    r <- out |> filter(cmp == c, align == "realigned")
    v <- if (r$mean_abs_z <= 1.15 * 0.846) {
        "PASS -- rebuild drift is no more than the sampler seed produces"
    } else {
        sprintf("EXCESS x%.2f over the null", r$mean_abs_z / 0.846)
    }
    cat(sprintf("  %s  %-28s mean |z| %.2f   %s\n", c, r$what, r$mean_abs_z, v))
}
