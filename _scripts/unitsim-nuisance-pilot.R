## PILOT, exploratory (2026-11-02, after the fixed-pin grid was read; Lucas
## asked whether any other model output could choose n_c / b_shape). Do the
## OTHER parameters differ across the 11 fixed-pin candidates of one dataset,
## and do they track |bias| in cycle_length better than elpd did (rho -0.23)?
## Only reps 1-3 of each window: 6 datasets, 66 fits. NO reading rule was set
## in advance and the correlations below are with |bias|, which is known only
## because this is simulated; a statistic is usable on real data only if it
## has an external target. Full-grid version would read all 120 datasets.
##
##   Rscript --vanilla _scripts/unitsim-nuisance-pilot.R > _data/unitsim-read-nuisance-pilot-<date>.txt

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
suppressPackageStartupMessages({library(rstan); library(dplyr)})
D <- "_data/unitsim/study/"
sfx <- c("", "-pinx2", "-pinhalf", "-f288", "-pinx2-f288", "-pinhalf-f288",
         "-pinmatch-f288", "-f512", "-pinx2-f512", "-pinhalf-f512", "-pinmatch-f512")
out <- list()
for (w in c(168L, 192L)) for (rep in 1:3) for (s in sfx) {
    nm <- sprintf("w%d-g384-fix-rep%d%s-ptruth", w, rep, s)
    f <- sprintf("%sfit-%s.rds", D, nm); if (!file.exists(f)) next
    fit <- readRDS(f); m <- as.matrix(fit)
    tr <- readRDS(sprintf("%schain-%s-c1.rds", D, nm))$truth
    out[[length(out) + 1]] <- tibble(
        window = w, rep = rep, cell = if (s == "") "truth-f384" else sub("^-", "", s),
        cl = mean(m[, "cycle_length[1]"]), bias = cl - tr$cycle_length,
        sd_iRBC = mean(m[, "sd_iRBC[1]"]), R = mean(m[, "R[1]"]),
        lt0 = mean(m[, "log10_total0[1]"]), b_off = mean(m[, "b_offset[1]"]),
        true_sd = unname(tr$sd_iRBC["MMV048_PIB.CH1"]), true_R = tr$R, true_lt0 = tr$log10_total0,
        true_boff = tr$b_offset)
}
x <- bind_rows(out)
cat("Truth: sd_iRBC", round(x$true_sd[1], 3), " R", round(x$true_R[1], 2),
    " log10_total0", round(x$true_lt0[1], 3), " b_offset", round(x$true_boff[1], 3), "\n\n")
x |> group_by(window, cell) |>
    summarise(n = n(), bias = mean(bias), sd_iRBC = mean(sd_iRBC), R = mean(R),
              lt0 = mean(lt0), b_off = mean(b_off), .groups = "drop") |>
    mutate(across(where(is.double), \(v) round(v, 3))) |> as.data.frame() |> print(row.names = FALSE)
cat("\nCorrelation over the 11 candidates within a dataset (mean over 6 datasets):\n")
x |> group_by(window, rep) |> summarise(
    r_sd = cor(sd_iRBC, abs(bias)), r_R = cor(R, abs(bias)),
    r_lt0 = cor(lt0, abs(bias)), r_boff = cor(b_off, abs(bias)), .groups = "drop") |>
    summarise(across(starts_with("r_"), \(v) round(mean(v), 3))) |> as.data.frame() |> print(row.names = FALSE)
