## What does fitting ABOVE the truth do to the OTHER parameters?
## (2026-11-09, Lucas, following "Choosing fit n_c from outside the data":
## the asymmetry argument was made on cycle_length alone.)
##
## The cleanest comparison in hand: the fixed-pin grid's truth-pin cells,
## gen n_c 384, b_shape FIXED at the truth 20.10, cycle_length prior centred
## at the truth, fitted at n_c 288 / 384 / 512 on identical data. b_shape
## cannot absorb anything (it is fixed), so whatever the wrong mesh costs has
## to appear in cycle_length, R, log10_total0, b_offset or sd_iRBC.
##
## Cells: posterior mean of each parameter, mean over the passing fits of a
## cell (58-60 of 60), against the truth. Paired: per replicate on identical
## data, both fits passing, mean (se) of the difference from the fit-384 fit.
## R and log10_total0 are known not to be identified by these data
## (2026-09-24), so read their columns as "where the misfit goes", not as
## recovery.
##
##   Rscript --vanilla _scripts/unitsim-ncother.R > _data/unitsim-read-ncother-<date>.txt

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
suppressPackageStartupMessages({library(rstan); library(dplyr); library(tidyr)})
options(width = 160)
D <- "_data/unitsim/study/"
PARS <- c(cycle_length = "cycle_length[1]", R = "R[1]", log10_total0 = "log10_total0[1]",
          b_offset = "b_offset[1]", sd_iRBC = "sd_iRBC[1]")

s <- bind_rows(lapply(Sys.glob(paste0(D, "sum-*.rds")), readRDS))
for (v in c("prior", "truth_set", "pin")) if (!v %in% names(s)) s[[v]] <- NA_character_
if (!"fit_nc" %in% names(s)) s$fit_nc <- NA_integer_
s <- s |> mutate(prior = coalesce(prior, "prod"), truth_set = coalesce(truth_set, "main"),
                 pin = coalesce(pin, "truth"), fit_nc = coalesce(fit_nc, 384L)) |>
    filter(truth_set == "main", gen_nc == 384L, method == "fix", prior == "truth",
           pin == "truth")
cat(nrow(s), "truth-pin fits (gen 384, b_shape fixed, prior at the truth);",
    sum(s$pass), "pass.\n\n")

tr <- readRDS(sprintf("%schain-%s-c1.rds", D, s$name[1]))$truth
TRUTH <- c(cycle_length = tr$cycle_length, R = tr$R, log10_total0 = tr$log10_total0,
           b_offset = tr$b_offset, sd_iRBC = unname(tr$sd_iRBC["MMV048_PIB.CH1"]))

rows <- lapply(seq_len(nrow(s)), function(i) {
    m <- as.matrix(readRDS(sprintf("%sfit-%s.rds", D, s$name[i])), pars = names(PARS))
    v <- colMeans(m[, PARS, drop = FALSE]); names(v) <- names(PARS)
    as_tibble(as.list(v))
})
x <- bind_cols(select(s, name, window, fit_nc, rep, pass), bind_rows(rows))

cat("Truth: ", paste(sprintf("%s %.3f", names(TRUTH), TRUTH), collapse = "; "), "\n\n", sep = "")
cat("1. PER CELL. Posterior mean of each parameter, mean over passing fits.\n\n")
xp <- filter(x, pass)
xp |> group_by(window, fit_nc) |>
    summarise(n = n(), across(all_of(names(PARS)), \(v) round(mean(v), 3)), .groups = "drop") |>
    as.data.frame() |> print(row.names = FALSE)

cat("\n2. BIAS PER CELL (posterior mean minus truth), same cells.\n\n")
xp |> group_by(window, fit_nc) |>
    summarise(n = n(), across(all_of(names(PARS)),
                              \(v) round(mean(v) - TRUTH[cur_column()], 3)), .groups = "drop") |>
    as.data.frame() |> print(row.names = FALSE)

cat("\n3. PAIRED, identical data: each fit minus the fit-384 fit of the same\n",
    "   replicate. mean (se) over pairs; both must pass.\n\n", sep = "")
ref <- filter(xp, fit_nc == 384L) |> select(window, rep, all_of(names(PARS)))
for (nc in c(288L, 512L)) {
    j <- inner_join(filter(xp, fit_nc == nc) |> select(window, rep, all_of(names(PARS))),
                    ref, by = c("window", "rep"), suffix = c("_a", "_b"))
    tab <- j |> group_by(window) |> summarise(pairs = n(), .groups = "drop")
    for (p in names(PARS)) {
        dd <- j[[paste0(p, "_a")]] - j[[paste0(p, "_b")]]
        agg <- tapply(dd, j$window, function(z) sprintf("%+.3f (%.3f)", mean(z), sd(z) / sqrt(length(z))))
        tab[[p]] <- unname(agg[as.character(tab$window)])
    }
    cat(sprintf("   fit %d - fit 384\n", nc))
    print(as.data.frame(tab), row.names = FALSE); cat("\n")
}
