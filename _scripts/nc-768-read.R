## Extend model A's n_c profile to 768 and apply the pre-registered rule.
##
## Thread 15. Reads the finished screen (_data/decay-law-unit*.rds, model A at
## n_c 96/192/384) and the extension (_data/decay-law-a768-unit*.rds, model A
## at 768) and reports the gain at each step.
##
## The rule was fixed before the job was submitted and is repeated in
## _scripts/decay-law-768.sh:
##   summed gain 384->768 < +3  AND  < 9/14 units positive -> plateau, no
##     structural rewrite is justified.
##   summed gain >= +8  AND  >= 10/14 units positive -> still climbing, the
##     Erlang family is fighting the data, an IPM is justified.
##   anything between -> ambiguous, decide on biology.
## Geometric decay of the earlier gains (+72.4 then +10.4, ratio 0.14)
## predicts +1.5, which falls in the plateau branch.
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   Rscript --vanilla _scripts/nc-768-read.R | tee _data/nc-768-read-$(date +%F).txt
##
## The saved output is tracked -- see the exceptions in _data/.gitignore --
## because numbers from it are quoted in the notes.

suppressPackageStartupMessages({library(dplyr); library(tidyr)})

f0 <- list.files("_data", "^decay-law-unit[0-9]+[.]rds$", full.names = TRUE)
f1 <- list.files("_data", "^decay-law-a768-unit[0-9]+[.]rds$", full.names = TRUE)
if (length(f0) == 0) stop("no decay-law-unit*.rds in _data/ -- run the screen first")
if (length(f1) == 0) stop("no decay-law-a768-unit*.rds in _data/ -- job not finished")

r <- bind_rows(lapply(c(f0, f1), readRDS)) |> filter(model == "A_sqrt")

## Same shape guard as decay-law-read.R, for the same reason: a short file is
## a leftover from an interrupted run and would silently drop a rung.
EXPECTED <- 4L
bad <- r |> count(unit, name = "rows") |> filter(rows != EXPECTED)
if (nrow(bad) > 0) {
    print(as.data.frame(bad))
    stop("unit(s) above do not have ", EXPECTED, " model-A rows (n_c 96, 192, ",
         "384, 768) -- interrupted or stale files.")
}
cat("units read:", n_distinct(r$unit), "of 14, all with", EXPECTED,
    "model-A rows\n\n")

lad <- r |> select(unit, n, n_c, ll) |>
    pivot_wider(names_from = n_c, values_from = ll, names_prefix = "nc") |>
    mutate(g96_192  = nc192 - nc96,
           g192_384 = nc384 - nc192,
           g384_768 = nc768 - nc384)

cat("=== model A as a likelihood profile over n_c, extended to 768 ===\n")
cat("Cells: maximised profile log-likelihood of model A, differenced across\n",
    "n_c WITHIN a unit, in log-likelihood units. Positive favours the higher\n",
    "n_c. Every rung has the same parameter count, so no complexity penalty\n",
    "applies. No priors anywhere. Aggregated over the 14 grp_init units that\n",
    "hold all 1130 observations; unpaired across units is meaningless here,\n",
    "so the summary is a sum and a sign count.\n\n", sep = "")
print(as.data.frame(lad), digits = 5)

s <- c(sum(lad$g96_192), sum(lad$g192_384), sum(lad$g384_768))
p <- c(sum(lad$g96_192 > 0), sum(lad$g192_384 > 0), sum(lad$g384_768 > 0))
cat(sprintf("\n  96 -> 192 : %+7.1f summed | %2d/%d units positive\n", s[1], p[1], nrow(lad)))
cat(sprintf(" 192 -> 384 : %+7.1f summed | %2d/%d units positive\n", s[2], p[2], nrow(lad)))
cat(sprintf(" 384 -> 768 : %+7.1f summed | %2d/%d units positive\n", s[3], p[3], nrow(lad)))
cat(sprintf("\nratio of successive gains: %.3f then %.3f\n", s[2]/s[1], s[3]/s[2]))
cat(sprintf("geometric prediction for 384 -> 768 was %+.1f\n", s[1] * (s[2]/s[1])^2))

cat("\n=== pre-registered verdict ===\n")
if (s[3] < 3 && p[3] < 9) {
    cat("PLATEAU. The gains are decaying as a geometric extrapolation predicts.\n",
        "The data's preferred dispersion is reachable inside the Erlang family,\n",
        "so the mesh/rate coupling is not actively distorting the fit. No\n",
        "structural rewrite is justified: pick n_c by elpd, and free the stage\n",
        "rates only if a dispersion ABOVE the Erlang floor is wanted.\n", sep = "")
} else if (s[3] >= 8 && p[3] >= 10) {
    cat("STILL CLIMBING. The gains are not decaying. The data want transit\n",
        "closer to deterministic than Erlang reaches at any affordable n_c,\n",
        "so the family is fighting them and the coupling is doing real damage.\n",
        "An IPM's free dispersion width is justified on evidence.\n", sep = "")
} else {
    cat("AMBIGUOUS by the pre-registered rule (needs < +3 and < 9/14 for a\n",
        "plateau, or >= +8 and >= 10/14 for still climbing). The decision goes\n",
        "back to biology. Report the numbers, do not pick a branch post hoc.\n", sep = "")
}
