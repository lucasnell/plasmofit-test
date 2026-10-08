## Aggregate the decay-law test across grp_init units.
##
## Model A is the sqrt law (Erlang chain alone, slope set by n_c). Model B
## adds a linear-in-cycles component via between-parasite heterogeneity in
## cycle duration. B nests A at sigma = 0, so B's log-likelihood can only be
## higher at the same n_c; the comparison is therefore made at MATCHED
## parameter count, with A spending one parameter on choosing its n_c from
## the grid exactly as B spends one on sigma.
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   Rscript --vanilla _scripts/decay-law-read.R | tee _data/decay-law-read-$(date +%F).txt
##
## The saved output is tracked -- see the exceptions in _data/.gitignore --
## because numbers from it are quoted in the notes.

suppressPackageStartupMessages({library(dplyr); library(tibble); library(readr)})

f <- list.files("_data", "^decay-law-unit[0-9]+[.]rds$", full.names = TRUE)
if (length(f) == 0) stop("no decay-law-unit*.rds in _data/")
r <- bind_rows(lapply(f, readRDS))

## Guard against a partial file. The job writes five rows per unit -- model A
## at n_c 96/192/384 and model B at 192/384 -- and only at the very end. A
## file with fewer rows is a leftover from an interrupted or differently
## configured run, and silently comparing it would drop whole rungs from that
## unit. This happened once: a killed smoke test left a 3-row unit07.
EXPECTED <- 5L
bad <- r |> count(unit, name = "rows") |> filter(rows != EXPECTED)
if (nrow(bad) > 0) {
    print(as.data.frame(bad))
    stop("unit(s) above do not have ", EXPECTED, " rows -- interrupted or ",
         "stale files. Delete them and re-run the task, or fix before trusting ",
         "the comparison.")
}
cat("units read:", length(unique(r$unit)), "of 14, all with",
    EXPECTED, "rows\n\n")

A <- r |> filter(model == "A_sqrt")    |> group_by(unit) |>
     slice_max(ll, n = 1, with_ties = FALSE) |> ungroup()
B <- r |> filter(model == "B_linear")  |> group_by(unit) |>
     slice_max(ll, n = 1, with_ties = FALSE) |> ungroup()

cmp <- A |> select(unit, n, ll_A = ll, A_n_c = n_c) |>
    left_join(B |> select(unit, ll_B = ll, B_n_c = n_c, sigma), by = "unit") |>
    mutate(d_ll = ll_B - ll_A)

cat("=== decay law, per grp_init unit ===\n")
cat("Cells: ll is the maximised profile log-likelihood with the observation\n",
    "sd profiled out; both models are matched at 6 parameters, A spending\n",
    "one on its choice of n_c and B one on sigma. d_ll > 0 favours the\n",
    "LINEAR law. sigma is the fitted log-sd of cycle duration between\n",
    "parasites, i.e. its coefficient of variation.\n\n", sep = "")
print(as.data.frame(cmp), digits = 4)

cat(sprintf("\ntotal d_ll (B - A): %+.2f over %d units | B wins %d, A wins %d\n",
            sum(cmp$d_ll), nrow(cmp), sum(cmp$d_ll > 0), sum(cmp$d_ll <= 0)))
cat(sprintf("median fitted sigma: %.4f | range %.4f-%.4f\n",
            median(cmp$sigma), min(cmp$sigma), max(cmp$sigma)))
cat(sprintf("A's chosen n_c: %s\n",
            paste(names(table(cmp$A_n_c)), table(cmp$A_n_c), sep = "x",
                  collapse = ", ")))

cat("\nReading it. Both models have the same parameter count, so d_ll is\n",
    "directly comparable and a couple of units of total d_ll is noise.\n",
    " - B clearly ahead AND consistent in sign -> the decay has a linear\n",
    "   component, and a 1-D IPM with a Gaussian kernel is the WRONG target\n",
    "   because it reproduces the sqrt law it would be replacing.\n",
    " - A ahead or a tie -> the sqrt law is adequate, and an IPM is then a\n",
    "   reasonable way to decouple the RATE from the mesh.\n",
    " - Mixed signs across units -> the design cannot tell, and the choice\n",
    "   has to be made on biology rather than on these data.\n", sep = "")

## ---------------------------------------------------------------------------
## Two diagnostics added 2026-10-08, after the first read of SLURM 29684.
##
## The first is a correction to the reading rule above, not a new test. B
## nests A at sigma = 0 (stated at the top of this file), and A's maximum is
## never at n_c = 96, so both models maximise over the same {192, 384} and
## ll_B >= ll_A must hold in exact arithmetic. A NEGATIVE d_ll is therefore
## impossible as evidence and can only be the optimiser stopping short. The
## rule as pre-registered counted those as "A wins" and read mixed signs as an
## uninformative design; in fact the negatives measure the noise floor, and a
## positive d_ll only counts as signal if it clears that floor.

m <- r |>
    select(unit, model, n_c, ll) |>
    tidyr::pivot_wider(names_from = model, values_from = ll) |>
    filter(!is.na(B_linear)) |>
    mutate(gap = B_linear - A_sqrt)

cat("\n\n=== nesting check: B - A at MATCHED n_c ===\n")
cat("Cells: difference in maximised profile log-likelihood, B minus A, at\n",
    "the same n_c, one row per unit x n_c. B contains A at sigma = 0, so\n",
    "every cell must be >= 0. Negative cells are optimiser shortfall and\n",
    "their magnitude is this test's noise floor.\n\n", sep = "")
print(as.data.frame(m |> arrange(gap) |> head(8)), digits = 4)

floor_ <- -min(m$gap)
cat(sprintf("\nnesting violations: %d of %d cells negative | worst %.3f | median violation %.4f\n",
            sum(m$gap < 0), nrow(m), min(m$gap),
            median(m$gap[m$gap < 0])))
cat(sprintf("noise floor (worst violation): %.3f log-likelihood units\n", floor_))

cat("\nPer-unit d_ll against that floor:\n")
sig <- cmp |> select(unit, n, d_ll) |> arrange(desc(d_ll)) |>
    mutate(clears_floor = d_ll > floor_)
print(as.data.frame(sig), digits = 4)
cat(sprintf("\nunits whose gain clears the noise floor: %d of %d\n",
            sum(sig$clears_floor), nrow(sig)))

## The second is a by-product, not what the job was built to measure: A on its
## own is a likelihood profile over n_c with no priors anywhere, so it is an
## independent check on the Bayesian n_c result.
cat("\n\n=== by-product: model A as a likelihood profile over n_c ===\n")
cat("Cells: maximised profile log-likelihood of model A at each n_c, and the\n",
    "gain of 192 and of 384 over 96, per unit, in log-likelihood units.\n",
    "Positive means the higher n_c fits better. No priors are involved.\n\n",
    sep = "")
lad <- r |> filter(model == "A_sqrt") |>
    select(unit, n, n_c, ll) |>
    tidyr::pivot_wider(names_from = n_c, values_from = ll,
                       names_prefix = "nc") |>
    mutate(g192 = nc192 - nc96, g384 = nc384 - nc96, g384_192 = nc384 - nc192)
print(as.data.frame(lad), digits = 5)
cat(sprintf("\nsummed over %d units: 192 over 96 = %+.1f | 384 over 96 = %+.1f | 384 over 192 = %+.1f\n",
            nrow(lad), sum(lad$g192), sum(lad$g384), sum(lad$g384_192)))
cat(sprintf("units where 192 beats 96: %d/%d | where 384 beats 192: %d/%d\n",
            sum(lad$g192 > 0), nrow(lad), sum(lad$g384_192 > 0), nrow(lad)))

## Power. The handoff pre-registered sigma ~ 0.033 as the between-parasite CV
## that would reproduce n_c = 192's spread under a linear law. The maximised
## gain is by definition >= the gain at any particular sigma, so each unit's
## d_ll is an UPPER BOUND on what a linear component of exactly that size
## could have bought. No refitting is needed to say this.
TARGET_SIGMA <- 0.033
cat("\n\n=== power: what a linear component of the predicted size could buy ===\n")
cat(sprintf(paste0("Cells: d_ll is the maximised gain of B over A, which bounds",
    " from above the\ngain at the pre-registered sigma = %.3f. sigma_hat is",
    " where B's optimiser\nactually landed.\n\n"), TARGET_SIGMA))
pw <- cmp |> select(unit, n, sigma_hat = sigma, d_ll) |>
    mutate(near_target = abs(sigma_hat - TARGET_SIGMA) < 0.01) |>
    arrange(desc(sigma_hat))
print(as.data.frame(pw), digits = 4)
cat(sprintf("\nunits landing within 0.01 of the predicted sigma: %d | their gains: %s\n",
            sum(pw$near_target),
            paste(sprintf("%+.3f", pw$d_ll[pw$near_target]), collapse = ", ")))
cat(sprintf("upper bound on the gain from sigma = %.3f, excluding the one unit\n",
            TARGET_SIGMA))
cat(sprintf("that clears the noise floor: at most %+.3f, median %+.4f\n",
            max(sig$d_ll[!sig$clears_floor]), median(sig$d_ll[!sig$clears_floor])))
