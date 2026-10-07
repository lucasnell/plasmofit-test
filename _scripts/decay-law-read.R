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
##   Rscript --vanilla _scripts/decay-law-read.R

suppressPackageStartupMessages({library(dplyr); library(tibble); library(readr)})

f <- list.files("_data", "^decay-law-unit[0-9]+[.]rds$", full.names = TRUE)
if (length(f) == 0) stop("no decay-law-unit*.rds in _data/")
r <- bind_rows(lapply(f, readRDS))
cat("units read:", length(unique(r$unit)), "of 14\n\n")

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
