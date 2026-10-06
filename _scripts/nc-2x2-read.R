## Thread 5: read the n_c misspecification 2x2.
##
## Generating n_c and fitting n_c are decoupled in wockner-schedule-sim.R, so
## the four cells separate estimator bias from misspecification bias:
##
##              fit 96              fit 192
##    gen  96   default             sim96_fit192
##    gen 192   sim192_fit96        sim192_fit192
##
## The diagonal is correctly specified. The off-diagonal is not.
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   Rscript --vanilla _scripts/nc-2x2-read.R

suppressPackageStartupMessages({library(dplyr); library(tibble); library(tidyr)})

CELLS <- c("default", "sim96_fit192", "sim192_fit96", "sim192_fit192")
GATE_RHAT <- 1.05

f <- list.files("_data", "^wock-schedsim-RES-.*[.]rds$", full.names = TRUE)
rows <- lapply(f, function(p) {
    x <- readRDS(p)
    if (!x$arm %in% CELLS) return(NULL)
    if (x$config != sprintf("%s-rep%d", x$arm, x$rep)) return(NULL)
    tibble(arm = x$arm, rep = x$rep,
           gen_n_c = as.integer(x$sim_n_c %||% 96L),
           fit_n_c = as.integer(x$fit_n_c %||% 96L),
           true_cl = x$true_cl,
           cl = mean(x$per_trial$cl_mean),
           bias = mean(x$per_trial$cl_mean) - x$true_cl,
           rhat = x$max_rhat, div = x$n_div, ess = x$min_ess)
})
r <- bind_rows(rows) |> filter(rep <= 3) |> arrange(arm, rep)
if (nrow(r) == 0) stop("no 2x2 results found")

cat("=== every cell, per replicate ===\n")
cat("Cells: cl is the mean posterior cycle_length over trials, in hours;\n",
    "bias is cl minus the simulated truth. gen/fit n_c are the n_c used to\n",
    "generate and to fit. converged is max R-hat < ", GATE_RHAT, ".\n\n", sep = "")
print(as.data.frame(r |> mutate(converged = rhat < GATE_RHAT)), digits = 4)

ok <- r |> filter(rhat < GATE_RHAT)
cat("\n=== bias by cell, converged replicates only ===\n")
summ <- ok |> group_by(arm, gen_n_c, fit_n_c) |>
    summarise(n = n(), mean_bias = mean(bias), min = min(bias), max = max(bias),
              .groups = "drop") |>
    arrange(gen_n_c, fit_n_c)
print(as.data.frame(summ), digits = 4)

cat("\n=== the two contrasts that matter ===\n")
g <- function(a) { v <- ok$bias[ok$arm == a]; if (length(v)) mean(v) else NA_real_ }
cat(sprintf("misspecification, UNDER (gen 192, fit 96): %s minus %s = %s h\n",
            fmt <- sprintf("%+.3f", g("sim192_fit96")),
            sprintf("%+.3f", g("sim192_fit192")),
            sprintf("%+.3f", g("sim192_fit96") - g("sim192_fit192"))))
cat(sprintf("misspecification, OVER  (gen 96, fit 192): %s minus %s = %s h\n",
            sprintf("%+.3f", g("sim96_fit192")),
            sprintf("%+.3f", g("default")),
            sprintf("%+.3f", g("sim96_fit192") - g("default"))))

cat("\nReading it:\n",
    " - OVER contrast NEGATIVE and near -3 h -> fitting above the true n_c\n",
    "   manufactures a downward cycle_length shift, and the real-data -3.20 h\n",
    "   is an artefact of over-large n_c rather than a correction.\n",
    " - OVER contrast near ZERO -> over-specifying n_c costs nothing, and the\n",
    "   real-data move is not explained this way.\n",
    " - UNDER contrast POSITIVE -> fitting below the true n_c inflates\n",
    "   cycle_length, which is the direction the real data would need if\n",
    "   n_c = 96 has been reading high all along.\n", sep = "")
