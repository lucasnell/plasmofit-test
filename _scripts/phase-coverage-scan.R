# How well does the observation schedule cover the cycle, as a function of the
# cycle length being considered?
#
# 87.8% of within-series intervals are 12 h. A period COMMENSURATE with 12 h
# -- 24 h, 48 h -- makes every observation land on one of only a few phases,
# because t/P mod 1 repeats. At such a period the data carry no information
# about where in the cycle anything sits: any b_offset fits equally well.
#
# That matters because the estimates are biased TOWARD 48 h and the prior is
# centred there. If the likelihood is degenerate in phase at commensurate
# periods, those periods carry extra parameter-space volume, and a posterior
# MEAN integrates over volume. It would be a marginalisation effect of the
# kind the MAP check could not test, arising from the sampling grid rather
# than from any prior -- which would explain why widening priors never moved
# it.
#
# Cells computed per candidate period P over the real observation times:
#   n_eff_phase  1 / sum(p_i^2) over 36 phase bins -- the effective number of
#                distinct phases visited, 36 if perfectly uniform, small if
#                the samples pile onto a few phases
#   max_gap      largest gap between consecutive sorted phases, in cycles;
#                0.028 if uniform, large if the cycle is poorly covered
#   R            mean resultant length of exp(2*pi*i*phase); 0 if uniform,
#                1 if every sample sits at one phase
#
# Deterministic, no fitting, no data values used -- only the observation times.
#
#     Rscript --vanilla _scripts/phase-coverage-scan.R \
#         2>&1 | tee _data/phase-coverage-scan.log

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
suppressPackageStartupMessages(library(tidyverse))

tv <- read_csv("_data/wockner-cleaned.csv", col_types = "cccdcdd")$time
NB <- 36

stat <- function(P) {
    ph <- (tv / P) %% 1
    b <- table(cut(ph, breaks = seq(0, 1, length.out = NB + 1),
                   include.lowest = TRUE))
    p <- as.numeric(b) / sum(b)
    s <- sort(unique(ph))
    gaps <- c(diff(s), 1 - max(s) + min(s))
    tibble(P = P, n_eff_phase = 1 / sum(p^2), max_gap = max(gaps),
           R = sqrt(mean(cos(2 * pi * ph))^2 + mean(sin(2 * pi * ph))^2))
}

grid <- map(seq(35, 50, by = 0.02), stat) |> list_rbind()

cat("=== commensurate periods vs their neighbourhood ===\n")
cat("Cells as in the header. 12 h grid, so P = 24 and 48 are exactly\n")
cat("commensurate (48/12 = 4 samples per cycle, repeating).\n\n")
grid |> filter(P %in% c(36, 40, 42, 44, 45, 45.02, 46, 47, 47.5, 47.9, 48,
                        48.1, 48.5, 49, 50)) |>
    mutate(across(where(is.numeric), \(x) round(x, 3))) |> print(n = Inf)

cat("\n=== where coverage is worst over 35-50 h ===\n")
cat("Lowest effective number of distinct phases:\n")
grid |> slice_min(n_eff_phase, n = 8) |>
    mutate(across(where(is.numeric), \(x) round(x, 3))) |> print(n = Inf)

cat("\nAt the simulated truth 45.012 h vs at the prior centre 48 h:\n")
a <- stat(45.012); b <- stat(48)
cat(sprintf("  45.012 h: %.1f effective phases, max gap %.3f, R %.3f\n",
            a$n_eff_phase, a$max_gap, a$R))
cat(sprintf("  48.000 h: %.1f effective phases, max gap %.3f, R %.3f\n",
            b$n_eff_phase, b$max_gap, b$R))
write_rds(grid, "_data/wock-phase-coverage.rds")
