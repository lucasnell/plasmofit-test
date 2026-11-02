## External anchors on n_c and b_shape (2026-11-02, Lucas's question: is there
## any other model output that could choose them?). Both are DERIVED
## QUANTITIES of the fitted model that have independently measurable
## counterparts, unlike elpd, which only scores fit to these data.
##
## n_c: in the Erlang chain a lineage traverses n_c exponential compartments
##   per cycle, so the cycle duration of one parasite has CV = 1 / sqrt(n_c)
##   (Erlang(k, rate) has mean k/rate, variance k/rate^2). Counterpart:
##   cell-to-cell variability in asexual cycle period.
## b_shape: the t = 0 age distribution is Beta(b, b) over one cycle, sd
##   1 / sqrt(4 (2b + 1)) cycles. Counterpart: the spread of the inoculum's
##   ages, e.g. the window over which 99% of the first wave bursts.
##
##   Rscript --vanilla _scripts/unitsim-external-anchors.R > _data/unitsim-read-anchors-<date>.txt
CL <- 42.844
cat("cycle_length 42.844 h (the single-unit truth).\n\n")
cat("n_c -> implied CV of one parasite's cycle duration, and its sd in hours\n\n")
for (nc in c(96, 192, 288, 384, 512, 768, 1024))
    cat(sprintf(" n_c %5d  CV %6.2f%%  sd %5.2f h\n", nc, 100 / sqrt(nc), CL / sqrt(nc)))
cat("\nInverse: the n_c implied by a given CV\n\n")
for (cv in c(0.02, 0.03, 0.04, 0.05, 0.0845, 0.10))
    cat(sprintf(" CV %5.2f%%  n_c %7.0f\n", 100 * cv, 1 / cv^2))
cat("\nb_shape -> sd of the t = 0 age distribution, and the window holding 99%\n\n")
for (b in c(5.86, 10.05, 14.00, 20.10, 40.20, 46.41, 74.60, 400))
    cat(sprintf(" b_shape %6.2f  sd %5.3f cycles = %5.2f h   99%% window %5.2f h\n",
                b, sqrt(1 / (4 * (2 * b + 1))), CL * sqrt(1 / (4 * (2 * b + 1))),
                CL * diff(qbeta(c(0.005, 0.995), b, b))))
cat("\nInverse: the b_shape implied by a 99% bursting window of a given width\n\n")
for (h in c(6, 9, 12, 17, 24)) {
    f <- function(b) CL * diff(qbeta(c(0.005, 0.995), b, b)) - h
    cat(sprintf(" window %5.1f h  b_shape %7.1f\n", h, uniroot(f, c(1.5, 5000))$root))
}
