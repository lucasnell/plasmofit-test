## External anchors on n_c and b_shape (2026-11-02; extended the same day
## once Lucas supplied the papers, _data/new-info/). Both are DERIVED
## QUANTITIES of the fitted model with independently measurable counterparts,
## unlike elpd, which only scores fit to these data.
##
## n_c: in the Erlang chain a lineage crosses n_c exponential compartments per
##   cycle, so ONE PARASITE'S cycle duration has CV = 1 / sqrt(n_c)
##   (Erlang(k, rate): mean k/rate, variance k/rate^2), redrawn independently
##   each cycle. Counterpart: WITHIN-LINEAGE, cycle-to-cycle variability.
## b_shape: the t = 0 age distribution is Beta(b, b) over one cycle, sd
##   1 / sqrt(4 (2b + 1)) cycles. Counterpart: the spread of the inoculum's
##   ages, e.g. the window over which 99% of the first wave bursts.
##
## WHICH MEASUREMENT MATCHES. Published variability in cycle period is of two
## kinds, and only one of them is what n_c represents:
##   persistent  each parasite keeps its own period, so the population's
##               spread grows LINEARLY in cycles. Smith et al. (2020) assume
##               exactly this ("each parasite has a constant natural IEC
##               frequency"; they state they do not model cycle-over-cycle
##               variability), and Park et al.'s founder variance sigma_i is
##               the same thing measured directly. THE MODEL HAS NO SUCH TERM.
##   per cycle   redrawn each cycle, so spread grows as sqrt(cycles). This is
##               the Erlang chain, and Park et al.'s sibling variance sigma_k
##               (progeny of ONE founder schizont) measures it.
## So Smith's CoV must NOT be read as 1 / sqrt(n_c); Park's sigma_k may be.
## Both are in vitro. Park's strain is 3D7HT-GFP, the trials' 3D7 background.
##
##   Rscript --vanilla _scripts/unitsim-external-anchors.R > _data/unitsim-read-anchors-<date>.txt

CL <- 42.844
nc_from_cv <- function(cv) 1 / cv^2
cat("Single-unit truth cycle_length 42.844 h. CV = sd / mean of one cycle's duration.\n\n")

cat("1. WHAT THE MODEL IMPLIES\n\n n_c -> CV of one parasite's cycle duration, and its sd in hours\n\n")
for (nc in c(96, 192, 288, 384, 512, 768, 1024))
    cat(sprintf("   n_c %5d  CV %6.2f%%  sd %5.2f h\n", nc, 100 / sqrt(nc), CL / sqrt(nc)))
cat("\n b_shape -> sd of the t = 0 age distribution, and the window holding 99%\n\n")
for (b in c(5.86, 10.05, 14.00, 20.10, 40.20, 46.41, 400))
    cat(sprintf("   b_shape %6.2f  sd %5.3f cycles = %5.2f h   99%% window %5.2f h\n",
                b, sqrt(1 / (4 * (2 * b + 1))), CL * sqrt(1 / (4 * (2 * b + 1))),
                CL * diff(qbeta(c(0.005, 0.995), b, b))))
cat("\n   b_shape implied by a 99% bursting window of a given width\n\n")
for (h in c(6, 9, 12, 17, 24)) {
    f <- function(b) CL * diff(qbeta(c(0.005, 0.995), b, b)) - h
    cat(sprintf("   window %5.1f h  b_shape %7.1f\n", h, uniroot(f, c(1.5, 5000))$root))
}

## Park et al. (2020) bioRxiv 344812, Table 1: time-lapse microscopy of
## 3D7HT-GFP, burst-to-burst interval of individual parasites. tau = mean
## period; sigma_i = founder (persistent, between lineages); sigma_k = sibling
## (per cycle, within a founder's progeny); sigma = sqrt(sigma_i^2+sigma_k^2).
PARK <- data.frame(
    set = c("synchronized rep 1", "unsynchronized rep 1", "synchronized rep 2"),
    tau = c(54.0, 57.5, 57.3), sd_i = c(2.9, 4.9, 2.7), sd_k = c(1.7, 5.0, 2.7),
    sd_tot = c(3.4, 7.0, 3.8))
cat("\n\n2. PARK ET AL. (2020), DIRECT SINGLE-CELL, 3D7HT-GFP\n",
    "n_c is implied by the PER-CYCLE part (sibling, sigma_k) alone. The founder\n",
    "part (sigma_i) is persistent lineage heterogeneity the model cannot express;\n",
    "it is shown only to say how much of the total it is.\n\n", sep = "")
PARK$cv_k <- PARK$sd_k / PARK$tau; PARK$cv_tot <- PARK$sd_tot / PARK$tau
PARK$n_c_from_k <- nc_from_cv(PARK$cv_k)
PARK$n_c_from_tot <- nc_from_cv(PARK$cv_tot)
for (i in seq_len(nrow(PARK))) with(PARK[i, ], cat(sprintf(
    "   %-20s tau %4.1f h  sigma_k %3.1f h (CV %5.2f%%) -> n_c %6.0f | sigma_i %3.1f h | total %3.1f h (CV %5.2f%%) -> n_c %5.0f\n",
    set, tau, sd_k, 100 * cv_k, n_c_from_k, sd_i, sd_tot, 100 * cv_tot, n_c_from_tot)))
cat(sprintf("\n   Their standard-culture control (Fig. S4, ring emergence after one cycle):\n   48.7 +/- 5.0 h, CV %4.2f%% -> n_c %5.0f, but that spread also carries the\n   synchronisation's own width, so it bounds n_c from below only.\n",
            100 * 5.0 / 48.7, nc_from_cv(5.0 / 48.7)))

## Smith et al. (2020) Science, Table S6: CoV of a log-logistic distribution of
## FIXED per-parasite periods, fitted to population staging data.
SMITH <- data.frame(strain = c("3D7", "SA250", "FVO", "D6", "HB3"),
                    period = c(39, 54, 43, 36, 46),
                    cov = c(2.76e-2, 1.02e-1, 8.17e-2, 2.27e-3, 8.27e-2),
                    model_error = c(1.16, 0.241, 0.307, 0.647, 0.168))
cat("\n\n3. SMITH ET AL. (2020), TABLE S6 — PERSISTENT heterogeneity, NOT n_c\n",
    "Each parasite keeps one period, so this spread grows linearly in cycles and\n",
    "is not 1 / sqrt(n_c). Listed for the record, with the n_c it WOULD imply if\n",
    "(wrongly) read as per-cycle. Model error is their own fit error, for scale.\n\n", sep = "")
for (i in seq_len(nrow(SMITH))) with(SMITH[i, ], cat(sprintf(
    "   %-6s mean period %2d h  CoV %6.3f%%  (model error %5.3f)  [would imply n_c %8.0f]\n",
    strain, period, 100 * cov, model_error, nc_from_cv(cov))))
cat(sprintf("   %-6s mean period %.2f h  CoV %6.3f%%  (fibroblast reference)\n",
            "circad", 24.93, 8.45))
cat("\n   3D7 at +/- 1 h on its mean period: CoV 3.65% (38 h), 2.76% (39 h), 4.20% (40 h).\n",
    "   3D7 has the LARGEST model error of the five strains (1.16-1.59 against\n",
    "   0.168-0.647), so it is the strain their model fits worst.\n", sep = "")

cat("\n\n4. WOCKNER ET AL. (2020), the source data\n",
    "   Sine-wave life cycle 38.8 h (95% CI 38.3-39.2) overall; by subject mean\n",
    "   39.5 h, sd 4.04 over 128 subjects, range 30.6-55.0 (Supplementary Table 5).\n",
    "   That sd is BETWEEN HOSTS, not between parasites, so it anchors neither.\n",
    "   Samples were run in duplicate or triplicate and GEOMETRICALLY AVERAGED on\n",
    "   the log10 scale before release, and Supplementary Table 2 holds only the\n",
    "   averages, so the replicate spread that would anchor sd_iRBC is not in it.\n", sep = "")
