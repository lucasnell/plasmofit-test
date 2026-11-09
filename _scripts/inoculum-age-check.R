## Is the fitted t = 0 age distribution consistent with what the inoculum is?
## (2026-11-09, Lucas.) Every QIMR-B IBSM inoculum came from ONE donor unit of
## blood drawn from a volunteer infected by mosquito bite (13 days after the
## bites, 6 h after onset of high fever), leucocyte-depleted, frozen in
## Glycerolyte 57; thawed, washed and injected within 60 min (McCarthy et al.
## 2011, PLoS One 6:e21914; Cheng et al. 1997; Wockner et al. 2020). So the
## injected iRBCs are the donor's CIRCULATING iRBCs: mature stages were
## sequestered in the donor and are absent from venous blood, and only rings
## survive Glycerolyte freezing in any case (ATCC handling sheets). They are
## also IDENTICAL across all 14 units, up to dilution.
##
## The model's own sequestration curve (package make_log_y_vals; Archer et al.
## 2018 fitted to Kriek et al. 2003) gives the probability that an iRBC of age a
## (hours) still circulates: y(a) = ((1 - p2) + (1 - p1) u) / (1 + u),
## u = 10^(p4 (p3 - a)), p1 = 11.3869 / 467.6209, p2 = 1, p3 = 18.5802 h,
## p4 = 0.2242 / h. Half of a cohort has sequestered by 18.58 h.
##
## This script:
##  1. the circulating window implied by y(a);
##  2. the widest initial age distribution a circulating-only inoculum can have
##     (donor asynchronous: uniform ages, thinned by y), as a b_shape floor;
##  3. for each of the 14 units of the real-data fit np_wide_total0_nc384_t
##     (no_pool, n_c 384, b_shape estimated under lognormal(2, 0.5)), the t = 0
##     age distribution implied by the posterior means of b_shape and b_offset
##     at its trial's cycle_length: modal age, sd, the share of mass at ages
##     the model itself treats as mostly sequestered, and the share the model
##     places in sequestered compartments at t = 0.
##
## The Beta parameterisation (package beta_starts): compartment i of n_c, at age
## (i - 0.5) / n_c cycles, gets weight Beta(b, b) density at
## x = ((i - 0.5) / n_c + b_offset) mod 1, so the modal age is
## (0.5 - b_offset) mod 1 cycles.
##
##   Rscript --vanilla _scripts/inoculum-age-check.R > _data/inoculum-age-check-<date>.txt

suppressPackageStartupMessages(library(rstan))
P1 <- 11.3869 / 467.6209; P2 <- 1; P3 <- 18.5802; P4 <- 0.2242
y_circ <- function(a) { u <- 10^(P4 * (P3 - a)); ((1 - P2) + (1 - P1) * u) / (1 + u) }
vi <- function(b) 1 / (4 * (2 * b + 1))
b_from_sd_cycles <- function(s) (1 / (4 * s^2) - 1) / 2

cat("1. CIRCULATING WINDOW (model's own sequestration curve)\n\n")
for (a in c(0, 6, 12, 15, 18.58, 21, 24, 30))
    cat(sprintf("   age %5.2f h: share still circulating %.3f\n", a, y_circ(a)))
q <- sapply(c(0.95, 0.5, 0.05), \(p) uniroot(\(a) y_circ(a) - p * (1 - P1), c(0, 60))$root)
cat(sprintf("\n   95%% / 50%% / 5%% of the age-0 circulating share left at %.1f / %.1f / %.1f h\n",
            q[1], q[2], q[3]))

cat("\n2. WIDEST POSSIBLE CIRCULATING INOCULUM\n",
    "   Donor fully asynchronous (uniform ages over the cycle), venous blood takes\n",
    "   only circulating iRBCs, so the inoculum's age density is proportional to\n",
    "   y(a). Its sd, and the b_shape with the same sd, by cycle length:\n\n", sep = "")
for (cl in c(39, 40.6, 42, 42.844, 44)) {
    a <- seq(0, cl, length.out = 20001); w <- y_circ(a); w <- w / sum(w)
    m <- sum(a * w); s <- sqrt(sum((a - m)^2 * w))
    cat(sprintf("   cycle_length %6.3f h: mean age %5.2f h, sd %4.2f h (%.4f cycles) -> b_shape %5.2f; 99%% of ages below %5.2f h\n",
                cl, m, s, s / cl, b_from_sd_cycles(s / cl), a[which(cumsum(w) >= 0.99)[1]]))
}
cat("   Any real donor is more synchronous than uniform, so b_shape should be at\n",
    "   least about this value and the modal age inside the circulating window.\n", sep = "")

cat("\n3. THE REAL-DATA FIT np_wide_total0_nc384_t, PER UNIT\n",
    "   Posterior means of b_shape and b_offset per grp_init unit, cycle_length of\n",
    "   its trial. mode = modal age at t = 0 (h); sd = sd of the Beta age\n",
    "   distribution (h); >18.6h = share of t = 0 mass at ages past the\n",
    "   half-sequestration age; seq0 = share the model puts in SEQUESTERED\n",
    "   compartments at t = 0 (1 - sum of y(age) x weight), which for an\n",
    "   inoculum of circulating iRBCs should be near 0.\n\n", sep = "")
d <- readRDS("_data/wock-data-np_wide_total0_nc384_t.rds")
lv <- attr(d, "levels")
f <- readRDS("_data/wock-fit-np_wide_total0_nc384_t.rds")
pm <- function(p) unname(colMeans(as.matrix(f, pars = p)))
bs <- pm("b_shape"); bo <- pm("b_offset"); cl <- pm("cycle_length")
unit_trial <- vapply(seq_along(lv$inoc), \(u) unique(d$grp_cl[d$grp_init == u]), integer(1))
n_c <- d$n_c
res <- do.call(rbind, lapply(seq_along(lv$inoc), \(u) {
    clu <- cl[unit_trial[u]]
    age_c <- (seq_len(n_c) - 0.5) / n_c
    x <- (age_c + bo[u]) %% 1
    w <- dbeta(x, bs[u], bs[u]); w <- w / sum(w)
    age_h <- age_c * clu
    data.frame(unit = lv$inoc[u], cycle_length = clu, b_shape = bs[u], b_offset = bo[u],
               mode_h = ((0.5 - bo[u]) %% 1) * clu,
               sd_h = sqrt(vi(bs[u])) * clu,
               past_18.6h = sum(w[age_h > P3]),
               seq0 = 1 - sum(w * y_circ(age_h)))
}))
num <- c("cycle_length", "b_shape", "b_offset", "mode_h", "sd_h", "past_18.6h", "seq0")
out <- res; out[num] <- lapply(out[num], \(v) round(v, 3))
options(width = 150); print(out, row.names = FALSE)
cat(sprintf("\n   Across the 14 units: modal age %.1f-%.1f h (median %.1f); sd %.1f-%.1f h;\n   sequestered share at t = 0 %.3f-%.3f (median %.3f).\n",
            min(res$mode_h), max(res$mode_h), median(res$mode_h), min(res$sd_h), max(res$sd_h),
            min(res$seq0), max(res$seq0), median(res$seq0)))
