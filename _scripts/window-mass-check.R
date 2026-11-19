## Would a 24 h window recover what entry 54's per-unit starts use? (2026-11-19,
## Lucas: "for 1, would it still only allow circulating?") For each unit of
## np_wide_total0_nc384_t (entry 54), the share of its t = 0 age distribution
## (posterior means of b_shape and b_offset, n_c 384, its trial's
## cycle_length; package beta_starts) in four age bands: inside the 13 h ring
## window; 13-24 h, which a 24 h window would add; 24 h to cycle_length minus
## 6 h; and the last 6 h of the cycle, i.e. "just before age 0" in wrapped
## terms: schizonts bursting within 6 h of injection. Also the share the
## model's own sequestration curve says would still circulate at each age,
## so the 13-24 h band can be read against a circulating-only inoculum.
##
##   Rscript --vanilla _scripts/window-mass-check.R > _data/window-mass-check-<date>.txt
.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
suppressPackageStartupMessages({library(plasmofit); library(rstan)})
d <- readRDS("_data/wock-data-np_wide_total0_nc384_t.rds"); lv <- attr(d, "levels")
f <- readRDS("_data/wock-fit-np_wide_total0_nc384_t.rds")
pm <- function(p) unname(colMeans(as.matrix(f, pars = p)))
bs <- pm("b_shape"); bo <- pm("b_offset"); cl <- pm("cycle_length")
ut <- vapply(seq_along(lv$inoc), \(u) unique(d$grp_cl[d$grp_init == u]), integer(1))
P1 <- 11.3869 / 467.6209; P3 <- 18.5802; P4 <- 0.2242
yc <- function(a) { u <- 10^(P4 * (P3 - a)); (1 - P1) * u / (1 + u) }
cat("Share of each unit's t = 0 age distribution by band (entry 54, posterior means).\n",
    "Bands in hours since invasion; 'last 6 h' = the final 6 h of the cycle, i.e.\n",
    "schizonts bursting within 6 h of injection. Rows sum to 1.\n\n", sep = "")
nc <- 384L
out <- t(vapply(seq_along(lv$inoc), \(u) {
    clu <- cl[ut[u]]; w <- plasmofit:::beta_starts(nc, bs[u], bo[u], 1)
    a <- (seq_len(nc) - 0.5) / nc * clu
    c(cycle = clu, `0-13h` = sum(w[a < 13]), `13-24h` = sum(w[a >= 13 & a < 24]),
      `24h-to-last6` = sum(w[a >= 24 & a < clu - 6]), `last 6h` = sum(w[a >= clu - 6]),
      circ_share_13_24 = sum(w[a >= 13 & a < 24] * yc(a[a >= 13 & a < 24])) / max(sum(w[a >= 13 & a < 24]), 1e-12))
}, numeric(6)))
tab <- data.frame(unit = lv$inoc, round(out, 3), check.names = FALSE)
print(tab, row.names = FALSE)
cat(sprintf("\nMedian over 14 units: 0-13 h %.2f, 13-24 h %.2f, 24 h to last 6 h %.2f, last 6 h %.2f\n",
            median(out[, 2]), median(out[, 3]), median(out[, 4]), median(out[, 5])))
cat(sprintf("Units with more mass in the last 6 h than in 13-24 h: %d of 14\n", sum(out[, 5] > out[, 3])))
cat(sprintf("Circulating share of a cohort, by the model's curve: %.2f at 13 h, %.2f at 18.6 h, %.2f at 24 h\n",
            yc(13), yc(18.58), yc(24)))
