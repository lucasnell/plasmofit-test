## How flat a logit-normal cycle_length prior can be made without a package
## change (2026-10-26, planning an uninformative-prior arm of the single-unit
## study). The package places logit_cl ~ normal(mean_logit_cl, sd_logit_cl)
## with cycle_length = 35 + 15 * inv_logit(logit_cl). For mean 0 and several
## sds: max / min of the implied density of cycle_length over 35.5-49.5 h
## (1 = flat), and Monte Carlo prior mean, sd, and mass within 1 h of a bound.
## Second table: prior width in hours by centre and sd_logit_cl.
##
##   Rscript --vanilla _scripts/unitsim-flatprior-shapes.R > _data/unitsim-read-flatprior-shapes-<date>.txt
cl <- seq(35.5, 49.5, by = 0.5); u <- (cl - 35) / 15; x <- qlogis(u)
for (s in c(1, 1.5, 1.6, 1.7, 1.8, 2, 3)) {
    d <- dnorm(x, 0, s) / (15 * u * (1 - u))
    set.seed(1); dr <- 35 + 15 * plogis(rnorm(1e6, 0, s))
    cat(sprintf("sd_logit_cl %.1f, mean 0: density over [35.5, 49.5] h max/min %.2f; prior mean %.2f h, sd %.2f h; P(within 1 h of a bound) %.3f\n",
                s, max(d) / min(d), mean(dr), sd(dr), mean(dr < 36 | dr > 49)))
}
cat(sprintf("uniform [35, 50]: mean 42.50, sd %.2f h, P(within 1 h of a bound) %.3f\n",
            15 / sqrt(12), 2 / 15))

## Added 2026-10-26 (weakly informative priors): the same logit sd is a
## different width in hours depending on the centre, because the logit map is
## steep near the bounds. Prior mean, median, and sd in hours (Monte Carlo,
## 1e6 draws) for centres 48 h (production), 42.844 h (the single-unit truth),
## and 42.5 h (the bounds' midpoint), by sd_logit_cl.
cat("\ncentre_h sd_logit_cl mean_h median_h sd_h P(>49h)\n")
for (ctr in c(48, 42.844, 42.5)) for (s in c(0.5, 0.75, 1, 1.25, 1.5)) {
    set.seed(1); dr <- 35 + 15 * plogis(rnorm(1e6, qlogis((ctr - 35) / 15), s))
    cat(sprintf("%8.3f %11.2f %6.2f %8.2f %4.2f %7.3f\n", ctr, s, mean(dr), median(dr), sd(dr), mean(dr > 49)))
}
