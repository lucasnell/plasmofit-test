## Write entry 59's MAXIMUM-LIKELIHOOD point as a Stan init list, for
## wockner-fit.R entry 61. 2026-11-23, Lucas ("1-3").
##
## WHY. Entry 59 (shared Beta on [0, 13] h plus a per-unit delay) failed the
## gate with its chains ~87 log-likelihood units below its own maximum
## (claude/findings.md, "Entries 58 and 59"). Starting chains AT that maximum
## asks whether a posterior exists there at all: if they stay and mix, the
## failure was the search; if they leave for the cycle_length bounds again,
## the region is not a posterior mode and the model is at fault.
##
## THIS IS A DIAGNOSTIC, NOT A POSTERIOR. Chains started in one place can
## agree with each other without the chain having explored anything, so a
## passing R-hat here does NOT license reporting entry 59's posterior. What
## it licenses is a statement about WHY the ordinary fit failed.
##
## SOURCE. _data/start-ml/task-59-*.rds (SLURM 35292), at the (alpha, beta)
## grid point that maximises the SUM over the 13 trials, since entry 59's
## start shape is shared. That is -766.28 summed over 1130 observations
## (_data/start-ml-2026-11-20.txt), at the model's own fixed sequestration
## midpoint, which is what the Stan program uses.
##
## The observation sds are not free in start-ml.R: they are profiled, i.e.
## set to the root mean squared residual of each trial x cohort group. Those
## values are inverted here to z_sd_iRBC = log(sd / 0.2), the parameter the
## Stan program actually declares.
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   Rscript --vanilla _scripts/start-ml-init.R > _data/start-ml-init-59-<date>.txt

.libPaths(c("/home2/lan68/plasmofit/.Rlib-dev", "/home/lan68/R/x86_64-pc-linux-gnu-library/4.6"))
suppressPackageStartupMessages(library(plasmofit))
OUT <- "_data/start-ml-init-59.rds"
MLDIR <- "_data/start-ml/"
WINDOW <- 13; DMAX <- 12
d <- readRDS("_data/wock-data-np_start13delay12_nc384_t.rds")
lv <- attr(d, "levels")
stopifnot(d$start_dist == 1L, d$start_delay == 1L, d$delay_max == DMAX,
          d$start_window == WINDOW, d$center_R == 1L, d$center_cl == 1L,
          d$b_shape_fixed == 0L, d$total0_anchor == 0L, d$sd_iRBC_fixed == 0L)
n_c <- d$n_c; nT <- d$n_grp_cl; nU <- d$n_grp_init
M <- as.integer(ceiling(max(d$ts) / d$min_cl) + 4); R_MAX <- as.integer(M * n_c + n_c)
LG <- lgamma(seq_len(R_MAX)); RL1 <- as.numeric(seq_len(R_MAX) - 1)
ilg <- function(x) 1 / (1 + exp(-pmin(pmax(x, -25), 25)))

obs_series <- rep(seq_len(d$n_ts), d$n_obs)
y_trans <- log10(d$y + 1)
obs_cl <- d$grp_cl[obs_series]; obs_init <- d$grp_init[obs_series]; obs_sd <- d$grp_sd[obs_series]

## ---- the shared (alpha, beta) that maximises the summed fit -------------------
per <- lapply(seq_len(nT), \(c) {
    x <- readRDS(sprintf("%stask-59-%02d.rds", MLDIR, c))
    list(units = x$units, runs = x$runs)
})
subs <- Reduce(intersect, lapply(per, \(p) names(p$runs)))
tot <- vapply(subs, \(sb) sum(vapply(per, \(p) max(vapply(p$runs[[sb]], `[[`, 0, "ll")), 0)), 0)
kb <- subs[which.max(tot)]
ab <- as.numeric(strsplit(kb, "_")[[1]])
cat(sprintf("best shared start over %d (alpha, beta) grid points: (%g, %g), summed log-lik %.2f\n",
            length(subs), ab[1], ab[2], max(tot)))

## ---- per-trial and per-unit parameters at that point --------------------------
cl <- numeric(nT); R <- numeric(nT); lt0 <- numeric(nU); dly <- numeric(nU)
for (c in seq_len(nT)) {
    r <- per[[c]]$runs[[kb]]
    th <- r[[which.max(vapply(r, `[[`, 0, "ll"))]]$par
    u <- per[[c]]$units; nu <- length(u)
    cl[c] <- 35 + 15 * ilg(th[1]); R[c] <- 50 * ilg(th[2])
    pm <- matrix(th[-(1:2)], nrow = nu, byrow = TRUE)
    lt0[u] <- pm[, 1]; dly[u] <- DMAX * ilg(pm[, 2])
}
cat(sprintf("cycle_length %.2f-%.2f h; R %.2f-%.2f; log10_total0 %.2f-%.2f; delay %.2f-%.2f h\n",
            min(cl), max(cl), min(R), max(R), min(lt0), max(lt0), min(dly), max(dly)))

## ---- residual sds, which start-ml.R profiled rather than estimated ------------
yh <- numeric(d$n_total_obs)
for (g in seq_len(nU)) {
    k <- which(obs_init == g)
    c <- unique(obs_cl[k]); stopifnot(length(c) == 1)
    ut <- sort(unique(d$ts[k]))
    tmin <- n_c * (ut - DMAX) / d$max_cl; tmax <- n_c * ut / d$min_cl
    rl <- as.integer(pmax(1, floor(tmin - 10 * sqrt(tmin)))); rh <- as.integer(ceiling(tmax + 10 * sqrt(tmax)))
    v <- plasmofit:::ew_poly_series_gen(cl[c], n_c, 0, WINDOW, ab[1], ab[2], dly[g], lt0[g],
                                        R[c], d$mu, length(ut), ut, M, R_MAX, LG, RL1, rl, rh,
                                        as.integer(pmax(1, ceiling(rl / n_c - 1))),
                                        as.integer(floor((rh - 2) / n_c + 1)))
    yh[k] <- v[match(d$ts[k], ut)]
}
resid <- y_trans - log10(yh + 1)
sd_g <- vapply(split(resid, obs_sd), \(x) sqrt(mean(x^2)), 0)
ll <- sum(dnorm(resid, 0, sd_g[obs_sd], log = TRUE))
cat(sprintf("summed log-likelihood rebuilt here: %.4f (start-ml: %.4f, difference %.2e)\n",
            ll, max(tot), ll - max(tot)))
stopifnot(abs(ll - max(tot)) < 1e-6)
cat(sprintf("sd_iRBC over %d trial x cohort groups: %.3f-%.3f (median %.3f)\n",
            length(sd_g), min(sd_g), max(sd_g), median(sd_g)))

## ---- as the Stan program's parameters -----------------------------------------
## alpha is at its lower bound 1 on the grid; nudged inside so the
## unconstrained value is finite.
logit_R <- log(R / (d$max_R - R))
logit_cl <- log((cl - d$min_cl) / (d$max_cl - cl))
init <- list(
    b_shape_free = numeric(0),
    b_off_vec    = array(numeric(0), c(0, 2)),
    start_alpha  = array(min(max(ab[1], 1.001), d$start_shape_cap - 0.001), 1),
    start_beta   = array(min(max(ab[2], 1.001), d$start_shape_cap - 0.001), 1),
    start_mode   = numeric(0),
    start_conc   = numeric(0),
    delay        = pmin(pmax(dly, 0.01), DMAX - 0.01),
    log10_total0 = lt0,
    delta_total0 = numeric(0),
    sigma_total0 = numeric(0),
    mu_logit_R   = mean(logit_R),
    sigma_logit_R = max(sd(logit_R), 0.05),
    eta_R        = logit_R,                 # center_R = 1: eta is logit_R itself
    mu_logit_cl  = mean(logit_cl),
    sigma_logit_cl = max(sd(logit_cl), 0.05),
    eta_cl       = logit_cl,                # center_cl = 1
    z_sd_iRBC    = log(sd_g / 0.2))
saveRDS(list(init = init, ll = ll, alpha_beta = ab, cycle_length = cl, R = R,
             log10_total0 = lt0, delay = dly, sd_iRBC = unname(sd_g)), OUT)
cat("written", OUT, "\n")
cat(sprintf("  sigma_logit_cl %.3f, sigma_logit_R %.3f; mu_logit_cl %.3f (%.2f h)\n",
            init$sigma_logit_cl, init$sigma_logit_R, init$mu_logit_cl,
            d$min_cl + (d$max_cl - d$min_cl) * ilg(init$mu_logit_cl)))

## ---- does Stan accept it, and agree about the likelihood? ----------------------
h <- suppressWarnings(rstan::sampling(plasmofit:::stanmodels$archer_fit, data = d, chains = 1,
                                      iter = 1, refresh = 0, seed = 1,
                                      algorithm = "Fixed_param", init = list(init)))
u <- rstan::unconstrain_pars(h, init)
cat(sprintf("unconstrain_pars round trip: %d coordinates (model has %d), all finite: %s\n",
            length(u), rstan::get_num_upars(h), all(is.finite(u))))
stopifnot(length(u) == rstan::get_num_upars(h), all(is.finite(u)))
cp <- rstan::constrain_pars(h, u)
cat(sprintf("  re-constrained: max |difference| in delay %.2e, log10_total0 %.2e, start_alpha %.2e\n",
            max(abs(cp$delay - init$delay)), max(abs(cp$log10_total0 - init$log10_total0)),
            abs(cp$start_alpha - init$start_alpha)))
hl <- suppressWarnings(rstan::sampling(plasmofit:::stanmodels$archer_fit, data = d, chains = 1,
                                       iter = 1, refresh = 0, seed = 1,
                                       algorithm = "Fixed_param", init = list(init)))
sll <- sum(as.matrix(hl, pars = "log_lik")[1, ])
cat(sprintf("Stan's summed log_lik at this init: %.4f (here %.4f, difference %.2e)\n",
            sll, ll, sll - ll))
