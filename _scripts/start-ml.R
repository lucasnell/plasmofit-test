## Maximum likelihood of the starting-distribution models, by optimiser from
## many starts, independent of whether MCMC chains can move between regions.
## 2026-11-20, Lucas ("do 2") after entries 58 and 59 failed the gate
## (claude/findings.md, "Entries 58 and 59").
##
## Models (wockner-fit.R entries; all no_pool, n_c 384, whole inoculum
## circulating at t = 0 for the windowed starts):
##   54  per-unit wrapped symmetric Beta over the whole cycle: b_shape in
##       [2, 250] and b_offset on the circle, per unit
##   55  one Beta(alpha, beta) on [0, 13] h, alpha, beta in [1, 3], shared
##   57  per-unit mode in [0, 1] and concentration in [2, 6] on [0, 13] h
##   58  57 on [-6, 13] h (ages below 0 = the cycle's last 6 h)
##   59  55 plus a per-unit delay in [0, 12] h
## Same parameter bounds as the Stan programs. PRIORS ARE IGNORED: this is the
## likelihood surface only. The likelihood is the Stan model's: normal on
## log10(y + 1), one sd per grp_sd group (trial x cohort), here PROFILED (set
## to its maximum-likelihood value, the root mean squared residual of the
## group). cycle_length and R are per trial and grp_sd nests in trial, so
## for 54, 57 and 58 the likelihood is a sum of 13 independent per-trial
## problems. 55 and 59 share (alpha, beta) across trials: they are profiled on
## an 8 x 8 grid of (alpha, beta), each grid point a sum of per-trial maxima,
## then (alpha, beta) is refined continuously from the best grid point.
##
## The forward map is the package's own (plasmofit:::ew_poly_series for 54,
## ew_poly_series_gen for the rest, scratch build cda47a0), with the
## Erlang-window truncation bounds computed exactly as the Stan programs'
## transformed data computes them. Check mode tests the whole likelihood
## against the fitted models' log_lik before anything runs.
##
## Starts, per trial and model (per grid point for 55 and 59):
##   54, 57, 58: 50 random starts, cycle_length stratified over [35.5, 49.5]
##       h, the rest uniform (log10_total0 on [-2.5, 2.5], R on [3, 40],
##       b_shape log-uniform on [2.5, 200], b_offset, mode on [0, 1],
##       concentration on [2, 6]); plus, from each of the fitted entry's 4
##       chains, the draw with the highest likelihood for that trial.
##   55, 59: 12 random starts per grid point, the 4 chain seeds with
##       (alpha, beta) replaced by the grid point's, and a warm start from the
##       best solution at the previous grid point.
##   Each start: Nelder-Mead (maxit 4000) then BFGS (maxit 300) on the
##   unconstrained scale. More starts can only raise a maximum, so every
##   number here is the best found, a lower bound on the true maximum.
##
## READING RULE, fixed before any result (2026-11-20). Free parameters
## including the 27 profiled sds: 54 95, 57 95, 58 95, 55 69, 59 83.
##   a. 57 or 58 against 54 (equal parameter counts): a difference of a few
##      units either way means the windowed start does not cost fit and the
##      sampling failures come from the parameterisation; a windowed model
##      far below 54 means the window itself costs fit.
##   b. 59 against 55 (nested, +14 parameters): a gain below 14
##      log-likelihood units (AIC's break-even, one unit per parameter) means
##      the delay does not earn its parameters.
##   c. 58 against 57 (equal counts): whether reaching below age 0 helps.
##   d. Per trial: whether optima with cycle_length at the 50 h (or 35 h)
##      bound are genuine local maxima, and how far below the best.
## These are maxima, not posteriors and not predictive comparisons; a
## difference between maxima says how well each model CAN fit, not which
## predicts better.
##
## MODES (SML_MODE):
##   check  R likelihood (with the draw's own sds) against the Stan fits'
##          summed log_lik at 3 draws per model
##   run    one task (SLURM_ARRAY_TASK_ID 1-65: model x trial), output
##          _data/start-ml/task-<model>-<trial>.rds
##   merge  combine, refine (alpha, beta) for 55 and 59, compare with the
##          MCMC draws, print the tables
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   SML_MODE=check Rscript --vanilla _scripts/start-ml.R > _data/start-ml-check-<date>.txt
##   sbatch _scripts/start-ml.sh
##   SML_MODE=merge Rscript --vanilla _scripts/start-ml.R > _data/start-ml-<date>.txt

.libPaths(c("/home2/lan68/plasmofit/.Rlib-dev", "/home/lan68/R/x86_64-pc-linux-gnu-library/4.6"))
suppressPackageStartupMessages({library(plasmofit); library(rstan)})
stopifnot(exists("ew_poly_series_gen", envir = asNamespace("plasmofit")))
MODE <- Sys.getenv("SML_MODE", "check")
OUTDIR <- "_data/start-ml/"
MODELS <- c("54", "55", "57", "58", "59")
CFG <- c(`54` = "np_wide_total0_nc384_t", `55` = "np_start13_nc384_t",
         `57` = "np_start13unit_nc384_t", `58` = "np_start13below6unit_nc384_t",
         `59` = "np_start13delay12_nc384_t")
WINDOW <- 13
BELOW <- c(`54` = NA, `55` = 0, `57` = 0, `58` = 6, `59` = 0)
DMAX <- c(`54` = 0, `55` = 0, `57` = 0, `58` = 0, `59` = 12)
GRID_A <- c(1, 1.1, 1.25, 1.5, 1.75, 2, 2.5, 3); GRID_B <- GRID_A
N_RAND <- 50; N_RAND_G <- 12
## SML_SMOKE=<dir>: a quick end-to-end test of run mode, 2 random starts and
## a 2 x 1 grid, written to <dir> instead of _data/start-ml/
if (nzchar(Sys.getenv("SML_SMOKE"))) {
    OUTDIR <- paste0(Sys.getenv("SML_SMOKE"), "/"); N_RAND <- 2; N_RAND_G <- 2
    GRID_A <- c(1, 1.1); GRID_B <- 2
}

d <- readRDS("_data/wock-data-np_start13unit_nc384_t.rds")
lv <- attr(d, "levels")
n_c <- 384L; MU <- d$mu
stopifnot(d$n_c == n_c, d$min_cl == 35, d$max_cl == 50, d$max_R == 50, d$max_shape == 250)
M <- as.integer(ceiling(max(d$ts) / d$min_cl) + 4); R_MAX <- as.integer(M * n_c + n_c)
LG <- lgamma(seq_len(R_MAX)); RL1 <- as.numeric(seq_len(R_MAX) - 1)
bounds <- function(ts, dmax) {
    tmin <- n_c * (ts - dmax) / d$max_cl; tmax <- n_c * ts / d$min_cl
    rh <- ceiling(tmax + 10 * sqrt(tmax)); rl <- pmax(1, floor(tmin - 10 * sqrt(tmin)))
    list(rl = as.integer(rl), rh = as.integer(rh), ml = as.integer(pmax(1, ceiling(rl / n_c - 1))),
         mh = as.integer(floor((rh - 2) / n_c + 1)))
}

## ---- data, per trial ---------------------------------------------------------
s_end <- cumsum(d$n_obs); s_start <- s_end - d$n_obs + 1
obs_series <- rep(seq_len(d$n_ts), d$n_obs)
y_trans <- log10(d$y + 1)
trial_data <- function(c, dmax) {
    k <- which(d$grp_cl[obs_series] == c)
    g_obs <- d$grp_init[obs_series[k]]
    units <- sort(unique(g_obs))
    U <- lapply(units, \(g) {
        ut <- sort(unique(d$ts[k[g_obs == g]]))
        c(list(g = g, ts = ut), bounds(ut, dmax))
    })
    list(c = c, k = k, units = units, U = U,
         uidx = match(g_obs, units),
         tidx = mapply(\(u, t) match(t, U[[u]]$ts), match(g_obs, units), d$ts[k]),
         y = y_trans[k], sdg = d$grp_sd[obs_series[k]])
}

## yhat for one unit; st = list(type, ...)
unit_yhat <- function(u, cl, R, lt0, st) {
    if (st$type == "wrapped") {
        plasmofit:::ew_poly_series(cl, n_c, st$bs, st$bo, lt0, R, MU, length(u$ts), u$ts,
                                   M, R_MAX, LG, RL1, u$rl, u$rh, u$ml, u$mh)
    } else {
        plasmofit:::ew_poly_series_gen(cl, n_c, st$below, WINDOW, st$a, st$b, st$delay, lt0, R, MU,
                                       length(u$ts), u$ts, M, R_MAX, LG, RL1, u$rl, u$rh, u$ml, u$mh)
    }
}
## log-likelihood of one trial: sd profiled (sd = NULL) or given per sd group
trial_ll <- function(td, cl, R, lt0, starts, sd = NULL) {
    yh <- numeric(length(td$k))
    for (j in seq_along(td$U)) {
        v <- unit_yhat(td$U[[j]], cl, R, lt0[j], starts[[j]])
        if (any(!is.finite(v)) || any(v < 0)) return(-Inf)
        sel <- td$uidx == j
        yh[sel] <- v[td$tidx[sel]]
    }
    r <- td$y - log10(yh + 1)
    if (is.null(sd)) {
        sum(vapply(split(r, td$sdg), \(x) { s2 <- mean(x^2); -length(x) / 2 * (log(2 * pi * s2) + 1) }, 0))
    } else sum(dnorm(r, 0, sd[td$sdg], log = TRUE))
}

## ---- parameterisation -----------------------------------------------------------
lgt <- function(p) log(p / (1 - p)); ilg <- function(x) 1 / (1 + exp(-pmin(pmax(x, -25), 25)))
## theta -> (cl, R, lt0, starts) for model m and nu units; ab = shared (a, b)
unpack <- function(th, m, nu, ab = NULL) {
    cl <- 35 + 15 * ilg(th[1]); R <- 50 * ilg(th[2])
    per <- matrix(th[-(1:2)], nrow = nu, byrow = TRUE)
    lt0 <- per[, 1]
    st <- lapply(seq_len(nu), \(j) switch(m,
        `54` = list(type = "wrapped", bs = 2 + 248 * ilg(per[j, 2]), bo = per[j, 3] %% 1),
        `57` = , `58` = { mo <- ilg(per[j, 2]); k <- 2 + 4 * ilg(per[j, 3])
                          list(type = "gen", below = BELOW[[m]], a = 1 + mo * (k - 2),
                               b = 1 + (1 - mo) * (k - 2), delay = 0) },
        `55` = list(type = "gen", below = 0, a = ab[1], b = ab[2], delay = 0),
        `59` = list(type = "gen", below = 0, a = ab[1], b = ab[2], delay = 12 * ilg(per[j, 2]))))
    list(cl = cl, R = R, lt0 = lt0, st = st)
}
n_per <- c(`54` = 3, `55` = 1, `57` = 3, `58` = 3, `59` = 2)
## natural parameters -> theta (for seeds from MCMC draws)
pack <- function(m, cl, R, lt0, extra) {   # extra: nu x (n_per - 1) matrix, natural scale
    per <- cbind(lt0, switch(m,
        `54` = cbind(lgt((extra[, 1] - 2) / 248), extra[, 2]),
        `57` = , `58` = cbind(lgt(pmin(pmax(extra[, 1], 1e-9), 1 - 1e-9)),
                              lgt(pmin(pmax((extra[, 2] - 2) / 4, 1e-9), 1 - 1e-9))),
        `55` = NULL,
        `59` = cbind(lgt(pmin(pmax(extra[, 1] / 12, 1e-9), 1 - 1e-9)))))
    c(lgt((cl - 35) / 15), lgt(R / 50), as.vector(t(per)))
}
## from a draw of a fitted entry, the natural parameters of trial c
draw_nat <- function(m, dr, i, td) {
    g <- td$units; c <- td$c
    extra <- switch(m,
        `54` = cbind(dr$b_shape[i, g], dr$b_offset[i, g]),
        `57` = , `58` = cbind(dr$start_mode[i, g], dr$start_conc[i, g]),
        `55` = matrix(numeric(0), length(g), 0),
        `59` = cbind(dr$delay[i, g]))
    list(cl = dr$cycle_length[i, c], R = dr$R[i, c], lt0 = dr$log10_total0[i, g], extra = extra,
         ab = if (m %in% c("55", "59")) c(dr$start_alpha[i], dr$start_beta[i]) else NULL)
}
draw_starts <- function(m, nat) {
    nu <- length(nat$lt0)
    lapply(seq_len(nu), \(j) switch(m,
        `54` = list(type = "wrapped", bs = nat$extra[j, 1], bo = nat$extra[j, 2]),
        `57` = , `58` = { mo <- nat$extra[j, 1]; k <- nat$extra[j, 2]
                          list(type = "gen", below = BELOW[[m]], a = 1 + mo * (k - 2),
                               b = 1 + (1 - mo) * (k - 2), delay = 0) },
        `55` = list(type = "gen", below = 0, a = nat$ab[1], b = nat$ab[2], delay = 0),
        `59` = list(type = "gen", below = 0, a = nat$ab[1], b = nat$ab[2], delay = nat$extra[j, 1])))
}
load_draws <- function(m) {
    f <- readRDS(sprintf("_data/wock-fit-%s.rds", CFG[[m]]))
    pars <- intersect(c("cycle_length", "R", "log10_total0", "b_shape", "b_offset", "start_mode",
                        "start_conc", "start_alpha", "start_beta", "delay", "sd_iRBC"), f@model_pars)
    dr <- rstan::extract(f, pars = pars, permuted = FALSE, inc_warmup = FALSE)  # iter x chain x par
    out <- lapply(pars, \(p) {
        cols <- grep(sprintf("^%s(\\[|$)", p), dimnames(dr)$parameters)
        matrix(dr[, , cols], nrow = dim(dr)[1] * dim(dr)[2])      # draw (chain-major blocks) x index
    })
    names(out) <- pars
    out$chain <- rep(seq_len(dim(dr)[2]), each = dim(dr)[1])
    out$log_lik <- apply(rstan::extract(f, "log_lik", permuted = FALSE), c(1, 2), sum) |> as.vector()
    out
}

## ---- check ----------------------------------------------------------------------
if (MODE == "check") {
    cat("CHECK. The R likelihood, with each draw's own sd_iRBC, summed over 13 trials,\n",
        "against the Stan fit's summed log_lik for the same draw (1130 observations).\n",
        "Cells: the two sums and their difference, 3 draws per entry.\n\n", sep = "")
    for (m in MODELS) {
        dr <- load_draws(m)
        TD <- lapply(seq_len(d$n_grp_cl), \(c) trial_data(c, DMAX[[m]]))
        for (i in c(1, 1500, 3999)) {
            tot <- sum(vapply(TD, \(td) { nat <- draw_nat(m, dr, i, td)
                trial_ll(td, nat$cl, nat$R, nat$lt0, draw_starts(m, nat), sd = dr$sd_iRBC[i, ]) }, 0))
            cat(sprintf("   entry %s draw %4d: R %.6f | Stan %.6f | difference %.2e\n",
                        m, i, tot, dr$log_lik[i], tot - dr$log_lik[i]))
        }
    }
    quit(save = "no")
}

## ---- run: one (model, trial) ----------------------------------------------------
opt_one <- function(th0, fn) {
    o1 <- tryCatch(optim(th0, fn, method = "Nelder-Mead", control = list(maxit = 4000, reltol = 1e-10)),
                   error = function(e) NULL)
    if (is.null(o1) || !is.finite(o1$value)) return(list(par = th0, value = Inf))
    o2 <- tryCatch(optim(o1$par, fn, method = "BFGS", control = list(maxit = 300, reltol = 1e-12)),
                   error = function(e) NULL)
    if (!is.null(o2) && is.finite(o2$value) && o2$value < o1$value) o2 else o1
}
rand_theta <- function(m, nu, i, n) {
    cl <- 35.5 + 14 * (i - 0.5) / n
    per <- t(vapply(seq_len(nu), \(j) c(runif(1, -2.5, 2.5), switch(m,
        `54` = c(lgt((exp(runif(1, log(2.5), log(200))) - 2) / 248), runif(1)),
        `57` = , `58` = c(lgt(runif(1, 0.02, 0.98)), lgt(runif(1, 0.02, 0.98))),
        `55` = NULL,
        `59` = lgt(runif(1, 0.02, 0.98)))), numeric(n_per[[m]])))
    c(lgt((cl - 35) / 15), lgt(runif(1, 3, 40) / 50), as.vector(t(per)))
}

if (MODE == "run") {
    task <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", NA))
    stopifnot(!is.na(task), task >= 1, task <= length(MODELS) * d$n_grp_cl)
    m <- MODELS[(task - 1) %/% d$n_grp_cl + 1]; c <- (task - 1) %% d$n_grp_cl + 1
    dir.create(OUTDIR, showWarnings = FALSE)
    OUT <- sprintf("%stask-%s-%02d.rds", OUTDIR, m, c)
    if (file.exists(OUT)) stop(OUT, " exists; refusing to overwrite")
    set.seed(20261120 + task)
    td <- trial_data(c, DMAX[[m]]); nu <- length(td$units)
    cat(sprintf("task %d: entry %s, trial %s, %d units, %d observations\n", task, m, lv$trial[c], nu, length(td$k)))
    ## MCMC seeds: per chain, the draw with the highest PROFILED likelihood for this trial
    dr <- load_draws(m)
    prof <- vapply(seq_along(dr$chain), \(i) { nat <- draw_nat(m, dr, i, td)
        trial_ll(td, nat$cl, nat$R, nat$lt0, draw_starts(m, nat)) }, 0)
    seeds <- lapply(1:4, \(ch) which(dr$chain == ch)[which.max(prof[dr$chain == ch])])
    seed_ll <- vapply(seeds, \(i) prof[i], 0)
    cat("   draws' profiled log-lik, best per chain:", sprintf("%.2f", seed_ll), "\n")
    runs <- list()
    solve_at <- function(ab, starts) {
        fn <- function(th) { p <- unpack(th, m, nu, ab); v <- -trial_ll(td, p$cl, p$R, p$lt0, p$st)
                             if (is.finite(v)) v else 1e10 }
        res <- lapply(starts, \(s) { o <- opt_one(s$th, fn); p <- unpack(o$par, m, nu, ab)
            list(src = s$src, ll = -o$value, par = o$par, cl = p$cl) })
        res
    }
    seed_theta <- function(i, ab = NULL) { nat <- draw_nat(m, dr, i, td)
        pack(m, nat$cl, nat$R, nat$lt0, nat$extra) }
    if (m %in% c("54", "57", "58")) {
        st <- c(lapply(seq_len(N_RAND), \(i) list(src = "random", th = rand_theta(m, nu, i, N_RAND))),
                lapply(1:4, \(ch) list(src = sprintf("chain%d", ch), th = seed_theta(seeds[[ch]]))))
        runs[["none"]] <- solve_at(NULL, st)
    } else {
        prev <- NULL
        for (a in GRID_A) for (b in GRID_B) {
            st <- c(lapply(seq_len(N_RAND_G), \(i) list(src = "random", th = rand_theta(m, nu, i, N_RAND_G))),
                    lapply(1:4, \(ch) list(src = sprintf("chain%d", ch), th = seed_theta(seeds[[ch]]))),
                    if (!is.null(prev)) list(list(src = "warm", th = prev)))
            r <- solve_at(c(a, b), st)
            prev <- r[[which.max(vapply(r, `[[`, 0, "ll"))]]$par
            runs[[sprintf("%g_%g", a, b)]] <- r
            cat(sprintf("   (a, b) = (%g, %g): best %.3f\n", a, b, max(vapply(r, `[[`, 0, "ll"))))
        }
    }
    saveRDS(list(model = m, trial = c, units = td$units, runs = runs, seed_ll = seed_ll,
                 draw_prof = prof, draw_chain = dr$chain), OUT)
    cat("written", OUT, "\n")
    quit(save = "no")
}

## ---- merge ---------------------------------------------------------------------
if (MODE == "merge") {
    TR <- lv$trial; nT <- d$n_grp_cl
    best_of <- function(r) max(vapply(r, `[[`, 0, "ll"))
    Rz <- lapply(MODELS, \(m) lapply(seq_len(nT), \(c) readRDS(sprintf("%stask-%s-%02d.rds", OUTDIR, m, c))))
    names(Rz) <- MODELS
    npar <- c(`54` = 95, `55` = 69, `57` = 95, `58` = 95, `59` = 83)
    ML <- c(); per_trial <- matrix(NA, nT, length(MODELS), dimnames = list(TR, MODELS)); ABbest <- list()
    for (m in MODELS) {
        if (m %in% c("54", "57", "58")) {
            per_trial[, m] <- vapply(Rz[[m]], \(x) best_of(x$runs$none), 0)
            ML[m] <- sum(per_trial[, m])
        } else {
            keys <- names(Rz[[m]][[1]]$runs)
            surf <- vapply(keys, \(k) sum(vapply(Rz[[m]], \(x) best_of(x$runs[[k]]), 0)), 0)
            kb <- keys[which.max(surf)]; ab0 <- as.numeric(strsplit(kb, "_")[[1]])
            ## refine (alpha, beta): per-trial optimisation warm-started from the
            ## best solution at the best grid point and at its grid neighbours
            TD <- lapply(seq_len(nT), \(c) trial_data(c, DMAX[[m]]))
            warm <- lapply(seq_len(nT), \(c) { r <- Rz[[m]][[c]]$runs
                ks <- keys[vapply(keys, \(k) { v <- as.numeric(strsplit(k, "_")[[1]])
                    abs(match(v[1], GRID_A) - match(ab0[1], GRID_A)) <= 1 &&
                    abs(match(v[2], GRID_B) - match(ab0[2], GRID_B)) <= 1 }, TRUE)]
                lapply(ks, \(k) r[[k]][[which.max(vapply(r[[k]], `[[`, 0, "ll"))]]$par) })
            prof_ab <- function(ab) {
                if (any(ab < 1) || any(ab > 3)) return(-1e10)
                sum(vapply(seq_len(nT), \(c) { td <- TD[[c]]; nu <- length(td$units)
                    fn <- function(th) { p <- unpack(th, m, nu, ab); v <- -trial_ll(td, p$cl, p$R, p$lt0, p$st)
                                         if (is.finite(v)) v else 1e10 }
                    -min(vapply(warm[[c]], \(th) optim(th, fn, method = "BFGS",
                                 control = list(maxit = 200))$value, 0)) }, 0))
            }
            ref <- optim(ab0, \(ab) -prof_ab(ab), method = "Nelder-Mead",
                         control = list(maxit = 60, reltol = 1e-8))
            ABbest[[m]] <- list(grid = ab0, grid_ll = max(surf), ref = ref$par, ref_ll = -ref$value,
                                surf = surf)
            kbest <- kb
            per_trial[, m] <- vapply(Rz[[m]], \(x) best_of(x$runs[[kbest]]), 0)
            ML[m] <- max(max(surf), -ref$value)
        }
    }
    cat("1. MAXIMUM LOG-LIKELIHOOD per model: the best found over all starts, summed\n",
        "   over 13 trials (1130 observations), sd per trial x cohort profiled, priors\n",
        "   ignored. npar = free parameters including the 27 sds. Higher is better;\n",
        "   every value is a lower bound on the true maximum. 'MCMC best' = the best\n",
        "   single draw of the fitted entry, re-scored with profiled sds (summed over\n",
        "   trials within that draw); the maximum must not be below it.\n\n", sep = "")
    for (m in MODELS) {
        dp <- Reduce(`+`, lapply(Rz[[m]], `[[`, "draw_prof"))
        cat(sprintf("   entry %s: max %8.2f | npar %3d | MCMC best %8.2f | max - MCMC best %+6.2f\n",
                    m, ML[m], npar[[m]], max(dp), ML[m] - max(dp)))
    }
    for (m in names(ABbest)) {
        x <- ABbest[[m]]
        cat(sprintf("   entry %s (alpha, beta): best grid point (%g, %g) at %.2f; refined (%.3f, %.3f) at %.2f\n",
                    m, x$grid[1], x$grid[2], x$grid_ll, x$ref[1], x$ref[2], x$ref_ll))
    }
    cat("\n2. DIFFERENCES between maxima (log-likelihood units; positive = first model\n",
        "   fits better), with the difference in free parameters.\n\n", sep = "")
    for (pr in list(c("57", "54"), c("58", "54"), c("58", "57"), c("59", "55"), c("55", "54"), c("59", "54")))
        cat(sprintf("   %s - %s: %+7.2f (parameters %+d)\n", pr[1], pr[2], ML[pr[1]] - ML[pr[2]],
                    npar[[pr[1]]] - npar[[pr[2]]]))
    cat("\n3. PER TRIAL, maximum log-likelihood (best found; for 55 and 59 at the\n",
        "   best (alpha, beta) grid point). Units: log-likelihood; cells are per\n",
        "   trial, summing to section 1 for 54, 57 and 58.\n\n", sep = "")
    print(round(per_trial, 2))
    cat("\n4. LOCAL OPTIMA AT THE cycle_length BOUNDS, per trial and model (54, 57, 58;\n",
        "   55 and 59 at their best grid point): best log-likelihood among starts that\n",
        "   ended with cycle_length >= 49.5 h, and <= 35.5 h, minus the trial's best\n",
        "   (NA = no start ended there); and how many of the 54 starts reached within\n",
        "   0.5 units of the trial's best.\n\n", sep = "")
    for (m in MODELS) {
        key <- if (m %in% c("54", "57", "58")) "none" else
            names(Rz[[m]][[1]]$runs)[which.max(vapply(names(Rz[[m]][[1]]$runs), \(k)
                sum(vapply(Rz[[m]], \(x) best_of(x$runs[[k]]), 0)), 0))]
        tab <- t(vapply(Rz[[m]], \(x) { r <- x$runs[[key]]; ll <- vapply(r, `[[`, 0, "ll"); cl <- vapply(r, `[[`, 0, "cl")
            b <- max(ll); c(hi = if (any(cl >= 49.5)) max(ll[cl >= 49.5]) - b else NA,
                            lo = if (any(cl <= 35.5)) max(ll[cl <= 35.5]) - b else NA,
                            n_best = sum(ll > b - 0.5), n = length(ll),
                            cl_best = cl[which.max(ll)]) }, numeric(5)))
        rownames(tab) <- TR
        cat(sprintf("   entry %s\n", m)); print(round(tab, 2))
    }
}
