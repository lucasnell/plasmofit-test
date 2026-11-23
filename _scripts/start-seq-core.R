## Shared core of the sequestration-schedule profiles: data, the forward map
## with a free schedule, the per-trial likelihood with profiled sds, the
## parameterisation, and the optimiser. Extracted from _scripts/start-seq.R on
## 2026-11-23 so that _scripts/start-seq2.R (the 2-D profile over midpoint and
## slope) cannot drift from it. seq_log_y() and trial_ll() gained a `p4`
## argument defaulting to the inherited 0.2242, so every earlier call is
## unchanged; start-seq.R's check and merge output is byte-identical before and
## after the extraction (verified 2026-11-23).
##
## Sourced, never run directly. The caller defines MODE, OUTDIR, the grids and
## the number of random starts.

.libPaths(c("/home2/lan68/plasmofit/.Rlib-dev", "/home/lan68/R/x86_64-pc-linux-gnu-library/4.6"))
suppressPackageStartupMessages({library(plasmofit); library(rstan)})
stopifnot(exists("ew_poly_series_gen", envir = asNamespace("plasmofit")))
WINDOW <- 13
DMAX <- c(`54` = 0, `57` = 0, `59` = 12)
P1 <- 11.3869 / 467.6209; P3_DEF <- 18.5802; P4 <- 0.2242
MODELS <- c("54", "57", "59")
CFG <- c(`54` = "np_wide_total0_nc384_t", `57` = "np_start13unit_nc384_t",
         `59` = "np_start13delay12_nc384_t")
MLDIR <- "_data/start-ml/"

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

## ---- the forward map with a free sequestration midpoint -----------------------
## log y(age) = log(1 - p1) - log1p(exp(-p4 (p3 - age) log 10)); compartment k
## is at age k * cycle_length / n_c, as make_log_y_vals() has it.
seq_log_y <- function(cl, p3, p4 = P4) {
    lu <- p4 * (p3 - seq_len(n_c) * (cl / n_c)) * log(10)
    log1p(-P1) - log1p(exp(-lu))
}
## D/C from a starting mass vector W and a sequestration curve.
## type "wrap": N_j(0) = y_j W_j, the composition of a population that aged
## into stage j in situ (entry 54). type "win": N_j(0) = W_j, the whole
## inoculum circulating (entries 55-59).
make_DC_free <- function(log_y, W, type) {
    c_k <- c(1, exp(log_y[-n_c] - log_y[1]))
    Nc <- if (type == "wrap") exp(log_y[1]) * (1 - c(0, -expm1(diff(log_y)))) * W else W / c_k
    DC <- numeric(3 * n_c - 1)
    for (dd in 0:(n_c - 1)) {
        DC[dd + 1]       <- sum(Nc[1:(n_c - dd)] * c_k[(dd + 1):n_c])
        DC[2 * n_c + dd] <- sum(W[1:(n_c - dd)]  * c_k[(dd + 1):n_c])
    }
    for (dd in 1:(n_c - 1)) DC[2 * n_c - dd] <- sum(W[(dd + 1):n_c] * c_k[1:(n_c - dd)])
    DC
}

## ---- data, per trial ---------------------------------------------------------
s_end <- cumsum(d$n_obs); s_start <- s_end - d$n_obs + 1
obs_series <- rep(seq_len(d$n_ts), d$n_obs)
y_trans <- log10(d$y + 1)
trial_data <- function(c, dmax) {
    k <- which(d$grp_cl[obs_series] == c)
    g_obs <- d$grp_init[obs_series[k]]
    units <- sort(unique(g_obs))
    U <- lapply(units, \(g) { ut <- sort(unique(d$ts[k[g_obs == g]]))
                              c(list(g = g, ts = ut), bounds(ut, dmax)) })
    list(c = c, k = k, units = units, U = U, uidx = match(g_obs, units),
         tidx = mapply(\(u, t) match(t, U[[u]]$ts), match(g_obs, units), d$ts[k]),
         y = y_trans[k], sdg = d$grp_sd[obs_series[k]])
}

unit_yhat <- function(u, cl, R, lt0, st, log_y) {
    total0 <- 10^lt0
    W <- if (st$type == "wrap") plasmofit:::beta_starts(n_c, st$bs, st$bo, total0)
         else plasmofit:::windowed_beta_starts_below(n_c, 0, WINDOW / cl, st$a, st$b, total0)
    tv <- if (st$type == "wrap") u$ts else u$ts - st$delay
    plasmofit:::ew_poly_eval(make_DC_free(log_y, W, if (st$type == "wrap") "wrap" else "win"),
                             n_c, cl, R, MU, length(tv), tv, M, R_MAX, LG, RL1,
                             u$rl, u$rh, u$ml, u$mh)
}
trial_ll <- function(td, cl, R, lt0, starts, p3, sd = NULL, p4 = P4) {
    log_y <- seq_log_y(cl, p3, p4)
    yh <- numeric(length(td$k))
    for (j in seq_along(td$U)) {
        v <- unit_yhat(td$U[[j]], cl, R, lt0[j], starts[[j]], log_y)
        if (any(!is.finite(v)) || any(v < 0)) return(-Inf)
        sel <- td$uidx == j
        yh[sel] <- v[td$tidx[sel]]
    }
    r <- td$y - log10(yh + 1)
    if (is.null(sd)) {
        sum(vapply(split(r, td$sdg), \(x) { s2 <- mean(x^2); -length(x) / 2 * (log(2 * pi * s2) + 1) }, 0))
    } else sum(dnorm(r, 0, sd[td$sdg], log = TRUE))
}

## ---- parameterisation (copied from start-ml.R) ----------------------------------
lgt <- function(p) log(p / (1 - p)); ilg <- function(x) 1 / (1 + exp(-pmin(pmax(x, -25), 25)))
unpack <- function(th, m, nu, ab = NULL) {
    cl <- 35 + 15 * ilg(th[1]); R <- 50 * ilg(th[2])
    per <- matrix(th[-(1:2)], nrow = nu, byrow = TRUE)
    lt0 <- per[, 1]
    st <- lapply(seq_len(nu), \(j) switch(m,
        `54` = list(type = "wrap", bs = 2 + 248 * ilg(per[j, 2]), bo = per[j, 3] %% 1),
        `57` = { mo <- ilg(per[j, 2]); k <- 2 + 4 * ilg(per[j, 3])
                 list(type = "gen", a = 1 + mo * (k - 2), b = 1 + (1 - mo) * (k - 2), delay = 0) },
        `59` = list(type = "gen", a = ab[1], b = ab[2], delay = 12 * ilg(per[j, 2]))))
    list(cl = cl, R = R, lt0 = lt0, st = st)
}
n_per <- c(`54` = 3, `57` = 3, `59` = 2)
pack <- function(m, cl, R, lt0, extra) {
    per <- cbind(lt0, switch(m,
        `54` = cbind(lgt((extra[, 1] - 2) / 248), extra[, 2]),
        `57` = cbind(lgt(pmin(pmax(extra[, 1], 1e-9), 1 - 1e-9)),
                     lgt(pmin(pmax((extra[, 2] - 2) / 4, 1e-9), 1 - 1e-9))),
        `59` = cbind(lgt(pmin(pmax(extra[, 1] / 12, 1e-9), 1 - 1e-9)))))
    c(lgt((cl - 35) / 15), lgt(R / 50), as.vector(t(per)))
}
draw_nat <- function(m, dr, i, td) {
    g <- td$units; c <- td$c
    extra <- switch(m,
        `54` = cbind(dr$b_shape[i, g], dr$b_offset[i, g]),
        `57` = cbind(dr$start_mode[i, g], dr$start_conc[i, g]),
        `59` = cbind(dr$delay[i, g]))
    list(cl = dr$cycle_length[i, c], R = dr$R[i, c], lt0 = dr$log10_total0[i, g], extra = extra,
         ab = if (m == "59") c(dr$start_alpha[i], dr$start_beta[i]) else NULL)
}
load_draws <- function(m) {
    f <- readRDS(sprintf("_data/wock-fit-%s.rds", CFG[[m]]))
    pars <- intersect(c("cycle_length", "R", "log10_total0", "b_shape", "b_offset", "start_mode",
                        "start_conc", "start_alpha", "start_beta", "delay", "sd_iRBC"), f@model_pars)
    dr <- rstan::extract(f, pars = pars, permuted = FALSE, inc_warmup = FALSE)
    out <- lapply(pars, \(p) {
        cols <- grep(sprintf("^%s(\\[|$)", p), dimnames(dr)$parameters)
        matrix(dr[, , cols], nrow = dim(dr)[1] * dim(dr)[2])
    })
    names(out) <- pars
    out$chain <- rep(seq_len(dim(dr)[2]), each = dim(dr)[1])
    out
}
## start-ml.R's best theta for this (model, trial), as a warm start
ml_best <- function(m, c, dir = MLDIR) {
    f <- sprintf("%stask-%s-%02d.rds", dir, m, c)
    if (!file.exists(f)) return(NULL)
    x <- readRDS(f)
    best <- NULL; bll <- -Inf
    for (k in names(x$runs)) for (r in x$runs[[k]])
        if (is.finite(r$ll) && r$ll > bll) { bll <- r$ll; best <- r$par }
    list(th = best, ll = bll)
}


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
        `57` = c(lgt(runif(1, 0.02, 0.98)), lgt(runif(1, 0.02, 0.98))),
        `59` = lgt(runif(1, 0.02, 0.98)))), numeric(n_per[[m]])))
    c(lgt((cl - 35) / 15), lgt(runif(1, 3, 40) / 50), as.vector(t(per)))
}

