## Validate the windowed, shared starting distribution (start_dist = 1)
## before anything samples with it. 2026-11-18, for
## claude/start-dist-design.md. The build lives in a scratch clone
## (/home2/lan68/plasmofit/plasmofit-dev, branch start-dist) installed into the
## scratch library /home2/lan68/plasmofit/.Rlib-dev; the live library still
## holds the old build, untouched.
##
## No sampling: every check evaluates the compiled model at fixed points
## (Fixed_param handles, log_prob, grad_log_prob), as unitsim.R's check mode
## does.
##
## MODES (START_VALIDATE):
##   old   in the LIVE library: log_prob and grad_log_prob of the production
##         data list (np_wide_total0_nc384_t, no_pool, n_c 384) at K fixed
##         unconstrained points, written to _data/start-dist-validate-old.rds.
##   new   in the SCRATCH library, the rest:
##     1. windowed_beta_starts() against exact pbeta differences, and its
##        continuity as cycle_length moves the window across a compartment
##        edge;
##     2. start_dist = 0 against the old build at the same K points: same
##        model, so log_prob and its gradient should agree to rounding;
##     3. start_dist = 1 on the production data, 13 h window, cap 3: the
##        parameter layout, finite log_prob and gradient, and the gradient
##        against central finite differences;
##     4. the matrix_exp cross-check (run_check) at n_c 96 under both
##        start_dist values: max_rel_diff compares the fast Erlang-window
##        likelihood path with an independent matrix exponential of the
##        same starting state.
##
##   ref1, ref23  references for start_dist 1 and for 2-3 from earlier
##         builds; see their blocks below.
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   START_VALIDATE=old Rscript --vanilla _scripts/start-dist-validate.R
##   START_VALIDATE=new Rscript --vanilla _scripts/start-dist-validate.R > _data/start-dist-validate-<date>.txt

MODE <- Sys.getenv("START_VALIDATE", "new")
LIVE <- "/home/lan68/R/x86_64-pc-linux-gnu-library/4.6"
DEV <- "/home2/lan68/plasmofit/.Rlib-dev"
## START_VALIDATE_LIB=live runs the `new` checks against the LIVE library,
## to verify the installed binary once the build is installed there; section
## 2 then compares it with the saved reference from the old build.
NEWLIB <- if (Sys.getenv("START_VALIDATE_LIB", "dev") == "live") LIVE else DEV
.libPaths(if (MODE %in% c("old", "ref1")) LIVE else NEWLIB)
suppressPackageStartupMessages({library(plasmofit); library(rstan)})
rstan_options(threads_per_chain = 1L)
REF <- "_data/start-dist-validate-old.rds"
K <- 8L
WINDOW <- 13; CAP <- 3

d0 <- readRDS("_data/wock-data-np_wide_total0_nc384_t.rds")
handle <- function(d, model = "archer_fit", init = "random")
    suppressWarnings(rstan::sampling(plasmofit:::stanmodels[[model]], data = d,
                                     chains = 1, iter = 1, refresh = 0, seed = 1,
                                     algorithm = "Fixed_param", init = init))

## ref1 (added 2026-11-19, before start_dist 2 and 3): in the LIVE library,
## which then held the start_dist 0/1 build (13db2bd), log_prob and gradient of
## start_dist = 1 (13 h, cap 3) at K fixed points, so the generalised build
## can be checked to leave start_dist = 1 unchanged.
REF1 <- "_data/start-dist-validate-ref1.rds"
if (MODE == "ref1") {
    cat("library:", find.package("plasmofit"), "\n")
    d <- d0; d$start_dist <- 1L; d$start_window <- WINDOW; d$start_shape_cap <- CAP
    f <- handle(d)
    np <- rstan::get_num_upars(f)
    set.seed(20261119)
    U <- matrix(rnorm(K * np, 0, 0.3), K, np)
    lp <- apply(U, 1, \(u) rstan::log_prob(f, u))
    gr <- t(apply(U, 1, \(u) rstan::grad_log_prob(f, u)))
    saveRDS(list(U = U, lp = lp, gr = gr, np = np), REF1)
    cat(sprintf("start_dist 1 reference: %d parameters, %d points, written to %s\n", np, K, REF1))
    quit(save = "no")
}

if (MODE == "old") {
    cat("library:", find.package("plasmofit"), "\n")
    f <- handle(d0)
    np <- rstan::get_num_upars(f)
    set.seed(20261118)
    U <- matrix(rnorm(K * np, 0, 0.3), K, np)
    lp <- apply(U, 1, \(u) rstan::log_prob(f, u))
    gr <- t(apply(U, 1, \(u) rstan::grad_log_prob(f, u)))
    saveRDS(list(U = U, lp = lp, gr = gr, np = np), REF)
    cat(sprintf("old build: %d unconstrained parameters; log_prob at %d points written to %s\n",
                np, K, REF))
    quit(save = "no")
}

cat("library:", find.package("plasmofit"), "\n")
stopifnot(startsWith(find.package("plasmofit"), NEWLIB),
          exists("windowed_beta_starts", envir = asNamespace("plasmofit")))
## dist >= 1, not dist == 1: until 2026-11-20 this read `dist == 1L`, so
## section 6's start_dist 2 and 3 checks ran with a 1 h window, not 13 h.
## The fitted entries 56-57 were unaffected (their data come from
## archer_stan_data() with start_window 13, checked at the end of section 6).
withS <- function(d, dist, window = WINDOW, cap = CAP) {
    d$start_dist <- as.integer(dist); d$start_window <- if (dist >= 1L) window else 1
    d$start_shape_cap <- cap
    ## the extended build (2026-11-20) also requires these three; the off
    ## values. Earlier builds ignore them (rstan reads only declared data).
    if (is.null(d$start_below)) { d$start_below <- 0; d$start_delay <- 0L; d$delay_max <- 0 }
    d
}

## ref23 (added 2026-11-20, before start_below and start_delay): in the
## SCRATCH library while it still held the fc7db77 build, log_prob and gradient
## of start_dist 2 and 3 (13 h, cap 3) at K fixed points, so the extended
## build can be checked to leave them unchanged (section 7).
REF23 <- "_data/start-dist-validate-ref23.rds"
if (MODE == "ref23") {
    out <- list()
    for (dist in 2:3) {
        f <- handle(withS(d0, dist))
        np <- rstan::get_num_upars(f)
        set.seed(20261120 + dist)
        U <- matrix(rnorm(K * np, 0, 0.3), K, np)
        out[[as.character(dist)]] <- list(U = U, np = np,
            lp = apply(U, 1, \(u) rstan::log_prob(f, u)),
            gr = t(apply(U, 1, \(u) rstan::grad_log_prob(f, u))))
        cat(sprintf("start_dist %d reference: %d parameters, %d points\n", dist, np, K))
    }
    saveRDS(out, REF23)
    cat("written to", REF23, "\n")
    quit(save = "no")
}

## ---- 1. windowed_beta_starts against exact -----------------------------------
cat("\n1. windowed_beta_starts() against exact pbeta differences over each\n",
    "   compartment's overlap with the window. Cells: max |difference| in\n",
    "   compartment weight (weights sum to 1).\n\n", sep = "")
exact <- function(n_c, w, a, b) {
    h <- 1 / n_c; out <- numeric(n_c)
    for (i in seq_len(n_c)) { lo <- (i - 1) * h; if (lo >= w) break
        hi <- min(i * h, w); out[i] <- pbeta(hi / w, a, b) - pbeta(lo / w, a, b) }
    out / sum(out)
}
worst <- 0
for (nc in c(96L, 384L, 512L)) for (cl in c(39.4, 42.84, 45.6))
    for (ab in list(c(1, 1), c(1.2, 1.2), c(1.5, 3), c(2, 2), c(3, 1), c(3, 3))) {
        w <- WINDOW / cl
        g <- plasmofit:::windowed_beta_starts(nc, w, ab[1], ab[2], 1)
        e <- exact(nc, w, ab[1], ab[2])
        dd <- max(abs(g - e)); worst <- max(worst, dd)
        stopifnot(abs(sum(g) - 1) < 1e-12, all(g >= 0), all(g[((seq_len(nc) - 1) / nc) >= w] == 0))
        if (nc == 384L && cl == 42.84)
            cat(sprintf("   n_c %d, cycle %.2f h, a %.1f b %.1f: max |diff| %.2e\n", nc, cl, ab[1], ab[2], dd))
    }
cat(sprintf("   worst over n_c {96, 384, 512} x 3 cycle lengths x 6 shapes: %.2e\n", worst))
## continuity: move the window edge across a compartment boundary
nc <- 384L; edge_cl <- WINDOW / (34 / nc)            # w exactly 34 compartments
tot <- sapply(edge_cl * (1 + c(-1e-7, 0, 1e-7)), \(cl)
    sum(abs(diff(rbind(plasmofit:::windowed_beta_starts(nc, WINDOW / cl, 2, 3, 1))))))
v <- lapply(edge_cl * (1 + c(-1e-7, 1e-7)), \(cl) plasmofit:::windowed_beta_starts(nc, WINDOW / cl, 2, 3, 1))
cat(sprintf("   continuity across a compartment edge (cycle_length +/- 1e-7 relative): max |change| %.2e\n",
            max(abs(v[[1]] - v[[2]]))))

## all circulating at t = 0 (added 2026-11-18 with that change): the windowed
## starting state has the whole inoculum in the circulating compartments
s0 <- plasmofit:::generate_starts_win(42.84, 384L, WINDOW, 2, 2, log10(1800))
cat(sprintf("   generate_starts_win(): circulating %.6g, sequestered %.3g, of total0 %.6g\n",
            sum(s0[1:384]), sum(s0[385:768]), 1800))
stopifnot(abs(sum(s0[1:384]) - 1800) < 1e-9, sum(s0[385:768]) == 0)
## how much that differs from splitting the start by y(age), as the old
## start does: matrix exponential of both states, circulating abundance at
## the observation times (72-216 h by 12 h), cycle 42.84 h, R 13.87, a = b = 2
yv <- plasmofit:::make_y_vals(42.84, 384L); Wv <- s0[1:384]
## From 12 h, to show WHERE it differs: only while the injected cells are in
## their first cycle. Observations start at 72 h.
ts <- seq(12, 216, by = 12)
yc_all <- plasmofit:::mat_exp_series(s0, 42.84, 384L, 13.87, 0, ts, 12)
yc_split <- plasmofit:::mat_exp_series(c(yv * Wv, (1 - yv) * Wv), 42.84, 384L, 13.87, 0, ts, 12)
cat("   all circulating / split by y(age), circulating abundance (cycle 42.84 h, R 13.87, a = b = 2):\n")
cat(sprintf("     t %3.0f h (%.2f cycles): %.8f\n", ts, ts / 42.84, yc_all / yc_split), sep = "")
cat(sprintf("   first-cycle Poisson mass at the first observation (72 h): %.2e\n",
            ppois(383, 384 * 72 / 42.84)))

## ---- 2. start_dist = 0 against the old build ---------------------------------
cat("\n2. start_dist = 0 against the OLD build's saved log density\n",
    "   (_data/start-dist-validate-old.rds, made in the live library BEFORE the\n",
    "   install), production data\n",
    "   (np_wide_total0_nc384_t: no_pool, n_c 384, 1130 observations), at ", K,
    " fixed\n   unconstrained points. Same model, so differences are rounding only.\n\n", sep = "")
ref <- readRDS(REF)
f0 <- handle(withS(d0, 0L))
stopifnot(rstan::get_num_upars(f0) == ref$np)
lp0 <- apply(ref$U, 1, \(u) rstan::log_prob(f0, u))
gr0 <- t(apply(ref$U, 1, \(u) rstan::grad_log_prob(f0, u)))
cat(sprintf("   parameters: %d, as the old build\n", ref$np))
cat(sprintf("   log_prob: max |new - old| %.2e (values %.1f to %.1f)\n",
            max(abs(lp0 - ref$lp)), min(ref$lp), max(ref$lp)))
cat(sprintf("   gradient: max |new - old| / max(1, |old|) %.2e over %d entries\n",
            max(abs(gr0 - ref$gr) / pmax(1, abs(ref$gr))), length(gr0)))

## ---- 3. start_dist = 1 on the production data --------------------------------
cat(sprintf("\n3. start_dist = 1, window %g h, cap %g, production data.\n\n", WINDOW, CAP))
f1 <- handle(withS(d0, 1L))
np1 <- rstan::get_num_upars(f1)
cat(sprintf("   parameters: %d (old %d: minus 14 b_shape_free and 14 x 2 b_off_vec, plus start_alpha and start_beta)\n",
            np1, ref$np))
stopifnot(np1 == ref$np - 14 - 28 + 2)
dims <- f1@par_dims
cat(sprintf("   output dims: b_shape %s, b_offset %s, start_alpha %s, start_beta %s\n",
            paste(dims$b_shape, collapse = "x"), paste(dims$b_offset, collapse = "x"),
            paste(dims$start_alpha, collapse = "x"), paste(dims$start_beta, collapse = "x")))
set.seed(7)
U1 <- matrix(rnorm(K * np1, 0, 0.3), K, np1)
lp1 <- apply(U1, 1, \(u) rstan::log_prob(f1, u))
gr1 <- t(apply(U1, 1, \(u) rstan::grad_log_prob(f1, u)))
cat(sprintf("   log_prob finite at %d of %d points; gradient finite in %d of %d entries\n",
            sum(is.finite(lp1)), K, sum(is.finite(gr1)), length(gr1)))
## where the start parameters sit in the unconstrained vector: they follow
## b_off_vec (now empty), so after b_shape_free (empty): the first two
skel <- rstan::unconstrain_pars(f1, rstan::constrain_pars(f1, U1[1, ]))
stopifnot(isTRUE(all.equal(skel, U1[1, ])))
fd_err <- c()
for (p in 1:2) for (q in 1:3) {
    u <- U1[q, ]; e <- 1e-5
    up <- u; up[p] <- up[p] + e; dn <- u; dn[p] <- dn[p] - e
    fd <- (rstan::log_prob(f1, up) - rstan::log_prob(f1, dn)) / (2 * e)
    fd_err <- c(fd_err, abs(fd - gr1[q, p]) / max(1, abs(gr1[q, p])))
}
cat(sprintf("   start_alpha, start_beta gradient against central differences (3 points each): max relative error %.2e\n",
            max(fd_err)))
cp <- rstan::constrain_pars(f1, U1[1, ])
cat(sprintf("   at point 1: start_alpha %.3f, start_beta %.3f (both inside [1, %g])\n",
            cp$start_alpha, cp$start_beta, CAP))

## ---- 4. matrix_exp cross-check ---------------------------------------------
cat("\n4. run_check at n_c 96: max_rel_diff between the Erlang-window\n",
    "   likelihood path and an independent matrix exponential of the same\n",
    "   starting state, at the Fixed_param initial point.\n\n", sep = "")
for (dist in 0:1) {
    d <- withS(d0, dist); d$n_c <- 96L; d$run_check <- 1L
    fc <- handle(d)
    mrd <- as.matrix(fc, pars = "max_rel_diff")[1, 1]
    cat(sprintf("   start_dist %d: max_rel_diff %.2e\n", dist, mrd))
}

## ---- 5. the production path: wockner-fit.R entry 55 -------------------------
## Added the same day, after entry 55 (np_start13_nc384_t) and the summary
## helpers' support for start_alpha/start_beta went into wockner-fit.R.
## Builds the entry's data list fresh through archer_stan_data(), compares it
## with the saved data of entry 54 (should differ only in the start fields),
## and runs wockner-fit.R's own draw_as_list() + unconstrain_pars() round trip
## on a start_dist = 1 handle and a start_dist = 0 one.
cat("\n5. wockner-fit.R entry 55 (np_start13_nc384_t) through archer_stan_data().\n\n")
suppressPackageStartupMessages({library(dplyr); library(readr)})
paras_df <- read_csv("_data/wockner-cleaned.csv", col_types = "cccdcdd") |>
    mutate(inoc = interaction(trial, inoc_size, drop = TRUE) |> paste(),
           obs_error = interaction(trial, cohort, drop = TRUE) |> paste())
d55 <- archer_stan_data(paras_df, series = "id", time = "time", abundance = "para",
                        grp_init = "inoc", grp_R = "trial", grp_cl = "trial",
                        grp_sd = "obs_error", calc_log_lik = 1L,
                        sd_log10_total0 = 1, n_c = 384L,
                        start_dist = 1L, start_window = 13, start_shape_cap = 3)
same <- vapply(names(d0), \(n) isTRUE(all.equal(d0[[n]], d55[[n]])), logical(1))
cat("   fields of entry 54's saved data that differ in entry 55's:",
    if (all(same)) "none" else paste(names(d0)[!same], collapse = ", "), "\n")
cat("   fields only in entry 55's:", paste(setdiff(names(d55), names(d0)), collapse = ", "), "\n")
stopifnot(all(same), identical(attr(d55, "levels"), attr(d0, "levels")))
bad <- tryCatch({ archer_stan_data(paras_df, series = "id", time = "time", abundance = "para",
                    grp_init = "inoc", grp_R = "trial", grp_cl = "trial", grp_sd = "obs_error",
                    start_dist = 1L); "no error" }, error = conditionMessage)
cat("   start_dist = 1 without a window is refused:", bad, "\n")
## wockner-fit.R's helpers, evaluated alone (not the whole script)
H <- new.env(); H$d <- d55
for (ex in parse("_scripts/wockner-fit.R"))
    if (is.call(ex) && identical(ex[[1]], as.name("<-")) &&
        as.character(ex[[2]]) %in% c("draw_as_list", "par_block_for", "start_pars_for",
                                     "bshape_par_for", "r_pars_for", "cl_pars_for",
                                     "anchor_pars_for")) eval(ex, H)
for (dist in 1:0) {
    dd <- if (dist == 1L) d55 else withS(d0, 0L)
    H$d <- dd
    fh <- handle(dd)
    init <- H$draw_as_list(fh, H$par_block_for("no_pool", fh, dd), 1L)
    u <- rstan::unconstrain_pars(fh, init)
    cat(sprintf("   start_dist %d: draw_as_list + unconstrain_pars round trip, %d coordinates; b_off_vec dims %s\n",
                dist, length(u), paste(dim(init$b_off_vec), collapse = "x")))
    stopifnot(length(u) == rstan::get_num_upars(fh))
}


## ---- 6. start_dist 2 and 3 (added 2026-11-19) ---------------------------------
## 2: a mode per unit and one shared concentration; 3: both per unit. Checks:
## start_dist 1 unchanged by the generalisation (against REF1 from the
## 13db2bd build), the parameter layout, the (mode, conc) -> (a, b) map, finite
## log_prob and gradient with spot finite differences, the matrix_exp
## cross-check, and wockner-fit.R entries 56-57's data lists and helpers.
if (file.exists(REF1)) {
    cat("\n6a. start_dist = 1 against the 13db2bd build (REF1), same 8 points.\n")
    r1 <- readRDS(REF1)
    f1b <- handle(withS(d0, 1L))
    stopifnot(rstan::get_num_upars(f1b) == r1$np)
    lpb <- apply(r1$U, 1, \(u) rstan::log_prob(f1b, u))
    grb <- t(apply(r1$U, 1, \(u) rstan::grad_log_prob(f1b, u)))
    cat(sprintf("   log_prob: max |new - ref| %.2e; gradient: max |new - ref| / max(1, |ref|) %.2e\n",
                max(abs(lpb - r1$lp)), max(abs(grb - r1$gr) / pmax(1, abs(r1$gr)))))
}
for (dist in 2:3) {
    cat(sprintf("\n6%s. start_dist = %d, window %g h, cap %g.\n", c("", "", "b", "c")[dist + 1], dist, WINDOW, CAP))
    fd <- handle(withS(d0, dist))
    npd <- rstan::get_num_upars(fd)
    expect <- ref$np - 42 + 14 + if (dist == 2) 1 else 14
    cat(sprintf("   parameters: %d (expected %d: old %d minus 42 start parameters, plus 14 modes and %d concentration%s)\n",
                npd, expect, ref$np, if (dist == 2) 1 else 14, if (dist == 2) "" else "s"))
    stopifnot(npd == expect)
    set.seed(10 + dist)
    Ud <- matrix(rnorm(K * npd, 0, 0.3), K, npd)
    lpd <- apply(Ud, 1, \(u) rstan::log_prob(fd, u))
    grd <- t(apply(Ud, 1, \(u) rstan::grad_log_prob(fd, u)))
    cat(sprintf("   log_prob finite at %d of %d points; gradient finite in %d of %d entries\n",
                sum(is.finite(lpd)), K, sum(is.finite(grd)), length(grd)))
    cp <- rstan::constrain_pars(fd, Ud[1, ])
    kk <- if (dist == 2) rep(cp$start_conc, 14) else cp$start_conc
    mode_back <- (cp$start_a - 1) / (cp$start_a + cp$start_b - 2)
    cat(sprintf("   map: max |mode from (a, b) - start_mode| %.2e; max |a + b - conc| %.2e; min(a, b) %.3f\n",
                max(abs(mode_back - cp$start_mode)), max(abs(cp$start_a + cp$start_b - kk)),
                min(c(cp$start_a, cp$start_b))))
    coords <- if (dist == 2) c(1, 7, 15) else c(1, 15, 28)
    err <- c()
    for (p in coords) for (q in 1:2) {
        u <- Ud[q, ]; e <- 1e-5; up <- u; up[p] <- up[p] + e; dn <- u; dn[p] <- dn[p] - e
        fdv <- (rstan::log_prob(fd, up) - rstan::log_prob(fd, dn)) / (2 * e)
        err <- c(err, abs(fdv - grd[q, p]) / max(1, abs(grd[q, p])))
    }
    cat(sprintf("   mode/conc gradient against central differences (coordinates %s, 2 points): max relative error %.2e\n",
                paste(coords, collapse = ", "), max(err)))
    dc <- withS(d0, dist); dc$n_c <- 96L; dc$run_check <- 1L
    cat(sprintf("   run_check at n_c 96: max_rel_diff %.2e\n",
                as.matrix(handle(dc), pars = "max_rel_diff")[1, 1]))
    dd <- archer_stan_data(paras_df, series = "id", time = "time", abundance = "para",
                           grp_init = "inoc", grp_R = "trial", grp_cl = "trial",
                           grp_sd = "obs_error", calc_log_lik = 1L, sd_log10_total0 = 1,
                           n_c = 384L, start_dist = dist, start_window = 13, start_shape_cap = 3)
    same <- vapply(names(d0), \(n) isTRUE(all.equal(d0[[n]], dd[[n]])), logical(1))
    H$d <- dd
    fh <- handle(dd)
    init <- H$draw_as_list(fh, H$par_block_for("no_pool", fh, dd), 1L)
    u <- rstan::unconstrain_pars(fh, init)
    cat(sprintf("   wockner-fit.R entry %d: data differs from entry 54's only in the start fields: %s; helper round trip %d coordinates\n",
                if (dist == 2) 56 else 57, all(same), length(u)))
    stopifnot(all(same), length(u) == npd)
}


## ---- 7. start_below and start_delay (added 2026-11-20) -----------------------
## For wockner-fit.R entries 58 (start_dist 3 on [-6, 13] h) and 59 (start_dist
## 1 on [0, 13] h plus a per-unit delay on [0, 12] h). claude/start-dist-design.md,
## "Below age 0 and a delay".
if (exists("windowed_beta_starts_below", envir = asNamespace("plasmofit"))) {
BELOW <- 6; DMAX <- 12
withX <- function(d, dist, below = 0, delay = FALSE, dmax = DMAX) {
    d <- withS(d, dist); d$start_below <- below
    d$start_delay <- as.integer(delay); d$delay_max <- if (delay) dmax else 0; d
}
cat("\n7a. start_dist 2 and 3 against the fc7db77 build (REF23, made in the scratch\n",
    "    library before this build replaced it), 8 points each, 13 h window.\n", sep = "")
r23 <- readRDS(REF23)
for (dist in 2:3) {
    r <- r23[[as.character(dist)]]
    f <- handle(withX(d0, dist))
    stopifnot(rstan::get_num_upars(f) == r$np)
    lp <- apply(r$U, 1, \(u) rstan::log_prob(f, u)); gr <- t(apply(r$U, 1, \(u) rstan::grad_log_prob(f, u)))
    cat(sprintf("   start_dist %d: log_prob max |new - ref| %.2e; gradient max |new - ref| / max(1, |ref|) %.2e\n",
                dist, max(abs(lp - r$lp)), max(abs(gr - r$gr) / pmax(1, abs(r$gr)))))
}

cat("\n7b. windowed_beta_starts_below() against exact pbeta differences over each\n",
    "    compartment's overlap with [0, w] (age = position) and [1 - wb, 1]\n",
    "    (age = position - 1), window [-6, 13] h. Cells: max |difference| in\n",
    "    compartment weight (weights sum to 1).\n\n", sep = "")
exact_below <- function(n_c, wb, w, a, b) {
    h <- 1 / n_c; sp <- wb + w; out <- numeric(n_c)
    P <- function(age_lo, age_hi) pbeta((age_hi + wb) / sp, a, b) - pbeta((age_lo + wb) / sp, a, b)
    for (i in seq_len(n_c)) { lo <- (i - 1) * h; hi <- i * h
        if (lo < w) out[i] <- out[i] + P(lo, min(hi, w))
        if (wb > 0 && hi > 1 - wb) out[i] <- out[i] + P(max(lo, 1 - wb) - 1, hi - 1) }
    out / sum(out)
}
worst <- 0; worst0 <- 0
for (nc in c(96L, 384L)) for (cl in c(39.4, 42.84, 45.6))
    for (ab in list(c(1, 1), c(1.2, 1.2), c(1.5, 3), c(2, 2), c(3, 1), c(3, 3))) {
        g <- plasmofit:::windowed_beta_starts_below(nc, BELOW / cl, WINDOW / cl, ab[1], ab[2], 1)
        e <- exact_below(nc, BELOW / cl, WINDOW / cl, ab[1], ab[2])
        pos <- (seq_len(nc) - 0.5) / nc
        stopifnot(abs(sum(g) - 1) < 1e-12, all(g >= 0),
                  all(g[pos > (WINDOW / cl) + 1 / nc & pos < 1 - BELOW / cl - 1 / nc] == 0))
        worst <- max(worst, max(abs(g - e)))
        ## with wb = 0 it must give windowed_beta_starts()
        g0 <- plasmofit:::windowed_beta_starts_below(nc, 0, WINDOW / cl, ab[1], ab[2], 1)
        worst0 <- max(worst0, max(abs(g0 - plasmofit:::windowed_beta_starts(nc, WINDOW / cl, ab[1], ab[2], 1))))
    }
cat(sprintf("   worst |below - exact| over n_c {96, 384} x 3 cycle lengths x 6 shapes: %.2e\n", worst))
cat(sprintf("   worst |below with wb = 0 - windowed_beta_starts()| over the same: %.2e\n", worst0))
nc <- 384L; edge_cl <- BELOW / (20 / nc)            # 1 - wb exactly on a compartment edge
v <- lapply(edge_cl * (1 + c(-1e-7, 1e-7)), \(cl) plasmofit:::windowed_beta_starts_below(nc, BELOW / cl, WINDOW / cl, 1.5, 3, 1))
cat(sprintf("   continuity as the lower edge crosses a compartment edge (cycle_length %.2f h +/- 1e-7 relative): max |change| %.2e\n",
            edge_cl, max(abs(v[[1]] - v[[2]]))))
gb <- plasmofit:::windowed_beta_starts_below(nc, BELOW / 42.84, WINDOW / 42.84, 1.5, 3, 1)
cat(sprintf("   share of the start below age 0 (last 6 h of the cycle), a 1.5 b 3, cycle 42.84 h: %.3f (exact pbeta(6/19) = %.3f)\n",
            sum(gb[((seq_len(nc) - 1) / nc) >= 0.5]), pbeta(6 / 19, 1.5, 3)))

cat("\n7c. All circulating against split by y(age) at t = 0 when part of the start\n",
    "    is late schizonts (window [-6, 13] h, a 1.5 b 3, cycle 42.84 h, R 13.87,\n",
    "    n_c 384): circulating abundance by matrix exponential, ratio all/split.\n\n", sep = "")
s0 <- plasmofit:::generate_starts_gen(42.84, 384L, BELOW, WINDOW, 1.5, 3, log10(1800))
yv <- plasmofit:::make_y_vals(42.84, 384L); Wv <- s0[1:384]
ts <- seq(12, 216, by = 12)
ra <- plasmofit:::mat_exp_series(s0, 42.84, 384L, 13.87, 0, ts, 12) /
      plasmofit:::mat_exp_series(c(yv * Wv, (1 - yv) * Wv), 42.84, 384L, 13.87, 0, ts, 12)
cat(sprintf("     t %3.0f h: %.10f\n", ts, ra), sep = "")
cat(sprintf("   max |ratio - 1| over the observation times 72-216 h: %.2e\n", max(abs(ra[ts >= 72] - 1))))

cat("\n7d. Entry 58's layout: start_dist 3, window [-6, 13] h, cap 3.\n")
f58 <- handle(withX(d0, 3L, below = BELOW))
np58 <- rstan::get_num_upars(f58)
cat(sprintf("   parameters: %d (start_dist 3 without start_below: %d)\n", np58, r23[["3"]]$np))
stopifnot(np58 == r23[["3"]]$np)
U <- r23[["3"]]$U
lp58 <- apply(U, 1, \(u) rstan::log_prob(f58, u)); gr58 <- t(apply(U, 1, \(u) rstan::grad_log_prob(f58, u)))
cat(sprintf("   log_prob finite at %d of %d points; gradient finite in %d of %d entries\n",
            sum(is.finite(lp58)), K, sum(is.finite(gr58)), length(gr58)))
cat(sprintf("   log_prob minus start_dist 3 on [0, 13] h at the same points (should be nonzero): %s\n",
            paste(sprintf("%.1f", lp58 - r23[["3"]]$lp), collapse = ", ")))
err <- c()
for (p in c(1, 15, 28)) for (q in 1:2) {
    u <- U[q, ]; e <- 1e-5; up <- u; up[p] <- up[p] + e; dn <- u; dn[p] <- dn[p] - e
    fdv <- (rstan::log_prob(f58, up) - rstan::log_prob(f58, dn)) / (2 * e)
    err <- c(err, abs(fdv - gr58[q, p]) / max(1, abs(gr58[q, p])))
}
cat(sprintf("   mode/conc gradient against central differences (coordinates 1, 15, 28; 2 points): max relative error %.2e\n", max(err)))

cat("\n7e. Entry 59's layout: start_dist 1 on [0, 13] h plus a delay per unit on [0, 12] h.\n")
d59x <- withX(d0, 1L, delay = TRUE)
f59 <- handle(d59x)
np59 <- rstan::get_num_upars(f59); r1 <- readRDS(REF1)
cat(sprintf("   parameters: %d (start_dist 1 without the delay %d, plus 14)\n", np59, r1$np))
stopifnot(np59 == r1$np + 14)
set.seed(59)
U59 <- matrix(rnorm(K * np59, 0, 0.3), K, np59)
lp59 <- apply(U59, 1, \(u) rstan::log_prob(f59, u)); gr59 <- t(apply(U59, 1, \(u) rstan::grad_log_prob(f59, u)))
cat(sprintf("   log_prob finite at %d of %d points; gradient finite in %d of %d entries\n",
            sum(is.finite(lp59)), K, sum(is.finite(gr59)), length(gr59)))
## the delay follows start_beta in the parameter order: coordinates 3-16
err <- c()
for (p in c(1, 3, 9, 16)) for (q in 1:2) {
    u <- U59[q, ]; e <- 1e-5; up <- u; up[p] <- up[p] + e; dn <- u; dn[p] <- dn[p] - e
    fdv <- (rstan::log_prob(f59, up) - rstan::log_prob(f59, dn)) / (2 * e)
    err <- c(err, abs(fdv - gr59[q, p]) / max(1, abs(gr59[q, p])))
}
cp <- rstan::constrain_pars(f59, U59[1, ])
cat(sprintf("   at point 1, delay (h): %s\n", paste(sprintf("%.2f", cp$delay), collapse = " ")))
cat(sprintf("   start_alpha and delay gradient against central differences (coordinates 1, 3, 9, 16; 2 points): max relative error %.2e\n", max(err)))
## delay -> 0 must give start_dist 1 without a delay: same other parameters,
## log density WITHOUT the Jacobian (the delay's flat prior is then constant 0)
f1n <- handle(withX(d0, 1L))
pl <- rstan::constrain_pars(f1n, r1$U[1, ])
dz <- c()
for (dl in c(1e-6, 1e-9)) {
    pl$delay <- rep(dl, 14)
    a <- rstan::log_prob(f59, rstan::unconstrain_pars(f59, pl), adjust_transform = FALSE)
    b <- rstan::log_prob(f1n, r1$U[1, ], adjust_transform = FALSE)
    dz <- c(dz, a - b)
}
cat(sprintf("   log density at delay 1e-6 h and 1e-9 h minus start_dist 1 without a delay (Jacobian off): %.2e, %.2e\n", dz[1], dz[2]))
## log_lik in generated quantities against the model block: two points
## differing only in the delays, where every prior is flat in what changed, so
## the change in log density (Jacobian off) must equal the change in sum(log_lik)
pA <- rstan::constrain_pars(f59, U59[1, ]); pB <- pA; pB$delay <- rev(pA$delay)
ll <- sapply(list(pA, pB), \(pp) sum(as.matrix(handle(d59x, init = list(pp)), pars = "log_lik")[1, ]))
lpd <- sapply(list(pA, pB), \(pp) rstan::log_prob(f59, rstan::unconstrain_pars(f59, pp), adjust_transform = FALSE))
cat(sprintf("   change in sum(log_lik) %.6f against change in log density %.6f (difference %.1e)\n",
            diff(ll), diff(lpd), diff(ll) - diff(lpd)))
## the same for entry 58: two points differing only in the start modes
qA <- rstan::constrain_pars(f58, U[1, ]); qB <- qA; qB$start_mode <- rev(qA$start_mode)
d58x <- withX(d0, 3L, below = BELOW)
ll <- sapply(list(qA, qB), \(pp) sum(as.matrix(handle(d58x, init = list(pp)), pars = "log_lik")[1, ]))
lpd <- sapply(list(qA, qB), \(pp) rstan::log_prob(f58, rstan::unconstrain_pars(f58, pp), adjust_transform = FALSE))
cat(sprintf("   entry 58: change in sum(log_lik) %.6f against change in log density %.6f (difference %.1e)\n",
            diff(ll), diff(lpd), diff(ll) - diff(lpd)))

## 7f was first the Stan run_check (mat_exp_at, one matrix-exponential action
## per observation); one Fixed_param handle at n_c 96 had not finished after
## 72 minutes, so the same comparison is made here per series in R with the
## exposed functions instead: the Erlang-window path ew_poly_series_gen(), with
## truncation bounds computed as transformed data computes them (from
## ts - delay_max, min_cl, max_cl), against mat_exp_series() of the same
## all-circulating start on a 0.5 h grid, read at ts - delay. Delays are
## multiples of 0.5 h so the read times are on that grid. Cycle lengths at
## and near both bounds stress the truncation.
cat("\n7f. Erlang-window path against matrix exponential, per series in R:\n",
    "    max over 7 read times (72-216 h by 24 h) of |fast - matrix_exp| / fast,\n",
    "    log10_total0 1, R 12, mu 0. delay_max 12 when a delay is on.\n\n", sep = "")
bounds <- function(ts, n_c, dmax) {
    M <- ceiling(max(d0$ts) / d0$min_cl) + 4; r_max <- M * n_c + n_c
    tmin <- n_c * (ts - dmax) / d0$max_cl; tmax <- n_c * ts / d0$min_cl
    r_hi <- ceiling(tmax + 10 * sqrt(tmax)); r_lo <- pmax(1, floor(tmin - 10 * sqrt(tmin)))
    list(M = as.integer(M), r_max = as.integer(r_max), lgam = lgamma(seq_len(r_max)),
         r_less1 = as.numeric(seq_len(r_max) - 1), r_lo = as.integer(r_lo), r_hi = as.integer(r_hi),
         m_lo = as.integer(pmax(1, ceiling(r_lo / n_c - 1))), m_hi = as.integer(floor((r_hi - 2) / n_c + 1)))
}
tsv <- seq(72, 216, by = 24); worst7f <- 0
for (lay in list(list("58: [-6, 13] h, no delay", BELOW, 0, 0, 3, 1.5),
                 list("59: [0, 13] h, delay 0.5 h", 0, 0.5, DMAX, 2, 2),
                 list("59: [0, 13] h, delay 6.5 h", 0, 6.5, DMAX, 2, 2),
                 list("59: [0, 13] h, delay 11.5 h", 0, 11.5, DMAX, 1.2, 3),
                 list("both: [-6, 13] h, delay 6.5 h", BELOW, 6.5, DMAX, 1.5, 3)))
    for (nc in c(96L, 384L)) {
        bb <- bounds(tsv, nc, lay[[4]])
        m <- sapply(c(35.2, 42.84, 49.8), \(cl) {
            fast <- plasmofit:::ew_poly_series_gen(cl, nc, lay[[2]], WINDOW, lay[[5]], lay[[6]], lay[[3]],
                        1, 12, 0, length(tsv), tsv, bb$M, bb$r_max, bb$lgam, bb$r_less1,
                        bb$r_lo, bb$r_hi, bb$m_lo, bb$m_hi)
            y0 <- plasmofit:::generate_starts_gen(cl, nc, lay[[2]], WINDOW, lay[[5]], lay[[6]], 1)
            mx <- plasmofit:::mat_exp_series(y0, cl, nc, 12, 0, tsv - lay[[3]], 0.5)
            max(abs(fast - mx) / fast)
        })
        worst7f <- max(worst7f, m)
        cat(sprintf("   %-32s n_c %3d, cycle 35.2 / 42.84 / 49.8 h: %s\n", lay[[1]], nc,
                    paste(sprintf("%.2e", m), collapse = " / ")))
    }
cat(sprintf("   worst: %.2e\n", worst7f))

cat("\n7g. wockner-fit.R entries 58 and 59 through archer_stan_data(), and refusals.\n")
mk <- function(...) archer_stan_data(paras_df, series = "id", time = "time", abundance = "para",
                                     grp_init = "inoc", grp_R = "trial", grp_cl = "trial",
                                     grp_sd = "obs_error", calc_log_lik = 1L, sd_log10_total0 = 1,
                                     n_c = 384L, start_window = 13, start_shape_cap = 3, ...)
for (en in list(list(58, mk(start_dist = 3L, start_below = 6), np58),
                list(59, mk(start_dist = 1L, start_delay = TRUE, delay_max = 12), np59))) {
    dd <- en[[2]]
    same <- vapply(names(d0), \(n) isTRUE(all.equal(d0[[n]], dd[[n]])), logical(1))
    H$d <- dd; fh <- handle(dd)
    init <- H$draw_as_list(fh, H$par_block_for("no_pool", fh, dd), 1L)
    u <- rstan::unconstrain_pars(fh, init)
    cat(sprintf("   entry %d: fields differing from entry 54's saved data: %s; helper round trip %d coordinates (expected %d)\n",
                en[[1]], paste(names(d0)[!same], collapse = ", "), length(u), en[[3]]))
    stopifnot(all(names(d0)[!same] %in% c("start_dist", "start_window")), length(u) == en[[3]])
}
for (bad in list(list(start_dist = 0L, start_below = 6), list(start_dist = 1L, start_delay = TRUE, delay_max = 72),
                 list(start_dist = 1L, start_delay = TRUE), list(start_dist = 3L, start_below = 30)))
    cat("   refused:", paste(names(bad), unlist(bad), sep = " = ", collapse = ", "), "->",
        tryCatch({ do.call(mk, bad); "NOT REFUSED" }, error = conditionMessage), "\n")
}
