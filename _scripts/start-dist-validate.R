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
##   cd /home2/lan68/plasmofit/plasmofit-test
##   START_VALIDATE=old Rscript --vanilla _scripts/start-dist-validate.R
##   START_VALIDATE=new Rscript --vanilla _scripts/start-dist-validate.R > _data/start-dist-validate-<date>.txt

MODE <- Sys.getenv("START_VALIDATE", "new")
LIVE <- "/home/lan68/R/x86_64-pc-linux-gnu-library/4.6"
DEV <- "/home2/lan68/plasmofit/.Rlib-dev"
.libPaths(if (MODE == "old") LIVE else DEV)
suppressPackageStartupMessages({library(plasmofit); library(rstan)})
rstan_options(threads_per_chain = 1L)
REF <- "_data/start-dist-validate-old.rds"
K <- 8L
WINDOW <- 13; CAP <- 3

d0 <- readRDS("_data/wock-data-np_wide_total0_nc384_t.rds")
handle <- function(d, model = "archer_fit")
    suppressWarnings(rstan::sampling(plasmofit:::stanmodels[[model]], data = d,
                                     chains = 1, iter = 1, refresh = 0, seed = 1,
                                     algorithm = "Fixed_param"))

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
stopifnot(startsWith(find.package("plasmofit"), DEV),
          exists("windowed_beta_starts", envir = asNamespace("plasmofit")))
withS <- function(d, dist, window = WINDOW, cap = CAP) {
    d$start_dist <- as.integer(dist); d$start_window <- if (dist == 1L) window else 1
    d$start_shape_cap <- cap; d
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

## ---- 2. start_dist = 0 against the old build ---------------------------------
cat("\n2. start_dist = 0 against the OLD build (live library), production data\n",
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
