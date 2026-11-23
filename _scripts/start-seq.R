## Profile likelihood of the SEQUESTRATION MIDPOINT p3, with the starting
## distribution of entries 54, 57 and 59 in turn. 2026-11-20, Lucas ("1-2")
## after the start-model maxima (claude/findings.md, "Maximum likelihood of
## the start models") showed the ring-only [0, 13] h window costs 29.8
## log-likelihood units against entry 54's whole-cycle start at equal
## parameter count.
##
## THE QUESTION. The model's circulating fraction is a fixed logistic in
## developmental age, y(age) = (1 - p1) / (1 + 10^(-p4 (p3 - age))), with
## p1 = 11.3869/467.6209, p3 = 18.5802 h and p4 = 0.2242 inherited from
## Archer et al. and NEVER estimated here. p3 sets how much of each cycle is
## visible. The starting distribution sets where in the cycle the population
## sits. Both bear on the phase and amplitude of the observed oscillation, and
## only one of them has ever been free: if p3 is wrong, the start has to move
## to compensate, and the ring-only start would be rejected for the wrong
## reason. This profiles p3 (shared by every trial and unit, as a parasite
## property) with each start model, which no run so far has done.
##
## They are separable in principle: moving the start translates the
## trajectory, while moving p3 changes the duty cycle -- mean y over one
## cycle is 0.278 at p3 = 12 h, 0.430 at 18.58 h, 0.603 at 26 h (cycle 42 h,
## n_c 384) -- so it changes the oscillation's shape, not only its position.
## Whether they are separable IN THESE DATA is what the flatness check below
## asks.
##
## READING RULE, fixed before any result (2026-11-20). "Maximum" means the
## best found over all starts, summed over 13 trials, at one p3; "profile"
## means that quantity as a function of p3; "peak" its argmax on the grid.
##   a. If entry 57's profile peak comes within 5 log-likelihood units of
##      entry 54's profile peak, the ring-only start is compatible with
##      these data once the sequestration schedule is free, and the 29.8-unit
##      deficit was an artefact of fixing p3 at 18.58 h.
##   b. If 57's peak is at a p3 more than 3 h from 18.58 h while 54's peak is
##      within 3 h of it, the start was absorbing schedule misspecification.
##   c. If both peak at about the same p3 and the gap stays near 29.8 units,
##      the deficit is about the start, not the schedule.
##   d. FLATNESS, checked first: if a model's profile varies by less than 2
##      units across the whole 12-26 h grid, p3 is not identified under that
##      start and nothing about the schedule follows from it. A start free to
##      move over the whole cycle (54) may well be flat where a constrained
##      one (57) is not; that contrast is itself the answer to the question.
##   e. p3's best value is a point estimate with no priors and no hierarchy.
##      It is NOT a reportable sequestration time, and p4 (the slope) stays
##      fixed throughout, so this profiles the schedule's POSITION only.
##
## METHOD. As _scripts/start-ml.R: priors ignored, each trial x cohort sd
## profiled (set to the root mean squared residual of its group), the same
## parameter bounds as the Stan programs, and the likelihood maximised by
## Nelder-Mead then BFGS from many starts. cycle_length and R are per trial
## and grp_sd nests in trial, so at a given p3 the problem is a sum of 13
## independent per-trial problems; p3 is shared, so the profile sums the
## per-trial maxima at each grid point. Every value is a lower bound on the
## true maximum.
##
## The forward map is the package's ew_poly_eval() fed a D/C pair built in R
## so that p3 is free (make_DC_free below). At the default p3 this reproduces
## plasmofit:::ew_poly_series() and ew_poly_series_gen() to 4e-16 relative,
## and seq_log_y() reproduces plasmofit:::make_log_y_vals() to 4e-15; check
## mode runs both, and also checks that the profile at the default p3 recovers
## start-ml.R's maximum for the same (model, trial).
##
## The parameterisation helpers (unpack, pack, rand_theta, draw_nat,
## load_draws) are COPIED VERBATIM from _scripts/start-ml.R, with the forward
## map replaced; the check above is what guarantees they still agree.
##
## Starts, per trial and p3 grid point: 12 random (cycle_length stratified
## over [35.5, 49.5] h, the rest uniform), the best draw of each of the fitted
## entry's 4 chains, a warm start from the best solution at the previous p3,
## and start-ml.R's own best solution. Entry 59 additionally profiles its
## SHARED (alpha, beta) on a 5 x 5 grid at each p3, with 4 random starts plus
## warm starts from the previous (alpha, beta) and the previous p3.
##
## MODES (SSQ_MODE):
##   check  the three verifications above
##   run    one task (SLURM_ARRAY_TASK_ID 1-39: model x trial), output
##          _data/start-seq/task-<model>-<trial>.rds
##   merge  combine and print the profiles
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   SSQ_MODE=check Rscript --vanilla _scripts/start-seq.R > _data/start-seq-check-<date>.txt
##   sbatch _scripts/start-seq.sh
##   SSQ_MODE=merge Rscript --vanilla _scripts/start-seq.R > _data/start-seq-<date>.txt

MODE <- Sys.getenv("SSQ_MODE", "check")
OUTDIR <- "_data/start-seq/"
P3_GRID <- sort(c(seq(12, 26, by = 1), 18.5802))
GRID_A <- c(1, 1.2, 1.5, 2, 3); GRID_B <- c(1, 1.5, 2, 2.5, 3)
N_RAND <- 12; N_RAND_G <- 4
source("_scripts/start-seq-core.R")
if (nzchar(Sys.getenv("SSQ_SMOKE"))) {          # quick end-to-end test of run mode
    OUTDIR <- paste0(Sys.getenv("SSQ_SMOKE"), "/"); N_RAND <- 2; N_RAND_G <- 2
    P3_GRID <- c(15, P3_DEF); GRID_A <- c(1, 1.5); GRID_B <- 2
}
## SSQ_P3ONLY=<i> restricts run mode to the i-th grid point and writes
## task-<model>-<trial>-p<p3>.rds instead of the whole-grid file, so one trial
## can be spread over the grid as parallel tasks (see _scripts/start-seq-p38.sh).
P3_ONLY <- Sys.getenv("SSQ_P3ONLY", "")
if (nzchar(P3_ONLY)) P3_GRID <- P3_GRID[as.integer(P3_ONLY)]

## ---- check ----------------------------------------------------------------------
if (MODE == "check") {
    cat("1. seq_log_y() against plasmofit:::make_log_y_vals() at the default p3.\n",
        "   Cells: max |difference| in log y over the n_c compartments.\n\n", sep = "")
    for (cl in c(35.2, 42.0, 49.8))
        cat(sprintf("   cycle %.1f h: %.2e\n", cl, max(abs(seq_log_y(cl, P3_DEF) - plasmofit:::make_log_y_vals(cl, n_c)))))
    cat("\n2. The R D/C path against the package's own series at the default p3.\n",
        "   Cells: max relative |difference| over t = 72-216 h by 12 h, cycle 42 h, R 12.\n\n", sep = "")
    ts <- seq(72, 216, by = 12); b <- bounds(ts, 0); cl <- 42; ly <- seq_log_y(cl, P3_DEF)
    ev <- function(W, ty) plasmofit:::ew_poly_eval(make_DC_free(ly, W, ty), n_c, cl, 12, 0,
              length(ts), ts, M, R_MAX, LG, RL1, b$rl, b$rh, b$ml, b$mh)
    f1 <- ev(plasmofit:::beta_starts(n_c, 25, 0.3, 10), "wrap")
    g1 <- plasmofit:::ew_poly_series(cl, n_c, 25, 0.3, 1, 12, 0, length(ts), ts, M, R_MAX, LG, RL1, b$rl, b$rh, b$ml, b$mh)
    cat(sprintf("   entry 54's start (b_shape 25, b_offset 0.3): %.2e\n", max(abs(f1 - g1) / g1)))
    f2 <- ev(plasmofit:::windowed_beta_starts_below(n_c, 0, WINDOW / cl, 1.2, 2.5, 10), "win")
    g2 <- plasmofit:::ew_poly_series_gen(cl, n_c, 0, WINDOW, 1.2, 2.5, 0, 1, 12, 0, length(ts), ts, M, R_MAX, LG, RL1, b$rl, b$rh, b$ml, b$mh)
    cat(sprintf("   windowed [0, 13] h start (a 1.2, b 2.5): %.2e\n", max(abs(f2 - g2) / g2)))
    cat("\n3. This script's likelihood at the DEFAULT p3 against start-ml.R's own best,\n",
        "   evaluated at start-ml's best parameter vector for the same (model, trial).\n",
        "   Cells: the two per-trial log-likelihoods and their difference.\n\n", sep = "")
    for (m in MODELS) {
        dr <- load_draws(m)
        for (c in c(1, 7, 13)) {
            td <- trial_data(c, DMAX[[m]]); nu <- length(td$units)
            x <- readRDS(sprintf("%stask-%s-%02d.rds", MLDIR, m, c))
            kb <- names(x$runs)[which.max(vapply(names(x$runs), \(k) max(vapply(x$runs[[k]], `[[`, 0, "ll")), 0))]
            r <- x$runs[[kb]][[which.max(vapply(x$runs[[kb]], `[[`, 0, "ll"))]]
            ab <- if (m == "59") as.numeric(strsplit(kb, "_")[[1]]) else NULL
            p <- unpack(r$par, m, nu, ab)
            v <- trial_ll(td, p$cl, p$R, p$lt0, p$st, P3_DEF)
            cat(sprintf("   entry %s trial %-13s: here %10.4f | start-ml %10.4f | difference %.2e\n",
                        m, lv$trial[c], v, r$ll, v - r$ll))
        }
    }
    cat("\n4. Duty cycle (mean y over one cycle, cycle 42 h) across the p3 grid.\n\n")
    cat(sprintf("   p3 %5.2f h: %.4f\n", P3_GRID, vapply(P3_GRID, \(p) mean(exp(seq_log_y(42, p))), 0)), sep = "")
    quit(save = "no")
}

## ---- run ------------------------------------------------------------------------
if (MODE == "run") {
    task <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", NA))
    stopifnot(!is.na(task), task >= 1, task <= length(MODELS) * d$n_grp_cl)
    m <- MODELS[(task - 1) %/% d$n_grp_cl + 1]; c <- (task - 1) %% d$n_grp_cl + 1
    dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
    OUT <- if (nzchar(P3_ONLY)) sprintf("%stask-%s-%02d-p%s.rds", OUTDIR, m, c, sprintf("%g", P3_GRID[1]))
           else sprintf("%stask-%s-%02d.rds", OUTDIR, m, c)
    if (file.exists(OUT)) stop(OUT, " exists; refusing to overwrite")
    set.seed(20261120 + 1000 + task)
    td <- trial_data(c, DMAX[[m]]); nu <- length(td$units)
    cat(sprintf("task %d: entry %s, trial %s, %d units, %d observations, %d p3 points\n",
                task, m, lv$trial[c], nu, length(td$k), length(P3_GRID)))
    dr <- load_draws(m)
    prof <- vapply(seq_along(dr$chain), \(i) { nat <- draw_nat(m, dr, i, td)
        st <- unpack(pack(m, nat$cl, nat$R, nat$lt0, nat$extra), m, nu, nat$ab)
        trial_ll(td, st$cl, st$R, st$lt0, st$st, P3_DEF) }, 0)
    seeds <- lapply(1:4, \(ch) which(dr$chain == ch)[which.max(prof[dr$chain == ch])])
    seed_th <- lapply(seeds, \(i) { nat <- draw_nat(m, dr, i, td); pack(m, nat$cl, nat$R, nat$lt0, nat$extra) })
    mlb <- ml_best(m, c)
    runs <- list()
    solve_at <- function(ab, starts, p3) {
        fn <- function(th) { p <- unpack(th, m, nu, ab)
            v <- -trial_ll(td, p$cl, p$R, p$lt0, p$st, p3); if (is.finite(v)) v else 1e10 }
        lapply(starts, \(s) { o <- opt_one(s$th, fn); p <- unpack(o$par, m, nu, ab)
                              list(src = s$src, ll = -o$value, par = o$par, cl = p$cl) })
    }
    prev_p3 <- NULL
    for (p3 in P3_GRID) {
        if (m != "59") {
            st <- c(lapply(seq_len(N_RAND), \(i) list(src = "random", th = rand_theta(m, nu, i, N_RAND))),
                    lapply(1:4, \(ch) list(src = sprintf("chain%d", ch), th = seed_th[[ch]])),
                    if (!is.null(mlb)) list(list(src = "startml", th = mlb$th)),
                    if (!is.null(prev_p3)) list(list(src = "warm_p3", th = prev_p3)))
            r <- solve_at(NULL, st, p3)
            prev_p3 <- r[[which.max(vapply(r, `[[`, 0, "ll"))]]$par
            runs[[sprintf("%g", p3)]] <- list(none = r)
            cat(sprintf("   p3 %6.3f: best %.3f\n", p3, max(vapply(r, `[[`, 0, "ll")))); flush.console()
        } else {
            sub <- list(); prev_ab <- NULL; bestp <- NULL; bestll <- -Inf
            for (a in GRID_A) for (b in GRID_B) {
                st <- c(lapply(seq_len(N_RAND_G), \(i) list(src = "random", th = rand_theta(m, nu, i, N_RAND_G))),
                        lapply(1:2, \(ch) list(src = sprintf("chain%d", ch), th = seed_th[[ch]])),
                        if (!is.null(prev_ab)) list(list(src = "warm_ab", th = prev_ab)),
                        if (!is.null(prev_p3)) list(list(src = "warm_p3", th = prev_p3)),
                        if (!is.null(mlb)) list(list(src = "startml", th = mlb$th)))
                r <- solve_at(c(a, b), st, p3)
                ib <- which.max(vapply(r, `[[`, 0, "ll")); prev_ab <- r[[ib]]$par
                if (r[[ib]]$ll > bestll) { bestll <- r[[ib]]$ll; bestp <- r[[ib]]$par }
                sub[[sprintf("%g_%g", a, b)]] <- r
            }
            prev_p3 <- bestp
            runs[[sprintf("%g", p3)]] <- sub
            cat(sprintf("   p3 %6.3f: best %.3f\n", p3, bestll)); flush.console()
        }
    }
    saveRDS(list(model = m, trial = c, units = td$units, p3_grid = P3_GRID, runs = runs,
                 ml_best = if (is.null(mlb)) NA else mlb$ll), OUT)
    cat("written", OUT, "\n")
    quit(save = "no")
}

## ---- merge ---------------------------------------------------------------------
## Entry 59's (alpha, beta) is SHARED by every trial, so its profile at a given
## p3 is max over (alpha, beta) of the SUM over trials, not the sum of each
## trial's own best (alpha, beta). Entries 54 and 57 share nothing across
## trials, so for them the two agree.
##
## A trial whose whole-grid file is absent is rebuilt from per-point files
## written by SSQ_P3ONLY runs; section 1 names any trial rebuilt that way and
## any grid point still missing.
if (MODE == "merge") {
    nT <- d$n_grp_cl; TR <- lv$trial
    kk <- function(x) sprintf("%g", P3_GRID[vapply(x, \(v) which.min(abs(P3_GRID - v)), 1L)])
    allk <- kk(P3_GRID); dkey <- kk(P3_DEF)
    best_starts <- function(r) max(vapply(r, `[[`, 0, "ll"))
    rebuilt <- character()
    ## runs_of(m, c): list by p3 key of (named vector over sub-keys of the best
    ## over starts). Sub-key is "none" for entries 54 and 57, "<a>_<b>" for 59.
    runs_of <- function(m, c) {
        f <- sprintf("%stask-%s-%02d.rds", OUTDIR, m, c)
        if (file.exists(f)) return(lapply(readRDS(f)$runs, \(sub) vapply(sub, best_starts, 0)))
        pf <- list.files(OUTDIR, sprintf("^task-%s-%02d-p.*\\.rds$", m, c), full.names = TRUE)
        if (!length(pf)) stop("no saved run for entry ", m, " trial ", c)
        rebuilt <<- c(rebuilt, sprintf("entry %s, trial %s: %d of %d points from per-point runs",
                                       m, TR[c], length(pf), length(P3_GRID)))
        out <- list()
        for (q in pf) { x <- readRDS(q)
                        for (k in names(x$runs)) out[[k]] <- vapply(x$runs[[k]], best_starts, 0) }
        out
    }
    LL <- lapply(MODELS, \(m) lapply(seq_len(nT), \(c) runs_of(m, c)))
    names(LL) <- MODELS
    ## profile: per p3 key present for every trial, max over shared sub-keys of
    ## the summed per-trial best. Also keep which sub-key attained it.
    prof_of <- function(m) {
        keys <- Reduce(intersect, lapply(LL[[m]], names))
        v <- numeric(length(keys)); sk <- character(length(keys)); names(v) <- keys
        for (j in seq_along(keys)) {
            subs <- Reduce(intersect, lapply(LL[[m]], \(x) names(x[[keys[j]]])))
            tot <- vapply(subs, \(sb) sum(vapply(LL[[m]], \(x) x[[keys[j]]][[sb]], 0)), 0)
            v[j] <- max(tot); sk[j] <- subs[which.max(tot)]
        }
        list(v = v, sub = setNames(sk, keys))
    }
    PR <- lapply(MODELS, prof_of); names(PR) <- MODELS
    PROF <- lapply(PR, `[[`, "v")
    ## start-ml.R's maximum, with (alpha, beta) shared the same way
    ml_shared <- function(m) {
        per <- lapply(seq_len(nT), \(c) { f <- sprintf("%stask-%s-%02d.rds", MLDIR, m, c)
            if (!file.exists(f)) return(NULL)
            vapply(readRDS(f)$runs, best_starts, 0) })
        if (any(vapply(per, is.null, TRUE))) return(NA_real_)
        subs <- Reduce(intersect, lapply(per, names))
        max(vapply(subs, \(sb) sum(vapply(per, \(x) x[[sb]], 0)), 0))
    }
    npar <- c(`54` = 95, `57` = 95, `59` = 83)
    cat("1. PROFILE of the sequestration midpoint p3. Cells: maximum log-likelihood\n",
        "   over 1130 observations at that p3, summed over 13 trials, sd per trial x\n",
        "   cohort profiled, priors ignored, p4 fixed at 0.2242. p3 is SHARED by every\n",
        "   trial and unit, and so is entry 59's start shape (alpha, beta), which is\n",
        "   therefore chosen once for all 13 trials at each p3. Higher is better; every\n",
        "   value is a lower bound on the true maximum. The row at 18.5802 h is the\n",
        "   value the model has always used. NA = that grid point is not available for\n",
        "   every trial in that model.\n\n", sep = "")
    if (length(rebuilt)) cat("   rebuilt from per-point runs: ", paste(unique(rebuilt), collapse = "; "), "\n\n", sep = "")
    tab <- data.frame(p3_h = P3_GRID)
    for (m in MODELS) tab[[m]] <- round(unname(PROF[[m]][match(allk, names(PROF[[m]]))]), 2)
    print(tab, row.names = FALSE)
    cat("\n2. PEAKS and FLATNESS (reading rule d first).\n\n")
    for (m in MODELS) {
        v <- PROF[[m]]; i <- which.max(v)
        cat(sprintf("   entry %s (npar %d): peak %.2f at p3 = %s h | at the default p3 %s | gain %s | range over %d points (%s-%s h) %.2f\n",
                    m, npar[[m]], v[i], names(v)[i],
                    if (dkey %in% names(v)) sprintf("%.2f", v[[dkey]]) else "NA",
                    if (dkey %in% names(v)) sprintf("%+.2f", v[i] - v[[dkey]]) else "NA",
                    length(v), names(v)[1], names(v)[length(v)], diff(range(v))))
    }
    cat("\n   FLAT (range < 2 units) means p3 is not identified under that start.\n")
    cat(sprintf("\n   entry 59's shared start at its peak p3: (alpha, beta) = %s\n",
                sub("_", ", ", PR[["59"]]$sub[[names(PROF[["59"]])[which.max(PROF[["59"]])]]])))
    cat("\n3. BETWEEN MODELS at each model's own peak, and at a COMMON p3 (reading\n",
        "   rule a and c). Differences in log-likelihood units; parameter difference\n",
        "   in brackets. 57 - 54 was -29.78 (0) with p3 fixed at 18.5802 h.\n\n", sep = "")
    pk <- vapply(MODELS, \(m) max(PROF[[m]]), 0)
    for (pr in list(c("57", "54"), c("59", "54"), c("59", "57")))
        cat(sprintf("   %s - %s at each peak: %+7.2f (%+d parameters)\n", pr[1], pr[2],
                    pk[[pr[1]]] - pk[[pr[2]]], npar[[pr[1]]] - npar[[pr[2]]]))
    for (m in MODELS) if (dkey %in% names(PROF[[m]]))
        cat(sprintf("   at the DEFAULT p3, %s - 54: %+7.2f (%+d parameters)\n", m,
                    PROF[[m]][[dkey]] - PROF[["54"]][[dkey]], npar[[m]] - npar[["54"]]))
    cat("\n4. PER TRIAL at each model's peak p3 (entry 59 at its peak shared start):\n",
        "   maximum log-likelihood, best found.\n\n", sep = "")
    pt <- vapply(MODELS, \(m) { k <- names(PROF[[m]])[which.max(PROF[[m]])]; sb <- PR[[m]]$sub[[k]]
        vapply(LL[[m]], \(x) x[[k]][[sb]], 0) }, numeric(nT))
    rownames(pt) <- TR
    print(round(pt, 2))
    cat("\n5. CONSISTENCY with start-ml.R at the default p3, both with (alpha, beta)\n",
        "   shared across trials. This profile has extra starts, so it should be equal\n",
        "   or higher; entry 59's shape grid here is 5 x 5 against start-ml's 8 x 8, so\n",
        "   a small shortfall there is the coarser grid, not a disagreement.\n\n", sep = "")
    for (m in MODELS) if (dkey %in% names(PROF[[m]])) {
        ms <- ml_shared(m)
        cat(sprintf("   entry %s: here %.2f | start-ml %.2f | difference %+.2f\n",
                    m, PROF[[m]][[dkey]], ms, PROF[[m]][[dkey]] - ms))
    }
}
