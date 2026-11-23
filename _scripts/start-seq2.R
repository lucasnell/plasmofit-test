## Profile likelihood over BOTH sequestration-schedule parameters: the
## midpoint p3 and the slope p4. 2026-11-23, Lucas ("1-3").
##
## WHY. _scripts/start-seq.R freed the schedule's POSITION only (p3 over
## 12-26 h, p4 fixed at the inherited 0.2242) and found that it recovers just
## 5.1 of the ring-only start's 29.8-unit deficit (claude/findings.md,
## "Profile of the sequestration midpoint"). The slope sets how sharply
## parasites leave circulation, which changes the oscillation's amplitude
## rather than its position, so it is a different degree of freedom and it
## has never been free either. This asks whether the gap that survived
## freeing the midpoint is still schedule-shaped.
##
## y(age) = (1 - p1) / (1 + 10^(-p4 (p3 - age))), p1 = 11.3869/467.6209
## inherited and still fixed. Grid: p3 in {15, 17, 18.5802, 20, 22} h and p4
## in {0.1121, 0.1682, 0.2242, 0.3363, 0.4484}, i.e. half, three quarters,
## one, one and a half, and twice the inherited slope. Both are SHARED by
## every trial and unit, so the profile sums the 13 per-trial maxima at each
## of the 25 grid points. Everything else is as start-seq.R: priors ignored,
## each trial x cohort sd profiled, the same bounds as the Stan programs, the
## shared core in _scripts/start-seq-core.R.
##
## ENTRY 59 IS CONDITIONED on a 3 x 3 shared start grid, (alpha, beta) in
## {1, 1.25, 1.5} x {1.25, 1.75, 2.25}, rather than start-seq.R's 5 x 5. It
## contains both the 1-D profile's best shape (1, 1.5) and start-ml.R's
## (1, 1.75), and it is what makes entry 59 affordable over a 2-D grid. The
## cost of that conditioning is measured at the reference point, below.
##
## READING RULE, fixed before any result (2026-11-23). "Peak" = the best grid
## point of the 25; the 1-D figures it is read against are 54 -744.25,
## 57 -768.93, 59 -765.81, with 59 - 54 = -21.55 at 12 fewer parameters.
##   f. If the 59 - 54 gap at the 2-D peaks differs from -21.55 by less than
##      5 units, the slope adds nothing and what remains of the deficit is
##      not schedule-shaped. That is the headline either way.
##   g. If freeing the slope moves either model's preferred MIDPOINT by more
##      than 3 h from its 1-D value, the 1-D midpoints were confounded with
##      the fixed slope and must not be quoted on their own again.
##   h. FLATNESS, checked first, as in start-seq.R: if a model's profile
##      varies by less than 2 units across the slope axis at its best
##      midpoint, the slope is not identified under that start.
##   i. If entry 59's value at the reference point (p3 = 18.5802, p4 =
##      0.2242) is more than 2 units below start-seq.R's -766.29 there, the
##      3 x 3 shape grid is costing too much and entry 59's numbers are
##      reported as conditional, not compared across models.
##   j. As before: no priors, no hierarchy, and a best (p3, p4) here is a
##      diagnostic, never a reportable sequestration schedule.
##
## TASKS (SLURM_ARRAY_TASK_ID 1-91), split so no task approaches the wall:
##   1-13   entry 54, one trial each, all 25 grid points
##   14-26  entry 57, one trial each, all 25 grid points
##   27-91  entry 59, one (trial, p3) each, all 5 slopes
##
## MODES (SSQ2_MODE): check / run / merge.
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   SSQ2_MODE=check Rscript --vanilla _scripts/start-seq2.R > _data/start-seq2-check-<date>.txt
##   sbatch _scripts/start-seq2.sh
##   SSQ2_MODE=merge Rscript --vanilla _scripts/start-seq2.R > _data/start-seq2-<date>.txt

MODE <- Sys.getenv("SSQ2_MODE", "check")
OUTDIR <- "_data/start-seq2/"
P3_GRID <- c(15, 17, 18.5802, 20, 22)
P4_GRID <- c(0.1121, 0.1682, 0.2242, 0.3363, 0.4484)
GRID_A <- c(1, 1.25, 1.5); GRID_B <- c(1.25, 1.75, 2.25)
N_RAND <- 12; N_RAND_G <- 4
source("_scripts/start-seq-core.R")
if (nzchar(Sys.getenv("SSQ2_SMOKE"))) {
    OUTDIR <- paste0(Sys.getenv("SSQ2_SMOKE"), "/"); N_RAND <- 2; N_RAND_G <- 2
    P3_GRID <- c(18.5802); P4_GRID <- c(0.1682, 0.2242); GRID_A <- c(1, 1.5); GRID_B <- 1.75
}
P4_DEF <- 0.2242
kk <- function(p3, p4) sprintf("%g_%g", p3, p4)
REF <- kk(P3_DEF, P4_DEF)
## start-seq.R's 1-D profile, for reading rules f, g and i
ONED <- c(`54` = -744.25, `57` = -768.93, `59` = -765.81)
ONED_P3 <- c(`54` = 12, `57` = 21, `59` = 18)
ONED_REF59 <- -766.29

## ---- check ----------------------------------------------------------------
if (MODE == "check") {
    cat("1. The shared core reproduces the package's forward map at the inherited\n",
        "   schedule, and the slope argument is wired through. Cells: max relative\n",
        "   |difference| in y_hat over t = 72-216 h by 12 h, cycle 42 h, R 12.\n\n", sep = "")
    ts <- seq(72, 216, by = 12); b <- bounds(ts, 0); cl <- 42
    W <- plasmofit:::windowed_beta_starts_below(n_c, 0, WINDOW / cl, 1.2, 2.5, 10)
    ev <- function(ly) plasmofit:::ew_poly_eval(make_DC_free(ly, W, "win"), n_c, cl, 12, 0,
              length(ts), ts, M, R_MAX, LG, RL1, b$rl, b$rh, b$ml, b$mh)
    g <- plasmofit:::ew_poly_series_gen(cl, n_c, 0, WINDOW, 1.2, 2.5, 0, 1, 12, 0, length(ts), ts,
                                        M, R_MAX, LG, RL1, b$rl, b$rh, b$ml, b$mh)
    cat(sprintf("   p3 default, p4 default (must match the package): %.2e\n",
                max(abs(ev(seq_log_y(cl, P3_DEF, P4_DEF)) - g) / g)))
    cat(sprintf("   p4 omitted vs p4 given explicitly (must be 0): %.2e\n",
                max(abs(ev(seq_log_y(cl, P3_DEF)) - ev(seq_log_y(cl, P3_DEF, P4_DEF))))))
    cat("\n2. The slope changes the curve, and in a different way from the midpoint.\n",
        "   Cells: age at which the circulating fraction falls to 0.25 and to 0.75 of\n",
        "   its maximum, and the duty cycle (mean over one 42 h cycle).\n\n", sep = "")
    for (p4 in P4_GRID) {
        y <- exp(seq_log_y(42, P3_DEF, p4)); age <- seq_len(n_c) * 42 / n_c; y <- y / max(y)
        cat(sprintf("   p4 %.4f: 75%% at %5.2f h, 25%% at %5.2f h (width %4.2f h), duty cycle %.4f\n",
                    p4, approx(rev(y), rev(age), 0.75)$y, approx(rev(y), rev(age), 0.25)$y,
                    approx(rev(y), rev(age), 0.25)$y - approx(rev(y), rev(age), 0.75)$y, mean(y * max(exp(seq_log_y(42, P3_DEF, p4))))))
    }
    cat("\n3. At the reference point this must reproduce start-ml.R's per-trial maxima,\n",
        "   as start-seq.R's check does. Cells: the two log-likelihoods and the difference.\n\n", sep = "")
    for (m in MODELS) {
        dr <- load_draws(m)
        for (cc in c(1, 13)) {
            td <- trial_data(cc, DMAX[[m]]); nu <- length(td$units)
            x <- readRDS(sprintf("%stask-%s-%02d.rds", MLDIR, m, cc))
            kb <- names(x$runs)[which.max(vapply(names(x$runs), \(k) max(vapply(x$runs[[k]], `[[`, 0, "ll")), 0))]
            r <- x$runs[[kb]][[which.max(vapply(x$runs[[kb]], `[[`, 0, "ll"))]]
            ab <- if (m == "59") as.numeric(strsplit(kb, "_")[[1]]) else NULL
            p <- unpack(r$par, m, nu, ab)
            v <- trial_ll(td, p$cl, p$R, p$lt0, p$st, P3_DEF, p4 = P4_DEF)
            cat(sprintf("   entry %s trial %-13s: here %10.4f | start-ml %10.4f | difference %.2e\n",
                        m, lv$trial[cc], v, r$ll, v - r$ll))
        }
    }
    quit(save = "no")
}

## ---- run ------------------------------------------------------------------
if (MODE == "run") {
    task <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", NA))
    nT <- d$n_grp_cl; nP <- length(P3_GRID)
    stopifnot(!is.na(task), task >= 1, task <= 2 * nT + nT * nP)
    dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
    if (task <= 2 * nT) {                       # entries 54 and 57: whole grid
        m <- MODELS[(task - 1) %/% nT + 1]
        cc <- (task - 1) %% nT + 1
        p3s <- P3_GRID
        OUT <- sprintf("%stask-%s-%02d.rds", OUTDIR, m, cc)
    } else {                                    # entry 59: one midpoint each
        m <- "59"
        j <- task - 2 * nT - 1
        cc <- j %/% nP + 1
        p3i <- j %% nP + 1
        p3s <- P3_GRID[p3i]
        OUT <- sprintf("%stask-59-%02d-p3%d.rds", OUTDIR, cc, p3i)
    }
    if (file.exists(OUT)) stop(OUT, " exists; refusing to overwrite")
    set.seed(20261123 + task)
    td <- trial_data(cc, DMAX[[m]]); nu <- length(td$units)
    cat(sprintf("task %d: entry %s, trial %s, %d units, %d observations, p3 %s, %d slopes\n",
                task, m, lv$trial[cc], nu, length(td$k),
                paste(sprintf("%g", p3s), collapse = ","), length(P4_GRID)))
    dr <- load_draws(m)
    prof <- vapply(seq_along(dr$chain), \(i) { nat <- draw_nat(m, dr, i, td)
        st <- unpack(pack(m, nat$cl, nat$R, nat$lt0, nat$extra), m, nu, nat$ab)
        trial_ll(td, st$cl, st$R, st$lt0, st$st, P3_DEF, p4 = P4_DEF) }, 0)
    seed_th <- lapply(1:4, \(ch) { i <- which(dr$chain == ch)[which.max(prof[dr$chain == ch])]
        nat <- draw_nat(m, dr, i, td); pack(m, nat$cl, nat$R, nat$lt0, nat$extra) })
    mlb <- ml_best(m, cc)
    solve_at <- function(ab, starts, p3, p4) {
        fn <- function(th) { pp <- unpack(th, m, nu, ab)
            v <- -trial_ll(td, pp$cl, pp$R, pp$lt0, pp$st, p3, p4 = p4)
            if (is.finite(v)) v else 1e10 }
        lapply(starts, \(st) { o <- opt_one(st$th, fn); pp <- unpack(o$par, m, nu, ab)
                               list(src = st$src, ll = -o$value, par = o$par, cl = pp$cl) })
    }
    runs <- list(); prev <- NULL
    for (p3 in p3s) for (p4 in P4_GRID) {
        key <- sprintf("%g_%g", p3, p4)
        if (m != "59") {
            st <- c(lapply(seq_len(N_RAND), \(i) list(src = "random", th = rand_theta(m, nu, i, N_RAND))),
                    lapply(1:4, \(ch) list(src = sprintf("chain%d", ch), th = seed_th[[ch]])),
                    if (!is.null(mlb)) list(list(src = "startml", th = mlb$th)),
                    if (!is.null(prev)) list(list(src = "warm", th = prev)))
            r <- solve_at(NULL, st, p3, p4)
            prev <- r[[which.max(vapply(r, `[[`, 0, "ll"))]]$par
            runs[[key]] <- list(none = r)
            cat(sprintf("   p3 %7g p4 %7.4f: best %.3f\n", p3, p4, max(vapply(r, `[[`, 0, "ll")))); flush.console()
        } else {
            sub <- list(); prev_ab <- NULL; bll <- -Inf; bp <- NULL
            for (a in GRID_A) for (b in GRID_B) {
                st <- c(lapply(seq_len(N_RAND_G), \(i) list(src = "random", th = rand_theta(m, nu, i, N_RAND_G))),
                        lapply(1:2, \(ch) list(src = sprintf("chain%d", ch), th = seed_th[[ch]])),
                        if (!is.null(prev_ab)) list(list(src = "warm_ab", th = prev_ab)),
                        if (!is.null(prev)) list(list(src = "warm", th = prev)),
                        if (!is.null(mlb)) list(list(src = "startml", th = mlb$th)))
                r <- solve_at(c(a, b), st, p3, p4)
                ib <- which.max(vapply(r, `[[`, 0, "ll")); prev_ab <- r[[ib]]$par
                if (r[[ib]]$ll > bll) { bll <- r[[ib]]$ll; bp <- r[[ib]]$par }
                sub[[sprintf("%g_%g", a, b)]] <- r
            }
            prev <- bp
            runs[[key]] <- sub
            cat(sprintf("   p3 %7g p4 %7.4f: best %.3f\n", p3, p4, bll)); flush.console()
        }
    }
    saveRDS(list(model = m, trial = cc, units = td$units, p3 = p3s, p4 = P4_GRID, runs = runs), OUT)
    cat("written", OUT, "\n")
    quit(save = "no")
}

## ---- merge ----------------------------------------------------------------
if (MODE == "merge") {
    nT <- d$n_grp_cl; TR <- lv$trial; nP <- length(P3_GRID)
    best_starts <- function(r) max(vapply(r, `[[`, 0, "ll"))
    runs_of <- function(m, cc) {
        if (m != "59") return(lapply(readRDS(sprintf("%stask-%s-%02d.rds", OUTDIR, m, cc))$runs,
                                     \(sub) vapply(sub, best_starts, 0)))
        out <- list()
        for (i in seq_len(nP)) { f <- sprintf("%stask-59-%02d-p3%d.rds", OUTDIR, cc, i)
            if (!file.exists(f)) next
            x <- readRDS(f)
            for (k in names(x$runs)) out[[k]] <- vapply(x$runs[[k]], best_starts, 0) }
        if (!length(out)) stop("no saved runs for entry 59 trial ", cc)
        out
    }
    LL <- lapply(MODELS, \(m) lapply(seq_len(nT), \(cc) runs_of(m, cc))); names(LL) <- MODELS
    prof_of <- function(m) {
        keys <- Reduce(intersect, lapply(LL[[m]], names))
        v <- vapply(keys, \(k) { subs <- Reduce(intersect, lapply(LL[[m]], \(x) names(x[[k]])))
            max(vapply(subs, \(sb) sum(vapply(LL[[m]], \(x) x[[k]][[sb]], 0)), 0)) }, 0)
        names(v) <- keys; v
    }
    PROF <- lapply(MODELS, prof_of); names(PROF) <- MODELS
    npar <- c(`54` = 95, `57` = 95, `59` = 83)
    cat("1. PROFILE over the sequestration midpoint p3 (rows) and slope p4 (columns).\n",
        "   Cells: maximum log-likelihood over 1130 observations at that (p3, p4),\n",
        "   summed over 13 trials, sd per trial x cohort profiled, priors ignored.\n",
        "   Both are SHARED by every trial and unit; so is entry 59's start shape,\n",
        "   here chosen once for all trials from a 3 x 3 grid. Higher is better; every\n",
        "   value is a lower bound on the true maximum. p4 = 0.2242 and p3 = 18.5802\n",
        "   are the inherited values. NA = grid point missing for some trial.\n\n", sep = "")
    for (m in MODELS) {
        cat(sprintf("   entry %s (npar %d)\n", m, npar[[m]]))
        tb <- outer(P3_GRID, P4_GRID, Vectorize(\(a, b) {
            k <- sprintf("%g_%g", a, b); if (k %in% names(PROF[[m]])) PROF[[m]][[k]] else NA_real_ }))
        dimnames(tb) <- list(sprintf("p3=%g", P3_GRID), sprintf("p4=%g", P4_GRID))
        print(round(tb, 2)); cat("\n")
    }
    cat("2. PEAKS, and reading rules g, h and i.\n\n")
    pk <- list()
    for (m in MODELS) {
        v <- PROF[[m]]; i <- which.max(v); kb <- strsplit(names(v)[i], "_")[[1]]
        pk[[m]] <- list(ll = v[[i]], p3 = as.numeric(kb[1]), p4 = as.numeric(kb[2]))
        rowk <- sprintf("%g_%g", pk[[m]]$p3, P4_GRID); rowk <- rowk[rowk %in% names(v)]
        cat(sprintf("   entry %s: peak %.2f at p3 = %g h, p4 = %g | 1-D peak %.2f at p3 = %g h | gain %+.2f | midpoint moved %+.2f h\n",
                    m, pk[[m]]$ll, pk[[m]]$p3, pk[[m]]$p4, ONED[[m]], ONED_P3[[m]],
                    pk[[m]]$ll - ONED[[m]], pk[[m]]$p3 - ONED_P3[[m]]))
        cat(sprintf("      across the slope axis at that midpoint: range %.2f over %d points -> %s\n",
                    diff(range(v[rowk])), length(rowk),
                    if (diff(range(v[rowk])) < 2) "FLAT, slope not identified" else "slope identified"))
    }
    cat(sprintf("\n   rule i, entry 59's 3 x 3 shape grid at the reference point: %.2f here against %.2f in start-seq.R, cost %.2f units -> %s\n",
                PROF[["59"]][[REF]], ONED_REF59, ONED_REF59 - PROF[["59"]][[REF]],
                if (ONED_REF59 - PROF[["59"]][[REF]] <= 2) "acceptable, comparable across models"
                else "TOO LARGE, entry 59's numbers are conditional"))
    cat("\n3. BETWEEN MODELS at each model's own 2-D peak (reading rule f). Differences\n",
        "   in log-likelihood units, parameter difference in brackets. With the slope\n",
        "   fixed the same comparisons were 57 - 54 -24.67, 59 - 54 -21.55, 59 - 57 +3.12.\n\n", sep = "")
    for (pr in list(c("57", "54"), c("59", "54"), c("59", "57")))
        cat(sprintf("   %s - %s: %+7.2f (%+d parameters)\n", pr[1], pr[2],
                    pk[[pr[1]]]$ll - pk[[pr[2]]]$ll, npar[[pr[1]]] - npar[[pr[2]]]))
    g59 <- pk[["59"]]$ll - pk[["54"]]$ll
    cat(sprintf("\n   rule f: 59 - 54 is %+.2f against -21.55 with the slope fixed, a change of %.2f -> %s\n",
                g59, abs(g59 + 21.55),
                if (abs(g59 + 21.55) < 5) "the slope adds nothing; what remains is NOT schedule-shaped"
                else "the slope matters; the remaining deficit is still schedule-shaped"))
}
