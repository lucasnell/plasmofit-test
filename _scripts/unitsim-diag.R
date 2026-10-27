## Can the data tell the less-biased single-unit fits apart? (2026-10-27,
## Lucas's question, asked after the main study and its extensions were
## read; exploratory, no reading rule was fixed in advance.)
##
## Two paired comparisons on IDENTICAL data, both fits passing the gate:
##   fix - est     b_shape fixed at the truth against estimated (same prior);
##   truth - prod  cycle_length prior centred at the truth against the
##                 production 48 h centre (same sd_logit_cl, so the width in
##                 hours differs too; findings.md correction 2026-10-26).
## The fixed / truth-centred fit is the less biased one on average. Measures:
##   sampler   divergences, treedepth, leapfrog steps, lp__ spread (from the
##             merged summaries). These describe the sampler, not the model.
##   elpd_loo  PSIS-LOO expected log predictive density over the unit's 34 or
##             46 observations, from log_lik recomputed with rstan::gqs on the
##             saved draws (calc_log_lik = 1; the fits were run with 0).
##             Higher = predicts held-out observations better. Paired
##             difference per replicate, with its se from the pointwise
##             differences (sd * sqrt(n), as loo_compare).
##   log BF    log Bayes factor (marginal likelihood ratio), positive = the
##             data favour the less-biased model.
##             truth : prod prior. Only the normal prior on logit_cl differs,
##               so BF = E over the production posterior of
##               N(logit_cl | m_truth, sd) / N(logit_cl | m_prod, sd);
##               estimated again from the truth-prior posterior as the
##               reciprocal, as a check.
##             fix : est (Savage-Dickey; the fixed model is the estimated one
##               at b_shape = truth, and the other priors are independent):
##               posterior density of log b_shape at log(truth) / its prior
##               density, the prior being lognormal(2, 0.5) truncated to
##               [2, max_shape]. Posterior density by kernel estimate.
##
## MODES (UNITSIM_DIAG):
##   check    log_lik from gqs against an independent R computation of one
##            draw (generate_starts + mat_exp_series, as unitsim.R simulates).
##   compute  one diag-<fit>.rds per merged study fit in _data/unitsim/study/
##            (skips existing; UNITSIM_DIAG_MATCH, a regex, restricts it).
##   read     the paired tables; write to _data/unitsim-read-diag-<date>.txt.
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   UNITSIM_DIAG=check   Rscript --vanilla _scripts/unitsim-diag.R
##   UNITSIM_DIAG=compute Rscript --vanilla _scripts/unitsim-diag.R
##   UNITSIM_DIAG=read    Rscript --vanilla _scripts/unitsim-diag.R > _data/unitsim-read-diag-<date>.txt

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
suppressPackageStartupMessages({
    library(plasmofit); library(rstan); library(loo)
    library(dplyr); library(tidyr)
})
options(width = 150)
DIR <- "_data/unitsim/study/"
MODE <- Sys.getenv("UNITSIM_DIAG", "read")
mod <- plasmofit:::stanmodels[["archer_fit_pooled_both"]]
diag_path <- function(nm) sprintf("%sdiag-%s.rds", DIR, nm)

log_lik_of <- function(nm) {
    f <- readRDS(sprintf("%sfit-%s.rds", DIR, nm))
    d <- readRDS(sprintf("%schain-%s-c1.rds", DIR, nm))$d
    d$calc_log_lik <- 1L
    dr <- as.matrix(f)
    g <- rstan::gqs(mod, data = d, draws = dr)
    ll <- rstan::extract(g, "log_lik", permuted = FALSE)    # iter x chain x obs
    list(f = f, d = d, dr = dr, ll = ll)
}

if (MODE == "check") {
    nm <- "w168-g384-est-rep1"
    x <- log_lik_of(nm)
    i <- 1L; p <- x$dr[i, ]; d <- x$d
    cl <- p[["cycle_length[1]"]]
    y0 <- plasmofit:::generate_starts(cl, d$n_c, p[["b_shape[1]"]], p[["b_offset[1]"]],
                                      p[["log10_total0[1]"]])
    u_ts <- sort(unique(d$ts))
    yh <- plasmofit:::mat_exp_series(y0, cl, d$n_c, p[["R[1]"]], d$mu, u_ts, d$dt_full)
    ll_r <- dnorm(log10(d$y + 1), log10(yh[match(d$ts, u_ts)] + 1), p[["sd_iRBC[1]"]], log = TRUE)
    ll_s <- x$ll[1, 1, ]
    cat(sprintf("%s draw 1: sum log_lik gqs %.6f, R %.6f; max |pointwise diff| %.2e over %d obs\n",
                nm, sum(ll_s), sum(ll_r), max(abs(ll_s - ll_r)), length(ll_r)))
    stopifnot(max(abs(ll_s - ll_r)) < 1e-6)
    quit(save = "no")
}

if (MODE == "compute") {
    sums <- Sys.glob(paste0(DIR, "sum-*.rds"))
    nms <- sub("^.*sum-(.*)\\.rds$", "\\1", sums)
    nms <- nms[!grepl("-pwide$", nms) | file.exists(sprintf("%sfit-%s.rds", DIR, nms))]
    ## UNITSIM_DIAG_MATCH (regex) restricts compute to matching fits, e.g. the
    ## fixed-pin grid (added 2026-10-27).
    mt <- Sys.getenv("UNITSIM_DIAG_MATCH", "")
    if (nzchar(mt)) nms <- nms[grepl(mt, nms)]
    n_done <- 0L
    for (nm in nms) {
        if (file.exists(diag_path(nm))) next
        x <- log_lik_of(nm)
        r_eff <- relative_eff(exp(x$ll))
        lo <- suppressWarnings(loo(x$ll, r_eff = r_eff))
        k <- lo$diagnostics$pareto_k
        ll_tot <- apply(x$ll, c(1, 2), sum)
        out <- list(name = nm, n_obs = dim(x$ll)[3],
                    elpd_loo = lo$estimates["elpd_loo", "Estimate"],
                    p_loo = lo$estimates["p_loo", "Estimate"],
                    k_max = max(k), k_bad = sum(k > 0.7),
                    pointwise = lo$pointwise[, "elpd_loo"],
                    ll_mean = mean(ll_tot),
                    logit_cl = x$dr[, "logit_cl"],
                    log_b = if ("b_shape_free[1]" %in% colnames(x$dr)) log(x$dr[, "b_shape_free[1]"]) else NULL,
                    mean_logit_cl = x$d$mean_logit_cl, sd_logit_cl = x$d$sd_logit_cl,
                    mean_log_b = x$d$mean_log_b_shape[1], sd_log_b = x$d$sd_log_b_shape[1],
                    max_shape = x$d$max_shape)
        saveRDS(out, diag_path(nm))
        n_done <- n_done + 1L
        if (n_done %% 50L == 0L) cat(format(Sys.time()), n_done, "done\n")
    }
    cat("computed", n_done, "this run;", length(Sys.glob(paste0(DIR, "diag-*.rds"))), "in all\n")
    quit(save = "no")
}

if (MODE != "read") stop("unknown UNITSIM_DIAG: ", MODE)

## ---- read ------------------------------------------------------------------
s <- bind_rows(lapply(Sys.glob(paste0(DIR, "sum-*.rds")), readRDS)) |>
    mutate(prior = coalesce(prior, "prod"), truth_set = coalesce(truth_set, "main"),
           lf_iter = vapply(per_chain, \(p) mean(p$leapfrog_per_iter), numeric(1)))
## The fixed-pin grid is read by _scripts/unitsim-pins-read.R.
if ("fit_nc" %in% names(s)) s <- filter(s, coalesce(fit_nc, 384L) == 384L)
if ("pin" %in% names(s)) s <- filter(s, coalesce(pin, "truth") == "truth")
dg <- lapply(Sys.glob(paste0(DIR, "diag-*.rds")), readRDS)
names(dg) <- vapply(dg, `[[`, character(1), "name")
s <- filter(s, name %in% names(dg)) |>
    mutate(elpd = vapply(dg[name], `[[`, numeric(1), "elpd_loo"),
           k_max = vapply(dg[name], `[[`, numeric(1), "k_max"),
           k_bad = vapply(dg[name], `[[`, numeric(1), "k_bad"),
           p_loo = vapply(dg[name], `[[`, numeric(1), "p_loo"))
cat(nrow(s), "fits with diagnostics;", sum(s$pass), "passing. Exploratory: no reading rule fixed in advance.\n")
cat(sprintf("Pareto k > 0.7 in %d of %d fits (max k over all fits %.2f).\n\n",
            sum(s$k_bad > 0), nrow(s), max(s$k_max)))

## Savage-Dickey log BF, fixed at truth : estimated, from an est fit.
sd_logbf <- function(g, b0) {
    lb <- g$log_b
    dens <- density(lb, n = 2048)
    post <- approx(dens$x, dens$y, xout = log(b0))$y
    if (is.na(post)) post <- 0
    Z <- pnorm(log(g$max_shape), g$mean_log_b, g$sd_log_b) - pnorm(log(2), g$mean_log_b, g$sd_log_b)
    prior <- dnorm(log(b0), g$mean_log_b, g$sd_log_b) / Z
    c(log(post) - log(prior), mean(lb > log(b0)))
}
## log BF of prior A : prior B (logit_cl normals, same sd) from B's posterior.
is_logbf <- function(g, mA) {
    lw <- dnorm(g$logit_cl, mA, g$sd_logit_cl, log = TRUE) -
          dnorm(g$logit_cl, g$mean_logit_cl, g$sd_logit_cl, log = TRUE)
    m <- max(lw); m + log(mean(exp(lw - m)))
}

pair_tab <- function(a, b, by, la, lb) {
    j <- inner_join(a, b, by = by, suffix = c("_a", "_b"))
    j$d_se <- as.numeric(mapply(\(x, y) sd(dg[[x]]$pointwise - dg[[y]]$pointwise) *
                         sqrt(length(dg[[x]]$pointwise)), j$name_a, j$name_b))
    j |> mutate(d_elpd = elpd_a - elpd_b, d_absbias = abs(bias_a) - abs(bias_b),
                d_bias = bias_a - bias_b)
}
summ <- function(j, grp) j |> group_by(across(all_of(grp))) |>
    summarise(pairs = n(),
              less_biased = sum(d_absbias < 0),
              d_elpd_mean = mean(d_elpd), d_elpd_se = sd(d_elpd) / sqrt(n()),
              favours = sum(d_elpd > 0), clear = sum(d_elpd > 2 * d_se),
              clear_against = sum(d_elpd < -2 * d_se),
              pair_se_med = median(d_se),
              agree = sum((d_elpd > 0) == (d_absbias < 0)),
              cor_elpd_absbias = cor(d_elpd, d_absbias),
              .groups = "drop") |>
    mutate(across(where(is.double), \(v) round(v, 3)))

cat("1. SAMPLER DIAGNOSTICS BY ARM. Cells: mean over passing fits of % divergent\n",
    "post-warmup iterations, % at max treedepth (12), leapfrog steps per iteration,\n",
    "lp__ spread across the 4 chain means, and posterior sd of cycle_length (h);\n",
    "mean bias (h) for reference. Unpaired. lp__ is not comparable between est\n",
    "and fix (different parameters) or between priors (different densities).\n\n", sep = "")
s |> filter(pass) |> group_by(truth_set, prior, window, gen_nc, method) |>
    summarise(n = n(), div_pct = mean(div_pct), td_cap_pct = mean(treedepth_cap_pct),
              lf_iter = mean(lf_iter), lp_spread = mean(lp_spread), cl_sd = mean(cl_sd),
              bias = mean(bias), .groups = "drop") |>
    mutate(across(where(is.double), \(v) round(v, 3))) |> as.data.frame() |> print(row.names = FALSE)

key <- c("truth_set", "prior", "window", "gen_nc", "rep")
sp <- filter(s, pass)
cat("\n2. FIX - EST (b_shape fixed at the truth minus estimated), identical data,\n",
    "both passing. less_biased = pairs where |bias| is smaller for fix. d_elpd =\n",
    "elpd_loo(fix) - elpd_loo(est), log predictive density units summed over the\n",
    "unit's observations, positive = fix predicts better: mean and se over pairs;\n",
    "favours = pairs with d_elpd > 0; clear / clear_against = pairs beyond +-2 se\n",
    "of their own pointwise difference (pair_se_med = median of that se); agree =\n",
    "pairs where the sign of d_elpd picks the less-biased fit; cor = correlation\n",
    "over pairs of d_elpd with the difference in |bias| (negative = elpd tracks bias).\n\n", sep = "")
j_fe <- pair_tab(filter(sp, method == "fix"), filter(sp, method == "est"),
                 key, "fix", "est")
print(as.data.frame(summ(j_fe, c("truth_set", "prior", "window", "gen_nc"))), row.names = FALSE)

cat("\n   Savage-Dickey log BF, fixed at truth : estimated, from each passing est\n",
    "fit (positive = the data favour fixing at the truth). Cells: mean, se, and\n",
    "count > 0 over fits; post_above = mean share of posterior b_shape draws above\n",
    "the truth. Kernel density on log b_shape, default bandwidth.\n\n", sep = "")
true_b <- function(nm) readRDS(sprintf("%schain-%s-c1.rds", DIR, nm))$truth$b_shape
est <- filter(sp, method == "est")
sdv <- t(vapply(est$name, \(nm) sd_logbf(dg[[nm]], true_b(nm)), numeric(2)))
est |> mutate(lbf = sdv[, 1], above = sdv[, 2]) |>
    group_by(truth_set, prior, window, gen_nc) |>
    summarise(n = n(), logbf_mean = mean(lbf), logbf_se = sd(lbf) / sqrt(n()),
              positive = sum(lbf > 0), post_above = mean(above), .groups = "drop") |>
    mutate(across(where(is.double), \(v) round(v, 3))) |> as.data.frame() |> print(row.names = FALSE)

if (any(sp$prior == "truth")) {
    key2 <- c("truth_set", "method", "window", "gen_nc", "rep")
    j_tp <- pair_tab(filter(sp, prior == "truth"), filter(sp, prior == "prod"),
                     key2, "truth", "prod")
    cat("\n3. TRUTH-CENTRED - PRODUCTION PRIOR, identical data, both passing. Columns\n",
        "as in 2, with 'truth-centred' in place of 'fix'.\n\n", sep = "")
    print(as.data.frame(summ(j_tp, c("truth_set", "method", "window", "gen_nc"))), row.names = FALSE)
    j_tp$lbf_from_prod <- as.numeric(mapply(\(a, b) is_logbf(dg[[b]], dg[[a]]$mean_logit_cl),
                                 j_tp$name_a, j_tp$name_b))
    j_tp$lbf_from_truth <- as.numeric(mapply(\(a, b) -is_logbf(dg[[a]], dg[[b]]$mean_logit_cl),
                                  j_tp$name_a, j_tp$name_b))
    cat("\n   log BF, truth-centred : production prior (positive = the data favour the\n",
        "truth-centred prior), importance-weighted from the production posterior;\n",
        "check = the same from the truth-centred posterior. Cells: mean, se, count > 0,\n",
        "and the largest |difference| between the two estimates.\n\n", sep = "")
    j_tp |> group_by(truth_set, method, window, gen_nc) |>
        summarise(pairs = n(), logbf_mean = mean(lbf_from_prod),
                  logbf_se = sd(lbf_from_prod) / sqrt(n()), positive = sum(lbf_from_prod > 0),
                  check_mean = mean(lbf_from_truth),
                  max_abs_diff = max(abs(lbf_from_prod - lbf_from_truth)), .groups = "drop") |>
        mutate(across(where(is.double), \(v) round(v, 3))) |> as.data.frame() |> print(row.names = FALSE)
}
