# Single-unit simulation study of cycle_length bias. Set up 2026-10-21 with
# Lucas, who chose it before the full-model study (_scripts/simstudy2.R, which
# is NOT to run unless he says so): if the concern is structural bias in the
# model, start with the smallest problem that has it.
#
# THE PROBLEM. One grp_init unit (trial x inoculum size). Inside a unit every
#   series shares every parameter, so it is one trajectory with its series as
#   replicates. Its series are simulated and fitted alone, so no other unit and
#   no hierarchy is involved: the bias is that of the likelihood, the nuisance
#   parameters, and the priors.
#
# WINDOW (Lucas, 2026-10-22, replacing a two-real-unit design that confounded
#   window length with series count, start times, cohorts, and truth): ONE
#   unit, MMV048_PIB.1800 (u7), ONE truth, two observation schedules.
#     w168  its six real series at their real times (96 or 120 to 168 h by
#           12 h, no 108; 34 observations);
#     w192  the same six series and times plus 180 and 192 h on each (46
#           observations). Only the end of the window differs (and with it,
#           12 more observations): in the real data every series ends at 168
#           or 192 h (Mefloquine.1800's one 216 h series aside).
#
# TRUTH: the posterior means of the real-data fit
#   np_wide_total0_nc384_t (wockner-fit.R entry 54: no_pool, n_c 384, b_shape
#   estimated under lognormal(2, 0.5), sd_log10_total0 1, adapt_delta 0.95,
#   max_treedepth 12; 41.815 h mean over 13 trials) for that unit's b_shape,
#   b_offset, log10_total0, its trial's R and cycle_length, and its cohorts'
#   sd_iRBC. At gen n_c 1024, b_shape is lowered so the total stage spread at
#   the unit's REAL last observation (168 h) equals the truth's at 384 (same
#   mapping and reason as simstudy2.R: 1 / (4 (2b' + 1)) = 1 / (4 (2b + 1)) +
#   k (1/384 - 1/1024), k = 168 / true cycle_length). The truth is the same
#   under both windows.
#
# SIMULATED DATA. The noise is drawn for all 1130 real observations with
#   set.seed(REP_SEEDS[rep]) exactly as simstudy2.R draws it, and the unit's
#   34 are taken from it, so w168's dataset is the same numbers as the unit's
#   part of the full-model dataset of the same (gen n_c, rep). The 12 added
#   observations of w192 take the next 12 draws of the same stream (series in
#   order of first appearance, then 180 before 192), so w192 = w168's data
#   plus 12 observations. Same seed at both windows, both gen n_c, and both
#   methods (paired on identical noise).
#
# FIT. model = "pooled_both" on the unit's data: one R and one cycle_length
#   shared by all series, with the production top-level priors placed
#   directly on them (cycle_length logit-normal centred 48 h, sd_logit_cl 1,
#   bounds 35-50 h; logit_R normal(-2, 1)); b_shape lognormal(2, 0.5) when
#   estimated (method "est") or fixed at the truth (method "fix");
#   log10_total0 normal(1, 1); n_c 384 always; adapt_delta 0.95,
#   max_treedepth 12; 4 chains; sampler seed 538065874, chain_id k.
#   Iterations per chain by method (Lucas, 2026-10-22, before any study fit):
#   est 700 warmup + 1000 draws (production); fix 2000 warmup + 2000 draws.
#   The probe (SLURM 32412) passed fixed b_shape in 3 of 4 fits at production
#   length, the failure marginal (max R-hat 1.170, chains' mean lp__ within
#   0.4). Sampler settings change exploration, not the posterior, so the two
#   arms are comparable; the longer setting applies to EVERY fixed fit, and
#   no failed fit is retried.
#   archer_stan_data() refuses grp_R or grp_cl with one level (it guards
#   no_pool's hierarchy). pooled_both broadcasts its single R and
#   cycle_length to every group, so the series are split into two labels
#   (alternating in order of first appearance) used for grp_R and grp_cl: the
#   likelihood is unchanged and the trajectory is computed twice per gradient
#   instead of once, on two threads. check mode verifies log_prob is the same
#   under a different split.
#
# DESIGN, fixed 2026-10-22 with Lucas before any study fit ran:
#   window   w168, w192 (above), one unit and one truth
#   gen_nc   384, 1024
#   method   est (b_shape estimated), fix (b_shape fixed at the truth)
#   rep      1-60 (REP_SEEDS; 1-16 shared with simstudy2.R, 17-60 drawn
#            2026-10-22 with set.seed(20261022)). n fixed in advance.
#   = 8 cells x 60 = 480 fits, 1920 chains, in _data/unitsim/study/.
#   Replicate 1 is refitted although the probes fitted some of its cells:
#   every study fit then comes from one frozen script and one setting.
#   Probes (2026-10-21, SLURM 32409, 32412) wrote to _data/unitsim/ and are
#   not part of the study.
#
# READING RULE, fixed 2026-10-22 before any study fit ran
#   (_scripts/unitsim-read.R):
#   gate      max R-hat < 1.05 over every monitored quantity (NA excluded),
#             4 chains. A failed fit is not used for bias and not retried.
#   bias      per fit: posterior mean of cycle_length - true cycle_length,
#             hours; positive = long.
#   per cell  n passing of 60 (an outcome in itself); mean and se (sd /
#             sqrt(n)) of bias over passing fits; READABLE if >= 45 of 60
#             pass. Display only: mean over all fits, and passing against
#             failing fits.
#   coverage  share of passing fits whose central 80% interval holds the
#             truth.
#   b_shape   est arm: log(posterior mean / true b_shape), mean over fits.
#   contrasts per rep, both fits passing:
#             w192 - w168 at each (gen_nc, method): same truth, same noise on
#               the 34 shared observations; the effect of the longer window;
#             est - fix (cost of estimating b_shape) at each (window, gen_nc):
#               identical data;
#             gen 1024 - 384 at each (window, method): same noise draws,
#               different noiseless trajectories, so partially paired.
#   failures  marginal 1.05-1.15, severe >= 1.15; ONE CHAIN APART if a
#             chain's mean lp__ is more than 3 units from the median of the
#             other three.
#
# PRIOR ARM, added 2026-10-23 (Lucas) after the main study's read, before any
#   prior-arm fit ran: the same 480 (window, gen_nc, method, rep) datasets
#   refitted with the cycle_length prior CENTRED AT THE TRUTH
#   (cl_prior_center = 42.844 h, so mean_logit_cl = logit((42.844 - 35) / 15);
#   sd_logit_cl 1 as production). Tokens end "-ptruth". The main study found
#   the w168 / estimated posterior close to the production prior (centred
#   48 h); this arm measures the bias when the prior's location is right.
#   Reading rule as above per cell, plus the paired contrast truth-centred -
#   production prior at each (window, gen_nc, method), identical data, both
#   passing. `UNITSIM_TASKS=prior` in tasks mode lists its 1920 tokens.
#
# LOW-SYNCHRONY TRUTH, added 2026-10-24 (Lucas) before any of its fits ran:
#   the main study's truth with ONLY b_shape changed, to 5.86, the lowest
#   posterior mean of the 14 units in the same real-data fit
#   (MMV048_PartB.2800); mapped to ~4.30 at gen n_c 1024 as above (k = 168 h /
#   true cycle_length). Same schedules, production prior, all 8 cells, 60
#   reps, same noise seeds (so paired on noise with the main study, on a
#   different noiseless trajectory). Tokens carry "-bslow" after the rep;
#   `UNITSIM_TASKS=bslow` lists its 1920. Reading rule as above per cell,
#   plus bslow - main bias at each (window, gen_nc, method), same noise.
#
# LAYOUT: pinned through _scripts/simstudy-pinned.sh. A benchmark of two
#   chains per L3 group (2 physical cores each) against one (4 cores) runs
#   first (2026-10-22; Lucas decides the layout from its numbers), in
#   _data/unitsim/bench-w4/ and bench-w2/; its fits are not part of the study.
#
# TASKS are tokens "w<168|192>-g<gen_nc>-<est|fix>-rep<r>-c<chain>", so file
#   names do not depend on a design grid. (The 2026-10-21/22 probes and the
#   layout benchmark used "u<unit>-..." tokens from earlier versions, frozen
#   as _data/unitsim/unitsim-run{,-2,-3}.R.) The pinned driver passes a token as
#   SLURM_ARRAY_TASK_ID.
#
# MODES (UNITSIM_MODE):
#   chain  (default) one chain for the token in SLURM_ARRAY_TASK_ID.
#   merge  merge every fit whose 4 chain files exist and that has no summary.
#   check  print the two schedules and the truth; for gen 384 and 1024, rep 1,
#          check w168's data against simstudy.R's, w192's shared observations
#          against w168's, and the split's log_prob invariance.
#   tasks  print the study's 1920 chain tokens, rep-major.
#   UNITSIM_ITER overrides iterations in chain mode (smoke test; "smoke-").
#   UNITSIM_DIR overrides the output directory (default _data/unitsim/study/).
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   UNITSIM_MODE=check Rscript --vanilla _scripts/unitsim.R
#   UNITSIM_MODE=tasks Rscript --vanilla _scripts/unitsim.R > _data/unitsim/study/tasks.txt
#   (submit through simstudy-pinned.sh; see TODO.md for the command)
#   UNITSIM_MODE=merge Rscript --vanilla _scripts/unitsim.R

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
suppressPackageStartupMessages({
    library(plasmofit)
    library(readr)
    library(dplyr)
    library(rstan)
})

N_CHAINS <- 4L
SETTINGS <- list(est = list(warmup = 700L, iter = 1700L),
                 fix = list(warmup = 2000L, iter = 4000L))
FIT_NC <- 384L
CONTROL <- list(adapt_delta = 0.95, max_treedepth = 12)
SEED_FIT <- 538065874L
REP_SEEDS <- c(742160183L, 1908445027L, 355721694L, 1414213562L, 1732050808L,
               273469144L, 923590678L, 832351853L, 1481098644L, 996583259L,
               1100299434L, 938347168L, 681577182L, 1986708724L, 1346550232L,
               1529979755L,                  # 1-16 as simstudy2.R
               504979627L, 1539660851L, 1627240572L, 610765287L, 567525470L,
               762245924L, 1968057179L, 992715895L, 1643895353L, 1163691978L,
               1888244126L, 1178995960L, 27758083L, 450201449L, 1211793044L,
               1227917299L, 265085477L, 1429626433L, 1635957329L, 1154687219L,
               942277168L, 205349015L, 1724489978L, 1988255862L, 634536402L,
               1519008358L, 434208822L, 1165502811L, 1554585708L, 756201702L,
               1695726929L, 843984591L, 749845090L, 712884944L, 2123819819L,
               1607165064L, 516892384L, 1034340295L, 647742746L, 1996862138L,
               628141401L, 1693362474L, 1341215164L, 536629132L)
stopifnot(length(REP_SEEDS) == 60L, !anyDuplicated(REP_SEEDS))
UNIT <- 7L                                   # MMV048_PIB.1800
STUDY <- expand.grid(window = c(168L, 192L), gen_nc = c(384L, 1024L),
                     method = c("est", "fix"), rep = seq_along(REP_SEEDS),
                     stringsAsFactors = FALSE)
STUDY$name <- with(STUDY, sprintf("w%d-g%d-%s-rep%d", window, gen_nc, method, rep))
TRUTH <- "np_wide_total0_nc384_t"

MODE <- Sys.getenv("UNITSIM_MODE", "chain")
iter_env <- Sys.getenv("UNITSIM_ITER", "")
SMOKE <- nzchar(iter_env)
DIR <- Sys.getenv("UNITSIM_DIR", "_data/unitsim/study/")
if (!endsWith(DIR, "/")) DIR <- paste0(DIR, "/")
PFX <- if (SMOKE) paste0(DIR, "smoke-") else DIR
chain_path <- function(nm, k) sprintf("%schain-%s-c%d.rds", PFX, nm, k)
fit_path <- function(nm) sprintf("%sfit-%s.rds", PFX, nm)
sum_path <- function(nm) sprintf("%ssum-%s.rds", PFX, nm)

parse_fit <- function(nm) {
    m <- regmatches(nm, regexec("^w(168|192)-g([0-9]+)-(est|fix)-rep([0-9]+)(-bslow)?(-ptruth)?$", nm))[[1]]
    if (length(m) != 7L) stop("bad fit name: ", nm)
    list(window = as.integer(m[2]), unit = UNIT, gen_nc = as.integer(m[3]),
         method = m[4], rep = as.integer(m[5]),
         truth_set = if (nzchar(m[6])) "bslow" else "main",
         prior = if (nzchar(m[7])) "truth" else "prod", name = nm)
}


# ------------------------------------------------------------------------ #
# Full-data structure, truth, and the unit subset
# ------------------------------------------------------------------------ #

paras_df <- read_csv("_data/wockner-cleaned.csv", col_types = "cccdcdd") |>
    mutate(inoc = interaction(trial, inoc_size, drop = TRUE) |> paste()) |>
    mutate(obs_error = interaction(trial, cohort, drop = TRUE) |> paste())

build_full <- function() archer_stan_data(
    data = paras_df, series = "id", time = "time", abundance = "para",
    grp_init = "inoc", grp_R = "trial", grp_cl = "trial", grp_sd = "obs_error",
    calc_log_lik = 0L, mean_log_b_shape = 2, sd_log_b_shape = 0.5,
    mean_log10_total0 = 1, sd_log10_total0 = 1)
d0 <- build_full()
lv0 <- attr(d0, "levels"); gc0 <- attr(d0, "grp_columns")
ARRAY_FIELDS <- names(d0)[lengths(d0) != 1L]
stopifnot(identical(attr(read_rds(sprintf("_data/wock-data-%s.rds", TRUTH)), "levels"), lv0))
UNITS <- lv0[[gc0$grp_init]]
TRIALS <- lv0[[gc0$grp_cl]]
COHORTS <- lv0[[gc0$grp_sd]]

f_truth <- read_rds(sprintf("_data/wock-fit-%s.rds", TRUTH))
pm <- function(p) unname(colMeans(as.matrix(f_truth, pars = p)))
T_FULL <- list(b_shape = pm("b_shape"), b_offset = pm("b_offset"),
               log10_total0 = pm("log10_total0"), R = pm("R"),
               sd_iRBC = pm("sd_iRBC"), cycle_length = pm("cycle_length"))
rm(f_truth)

obs_key0 <- paste(rep(lv0$series, d0$n_obs), d0$ts)
sd_obs0 <- T_FULL$sd_iRBC[rep(d0$grp_sd, d0$n_obs)]

## Noise for all 1130 observations, drawn exactly as simstudy2.R's
## add_noise(): rnorm(n, 0, sd) is sd * the same standard normal draws.
noise_full <- function(rep) { set.seed(REP_SEEDS[rep]); rnorm(d0$n_total_obs, 0, sd_obs0) }

var_init <- function(b) 1 / (4 * (2 * b + 1))
b_from_var <- function(v) (1 / (4 * v) - 1) / 2

## The unit's rows, its Stan data under a given two-label split, and its truth.
unit_rows <- function(u) paras_df |> filter(inoc == UNITS[u])

## The schedule: the unit's real rows (w168), or those plus 180 and 192 h on
## every series (w192). `para` of an added row is a placeholder; y is replaced
## by simulated data before any fit.
window_rows <- function(window) {
    df <- unit_rows(UNIT)
    stopifnot(max(df$time) == 168)
    if (window == 168L) return(df)
    ids <- unique(df$id)
    add <- df[match(rep(ids, each = 2L), df$id), ] |>
        mutate(time = rep(c(180, 192), length(ids)), para = 1)
    bind_rows(df, add) |> arrange(match(id, ids), time)
}
ADDED_KEYS <- { ids <- unique(unit_rows(UNIT)$id)
                paste(rep(ids, each = 2L), rep(c(180, 192), length(ids))) }

build_unit <- function(window, split = "alt", ...) {
    df <- window_rows(window)
    ids <- unique(df$id)
    lab <- if (split == "alt") rep_len(c("A", "B"), length(ids))
           else c(rep("A", ceiling(length(ids) / 2)), rep("B", floor(length(ids) / 2)))
    df$half <- lab[match(df$id, ids)]
    args <- modifyList(list(data = df, series = "id", time = "time",
                            abundance = "para", grp_init = "inoc",
                            grp_R = "half", grp_cl = "half",
                            grp_sd = "obs_error", calc_log_lik = 0L, n_c = FIT_NC,
                            mean_log_b_shape = 2, sd_log_b_shape = 0.5,
                            mean_log10_total0 = 1, sd_log10_total0 = 1),
                       list(...))
    d <- do.call(archer_stan_data, args)
    ## rstan passes a length-1 R vector as a scalar, but per-unit fields are
    ## Stan vectors/arrays (n_grp_init = 1 here). Fields that are not length
    ## 1 in the full data are arrays in Stan; keep their dims.
    for (nm in union(ARRAY_FIELDS, "b_shape_data"))
        if (!is.null(d[[nm]]) && length(d[[nm]]) == 1L) d[[nm]] <- as.array(d[[nm]])
    d
}

## The low-synchrony truth's b_shape: the lowest unit posterior mean.
B_LOW <- T_FULL$b_shape[UNITS == "MMV048_PartB.2800"]
stopifnot(length(B_LOW) == 1L, B_LOW == min(T_FULL$b_shape))

unit_truth <- function(u, gen_nc, b384 = T_FULL$b_shape[u]) {
    tr <- unique(d0$grp_cl[d0$grp_init == u]); stopifnot(length(tr) == 1L)
    df <- unit_rows(u)
    k_last <- max(df$time) / T_FULL$cycle_length[tr]
    b <- b_from_var(var_init(b384) + k_last * (1 / FIT_NC - 1 / gen_nc))
    stopifnot(b > 2, b < d0$max_shape)
    list(unit = UNITS[u], trial = TRIALS[tr], b_shape = b, b_shape_384 = b384,
         k_last = k_last, b_offset = T_FULL$b_offset[u],
         log10_total0 = T_FULL$log10_total0[u], R = T_FULL$R[tr],
         cycle_length = T_FULL$cycle_length[tr],
         sd_iRBC = setNames(T_FULL$sd_iRBC, COHORTS))
}

make_fit_inputs <- function(nm, split = "alt") {
    cf <- parse_fit(nm)
    truth <- if (cf$truth_set == "bslow") unit_truth(cf$unit, cf$gen_nc, B_LOW)
             else unit_truth(cf$unit, cf$gen_nc)
    args <- list(window = cf$window, split = split)
    if (cf$method == "fix") args$b_shape <- truth$b_shape
    if (cf$prior == "truth") args$cl_prior_center <- truth$cycle_length
    d <- do.call(build_unit, args)
    stopifnot(d$n_grp_init == 1L, d$n_c == FIT_NC,
              d$b_shape_fixed == as.integer(cf$method == "fix"))
    lv <- attr(d, "levels")
    y0 <- plasmofit:::generate_starts(truth$cycle_length, cf$gen_nc, truth$b_shape,
                                      truth$b_offset, truth$log10_total0)
    ## Every series in the unit shares every parameter, so they share one
    ## trajectory: evaluate it once at the unit's distinct times (one matrix
    ## exponential, ~22 s at n_c 1024; gotchas.md, "Trajectory cost is roughly
    ## CUBIC"). mat_exp_series steps from t = 0, so a value at time t does not
    ## depend on which other times are requested: identical to per-series calls.
    u_ts <- sort(unique(d$ts))
    y_u <- plasmofit:::mat_exp_series(y0, truth$cycle_length, cf$gen_nc,
                                      truth$R, d$mu, u_ts, d$dt_full)
    y_hat <- y_u[match(d$ts, u_ts)]
    stopifnot(all(is.finite(y_hat)), all(y_hat >= 0))
    ## Noise: real observations from the full noise vector, by (series, time);
    ## added ones from the next draws of the same stream (ADDED_KEYS order).
    key <- paste(rep(lv$series, d$n_obs), d$ts)
    sd_d <- unname(truth$sd_iRBC[lv[[attr(d, "grp_columns")$grp_sd]]][rep(d$grp_sd, d$n_obs)])
    pos <- match(key, obs_key0)
    added <- is.na(pos)
    stopifnot(!anyDuplicated(pos[!added]), all(key[added] %in% ADDED_KEYS),
              sum(added) == if (cf$window == 192L) length(ADDED_KEYS) else 0L,
              isTRUE(all.equal(sd_obs0[pos[!added]], sd_d[!added])))
    set.seed(REP_SEEDS[cf$rep])
    z_full <- rnorm(d0$n_total_obs, 0, sd_obs0)
    z_add <- rnorm(length(ADDED_KEYS))
    eps <- numeric(d$n_total_obs)
    eps[!added] <- z_full[pos[!added]]
    eps[added] <- z_add[match(key[added], ADDED_KEYS)] * sd_d[added]
    d$y <- pmax(0, 10^(log10(y_hat + 1) + eps) - 1)
    list(cf = cf, truth = truth, d = d, y_hat = y_hat, pos = pos, key = key)
}


# ------------------------------------------------------------------------ #
# Modes
# ------------------------------------------------------------------------ #

dir.create(DIR, showWarnings = FALSE)

if (MODE == "check") {
    t384 <- unit_truth(UNIT, 384L); t1024 <- unit_truth(UNIT, 1024L)
    cat(sprintf("Unit %s, truth: cycle_length %.3f h, R %.2f, b_shape %.2f (gen 384) / %.2f (gen 1024, mapped at k = %.3f), log10_total0 %.3f, b_offset %.3f\n",
                UNITS[UNIT], t384$cycle_length, t384$R, t384$b_shape, t1024$b_shape,
                t384$k_last, t384$log10_total0, t384$b_offset))
    for (w in c(168L, 192L)) {
        df <- window_rows(w)
        cat(sprintf("\nw%d: %d series, %d observations\n", w, n_distinct(df$id), nrow(df)))
        df |> group_by(id) |> summarise(times = paste(time, collapse = ","), .groups = "drop") |>
            as.data.frame() |> print(row.names = FALSE)
    }
    old_y <- read_rds("_data/simstudy/chain-B-g384-f384-bs1.0-wide-rep1-c1.rds")$d$y
    mod <- plasmofit:::stanmodels[["archer_fit_pooled_both"]]
    for (g in c(384L, 1024L)) {
        a <- make_fit_inputs(sprintf("w168-g%d-est-rep1", g))
        b <- make_fit_inputs(sprintf("w192-g%d-est-rep1", g))
        cat(sprintf("\ngen %d, rep 1: log10(y + 1) mean w168 %.4f, w192 %.4f\n", g,
                    mean(log10(a$d$y + 1)), mean(log10(b$d$y + 1))))
        stopifnot(identical(a$truth, b$truth))
        shared <- match(a$key, b$key); stopifnot(!anyNA(shared))
        cat(sprintf("  w192 at the 34 shared observations vs w168: max |diff| %.2e\n",
                    max(abs(b$d$y[shared] - a$d$y))))
        stopifnot(isTRUE(all.equal(b$d$y[shared], a$d$y, tolerance = 1e-12)))
        if (g == 384L) {
            cat(sprintf("  w168 vs simstudy.R B-g384 rep 1 at these observations: max |diff| %.2e\n",
                        max(abs(a$d$y - old_y[a$pos]))))
            stopifnot(isTRUE(all.equal(a$d$y, old_y[a$pos], tolerance = 1e-10)))
        }
        ## Split invariance of log_prob (alternating vs first-half / second-half).
        for (x in list(a, b)) {
            xb <- make_fit_inputs(x$cf$name, split = "halves")
            fa <- rstan::sampling(mod, data = x$d, chains = 1, iter = 1, refresh = 0,
                                  seed = 1, algorithm = "Fixed_param")
            fb <- rstan::sampling(mod, data = xb$d, chains = 1, iter = 1, refresh = 0,
                                  seed = 1, algorithm = "Fixed_param")
            set.seed(7); dmax <- 0
            for (j in 1:3) {
                up <- rnorm(rstan::get_num_upars(fa), 0, 0.5)
                dmax <- max(dmax, abs(rstan::log_prob(fa, up) - rstan::log_prob(fb, up)))
            }
            cat(sprintf("  %s split check: max |diff| in log_prob over 3 points %.2e\n",
                        x$cf$name, dmax))
            stopifnot(dmax < 1e-6)
        }
        ## Prior arm: identical data, only the cycle_length prior's centre moves.
        pa <- make_fit_inputs(sprintf("w168-g%d-est-rep1-ptruth", g))
        stopifnot(identical(pa$d$y, a$d$y))
        nd <- setdiff(names(a$d)[!vapply(names(a$d), \(n) isTRUE(all.equal(a$d[[n]], pa$d[[n]])), logical(1))], character(0))
        cat(sprintf("  prior arm: same y; fields differing from production: %s; mean_logit_cl %.4f -> %.4f (implied centre %.3f h)\n",
                    paste(nd, collapse = ", "), a$d$mean_logit_cl, pa$d$mean_logit_cl,
                    35 + 15 * plogis(pa$d$mean_logit_cl)))
        stopifnot(identical(nd, "mean_logit_cl"))
        ## Low-synchrony truth: only b_shape differs from the main truth.
        lo <- make_fit_inputs(sprintf("w168-g%d-fix-rep1-bslow", g))
        mn <- make_fit_inputs(sprintf("w168-g%d-fix-rep1", g))
        dt <- names(mn$truth)[!vapply(names(mn$truth), \(n) identical(mn$truth[[n]], lo$truth[[n]]), logical(1))]
        cat(sprintf("  bslow truth: fields differing from main: %s; b_shape %.3f (main %.3f); fixed b_shape in fit data %.3f\n",
                    paste(dt, collapse = ", "), lo$truth$b_shape, mn$truth$b_shape, lo$d$b_shape_data))
        stopifnot(setequal(dt, c("b_shape", "b_shape_384")), lo$d$b_shape_data == lo$truth$b_shape)
    }
    quit(save = "no")
}

if (MODE == "tasks") {
    ## rep-major, chains innermost
    nms <- STUDY$name
    which_tasks <- Sys.getenv("UNITSIM_TASKS", "main")
    stopifnot(which_tasks %in% c("main", "prior", "bslow"))
    if (which_tasks == "prior") nms <- paste0(nms, "-ptruth")
    if (which_tasks == "bslow") nms <- paste0(nms, "-bslow")
    cat(paste0(rep(nms, each = N_CHAINS), "-c", seq_len(N_CHAINS)), sep = "\n")
    quit(save = "no")
}

if (MODE == "chain") {
    tok <- Sys.getenv("SLURM_ARRAY_TASK_ID", "")
    m <- regmatches(tok, regexec("^(.*)-c([0-9]+)$", tok))[[1]]
    stopifnot(length(m) == 3L)
    nm <- m[2]; k <- as.integer(m[3]); stopifnot(k >= 1L, k <= N_CHAINS)
    out <- chain_path(nm, k)
    st <- SETTINGS[[parse_fit(nm)$method]]
    ITER <- st$iter; WARMUP <- st$warmup
    if (SMOKE) { ITER <- as.integer(iter_env); WARMUP <- ITER %/% 2L }
    cat(sprintf("task %s: fit %s, chain %d of %d, %d iter (%d warmup), cpus %s\n",
                tok, nm, k, N_CHAINS, ITER, WARMUP,
                sub(".*:\\s*", "", grep("Cpus_allowed_list",
                                        readLines("/proc/self/status"), value = TRUE))))
    if (file.exists(out)) { cat("exists, skipping:", out, "\n"); quit(save = "no") }
    x <- make_fit_inputs(nm)
    fk <- archer_fit(x$d, model = "pooled_both", chains = 1L, iter = ITER,
                     warmup = WARMUP, seed = SEED_FIT, chain_id = k,
                     threads_per_chain = 2L, control = CONTROL)
    write_rds(list(fit = fk, d = x$d, truth = x$truth, cf = x$cf, chain = k,
                   seed = SEED_FIT, iter = ITER, warmup = WARMUP,
                   control = CONTROL), out)
    print(rstan::get_elapsed_time(fk))
    cat("wrote", out, "\n")
    quit(save = "no")
}

if (MODE == "merge") {
    files <- Sys.glob(paste0(PFX, "chain-*-c*.rds"))
    names_all <- unique(sub(paste0("^", PFX, "chain-(.*)-c[0-9]+\\.rds$"), "\\1", files))
    n_done <- 0L
    for (nm in names_all) {
        paths <- vapply(seq_len(N_CHAINS), \(k) chain_path(nm, k), character(1))
        if (file.exists(sum_path(nm)) || !all(file.exists(paths))) next
        parts <- lapply(paths, read_rds)
        for (pt in parts) stopifnot(identical(pt$d, parts[[1]]$d),
                                    identical(pt$truth, parts[[1]]$truth),
                                    pt$seed == SEED_FIT)
        stopifnot(identical(vapply(parts, `[[`, integer(1), "chain"), seq_len(N_CHAINS)))
        f <- rstan::sflist2stanfit(lapply(parts, `[[`, "fit"))
        tru <- parts[[1]]$truth; cf <- parts[[1]]$cf
        n_draws <- parts[[1]]$iter - parts[[1]]$warmup
        el <- do.call(rbind, lapply(parts, \(p) rstan::get_elapsed_time(p$fit)))
        rm(parts)
        write_rds(f, fit_path(nm))

        s <- rstan::summary(f)$summary
        rh <- s[, "Rhat"]; rh <- rh[is.finite(rh)]
        cl <- as.matrix(f, pars = "cycle_length")[, 1]
        bs <- as.matrix(f, pars = "b_shape")[, 1]
        lp_chain <- apply(rstan::extract(f, "lp__", permuted = FALSE), 2, mean)
        lp_apart <- max(vapply(seq_along(lp_chain), \(k)
                                   abs(lp_chain[k] - median(lp_chain[-k])), numeric(1)))
        sp <- rstan::get_sampler_params(f, inc_warmup = FALSE)
        lf <- vapply(sp, \(x) sum(x[, "n_leapfrog__"]), numeric(1))
        q <- quantile(cl, c(0.1, 0.9))
        row <- tibble(name = nm, unit = tru$unit, window = cf$window,
                      prior = if (is.null(cf$prior)) "prod" else cf$prior,
                      truth_set = if (is.null(cf$truth_set)) "main" else cf$truth_set,
                      gen_nc = cf$gen_nc, method = cf$method,
                      rep = cf$rep, max_rhat = max(rh), worst = names(which.max(rh)),
                      pass = max(rh) < 1.05,
                      true_cl = tru$cycle_length, cl_mean = mean(cl), cl_sd = sd(cl),
                      bias = mean(cl) - tru$cycle_length,
                      in_80 = tru$cycle_length >= q[1] && tru$cycle_length <= q[2],
                      bs_logratio_gen = log(mean(bs) / tru$b_shape),
                      bs_logratio_384 = log(mean(bs) / tru$b_shape_384),
                      lp_spread = diff(range(lp_chain)), lp_apart = lp_apart,
                      div_pct = 100 * mean(unlist(lapply(sp, \(x) x[, "divergent__"]))),
                      treedepth_cap_pct = 100 * mean(unlist(lapply(sp, \(x)
                                              x[, "treedepth__"] >= CONTROL$max_treedepth))),
                      hours_max_chain = max(rowSums(el)) / 3600)
        row$per_chain <- list(tibble(chain = seq_len(N_CHAINS), lp_mean = lp_chain,
                                     cl_mean = apply(rstan::extract(f, "cycle_length", permuted = FALSE)[, , 1], 2, mean),
                                     hours = rowSums(el) / 3600,
                                     warmup_hours = el[, "warmup"] / 3600,
                                     leapfrog_per_iter = lf / n_draws,
                                     ms_per_grad = 1e3 * el[, "sample"] / lf))
        write_rds(row, sum_path(nm))
        n_done <- n_done + 1L
        cat(sprintf("merged %s: max R-hat %.3f %s, bias %+.3f h, %.2f h longest chain\n",
                    nm, row$max_rhat, if (row$pass) "PASS" else "FAIL",
                    row$bias, row$hours_max_chain))
    }
    cat("merged", n_done, "fits this run;",
        length(Sys.glob(paste0(PFX, "sum-*.rds"))), "fits merged in all\n")
    quit(save = "no")
}

stop("unknown UNITSIM_MODE: ", MODE)
