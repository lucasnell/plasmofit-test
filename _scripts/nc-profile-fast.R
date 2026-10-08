## The dispersion profile, densely, using the convolution forward map.
## Thread 15; the fast replacement for the discrete ladder, and the IPM's own
## identifiability check. See claude/ipm-decision.md.
##
## WHY THIS EXISTS. The convolution in _scripts/convolution-forward-map.R is
## the SAME model as mat_exp_series, agreeing to 1e-12, and a whole unit fit
## takes 12 s at n_c = 384 where SLURM 29684 took 1777 s. That buys two things
## the discrete ladder could not afford:
##   A. the exact chain at nine integer rungs instead of three, and
##   B. the gamma-kernel IPM at a FIXED mesh with CONTINUOUS n_eff, which is
##      the forward map an IPM would actually use and the only way to ask
##      whether the dispersion is identified rather than quantised.
##
## THE READING RULE IS NOT NEW. It is the one pre-registered in
## _scripts/nc-dispersion-profile-read.R before any of this existed: turns
## over / saturates / still climbing / flat, judged on whether the
## 2-log-likelihood interval is bounded in the fine direction.
##
## BUILT-IN REGRESSION CHECK. The chain rungs at 96, 192 and 384 must
## reproduce SLURM 29684's saved log-likelihoods. The script compares them and
## stops if any differs by more than TOL, because a silent drift here would
## make every number below incomparable with the rest of thread 15.
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   sbatch _scripts/nc-profile-fast.sh
##   Rscript --vanilla _scripts/nc-profile-fast-read.R | tee _data/nc-profile-fast-$(date +%F).txt

suppressPackageStartupMessages({library(readr); library(dplyr); library(tibble)})
source("_scripts/convolution-forward-map.R")

NC_CHAIN <- c(64L, 96L, 128L, 192L, 256L, 384L, 512L, 768L, 1024L)
M_IPM    <- 192L                       # mesh held FIXED, which is the point
NE_IPM   <- round(exp(seq(log(48), log(4096), length.out = 16)), 2)
M_CHECK  <- 384L                       # mesh-convergence check
## MUST be values that appear in NE_IPM, or the reader pairs them against
## nothing and reports NA. They were not, in the 2026-10-08 run; fixed here.
## Chosen at the fine end, because that is where the verdict is decided.
NE_CHECK <- c(514.25, 2264.01, 4096.00)
MU       <- 0
TOL      <- 0.01                       # log-likelihood units

## Start count, 2026-10-08. Two starts were MEASURED to be too few at fine
## rungs: _scripts/profile-noise-check.R refit n_eff = 4096 with eight and
## DSM265|1800 alone gained +11.29 log-likelihood units, which was the whole
## of the pooled drop that had looked like a turnover. The first two entries
## are the original pair, so NCPF_STARTS=2 reproduces the earlier run exactly.
STARTS_ALL <- list(
    c( 0.0, log(15),   0.25, 1.0, log(5.5)),
    c( 0.5, log(60),   0.60, 0.5, log(8)),
    c(-0.5, log(8),    0.05, 2.0, log(3)),
    c( 0.9, log(300),  0.85, 0.0, log(15)),
    c( 0.2, log(1500), 0.40, 1.5, log(6)),
    c(-0.2, log(40),   0.70, 0.8, log(10)),
    c( 0.7, log(5),    0.15, 1.2, log(20)),
    c( 0.0, log(900),  0.95, 0.3, log(4)))
N_STARTS <- as.integer(Sys.getenv("NCPF_STARTS", "2"))
if (is.na(N_STARTS) || N_STARTS < 1 || N_STARTS > length(STARTS_ALL))
    stop("NCPF_STARTS must be 1-", length(STARTS_ALL))
STARTS <- STARTS_ALL[seq_len(N_STARTS)]

## NCPF_TAG appends to the output name so a re-run lands BESIDE the earlier
## one instead of destroying the thing it is being compared against.
TAG <- Sys.getenv("NCPF_TAG", "")
OUT <- sprintf("_data/nc-profile-fast-unit%%02d%s.rds", TAG)
cat("starts:", N_STARTS, "| output:", sprintf(OUT, 0), "\n")

## Refuse to destroy an earlier run. The untagged name is the 2026-10-08
## two-start run, which is the baseline every later run is compared against,
## and the default TAG is empty so that NCPF_STARTS=2 can reproduce it. Those
## two facts together make an accidental overwrite easy, so it is blocked.
check_out <- function(path) {
    if (file.exists(path) && !nzchar(Sys.getenv("NCPF_FORCE")))
        stop(path, " already exists. Set NCPF_TAG to write beside it, or ",
             "NCPF_FORCE=1 to overwrite deliberately.")
}

d <- read_csv("_data/wockner-cleaned.csv", col_types = "cccdcdd") |>
    mutate(unit = paste(trial, inoc_size, sep = "|"))
units <- sort(unique(d$unit))
task <- suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", NA)))
if (is.na(task)) stop("set SLURM_ARRAY_TASK_ID to 1-", length(units))
if (task < 1 || task > length(units))
    stop("SLURM_ARRAY_TASK_ID must be 1-", length(units), ", got ", task)

u  <- units[task]; du <- filter(d, unit == u)
yo <- du$para
## mat_exp_series needed strictly increasing times; the convolution does not,
## but the trajectory is still evaluated once per distinct time and indexed
## back out, because that is where the cost is.
ut  <- sort(unique(du$time)); idx <- match(du$time, ut)
stopifnot(!is.unsorted(ut, strictly = TRUE), !anyNA(idx))
cat("unit", task, "of", length(units), ":", u, "|", nrow(du), "observations,",
    length(ut), "distinct times\n"); flush.console()
check_out(sprintf(OUT, task))

## Same parameterisation and same profiled normal likelihood on log10(y + 1)
## as _scripts/decay-law-test.R, so the numbers are comparable to 29684.
nll <- function(par, fwd) {
    cl <- 35 + 15 / (1 + exp(-par[1])); bs <- exp(par[2]); bo <- par[3] %% 1
    lt0 <- par[4]; R <- 1 + exp(par[5])
    if (!is.finite(bs) || bs > 5000 || bs <= 2 || !is.finite(R) || R > 200)
        return(1e10)
    yu <- tryCatch(fwd(cl, bs, bo, lt0, R), error = function(e) NULL)
    if (is.null(yu) || any(!is.finite(yu)) || any(yu < 0)) return(1e10)
    r <- log10(yu[idx] + 1) - log10(yo + 1)
    n <- length(r); s2 <- sum(r^2) / n
    if (!is.finite(s2) || s2 <= 0) return(1e10)
    0.5 * n * (log(2 * pi * s2) + 1)
}
fit1 <- function(fwd) {
    best <- NULL
    for (s0 in STARTS) {
        o <- tryCatch(optim(s0, nll, fwd = fwd, method = "Nelder-Mead",
                            control = list(maxit = 2000, reltol = 1e-9)),
                      error = function(e) NULL)
        if (!is.null(o) && (is.null(best) || o$value < best$value)) best <- o
    }
    best
}
row1 <- function(model, knot, M, o, secs) tibble(
    unit = u, model = model, knot = knot, M = M, ll = -o$value,
    cl_hat = 35 + 15 / (1 + exp(-o$par[1])), bs_hat = exp(o$par[2]),
    bo_hat = o$par[3] %% 1, lt0_hat = o$par[4], R_hat = 1 + exp(o$par[5]),
    n = length(yo), secs = secs)

res <- list()
for (nc in NC_CHAIN) {
    t0 <- proc.time()[["elapsed"]]
    o <- fit1(function(cl, bs, bo, lt0, R) cfm_chain(ut, cl, nc, R, MU, bs, bo, lt0))
    res[[length(res) + 1]] <- row1("chain", nc, nc, o, proc.time()[["elapsed"]] - t0)
    cat(sprintf("  chain  n_c=%5d  ll %+.4f  (%.0f s)\n", nc, -o$value,
                proc.time()[["elapsed"]] - t0)); flush.console()
}
for (ne in NE_IPM) {
    t0 <- proc.time()[["elapsed"]]
    o <- fit1(function(cl, bs, bo, lt0, R) cfm_gamma(ut, cl, M_IPM, ne, R, MU, bs, bo, lt0))
    res[[length(res) + 1]] <- row1("ipm_gamma", ne, M_IPM, o, proc.time()[["elapsed"]] - t0)
    cat(sprintf("  gamma  n_eff=%7.1f  ll %+.4f  (%.0f s)\n", ne, -o$value,
                proc.time()[["elapsed"]] - t0)); flush.console()
}
## Is the fixed mesh fine enough at the extremes? Repeat three rungs on a
## mesh twice as fine; a large shift would mean the IPM rows are mesh artefacts.
for (ne in NE_CHECK) {
    t0 <- proc.time()[["elapsed"]]
    o <- fit1(function(cl, bs, bo, lt0, R) cfm_gamma(ut, cl, M_CHECK, ne, R, MU, bs, bo, lt0))
    res[[length(res) + 1]] <- row1("ipm_gamma_meshcheck", ne, M_CHECK, o,
                                   proc.time()[["elapsed"]] - t0)
    cat(sprintf("  mesh   n_eff=%7.1f M=%d  ll %+.4f  (%.0f s)\n", ne, M_CHECK,
                -o$value, proc.time()[["elapsed"]] - t0)); flush.console()
}
r <- bind_rows(res)

## Regression check against the validated mat_exp_series path.
old <- sprintf("_data/decay-law-unit%02d.rds", task)
if (file.exists(old)) {
    o <- readRDS(old) |> filter(model == "A_sqrt") |> select(knot = n_c, ll_old = ll)
    cmp <- r |> filter(model == "chain") |> select(knot, ll_new = ll) |>
        inner_join(o, by = "knot") |> mutate(d = ll_new - ll_old)
    print(as.data.frame(cmp), digits = 7)
    if (nrow(cmp) > 0 && max(abs(cmp$d)) > TOL)
        stop("chain rungs disagree with SLURM 29684 by more than ", TOL,
             " log-likelihood units -- the convolution harness has drifted ",
             "from mat_exp_series and nothing here is comparable.")
    cat("regression check vs 29684: max |difference| =",
        sprintf("%.2e", max(abs(cmp$d))), "log-likelihood units, PASS\n")
} else {
    cat("NOTE: no", old, "to check against; regression check skipped.\n")
}

saveRDS(r, sprintf(OUT, task))
cat("\nwrote ", sprintf(OUT, task), "\n", sep = "")
