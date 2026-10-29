## Does the total-spread formula used to map b_shape across n_c hold for the
## model itself? (2026-10-29, before the fixed-pin grid runs.)
##
## The formula (unitsim.R header; simstudy2.R): stage-position variance in
## cycles^2 at time t is var_init(b) + k / n_c, with var_init(b) = 1 / (4 (2b
## + 1)) (a Beta(b, b) start) and k = t / cycle_length (Poisson progress
## through n_c compartments). It treats each parasite separately; the
## model's population also multiplies R-fold at each cycle end, which favours
## lineages that have advanced further, and the positions wrap at 1.
##
## Here the model's own state is propagated exactly: y0 = generate_starts(),
## state(t) = expm(A t) y0 with A = build_A() (the package's linear chain,
## circulating and sequestered compartments), at the main truth's parameters
## (cycle_length 42.844 h, R 13.87, b_offset 0.226, mu 0). The stage
## distribution at t is circulating + sequestered abundance per stage,
## normalised. Its spread is reported as a wrapped-normal sd in cycles,
## sqrt(-2 log rho) / (2 pi), with rho the mean resultant length of the
## positions (i - 0.5) / n_c; this equals the sd for a wrapped normal and is
## close to it for any distribution as narrow as these.
##
## Pairs: the reference (b_shape 20.10, n_c 384); the main study's gen 1024
## mapping (9.54, 1024); the fixed-pin grid's spread-matched pins (46.41,
## 288) and (14.00, 512); and the truth pin at the other fit n_c.
##
##   Rscript --vanilla _scripts/unitsim-spread-check.R > _data/unitsim-check-spread-<date>.txt

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
suppressPackageStartupMessages({library(plasmofit); library(Matrix)})
CL <- 42.844; R0 <- 13.87; B_OFF <- 0.226; MU <- 0
TIMES <- c(0, 72, 120, 168, 192)
STEP <- 12
vi <- function(b) 1 / (4 * (2 * b + 1))
formula_sd <- function(b, n_c, t) sqrt(vi(b) + (t / CL) / n_c)
wrapped_sd <- function(p, n_c) {
    x <- (seq_len(n_c) - 0.5) / n_c
    z <- sum(p * exp(2i * pi * x)) / sum(p)
    sqrt(-2 * log(Mod(z))) / (2 * pi)
}

model_sd <- function(b, n_c, R = R0) {
    y0 <- plasmofit:::generate_starts(CL, n_c, b, B_OFF, 1)
    q <- plasmofit:::make_q_vals(CL, n_c)
    A <- plasmofit:::build_A(CL, MU, MU, q, R, n_c)
    P <- as.matrix(expm(Matrix(A * STEP)))
    st <- y0; out <- numeric(length(TIMES)); tt <- 0
    for (j in seq_along(TIMES)) {
        while (tt < TIMES[j]) { st <- P %*% st; tt <- tt + STEP }
        p <- st[1:n_c] + st[(n_c + 1):(2 * n_c)]
        stopifnot(all(p >= -1e-8 * max(p)))
        out[j] <- wrapped_sd(pmax(p, 0), n_c)
    }
    out
}

PAIRS <- data.frame(label = c("reference", "gen-1024 mapping", "match pin, fit 288",
                              "match pin, fit 512", "truth pin, fit 288", "truth pin, fit 512"),
                    b = c(20.098, 9.540, 46.41, 14.00, 20.098, 20.098),
                    n_c = c(384L, 1024L, 288L, 512L, 288L, 512L))
res <- do.call(rbind, lapply(seq_len(nrow(PAIRS)), \(i) {
    m <- model_sd(PAIRS$b[i], PAIRS$n_c[i])
    data.frame(label = PAIRS$label[i], b_shape = PAIRS$b[i], n_c = PAIRS$n_c[i], t_h = TIMES,
               formula_sd = formula_sd(PAIRS$b[i], PAIRS$n_c[i], TIMES),
               model_sd_R1 = model_sd(PAIRS$b[i], PAIRS$n_c[i], R = 1), model_sd = m)
}))
res$model_minus_formula <- res$model_sd - res$formula_sd
ref <- res[res$label == "reference", c("t_h", "model_sd")]
res$model_minus_reference <- res$model_sd - ref$model_sd[match(res$t_h, ref$t_h)]

cat("Stage-position spread (sd, cycles) by time. formula_sd = sqrt(1/(4(2b+1)) + (t/cl)/n_c);\n",
    "model_sd = wrapped-normal sd of the model's own stage distribution (circulating +\n",
    "sequestered) at the truth's R 13.87; model_sd_R1 = the same with R = 1 (no\n",
    "multiplication, so every lineage counts once: what the formula describes);\n",
    "model_minus_formula = model_sd - formula_sd; model_minus_reference = model_sd\n",
    "minus the reference's (b_shape 20.10,\n",
    "n_c 384) at the same time. The b_shape matching targets equal spread at 168 h.\n\n", sep = "")
cols <- c("formula_sd", "model_sd_R1", "model_sd", "model_minus_formula", "model_minus_reference")
res[, cols] <- round(res[, cols], 4)
options(width = 150)
print(res, row.names = FALSE)
