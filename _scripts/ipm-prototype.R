## Prototype and validate a transport-with-dispersion (IPM) forward map
## against the Erlang chain. Thread 15; the de-risking step named in
## claude/ipm-decision.md.
##
## THE DERIVATION, which is what makes this cheap. Work in ABSOLUTE
## developmental age rather than age modulo the cycle. Then:
##   - Transport is a pure-birth process: a parasite's absolute stage after
##     time t is its starting stage plus Poisson(lambda * t), lambda = n_c/cl.
##     No wrap, no boundary condition.
##   - Growth is a weight, not a flux: a parasite at absolute stage j has
##     divided d = floor((j-1)/n_c) times, so it carries R^d.
##   - Sequestration is a weight too, because circulating status is a function
##     of position within the cycle and RESETS at division (the ODE sends both
##     N[n_c] and S[n_c] into N[1]).
## So the whole trajectory is a convolution, evaluated once per observation
## time. There is no time stepping, hence NO NUMERICAL DIFFUSION -- the hazard
## flagged in claude/ipm-decision.md does not arise, because the kernel is
## applied analytically rather than accumulated over steps.
##
## TWO SUBTLETIES, both found by validating rather than by reading the source.
##  1. The chain applies the hazard with a ONE-STAGE LAG: the transition from
##     stage k-1 to k uses q[k-1], so P(circulating at stage k) is G[k-1], not
##     y[k]. This is the O(1/n_c) delay the Stan source documents.
##  2. Parasites in their FIRST cycle have not been reset, so their weight
##     depends on where they started: P(circ at k | start i) = y[i]*G[k-1]/G[i-1].
##     That is separable, so it costs one extra convolution rather than
##     breaking the method.
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   Rscript --vanilla _scripts/ipm-prototype.R | tee _data/ipm-prototype-$(date +%F).txt

suppressPackageStartupMessages(library(plasmofit))

CL <- 45; BS <- 400; BO <- 0.3; LT0 <- 5.5; RR <- 8; MU <- 0
TS <- c(72, 96, 120, 144, 168, 192, 216); DTF <- 12

## The forward maps live in one place, so this script and the fitting script
## cannot drift apart. The derivation and its two load-bearing subtleties are
## documented there.
source("_scripts/convolution-forward-map.R")

mesh_weights <- function(M, cl = CL) cfm_weights(M, cl, BS, BO, LT0)
chain_conv <- function(ts, cl = CL, n_c, R = RR, mu = MU)
    cfm_chain(ts, cl, n_c, R, mu, BS, BO, LT0)
ipm_gauss <- function(ts, cl = CL, M, sigma_c, R = RR, mu = MU)
    cfm_gauss(ts, cl, M, sigma_c, R, mu, BS, BO, LT0)
ipm_gamma <- function(ts, cl = CL, M, n_eff, R = RR, mu = MU)
    cfm_gamma(ts, cl, M, n_eff, R, mu, BS, BO, LT0)

## The likelihood is normal on log10(y + 1), so discrepancies are measured in
## log10 units, not as relative error. A relative error is dominated by the
## troughs, where the observable is four orders of magnitude below the peaks.
d10 <- function(a, b) max(abs(log10(a) - log10(b)))

cat("=== 1. the Poisson convolution reproduces mat_exp_series ===\n")
cat("Cells: max over the 7 observation times of |convolution - mat_exp_series|\n",
    "/ mat_exp_series, at cycle_length 45, b_shape 400, R 8, mu 0. This is a\n",
    "reimplementation from the ODE's structure, not a refactor, so agreement\n",
    "at 1e-12 is an independent check of both.\n\n", sep = "")
for (n_c in c(96L, 192L, 384L)) {
    y0  <- plasmofit:::generate_starts(CL, n_c, BS, BO, LT0)
    ref <- plasmofit:::mat_exp_series(y0, CL, n_c, RR, MU, TS, DTF)
    got <- chain_conv(TS, n_c = n_c)
    cat(sprintf("  n_c = %3d : max rel err = %.3e\n", n_c,
                max(abs(got - ref) / ref)))
}

cat("\n=== 2. a continuous kernel is NOT a reparameterisation of the chain ===\n")
cat("Cells: max over the 7 times of |log10(IPM) - log10(chain)|, in log10\n",
    "units, at mesh M = n_c with the kernel's variance matched to the chain's.\n",
    "The only difference is the kernel: the chain's stage at time t is exactly\n",
    "POISSON, a lattice distribution, while both IPM kernels are continuous.\n",
    "For scale, the residual sd of these ML fits is about 0.48 log10 units.\n\n",
    sep = "")
for (n_c in c(96L, 192L, 384L)) {
    ch <- chain_conv(TS, n_c = n_c)
    cat(sprintf("  n_c = %3d : gaussian %.4f | gamma %.4f\n", n_c,
                d10(ipm_gauss(TS, M = n_c, sigma_c = CL / sqrt(n_c)), ch),
                d10(ipm_gamma(TS, M = n_c, n_eff = n_c), ch)))
}
cat("\nThe gamma keeps the chain's right skew and the gaussian does not, yet\n",
    "they differ by under 0.003 log10 units, so SKEW IS NOT THE GAP. The gap\n",
    "grows with n_c, which is the opposite of the kernel-shape argument, and\n",
    "the next table says why.\n", sep = "")

cat("\n=== 3. the gap is one point, and it is the deepest trough ===\n")
cat("Cells: the chain's predicted circulating parasitaemia and the gamma IPM's,\n",
    "per observation time, at n_c = 384; dlog10 is log10(IPM) - log10(chain).\n",
    "The observable spans four orders of magnitude across the cycle because\n",
    "sequestration hides ~60% of it, and at a trough the value is set by the\n",
    "TAIL of the age distribution -- exactly where a lattice Poisson and a\n",
    "continuous gamma differ.\n\n", sep = "")
ch <- chain_conv(TS, n_c = 384L); ip <- ipm_gamma(TS, M = 384L, n_eff = 384)
print(data.frame(t = TS, cycles = round(TS / CL, 2), chain = ch, ipm = ip,
                 dlog10 = log10(ip) - log10(ch)), digits = 5)
cat(sprintf("\ntrough / peak ratio in the chain: %.0e. All of the discrepancy is at\n",
            min(ch) / max(ch)))
cat("t = 72 h; the other six times agree to 0.006 log10 units or better.\n")
cat("Troughs are also where low-end censoring is an open question (TODO), so\n",
    "this is the one place the kernel choice and a known data issue coincide.\n",
    sep = "")

cat("\n=== 4. the decoupling, which is the whole point ===\n")
cat("Cells: max |log10 difference| from the chain at n_c = 384, the dispersion\n",
    "the data prefer. Row 1 is the chain forced onto a coarse mesh, which also\n",
    "coarsens its dispersion because it has only one knob. The IPM rows use\n",
    "the SAME coarse mesh with the fine dispersion dialled in.\n\n", sep = "")
tgt <- chain_conv(TS, n_c = 384L)
cat(sprintf("  chain, n_c = 96 (mesh and rate locked together) : %.4f\n",
            d10(chain_conv(TS, n_c = 96L), tgt)))
for (M in c(96L, 192L))
    cat(sprintf("  IPM gamma, M = %3d, n_eff = 384                 : %.4f\n",
                M, d10(ipm_gamma(TS, M = M, n_eff = 384), tgt)))
cat("\nA factor of about 25 in log10 units. Mesh and rate are separated.\n")

cat("\n=== 5. dispersion is continuous for the IPM and quantised for the chain ===\n")
cat("Cells: sigma_c is the sd in hours of a parasite's age advance after one\n",
    "cycle; the last column is the stage distribution's sd at the final\n",
    "observation. The chain reaches only these rows, at cubic cost to move\n",
    "between them; the IPM reaches any value at a fixed mesh and fixed cost.\n\n",
    sep = "")
nn <- c(96, 192, 384, 768)
print(data.frame(n_c = nn, sigma_c_h = CL / sqrt(nn),
                 sd_at_4.8_cycles_h = CL * sqrt(4.8 / nn)), digits = 4)

cat("\n=== 6. cost ===\n")
cat("Cells: seconds for one 7-point trajectory, median of 5 runs.\n",
    "mat_exp_series exponentiates a 2*n_c square matrix, which is cubic; the\n",
    "convolutions are one FFT pair per time point.\n\n", sep = "")
tm <- function(f, n = 5) median(replicate(n, system.time(f())[["elapsed"]]))
for (n_c in c(96L, 192L, 384L)) {
    y0 <- plasmofit:::generate_starts(CL, n_c, BS, BO, LT0)

    a <- tm(function() plasmofit:::mat_exp_series(y0, CL, n_c, RR, MU, TS, DTF))
    b <- tm(function() chain_conv(TS, n_c = n_c))
    cc <- tm(function() ipm_gamma(TS, M = n_c, n_eff = n_c))
    cat(sprintf("  n_c = %3d : mat_exp %7.4f | poisson-conv %7.4f | ipm-gamma %7.4f | speedup %5.0fx\n",
                n_c, a, b, cc, a / cc))
}
cat("\nThe speedup is available WITHOUT any rewrite: the Poisson convolution is\n",
    "the same model as mat_exp_series, agreeing to 1e-12, and it is the column\n",
    "to port if cost is the motive rather than decoupling.\n", sep = "")
