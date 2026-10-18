## Paired PSIS-LOO for the n_c ladder read by nc-ladder-read.R.
##
## np_wide_both (n_c = 96), np_nc192, np_nc384 differ in n_c and nothing else.
## The elpd figures quoted for this ladder (findings.md, "`n_c` is a biological
## assumption") had no saved output behind them; this writes one. Reads the
## saved LOO objects only (their observation-level element), no refitting.
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   Rscript --vanilla _scripts/nc-ladder-loo.R > _data/nc-ladder-loo-<date>.txt

suppressPackageStartupMessages(library(loo))

ARMS <- c(np_wide_both = 96L, np_nc192 = 192L, np_nc384 = 384L)

loos <- lapply(names(ARMS), function(nm) {
    readRDS(sprintf("_data/wock-fit-LOO-%s.rds", nm))$observation
})
names(loos) <- sprintf("n_c %d (%s)", ARMS, names(ARMS))

cat("=== elpd by rung ===\n")
cat("Cells: elpd_loo and its se, summed over observations; n_obs is the\n",
    "number of pointwise terms; pareto_k_gt_0.7 counts observations with\n",
    "Pareto k above 0.7. Higher elpd is better.\n\n", sep = "")
tab <- do.call(rbind, lapply(names(loos), function(nm) {
    l <- loos[[nm]]
    e <- l$estimates
    data.frame(rung = nm, n_obs = nrow(l$pointwise),
               elpd = round(e["elpd_loo", "Estimate"], 1),
               se = round(e["elpd_loo", "SE"], 1),
               pareto_k_gt_0.7 = sum(l$diagnostics$pareto_k > 0.7))
}))
print(tab, row.names = FALSE)

cat("\n=== paired differences ===\n")
cat("Cells: elpd of the higher n_c minus the lower, summed over the shared\n",
    "observations; se is the standard error of the paired difference.\n",
    "Positive favours the higher n_c.\n\n", sep = "")
pair <- function(lo, hi) {
    d <- loos[[hi]]$pointwise[, "elpd_loo"] - loos[[lo]]$pointwise[, "elpd_loo"]
    stopifnot(length(d) == nrow(loos[[lo]]$pointwise))
    data.frame(comparison = sprintf("%s over %s", sub(" \\(.*", "", hi),
                                    sub(" \\(.*", "", lo)),
               elpd_diff = round(sum(d), 1),
               se_diff = round(sqrt(length(d)) * sd(d), 1))
}
nm <- names(loos)
print(rbind(pair(nm[1], nm[2]), pair(nm[1], nm[3]), pair(nm[2], nm[3])),
      row.names = FALSE)
