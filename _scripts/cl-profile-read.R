## Read the cycle_length profile likelihood (SLURM 30869, _scripts/cl-profile.R).
## Reading rule, fixed in cl-profile.R before any profile existed: TWO PEAKS if
## two local maxima are separated by a dip of at least 2 log-likelihood units
## below the lower of them; otherwise ONE.
##
## Cells: per unit and b_shape mode. ll is the profile log-likelihood (normal
## on log10(y + 1), one residual variance per unit, prior-free), d_ll = ll
## minus the unit-mode maximum (0 is best, more negative is worse). Peaks are
## local maxima of the grid profile; "dip" is the lowest ll between two
## adjacent peaks, reported as how far it falls below the lower peak.
## sweep_gain = how much the warm-start sweeps raised the best grid point over
## the 8 cold starts, the optimiser's own error on this surface.
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   Rscript --vanilla _scripts/cl-profile-read.R > _data/cl-profile-read-<date>.txt

fs <- Sys.glob("_data/cl-profile-*-*.rds")
cat("profiles present:", length(fs), "of 6\n")
for (f in fs) {
    r <- readRDS(f); r <- r[order(r$cl), ]
    r$d_ll <- r$ll - max(r$ll)
    n <- nrow(r); ll <- r$ll
    pk <- which(sapply(seq_len(n), function(i)
        (i == 1 || ll[i] > ll[i - 1]) && (i == n || ll[i] > ll[i + 1])))
    cat(sprintf("\n=== %s, b_shape %s (%d observations) ===\n", r$unit[1], r$mode[1], r$n[1]))
    cat(sprintf("max ll %.3f at cl %.1f h; largest sweep gain over cold starts %.4f\n",
                max(ll), r$cl[which.max(ll)], max(r$ll - r$cold_ll)))
    cat("peaks (cl h, ll, d_ll):", paste(sprintf("%.1f (%.2f, %.2f)", r$cl[pk], ll[pk], r$d_ll[pk]),
                                         collapse = "; "), "\n")
    two <- FALSE
    if (length(pk) >= 2) for (k in seq_len(length(pk) - 1)) {
        a <- pk[k]; b <- pk[k + 1]
        dip <- min(ll[a:b]); depth <- min(ll[a], ll[b]) - dip
        cat(sprintf("  between %.1f and %.1f h: dip to %.2f, %.2f below the lower peak\n",
                    r$cl[a], r$cl[b], dip, depth))
        if (depth >= 2) two <- TRUE
    }
    cat("VERDICT:", if (two) "TWO PEAKS" else "ONE", "\n")
    print(data.frame(cl = r$cl, ll = round(r$ll, 3), d_ll = round(r$d_ll, 2),
                     bs = round(r$bs, 1), bo = round(r$bo, 3), R = round(r$R, 2)),
          row.names = FALSE)
}
