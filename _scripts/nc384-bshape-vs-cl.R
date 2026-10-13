## Does a trial's cycle_length shift between the pinned-15 and free-b_shape
## fits at n_c = 384 track how far its free b_shape sits from 15?
## Inputs: SLURM 30830 entry 43 (np_bs15_nc384) and 30842 entry 54
## (np_wide_total0_nc384_t), both converged.
##
## Cells of the printed table: b_free is the posterior mean of b_shape for
## that grp_init unit in the free fit; d_cl is free minus pinned-15 posterior
## mean cycle_length (h) for the unit's trial (grp_cl). SJ733IBSMCS has two
## units and so appears twice with one d_cl.
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   Rscript --vanilla _scripts/nc384-bshape-vs-cl.R \
##       > _data/nc384-bshape-vs-cl-<date>.txt

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
suppressPackageStartupMessages(library(rstan))

d <- readRDS("_data/wock-data-np_wide_total0_nc384_t.rds"); lv <- attr(d, "levels")
stopifnot(identical(lv, attr(readRDS("_data/wock-data-np_bs15_nc384.rds"), "levels")))
b    <- colMeans(extract(readRDS("_data/wock-fit-np_wide_total0_nc384_t.rds"), "b_shape")[[1]])
cl15 <- colMeans(extract(readRDS("_data/wock-fit-np_bs15_nc384.rds"), "cycle_length")[[1]])
clf  <- colMeans(extract(readRDS("_data/wock-fit-np_wide_total0_nc384_t.rds"), "cycle_length")[[1]])
tr <- sub("\\.[0-9]+$", "", lv$inoc)
stopifnot(all(tr %in% lv$trial))
x <- data.frame(unit = lv$inoc, b_free = round(b, 2),
                d_cl = round((clf - cl15)[match(tr, lv$trial)], 2))
print(x[order(-abs(x$d_cl)), ], row.names = FALSE)
cat("\nSpearman, |log(b_free/15)| vs |d_cl|, 14 units:",
    round(cor(abs(log(x$b_free / 15)), abs(x$d_cl), method = "spearman"), 3), "\n")
cat("Spearman, log(b_free) vs d_cl, 14 units:",
    round(cor(log(x$b_free), x$d_cl, method = "spearman"), 3), "\n")
