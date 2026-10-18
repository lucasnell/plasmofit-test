## Estimate when the simulation study's chain array (SLURM 30875) finishes.
## Replays the remaining queue in array order through 96 slots, drawing each
## chain's run time from the completed chains of the same class (fit n_c x
## b_shape method). Running chains get their class time minus time already
## run, resampled until positive. 1000 replays; reports quantiles of the end.
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   Rscript --vanilla _scripts/simstudy-eta.R > _data/simstudy-eta-<date>.txt

set.seed(1)
SLOTS <- 96L
cls <- function(id) {
    f <- (id - 1L) %/% 8L + 1L; j <- (f - 1L) %% 48L + 1L
    paste0(ifelse(j <= 24L, "f384", "f192"), "-", ifelse(j %% 2L == 1L, "wide", "fixed"))
}
secs <- function(x) vapply(strsplit(x, "[-:]"), function(p) {
    p <- as.numeric(p)
    if (length(p) == 4) sum(p * c(86400, 3600, 60, 1))
    else if (length(p) == 3) sum(p * c(3600, 60, 1)) else sum(p * c(60, 1))
}, numeric(1)) / 3600

done <- read.table(pipe("sacct -j 30875 -X -n -P --format=JobID,Elapsed,State"),
                   sep = "|", col.names = c("job", "el", "state"))
done <- done[done$state == "COMPLETED", ]
done$id <- as.integer(sub(".*_", "", done$job)); done$h <- secs(done$el)
pool <- split(done$h, cls(done$id))

q <- read.table(pipe("squeue -j 30875 -h -r -t PD,R -o '%i %T %M'"),
                col.names = c("job", "state", "el"))
q$id <- as.integer(sub(".*_", "", q$job)); q$cl <- cls(q$id)
run <- q[q$state == "RUNNING", ]; pend <- q[q$state == "PENDING", ]
pend <- pend[order(pend$id), ]
run$done_h <- secs(run$el)

cat("completed chains by class: mean hours (n)\n")
for (k in names(pool)) cat(sprintf("  %-11s %.2f (%d)\n", k, mean(pool[[k]]), length(pool[[k]])))
cat(sprintf("remaining: %d running, %d pending\n", nrow(run), nrow(pend)))
print(table(class = q$cl, state = q$state))

ends <- replicate(1000, {
    free <- vapply(seq_len(nrow(run)), function(r) {
        p <- pool[[run$cl[r]]]; left <- p[p > run$done_h[r]]
        if (length(left) == 0) 0.25 else sample(left, 1) - run$done_h[r]
    }, numeric(1))
    slots <- c(free, rep(0, SLOTS - length(free)))
    for (r in seq_len(nrow(pend))) {
        s <- which.min(slots)
        slots[s] <- slots[s] + sample(pool[[pend$cl[r]]], 1)
    }
    max(slots)
})
now <- Sys.time()
qs <- quantile(ends, c(0.1, 0.5, 0.9))
cat(sprintf("\nnow %s\nchain array ends, hours from now: 10%% %.1f, median %.1f, 90%% %.1f\n",
            format(now, "%Y-%m-%d %H:%M"), qs[1], qs[2], qs[3]))
cat("as clock times:", paste(format(now + qs * 3600, "%a %Y-%m-%d %H:%M"), collapse = " | "), "\n")
