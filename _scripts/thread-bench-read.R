## Read the thread-layout benchmark, SLURM 30846-30853 (2026-10-13), from
## _data/thread-bench-<jobid>.out as listed in
## _data/thread-bench-jobs-2026-10-13.txt. Layouts ran one at a time
## (--dependency=singleton); entry 42 (SLURM 30830_42, OS cores 0-1 and
## their siblings 128-129) ran throughout and is the only other load of ours.
##
## On cbsugreischar (2 x EPYC 7702) logical CPU c and c+128 are one physical
## core, and physical cores 4k..4k+3 share one 16 MB L3 cache ("L3 group" k).
## Each chain is labelled with the L3 groups its cores touch and how many
## OTHER busy chains (or entry 42) touch the same groups.
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   Rscript --vanilla _scripts/thread-bench-read.R \
##       > _data/thread-bench-read-2026-10-13.txt

jobs <- read.table("_data/thread-bench-jobs-2026-10-13.txt", col.names = c("tag", "job"))
expand <- function(s) unlist(lapply(strsplit(s, ",")[[1]], function(p) {
    r <- as.integer(strsplit(p, "-")[[1]]); if (length(r) == 1) r else r[1]:r[2] }))
phys <- function(cpus) unique(ifelse(cpus >= 128, cpus - 128, cpus))
ENTRY42 <- 0:1   # physical cores of the running production fit

rows <- list()
for (k in seq_len(nrow(jobs))) {
    L <- readLines(sprintf("_data/thread-bench-%s.out", jobs$job[k]))
    lay <- L[grepl("^# layout", L)]
    get <- function(key, s) sub(sprintf(".*%s=([^ ]+).*", key), "\\1", s)
    cpuset <- expand(get("cpuset", lay)); pin <- get("pin", lay)
    res <- L[grepl("^tag=", L)]
    ch <- L[grepl("^# chain", L)]
    n_ch <- length(res)
    ## cores per chain: pinned -> from "# chain i cpus"; unpinned -> whole set
    chain_cores <- if (length(ch)) lapply(sub(".*cpus ", "", ch), function(s) phys(expand(s)))
                   else rep(list(phys(cpuset)), n_ch)
    draw <- as.integer(get("draw", res))
    chain_of <- draw / 60       # thread-bench.sh: BENCH_DRAW = chain * 60
    groups <- lapply(chain_cores, function(cc) sort(unique(cc %/% 4)))
    for (j in seq_along(res)) {
        i <- chain_of[j]; g <- groups[[i]]
        others <- sum(vapply(setdiff(seq_along(groups), i), function(o)
            any(groups[[o]] %in% g), logical(1)))
        if (length(ch) == 0) others <- n_ch - 1   # unpinned: all chains share all cores
        rows[[length(rows) + 1]] <- data.frame(
            layout = jobs$tag[k], task_cpus = length(cpuset), pin = pin,
            chain = i, phys_cores = length(chain_cores[[i]]),
            threads = as.integer(get("active", res[j])),
            l3_groups = paste(g, collapse = "+"),
            shares_l3_with = others + any(ENTRY42 %/% 4 %in% g),
            ms_per_grad = as.numeric(get("ms_per_grad", res[j])),
            cpu_per_wall = as.numeric(get("cpu_per_wall", res[j])))
    }
}
R <- do.call(rbind, rows)
cat("Cells: one row per chain process. ms_per_grad is wall time per\n",
    "grad_log_prob call over 300 calls at a saved posterior draw of entry 43\n",
    "(np_bs15_nc384, n_c 384), milliseconds, lower is faster. threads = threads\n",
    "that did work; cpu_per_wall = CPU-seconds per wall-second for the process.\n",
    "shares_l3_with = number of other busy chains, plus entry 42, whose cores\n",
    "touch one of this chain's L3 groups. Unpinned layouts share everything.\n\n", sep = "")
print(R, row.names = FALSE)

cat("\n=== by layout ===\n")
cat("Cells: per layout, mean and range of ms_per_grad over its chains; grads_per_s_per_core\n",
    "= total gradients per second across the layout's chains divided by the physical\n",
    "cores the task held (throughput per core; higher is better).\n\n", sep = "")
S <- do.call(rbind, lapply(split(R, R$layout), function(x) data.frame(
    layout = x$layout[1], task_cpus = x$task_cpus[1], pin = x$pin[1],
    chains = nrow(x), threads_per_chain = x$threads[1],
    ms_mean = round(mean(x$ms_per_grad), 1),
    ms_range = sprintf("%.1f-%.1f", min(x$ms_per_grad), max(x$ms_per_grad)),
    grads_per_s_per_core = round(sum(1000 / x$ms_per_grad) / (x$task_cpus[1] / 2), 1))))
print(S[order(S$task_cpus, S$layout), ], row.names = FALSE)

cat("\n=== isolating L3 sharing: same layout, chains alone vs not alone in their L3 group ===\n")
for (lay in c("B", "F")) {
    x <- R[R$layout == lay, ]
    cat(sprintf("%s: alone %s ms | shared %s ms\n", lay,
                paste(x$ms_per_grad[x$shares_l3_with == 0], collapse = ", "),
                paste(x$ms_per_grad[x$shares_l3_with > 0], collapse = ", ")))
}

## The probe (SLURM 30843-30845, one chain, 4 CPUs, 30 gradients): neither
## STAN_NUM_THREADS=1 nor RcppParallel::setThreadOptions(1) changed how many
## threads did work. Raw lines, so that claim is re-derivable.
cat("\n=== probe: thread-count settings (30843-30845) ===\n")
for (j in readLines("_data/thread-bench-probe-jobs-2026-10-13.txt"))
    cat(grep("^tag=", readLines(sprintf("_data/thread-bench-%s.out", j)), value = TRUE), sep = "\n")
