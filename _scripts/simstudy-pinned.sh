#!/bin/bash
# Pinned-layout workers for the simulation study (_scripts/simstudy.R): one
# job holds whole physical cores, splits them into workers of 4 physical
# cores (whole L3 groups first), and each worker runs chain tasks from a list
# one after another, each chain pinned with `taskset` to its 4 cores (4 Stan
# threads; the thread pool follows the CPUs a process may use).
#
# WHY. In the study's original layout (simstudy.sh: one single-threaded
# chain per core, 96 at once, 4 chains per L3 group) fit n_c 384 chains ran
# at ~260 ms per gradient, ~6.8x the 38 ms of one chain per core with few
# running (claude/findings.md, "Where the simulation study's time goes").
# The pinned layout ran entry 43 at 28 ms. Added 2026-10-18 to re-run the 192
# pending fit n_c 384 chain tasks cancelled from SLURM 30875.
#
# HOW. Workers claim task IDs from SIMSTUDY_PIN_LIST (one per line) with an
# atomic `mkdir _data/simstudy/claims/<id>`, so several jobs can share one
# list. Each claimed ID runs the frozen script exactly as simstudy.sh would,
# with SLURM_ARRAY_TASK_ID set to that ID, so file names, fits, chains, and
# seeds are those of the original array. Thread count changes the order of
# floating-point sums, so a chain is not bit-identical to a 1-thread run of
# the same seed; it is an equally valid draw. A claim whose chain file never
# appears (worker died, walltime) is listed at the end; remove the claim
# directory and rerun to retry it.
#
# FOOTPRINT. Counts against the study's 75% allowance (192 CPUs, 773,505 MB)
# together with 30875; see TODO.md for the throttle split.
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   sbatch --cpus-per-task=64 --mem=40G \
#       --export=ALL,SIMSTUDY_PIN_LIST=_data/simstudy/pinned-tasks.txt \
#       _scripts/simstudy-pinned.sh
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=3-00:00:00
#SBATCH --job-name=simstudy-pin
#SBATCH --chdir=/home2/lan68/plasmofit/plasmofit-test
#SBATCH --output=_data/simstudy/logs/pin-%j.out
#SBATCH --error=_data/simstudy/logs/pin-%j.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

set -u
LIST=${SIMSTUDY_PIN_LIST:?set SIMSTUDY_PIN_LIST}
CLAIMS=_data/simstudy/claims
mkdir -p "$CLAIMS"
JOB=${SLURM_JOB_ID:-local}

CPUSET=$(grep Cpus_allowed_list /proc/self/status | awk '{print $2}')
expand() { tr ',' '\n' <<< "$1" | awk -F- '{if (NF==1) print $1; else for (i=$1;i<=$2;i++) print i}'; }
ALLOWED=" $(expand "$CPUSET" | tr '\n' ' ') "

## Physical cores held (first hardware thread of each), tagged with L3 group.
CORES=()
for c in $(expand "$CPUSET" | sort -n); do
    first=$(sed 's/-/,/' /sys/devices/system/cpu/cpu$c/topology/thread_siblings_list | cut -d, -f1)
    [ "$c" -eq "$first" ] || continue
    g=$(cut -d, -f1 <<< "$(sed 's/-/,/' /sys/devices/system/cpu/cpu$c/cache/index3/shared_cpu_list)")
    CORES+=("$g:$c")
done

## Workers: whole L3 groups (4 held cores) first, then the leftover cores in
## group order, 4 at a time. Fewer than 4 left over are left idle.
declare -A BYG
for gc in "${CORES[@]}"; do BYG[${gc%%:*}]+="${gc##*:} "; done
WORKERS=(); LEFT=()
for g in $(printf '%s\n' "${!BYG[@]}" | sort -n); do
    set -- ${BYG[$g]}
    if [ $# -ge 4 ]; then WORKERS+=("$(IFS=,; echo "$1,$2,$3,$4") whole"); shift 4; fi
    LEFT+=("$@")
done
while [ ${#LEFT[@]} -ge 4 ]; do
    WORKERS+=("$(IFS=,; echo "${LEFT[*]:0:4}") split"); LEFT=("${LEFT[@]:4}")
done
echo "job ${JOB} cpuset ${CPUSET}; ${#CORES[@]} physical cores, ${#WORKERS[@]} workers; idle cores: ${LEFT[*]:-none}"
for w in "${WORKERS[@]}"; do echo "  worker cores ${w% *} (${w#* } L3 group)"; done
[ ${#WORKERS[@]} -ge 1 ] || { echo "ERROR: no worker can be formed" >&2; exit 1; }
[ "${SIMSTUDY_PIN_DRYRUN:-0}" = "1" ] && exit 0

worker() {
    local cores=$1 n=0
    for id in $(cat "$LIST"); do
        mkdir "$CLAIMS/$id" 2>/dev/null || continue
        echo "$(date '+%F %T') job ${JOB} cores ${cores}: task ${id}"
        SLURM_ARRAY_TASK_ID=$id SIMSTUDY_MODE=chain taskset -c "$cores" \
            /programs/R-4.6.1/bin/Rscript --vanilla _data/simstudy/simstudy-run.R \
            > "_data/simstudy/logs/${id}.out" 2> "_data/simstudy/logs/${id}.err" &
        local pid=$!
        if [ $n -eq 0 ]; then   # once per worker: confirm the thread count
            ( sleep 600; ps -o nlwp=,pcpu= -p $pid | awk -v c="$cores" -v t="$id" \
                '{print "task " t " on " c ": " $1 " threads, " $2 "% cpu after 10 min"}' ) &
        fi
        wait $pid || echo "$(date '+%F %T') task ${id} exited non-zero"
        n=$((n + 1))
    done
    echo "$(date '+%F %T') worker ${cores} done after ${n} tasks"
}

for w in "${WORKERS[@]}"; do worker "${w% *}" & done
wait

echo "claims without a finished chain, including any still running in other jobs (retry by removing the claim):"
for d in "$CLAIMS"/*; do
    id=$(basename "$d")
    grep -q "wrote" "_data/simstudy/logs/${id}.out" 2>/dev/null || echo "  $id"
done
