#!/bin/bash -l

# Run a wockner-fit.R config with each chain as its own process, pinned to one
# whole L3 cache group (4 physical cores, one thread per core), then merge.
#
# WHY. Benchmark 30846-30853 (claude/findings.md, "Thread-layout benchmark"):
# production runs 4 chains on 2 shared physical cores at ~83 ms per gradient;
# a chain alone in its L3 group with one thread per physical core runs at
# ~13-17 ms. A chain that shares its L3 group with another busy chain runs at
# half speed, and hyperthread siblings slow it, so each chain gets a whole
# group and the siblings are left idle.
#
# HOW. 40 logical CPUs = 20 physical cores, which always contain at least
# four whole 4-core L3 groups however SLURM aligns the allocation. Groups are
# read from /sys, not assumed. Each chain runs WOCKFIT_MODE=chain under
# `taskset` with WOCKFIT_THREADS=4; once all four have saved their files,
# WOCKFIT_MODE=merge combines them and writes the usual four outputs
# (wock-fit-, wock-data-, wock-fit-LOO-, wock-fit-RES-<cfg>.rds). Chain files
# wock-chain-<cfg>-<k>.rds stay on disk, so a failed merge reruns in minutes:
#   cd _data && WOCKFIT_MODE=merge SLURM_ARRAY_TASK_ID=<i> SLURM_CPUS_PER_TASK=4 \
#       WOCKFIT_SEED=<seed> WOCKFIT_SUFFIX=<suffix> Rscript --vanilla ../_scripts/wockner-fit.R
#
# FOOTPRINT. 40 CPUs and 32G per task: at most 3 at once inside the
# 128 CPU / 515,670 MB half-node budget (M = floor(min(128/40, 515670/32768)) = 3).
# Pass --array=<i>%3 or fewer for several configs.
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   sbatch --array=43 --export=ALL,WOCKFIT_SEED=20261010,WOCKFIT_SUFFIX=-pinned \
#       _scripts/wockner-fit-pinned.sh

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=32G
#SBATCH --time=3-00:00:00
#SBATCH --job-name=wock-pin
#SBATCH --chdir=/home2/lan68/plasmofit/plasmofit-test/_data
#SBATCH --output=wock-pin-%a.out
#SBATCH --error=wock-pin-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

set -u
N_CHAINS=4
CPUSET=$(grep Cpus_allowed_list /proc/self/status | awk '{print $2}')
expand() { tr ',' '\n' <<< "$1" | awk -F- '{if (NF==1) print $1; else for (i=$1;i<=$2;i++) print i}'; }
ALLOWED=" $(expand "$CPUSET" | tr '\n' ' ') "

## One line per whole L3 group inside the allocation: its physical cores
## (first hardware thread of each), comma-separated.
L3GROUPS=()
declare -A SEEN
for c in $(expand "$CPUSET" | sort -n); do
    shared=$(cat /sys/devices/system/cpu/cpu$c/cache/index3/shared_cpu_list)
    [ -n "${SEEN[$shared]:-}" ] && continue
    SEEN[$shared]=1
    ## physical cores of this group = members that are the lower sibling
    phys=()
    whole=1
    for m in $(expand "$shared"); do
        first=$(cut -d, -f1 <<< "$(sed 's/-/,/' /sys/devices/system/cpu/cpu$m/topology/thread_siblings_list)")
        [ "$m" -eq "$first" ] || continue
        phys+=("$m")
        [[ "$ALLOWED" == *" $m "* ]] || whole=0
    done
    [ "$whole" -eq 1 ] && L3GROUPS+=("$(IFS=,; echo "${phys[*]}")")
done

echo "task ${SLURM_ARRAY_TASK_ID} cpuset ${CPUSET}; whole L3 groups: ${L3GROUPS[*]}"
if [ "${#L3GROUPS[@]}" -lt "$N_CHAINS" ]; then
    echo "ERROR: only ${#L3GROUPS[@]} whole L3 groups in the allocation, need ${N_CHAINS}" >&2
    exit 1
fi

if [ "${WOCKFIT_PIN_DRYRUN:-0}" = "1" ]; then
    for k in $(seq 1 "$N_CHAINS"); do echo "dry run: chain $k -> cores ${L3GROUPS[$((k - 1))]}"; done
    exit 0
fi

PIDS=()
for k in $(seq 1 "$N_CHAINS"); do
    cores=${L3GROUPS[$((k - 1))]}
    nthr=$(tr ',' '\n' <<< "$cores" | wc -l)
    echo "chain $k -> cores $cores ($nthr threads)"
    WOCKFIT_MODE=chain WOCKFIT_CHAIN=$k WOCKFIT_THREADS=$nthr \
        taskset -c "$cores" Rscript --vanilla ../_scripts/wockner-fit.R \
        > "wock-pin-${SLURM_ARRAY_TASK_ID}-chain${k}.log" 2>&1 &
    PIDS+=($!)
done

FAIL=0
for k in $(seq 1 "$N_CHAINS"); do
    wait "${PIDS[$((k - 1))]}" || { echo "chain $k exited non-zero" >&2; FAIL=1; }
done
[ "$FAIL" -eq 0 ] || { echo "ERROR: a chain failed; see wock-pin-${SLURM_ARRAY_TASK_ID}-chain*.log" >&2; exit 1; }

echo "all chains done; merging"
WOCKFIT_MODE=merge Rscript --vanilla ../_scripts/wockner-fit.R
