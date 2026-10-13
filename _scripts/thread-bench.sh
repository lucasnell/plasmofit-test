#!/bin/bash -l

# Launch BENCH_CHAINS copies of _scripts/thread-bench.R at once inside one
# SLURM task, to reproduce or vary a production task's layout.
#
# BENCH_PIN controls which of the task's CPUs each chain process may use.
# The thread pool sizes itself to the CPUs a process is allowed, and neither
# STAN_NUM_THREADS nor RcppParallel::setThreadOptions changed that
# (probe 30843-30845), so pinning is how threads per chain are set:
#   none   every chain sees all the task's CPUs (production as it runs now)
#   split  physical cores divided evenly, each chain gets its cores AND their
#          hyperthread siblings
#   phys   physical cores divided evenly, first hardware thread only
# On cbsugreischar logical CPU c and c+128 are the two threads of one core.
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   sbatch --cpus-per-task=16 \
#       --export=ALL,BENCH_KNOB=default,BENCH_PIN=split,BENCH_CHAINS=4,BENCH_TAG=B \
#       _scripts/thread-bench.sh
#
# Output lands in _data/thread-bench-<jobid>.out, one line per chain.

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --job-name=thread-bench
#SBATCH --chdir=/home2/lan68/plasmofit/plasmofit-test
#SBATCH --output=_data/thread-bench-%j.out
#SBATCH --error=_data/thread-bench-%j.err

CHAINS=${BENCH_CHAINS:-4}
PIN=${BENCH_PIN:-none}
CPUSET=$(grep Cpus_allowed_list /proc/self/status | awk '{print $2}')

## Expand "2-3,130-131" to one CPU per line, keep physical cores (< 128).
expand() { tr ',' '\n' <<< "$1" | awk -F- '{if (NF==1) print $1; else for (i=$1;i<=$2;i++) print i}'; }
PHYS=($(expand "$CPUSET" | awk '$1 < 128' | sort -n))
PER=$(( ${#PHYS[@]} / CHAINS ))

echo "# layout tag=${BENCH_TAG} pin=${PIN} chains=${CHAINS} cpus=${SLURM_CPUS_PER_TASK} cpuset=${CPUSET} phys_cores=${#PHYS[@]} per_chain=${PER}"
if [ "$PIN" != "none" ] && [ "$PER" -lt 1 ]; then echo "# ERROR fewer physical cores than chains"; exit 1; fi

for i in $(seq 0 $((CHAINS - 1))); do
    CMD="Rscript --vanilla _scripts/thread-bench.R"
    if [ "$PIN" != "none" ]; then
        MINE=("${PHYS[@]:$((i * PER)):$PER}")
        LIST=$(IFS=,; echo "${MINE[*]}")
        if [ "$PIN" = "split" ]; then
            SIB=$(for c in "${MINE[@]}"; do echo $((c + 128)); done | paste -sd,)
            LIST="${LIST},${SIB}"
        fi
        CMD="taskset -c ${LIST} ${CMD}"
        echo "# chain $((i + 1)) cpus ${LIST}"
    fi
    BENCH_DRAW=$(( (i + 1) * 60 )) $CMD &
done
wait
