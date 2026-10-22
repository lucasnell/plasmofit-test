#!/bin/bash
# Runs _scripts/loadbalance-test.R (header there): parts full1 and unit on one
# physical core, full4 on four physical cores, all pinned with taskset, hyper-
# thread siblings idle. Asks for 8 CPUs (4 physical cores) and reports whether
# they form one L3 group; full4 on a split group is slower and is flagged.
# Added 2026-10-22.
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   sbatch _scripts/loadbalance-test.sh
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --job-name=loadbal
#SBATCH --chdir=/home2/lan68/plasmofit/plasmofit-test
#SBATCH --output=_data/loadbalance-%j.out
#SBATCH --error=_data/loadbalance-%j.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

set -u
expand() { tr ',' '\n' <<< "$1" | awk -F- '{if (NF==1) print $1; else for (i=$1;i<=$2;i++) print i}'; }
CPUSET=$(grep Cpus_allowed_list /proc/self/status | awk '{print $2}')
CORES=(); GROUPS_SEEN=()
for c in $(expand "$CPUSET" | sort -n); do
    first=$(sed 's/-/,/' /sys/devices/system/cpu/cpu$c/topology/thread_siblings_list | cut -d, -f1)
    [ "$c" -eq "$first" ] || continue
    CORES+=("$c")
    GROUPS_SEEN+=("$(cut -d, -f1 <<< "$(sed 's/-/,/' /sys/devices/system/cpu/cpu$c/cache/index3/shared_cpu_list)")")
done
NG=$(printf '%s\n' "${GROUPS_SEEN[@]:0:4}" | sort -u | wc -l)
FOUR=$(IFS=,; echo "${CORES[*]:0:4}")
D=$(date +%F)
echo "cpuset $CPUSET; physical cores ${CORES[*]}; first four ($FOUR) span $NG L3 group(s)"
R=/programs/R-4.6.1/bin/Rscript
taskset -c "${CORES[0]}" $R --vanilla _scripts/loadbalance-test.R full1 > _data/loadbalance-full1-$D.txt 2>&1
taskset -c "${CORES[0]}" $R --vanilla _scripts/loadbalance-test.R unit  > _data/loadbalance-unit-$D.txt 2>&1
{ echo "4 physical cores $FOUR span $NG L3 group(s) (1 = one whole group)"
  taskset -c "$FOUR" $R --vanilla _scripts/loadbalance-test.R full4; } > _data/loadbalance-full4-$D.txt 2>&1
echo done
