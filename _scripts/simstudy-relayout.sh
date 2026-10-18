#!/bin/bash
# One-off, 2026-10-18: move the simulation study's pending fit n_c 384 chain
# tasks out of SLURM 30875 (packed single-threaded layout) into two pinned
# jobs (_scripts/simstudy-pinned.sh), inside the study's 75% allowance of 192
# CPUs, and make the merge job 30876 wait for them too. Run from a login
# shell; it prints and records everything it does in
# _data/simstudy/relayout-<date>.txt.
#
# BUDGET. 30875's throttle drops from 96 to 48 tasks (96 CPUs). The pinned
# jobs take 64 CPUs (8 workers) and 32 CPUs (4 workers). 30875 has R tasks
# running now and starts no new one until fewer than 48 run, so:
#   P1 may start once R - 64 of them have ended:  2 * 64 + 64 = 192
#   P2 may start once R - 48 have ended:          2 * 48 + 64 + 32 = 192
# Each pinned job depends (afterany) on that many running 30875 tasks,
# chosen as those expected to end soonest: fit n_c 192 chains first, then
# fit n_c 384, longest-running first within each.
#
#   cd /home2/lan68/plasmofit/plasmofit-test && bash _scripts/simstudy-relayout.sh

set -eu
cd /home2/lan68/plasmofit/plasmofit-test
LOG=_data/simstudy/relayout-$(date +%F).txt
exec > >(tee -a "$LOG") 2>&1
echo "== $(date '+%F %T')"

cls() { awk '{split($1,a,"_"); id=a[2]; f=int((id-1)/8)+1; j=(f-1)%48+1;
              print id, (j<=24?384:192), (j%2==1?"wide":"fixed"), $2}'; }
hrs() { awk '{n=split($4,p,/[-:]/); s=0;
              if (n==4) s=p[1]*86400+p[2]*3600+p[3]*60+p[4];
              else if (n==3) s=p[1]*3600+p[2]*60+p[3]; else s=p[1]*60+p[2];
              print $1, $2, $3, s}'; }

scontrol update JobId=30875 ArrayTaskThrottle=48
echo "30875 throttle: $(scontrol show job 30875 | grep -o 'ArrayTaskThrottle=[0-9]*' | head -1)"

LIST=_data/simstudy/pinned-tasks.txt
squeue -j 30875 -h -r -t PD -o "%i %M" | cls | awk '$2==384 && $3=="wide" {print $1}' | sort -n > "$LIST"
echo "pending fit n_c 384 tasks to move: $(wc -l < "$LIST") ($(head -1 "$LIST")-$(tail -1 "$LIST"))"
scancel $(sed 's/^/30875_/' "$LIST")
sleep 5
LEFT=$(squeue -j 30875 -h -r -t PD -o "%i" | sed 's/.*_//' | grep -Fxf "$LIST" | wc -l || true)
echo "still pending after scancel: $LEFT"
[ "$LEFT" -eq 0 ] || { echo "ERROR: some tasks were not cancelled; stopping before submitting" >&2; exit 1; }

## MOVE_RUNNING="id id ...": running fit n_c 384 tasks to cancel and re-run
## pinned as well (Lucas, 2026-10-18: the 16 that had run under ~1 h).
if [ -n "${MOVE_RUNNING:-}" ]; then
    echo "cancelling running tasks to re-run pinned: $MOVE_RUNNING"
    scancel $(tr ' ' '\n' <<< "$MOVE_RUNNING" | grep -v '^$' | sed 's/^/30875_/')
    tr ' ' '\n' <<< "$MOVE_RUNNING" | grep -v '^$' >> "$LIST"
    sort -n -u -o "$LIST" "$LIST"
    for i in $(seq 1 24); do
        STILL=$(squeue -j 30875 -h -r -o "%i" | sed 's/.*_//' | grep -Fxf <(tr ' ' '\n' <<< "$MOVE_RUNNING") | wc -l || true)
        [ "$STILL" -eq 0 ] && break; sleep 5
    done
    echo "of those, still in the queue: $STILL; list now $(wc -l < "$LIST") tasks"
    [ "$STILL" -eq 0 ] || { echo "ERROR: running tasks did not end; stopping before submitting" >&2; exit 1; }
fi

RUN=$(squeue -j 30875 -h -r -t R -o "%i %M" | cls | hrs | sort -k2,2n -k4,4nr | awk '{print $1}')
R=$(wc -w <<< "$RUN")
deps() { local n=$1; [ "$n" -gt 0 ] || return 0
         echo "--dependency=afterany:$(tr ' ' '\n' <<< "$RUN" | head -n "$n" | sed 's/^/30875_/' | paste -sd:)"; }
N1=$(( R > 64 ? R - 64 : 0 )); N2=$(( R > 48 ? R - 48 : 0 ))
echo "30875 running: $R; P1 waits for $N1 of them, P2 for $N2"

P1=$(sbatch --parsable --cpus-per-task=64 --mem=40G $(deps $N1) \
        --export=ALL,SIMSTUDY_PIN_LIST=$LIST _scripts/simstudy-pinned.sh)
P2=$(sbatch --parsable --cpus-per-task=32 --mem=20G $(deps $N2) \
        --export=ALL,SIMSTUDY_PIN_LIST=$LIST _scripts/simstudy-pinned.sh)
echo "pinned jobs: P1 $P1 (64 CPUs), P2 $P2 (32 CPUs)"

scontrol update JobId=30876 Dependency=afterany:30875:$P1:$P2
echo "30876 now: $(squeue -j 30876 -h -o '%E')"
squeue -u "$USER" -h -o "%i %T %C %m" | awk '$2=="RUNNING"{c+=$3} END{print "our running CPUs now:", c}'
