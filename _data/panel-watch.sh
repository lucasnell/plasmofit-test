#!/bin/bash
# Scaffolding: waits for SLURM 28931 to drain, then runs the thread 2/3 panel.
# Deleted once its log reads clean; the panel's own output is _data/prior-panel.log.
cd /home2/lan68/plasmofit/plasmofit-test || exit 1
echo "=== waiting for SLURM 28931 (wock-fit tasks 10-13) ==="
while squeue -u lan68 -h -j 28931 2>/dev/null | grep -q .; do sleep 60; done
echo "=== 28931 drained at $(date) ==="
sacct -j 28931 --format=JobID,State,Elapsed,ExitCode -P | grep -v '\.batch'
echo
echo "=== fits on disk ==="
for c in no_pool pooled_cl np_anchor np_wide_total0 np_wide_bshape \
         np_wide_both pl_wide_total0 pl_wide_both; do
    f="_data/wock-fit-${c}.rds"
    if [ -f "$f" ]; then printf "  PRESENT %-16s %s\n" "$c" "$(stat -c %y "$f" | cut -d. -f1)"
    else printf "  MISSING %-16s\n" "$c"; fi
done
echo
echo "=== running wockner-prior-panel.R ==="
R_ENVIRON_USER=/dev/null srun -N 1 -n 1 -c 4 --mem=48G Rscript --vanilla \
    _scripts/wockner-prior-panel.R > _data/prior-panel.log 2>&1
echo "panel exited $? -- output in _data/prior-panel.log"
