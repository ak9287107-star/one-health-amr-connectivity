#!/usr/bin/env bash
# Stage C: Kleborate (kpsc preset) for all K. pneumoniae cohort genomes. Resumable.
set -uo pipefail
BASE="$HOME/workspace/one_health_amr_pipeline/full_protocol_run"
export PATH="$HOME/workspace/tools/envs/oh_tools/bin:$HOME/workspace/tools/envs/oh_amr/bin:$PATH"
mkdir -p "$BASE/stageC/kleborate"
python3 - "$BASE/stageB/cohort_manifest.csv" <<'EOF' > "$BASE/stageC/kp_list.txt"
import csv, sys
for r in csv.DictReader(open(sys.argv[1])):
    if 'Klebsiella' in r['taxgroup_name']:
        print(r['asm_acc'])
EOF
while read -r acc; do
  outdir="$BASE/stageC/kleborate/$acc"
  if ls "$outdir"/klebsiella_pneumo_complex_output.txt >/dev/null 2>&1; then continue; fi
  kleborate -a "$BASE/stageB/genomes/${acc}.fna" -o "$outdir" -p kpsc -t 2 > "$BASE/stageC/kleborate/${acc}.log" 2>&1 \
    || echo "$acc KLEBORATE_FAIL" >> "$BASE/stageC/failures.txt"
  echo "kleborate $acc $(date +%H:%M:%S)" >> "$BASE/stageC/kleborate_progress.log"
done < "$BASE/stageC/kp_list.txt"
echo KLEBORATE_BATCH_COMPLETE >> "$BASE/stageC/kleborate_progress.log"
