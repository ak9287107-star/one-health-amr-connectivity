#!/usr/bin/env bash
# Stage C: ClermonTyping (--minimal) for all E. coli cohort genomes. Resumable.
set -uo pipefail
BASE="$HOME/workspace/one_health_amr_pipeline/full_protocol_run"
export PATH="$HOME/workspace/tools/envs/oh_tools/bin:$HOME/workspace/tools/envs/oh_amr/bin:$PATH"
mkdir -p "$BASE/stageC/clermont"
WORK="$(mktemp -d)"; trap 'rm -rf "$WORK"' EXIT
cd "$WORK"
python3 - "$BASE/stageB/cohort_manifest.csv" <<'EOF' > "$WORK/ec_list.txt"
import csv, sys
for r in csv.DictReader(open(sys.argv[1])):
    if 'Klebsiella' not in r['taxgroup_name']:
        print(r['asm_acc'])
EOF
while read -r acc; do
  outphy="$BASE/stageC/clermont/${acc}_phylogroups.txt"
  [[ -s "$outphy" ]] && continue
  rm -rf "ct_$acc"
  if bash "$HOME/workspace/tools/clermontyping_src/clermonTyping.sh" --fasta "$BASE/stageB/genomes/${acc}.fna" --name "ct_$acc" --minimal > "$BASE/stageC/clermont/${acc}.log" 2>&1; then
    cp ct_$acc/*_phylogroups.txt "$outphy" 2>/dev/null || echo "$acc CLERMONT_OUTPUT_MISSING" >> "$BASE/stageC/failures.txt"
  else
    echo "$acc CLERMONT_FAIL" >> "$BASE/stageC/failures.txt"
  fi
  rm -rf "ct_$acc"
  echo "clermont $acc $(date +%H:%M:%S)" >> "$BASE/stageC/clermont_progress.log"
done < "$WORK/ec_list.txt"
echo CLERMONT_BATCH_COMPLETE >> "$BASE/stageC/clermont_progress.log"
