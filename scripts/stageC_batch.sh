#!/usr/bin/env bash
# Stage C batch: ClermonTyping for E. coli, Kleborate for K. pneumoniae.
# Resumable via skip logic. Usage: stageC_batch.sh <cohort_manifest> <genome_dir> <out_dir>
set -uo pipefail
MANIFEST="$(readlink -f "${1:?}")"; GDIR="$(readlink -f "${2:?}")"; OUT="$(readlink -f "${3:?}")"
mkdir -p "$OUT/clermont" "$OUT/kleborate"
export PATH="$HOME/workspace/tools/envs/oh_tools/bin:$HOME/workspace/tools/envs/oh_amr/bin:$PATH"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
cd "$WORK"

python3 - "$MANIFEST" <<'EOF' > "$WORK/lists.txt"
import csv, sys
for r in csv.DictReader(open(sys.argv[1])):
    sp = 'kp' if 'Klebsiella' in r['taxgroup_name'] else 'ec'
    print(f"{sp}\t{r['asm_acc']}")
EOF

while IFS=$'\t' read -r sp acc; do
  fna="$GDIR/${acc}.fna"
  [[ -s "$fna" ]] || continue
  if [[ "$sp" == "ec" ]]; then
    outphy="$OUT/clermont/${acc}_phylogroups.txt"
    [[ -s "$outphy" ]] && continue
    rm -rf "ct_$acc"
    if bash "$HOME/workspace/tools/clermontyping_src/clermonTyping.sh" --fasta "$fna" --name "ct_$acc" --minimal > "$OUT/clermont/${acc}.log" 2>&1; then
      cp "ct_$acc/ct_${acc}_phylogroups.txt" "$outphy" 2>/dev/null || cp ct_$acc/*_phylogroups.txt "$outphy" 2>/dev/null || true
    fi
    rm -rf "ct_$acc"
  else
    outdir="$OUT/kleborate/$acc"
    if [[ ! -s "$outdir/kleborate.tsv" && ! -s "$outdir"/*kleborate*.tsv ]]; then
      kleborate -a "$fna" -o "$outdir" -p kpsc -t 2 > "$OUT/kleborate/${acc}.log" 2>&1 \
        || echo "$acc KLEBORATE_FAIL" >> "$OUT/failures.txt"
    fi
  fi
  echo "typed $acc $(date +%H:%M:%S)" >> "$OUT/progress.log"
done < "$WORK/lists.txt"
echo "STAGE_C_BATCH_COMPLETE" >> "$OUT/progress.log"
