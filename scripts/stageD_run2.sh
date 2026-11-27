#!/usr/bin/env bash
# Stage D v2: sequential AMRFinderPlus (--plus) + MOB-recon, one genome at a time
# (all 4 threads to a single genome; measured faster than parallel split on this
# 1.5-effective-CPU VM). Processing order comes from an ordered accession list
# (stratified round-robin, seed-fixed) so any completed prefix is an approximately
# balanced subsample. Resumable via skip logic.
# Usage: stageD_run2.sh <ordered_acc_list> <genome_dir> <out_dir>
set -uo pipefail
LIST="${1:?ordered accession list}"
GDIR="${2:?genome dir}"
OUT="${3:?out dir}"
export PATH="$HOME/workspace/tools/envs/oh_amr/bin:$PATH"
mkdir -p "$OUT/amrfinder" "$OUT/mob_suite" "$OUT/logs"
: >> "$OUT/failures.txt"

while read -r acc; do
  [[ -n "$acc" ]] || continue
  fna="$GDIR/${acc}.fna"
  [[ -s "$fna" ]] || { echo "$acc NO_FASTA" >> "$OUT/failures.txt"; continue; }
  t0=$(date +%s)
  amr_out="$OUT/amrfinder/${acc}_amrfinder.tsv"
  if [[ ! -s "$amr_out" ]]; then
    org="$(awk -v a="$acc" '$1==a{print $2}' "$OUT/org_map.tsv" 2>/dev/null || true)"
    if [[ -n "$org" ]]; then
      amrfinder -n "$fna" -O "$org" --plus --threads 4 -o "$amr_out" > "$OUT/logs/${acc}_amrfinder.log" 2>&1 \
        || echo "$acc AMRFINDER_FAIL" >> "$OUT/failures.txt"
    else
      amrfinder -n "$fna" --plus --threads 4 -o "$amr_out" > "$OUT/logs/${acc}_amrfinder.log" 2>&1 \
        || echo "$acc AMRFINDER_FAIL" >> "$OUT/failures.txt"
    fi
  fi
  mob_dir="$OUT/mob_suite/${acc}"
  if [[ ! -s "$mob_dir/contig_report.txt" ]]; then
    mob_recon --infile "$fna" --outdir "$mob_dir" --num_threads 4 --force > "$OUT/logs/${acc}_mob.log" 2>&1 \
      || echo "$acc MOB_FAIL" >> "$OUT/failures.txt"
  fi
  echo "done $acc $(( $(date +%s) - t0 ))s $(date +%H:%M:%S)" >> "$OUT/progress.log"
done < "$LIST"
echo "STAGE_D_V2_COMPLETE" >> "$OUT/progress.log"
