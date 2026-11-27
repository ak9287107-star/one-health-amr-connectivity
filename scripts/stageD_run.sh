#!/usr/bin/env bash
# Stage D: AMRFinderPlus (--plus) + MOB-suite mob_recon per cohort genome.
# Resumable: skips genomes whose outputs already exist and are non-empty.
# Usage: stageD_run.sh <genome_dir> <out_dir> [parallel_jobs] [threads_each]
set -uo pipefail
GDIR="${1:?genome dir}"
OUT="${2:?out dir}"
JOBS="${3:-2}"
THREADS="${4:-2}"
export PATH="$HOME/workspace/tools/envs/oh_amr/bin:$PATH"
mkdir -p "$OUT/amrfinder" "$OUT/mob_suite" "$OUT/logs"

run_one() {
  fna="$1"
  acc="$(basename "$fna" .fna)"
  amr_out="$OUT/amrfinder/${acc}_amrfinder.tsv"
  if [[ ! -s "$amr_out" ]]; then
    org=""
    # organism flag decided by caller via sidecar file org_map.tsv (acc \t org)
    org="$(awk -v a="$acc" '$1==a{print $2}' "$OUT/org_map.tsv" 2>/dev/null || true)"
    if [[ -n "$org" ]]; then
      amrfinder -n "$fna" -O "$org" --plus --threads "$THREADS" -o "$amr_out" > "$OUT/logs/${acc}_amrfinder.log" 2>&1 \
        || { echo "$acc AMRFINDER_FAIL" >> "$OUT/failures.txt"; }
    else
      amrfinder -n "$fna" --plus --threads "$THREADS" -o "$amr_out" > "$OUT/logs/${acc}_amrfinder.log" 2>&1 \
        || { echo "$acc AMRFINDER_FAIL" >> "$OUT/failures.txt"; }
    fi
  fi
  mob_dir="$OUT/mob_suite/${acc}"
  if [[ ! -s "$mob_dir/contig_report.txt" ]]; then
    mob_recon --infile "$fna" --outdir "$mob_dir" --num_threads "$THREADS" --force > "$OUT/logs/${acc}_mob.log" 2>&1 \
      || { echo "$acc MOB_FAIL" >> "$OUT/failures.txt"; }
  fi
  echo "done $acc $(date +%H:%M:%S)"
}
export OUT THREADS
export -f run_one

ls "$GDIR"/*.fna | xargs -P "$JOBS" -I{} bash -c 'run_one "$@"' _ {}
echo "STAGE_D_BATCH_COMPLETE amrfinder=$(ls "$OUT"/amrfinder/*_amrfinder.tsv 2>/dev/null | wc -l) mob=$(ls -d "$OUT"/mob_suite/*/ 2>/dev/null | wc -l)"
