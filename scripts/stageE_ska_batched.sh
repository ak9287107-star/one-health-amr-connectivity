#!/usr/bin/env bash
# Stage E SKA, batched: whole-cohort `ska build` was OOM-killed on this 7.7 GB VM
# (93 E. coli genomes in one build). Build in chunks of 20 (1 thread), merge, distance.
set -uo pipefail
BASE="$HOME/workspace/one_health_amr_pipeline/full_protocol_run"
export PATH="$HOME/workspace/tools/envs/oh_amr/bin:$PATH"
cd "$BASE/stageE/ska"
CHUNK=20

for grp in ecoli kpneumoniae; do
  [[ -s "distances_${grp}.tsv" ]] && { echo "$grp already done"; continue; }
  rm -rf "chunks_$grp"; mkdir -p "chunks_$grp"
  split -l $CHUNK -d "samples_${grp}.tsv" "chunks_$grp/part_"
  skfs=()
  for part in chunks_$grp/part_*; do
    name="chunks_$grp/$(basename "$part")"
    if [[ ! -s "${name}.skf" ]]; then
      ska build -f "$part" -o "$name" --threads 1 > "${name}.log" 2>&1 || { echo "BUILD_FAIL $name"; continue; }
    fi
    skfs+=("${name}.skf")
  done
  echo "$grp chunks built: ${#skfs[@]}"
  if [[ ${#skfs[@]} -gt 1 ]]; then
    ska merge -o "genomes_${grp}" "${skfs[@]}" > "merge_${grp}.log" 2>&1
  elif [[ ${#skfs[@]} -eq 1 ]]; then
    cp "${skfs[0]}" "genomes_${grp}.skf"
  fi
  if [[ -s "genomes_${grp}.skf" ]]; then
    ska distance "genomes_${grp}.skf" --threads 2 > "distances_${grp}.tsv" 2> "ska_distance_${grp}.log"
    echo "$grp distances: $(grep -c "" "distances_${grp}.tsv") lines"
  fi
done
echo SKA_BATCHED_DONE
