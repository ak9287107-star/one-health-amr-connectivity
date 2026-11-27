#!/usr/bin/env bash
# Stage E (SKA part): per-species SKA2 build + pairwise distances for the core cohort.
set -uo pipefail
BASE="$HOME/workspace/one_health_amr_pipeline/full_protocol_run"
export PATH="$HOME/workspace/tools/envs/oh_amr/bin:$PATH"
mkdir -p "$BASE/stageE/ska"
cd "$BASE/stageE/ska"

python3 - "$BASE/stageB/cohort_manifest.csv" "$BASE/stageB/genomes" <<'EOF'
import csv, os, sys
man, gdir = sys.argv[1], sys.argv[2]
groups = {}
for r in csv.DictReader(open(man)):
    key = 'kpneumoniae' if 'Klebsiella' in r['taxgroup_name'] else 'ecoli'
    fna = os.path.join(gdir, r['asm_acc'] + '.fna')
    if os.path.exists(fna):
        groups.setdefault(key, []).append((r['asm_acc'], os.path.abspath(fna)))
for key, rows in groups.items():
    with open(f'samples_{key}.tsv', 'w') as fh:
        for acc, fna in sorted(rows):
            fh.write(f'{acc}\t{fna}\n')
    print(key, len(rows))
EOF

for grp in ecoli kpneumoniae; do
  if [[ ! -s "genomes_${grp}.skf" ]]; then
    ska build -f "samples_${grp}.tsv" -o "genomes_${grp}" --threads 2 > "ska_build_${grp}.log" 2>&1
  fi
  if [[ -s "genomes_${grp}.skf" && ! -s "distances_${grp}.tsv" ]]; then
    ska distance "genomes_${grp}.skf" --threads 2 > "distances_${grp}.tsv" 2> "ska_distance_${grp}.log"
  fi
  echo "ska $grp done $(date +%H:%M:%S)"
done
echo SKA_ALL_DONE
