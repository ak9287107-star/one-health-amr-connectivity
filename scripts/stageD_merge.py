#!/usr/bin/env python3
"""Stage D merge: AMRFinderPlus hits x MOB-suite contig calls -> localization table
+ Plasmid Localization Index (mirrors run_amr_connectivity_pipeline.sh stage 5).
Runs on the subset of cohort genomes having BOTH outputs; N is printed and written.
"""
import csv, glob, os, sys
import pandas as pd

BASE = os.path.expanduser('~/workspace/one_health_amr_pipeline/full_protocol_run')
man = pd.read_csv(os.path.join(BASE, 'stageB', 'cohort_manifest.csv'))
amr_dir, mob_dir = os.path.join(BASE, 'stageD', 'amrfinder'), os.path.join(BASE, 'stageD', 'mob_suite')
sum_dir = os.path.join(BASE, 'stageD', 'summary')
os.makedirs(sum_dir, exist_ok=True)

records, used = [], []
for amr_file in sorted(glob.glob(os.path.join(amr_dir, '*_amrfinder.tsv'))):
    acc = os.path.basename(amr_file).replace('_amrfinder.tsv', '')
    contig_rep = os.path.join(mob_dir, acc, 'contig_report.txt')
    mob_typer = os.path.join(mob_dir, acc, 'mobtyper_results.txt')
    if not os.path.exists(contig_rep):
        continue
    amr_df = pd.read_csv(amr_file, sep='\t')
    if amr_df.empty:
        continue
    mob_c = pd.read_csv(contig_rep, sep='\t')
    mob_c['contig_key'] = mob_c['contig_id'].astype(str).str.split().str[0]
    amr_df['contig_key'] = amr_df['Contig id'].astype(str).str.split().str[0]
    merged = amr_df.merge(mob_c[['contig_key', 'molecule_type', 'primary_cluster_id', 'secondary_cluster_id']],
                          on='contig_key', how='left')
    if os.path.exists(mob_typer):
        mt = pd.read_csv(mob_typer, sep='\t')
        if 'primary_cluster_id' in mt.columns:
            keep = [c for c in ['primary_cluster_id', 'rep_type(s)', 'predicted_mobility', 'size'] if c in mt.columns]
            merged = merged.merge(mt[keep].drop_duplicates(subset=['primary_cluster_id']),
                                  on='primary_cluster_id', how='left')
    merged['asm_acc'] = acc
    records.append(merged)
    used.append(acc)

print(f'genomes with both AMRFinder+MOB: {len(used)} -> {used}')
if not records:
    print('No merged records yet.')
    sys.exit(0)
full = pd.concat(records, ignore_index=True)
meta_cols = [c for c in ['asm_acc', 'country', 'income_group', 'one_health_sector', 'collection_year', 'taxgroup_name'] if c in man.columns]
full = full.merge(man[meta_cols].drop_duplicates(subset=['asm_acc']), on='asm_acc', how='left')
out = os.path.join(sum_dir, 'amr_plasmid_localization_table.tsv')
full.to_csv(out, sep='\t', index=False)
gene_col = 'Element symbol' if 'Element symbol' in full.columns else 'Gene symbol'
summary = (full.groupby(['country', 'one_health_sector', gene_col], dropna=False)
           .agg(n_hits=(gene_col, 'count'),
                n_plasmid=('molecule_type', lambda s: (s == 'plasmid').sum()),
                n_chromosome=('molecule_type', lambda s: (s == 'chromosome').sum()),
                unique_plasmid_clusters=('primary_cluster_id', lambda s: s.dropna().nunique()))
           .reset_index())
summary['plasmid_localization_index'] = (summary['n_plasmid'] / summary['n_hits']).round(4)
summary.to_csv(os.path.join(sum_dir, 'plasmid_localization_index_summary.tsv'), sep='\t', index=False)
print(f'wrote {out}: {len(full)} hits from {len(used)} genomes')
print(full['molecule_type'].value_counts(dropna=False).to_string())
