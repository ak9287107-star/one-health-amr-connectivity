#!/usr/bin/env python3
"""Stage E: connectivity analyses (Protocol S4 Stage 5).
Inputs (whatever exists at run time; Ns reported per analysis):
  - stageE/ska/distances_{ecoli,kpneumoniae}.tsv  (SKA2 pairwise SNP distances)
  - stageB/cohort_manifest.csv                    (metadata + NCBI AMRFinder genotypes)
  - stageC/all_mlst.tsv                           (MLST, if present)
  - stageD/mob_suite/<acc>/mobtyper_results.txt   (MOB-suite plasmid clusters)
  - pilot mob outputs for the sensitivity cohort (separate, labelled)
Outputs in stageE/: ska_cross_sector_pairs.tsv, bray_curtis summaries,
  plasmid_cluster_sharing.tsv, L1L4_scores.tsv, stageE_summary.md
Bray-Curtis here is computed directly (no scikit-bio dependency).
SKA distances are split-kmer approximations of core-SNP distances -- reported as such.
"""
import csv, os, re, sys, json
from collections import defaultdict, Counter
import numpy as np

BASE = os.path.expanduser('~/workspace/one_health_amr_pipeline/full_protocol_run')
OUT = os.path.join(BASE, 'stageE')
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(20261124)

man = {}
with open(os.path.join(BASE, 'stageB', 'cohort_manifest.csv'), newline='') as fh:
    for r in csv.DictReader(fh):
        man[r['asm_acc']] = r

def genes_of(r):
    return set(g for g in r['amr_genes'].split(';') if g)

notes = []
summary = ['# Stage E — Connectivity Results\n']

# ---------- 1. SKA cross-sector pairs ----------
ska_stats = {}
for grp, label in (('ecoli', 'E.coli and Shigella'), ('kpneumoniae', 'Klebsiella pneumoniae')):
    path = os.path.join(OUT, 'ska', f'distances_{grp}.tsv')
    if not os.path.exists(path):
        ska_stats[grp] = None
        notes.append(f'SKA distances for {grp}: NOT AVAILABLE at analysis time')
        continue
    pairs = []
    with open(path, newline='') as fh:
        rd = csv.reader(fh, delimiter='\t')
        header = next(rd)
        for row in rd:
            if len(row) < 3:
                continue
            a, b = row[0], row[1]
            try:
                snps = float(row[2])
            except ValueError:
                continue
            if a in man and b in man and snps <= 100:
                pairs.append({'acc1': a, 'acc2': b, 'ska_snps': snps,
                              'species': label,
                              'sector1': man[a]['one_health_sector'], 'sector2': man[b]['one_health_sector'],
                              'country1': man[a]['country'], 'country2': man[b]['country'],
                              'cross_sector': man[a]['one_health_sector'] != man[b]['one_health_sector'],
                              'same_country': man[a]['country'] == man[b]['country'],
                              'tier': '<=20' if snps <= 20 else '21-100'})
    ska_stats[grp] = pairs
    with open(os.path.join(OUT, f'ska_pairs_{grp}.tsv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, delimiter='\t', fieldnames=['acc1', 'acc2', 'species', 'ska_snps', 'tier', 'sector1', 'sector2', 'cross_sector', 'country1', 'country2', 'same_country'])
        w.writeheader(); w.writerows(pairs)
    n_cross = sum(p['cross_sector'] for p in pairs)
    n20 = sum(p['tier'] == '<=20' for p in pairs)
    n20_cross = sum(p['cross_sector'] and p['tier'] == '<=20' for p in pairs)
    summary.append(f'## SKA2 ({label}): {len(pairs)} pairs <=100 SNP; {n20} at <=20; '
                   f'cross-sector {n_cross} (<=20: {n20_cross})')

# ---------- 2. Bray-Curtis on ARG profiles (cohort, NCBI genotypes) ----------
accs = sorted(man)
all_genes = sorted(set().union(*[genes_of(man[a]) for a in accs])) if accs else []
gidx = {g: i for i, g in enumerate(all_genes)}
X = np.zeros((len(accs), len(all_genes)))
for i, a in enumerate(accs):
    for g in genes_of(man[a]):
        X[i, gidx[g]] = 1.0

def bray_curtis(u, v):
    d = np.abs(u - v).sum()
    s = (u + v).sum()
    return d / s if s > 0 else 0.0

# group-mean profiles per species x sector and species x income
groups = defaultdict(list)
for i, a in enumerate(accs):
    r = man[a]
    sp = 'Kp' if 'Klebsiella' in r['taxgroup_name'] else 'Ec'
    groups[(sp, r['one_health_sector'])].append(i)
    groups[(sp, r['income_group'])].append(i)
bc_rows = []
keys = sorted(groups)
for i in range(len(keys)):
    for j in range(i + 1, len(keys)):
        k1, k2 = keys[i], keys[j]
        if k1[0] != k2[0]:
            continue
        p1 = X[groups[k1]].mean(axis=0)
        p2 = X[groups[k2]].mean(axis=0)
        bc_rows.append({'species': k1[0], 'group1': k1[1], 'group2': k2[1],
                        'bray_curtis_ARG_profile': round(bray_curtis(p1, p2), 4),
                        'n1': len(groups[k1]), 'n2': len(groups[k2])})
with open(os.path.join(OUT, 'bray_curtis_ARG_groups.tsv'), 'w', newline='') as fh:
    w = csv.DictWriter(fh, delimiter='\t', fieldnames=['species', 'group1', 'group2', 'bray_curtis_ARG_profile', 'n1', 'n2'])
    w.writeheader(); w.writerows(bc_rows)

# per-genome pairwise BC: cross-sector vs within-sector (subsample pairs for speed)
pair_idx = [(i, j) for i in range(len(accs)) for j in range(i + 1, len(accs))]
if len(pair_idx) > 20000:
    sel = rng.choice(len(pair_idx), 20000, replace=False)
    pair_idx = [pair_idx[s] for s in sel]
bc_cross, bc_within = [], []
for i, j in pair_idx:
    a, b = accs[i], accs[j]
    if ('Klebsiella' in man[a]['taxgroup_name']) != ('Klebsiella' in man[b]['taxgroup_name']):
        continue
    v = bray_curtis(X[i], X[j])
    if man[a]['one_health_sector'] != man[b]['one_health_sector']:
        bc_cross.append(v)
    else:
        bc_within.append(v)
summary.append(f'## Bray-Curtis ARG profiles (cohort, per-genome pairs, n_cross={len(bc_cross)}, n_within={len(bc_within)}): '
               f'mean cross-sector {np.mean(bc_cross):.3f} vs within-sector {np.mean(bc_within):.3f}'
               if bc_cross and bc_within else '## Bray-Curtis: insufficient pairs')

# rarefaction: equal-N (min group size) subsampling of sector profiles, 100 reps
sector_groups = {k: v for k, v in groups.items() if k[1] in ('Human', 'Livestock_Poultry', 'Wastewater_Aquatic')}
rare_rows = []
for sp in ('Ec', 'Kp'):
    sg = {k[1]: v for k, v in sector_groups.items() if k[0] == sp}
    if len(sg) < 2:
        continue
    nmin = min(len(v) for v in sg.values())
    secs = sorted(sg)
    for x in range(len(secs)):
        for y in range(x + 1, len(secs)):
            vals = []
            for rep in range(100):
                i1 = rng.choice(sg[secs[x]], nmin, replace=False)
                i2 = rng.choice(sg[secs[y]], nmin, replace=False)
                vals.append(bray_curtis(X[i1].mean(axis=0), X[i2].mean(axis=0)))
            rare_rows.append({'species': sp, 'sector1': secs[x], 'sector2': secs[y], 'n_per_group': nmin,
                              'bc_mean': round(float(np.mean(vals)), 4), 'bc_sd': round(float(np.std(vals)), 4)})
with open(os.path.join(OUT, 'bray_curtis_rarefaction.tsv'), 'w', newline='') as fh:
    w = csv.DictWriter(fh, delimiter='\t', fieldnames=['species', 'sector1', 'sector2', 'n_per_group', 'bc_mean', 'bc_sd'])
    w.writeheader(); w.writerows(rare_rows)

# ---------- 3. Shared STs across sectors (MLST) ----------
st_by_acc = {}
mlst_path = os.path.join(BASE, 'stageC', 'all_mlst.tsv')
if os.path.exists(mlst_path):
    with open(mlst_path, newline='') as fh:
        for row in csv.reader(fh, delimiter='\t'):
            if len(row) >= 3:
                acc = os.path.basename(row[0]).replace('.fna', '')
                st_by_acc[acc] = (row[1], row[2])
st_sectors = defaultdict(set)
st_counts = Counter()
for a, (scheme, st) in st_by_acc.items():
    if a in man and st not in ('-', ''):
        key = (scheme, st)
        st_sectors[key].add(man[a]['one_health_sector'])
        st_counts[key] += 1
shared_sts = {k: v for k, v in st_sectors.items() if len(v) >= 2}
with open(os.path.join(OUT, 'shared_STs_across_sectors.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['scheme', 'ST', 'n_genomes', 'sectors'])
    for k in sorted(shared_sts, key=lambda k: -st_counts[k]):
        w.writerow([k[0], k[1], st_counts[k], ';'.join(sorted(shared_sts[k]))])
summary.append(f'## MLST: {len(st_by_acc)} genomes typed; {len(shared_sts)} STs found in >=2 sectors '
               f'of {len(st_sectors)} distinct STs')

# ---------- 4. MOB-suite plasmid cluster sharing (core cohort, completed subset) ----------
mob_dir = os.path.join(BASE, 'stageD', 'mob_suite')
cluster_sectors = defaultdict(set)
cluster_species = {}
n_mob = 0
if os.path.isdir(mob_dir):
    for acc in sorted(os.listdir(mob_dir)):
        mt = os.path.join(mob_dir, acc, 'mobtyper_results.txt')
        if os.path.exists(mt) and acc in man:
            n_mob += 1
            with open(mt, newline='') as fh:
                rd = csv.DictReader(fh, delimiter='\t')
                for row in rd:
                    cid = row.get('primary_cluster_id', '')
                    if cid and cid != '-':
                        cluster_sectors[cid].add(man[acc]['one_health_sector'])
                        cluster_species[cid] = 'Kp' if 'Klebsiella' in man[acc]['taxgroup_name'] else 'Ec'
shared_clusters = {c: s for c, s in cluster_sectors.items() if len(s) >= 2}
with open(os.path.join(OUT, 'plasmid_cluster_sharing.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['primary_cluster_id', 'species', 'sectors', 'n_sectors'])
    for c in sorted(shared_clusters):
        w.writerow([c, cluster_species.get(c, ''), ';'.join(sorted(shared_clusters[c])), len(shared_clusters[c])])
    w.writerow(['# MOB-processed genomes in core cohort', n_mob, '', ''])
summary.append(f'## MOB-suite plasmid clusters (core cohort subset, n={n_mob} genomes): '
               f'{len(cluster_sectors)} clusters, {len(shared_clusters)} in >=2 sectors')

# ---------- 5. L1-L4 connectivity scores (operationalized from protocol S2/S4) ----------
# Definitions (documented operationalization of Li et al. 2025 as described in the protocol):
#  L1 ecological: fraction of (country) cells with the species present in >=2 sectors (audit-level, per income group)
#  L2 microbial: shared-ST fraction + cross-sector SKA-pair fraction (sequence sharing)
#  L3 genetic: 1 - Bray-Curtis(ARG profiles) between sectors + shared plasmid-cluster fraction
#  L4 risk: critical-gene (blaNDM/blaKPC/blaOXA-48-like/mcr) carriage shared across sectors
# Composite connectivity = (2 * sequence_sharing + profile_similarity) / 3  (protocol: 2x weighting)
crit = re.compile(r'^(blaNDM|blaKPC|blaOXA-48|blaOXA-181|blaOXA-232|mcr-)')
l4_rows = []
for sp_label, sp_key in (('E.coli and Shigella', 'Ec'), ('Klebsiella pneumoniae', 'Kp')):
    for grp_name, filt in (('LMIC', lambda r: r['income_group'] == 'LMIC'),
                           ('HIC', lambda r: r['income_group'] == 'HIC')):
        sub = [a for a in accs if man[a]['taxgroup_name'] == sp_label and filt(man[a])]
        if len(sub) < 6:
            continue
        sec_genes = defaultdict(set)
        sec_crit = defaultdict(set)
        for a in sub:
            sec_genes[man[a]['one_health_sector']] |= genes_of(man[a])
            sec_crit[man[a]['one_health_sector']] |= {g for g in genes_of(man[a]) if crit.match(g)}
        secs = sorted(sec_genes)
        if len(secs) < 2:
            continue
        # L3: mean pairwise 1-BC between sector mean profiles
        sims = []
        for x in range(len(secs)):
            for y in range(x + 1, len(secs)):
                i1 = [accs.index(a) for a in sub if man[a]['one_health_sector'] == secs[x]]
                i2 = [accs.index(a) for a in sub if man[a]['one_health_sector'] == secs[y]]
                sims.append(1 - bray_curtis(X[i1].mean(axis=0), X[i2].mean(axis=0)))
        l3 = float(np.mean(sims)) if sims else 0.0
        # L4: fraction of critical genes (union) found in >=2 sectors
        all_crit = set().union(*sec_crit.values()) if sec_crit else set()
        shared_crit = {g for g in all_crit if sum(g in sec_crit[s] for s in secs) >= 2}
        l4 = len(shared_crit) / len(all_crit) if all_crit else 0.0
        # L2 components for this subset
        sts_sub = defaultdict(set)
        for a in sub:
            if a in st_by_acc and st_by_acc[a][1] not in ('-', ''):
                sts_sub[st_by_acc[a][1]].add(man[a]['one_health_sector'])
        l2_st = (sum(1 for v in sts_sub.values() if len(v) >= 2) / len(sts_sub)) if sts_sub else 0.0
        ska_pairs = ska_stats.get('ecoli' if sp_key == 'Ec' else 'kpneumoniae') or []
        ska_sub = [p for p in ska_pairs if p['acc1'] in set(sub) and p['acc2'] in set(sub)]
        l2_ska = (sum(p['cross_sector'] for p in ska_sub) / len(ska_sub)) if ska_sub else 0.0
        seq_share = float(np.mean([l2_st, l2_ska])) if (sts_sub or ska_sub) else 0.0
        composite = (2 * seq_share + l3) / 3
        l4_rows.append({'species': sp_label, 'income_group': grp_name, 'n_genomes': len(sub),
                        'L2_shared_ST_fraction': round(l2_st, 4), 'L2_cross_sector_SKA_fraction': round(l2_ska, 4),
                        'L3_ARG_profile_similarity': round(l3, 4), 'L4_critical_gene_sharing': round(l4, 4),
                        'composite_connectivity_2x_seq': round(composite, 4)})
with open(os.path.join(OUT, 'L1L4_scores.tsv'), 'w', newline='') as fh:
    if l4_rows:
        w = csv.DictWriter(fh, delimiter='\t', fieldnames=list(l4_rows[0].keys()))
        w.writeheader(); w.writerows(l4_rows)
summary.append('## L1-L4 scores: see L1L4_scores.tsv (L1 ecological co-occurrence is in the Stage A audit: '
               '65 country x species cells with >=2 sectors in the Primary Cohort)')

with open(os.path.join(OUT, 'stageE_summary.md'), 'w') as fh:
    fh.write('\n\n'.join(summary) + '\n\nNotes:\n' + '\n'.join('- ' + n for n in notes) + '\n')
print('\n\n'.join(summary))
print('NOTES:', notes)
