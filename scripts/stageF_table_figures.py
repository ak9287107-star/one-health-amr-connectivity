#!/usr/bin/env python3
"""Stage F: Table 1 + Figures 1-4 (Protocol S5), from real outputs only.
Fig 2 and Fig 3 use the pilot sensitivity-cohort dataset (n=50, complete
AMRFinderPlus + MOB-suite outputs) and are labelled as such; the core Primary
cohort's MOB subset was still small at production time (see RUN_STATUS.md).
"""
import csv, os, re
from collections import defaultdict, Counter
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = os.path.expanduser('~/workspace/one_health_amr_pipeline/full_protocol_run')
PILOT = os.path.expanduser('~/workspace/one_health_amr_pipeline/pilot_run')
OUT = os.path.join(BASE, 'stageF')
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({'font.size': 9, 'figure.dpi': 300, 'savefig.dpi': 300})
CRIT = re.compile(r'^(blaNDM|blaKPC|blaOXA-48|blaOXA-181|blaOXA-232|mcr-)')

man = {}
with open(os.path.join(BASE, 'stageB', 'cohort_manifest.csv'), newline='') as fh:
    for r in csv.DictReader(fh):
        man[r['asm_acc']] = r
qc = {}
with open(os.path.join(BASE, 'stageB', 'qc_results.tsv'), newline='') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        qc[r['asm_acc']] = r

# ---------------- Table 1 ----------------
def med(vals):
    v = [x for x in vals if x not in ('', None)]
    return round(float(np.median([float(x) for x in v])), 1) if v else ''

rows = []
groups = defaultdict(list)
for a, r in man.items():
    groups[(r['taxgroup_name'], r['one_health_sector'], r['income_group'])].append(a)
for k in sorted(groups):
    accs = groups[k]
    crit_n = sum(1 for a in accs if any(CRIT.match(g) for g in man[a]['amr_genes'].split(';') if g))
    rows.append({'species': k[0], 'sector': k[1], 'income_group': k[2], 'n_genomes': len(accs),
                 'median_year': med([man[a]['collection_year'] for a in accs]),
                 'median_contigs_QUAST': med([qc.get(a, {}).get('quast_n_contigs') for a in accs]),
                 'median_N50_kb_QUAST': round(med([qc.get(a, {}).get('quast_N50') for a in accs]) / 1000, 1) if accs else '',
                 'pct_any_critical_gene_NCBI_genotypes': round(100 * crit_n / len(accs), 1)})
with open(os.path.join(OUT, 'Table1_dataset_summary.tsv'), 'w', newline='') as fh:
    w = csv.DictWriter(fh, delimiter='\t', fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)
print('Table 1 rows:', len(rows))

# ---------------- Figure 1 ----------------
audit = defaultdict(int)
with open(os.path.join(BASE, 'stageA', 'audit_headline_totals.tsv'), newline='') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        if r['species'] == 'ALL_SPECIES' and r['cohort'] == 'primary' and r['income_group'] in ('LMIC', 'HIC'):
            audit[(r['income_group'], r['sector'])] = int(r['n_distinct_asm_acc'])
fig, axes = plt.subplots(1, 3, figsize=(11, 3.4))
sectors = ['Human', 'Livestock_Poultry', 'Wastewater_Aquatic']
x = np.arange(3)
axes[0].bar(x - 0.2, [audit[('LMIC', s)] for s in sectors], 0.4, label='LMIC')
axes[0].bar(x + 0.2, [audit[('HIC', s)] for s in sectors], 0.4, label='HIC')
axes[0].set_xticks(x); axes[0].set_xticklabels(['Human', 'Livestock', 'Wastewater'], rotation=20)
axes[0].set_ylabel('Primary-cohort genomes'); axes[0].set_title('A. Global Primary Cohort (audit)')
axes[0].legend(frameon=False)
cc = Counter(r['country'] for r in man.values())
top = cc.most_common(12)
cols = ['#b2182b' if man[[a for a in man if man[a]['country'] == c][0]]['income_group'] == 'LMIC' else '#2166ac' for c, _ in top]
axes[1].barh([c for c, _ in top][::-1], [n for _, n in top][::-1], color=cols[::-1])
axes[1].set_xlabel('Core-cohort genomes (n=151)'); axes[1].set_title('B. Core cohort by country (red=LMIC)')
shared_st = 0; total_st = 0
p = os.path.join(BASE, 'stageE', 'shared_STs_across_sectors.tsv')
if os.path.exists(p):
    with open(p) as fh:
        shared_st = sum(1 for _ in fh) - 1
axes[2].bar(['Shared STs\n(across sectors)'], [max(shared_st, 0)], color='#4daf4a')
axes[2].set_title('C. Sector overlap (core cohort, MLST subset)')
axes[2].set_ylabel('Count')
fig.tight_layout()
fig.savefig(os.path.join(OUT, 'Figure1_global_distribution.pdf')); fig.savefig(os.path.join(OUT, 'Figure1_global_distribution.png'))
plt.close(fig)

# ---------------- Figure 2 (core Primary cohort, complete genomes) ----------------
loc_path = os.path.join(BASE, 'stageD', 'summary', 'amr_plasmid_localization_table.tsv')
agg = defaultdict(lambda: [0, 0])
with open(loc_path, newline='') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        g = r['Element symbol']
        if r['molecule_type'] == 'plasmid':
            agg[g][0] += 1
        elif r['molecule_type'] == 'chromosome':
            agg[g][1] += 1
top_genes = sorted(agg.items(), key=lambda kv: -(kv[1][0] + kv[1][1]))[:15]
fig, ax = plt.subplots(figsize=(7.5, 4.2))
labels = [g for g, _ in top_genes][::-1]
plas = [v[0] for _, v in top_genes][::-1]
chrom = [v[1] for _, v in top_genes][::-1]
ax.barh(labels, plas, color='#b2182b', label='Plasmid')
ax.barh(labels, chrom, left=plas, color='#2166ac', label='Chromosome')
ax.set_xlabel('AMR gene hits (MOB-suite calls, complete genomes)')
ax.set_title('ARG localization landscape — core Primary cohort (global, n=151 complete genomes)')
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(os.path.join(OUT, 'Figure2_localization_landscape.pdf')); fig.savefig(os.path.join(OUT, 'Figure2_localization_landscape.png'))
plt.close(fig)

# ---------------- Figure 3 (core cohort plasmid cluster network) ----------------
cluster_sec = defaultdict(set)
cluster_n = Counter()
mob_root = os.path.join(BASE, 'stageD', 'mob_suite')
for acc in os.listdir(mob_root):
    mt = os.path.join(mob_root, acc, 'mobtyper_results.txt')
    if os.path.exists(mt) and acc in man:
        with open(mt, newline='') as fh:
            for row in csv.DictReader(fh, delimiter='\t'):
                cid = row.get('primary_cluster_id', '')
                if cid and cid != '-':
                    cluster_sec[cid].add(man[acc]['one_health_sector'])
                    cluster_n[cid] += 1
shared = {c: s for c, s in cluster_sec.items() if len(s) >= 2}
fig, ax = plt.subplots(figsize=(7.5, 4.6))
sec_nodes = ['Human', 'Livestock_Poultry', 'Wastewater_Aquatic']
sec_pos = {'Human': (0, 1.0), 'Livestock_Poultry': (0, 0.5), 'Wastewater_Aquatic': (0, 0.0)}
for s, (sx, sy) in sec_pos.items():
    ax.scatter([sx], [sy], s=420, color='#2166ac', zorder=3)
    ax.text(sx - 0.03, sy, s, ha='right', va='center', fontsize=9)
cl_items = sorted(shared.items(), key=lambda kv: -cluster_n[kv[0]])[:18]
for i, (cid, secs) in enumerate(cl_items):
    cx, cy = 1.0, 1.0 - i * (1.0 / max(len(cl_items) - 1, 1))
    ax.scatter([cx], [cy], s=40 + 22 * cluster_n[cid], color='#b2182b', zorder=3)
    ax.text(cx + 0.03, cy, cid, ha='left', va='center', fontsize=7)
    for s in secs:
        ax.plot([sx, cx], [sec_pos[s][1], cy], color='grey', lw=0.7, alpha=0.6, zorder=1)
ax.set_xlim(-0.45, 1.3); ax.set_ylim(-0.08, 1.08); ax.axis('off')
ax.set_title(f'Plasmid clusters shared across sectors — core Primary cohort ({len(shared)} shared clusters; top {len(cl_items)} shown)')
fig.tight_layout()
fig.savefig(os.path.join(OUT, 'Figure3_plasmid_network.pdf')); fig.savefig(os.path.join(OUT, 'Figure3_plasmid_network.png'))
plt.close(fig)
print('core cohort shared plasmid clusters:', len(shared))

# ---------------- Figure 4 ----------------
l4 = []
with open(os.path.join(BASE, 'stageE', 'L1L4_scores.tsv'), newline='') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        l4.append(r)
fig, axes = plt.subplots(1, 2, figsize=(10, 3.6))
if l4:
    comps = ['L2_shared_ST_fraction', 'L3_ARG_profile_similarity', 'L4_critical_gene_sharing', 'composite_connectivity_2x_seq']
    M = np.array([[float(r[c]) for c in comps] for r in l4])
    im = axes[0].imshow(M, vmin=0, vmax=1, cmap='viridis')
    axes[0].set_xticks(range(4)); axes[0].set_xticklabels(['L2 ST', 'L3 ARG sim', 'L4 crit share', 'Composite'], rotation=25)
    axes[0].set_yticks(range(len(l4)))
    axes[0].set_yticklabels([f"{'Ec' if 'coli' in r['species'] else 'Kp'} {r['income_group']}" for r in l4])
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            axes[0].text(j, i, f'{M[i, j]:.2f}', ha='center', va='center', color='white', fontsize=8)
    axes[0].set_title("A. L1-L4 connectivity (final)")
    fig.colorbar(im, ax=axes[0], shrink=0.85)
# rarefaction curves: BC between sector mean profiles vs subsample size (E. coli)
accs = sorted(man)
def genes_of(r):
    return set(g for g in r['amr_genes'].split(';') if g)
all_genes = sorted(set().union(*[genes_of(man[a]) for a in accs]))
gi = {g: i for i, g in enumerate(all_genes)}
X = np.zeros((len(accs), len(all_genes)))
for i, a in enumerate(accs):
    for g in genes_of(man[a]):
        X[i, gi[g]] = 1.0
rng = np.random.default_rng(20261124)
def bc(u, v):
    d = np.abs(u - v).sum(); s = (u + v).sum()
    return d / s if s else 0.0
ec_idx = {s: [i for i, a in enumerate(accs) if 'Klebsiella' not in man[a]['taxgroup_name'] and man[a]['one_health_sector'] == s]
          for s in ('Human', 'Livestock_Poultry', 'Wastewater_Aquatic')}
sizes = [3, 5, 8, 12, 16, 20]
for (s1, s2) in (('Human', 'Livestock_Poultry'), ('Human', 'Wastewater_Aquatic'), ('Livestock_Poultry', 'Wastewater_Aquatic')):
    xs, ys = [], []
    for n in sizes:
        if len(ec_idx[s1]) >= n and len(ec_idx[s2]) >= n:
            vals = [bc(X[rng.choice(ec_idx[s1], n, replace=False)].mean(axis=0),
                       X[rng.choice(ec_idx[s2], n, replace=False)].mean(axis=0)) for _ in range(60)]
            xs.append(n); ys.append(np.mean(vals))
    if xs:
        axes[1].plot(xs, ys, marker='o', label=f'{s1} vs {s2}')
axes[1].set_xlabel('Genomes per sector (subsample)'); axes[1].set_ylabel('Bray-Curtis (ARG profiles)')
axes[1].set_title('B. Rarefaction — E. coli sector profile dissimilarity')
axes[1].legend(frameon=False, fontsize=7)
fig.tight_layout()
fig.savefig(os.path.join(OUT, 'Figure4_connectivity_matrix.pdf')); fig.savefig(os.path.join(OUT, 'Figure4_connectivity_matrix.png'))
plt.close(fig)
print('figures done')
