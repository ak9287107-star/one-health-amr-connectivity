#!/usr/bin/env python3
"""Stage B (part 1): Sampling frame + core cohort manifest (Protocol S3, S4 Stage 1, S7.1).

Population: global Primary Cohort (asm_level Complete Genome / Chromosome) restricted to
country x species combinations meeting the protocol geographic rule (>=2 sectors
represented in the Primary Cohort; stageA/shortlist_primary.tsv).
Frame: every eligible Primary genome in those cells, stratified by
country x species x sector. Within each stratum, genomes are shuffled with a fixed
seed (20261124) and ranked; rank 1 per stratum forms the CORE analysis cohort
(one genome per country x species x sector cell; 151 cells). This is the protocol's
S7.1 subsampling mitigation applied prospectively: the full Primary Cohort is 10,073
genomes globally (8,987 in shortlist cells) and at the pilot-measured ~6 min/genome
for AMRFinderPlus+MOB-suite could not be processed on available compute (~1,000 h).
The full frame is saved so the sample is reproducible and extensible (rank 2, 3...
can be processed later under the same design).
"""
import csv, json, os, random
from collections import defaultdict

BASE = os.path.expanduser('~/workspace/one_health_amr_pipeline/full_protocol_run')
SEED = 20261124

# income lookup (same as Stage A)
wb = {}
with open(os.path.join(BASE, 'stageA', 'worldbank_countries.json'), encoding='utf-8') as fh:
    data = json.load(fh)
for c in data[1]:
    if c.get('region', {}).get('id') == 'NA':
        continue
    wb[c['name']] = (c['id'], c.get('incomeLevel', {}).get('id', ''))
ALIASES = {'USA': 'United States', 'United States': 'United States',
           'Viet Nam': 'Viet Nam', 'Vietnam': 'Viet Nam', 'Gaza Strip': 'West Bank and Gaza',
           'West Bank': 'West Bank and Gaza', 'Puerto Rico': 'Puerto Rico (US)',
           'South Korea': 'Korea, Rep.', 'Russia': 'Russian Federation', 'Iran': 'Iran, Islamic Rep.',
           'Egypt': 'Egypt, Arab Rep.', 'Hong Kong': 'Hong Kong SAR, China', 'Laos': 'Lao PDR',
           'Syria': 'Syrian Arab Republic', 'Czech Republic': 'Czechia', 'Turkiye': 'Turkiye',
           'Turkey': 'Turkiye', 'Congo': 'Congo, Rep.', 'Democratic Republic of the Congo': 'Congo, Dem. Rep.',
           'Venezuela': 'Venezuela, RB', 'Yemen': 'Yemen, Rep.', 'Kyrgyzstan': 'Kyrgyz Republic',
           'Gambia': 'Gambia, The', 'Bahamas': 'Bahamas, The', 'Cape Verde': 'Cabo Verde',
           'Swaziland': 'Eswatini', 'East Timor': 'Timor-Leste', 'Palestine': 'West Bank and Gaza',
           'Slovakia': 'Slovak Republic', 'Moldova': 'Moldova', 'Brunei': 'Brunei Darussalam'}
def income_for(country):
    name = ALIASES.get(country, country)
    if name in wb:
        iso3, iid = wb[name]
        return iso3, iid, ('HIC' if iid == 'HIC' else ('LMIC' if iid in ('LIC', 'LMC', 'UMC') else 'Unclassified'))
    return '', '', 'Unclassified'

shortlist = set()
with open(os.path.join(BASE, 'stageA', 'shortlist_primary.tsv'), newline='') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        shortlist.add((r['country'], r['species']))

rows = []
with open(os.path.join(BASE, 'stageA', 'eligible_rows.jsonl'), encoding='utf-8') as fh:
    for line in fh:
        r = json.loads(line)
        if r['asm_level'] not in ('Complete Genome', 'Chromosome'):
            continue
        if (r['country'], r['taxgroup_name']) not in shortlist:
            continue
        iso3, iid, grp = income_for(r['country'])
        r['iso3'], r['wb_income_id'], r['income_group'] = iso3, iid, grp
        rows.append(r)
print(f'frame population (primary, shortlist cells): {len(rows)}')

strata = defaultdict(list)
for r in rows:
    strata[(r['country'], r['taxgroup_name'], r['one_health_sector'])].append(r)
rng = random.Random(SEED)
for k in strata:
    strata[k].sort(key=lambda r: r['asm_acc'])
    rng.shuffle(strata[k])
    for i, r in enumerate(strata[k], 1):
        r['stratum'] = '|'.join(k)
        r['within_stratum_rank'] = i
        r['sampled_core'] = (i == 1)

frame_cols = ['asm_acc', 'target_acc', 'biosample_acc', 'taxgroup_name', 'scientific_name', 'country',
              'iso3', 'wb_income_id', 'income_group', 'geo_loc_name', 'collection_year', 'collection_date',
              'one_health_sector', 'epi_type', 'host', 'isolation_source', 'asm_level',
              'asm_stats_contig_n50', 'asm_stats_n_contig', 'asm_stats_length_bp', 'stratum',
              'within_stratum_rank', 'sampled_core', 'amr_genes', 'stress_metal_biocide_genes', 'AST_phenotypes']
os.makedirs(os.path.join(BASE, 'stageB'), exist_ok=True)
with open(os.path.join(BASE, 'stageB', 'sampling_frame.csv'), 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=frame_cols, extrasaction='ignore')
    w.writeheader()
    for r in sorted(rows, key=lambda r: (r['stratum'], r['within_stratum_rank'])):
        w.writerow(r)

core = [r for r in rows if r['sampled_core']]
core.sort(key=lambda r: (r['taxgroup_name'], r['country'], r['one_health_sector']))
with open(os.path.join(BASE, 'stageB', 'cohort_manifest.csv'), 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=frame_cols, extrasaction='ignore')
    w.writeheader()
    w.writerows(core)
print(f'core cohort: {len(core)} genomes from {len(strata)} strata')
from collections import Counter
print('core by species:', Counter(r['taxgroup_name'] for r in core))
print('core by sector:', Counter(r['one_health_sector'] for r in core))
print('core by income:', Counter(r['income_group'] for r in core))
print('core by country (top):', Counter(r['country'] for r in core).most_common(10))
