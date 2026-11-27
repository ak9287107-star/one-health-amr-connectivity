#!/usr/bin/env python3
"""Stage A: Global coverage audit (Protocol S3).
Source: NCBI Pathogen Detection FTP AMR metadata (local copies downloaded 2026-11-23
from the same URLs the pilot used; see PILOT_REPORT.md S1).
Filters mirror Query 2 / build_manifest.py exactly (patched sector logic, no bare 'water'):
  asm_acc + geo_loc_name present, amrfinder_version present, collection_year 2015-2026,
  sector in {Human, Livestock_Poultry, Wastewater_Aquatic}; Klebsiella filtered to
  scientific_name='Klebsiella pneumoniae'.
Cohorts: Primary = asm_level in {Complete Genome, Chromosome};
         Sensitivity = N50>=50000 AND n_contig<=200 (any level);
         All = all eligible levels.
Distinct asm_acc per (country, species, sector); an asm_acc in >1 sector/country keeps
its first occurrence for cohort-row purposes but is counted once per cell it appears in
(cells are the audit unit, matching Query 1 semantics).
Income: World Bank API country list (incomeLevel), fetched 2026-11-24.
"""
import csv, re, sys, os, json
from collections import Counter, defaultdict

BASE = os.path.expanduser('~/workspace/one_health_amr_pipeline')
SRC = os.path.join(BASE, 'pilot_run', 'source_data')
OUT = os.path.join(BASE, 'full_protocol_run', 'stageA')
os.makedirs(OUT, exist_ok=True)

WASTE = re.compile(r'wastewater|sewage|sludge|effluent|influent|river|surface water|pond|lake|canal|stream|aquatic|sediment|drain')
LIVE = re.compile(r'gallus|chicken|broiler|poultry|turkey|duck|sus scrofa|pig|swine|porcine|bos taurus|cattle|cow|bovine|calf|sheep|goat|livestock|meat|caec|ceca|abattoir|slaughter')
HUMAN_HOST = re.compile(r'homo sapiens|human')
YEAR = re.compile(r'^(20[12][0-9])')

def sector(host, source, epi):
    if WASTE.search((source or '').lower()): return 'Wastewater_Aquatic'
    if LIVE.search(((host or '') + ' ' + (source or '')).lower()): return 'Livestock_Poultry'
    if HUMAN_HOST.search((host or '').lower()) or (epi or '').lower() == 'clinical': return 'Human'
    return 'Other_Unclassified'

def clean(v):
    if v is None: return ''
    v = v.strip()
    return '' if v == 'NULL' else v

def parse_genotypes(s):
    s = clean(s).strip('"')
    if not s: return []
    return [g.strip().strip('"') for g in s.split(',') if g.strip()]

# ---------- World Bank income ----------
wb = {}
with open(os.path.join(OUT, 'worldbank_countries.json'), encoding='utf-8') as fh:
    data = json.load(fh)
for c in data[1]:
    if c.get('region', {}).get('id') == 'NA':  # aggregates
        continue
    wb[c['name']] = (c['id'], c.get('incomeLevel', {}).get('id', ''), c.get('incomeLevel', {}).get('value', ''))
ALIASES = {
    'United States': 'United States', 'USA': 'United States',
    'Viet Nam': 'Viet Nam', 'Vietnam': 'Viet Nam',
    'Gaza Strip': 'West Bank and Gaza', 'West Bank': 'West Bank and Gaza',
    'Puerto Rico': 'Puerto Rico (US)',
    'South Korea': 'Korea, Rep.', 'Korea, South': 'Korea, Rep.', 'Republic of Korea': 'Korea, Rep.',
    'North Korea': "Korea, Dem. People's Rep.",
    'Russia': 'Russian Federation', 'Russian Federation': 'Russian Federation',
    'Iran': 'Iran, Islamic Rep.', 'Egypt': 'Egypt, Arab Rep.', 'Yemen': 'Yemen, Rep.',
    'Venezuela': 'Venezuela, RB', 'Bolivia': 'Bolivia', 'Tanzania': 'Tanzania',
    'United Kingdom': 'United Kingdom', 'UK': 'United Kingdom',
    'Hong Kong': 'Hong Kong SAR, China', 'Laos': 'Lao PDR', "Lao People's Democratic Republic": 'Lao PDR',
    'Syria': 'Syrian Arab Republic', 'Moldova': 'Moldova', 'Brunei': 'Brunei Darussalam',
    'Congo': 'Congo, Rep.', 'Democratic Republic of the Congo': 'Congo, Dem. Rep.',
    'Ivory Coast': "Cote d'Ivoire", "Côte d'Ivoire": "Cote d'Ivoire",
    'Czech Republic': 'Czechia', 'Czechia': 'Czechia', 'Slovakia': 'Slovak Republic',
    'Kyrgyzstan': 'Kyrgyz Republic', 'Gambia': 'Gambia, The', 'Bahamas': 'Bahamas, The',
    'Micronesia': 'Micronesia, Fed. Sts.', 'Cape Verde': 'Cabo Verde',
    'Swaziland': 'Eswatini', 'East Timor': 'Timor-Leste', 'Palestine': 'West Bank and Gaza',
    'Macau': 'Macao SAR, China', 'Saint Kitts and Nevis': 'St. Kitts and Nevis',
    'Saint Lucia': 'St. Lucia', 'Saint Vincent and the Grenadines': 'St. Vincent and the Grenadines',
    'Sao Tome and Principe': 'Sao Tome and Principe', 'Turkiye': 'Turkiye', 'Turkey': 'Turkiye',
    'Netherlands': 'Netherlands', 'France': 'France', 'Germany': 'Germany',
}
def income_for(country):
    name = ALIASES.get(country, country)
    if name in wb:
        iso3, iid, ival = wb[name]
        group = 'HIC' if iid == 'HIC' else ('LMIC' if iid in ('LIC', 'LMC', 'UMC') else 'Unclassified')
        return iso3, iid, ival, group
    return '', '', '', 'Unclassified'

# ---------- stream metadata ----------
# cell[(country,species,sector)][cohort] = set(asm_acc)
cells = defaultdict(lambda: {'primary': set(), 'sensitivity': set(), 'all': set()})
row_store = {}  # asm_acc -> full row (first occurrence), for Stage B reuse

def load_group(path, label, species_filter=None):
    n = 0
    with open(path, newline='', encoding='utf-8', errors='replace') as fh:
        r = csv.DictReader(fh, delimiter='\t')
        for row in r:
            if species_filter and clean(row.get('scientific_name')) != species_filter:
                continue
            asm = clean(row.get('asm_acc'))
            geo = clean(row.get('geo_loc_name'))
            if not asm or not geo: continue
            if not clean(row.get('amrfinder_version')): continue
            cdate = clean(row.get('collection_date'))
            ym = YEAR.match(cdate)
            if not ym: continue
            year = int(ym.group(1))
            if not (2015 <= year <= 2026): continue
            sec = sector(clean(row.get('host')), clean(row.get('isolation_source')), clean(row.get('epi_type')))
            if sec not in ('Human', 'Livestock_Poultry', 'Wastewater_Aquatic'): continue
            country = geo.split(':')[0].strip()
            lvl = clean(row.get('asm_level'))
            try: n50 = int(float(clean(row.get('asm_stats_contig_n50')) or 0))
            except ValueError: n50 = 0
            try: nctg = int(float(clean(row.get('asm_stats_n_contig')) or 10**9))
            except ValueError: nctg = 10**9
            cell = cells[(country, label, sec)]
            cell['all'].add(asm)
            if lvl in ('Complete Genome', 'Chromosome'): cell['primary'].add(asm)
            if n50 >= 50000 and nctg <= 200: cell['sensitivity'].add(asm)
            if asm not in row_store:
                amr = parse_genotypes(row.get('AMR_genotypes'))
                stress = parse_genotypes(row.get('stress_genotypes'))
                row_store[asm] = {
                    'asm_acc': asm, 'target_acc': clean(row.get('target_acc')),
                    'biosample_acc': clean(row.get('biosample_acc')),
                    'taxgroup_name': label, 'country': country, 'geo_loc_name': geo,
                    'collection_year': year, 'collection_date': cdate,
                    'one_health_sector': sec, 'epi_type': clean(row.get('epi_type')),
                    'host': clean(row.get('host')), 'isolation_source': clean(row.get('isolation_source')),
                    'asm_level': lvl, 'asm_stats_contig_n50': n50, 'asm_stats_n_contig': nctg,
                    'asm_stats_length_bp': clean(row.get('asm_stats_length_bp')),
                    'scientific_name': clean(row.get('scientific_name')),
                    'amr_genes': ';'.join(sorted(set(amr))),
                    'stress_metal_biocide_genes': ';'.join(sorted(set(stress))),
                    'AST_phenotypes': clean(row.get('AST_phenotypes')),
                }
            n += 1
    print(f'{label}: {n} eligible rows', file=sys.stderr)

load_group(os.path.join(SRC, 'ecoli.amr.metadata.tsv'), 'E.coli and Shigella')
load_group(os.path.join(SRC, 'klebsiella.amr.metadata.tsv'), 'Klebsiella pneumoniae', species_filter='Klebsiella pneumoniae')

# ---------- write audit ----------
audit_path = os.path.join(OUT, 'global_coverage_audit.tsv')
with open(audit_path, 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['country', 'iso3', 'wb_income_id', 'wb_income_value', 'income_group', 'species', 'sector',
                'n_primary', 'n_sensitivity', 'n_all_levels'])
    for (country, species, sec) in sorted(cells):
        iso3, iid, ival, grp = income_for(country)
        c = cells[(country, species, sec)]
        w.writerow([country, iso3, iid, ival, grp, species, sec,
                    len(c['primary']), len(c['sensitivity']), len(c['all'])])

# shortlists: countries with >=2 sectors represented (n>0) per cohort
def shortlist(cohort):
    by_country = defaultdict(set)
    celln = {}
    for (country, species, sec), c in cells.items():
        if len(c[cohort]) > 0:
            by_country[(country, species)].add(sec)
            celln[(country, species, sec)] = len(c[cohort])
    out = []
    for (country, species), secs in sorted(by_country.items()):
        if len(secs) >= 2:
            iso3, iid, ival, grp = income_for(country)
            row = {'country': country, 'income_group': grp, 'species': species,
                   'sectors_represented': len(secs)}
            tot = 0
            for s in ('Human', 'Livestock_Poultry', 'Wastewater_Aquatic'):
                row['n_' + s] = celln.get((country, species, s), 0)
                tot += row['n_' + s]
            row['n_total'] = tot
            out.append(row)
    return out

for cohort in ('primary', 'sensitivity'):
    sl = shortlist(cohort)
    with open(os.path.join(OUT, f'shortlist_{cohort}.tsv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, delimiter='\t', fieldnames=['country', 'income_group', 'species', 'sectors_represented', 'n_Human', 'n_Livestock_Poultry', 'n_Wastewater_Aquatic', 'n_total'])
        w.writeheader(); w.writerows(sl)
    print(f'shortlist {cohort}: {len(sl)} country x species entries', file=sys.stderr)

# headline totals per species x sector x income x cohort (distinct asm_acc globally)
tot = defaultdict(set)
for (country, species, sec), c in cells.items():
    _, _, _, grp = income_for(country)
    for cohort in ('primary', 'sensitivity', 'all'):
        tot[(species, sec, grp, cohort)] |= c[cohort]
        tot[(species, sec, 'ALL_INCOME', cohort)] |= c[cohort]
        tot[('ALL_SPECIES', sec, grp, cohort)] |= c[cohort]
        tot[('ALL_SPECIES', sec, 'ALL_INCOME', cohort)] |= c[cohort]
with open(os.path.join(OUT, 'audit_headline_totals.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['species', 'sector', 'income_group', 'cohort', 'n_distinct_asm_acc'])
    for k in sorted(tot): w.writerow([*k, len(tot[k])])

unmatched = sorted({ct for (ct, _, _) in cells if income_for(ct)[3] == 'Unclassified'})
print('countries total:', len({ct for (ct, _, _) in cells}), '| income-unmatched:', unmatched, file=sys.stderr)
print('stored rows (distinct asm_acc):', len(row_store), file=sys.stderr)

# save row store for Stage B
with open(os.path.join(OUT, 'eligible_rows.jsonl'), 'w', encoding='utf-8') as fh:
    for asm, r in row_store.items():
        fh.write(json.dumps(r) + '\n')
print('Stage A done.')
