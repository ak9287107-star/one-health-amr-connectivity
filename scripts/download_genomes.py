#!/usr/bin/env python3
"""Download cohort genomes from GenBank FTP (canonical layout), pilot PATCH 1 route.
Skip logic: existing non-empty .fna files are kept. Pilot genomes are reused by
hardlink/copy if the accession was already downloaded for the pilot.
Usage: download_genomes.py <manifest_csv> <genome_dir> <failures_txt>
"""
import csv, os, re, subprocess, sys

manifest, gdir, failpath = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(gdir, exist_ok=True)
PILOT = os.path.expanduser('~/workspace/one_health_amr_pipeline/pilot_run/oh_amr_results/genomes')

accs = []
with open(manifest, newline='') as fh:
    for r in csv.DictReader(fh):
        accs.append(r['asm_acc'])
print(f'{len(accs)} accessions', flush=True)

fails = []
done = 0
for acc in accs:
    out = os.path.join(gdir, acc + '.fna')
    if os.path.exists(out) and os.path.getsize(out) > 0:
        done += 1
        continue
    pilot_f = os.path.join(PILOT, acc + '.fna')
    if os.path.exists(pilot_f) and os.path.getsize(pilot_f) > 0:
        subprocess.run(['cp', pilot_f, out], check=True)
        done += 1
        print(f'{acc} reused from pilot', flush=True)
        continue
    digits = re.sub(r'^GCA_([0-9]{9})\..*', r'\1', acc)
    if digits == acc:
        # GCF or other prefix: try GCA layout anyway with digits from acc
        m = re.search(r'([0-9]{9})', acc)
        if not m:
            fails.append((acc, 'UNPARSEABLE_ACC'))
            continue
        digits = m.group(1)
    base = f'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/{digits[0:3]}/{digits[3:6]}/{digits[6:9]}/'
    try:
        listing = subprocess.run(['curl', '-s', '--max-time', '60', base],
                                 capture_output=True, text=True).stdout
    except Exception as e:
        fails.append((acc, f'LISTING_ERROR {e}'))
        continue
    m = re.search(re.escape(acc) + r'_[A-Za-z0-9._-]+', listing)
    if not m:
        fails.append((acc, 'NO_FOLDER'))
        continue
    folder = m.group(0)
    url = f'{base}{folder}/{folder}_genomic.fna.gz'
    gz = out + '.gz'
    rc = subprocess.run(['curl', '-s', '--fail', '--max-time', '600', '-o', gz, url]).returncode
    if rc != 0:
        fails.append((acc, 'DOWNLOAD_FAILED'))
        if os.path.exists(gz):
            os.remove(gz)
        continue
    subprocess.run(['gunzip', '-f', gz], check=True)
    done += 1
    print(f'{acc} downloaded ({done})', flush=True)

with open(failpath, 'w') as fh:
    for acc, why in fails:
        fh.write(f'{acc}\t{why}\n')
print(f'DONE downloaded/staged={done} failures={len(fails)}')
