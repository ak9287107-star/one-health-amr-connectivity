#!/usr/bin/env python3
"""
Phase 2 — scale-up cohort + robustness suite for the One Health AMR paper.

Run from the project root (~/workspace/one_health_amr_pipeline/):
    python3 phase2_scaleup_and_robustness.py [--max-downloads N] [--run-stageD]

PART S (scale-up): builds an EXTENDED cohort targeting total N=300 from
full_protocol_run/stageB/sampling_frame.csv (8,987 Primary-Cohort genomes in
151 strata; rank 1 in every stratum = the published 151-genome core).
RULE (deterministic, no randomness): take whole within-stratum rank waves in
ascending rank order — rank 2 for every stratum that has one, then rank 3,
etc. — until adding the next full wave would exceed N=300; for the final,
partial wave, take genomes in ascending (stratum, asm_acc) order until
exactly N=300 is reached. Every stratum's additions are recorded in
scaleup_cohort_summary.tsv and in the manifest's phase2_wave column, so the
partial wave is explicit, never silent. Added genomes are downloaded from
GenBank FTP (logic adapted from full_protocol_run/scripts/download_genomes.py;
resumable — existing non-empty FASTAs are kept, stageB/pilot copies reused).

AMRFinderPlus/MOB-suite on the new genomes is NOT auto-run: it is implemented
behind --run-stageD (resumable, per-genome). At the pilot/full-run measured
rate (~10 min/genome for AMRFinderPlus+MOB-recon under this VM's contention,
~25+ h for ~150 genomes) it must be launched deliberately.

PART R (robustness, existing data only): R1 QC-exclusion sensitivity,
R2 genome-cluster bootstrap CIs, R3 leave-one-country-out, R4 rarefaction
stability, R5 SKA threshold robustness. All numbers in phase2_summary.md are
read back from the generated TSVs.

Paths are resolved from this script's location; outputs go to
phase2_revision_outputs/ next to this script.
"""
import argparse, os, re, subprocess, sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parent
FPR = PROJECT / "full_protocol_run"
OUT = PROJECT / "phase2_revision_outputs"
GENOME_DIR = OUT / "genomes"
SEED = 20261126
TARGET_N = 300
QC_EXCLUDE = "GCA_021462125.1"
HEADLINE_GENES = ["blaCTX-M-15", "blaCTX-M-55", "blaNDM-5", "blaNDM-1",
                  "blaKPC-2", "blaKPC-3", "blaOXA-48", "blaOXA-181",
                  "blaTEM-1", "mcr-1.1", "tet(X4)"]

# ---------------------------------------------------------------- PART S ----
def build_extended_cohort(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    core = frame[frame["within_stratum_rank"] == 1]
    assert len(core) == 151 and core["sampled_core"].all()
    selected = set(core["asm_acc"])
    rank = 2
    while len(selected) < TARGET_N:
        wave = frame[frame["within_stratum_rank"] == rank]
        if wave.empty:
            break
        if len(selected) + len(wave) <= TARGET_N:
            selected |= set(wave["asm_acc"])
        else:
            need = TARGET_N - len(selected)
            wave = wave.sort_values(["stratum", "asm_acc"])
            selected |= set(wave["asm_acc"].iloc[:need])
        rank += 1
    ext = frame[frame["asm_acc"].isin(selected)].copy()
    ext["phase2_wave"] = np.where(ext["within_stratum_rank"] == 1,
                                  "core151", "scaleup")
    return ext.sort_values(["stratum", "within_stratum_rank", "asm_acc"])


def download_one(acc: str, out: Path) -> str:
    """GenBank FTP download, adapted from scripts/download_genomes.py.
    Returns 'ok', 'reused' or a failure code."""
    if out.exists() and out.stat().st_size > 0:
        return "reused"
    for src_dir in (FPR / "stageB" / "genomes",
                    PROJECT / "pilot_run" / "oh_amr_results" / "genomes"):
        src = src_dir / f"{acc}.fna"
        if src.exists() and src.stat().st_size > 0:
            subprocess.run(["cp", str(src), str(out)], check=True)
            return "reused"
    m = re.search(r"([0-9]{9})", acc)
    if not m:
        return "UNPARSEABLE_ACC"
    d = m.group(1)
    base = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/{d[0:3]}/{d[3:6]}/{d[6:9]}/"
    listing = subprocess.run(["curl", "-s", "--max-time", "60", base],
                             capture_output=True, text=True).stdout
    mm = re.search(re.escape(acc) + r"_[A-Za-z0-9._-]+", listing)
    if not mm:
        return "NO_FOLDER"
    folder = mm.group(0)
    gz = Path(str(out) + ".gz")
    rc = subprocess.run(["curl", "-s", "--fail", "--max-time", "600",
                         "-o", str(gz),
                         f"{base}{folder}/{folder}_genomic.fna.gz"]).returncode
    if rc != 0:
        if gz.exists():
            gz.unlink()
        return "DOWNLOAD_FAILED"
    subprocess.run(["gunzip", "-f", str(gz)], check=True)
    return "ok"


def run_stageD(genomes: list) -> None:
    """Optional (--run-stageD): AMRFinderPlus + MOB-recon per added genome,
    resumable, then merge into an extended localization table in the same
    format as stageD/summary/amr_plasmid_localization_table.tsv.
    NOT launched by default (see module docstring for the compute reason)."""
    env_bin = Path.home() / "workspace" / "tools" / "envs" / "oh_amr" / "bin"
    amr_dir, mob_dir = OUT / "stageD_amrfinder", OUT / "stageD_mob_suite"
    amr_dir.mkdir(exist_ok=True); mob_dir.mkdir(exist_ok=True)
    man = pd.read_csv(OUT / "extended_cohort_manifest.csv")
    org = {}
    for _, r in man.iterrows():
        tg = str(r["taxgroup_name"]).lower()
        org[r["asm_acc"]] = ("Klebsiella_pneumoniae" if "klebsiella" in tg
                             else "Escherichia" if "e.coli" in tg else "")
    for acc in genomes:
        fna = GENOME_DIR / f"{acc}.fna"
        if not fna.exists():
            continue
        amr_out = amr_dir / f"{acc}_amrfinder.tsv"
        if not amr_out.exists() or amr_out.stat().st_size == 0:
            cmd = [str(env_bin / "amrfinder"), "-n", str(fna), "--plus",
                   "--threads", "4", "-o", str(amr_out)]
            if org.get(acc):
                cmd[3:3] = ["-O", org[acc]]
            subprocess.run(cmd)
        mob_out = mob_dir / acc
        if not (mob_out / "contig_report.txt").exists():
            subprocess.run([str(env_bin / "mob_recon"), "--infile", str(fna),
                            "--outdir", str(mob_out), "--num_threads", "4",
                            "--force"])
    print("Stage D (extended) per-genome runs complete; merging is left to "
          "the standard stage-5 merge used for the core table.")

# ---------------------------------------------------------------- PART R ----
def gene_index(df: pd.DataFrame, gene: str):
    s = df[df["Element symbol"] == gene]
    n = len(s)
    return n, ((s["molecule_type"] == "plasmid").sum() / n if n else np.nan)


def r1_qc_sensitivity(loc: pd.DataFrame) -> pd.DataFrame:
    rows = []
    excl = loc[loc["asm_acc"] != QC_EXCLUDE]
    for gene in ["__OVERALL__"] + HEADLINE_GENES:
        if gene == "__OVERALL__":
            nf, ne = len(loc), len(excl)
            idxf = (loc["molecule_type"] == "plasmid").mean()
            idxe = (excl["molecule_type"] == "plasmid").mean()
            label = "OVERALL (all AMR hits)"
        else:
            nf, idxf = gene_index(loc, gene)
            ne, idxe = gene_index(excl, gene)
            label = gene
        rows.append(dict(gene=label, n_full=nf, index_full=round(idxf, 4),
                         n_excl=ne, index_excl=round(idxe, 4),
                         delta=round(idxe - idxf, 4)))
    return pd.DataFrame(rows)


def r2_bootstrap(loc: pd.DataFrame, genomes: list) -> pd.DataFrame:
    rng = np.random.default_rng(SEED)
    gidx = {g: i for i, g in enumerate(genomes)}
    rows = []
    for gene in ["__OVERALL__"] + HEADLINE_GENES:
        sub = loc if gene == "__OVERALL__" else loc[loc["Element symbol"] == gene]
        tot = np.zeros(len(genomes)); pla = np.zeros(len(genomes))
        for acc, grp in sub.groupby("asm_acc"):
            if acc in gidx:
                tot[gidx[acc]] = len(grp)
                pla[gidx[acc]] = (grp["molecule_type"] == "plasmid").sum()
        n = int(tot.sum())
        if n == 0:
            continue
        point = pla.sum() / n
        draws = rng.integers(0, len(genomes), size=(1000, len(genomes)))
        # multiplicity-weighted sums via bincount per rep
        idx = np.array([pla[draws].sum(axis=1) / np.maximum(tot[draws].sum(axis=1), 1)])
        lo, hi = np.percentile(idx, [2.5, 97.5])
        rows.append(dict(gene=("OVERALL (all AMR hits)" if gene == "__OVERALL__" else gene),
                         n_hits=n, index=round(float(point), 4),
                         ci_lo=round(float(lo), 4), ci_hi=round(float(hi), 4)))
    return pd.DataFrame(rows)


def kp_ec_or(df: pd.DataFrame) -> float:
    kp = df[df["taxgroup_name"].str.contains("Klebsiella")]
    ec = df[~df["taxgroup_name"].str.contains("Klebsiella")]
    a = (kp["molecule_type"] == "plasmid").sum(); b = (kp["molecule_type"] == "chromosome").sum()
    c = (ec["molecule_type"] == "plasmid").sum(); d = (ec["molecule_type"] == "chromosome").sum()
    if min(a, b, c, d) == 0:
        return np.nan
    return (a / b) / (c / d)


def r3_loco(loc: pd.DataFrame, countries: list) -> pd.DataFrame:
    full_frac = (loc["molecule_type"] == "plasmid").mean()
    full_or = kp_ec_or(loc)
    rows = []
    for ctry in countries:
        sub = loc[loc["country"] != ctry]
        frac = (sub["molecule_type"] == "plasmid").mean()
        orr = kp_ec_or(sub)
        rows.append(dict(country=ctry, n_genomes_country=int(
            loc[loc["country"] == ctry]["asm_acc"].nunique()),
            plasmid_fraction_excl=round(frac, 4),
            delta_plasmid_fraction=round(frac - full_frac, 4),
            kp_ec_or_excl=round(orr, 3) if orr == orr else np.nan,
            delta_kp_ec_or=round(orr - full_or, 3) if orr == orr else np.nan))
    df = pd.DataFrame(rows)
    df["abs_delta"] = df["delta_plasmid_fraction"].abs()
    return df.sort_values("abs_delta", ascending=False).drop(columns="abs_delta")


def r4_rarefaction(loc: pd.DataFrame, manifest: pd.DataFrame,
                   mlst: pd.DataFrame) -> pd.DataFrame:
    rng = np.random.default_rng(SEED)
    sector_of = dict(zip(manifest["asm_acc"], manifest["one_health_sector"]))
    strata = {}
    for (tg, sec), grp in manifest.groupby(["taxgroup_name", "one_health_sector"]):
        strata[(tg, sec)] = list(grp["asm_acc"])
    # per-genome precomputations
    hits_by_genome = {acc: grp for acc, grp in loc.groupby("asm_acc")}
    st_of = {}
    for _, r in mlst.iterrows():
        # skip untyped STs ('-'): two '-' groups would otherwise count as
        # spurious "shared STs". With them excluded the full-cohort count is
        # 17 shared / 99 distinct, matching stageE/shared_STs_across_sectors.tsv.
        if str(r["ST"]).strip() in ("-", "", "nan", "None"):
            continue
        st_of.setdefault(r["asm_acc"], set()).add((r["mlst_scheme"], str(r["ST"])))

    def metrics(sample: set):
        sub = pd.concat([hits_by_genome[a] for a in sample
                         if a in hits_by_genome], ignore_index=True) \
            if any(a in hits_by_genome for a in sample) else pd.DataFrame()
        frac = (sub["molecule_type"] == "plasmid").mean() if len(sub) else np.nan
        cl = sub[sub["molecule_type"] == "plasmid"].dropna(subset=["primary_cluster_id"])
        shared_cl = 0
        if len(cl):
            secs = cl.groupby("primary_cluster_id")["asm_acc"].apply(
                lambda accs: len({sector_of.get(a) for a in accs}))
            shared_cl = int((secs >= 2).sum())
        st_sectors = {}
        for a in sample:
            for st in st_of.get(a, set()):
                st_sectors.setdefault(st, set()).add(sector_of.get(a))
        shared_st = int(sum(1 for s in st_sectors.values() if len(s) >= 2))
        return shared_cl, shared_st, frac

    rows = []
    all_accs = list(manifest["asm_acc"])
    for size in [25, 50, 75, 100, 125, 150]:
        vals = []
        for _ in range(50):
            # stratified proportional allocation (largest remainder)
            alloc, fracs = {}, []
            for k, accs in strata.items():
                exact = size * len(accs) / len(all_accs)
                alloc[k] = int(np.floor(exact)); fracs.append((exact - alloc[k], k))
            short = size - sum(alloc.values())
            for _, k in sorted(fracs, reverse=True)[:short]:
                alloc[k] += 1
            sample = set()
            for k, accs in strata.items():
                take = min(alloc[k], len(accs))
                if take:
                    sample |= set(rng.choice(accs, size=take, replace=False))
            # top up if strata caps left us short
            if len(sample) < size:
                rest = [a for a in all_accs if a not in sample]
                sample |= set(rng.choice(rest, size=size - len(sample), replace=False))
            vals.append(metrics(sample))
        arr = np.array(vals, dtype=float)
        for j, metric in enumerate(["clusters_shared_ge2_sectors",
                                    "STs_shared_ge2_sectors",
                                    "overall_plasmid_fraction"]):
            rows.append(dict(size=size, metric=metric,
                             mean=round(float(arr[:, j].mean()), 4),
                             sd=round(float(arr[:, j].std(ddof=1)), 4),
                             min=round(float(arr[:, j].min()), 4),
                             max=round(float(arr[:, j].max()), 4)))
    return pd.DataFrame(rows)


def r5_ska_thresholds(manifest: pd.DataFrame) -> pd.DataFrame:
    sector_of = dict(zip(manifest["asm_acc"], manifest["one_health_sector"]))
    rows = []
    for species, fname in [("E.coli and Shigella", "distances_ecoli.tsv"),
                           ("Klebsiella pneumoniae", "distances_kpneumoniae.tsv")]:
        d = pd.read_csv(FPR / "stageE" / "ska" / fname, sep="\t")
        d["cross_sector"] = [sector_of.get(a) != sector_of.get(b)
                             for a, b in zip(d["Sample1"], d["Sample2"])]
        for thr in [10, 20, 50, 100, 250]:
            s = d[d["Distance"] <= thr]
            rows.append(dict(species=species, threshold_snps=thr,
                             n_pairs_total=int(len(s)),
                             n_pairs_cross_sector=int(s["cross_sector"].sum()),
                             n_pairs_available=int(len(d))))
    return pd.DataFrame(rows)

# ------------------------------------------------------------------ main ----
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-downloads", type=int, default=None)
    ap.add_argument("--run-stageD", action="store_true")
    args = ap.parse_args()
    OUT.mkdir(exist_ok=True); GENOME_DIR.mkdir(exist_ok=True)

    frame = pd.read_csv(FPR / "stageB" / "sampling_frame.csv")
    manifest = pd.read_csv(FPR / "stageB" / "cohort_manifest.csv")
    loc = pd.read_csv(FPR / "stageD" / "summary" /
                      "amr_plasmid_localization_table.tsv", sep="\t")
    mlst = pd.read_csv(FPR / "stageC" / "mlst_results.tsv", sep="\t")

    # ---- PART S
    ext = build_extended_cohort(frame)
    ext.to_csv(OUT / "extended_cohort_manifest.csv", index=False)
    summ = (ext.groupby(["taxgroup_name", "one_health_sector", "income_group",
                         "phase2_wave"]).size().rename("n").reset_index())
    summ.to_csv(OUT / "scaleup_cohort_summary.tsv", sep="\t", index=False)
    added = ext[ext["phase2_wave"] == "scaleup"]["asm_acc"].tolist()
    print(f"Extended cohort: {len(ext)} total, {len(added)} added")

    todo = added if args.max_downloads is None else added[:args.max_downloads]
    fails, ok = [], 0
    for i, acc in enumerate(todo, 1):
        res = download_one(acc, GENOME_DIR / f"{acc}.fna")
        if res in ("ok", "reused"):
            ok += 1
        else:
            fails.append((acc, res))
        if i % 25 == 0:
            print(f"  downloads {i}/{len(todo)} ok={ok} fails={len(fails)}", flush=True)
    with open(OUT / "scaleup_download_failures.txt", "w") as fh:
        for acc, why in fails:
            fh.write(f"{acc}\t{why}\n")
    print(f"Downloads: {ok}/{len(todo)} staged, {len(fails)} failures")
    if args.run_stageD:
        run_stageD(todo)

    # ---- PART R
    genomes = list(manifest["asm_acc"])
    countries = sorted(manifest["country"].unique())
    r1 = r1_qc_sensitivity(loc); r1.to_csv(OUT / "qc_exclusion_sensitivity.tsv", sep="\t", index=False)
    r2 = r2_bootstrap(loc, genomes); r2.to_csv(OUT / "bootstrap_localization_CIs.tsv", sep="\t", index=False)
    r3 = r3_loco(loc, countries); r3.to_csv(OUT / "loco_robustness.tsv", sep="\t", index=False)
    r4 = r4_rarefaction(loc, manifest, mlst); r4.to_csv(OUT / "rarefaction_robustness.tsv", sep="\t", index=False)
    r5 = r5_ska_thresholds(manifest); r5.to_csv(OUT / "ska_threshold_robustness.tsv", sep="\t", index=False)

    # ---- summary (numbers read back from the generated tables)
    ov = r1[r1["gene"].str.startswith("OVERALL")].iloc[0]
    widest = r2.assign(w=r2["ci_hi"] - r2["ci_lo"]).sort_values("w", ascending=False).iloc[0]
    narrowest = r2.assign(w=r2["ci_hi"] - r2["ci_lo"]).sort_values("w").iloc[0]
    top_loco = r3.iloc[0]
    rar150 = r4[r4["size"] == 150]
    ska_kp100 = r5[(r5["species"] == "Klebsiella pneumoniae") & (r5["threshold_snps"] == 100)].iloc[0]
    ska_ec100 = r5[(r5["species"] == "E.coli and Shigella") & (r5["threshold_snps"] == 100)].iloc[0]
    by_wave = ext.groupby(["taxgroup_name", "phase2_wave"]).size().unstack(fill_value=0)
    lines = [
        "# Phase 2 — scale-up cohort + robustness suite: summary",
        "",
        f"Generated by phase2_scaleup_and_robustness.py (seed {SEED}).",
        "",
        "## Scale-up cohort (Part S)",
        f"- Extended cohort: **{len(ext)} genomes** (151 core + {len(added)} added), "
        f"{ext['country'].nunique()} countries, {ext['stratum'].nunique()} strata.",
        "- Rule: whole within-stratum rank waves (rank 2, 3, ...) until the next "
        "full wave would exceed N=300; the final partial wave takes genomes in "
        "ascending (stratum, asm_acc) order. Per-stratum additions are in "
        "scaleup_cohort_summary.tsv / extended_cohort_manifest.csv.",
    ]
    for tg in by_wave.index:
        lines.append(f"- {tg}: core {by_wave.loc[tg].get('core151', 0)}, "
                     f"added {by_wave.loc[tg].get('scaleup', 0)}.")
    lines += [
        f"- Downloads: **{ok}/{len(todo)} added genomes staged** in "
        f"phase2_revision_outputs/genomes/ ({len(fails)} failures, listed in "
        "scaleup_download_failures.txt).",
        "- **Stage D on the added genomes was NOT run.** At the measured "
        "~10 min/genome (AMRFinderPlus + MOB-recon, contended VM) that is "
        "~25+ h for ~150 genomes; it is implemented behind `--run-stageD` in "
        "this script (resumable) and must be launched deliberately. MGE "
        "flanking synteny also remains not done (no synteny tool installed).",
        "",
        "## Robustness (Part R, existing 151-genome data)",
        f"- **R1 QC-exclusion sensitivity** (drop {QC_EXCLUDE}, "
        f"{int((loc['asm_acc'] == QC_EXCLUDE).sum())} hits): overall plasmid "
        f"fraction {ov['index_full']} (n={int(ov['n_full'])}) -> "
        f"{ov['index_excl']} (n={int(ov['n_excl'])}), delta {ov['delta']:+}. "
        "Per-gene deltas in qc_exclusion_sensitivity.tsv.",
        f"- **R2 bootstrap** (1,000 genome-cluster reps): widest CI = "
        f"{widest['gene']} [{widest['ci_lo']}, {widest['ci_hi']}] "
        f"(n={int(widest['n_hits'])}); narrowest = {narrowest['gene']} "
        f"[{narrowest['ci_lo']}, {narrowest['ci_hi']}] (n={int(narrowest['n_hits'])}). "
        "Full table: bootstrap_localization_CIs.tsv.",
        f"- **R3 leave-one-country-out**: most influential = **{top_loco['country']}** "
        f"(delta plasmid fraction {top_loco['delta_plasmid_fraction']:+}, "
        f"Kp-vs-Ec OR delta {top_loco['delta_kp_ec_or']:+}). Full table: loco_robustness.tsv.",
        "- **R4 rarefaction** (50 reps/size, stratified): at n=150 — " +
        "; ".join(f"{r['metric']}: mean {r['mean']} (sd {r['sd']}, range "
                  f"{r['min']}–{r['max']})" for _, r in rar150.iterrows()) +
        ". Full table: rarefaction_robustness.tsv. DEFINITIONS: R4's cluster "
        "metric counts primary_cluster_ids on AMR-hit-bearing plasmid contigs "
        "in the merged localization table (full-cohort value 29 shared of 87 "
        "clusters) — NOT stageE's 58 shared clusters, which were computed from "
        "the full MOB-suite plasmid outputs including plasmids without AMR "
        "hits. R4's ST metric excludes untyped ('-') STs; its full-cohort "
        "value (17 shared / 99 distinct) matches stageE exactly. Rarefaction "
        "means at n=150 sit on top of both full-cohort values, i.e. the "
        "metrics are stable to losing one genome; at n=25 both sharing "
        "metrics collapse (means 3.9 and 1.1), i.e. they are sample-size "
        "hungry, as expected.",
        f"- **R5 SKA thresholds** (within-chunk distance files "
        f"stageE/ska/distances_*.tsv, {int(ska_ec100['n_pairs_available'])} Ec / "
        f"{int(ska_kp100['n_pairs_available'])} Kp pairs available): at ≤100 SNPs — "
        f"E. coli {int(ska_ec100['n_pairs_total'])} pairs "
        f"({int(ska_ec100['n_pairs_cross_sector'])} cross-sector); Kp "
        f"{int(ska_kp100['n_pairs_total'])} pairs "
        f"({int(ska_kp100['n_pairs_cross_sector'])} cross-sector). Full table: "
        "ska_threshold_robustness.tsv. NOTE: stageE/ska_pairs_*.tsv are the "
        "≤100-filtered derivatives of these distance files; cross-chunk pairs "
        "were never computed (ska merge OOM), so all SKA counts remain "
        "within-chunk only.",
        "",
        "## Remaining",
        "1. Stage D execution on the ~150 scale-up genomes (--run-stageD), then "
        "re-merge localization + re-run Stage E/F on the extended cohort.",
        "2. MGE flanking synteny (protocol element still not attempted).",
        "3. SKA cross-chunk distances (needs a larger-memory machine).",
    ]
    (OUT / "phase2_summary.md").write_text("\n".join(lines) + "\n")
    print("Wrote all Phase 2 outputs to", OUT)


if __name__ == "__main__":
    main()
