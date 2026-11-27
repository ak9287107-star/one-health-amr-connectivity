# Pipeline run commands (as executed, Nov 2026)

Environment: micromamba env `oh_amr` (NCBI datasets 18.41.1, AMRFinderPlus 4.2.7,
MOB-suite 3.1.9, mlst 2.35.0, SKA 0.5.1); BUSCO used for assembly QC (CheckM2
unavailable — documented substitution in the manuscript).

## Stage A — cohort frame (NCBI Pathogen Detection snapshots)
- Snapshots: PDG000000004.6463 (E. coli/Shigella), PDG000000012.2587 (Klebsiella)
- `stageA_audit.py` — eligibility audit → 433,603 eligible assemblies, 151 countries;
  Primary Cohort (Complete/Chromosome only) = 10,073 genomes.

## Stage B — stratified sampling frame + core cohort
- `stageB_frame.py` — builds the 8,987-genome sampling frame in shortlist cells and the
  stratified core cohort with fixed seed **20261124**; rank 1 per stratum = the 151-genome
  core cohort (`data/accessions_core151.tsv`).
- `download_genomes.py` — per-accession GenBank FTP download of genome FASTAs.

## Stage C — typing
- `stageC_batch.sh` (mlst), `stageC_clermont.sh` (ClermonTyping, E. coli),
  `stageC_kleborate.sh` (Kleborate, K. pneumoniae).

## Stage D — AMR annotation + plasmid reconstruction (per genome)
Executed by `stageD_run.sh` / `stageD_run2.sh`; core commands:

```
amrfinder -n <genome.fna> -O <organism> --plus --threads 4 -o <acc>_amrfinder.tsv
mob_recon --infile <genome.fna> --outdir <acc>_mob --num_threads 4 --force
```

- `stageD_merge.py` — merges AMRFinderPlus hits with MOB-suite molecule typing into
  `amr_plasmid_localization_table.tsv` (4,948 hits; plasmid/chromosome assignment).

## Stage E — connectivity
- `stageE_ska.sh` / `stageE_ska_batched.sh` (SKA per species, within-chunk),
  `stageE_analysis.py` (Bray–Curtis, L1–L4 connectivity scores, plasmid-cluster sharing).

## Stage F — statistics (R, lme4)
- `stageF_frame_models.R` — frame-level glmer models (n = 8,987).
- `stageF_cohort_model.R` — cohort-level glmer on the 4,948-hit localization table.
- `stageF_table_figures.py` — Table 1 and Figures 1–4.

## Phase 2 — scale-up + robustness (revision)
- `phase2_scaleup_and_robustness.py` — extended 300-genome cohort
  (`data/accessions_scaleup149.tsv` = the 149 staged additions), QC-exclusion
  sensitivity, bootstrap CIs, leave-one-country-out, rarefaction, SKA thresholds.
- `m6_frame_robustness_models.R` — M6a–M6g frame-model robustness suite.
