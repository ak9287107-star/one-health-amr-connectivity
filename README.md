# Plasmid Localization and Last-Resort Antimicrobial Resistance Across Human, Livestock, and Wastewater Sectors: A Complete-Genome Study from 50 Countries

**Authors:** Muhammad Ahmad Khan, Malik Muhammad Osama, Abdul Hakim
Centre of Biotechnology and Microbiology (COBAM), University of Peshawar, Pakistan

LaTeX source of the preprint-ready manuscript, with figures, supplementary Table S1,
genome accession lists, and the full analysis scripts.

**DOI (v1.0.0):** [10.5281/zenodo.24589690](https://doi.org/10.5281/zenodo.24589690) — Zenodo record: https://zenodo.org/records/24589690
**Concept DOI (all versions):** [10.5281/zenodo.24589689](https://doi.org/10.5281/zenodo.24589689)

## Contents

- `main.tex`, `references.bib` — manuscript source
- `figures/fig1.png` – `fig3.png` — manuscript figures
- `TableS1_per_plasmid_cargo_n150.tsv` — Supplementary Table S1
- `data/accessions_core151.tsv` — the 151-genome core analysis cohort (accessions + metadata)
- `data/accessions_scaleup149.tsv` — the 149 staged scale-up genomes (Phase 2 revision cohort, total 300)
- `scripts/` — all analysis code, with `scripts/RUN_COMMANDS.md` documenting the exact
  run commands (sampling seed 20261124, AMRFinderPlus and MOB-suite invocations, R models)

## Data sources

NCBI Pathogen Detection snapshots PDG000000004.6463 and PDG000000012.2587
(retrieved 23 Nov 2026); World Bank country income classifications (API, 24 Nov 2026).
All genomes are publicly available under the accessions listed in `data/`.

## License

MIT — see `LICENSE`.
