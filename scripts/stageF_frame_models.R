#!/usr/bin/env Rscript
# Stage F (metadata level): mixed-effects models on the Stage B sampling frame.
# Data: NCBI Pathogen Detection AMRFinderPlus genotypes (FTP AMR_genotypes) for ALL
# Primary-Cohort genomes in shortlist cells (n=8,987) -- no local compute needed.
# Outcome: carriage of critical-priority genes. Predictors: sector, income group,
# species, collection year; random intercept: country (protocol S4 Stage 6 design).
suppressMessages({library(lme4)})

frame <- read.csv("~/workspace/one_health_amr_pipeline/full_protocol_run/stageB/sampling_frame.csv",
                  stringsAsFactors = FALSE)
has <- function(frame, pattern) grepl(pattern, frame$amr_genes, fixed = FALSE)
frame$any_critical <- as.integer(has(frame, "^blaNDM|^blaKPC|^blaOXA-48|^blaOXA-181|^blaOXA-232|mcr-"))
frame$blaNDM <- as.integer(has(frame, "^blaNDM"))
frame$mcr <- as.integer(has(frame, "mcr-"))
frame$sector <- factor(frame$one_health_sector, levels = c("Human", "Livestock_Poultry", "Wastewater_Aquatic"))
frame$income_group <- factor(frame$income_group, levels = c("HIC", "LMIC", "Unclassified"))
frame$species <- factor(ifelse(grepl("Klebsiella", frame$taxgroup_name), "Kp", "Ec"))
frame$year_c <- frame$collection_year - 2020

out <- file("~/workspace/one_health_amr_pipeline/full_protocol_run/stageF/frame_models_summary.txt", open = "wt")
sink(out, type = "output")
cat("Stage F frame-level mixed models (Primary Cohort frame, NCBI genotypes)\n")
cat("N genomes:", nrow(frame), "\n\n")
for (resp in c("any_critical", "blaNDM", "mcr")) {
  cat("==== outcome:", resp, " prevalence:", mean(frame[[resp]]), "====\n")
  f <- as.formula(paste(resp, "~ sector + income_group + species + year_c + (1|country)"))
  m <- try(glmer(f, data = frame, family = binomial,
                 control = glmerControl(optimizer = "bobyqa")), silent = TRUE)
  if (inherits(m, "try-error")) {
    cat("glmer failed; falling back to glm without random intercept\n")
    m <- glm(update(f, . ~ . - (1|country)), data = frame, family = binomial)
  }
  print(summary(m))
  cat("\nOdds ratios (fixed effects):\n")
  print(round(exp(fixef(m)), 3))
  cat("\n")
}
sink()
cat("done\n")
