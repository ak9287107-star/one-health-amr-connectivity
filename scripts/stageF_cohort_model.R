#!/usr/bin/env Rscript
# Stage F cohort-level mixed model (Protocol S4 Stage 6):
# P(AMR hit is plasmid-borne) ~ sector + income + species + year + (1|country) + (1|gene)
# Data: stageD/summary/amr_plasmid_localization_table.tsv (all 151 core genomes).
suppressMessages(library(lme4))
d <- read.delim('stageD/summary/amr_plasmid_localization_table.tsv', stringsAsFactors = FALSE)
d <- d[d$molecule_type %in% c('plasmid', 'chromosome'), ]
d$plasmid <- as.integer(d$molecule_type == 'plasmid')
gene_col <- if ('Element.symbol' %in% names(d)) 'Element.symbol' else 'Element symbol'
d$gene <- d[[gene_col]]
d$year_c <- suppressWarnings(as.numeric(d$collection_year)) - 2020
d <- d[!is.na(d$year_c) & d$gene != '' & !is.na(d$gene), ]
cat('hits analysed:', nrow(d), ' from ', length(unique(d$asm_acc)), ' genomes\n')
sink('stageF/cohort_localization_model.txt')
cat('Stage F cohort-level mixed model: plasmid localization (151 core genomes, AMRFinderPlus x MOB-suite)\n')
cat('N hits:', nrow(d), '\n\n')
m <- glmer(plasmid ~ one_health_sector + income_group + taxgroup_name + year_c +
             (1 | country) + (1 | gene),
           data = d, family = binomial,
           control = glmerControl(optimizer = 'bobyqa'))
print(summary(m))
cat('\nOdds ratios (fixed effects):\n')
print(round(exp(fixef(m)), 3))
sink()
cat('model written to stageF/cohort_localization_model.txt\n')
