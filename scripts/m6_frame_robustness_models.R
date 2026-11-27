#!/usr/bin/env Rscript
# M6 — Frame-level robustness models for the One Health AMR paper (Phase 2).
#
# Run from the project root:
#   Rscript phase2_revision_outputs/m6_frame_robustness_models.R
# Paths: the project root is resolved as the PARENT of this script's own
# directory (via commandArgs --file), so the script also works from any cwd.
# Data: full_protocol_run/stageB/sampling_frame.csv (n=8,987).
# Outcome definitions and factor codings are copied EXACTLY from
# full_protocol_run/scripts/stageF_frame_models.R for comparability:
#   any_critical = grepl("^blaNDM|^blaKPC|^blaOXA-48|^blaOXA-181|^blaOXA-232|mcr-", amr_genes)
#   blaNDM = grepl("^blaNDM", amr_genes);  mcr = grepl("mcr-", amr_genes)
#   sector: Human (ref) / Livestock_Poultry / Wastewater_Aquatic
#   income_group: HIC (ref) / LMIC / Unclassified; species: Ec (ref) / Kp
#   year_c = collection_year - 2020; random intercept: country.
# Models: M6a base glmer (replication of Stage F) | M6b + sector:income_group
#   M6c LMIC-only | M6d HIC-only (income term dropped in both subsets — it is
#   constant there) | M6e 2020-2026 only | M6f drop top-3 countries by frame N
#   M6g glm, no random intercept. blaNDM/mcr: M6a spec only; separation or
#   non-convergence is reported plainly and no ORs are forced.
# OR CIs are computed manually from fixef + vcov (lme4 only, no broom).
suppressMessages(library(lme4))

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
script_dir <- if (length(file_arg)) dirname(normalizePath(file_arg)) else getwd()
project_root <- dirname(script_dir)
out_dir <- file.path(project_root, "phase2_revision_outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

frame <- read.csv(file.path(project_root, "full_protocol_run", "stageB",
                            "sampling_frame.csv"), stringsAsFactors = FALSE)
has <- function(p) grepl(p, frame$amr_genes, fixed = FALSE)
frame$any_critical <- as.integer(has("^blaNDM|^blaKPC|^blaOXA-48|^blaOXA-181|^blaOXA-232|mcr-"))
frame$blaNDM <- as.integer(has("^blaNDM"))
frame$mcr <- as.integer(has("mcr-"))
frame$sector <- factor(frame$one_health_sector,
                       levels = c("Human", "Livestock_Poultry", "Wastewater_Aquatic"))
frame$income_group <- factor(frame$income_group,
                             levels = c("HIC", "LMIC", "Unclassified"))
frame$species <- factor(ifelse(grepl("Klebsiella", frame$taxgroup_name), "Kp", "Ec"))
frame$year_c <- frame$collection_year - 2020

tidy_rows <- list()
notes <- character()
add_tidy <- function(model_label, outcome, fit) {
  fe <- if (inherits(fit, "merMod")) fixef(fit) else coef(fit)
  V <- as.matrix(vcov(fit))
  se <- sqrt(diag(V))
  for (nm in names(fe)) {
    if (is.na(fe[nm]) || is.na(se[nm])) next
    z <- fe[nm] / se[nm]
    tidy_rows[[length(tidy_rows) + 1]] <<- data.frame(
      model = model_label, outcome = outcome, term = nm,
      OR = round(exp(fe[nm]), 4),
      CI_lo = round(exp(fe[nm] - 1.96 * se[nm]), 4),
      CI_hi = round(exp(fe[nm] + 1.96 * se[nm]), 4),
      p = signif(2 * pnorm(-abs(z)), 4), stringsAsFactors = FALSE)
  }
}

fit_one <- function(label, outcome, formula, data, mixed = TRUE,
                    sink_summary = TRUE) {
  fit <- tryCatch(
    if (mixed) glmer(formula, data = data, family = binomial,
                     control = glmerControl(optimizer = "bobyqa"))
    else glm(formula, data = data, family = binomial),
    error = function(e) e)
  if (inherits(fit, "error")) {
    notes <<- c(notes, sprintf("%s [%s]: FIT FAILED — %s", label, outcome,
                               conditionMessage(fit)))
    cat(sprintf("\n#### %s — outcome %s: FIT FAILED — %s\n", label, outcome,
                conditionMessage(fit)))
    return(invisible(NULL))
  }
  # separation / non-identifiability check: huge coefficients or SEs
  fe <- if (inherits(fit, "merMod")) fixef(fit) else coef(fit)
  se <- sqrt(diag(as.matrix(vcov(fit))))
  separated <- any(abs(fe) > 10, na.rm = TRUE) || any(se > 10, na.rm = TRUE)
  conv_ok <- TRUE
  if (inherits(fit, "merMod")) {
    msgs <- fit@optinfo$conv$lme4$messages
    conv_ok <- is.null(msgs) || !any(grepl("failed to converge|degenerate",
                                          msgs))
  }
  cat(sprintf("\n#### %s — outcome %s (n=%d, events=%d)\n",
              label, outcome, nrow(data), sum(data[[outcome]])))
  if (sink_summary) print(summary(fit))
  if (separated || !conv_ok) {
    notes <<- c(notes, sprintf(
      "%s [%s]: fitted but NOT INTERPRETABLE (separation=%s, converged=%s); ORs withheld.",
      label, outcome, separated, conv_ok))
    cat("** NOT INTERPRETABLE (separation / non-convergence) — ORs withheld. **\n")
    return(invisible(NULL))
  }
  add_tidy(label, outcome, fit)
  invisible(fit)
}

base_rhs <- "sector + income_group + species + year_c"
base_rhs_noinc <- "sector + species + year_c"
top3 <- names(sort(table(frame$country), decreasing = TRUE))[1:3]

out_con <- file(file.path(out_dir, "m6_frame_robustness_models.txt"),
                open = "wt")
sink(out_con, type = "output")
cat("M6 frame-level robustness models (Phase 2 revision)\n")
cat("Frame N:", nrow(frame), " genomes;", length(unique(frame$country)),
    "countries\n")
cat("Outcome prevalences: any_critical", mean(frame$any_critical),
    "| blaNDM", mean(frame$blaNDM), "| mcr", mean(frame$mcr), "\n")
cat("M6f drops top-3 countries by frame N:", paste(top3, collapse = ", "), "\n")

f_base <- as.formula(paste("any_critical ~", base_rhs, "+ (1|country)"))
fit_one("M6a", "any_critical", f_base, frame)
fit_one("M6b", "any_critical",
        as.formula(paste("any_critical ~ sector * income_group + species +",
                         "year_c + (1|country)")), frame)
fit_one("M6c", "any_critical",
        as.formula(paste("any_critical ~", base_rhs_noinc, "+ (1|country)")),
        frame[frame$income_group == "LMIC", ])
fit_one("M6d", "any_critical",
        as.formula(paste("any_critical ~", base_rhs_noinc, "+ (1|country)")),
        frame[frame$income_group == "HIC", ])
fit_one("M6e", "any_critical", f_base,
        frame[frame$collection_year >= 2020 & frame$collection_year <= 2026, ])
fit_one("M6f", "any_critical", f_base, frame[!(frame$country %in% top3), ])
fit_one("M6g", "any_critical",
        as.formula(paste("any_critical ~", base_rhs)), frame, mixed = FALSE)

for (resp in c("blaNDM", "mcr")) {
  fit_one("M6a", resp,
          as.formula(paste(resp, "~", base_rhs, "+ (1|country)")), frame)
}
cat("\n==== NOTES ====\n")
for (n in notes) cat("-", n, "\n")
sink()
close(out_con)

tidy <- do.call(rbind, tidy_rows)
write.table(tidy, file.path(out_dir, "m6_or_comparison.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("M6 done:", nrow(tidy), "tidy OR rows written\n")
