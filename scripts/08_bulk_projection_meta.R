# 08_bulk_projection_meta.R
# Cross-cohort projection of the epithelial-derived core60 as harmonized core57,
# followed by REML random-effects meta-analysis.
#
# PROVENANCE STATUS:
# This is the strongest recovered script section. The GSE53622, TCGA-ESCC,
# standardization, regression, and metafor workflow are reconstructed directly
# from the executed R history saved on 2026-10-08. The parallel GSE53624 block
# is reconstructed from the same executed object/formula and final Methods.

source("scripts/_helpers.R")
if (!requireNamespace("metafor", quietly = TRUE)) stop("Package 'metafor' is required.")

out_dir <- p_out("08_bulk_projection")
ensure_dir(out_dir)
core60 <- readRDS(p_out("05_core60", "core60_genes_reconstructed.rds"))

expr24 <- read_geo_gene_matrix(p_data("03_bulk_analysis", "GSE53624_gene_expression_119tumors.rds"))
expr22 <- read_geo_gene_matrix(p_data("03_bulk_analysis", "GSE53622_gene_expression_60tumors.rds"))
clin24 <- as.data.frame(readRDS(p_data("03_bulk_analysis", "GSE53624_clinical_main.rds")))
clin22 <- as.data.frame(readRDS(p_data("03_bulk_analysis", "GSE53622_clinical_main.rds")))
tcga83 <- readRDS(p_data("06_TCGA_ESCC", "TCGA_ESCC_voom_83patients.rds"))

common57 <- Reduce(intersect, list(core60, rownames(expr24), rownames(expr22), rownames(tcga83$E)))
missing3 <- setdiff(core60, common57)
cat("Jointly measurable core genes:", length(common57), "\n")
cat("Not jointly measurable:", paste(missing3, collapse = ", "), "\n")
stopifnot(length(common57) == 57L)
stopifnot(setequal(missing3, c("DHRS2", "ESRP2", "SULT2B1")))

fit_geo <- function(expr, clin, cohort) {
  score_obj <- patient_mean_z_score(expr, common57)
  clin <- align_clinical(expr, clin)
  clin$core_score <- as.numeric(score_obj$score)
  clin$core_score_SD <- as.numeric(scale(clin$core_score))
  if (!"N_group" %in% names(clin)) clin$N_group <- make_n_group(clin$N_status)
  clin$N_group <- factor(as.character(clin$N_group), levels = c("N0", "Nplus"))
  required <- c("T_stage", "grade", "age", "sex")
  miss <- setdiff(required, names(clin))
  if (length(miss)) stop(cohort, ": missing covariates: ", paste(miss, collapse = ", "))
  fit <- stats::lm(core_score_SD ~ N_group + T_stage + grade + age + sex, data = clin)
  a <- extract_coef(fit)
  data.frame(cohort = cohort, beta = a$beta, se = a$se, p_value = a$p_value,
             ci_low = a$ci_low, ci_high = a$ci_high, stringsAsFactors = FALSE)
}

r24 <- fit_geo(expr24, clin24, "GSE53624")
r22 <- fit_geo(expr22, clin22, "GSE53622")

# TCGA block mirrors the executed 2026-10-08 history and uses the saved design.
score_tc <- patient_mean_z_score(tcga83$E, common57)$score
N_group_tcga <- factor(ifelse(tcga83$design[, "N_groupNplus"] == 1, "Nplus", "N0"),
                       levels = c("N0", "Nplus"))
analysis_tcga <- data.frame(
  core_score = as.numeric(score_tc),
  N_group = N_group_tcga,
  T_stageT2 = tcga83$design[, "T_stageT2"],
  T_stageT3 = tcga83$design[, "T_stageT3"],
  T_stageT4 = tcga83$design[, "T_stageT4"],
  grademoderately = tcga83$design[, "grademoderately"],
  gradepoorly = tcga83$design[, "gradepoorly"],
  age_years = tcga83$design[, "age_years"],
  sexmale = tcga83$design[, "sexmale"]
)
analysis_tcga$core_score_SD <- as.numeric(scale(analysis_tcga$core_score))
fit_tc <- stats::lm(core_score_SD ~ N_group + T_stageT2 + T_stageT3 + T_stageT4 +
                      grademoderately + gradepoorly + age_years + sexmale,
                    data = analysis_tcga)
a_tc <- extract_coef(fit_tc)
rTC <- data.frame(cohort = "TCGA-ESCC", beta = a_tc$beta, se = a_tc$se, p_value = a_tc$p_value,
                  ci_low = a_tc$ci_low, ci_high = a_tc$ci_high, stringsAsFactors = FALSE)

cohort_effects <- rbind(r24, r22, rTC)
utils::write.csv(cohort_effects, file.path(out_dir, "core57_standardized_bulk_effects_reconstructed.csv"), row.names = FALSE)

meta_fit <- metafor::rma(yi = beta, sei = se, data = cohort_effects, method = "REML")
saveRDS(meta_fit, file.path(out_dir, "core57_random_effects_meta_reconstructed.rds"))
capture.output(meta_fit, file = file.path(out_dir, "core57_random_effects_meta_reconstructed.txt"))

meta_summary <- data.frame(
  estimate = as.numeric(meta_fit$beta),
  se = meta_fit$se,
  p_value = meta_fit$pval,
  ci_low = meta_fit$ci.lb,
  ci_high = meta_fit$ci.ub,
  tau2 = meta_fit$tau2,
  I2 = meta_fit$I2,
  Q = meta_fit$QE,
  Q_p = meta_fit$QEp
)
utils::write.csv(meta_summary, file.path(out_dir, "core57_random_effects_meta_summary.csv"), row.names = FALSE)

print(cohort_effects)
print(meta_summary)
