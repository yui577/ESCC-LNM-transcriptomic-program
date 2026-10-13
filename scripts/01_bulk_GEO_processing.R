# 01_bulk_GEO_processing.R
# Reconstructs the reported pN-associated limma analyses for GSE53624/GSE53622
# from the saved gene-level expression and clinical RDS objects.
#
# PROVENANCE STATUS:
# - Final adjusted limma model commands were recovered from the preserved execution transcript.
# - The original raw-to-gene-level GEO processing script was not preserved.
# - The inputs used here are saved products of the executed analysis.
# - The model specification is locked to the final manuscript Methods.
# - This script is therefore an audited downstream reconstruction, not an
#   original historical source file.

source("scripts/_helpers.R")
if (!requireNamespace("limma", quietly = TRUE)) stop("Package 'limma' is required.")

out_dir <- p_out("01_bulk_GEO")
ensure_dir(out_dir)

run_geo_cohort <- function(cohort) {
  expr_path <- p_data("03_bulk_analysis", paste0(cohort, "_gene_expression_", ifelse(cohort == "GSE53624", "119", "60"), "tumors.rds"))
  clin_path <- p_data("03_bulk_analysis", paste0(cohort, "_clinical_main.rds"))
  assert_file(expr_path); assert_file(clin_path)

  expr <- read_geo_gene_matrix(expr_path)
  clin <- readRDS(clin_path)
  clin <- as.data.frame(clin)
  clin <- align_clinical(expr, clin)

  if (!"N_group" %in% names(clin)) {
    if (!"N_status" %in% names(clin)) stop(cohort, ": no N_group/N_status variable.")
    clin$N_group <- make_n_group(clin$N_status)
  } else {
    clin$N_group <- factor(as.character(clin$N_group), levels = c("N0", "Nplus"))
  }

  required <- c("N_group", "T_stage", "grade", "age", "sex")
  missing_vars <- setdiff(required, names(clin))
  if (length(missing_vars)) stop(cohort, ": missing clinical columns: ", paste(missing_vars, collapse = ", "))

  design <- stats::model.matrix(~ N_group + T_stage + grade + age + sex, data = clin)
  if (!"N_groupNplus" %in% colnames(design)) stop(cohort, ": N_groupNplus coefficient not created.")

  fit <- limma::lmFit(expr, design)
  fit <- limma::eBayes(fit, robust = TRUE)
  tt <- limma::topTable(fit, coef = "N_groupNplus", number = Inf, adjust.method = "BH", sort.by = "P")
  tt$gene <- rownames(tt)

  saveRDS(list(expr = expr, clinical = clin, design = design, fit = fit, table = tt),
          file.path(out_dir, paste0(cohort, "_limma_reconstructed.rds")))
  utils::write.table(tt,
                     file.path(out_dir, paste0(cohort, "_limma_reconstructed.tsv")),
                     sep = "\t", quote = FALSE, row.names = FALSE)

  list(expr = expr, clinical = clin, design = design, fit = fit, table = tt)
}

res_53624 <- run_geo_cohort("GSE53624")
res_53622 <- run_geo_cohort("GSE53622")

cat("GSE53624 samples:", ncol(res_53624$expr), "genes:", nrow(res_53624$expr), "\n")
cat("GSE53622 samples:", ncol(res_53622$expr), "genes:", nrow(res_53622$expr), "\n")
cat("Expected final cohort sizes: 119 and 60.\n")
