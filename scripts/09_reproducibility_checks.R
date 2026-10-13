# 09_reproducibility_checks.R
# Final audit checks and core60 functional over-representation analysis.
#
# This script intentionally focuses on numerical/statistical reproducibility.
# Exact historical figure-drawing source files were not preserved and are not
# fabricated here.

source("scripts/_helpers.R")

out_dir <- p_out("09_repro_checks")
ensure_dir(out_dir)
checks <- list()
add_check <- function(name, pass, observed = NA_character_, expected = NA_character_) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = name, pass = isTRUE(pass), observed = as.character(observed), expected = as.character(expected),
    stringsAsFactors = FALSE
  )
}

# 1. Bulk candidate pathways and positive control from saved executed outputs.
hall_saved <- utils::read.delim(p_data("03_bulk_analysis", "Hallmark_three_cohort_comparison.tsv"),
                                check.names = FALSE, stringsAsFactors = FALSE)
cand_saved <- subset(hall_saved, NES_53624 > 0 & NES_53622 > 0 & NES_TCGA > 0 & n_FDR05 >= 2)$pathway
expected_cand <- c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_ADIPOGENESIS",
                   "HALLMARK_GLYCOLYSIS", "HALLMARK_ESTROGEN_RESPONSE_EARLY")
add_check("bulk_candidate4", setequal(cand_saved, expected_cand), paste(sort(cand_saved), collapse = ";"), paste(sort(expected_cand), collapse = ";"))

pc <- utils::read.csv(p_data("03_bulk_analysis", "tumor_normal_positive_control_summary.csv"))
rho <- pc$value[pc$metric == "Spearman_rho"]
dir <- pc$value[pc$metric == "Direction_concordance"]
add_check("positive_control_rho", abs(rho - 0.921761) < 1e-5, sprintf("%.6f", rho), "0.921761")
add_check("positive_control_direction", abs(dir - 0.873147) < 1e-5, sprintf("%.6f", dir), "0.873147")

# 2. Main epithelial cohort composition.
m52 <- utils::read.csv(p_data("05_scRNA_GSE160269", "GSE160269_Epithelial_main52_analysis_metadata.csv"), stringsAsFactors = FALSE)
add_check("epi52_n", nrow(m52) == 52, nrow(m52), 52)
add_check("epi52_pN0", sum(m52$N_group == "N0") == 28, sum(m52$N_group == "N0"), 28)
add_check("epi52_pNplus", sum(m52$N_group == "Nplus") == 24, sum(m52$N_group == "Nplus"), 24)
add_check("epi52_cells", sum(m52$n_cells) == 44122, sum(m52$n_cells), 44122)

# 3. Leading-edge recurrence/core60.
rob <- utils::read.csv(p_data("05_scRNA_GSE160269", "GSE160269_EstrogenEarly_leading_edge_robustness.csv"))
rec_counts <- table(factor(rob$n_models, levels = 1:4))
add_check("leading_edge_union", nrow(rob) == 77, nrow(rob), 77)
add_check("leading_edge_recurrence_1_2_3_4",
          identical(as.integer(rec_counts), c(7L,4L,6L,60L)),
          paste(as.integer(rec_counts), collapse = ","), "7,4,6,60")
core_saved <- utils::read.csv(p_data("05_scRNA_GSE160269", "GSE160269_pNplus_Epithelial_core60_genes.csv"))$gene
add_check("core60_size", length(core_saved) == 60, length(core_saved), 60)

# 4. Saved executed core60 score model.
fit_core <- readRDS(p_data("05_scRNA_GSE160269", "GSE160269_Epithelial_core60_score_lm_main.rds"))
a_core <- extract_coef(fit_core)
add_check("core60_main_beta", abs(a_core$beta - 0.58) < 0.02, sprintf("%.6f", a_core$beta), "~0.58")
add_check("core60_main_p", abs(log10(a_core$p_value) - log10(7.8e-6)) < 0.2, format(a_core$p_value, scientific = TRUE), "~7.8e-6")

# 5. Saved residual model.
fit_resid <- readRDS(p_data("05_scRNA_GSE160269", "GSE160269_core60_adjusted_for_published_programs.rds"))
a_resid <- extract_coef(fit_resid)
add_check("core60_residual_beta", abs(a_resid$beta - 0.52) < 0.03, sprintf("%.6f", a_resid$beta), "~0.52")
add_check("core60_residual_p", abs(a_resid$p_value - 0.017) < 0.005, sprintf("%.6f", a_resid$p_value), "~0.017")

# 6. Myeloid ESTROGEN_RESPONSE_EARLY saved GSEA.
extract_gsea <- function(path, pathway) {
  x <- as.data.frame(readRDS(path), stringsAsFactors = FALSE)
  idx <- which(trimws(as.character(x$pathway)) == trimws(pathway))
  if (length(idx) != 1L) {
    stop("Pathway not found uniquely in ", path, " (matches=", length(idx), ")")
  }
  x[idx, , drop = FALSE]
}
my_main <- extract_gsea(p_data("05_scRNA_GSE160269", "GSE160269_Myeloid_Hallmark_GSEA_main.rds"), "HALLMARK_ESTROGEN_RESPONSE_EARLY")
my_t4 <- extract_gsea(p_data("05_scRNA_GSE160269", "GSE160269_Myeloid_Hallmark_GSEA_T4fixed.rds"), "HALLMARK_ESTROGEN_RESPONSE_EARLY")
add_check("myeloid_EstrogenEarly_main_NES", round(my_main$NES, 2) == 1.51, round(my_main$NES, 2), 1.51)
add_check("myeloid_EstrogenEarly_main_FDR", round(my_main$padj, 3) == 0.023, round(my_main$padj, 3), 0.023)
add_check("myeloid_EstrogenEarly_T4_NES", round(my_t4$NES, 2) == 1.58, round(my_t4$NES, 2), 1.58)
add_check("myeloid_EstrogenEarly_T4_FDR", round(my_t4$padj, 3) == 0.028, round(my_t4$padj, 3), 0.028)

# 6b. Reconstructed cross-cell-type results, if script 07 has been run.
recon_cross <- p_out("07_cross_celltype", "cross_celltype_EstrogenEarly_summary.csv")
if (file.exists(recon_cross)) {
  rc <- utils::read.csv(recon_cross, stringsAsFactors = FALSE)
  get_rc <- function(ct, model) rc[rc$cell_type == ct & rc$model == model, , drop = FALSE]
  r1 <- get_rc("Myeloid", "T_binary")
  r2 <- get_rc("Myeloid", "T1-T4")
  add_check("reconstructed_myeloid_main_NES", nrow(r1) == 1 && abs(r1$NES - 1.508108) < 1e-5,
            if (nrow(r1)) sprintf("%.6f", r1$NES) else "missing", "1.508108")
  add_check("reconstructed_myeloid_main_FDR", nrow(r1) == 1 && abs(r1$padj - 0.02298321) < 1e-7,
            if (nrow(r1)) sprintf("%.8f", r1$padj) else "missing", "0.02298321")
  add_check("reconstructed_myeloid_T4_NES", nrow(r2) == 1 && abs(r2$NES - 1.580379) < 1e-5,
            if (nrow(r2)) sprintf("%.6f", r2$NES) else "missing", "1.580379")
  add_check("reconstructed_myeloid_T4_FDR", nrow(r2) == 1 && abs(r2$padj - 0.0281) < 1e-7,
            if (nrow(r2)) sprintf("%.7f", r2$padj) else "missing", "0.0281000")
}

# 7. Bulk core57 saved executed effects and REML summary.
eff <- utils::read.csv(p_data("03_bulk_analysis", "core57_standardized_bulk_effects.csv"))
expected_beta <- c(0.007336104, 0.2772910, 0.652862860)
expected_se <- c(0.176069518, 0.2671785, 0.241733479)
add_check("core57_cohort_betas", max(abs(eff$beta - expected_beta)) < 1e-6,
          paste(round(eff$beta, 6), collapse = ","), paste(round(expected_beta, 6), collapse = ","))
add_check("core57_cohort_SE", max(abs(eff$se - expected_se)) < 1e-6,
          paste(round(eff$se, 6), collapse = ","), paste(round(expected_se, 6), collapse = ","))

if (requireNamespace("metafor", quietly = TRUE)) {
  mf <- metafor::rma(yi = beta, sei = se, data = eff, method = "REML")
  add_check("core57_meta_beta", abs(as.numeric(mf$beta) - 0.2866) < 5e-4, sprintf("%.4f", mf$beta), "0.2866")
  add_check("core57_meta_p", abs(mf$pval - 0.1476) < 5e-4, sprintf("%.4f", mf$pval), "0.1476")
  add_check("core57_meta_I2", abs(mf$I2 - 56.58) < 0.1, sprintf("%.2f", mf$I2), "56.58")
}

# 8. Corrected S9 values, if reconstructed script 06 has been run.
s9_new <- p_out("06_published_programs", "Published_program_pN_associations_8programs_BH.csv")
if (file.exists(s9_new)) {
  s9 <- utils::read.csv(s9_new, stringsAsFactors = FALSE)
  expected_fdr <- c(
    Epi1Epi2_set2 = 0.00458520330415088,
    Stress = 0.00795401247132884,
    Mucosal = 0.1659752227915688,
    Epi1Epi2_set1 = 0.1659752227915688,
    AP = 0.3939322644906208,
    Mes = 0.4919916902907933,
    Oxd = 0.6468973439855371,
    Cycling = 0.782895454604577
  )
  obs <- setNames(s9$FDR, s9$Program)[names(expected_fdr)]
  add_check("corrected_S9_BH8", max(abs(obs - expected_fdr), na.rm = TRUE) < 1e-8,
            paste(format(obs, digits = 7), collapse = ","), "final corrected 8-program BH values")
}

# 9. Functional ORA, only if the unfiltered 17,012-gene pseudobulk matrix is available.
counts_file <- p_data("05_scRNA_GSE160269", "GSE160269_Epithelial_main52_pseudobulk_counts.tsv.gz")
if (file.exists(counts_file) && requireNamespace("gprofiler2", quietly = TRUE)) {
  counts <- read_pseudobulk_tsv(counts_file)
  background <- rownames(counts)
  if (length(background) != 17012L) warning("Expected 17,012-gene unfiltered background; observed ", length(background))
  gp <- gprofiler2::gost(
    query = core_saved,
    organism = "hsapiens",
    ordered_query = FALSE,
    significant = TRUE,
    user_threshold = 0.05,
    correction_method = "g_SCS",
    domain_scope = "custom",
    custom_bg = background,
    sources = c("GO:BP", "REAC")
  )
  # Preserve the complete g:Profiler response losslessly as RDS.
  saveRDS(gp, file.path(out_dir, "core60_GO_REACTOME_gprofiler.rds"))

  if (!is.null(gp$result)) {
    # gprofiler2 returns list-columns (for example parents/intersections).
    # Base write.csv() cannot encode those columns; flatten only for the
    # human-readable CSV while retaining the complete object in the RDS above.
    write_flat_csv(
      gp$result,
      file.path(out_dir, "core60_GO_REACTOME_gprofiler.csv"),
      row.names = FALSE
    )
  }
}

check_df <- do.call(rbind, checks)
utils::write.csv(check_df, file.path(out_dir, "reproducibility_checklist.csv"), row.names = FALSE)
print(check_df)
if (any(!check_df$pass)) warning("One or more reproducibility checks failed; inspect output before release.")
