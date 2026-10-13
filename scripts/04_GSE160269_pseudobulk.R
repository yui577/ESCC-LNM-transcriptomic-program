# 04_GSE160269_pseudobulk.R
# GSE160269 epithelial patient-level pseudobulk sensitivity analyses.
#
# PROVENANCE STATUS (recovered workflow):
# Key model/GSEA commands were recovered from the user's prior 2026-10-08 chat,
# not inferred only from the manuscript. In particular, the original epithelial
# limma fits used robust empirical Bayes, and the Hallmark sets were read from
# the local MSigDB 2025.1 GMT file.

source("scripts/_helpers.R")

out_dir <- p_out("04_GSE160269_epithelial")
ensure_dir(out_dir)

# Recovered original command:
# hallmark_list_160269 <- fgsea::gmtPathways(
#   "05_scRNA_GSE160269/h.all.v2025.1.Hs.symbols.gmt"
# )
hallmark_list_160269 <- hallmark_sets_gse160269()
stopifnot(length(hallmark_list_160269) == 50L)

prepare_meta <- function(path) {
  assert_file(path)
  m <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"N_group" %in% names(m)) {
    if (!"N_stage" %in% names(m)) stop("Metadata lacks N_group/N_stage.")
    m$N_group <- make_n_group(m$N_stage)
  } else {
    m$N_group <- factor(m$N_group, levels = c("N0", "Nplus"))
  }
  m$T_stage4 <- make_t_stage4(m$T_stage)
  m$T_group <- make_t_group(m$T_stage)
  m
}

# Build one filtered/TMM DGEList per patient cohort using the collapsed-T design.
# This matches the recovered original 52-patient code:
#   design_epi52 <- model.matrix(~ T_group + N_group, data = meta_epi52)
#   y_epi52 <- DGEList(counts = counts_epi52)
#   keep_epi52 <- filterByExpr(y_epi52, design = design_epi52)
#   y_epi52 <- y_epi52[keep_epi52, , keep.lib.sizes = FALSE]
#   y_epi52 <- calcNormFactors(y_epi52, method = "TMM")
prepare_one_cohort <- function(count_path, meta_path) {
  counts <- read_pseudobulk_tsv(count_path)
  meta <- prepare_meta(meta_path)
  prep <- prepare_pseudobulk_dge(counts, meta, ~ T_group + N_group)
  prep
}

rank_from_fit <- function(fit_obj) {
  r <- fit_obj$fit$t[, "N_groupNplus"]
  names(r) <- rownames(fit_obj$fit$t)
  r <- r[is.finite(r)]
  sort(r, decreasing = TRUE)
}

save_gsea <- function(gsea, prefix) {
  saveRDS(gsea, file.path(out_dir, paste0(prefix, ".rds")))
  write_flat_csv(gsea, file.path(out_dir, paste0(prefix, ".csv")), row.names = FALSE)
}

run_epi52 <- function(count_path, meta_path) {
  prep <- prepare_one_cohort(count_path, meta_path)
  y_epi52 <- prep$y
  meta_epi52 <- prep$meta

  # Recovered original primary 52-patient commands:
  # design_epi52 <- model.matrix(~ T_group + N_group, data = meta_epi52)
  # v_epi52 <- voom(y_epi52, design_epi52, plot = TRUE)
  # fit_epi52 <- lmFit(v_epi52, design_epi52)
  # fit_epi52 <- eBayes(fit_epi52, trend = FALSE, robust = TRUE)
  fit_epi52 <- run_limma_from_dge(
    y_epi52, meta_epi52, ~ T_group + N_group,
    robust = TRUE, trend = FALSE, voom_plot = FALSE
  )
  rank_epi52 <- rank_from_fit(fit_epi52)

  # Recovered original code had set.seed(123) immediately before this call.
  set.seed(123)
  gsea_epi52 <- fgsea::fgseaMultilevel(
    pathways = hallmark_list_160269,
    stats = rank_epi52,
    minSize = 10,
    maxSize = 500
  )
  gsea_epi52 <- gsea_epi52[order(gsea_epi52$padj), ]

  # Recovered original T4 sensitivity commands reused the SAME y_epi52:
  # design_epi52_T4 <- model.matrix(~ T_stage4 + N_group, data=meta_epi52)
  # v_epi52_T4 <- voom(y_epi52, design_epi52_T4)
  # fit_epi52_T4 <- lmFit(v_epi52_T4, design_epi52_T4)
  # fit_epi52_T4 <- eBayes(fit_epi52_T4, robust = TRUE)
  fit_epi52_T4 <- run_limma_from_dge(
    y_epi52, meta_epi52, ~ T_stage4 + N_group,
    robust = TRUE, trend = FALSE, voom_plot = FALSE
  )
  rank_epi52_T4 <- rank_from_fit(fit_epi52_T4)

  # In the recovered preserved execution transcript there was no new set.seed() immediately before
  # gsea_epi52_T4; preserve the RNG stream after gsea_epi52.
  gsea_epi52_T4 <- fgsea::fgseaMultilevel(
    pathways = hallmark_list_160269,
    stats = rank_epi52_T4,
    minSize = 10,
    maxSize = 500
  )
  gsea_epi52_T4 <- gsea_epi52_T4[order(gsea_epi52_T4$padj), ]

  saveRDS(fit_epi52, file.path(out_dir, "GSE160269_Epithelial_main52_limma.rds"))
  saveRDS(fit_epi52_T4, file.path(out_dir, "GSE160269_Epithelial_main52_limma_T4fixed.rds"))
  save_gsea(gsea_epi52, "GSE160269_Epithelial_main52_Hallmark_GSEA")
  save_gsea(gsea_epi52_T4, "GSE160269_Epithelial_main52_Hallmark_GSEA_T4fixed")

  list(y = y_epi52, meta = meta_epi52,
       fit = fit_epi52, fit_T4 = fit_epi52_T4,
       gsea = gsea_epi52, gsea_T4 = gsea_epi52_T4)
}

run_epi55 <- function(count_path, meta_path) {
  prep <- prepare_one_cohort(count_path, meta_path)
  y_epi55 <- prep$y
  meta_epi55 <- prep$meta

  # The preserved execution transcript confirms the 55-patient analyses used the same core workflow
  # and set.seed(123) before both fgsea calls. Robust eBayes is applied here to
  # mirror the recovered epithelial analysis convention.
  fit_epi55 <- run_limma_from_dge(
    y_epi55, meta_epi55, ~ T_group + N_group,
    robust = TRUE, trend = FALSE, voom_plot = FALSE
  )
  rank_epi55 <- rank_from_fit(fit_epi55)
  set.seed(123)
  gsea_epi55 <- fgsea::fgseaMultilevel(
    pathways = hallmark_list_160269,
    stats = rank_epi55,
    minSize = 10,
    maxSize = 500
  )
  gsea_epi55 <- gsea_epi55[order(gsea_epi55$padj), ]

  fit_epi55_T4 <- run_limma_from_dge(
    y_epi55, meta_epi55, ~ T_stage4 + N_group,
    robust = TRUE, trend = FALSE, voom_plot = FALSE
  )
  rank_epi55_T4 <- rank_from_fit(fit_epi55_T4)
  set.seed(123)
  gsea_epi55_T4 <- fgsea::fgseaMultilevel(
    pathways = hallmark_list_160269,
    stats = rank_epi55_T4,
    minSize = 10,
    maxSize = 500
  )
  gsea_epi55_T4 <- gsea_epi55_T4[order(gsea_epi55_T4$padj), ]

  saveRDS(fit_epi55, file.path(out_dir, "GSE160269_Epithelial_sensitivity55_limma.rds"))
  saveRDS(fit_epi55_T4, file.path(out_dir, "GSE160269_Epithelial_sensitivity55_limma_T4fixed.rds"))
  save_gsea(gsea_epi55, "GSE160269_Epithelial_sensitivity55_Hallmark_GSEA")
  save_gsea(gsea_epi55_T4, "GSE160269_Epithelial_sensitivity55_Hallmark_GSEA_T4fixed")

  list(y = y_epi55, meta = meta_epi55,
       fit = fit_epi55, fit_T4 = fit_epi55_T4,
       gsea = gsea_epi55, gsea_T4 = gsea_epi55_T4)
}

main_counts <- p_data("05_scRNA_GSE160269", "GSE160269_Epithelial_main52_pseudobulk_counts.tsv.gz")
main_meta   <- p_data("05_scRNA_GSE160269", "GSE160269_Epithelial_main52_patient_metadata.csv")
sens_counts <- p_data("05_scRNA_GSE160269", "GSE160269_Epithelial_sensitivity55_pseudobulk_counts.tsv.gz")
sens_meta   <- p_data("05_scRNA_GSE160269", "GSE160269_Epithelial_sensitivity55_patient_metadata.csv")

if (all(file.exists(c(main_counts, main_meta, sens_counts, sens_meta)))) {
  epi52 <- run_epi52(main_counts, main_meta)
  epi55 <- run_epi55(sens_counts, sens_meta)

  # Quick locked-result printout for direct comparison with the executed study.
  target <- "HALLMARK_ESTROGEN_RESPONSE_EARLY"
  pick <- function(x) x[x$pathway == target, c("NES", "pval", "padj", "size", "leadingEdge")]
  cat("\n52 / T-binary:\n"); print(pick(epi52$gsea))
  cat("\n52 / T1-T4:\n"); print(pick(epi52$gsea_T4))
  cat("\n55 / T-binary:\n"); print(pick(epi55$gsea))
  cat("\n55 / T1-T4:\n"); print(pick(epi55$gsea_T4))

  cat("Recomputed four epithelial analyses using recovered original model details.\n")
} else {
  message("Patient-level pseudobulk count files were not all found. Staging saved executed GSEA objects instead.")
  saved <- c(
    "GSE160269_Epithelial_main52_Hallmark_GSEA.rds",
    "GSE160269_Epithelial_main52_Hallmark_GSEA_T4fixed.rds",
    "GSE160269_Epithelial_sensitivity55_Hallmark_GSEA.rds",
    "GSE160269_Epithelial_sensitivity55_Hallmark_GSEA_T4fixed.rds"
  )
  for (nm in saved) {
    src <- p_data("05_scRNA_GSE160269", nm)
    assert_file(src)
    file.copy(src, file.path(out_dir, nm), overwrite = TRUE)
  }
  cat("Fallback completed: saved executed GSEA objects staged; no upstream re-analysis performed.\n")
}
