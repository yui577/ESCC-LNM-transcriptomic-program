# 07_cross_celltype_analysis.R
# Cross-cell-type GSE160269 pseudobulk analysis.
#
# PROVENANCE STATUS: command-level reconstruction from the prior executed chat
# workflow, checked against the saved original GSEA RDS objects.
# Key recovered details:
#   * main design: ~ T_group + N_group
#   * filterByExpr performed once with the main design
#   * same filtered/TMM-normalized DGEList reused for T-stage-4 model
#   * voom -> lmFit -> eBayes() (default robust = FALSE)
#   * rank taken directly from fit$t[, "N_groupNplus"]
#   * local Hallmark GMT h.all.v2025.1.Hs.symbols.gmt
#   * set.seed(20261011)
#   * Hallmark fgseaMultilevel minSize=15, maxSize=500, eps=0
#   * core60 fgseaMultilevel minSize=10, maxSize=500, eps=0

source("scripts/_helpers.R")

if (!requireNamespace("edgeR", quietly = TRUE)) stop("Package 'edgeR' is required.")
if (!requireNamespace("limma", quietly = TRUE)) stop("Package 'limma' is required.")
if (!requireNamespace("fgsea", quietly = TRUE)) stop("Package 'fgsea' is required.")

out_dir <- p_out("07_cross_celltype")
ensure_dir(out_dir)

hallmark_gmt <- fgsea::gmtPathways(
  p_data("05_scRNA_GSE160269", "h.all.v2025.1.Hs.symbols.gmt")
)
core60 <- readRDS(p_out("05_core60", "core60_genes_reconstructed.rds"))

clinical <- utils::read.csv(
  p_data("05_scRNA_GSE160269", "GSE160269_patient_clinical_Nstage.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
clinical$N_group <- factor(clinical$N_group, levels = c("N0", "Nplus"))
clinical$T_stage4 <- make_t_stage4(clinical$T_stage)
clinical$T_group <- factor(
  ifelse(clinical$T_stage %in% c("T1", "T2"), "Early_T", "Advanced_T"),
  levels = c("Early_T", "Advanced_T")
)

as_count_matrix <- function(obj) {
  if (inherits(obj, "DGEList")) obj <- obj$counts
  if (is.list(obj) && !is.null(obj$counts)) obj <- obj$counts
  m <- as.matrix(obj)
  storage.mode(m) <- "numeric"
  m
}

flat_write <- function(x, path) {
  write_flat_csv(x, path, row.names = FALSE)
}

extract_one <- function(x, pathway) {
  d <- as.data.frame(x, stringsAsFactors = FALSE)
  idx <- which(trimws(as.character(d$pathway)) == trimws(pathway))
  if (length(idx) != 1L) {
    stop("Pathway not found uniquely: ", pathway,
         " (matches=", length(idx), ")")
  }
  d[idx, , drop = FALSE]
}

run_celltype_exact <- function(label, filename) {
  counts <- as_count_matrix(readRDS(p_data("05_scRNA_GSE160269", filename)))
  meta <- clinical[match(colnames(counts), clinical$patient_id), , drop = FALSE]
  if (anyNA(meta$patient_id)) stop(label, ": unable to match patient IDs")
  stopifnot(identical(as.character(meta$patient_id), colnames(counts)))

  # ---------- main collapsed-T model ----------
  design_main <- stats::model.matrix(~ T_group + N_group, data = meta)

  dge <- edgeR::DGEList(counts = counts)
  keep <- edgeR::filterByExpr(dge, design = design_main)
  dge <- dge[keep, , keep.lib.sizes = FALSE]

  if ("normLibSizes" %in% getNamespaceExports("edgeR")) {
    dge <- edgeR::normLibSizes(dge, method = "TMM")
  } else {
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
  }

  v_main <- limma::voom(dge, design_main, plot = FALSE)
  fit_main <- limma::lmFit(v_main, design_main)
  fit_main <- limma::eBayes(fit_main)

  rank_main <- sort(
    stats::setNames(
      fit_main$t[, "N_groupNplus"],
      rownames(fit_main$t)
    ),
    decreasing = TRUE
  )

  set.seed(20261011)
  gsea_main <- fgsea::fgseaMultilevel(
    pathways = hallmark_gmt,
    stats = rank_main,
    minSize = 15,
    maxSize = 500,
    eps = 0
  )

  # ---------- four-level T model ----------
  # IMPORTANT: reuse the same already-filtered/TMM-normalized DGEList.
  design_t4 <- stats::model.matrix(~ T_stage4 + N_group, data = meta)
  v_t4 <- limma::voom(dge, design_t4, plot = FALSE)
  fit_t4 <- limma::lmFit(v_t4, design_t4)
  fit_t4 <- limma::eBayes(fit_t4)

  rank_t4 <- sort(
    stats::setNames(
      fit_t4$t[, "N_groupNplus"],
      rownames(fit_t4$t)
    ),
    decreasing = TRUE
  )

  set.seed(20261011)
  gsea_t4 <- fgsea::fgseaMultilevel(
    pathways = hallmark_gmt,
    stats = rank_t4,
    minSize = 15,
    maxSize = 500,
    eps = 0
  )

  # ---------- focused core60 projection ----------
  core_main <- list(core60 = intersect(core60, names(rank_main)))
  core_t4 <- list(core60 = intersect(core60, names(rank_t4)))

  set.seed(20261011)
  core_gsea_main <- fgsea::fgseaMultilevel(
    pathways = core_main,
    stats = rank_main,
    minSize = 10,
    maxSize = 500,
    eps = 0
  )

  set.seed(20261011)
  core_gsea_t4 <- fgsea::fgseaMultilevel(
    pathways = core_t4,
    stats = rank_t4,
    minSize = 10,
    maxSize = 500,
    eps = 0
  )

  # Save reconstructed objects.
  saveRDS(fit_main, file.path(out_dir, paste0(label, "_limma_main_fit.rds")))
  saveRDS(fit_t4, file.path(out_dir, paste0(label, "_limma_T4fixed_fit.rds")))
  saveRDS(gsea_main, file.path(out_dir, paste0(label, "_Hallmark_GSEA_main.rds")))
  saveRDS(gsea_t4, file.path(out_dir, paste0(label, "_Hallmark_GSEA_T4fixed.rds")))
  saveRDS(core_gsea_main, file.path(out_dir, paste0(label, "_core60_GSEA_main.rds")))
  saveRDS(core_gsea_t4, file.path(out_dir, paste0(label, "_core60_GSEA_T4fixed.rds")))
  flat_write(gsea_main, file.path(out_dir, paste0(label, "_Hallmark_GSEA_main.csv")))
  flat_write(gsea_t4, file.path(out_dir, paste0(label, "_Hallmark_GSEA_T4fixed.csv")))
  flat_write(core_gsea_main, file.path(out_dir, paste0(label, "_core60_GSEA_main.csv")))
  flat_write(core_gsea_t4, file.path(out_dir, paste0(label, "_core60_GSEA_T4fixed.csv")))

  target <- "HALLMARK_ESTROGEN_RESPONSE_EARLY"
  target_main <- extract_one(gsea_main, target)
  target_t4 <- extract_one(gsea_t4, target)

  diagnostics <- data.frame(
    cell_type = label,
    raw_genes = nrow(counts),
    filtered_genes = nrow(dge),
    estrogen_total_genes = length(hallmark_gmt[[target]]),
    estrogen_measurable = sum(hallmark_gmt[[target]] %in% names(rank_main)),
    n_patients = ncol(counts),
    stringsAsFactors = FALSE
  )

  # Hard sanity check recovered from the original Myeloid execution.
  if (label == "Myeloid") {
    stopifnot(nrow(counts) == 15175L)
    stopifnot(nrow(dge) == 11084L)
    stopifnot(sum(hallmark_gmt[[target]] %in% names(rank_main)) == 120L)
  }

  list(
    label = label,
    n = ncol(counts),
    diagnostics = diagnostics,
    g_main = gsea_main,
    g_t4 = gsea_t4,
    c_main = core_gsea_main,
    c_t4 = core_gsea_t4,
    n_core_main = length(core_main$core60),
    n_core_t4 = length(core_t4$core60),
    target_main = target_main,
    target_t4 = target_t4
  )
}

results <- list(
  Myeloid = run_celltype_exact("Myeloid", "GSE160269_Myeloid_60patients_pseudobulk_counts.rds"),
  Tcell = run_celltype_exact("Tcell", "GSE160269_Tcell_60patients_pseudobulk_counts.rds"),
  Bcell = run_celltype_exact("Bcell", "GSE160269_Bcell_56patients_pseudobulk_counts.rds"),
  Fibroblast = run_celltype_exact("Fibroblast", "GSE160269_Fibroblast_49patients_pseudobulk_counts.rds")
)

diagnostics <- do.call(rbind, lapply(results, `[[`, "diagnostics"))
utils::write.csv(diagnostics, file.path(out_dir, "cross_celltype_filter_diagnostics.csv"), row.names = FALSE)

hallmark_summary <- do.call(rbind, lapply(results, function(x) {
  a <- x$target_main
  b <- x$target_t4
  rbind(
    data.frame(cell_type = x$label, model = "T_binary", n_patients = x$n,
               NES = a$NES, pval = a$pval, padj = a$padj, size = a$size),
    data.frame(cell_type = x$label, model = "T1-T4", n_patients = x$n,
               NES = b$NES, pval = b$pval, padj = b$padj, size = b$size)
  )
}))
utils::write.csv(hallmark_summary, file.path(out_dir, "cross_celltype_EstrogenEarly_summary.csv"), row.names = FALSE)

core_summary <- do.call(rbind, lapply(results, function(x) {
  a <- extract_one(x$c_main, "core60")
  b <- extract_one(x$c_t4, "core60")
  rbind(
    data.frame(cell_type = x$label, model = "T_binary", n_patients = x$n,
               n_core_measurable = x$n_core_main,
               NES = a$NES, pval = a$pval, padj = a$padj, size = a$size),
    data.frame(cell_type = x$label, model = "T1-T4", n_patients = x$n,
               n_core_measurable = x$n_core_t4,
               NES = b$NES, pval = b$pval, padj = b$padj, size = b$size)
  )
}))
utils::write.csv(core_summary, file.path(out_dir, "cross_celltype_core60_summary.csv"), row.names = FALSE)

# Optional local audit against the preserved original saved outputs.
audit_rows <- list()
for (label in names(results)) {
  original_prefix <- switch(
    label,
    Myeloid = "GSE160269_Myeloid",
    Tcell = "GSE160269_Tcell",
    Bcell = "GSE160269_Bcell",
    Fibroblast = "GSE160269_Fibroblast"
  )
  for (model in c("main", "T4fixed")) {
    original_file <- p_data(
      "05_scRNA_GSE160269",
      paste0(original_prefix, "_Hallmark_GSEA_", model, ".rds")
    )
    if (file.exists(original_file)) {
      old <- extract_one(readRDS(original_file), "HALLMARK_ESTROGEN_RESPONSE_EARLY")
      new <- if (model == "main") results[[label]]$target_main else results[[label]]$target_t4
      audit_rows[[length(audit_rows) + 1L]] <- data.frame(
        cell_type = label,
        model = model,
        old_NES = old$NES,
        new_NES = new$NES,
        old_padj = old$padj,
        new_padj = new$padj,
        old_size = old$size,
        new_size = new$size,
        delta_NES = new$NES - old$NES,
        delta_padj = new$padj - old$padj,
        stringsAsFactors = FALSE
      )
    }
  }
}
if (length(audit_rows)) {
  audit_compare <- do.call(rbind, audit_rows)
  utils::write.csv(audit_compare, file.path(out_dir, "saved_vs_reconstructed_EstrogenEarly.csv"), row.names = FALSE)
  cat("\n===== FILTER DIAGNOSTICS =====\n")
  print(diagnostics)
  cat("\n===== ESTROGEN_RESPONSE_EARLY =====\n")
  print(hallmark_summary)
  cat("\n===== SAVED vs RECONSTRUCTED =====\n")
  print(audit_compare)
} else {
  print(diagnostics)
  print(hallmark_summary)
}

cat("\n===== CORE60 CROSS-CELLTYPE =====\n")
print(core_summary)
