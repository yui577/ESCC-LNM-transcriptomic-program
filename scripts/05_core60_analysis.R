# 05_core60_analysis.R
# Derives the recurrent ESTROGEN_RESPONSE_EARLY leading-edge core60 and
# reconstructs the patient-level epithelial core60 score analyses.
#
# PROVENANCE STATUS:
# - Core60 is reconstructed directly from the four saved/recomputed fgsea objects.
# - The patient score uses the recovered original unfiltered-pseudobulk TMM/logCPM workflow.
# - If a newly reconstructed main52 limma table is unavailable, gene-level
#   core60 statistics are recovered from the executed checkpoint object.

source("scripts/_helpers.R")
if (!requireNamespace("edgeR", quietly = TRUE)) stop("Package 'edgeR' is required.")

in_dir <- p_out("04_GSE160269_epithelial")
out_dir <- p_out("05_core60")
ensure_dir(out_dir)

gsea_files <- c(
  main52 = file.path(in_dir, "GSE160269_Epithelial_main52_Hallmark_GSEA.rds"),
  main52_T4 = file.path(in_dir, "GSE160269_Epithelial_main52_Hallmark_GSEA_T4fixed.rds"),
  sensitivity55 = file.path(in_dir, "GSE160269_Epithelial_sensitivity55_Hallmark_GSEA.rds"),
  sensitivity55_T4 = file.path(in_dir, "GSE160269_Epithelial_sensitivity55_Hallmark_GSEA_T4fixed.rds")
)
for (f in gsea_files) assert_file(f)

target <- "HALLMARK_ESTROGEN_RESPONSE_EARLY"
get_leading <- function(path) {
  x <- readRDS(path)
  i <- which(x$pathway == target)
  if (length(i) != 1) stop("Expected one ", target, " row in ", path)
  as.character(x$leadingEdge[[i]])
}

leading <- lapply(gsea_files, get_leading)
leading_sizes <- vapply(leading, length, integer(1))
all_genes <- sort(unique(unlist(leading)))
recurrence <- vapply(all_genes, function(g) sum(vapply(leading, function(x) g %in% x, logical(1))), integer(1))
robustness <- data.frame(gene = all_genes, n_models = recurrence)
robustness <- robustness[order(-robustness$n_models, robustness$gene), ]
core60 <- robustness$gene[robustness$n_models == 4]

utils::write.csv(robustness, file.path(out_dir, "GSE160269_EstrogenEarly_leading_edge_robustness_reconstructed.csv"), row.names = FALSE)
utils::write.csv(data.frame(gene = core60), file.path(out_dir, "core60_genes_reconstructed.csv"), row.names = FALSE)
saveRDS(core60, file.path(out_dir, "core60_genes_reconstructed.rds"))

cat("Leading-edge sizes:", paste(leading_sizes, collapse = ", "), "\n")
cat("Union size:", length(all_genes), "core60 size:", length(core60), "\n")
stopifnot(length(all_genes) == 77L, length(core60) == 60L)
stopifnot(identical(as.integer(table(factor(recurrence, levels = 1:4))), c(7L, 4L, 6L, 60L)))

# Patient-level score: recovered old-chat code used the UNFILTERED 17,012-gene
# 52-patient pseudobulk matrix, created a fresh DGEList, applied TMM library
# normalization (without filterByExpr), then cpm(log=TRUE, prior.count=2).
counts_file <- p_data("05_scRNA_GSE160269", "GSE160269_Epithelial_main52_pseudobulk_counts.tsv.gz")
meta_file <- p_data("05_scRNA_GSE160269", "GSE160269_Epithelial_main52_analysis_metadata.csv")
assert_file(counts_file); assert_file(meta_file)

epi_pb <- read_pseudobulk_tsv(counts_file)
if (nrow(epi_pb) != 17012L) warning("Expected 17,012 unfiltered genes; observed ", nrow(epi_pb))
dge_epi_score <- edgeR::DGEList(counts = epi_pb)
if ("normLibSizes" %in% getNamespaceExports("edgeR")) {
  dge_epi_score <- edgeR::normLibSizes(dge_epi_score, method = "TMM")
} else {
  dge_epi_score <- edgeR::calcNormFactors(dge_epi_score, method = "TMM")
}
logcpm <- edgeR::cpm(dge_epi_score, log = TRUE, prior.count = 2)

meta <- utils::read.csv(meta_file, stringsAsFactors = FALSE, check.names = FALSE)
idx <- match(colnames(logcpm), meta$patient_id)
if (anyNA(idx)) stop("Unable to align main52 pseudobulk columns to metadata.")
meta <- meta[idx, , drop = FALSE]
meta$N_group <- factor(meta$N_group, levels = c("N0", "Nplus"))
meta$T_stage4 <- make_t_stage4(meta$T_stage)
meta$T_group <- make_t_group(meta$T_stage)

score_obj <- patient_mean_z_score(logcpm, core60)
if (length(score_obj$genes) != 60L) stop("Not all core60 genes were present in unfiltered main52 logCPM matrix.")
meta$core60_score <- as.numeric(score_obj$score)
meta$core60_score_z <- as.numeric(scale(meta$core60_score))

fit_main <- stats::lm(core60_score ~ T_group + N_group, data = meta)
fit_t4   <- stats::lm(core60_score ~ T_stage4 + N_group, data = meta)

saveRDS(fit_main, file.path(out_dir, "GSE160269_Epithelial_core60_score_lm_main_reconstructed.rds"))
saveRDS(fit_t4,   file.path(out_dir, "GSE160269_Epithelial_core60_score_lm_T4fixed_reconstructed.rds"))
utils::write.csv(meta, file.path(out_dir, "GSE160269_Epithelial_core60_patient_scores.csv"), row.names = FALSE)

coef_main <- extract_coef(fit_main)
utils::write.csv(coef_main, file.path(out_dir, "core60_score_main_association.csv"), row.names = FALSE)
print(coef_main)

# Gene-level statistics: prefer newly reconstructed limma output; otherwise use checkpoint.
limma_new <- file.path(in_dir, "GSE160269_Epithelial_main52_limma.rds")
if (file.exists(limma_new)) {
  lim <- readRDS(limma_new)$table
  core_stats <- lim[match(core60, lim$gene), , drop = FALSE]
} else {
  ck_file <- p_data("09_scripts", "checkpoint_core57_2026-10-08.RData")
  assert_file(ck_file)
  ck <- new.env()
  load(ck_file, envir = ck)
  if (!exists("core60_result", envir = ck, inherits = FALSE)) stop("Checkpoint lacks core60_result.")
  core_stats <- get("core60_result", envir = ck)
}
utils::write.csv(core_stats, file.path(out_dir, "core60_gene_statistics_reconstructed.csv"), row.names = FALSE)
