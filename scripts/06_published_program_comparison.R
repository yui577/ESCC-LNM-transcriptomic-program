# 06_published_program_comparison.R
# Reconstructs comparison of core60 with eight published GSE160269 epithelial
# marker-based programs (overlap, correlations, pN associations, residual model).
#
# PROVENANCE STATUS:
# - Marker-table parsing (MOESM8, sheet=Genes, rows 4:33, eight columns), the
#   residual-model formula, and key object names were recovered from the preserved execution transcript.
# - Patient-level marker-score construction and the eight-program BH analysis are
#   curated from the recovered workflow and saved outputs; this is not claimed to
#   be a byte-for-byte historical source file.

source("scripts/_helpers.R")
if (!requireNamespace("edgeR", quietly = TRUE)) stop("Package 'edgeR' is required.")

out_dir <- p_out("06_published_programs")
ensure_dir(out_dir)

core60 <- readRDS(p_out("05_core60", "core60_genes_reconstructed.rds"))
score_meta <- utils::read.csv(p_out("05_core60", "GSE160269_Epithelial_core60_patient_scores.csv"),
                              stringsAsFactors = FALSE, check.names = FALSE)
score_meta$N_group <- factor(score_meta$N_group, levels = c("N0", "Nplus"))
score_meta$T_stage4 <- make_t_stage4(score_meta$T_stage)

# Recovered old-chat scoring workflow: use the UNFILTERED 17,012-gene main52
# pseudobulk count matrix, TMM-normalize without filterByExpr, and transform to
# logCPM with prior.count=2. This is the same expression basis used for core60
# and the reconstructed published marker-based patient scores.
epi_pb <- read_pseudobulk_tsv(
  p_data("05_scRNA_GSE160269", "GSE160269_Epithelial_main52_pseudobulk_counts.tsv.gz")
)
dge_epi_score <- edgeR::DGEList(counts = epi_pb)
if ("normLibSizes" %in% getNamespaceExports("edgeR")) {
  dge_epi_score <- edgeR::normLibSizes(dge_epi_score, method = "TMM")
} else {
  dge_epi_score <- edgeR::calcNormFactors(dge_epi_score, method = "TMM")
}
logcpm <- edgeR::cpm(dge_epi_score, log = TRUE, prior.count = 2)
idx <- match(colnames(logcpm), score_meta$patient_id)
if (anyNA(idx)) stop("Unable to align main52 logCPM to score metadata.")
score_meta <- score_meta[idx, , drop = FALSE]

marker_file <- p_data("05_scRNA_GSE160269", "41467_2021_25539_MOESM8_ESM.xlsx")
assert_file(marker_file)

# Recovered old-chat code read the worksheet without headers, then selected
# rows 4:33 (the 30 marker genes) and the first eight program columns.
if (requireNamespace("readxl", quietly = TRUE)) {
  prog_raw <- readxl::read_excel(marker_file, sheet = "Genes", col_names = FALSE)
  program_mat <- as.data.frame(prog_raw[4:33, 1:8])
} else if (requireNamespace("openxlsx", quietly = TRUE)) {
  prog_raw <- openxlsx::read.xlsx(marker_file, sheet = "Genes", colNames = FALSE)
  program_mat <- as.data.frame(prog_raw[4:33, 1:8])
} else {
  stop("Package 'readxl' or 'openxlsx' is required to read published marker spreadsheet.")
}

program_names <- c("Mucosal", "Stress", "AP", "Cycling",
                   "Epi1Epi2_set1", "Epi1Epi2_set2", "Mes", "Oxd")
names(program_mat) <- program_names
# Recovered pattern: programs8 <- lapply(program_mat, function(x)
#   unique(trimws(as.character(x)[!is.na(x) & nzchar(x)])))
markers <- lapply(program_mat, function(x) {
  x <- as.character(x)
  unique(trimws(x[!is.na(x) & nzchar(x)]))
})

# Direct overlap.
overlap_tab <- do.call(rbind, lapply(program_names, function(nm) {
  ov <- intersect(core60, markers[[nm]])
  data.frame(program = nm,
             n_markers = length(markers[[nm]]),
             n_overlap = length(ov),
             overlap_genes = paste(ov, collapse = ";"),
             stringsAsFactors = FALSE)
}))
utils::write.csv(overlap_tab, file.path(out_dir, "core60_published_program_overlap.csv"), row.names = FALSE)
unique_overlap <- unique(unlist(strsplit(paste(overlap_tab$overlap_genes, collapse = ";"), ";", fixed = TRUE)))
unique_overlap <- unique_overlap[nzchar(unique_overlap)]
cat("Unique core60 genes overlapping any published program:", length(unique_overlap), "\n")

# Reconstructed patient-level marker-based program scores.
for (nm in program_names) {
  use <- intersect(markers[[nm]], rownames(logcpm))
  if (length(use) == 0) stop("No measurable markers for program: ", nm)
  z <- z_by_gene(logcpm[use, , drop = FALSE])
  raw_score <- colMeans(z, na.rm = TRUE)
  score_meta[[paste0(nm, "_score")]] <- as.numeric(raw_score)
  score_meta[[paste0(nm, "_z")]] <- as.numeric(scale(raw_score))
}
# Preserve the historical residual-model term name recorded in the saved model
# object and locked Methods while retaining the explicit S9 program label above.
score_meta$Epi_set2_z <- score_meta$Epi1Epi2_set2_z
utils::write.csv(score_meta, file.path(out_dir, "published_program_patient_scores.csv"), row.names = FALSE)

# Spearman correlations with core60 score.
cor_tab <- do.call(rbind, lapply(program_names, function(nm) {
  ct <- suppressWarnings(stats::cor.test(score_meta$core60_score,
                                         score_meta[[paste0(nm, "_score")]],
                                         method = "spearman", exact = FALSE))
  data.frame(program = nm,
             spearman_rho = unname(ct$estimate),
             p_value = ct$p.value)
}))
cor_tab$FDR <- p.adjust(cor_tab$p_value, method = "BH")
utils::write.csv(cor_tab, file.path(out_dir, "core60_program_correlations.csv"), row.names = FALSE)

# Four-level T-stage-adjusted pN associations of standardized program scores.
assoc <- do.call(rbind, lapply(program_names, function(nm) {
  zcol <- paste0(nm, "_z")
  fit <- stats::lm(stats::reformulate(c("T_stage4", "N_group"), response = zcol), data = score_meta)
  a <- extract_coef(fit)
  data.frame(Program = nm,
             Beta_SD = a$beta,
             SE = a$se,
             P_value = a$p_value,
             CI_low = a$ci_low,
             CI_high = a$ci_high,
             stringsAsFactors = FALSE)
}))
# IMPORTANT: final corrected S9 applies BH across the eight published programs only.
assoc$FDR <- p.adjust(assoc$P_value, method = "BH")
assoc <- assoc[order(assoc$FDR, assoc$P_value), ]
utils::write.csv(assoc, file.path(out_dir, "Published_program_pN_associations_8programs_BH.csv"), row.names = FALSE)

# Baseline standardized core60 model and exploratory residual model.
fit_core_base <- stats::lm(core60_score_z ~ T_stage4 + N_group, data = score_meta)
fit_core_resid <- stats::lm(core60_score_z ~ T_stage4 + N_group + Stress_z + Epi_set2_z,
                            data = score_meta)
saveRDS(fit_core_base, file.path(out_dir, "core60_T4_standardized_model.rds"))
saveRDS(fit_core_resid, file.path(out_dir, "core60_adjusted_Stress_Epi_set2_model.rds"))

resid_tab <- rbind(
  data.frame(model = "T_stage4 + N_group", extract_coef(fit_core_base)),
  data.frame(model = "T_stage4 + N_group + Stress + Epi_set2", extract_coef(fit_core_resid))
)
utils::write.csv(resid_tab, file.path(out_dir, "core60_adjusted_program_model.csv"), row.names = FALSE)

print(cor_tab[order(-cor_tab$spearman_rho), ])
print(assoc)
print(resid_tab)
