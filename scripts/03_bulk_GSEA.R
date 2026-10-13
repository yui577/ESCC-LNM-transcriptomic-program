# 03_bulk_GSEA.R
# Hallmark GSEA across GSE53624, GSE53622 and TCGA-ESCC, plus the
# tumor-versus-normal GEO positive-control concordance calculation.
#
# PROVENANCE STATUS:
# The bulk Hallmark ranking/fgsea commands (msigdbr Hallmark, moderated-t rank,
# set.seed(123), fgseaMultilevel minSize=15/maxSize=500) were recovered from the
# preserved execution transcript. The positive-control calculation is an audited reconstruction from
# the saved executed gene-level concordance object.

source("scripts/_helpers.R")

geo24_file <- p_out("01_bulk_GEO", "GSE53624_limma_reconstructed.rds")
geo22_file <- p_out("01_bulk_GEO", "GSE53622_limma_reconstructed.rds")
tcga_file  <- p_out("02_TCGA_ESCC", "TCGA_ESCC_limma_reconstructed.rds")
assert_file(geo24_file); assert_file(geo22_file); assert_file(tcga_file)

geo24 <- readRDS(geo24_file)
geo22 <- readRDS(geo22_file)
tcga  <- readRDS(tcga_file)

hallmarks <- hallmark_sets()

stats_from_tt <- function(tt) {
  stat <- tt$t
  names(stat) <- if ("gene" %in% names(tt)) tt$gene else rownames(tt)
  stat
}

gsea24 <- run_fgsea(stats_from_tt(geo24$table), hallmarks, min_size = 15, max_size = 500, seed = 123)
gsea22 <- run_fgsea(stats_from_tt(geo22$table), hallmarks, min_size = 15, max_size = 500, seed = 123)
gseaTC <- run_fgsea(stats_from_tt(tcga$table),  hallmarks, min_size = 15, max_size = 500, seed = 123)

out_dir <- p_out("03_bulk_GSEA")
ensure_dir(out_dir)
saveRDS(gsea24, file.path(out_dir, "GSE53624_Hallmark_GSEA_reconstructed.rds"))
saveRDS(gsea22, file.path(out_dir, "GSE53622_Hallmark_GSEA_reconstructed.rds"))
saveRDS(gseaTC, file.path(out_dir, "TCGA_ESCC_Hallmark_GSEA_reconstructed.rds"))

pick <- function(x, suffix) {
  data.frame(
    pathway = x$pathway,
    NES = x$NES,
    P = x$pval,
    FDR = x$padj,
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) |>
    stats::setNames(c("pathway", paste0("NES_", suffix), paste0("P_", suffix), paste0("FDR_", suffix)))
}

cmp <- Reduce(function(x, y) merge(x, y, by = "pathway", all = TRUE),
              list(pick(gsea24, "53624"), pick(gsea22, "53622"), pick(gseaTC, "TCGA")))

cmp$same_direction_3 <- with(cmp,
  (NES_53624 > 0 & NES_53622 > 0 & NES_TCGA > 0) |
  (NES_53624 < 0 & NES_53622 < 0 & NES_TCGA < 0))
cmp$significant_53624 <- cmp$FDR_53624 < 0.05
cmp$significant_53622 <- cmp$FDR_53622 < 0.05
cmp$significant_TCGA  <- cmp$FDR_TCGA  < 0.05
cmp$n_FDR05 <- rowSums(cbind(cmp$significant_53624, cmp$significant_53622, cmp$significant_TCGA), na.rm = TRUE)

cmp <- cmp[order(-cmp$n_FDR05, pmin(cmp$FDR_53624, cmp$FDR_53622, cmp$FDR_TCGA, na.rm = TRUE)), ]
utils::write.table(cmp, file.path(out_dir, "Hallmark_three_cohort_comparison_reconstructed.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)

candidate4 <- subset(cmp,
                     NES_53624 > 0 & NES_53622 > 0 & NES_TCGA > 0 & n_FDR05 >= 2)
utils::write.csv(candidate4, file.path(out_dir, "bulk_prioritized_candidate4.csv"), row.names = FALSE)

expected_candidate4 <- sort(c(
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_ADIPOGENESIS",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_ESTROGEN_RESPONSE_EARLY"
))
stopifnot(identical(sort(candidate4$pathway), expected_candidate4))

# Positive control: saved gene-level tumor-vs-normal moderated t statistics.
concord_file <- p_data("03_bulk_analysis", "GSE53624_GSE53622_tumor_normal_gene_concordance.rds")
assert_file(concord_file)
concord <- readRDS(concord_file)
required_cols <- c("gene_symbol", "t_53624", "t_53622")
if (!all(required_cols %in% names(concord))) stop("Unexpected positive-control object structure.")

rho <- stats::cor(concord$t_53624, concord$t_53622, method = "spearman", use = "complete.obs")
dir_concordance <- mean(sign(concord$t_53624) == sign(concord$t_53622), na.rm = TRUE)
pc <- data.frame(metric = c("Spearman_rho", "Direction_concordance"),
                 value = c(rho, dir_concordance))
utils::write.csv(pc, file.path(out_dir, "tumor_normal_positive_control_reconstructed.csv"), row.names = FALSE)

cat("Prioritized pathways (must be 4):\n")
print(candidate4[, c("pathway", "NES_53624", "FDR_53624", "NES_53622", "FDR_53622", "NES_TCGA", "FDR_TCGA")])
cat(sprintf("Positive control: rho=%.6f; directional concordance=%.6f\n", rho, dir_concordance))
