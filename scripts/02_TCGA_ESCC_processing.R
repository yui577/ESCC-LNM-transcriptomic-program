# 02_TCGA_ESCC_processing.R
# Reconstructs the reported TCGA-ESCC pN limma analysis from the saved voom EList.
#
# PROVENANCE STATUS:
# - The final TCGA edgeR/voom/limma commands were recovered from the preserved execution transcript.
# - GDC acquisition, histology restriction, clinical matching, filterByExpr,
#   TMM normalization, and voom transformation were executed previously and
#   preserved in TCGA_ESCC_voom_83patients.rds.
# - The original upstream script was not preserved.
# - This script re-fits the locked final design to the saved voom object.

source("scripts/_helpers.R")
if (!requireNamespace("limma", quietly = TRUE)) stop("Package 'limma' is required.")

infile <- p_data("06_TCGA_ESCC", "TCGA_ESCC_voom_83patients.rds")
assert_file(infile)
out_dir <- p_out("02_TCGA_ESCC")
ensure_dir(out_dir)

v <- readRDS(infile)
if (is.null(v$E) || is.null(v$design)) stop("Saved TCGA object must contain E and design.")

expected_terms <- c("N_groupNplus", "T_stageT2", "T_stageT3", "T_stageT4",
                    "grademoderately", "gradepoorly", "age_years", "sexmale")
missing_terms <- setdiff(expected_terms, colnames(v$design))
if (length(missing_terms)) stop("TCGA design missing terms: ", paste(missing_terms, collapse = ", "))

# lmFit on the EList retains saved voom weights.
fit <- limma::lmFit(v, design = v$design)
fit <- limma::eBayes(fit, robust = TRUE)
tt <- limma::topTable(fit, coef = "N_groupNplus", number = Inf, adjust.method = "BH", sort.by = "P")
tt$gene <- rownames(tt)

saveRDS(list(voom = v, fit = fit, table = tt),
        file.path(out_dir, "TCGA_ESCC_limma_reconstructed.rds"))
utils::write.table(tt,
                   file.path(out_dir, "TCGA_ESCC_limma_reconstructed.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)

cat("TCGA-ESCC samples:", ncol(v$E), "genes:", nrow(v$E), "\n")
cat("Expected final dimensions: 83 patients and 16,825 genes.\n")
