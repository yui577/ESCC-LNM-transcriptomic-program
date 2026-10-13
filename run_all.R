# Run from the repository root.
# Set ESCC_DATA_ROOT to the original analysis project directory, e.g. on Windows:
# Sys.setenv(ESCC_DATA_ROOT = "<path-to-local-ESCC-project>")
# Optionally set a separate output location:
# Sys.setenv(ESCC_OUT_ROOT = "<path-to-output-directory>")

scripts <- c(
  "01_bulk_GEO_processing.R",
  "02_TCGA_ESCC_processing.R",
  "03_bulk_GSEA.R",
  "04_GSE160269_pseudobulk.R",
  "05_core60_analysis.R",
  "06_published_program_comparison.R",
  "07_cross_celltype_analysis.R",
  "08_bulk_projection_meta.R",
  "09_reproducibility_checks.R"
)
for (s in scripts) {
  message("\n===== Running ", s, " =====")
  source(file.path("scripts", s), echo = FALSE)
}
