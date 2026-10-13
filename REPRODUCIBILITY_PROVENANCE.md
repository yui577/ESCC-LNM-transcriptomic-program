# Reproducibility provenance and audit status

## Evidence recovered

The original analysis project directory contained a structured workflow with folders for bulk analysis, GSE160269, TCGA-ESCC, figures, tables, and scripts. Only one historical R command log was preserved as a source-like artifact: `09_scripts/R_history_2026-10-08.Rhistory`.

That history documents the downstream core57 bulk-projection work, including:

- use of the same 57 jointly measurable genes in GSE53624, GSE53622 and TCGA-ESCC;
- gene-wise standardization across patients;
- patient-level mean core score construction;
- within-cohort score standardization;
- adjustment for T stage, grade, age and sex;
- cohort effects used in the final analysis; and
- `metafor::rma(..., method = "REML")` random-effects meta-analysis.

A checkpoint file saved in the same session contains executed objects including `core60`, `core60_result`, four epithelial GSEA objects, three bulk analysis objects, and the core57 meta-analysis objects.

## Preserved analytical inputs / intermediate objects

The reconstructed scripts depend on saved objects from the executed analysis, including:

### Bulk GEO

- `03_bulk_analysis/GSE53624_gene_expression_119tumors.rds`
- `03_bulk_analysis/GSE53624_clinical_main.rds`
- `03_bulk_analysis/GSE53622_gene_expression_60tumors.rds`
- `03_bulk_analysis/GSE53622_clinical_main.rds`
- `03_bulk_analysis/GSE53624_GSE53622_tumor_normal_gene_concordance.rds`

### TCGA

- `06_TCGA_ESCC/TCGA_ESCC_voom_83patients.rds`

### GSE160269

- `05_scRNA_GSE160269/GSE160269_Epithelial_main52_pseudobulk_counts.tsv.gz`
- `05_scRNA_GSE160269/GSE160269_Epithelial_sensitivity55_pseudobulk_counts.tsv.gz`
- `05_scRNA_GSE160269/GSE160269_Epithelial_main52_edgeR_TMM.rds`
- `05_scRNA_GSE160269/GSE160269_Epithelial_main52_analysis_metadata.csv`
- `05_scRNA_GSE160269/GSE160269_Epithelial_main52_patient_metadata.csv`
- `05_scRNA_GSE160269/GSE160269_Epithelial_sensitivity55_patient_metadata.csv`
- `05_scRNA_GSE160269/GSE160269_Myeloid_60patients_pseudobulk_counts.rds`
- `05_scRNA_GSE160269/GSE160269_Tcell_60patients_pseudobulk_counts.rds`
- `05_scRNA_GSE160269/GSE160269_Bcell_56patients_pseudobulk_counts.rds`
- `05_scRNA_GSE160269/GSE160269_Fibroblast_49patients_pseudobulk_counts.rds`
- `05_scRNA_GSE160269/GSE160269_patient_clinical_Nstage.csv`
- original-study supplementary marker/clinical workbooks.

## Saved results used as audit references

The following executed outputs provide independent checks on reconstructed calculations:

- three-cohort Hallmark comparison;
- tumor-versus-normal positive-control summary;
- four epithelial Hallmark GSEA result objects;
- leading-edge recurrence table and core60 gene list;
- core60 patient-score model objects;
- residual model adjusted for Stress and Epi1/Epi2 set2;
- cross-cell-type GSEA objects;
- core57 cohort effects and REML meta-analysis output.

## Upstream steps not recoverable as original code

The following historical source scripts were not found and should not be represented as original executed code:

1. initial GEO download and GPL18109/GENCODE-v40 reannotation code;
2. initial TCGA/GDC download, ESCC restriction, clinical harmonization and raw-count preprocessing code;
3. initial GSE160269 cell-to-patient aggregation code;
4. exact final figure-layout/styling scripts.

The manuscript Methods describe these procedures, but the repository distinguishes that methodological description from recovered executable provenance.

## Validation and release status

The reconstructed workflow was executed in the original Windows/R environment and `09_reproducibility_checks.R` returned **28/28 PASS** on 2026-10-13. This provenance statement is retained in the public release, and the reconstructed scripts are not labeled as original historical source files.


## Additional recovery from the preserved interactive analysis transcript

The preserved execution transcript containing the GSE160269 interactive workflow was re-examined during repository reconstruction. It confirmed two details that were not preserved in the late-session `.Rhistory`:

- epithelial 52/55-patient limma fits used `eBayes(..., robust = TRUE)`;
- GSE160269 Hallmark analyses used the project-local `05_scRNA_GSE160269/h.all.v2025.1.Hs.symbols.gmt`;
- cross-cell-type Myeloid/T/B/Fibroblast fits used ordinary `eBayes(fit)` (not robust=TRUE), while using the same local Hallmark GMT.

These recovered details are encoded in scripts 04 and 07.


## Additional command-level recovery

The preserved interactive analysis transcript was sufficient to restore command-level details for bulk GEO final limma fits, TCGA-ESCC limma, bulk Hallmark GSEA, GSE160269 epithelial 52/55-patient limma/GSEA, core60 leading-edge recurrence, and the myeloid cross-cell-type pattern. Those recovered commands are encoded in scripts 01–05 and 07 and summarized in `CODE_PROVENANCE.tsv`.  This improves provenance but does not convert the curated repository files into byte-for-byte original historical source scripts.
