# Analysis scripts

This directory contains the **reconstructed reproducibility workflow** for the ESCC nodal-involvement manuscript.

The original `.R` source files were not preserved. These scripts are therefore not presented as byte-for-byte historical source code. They were reconstructed from preserved analysis objects, the saved R command history for the downstream core57 analysis, the finalized Methods, and saved numerical outputs.

| Script | Component | Provenance |
|---|---|---|
| `01_bulk_GEO_processing.R` | Re-fit GSE53624/GSE53622 pN limma models from preserved gene-level expression/clinical objects | reconstructed from saved executed inputs + Methods |
| `02_TCGA_ESCC_processing.R` | Re-fit TCGA-ESCC model from preserved voom EList/design | reconstructed from saved executed input + Methods |
| `03_bulk_GSEA.R` | Hallmark GSEA, three-cohort comparison, four-pathway prioritization, GEO positive control | reconstructed + checked against saved outputs |
| `04_GSE160269_pseudobulk.R` | Four epithelial pseudobulk sensitivity analyses | reconstructed from preserved pseudobulk objects/counts + Methods |
| `05_core60_analysis.R` | Leading-edge recurrence, core60 definition, patient score | reconstructed from saved GSEA/DGEList/checkpoint objects |
| `06_published_program_comparison.R` | Published marker overlap, correlations, corrected eight-program BH analysis, residual model | reconstructed from source marker table + saved expression/model objects |
| `07_cross_celltype_analysis.R` | Myeloid/T/B/fibroblast Hallmark and core60 analyses | reconstructed from preserved pseudobulk counts + Methods |
| `08_bulk_projection_meta.R` | core57 projection and REML random-effects meta-analysis | strongest provenance; largely recovered from executed R history |
| `09_reproducibility_checks.R` | Locked numerical audit checks and optional core60 ORA | audit script |

Run `../run_all.R` from the repository root after setting `ESCC_DATA_ROOT` and `ESCC_OUT_ROOT` as described in the repository README.

See `../REPRODUCIBILITY_PROVENANCE.md` for the exact recovery limitations and release policy.
