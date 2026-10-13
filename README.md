# ESCC nodal-involvement transcriptomic analysis — audited reproducibility reconstruction

Reproducibility repository for the manuscript:

> *Cross-cohort transcriptomic analysis reveals heterogeneous bulk signatures and a recurrent epithelial transcriptional program associated with pathological nodal involvement in esophageal squamous cell carcinoma*

## What this repository is

The analysis was originally performed interactively in RStudio. Most standalone historical `.R` source files were not saved. This repository is therefore a **curated and audited reproducibility reconstruction**, not a claim of byte-for-byte recovery of the original scripts.

The reconstruction uses four evidence layers:

1. executed R history (`provenance/R_history_2026-10-08.Rhistory`);
2. saved intermediate/model/result objects from the executed analysis;
3. command blocks preserved in an interactive analysis transcript from the original RStudio workflow; and
4. locked manuscript/supplementary numerical outputs used as audit targets.

Stage-specific evidence is documented in `CODE_PROVENANCE.tsv`, `RECOVERY_STATUS.md`, and `REPRODUCIBILITY_PROVENANCE.md`. Data inputs are summarized in `DATA_MANIFEST.tsv` and `DATA_AVAILABILITY.md`. The completed local validation is reported in `VALIDATION_STATUS.md` and `VALIDATION_RESULTS.tsv`.

## Validation status

**PASS: 28/28 locked reproducibility checks passed on 2026-10-13** against the preserved ESCC project inputs/intermediate objects used in the manuscript analysis. See `VALIDATION_STATUS.md` and `VALIDATION_RESULTS.tsv`.

This validation establishes numerical reproducibility of the locked manuscript-result chain from those preserved inputs. It does **not** claim that every public raw dataset can be regenerated end-to-end using only the files committed to this repository.

## Analysis modules

| Script | Main task |
|---|---|
| `01_bulk_GEO_processing.R` | Re-fit GSE53624/GSE53622 adjusted pN limma models from saved gene-level inputs |
| `02_TCGA_ESCC_processing.R` | Re-fit the 83-patient TCGA-ESCC adjusted limma model from the saved voom object |
| `03_bulk_GSEA.R` | Three-cohort Hallmark GSEA, four-pathway prioritization, tumor-vs-normal positive control |
| `04_GSE160269_pseudobulk.R` | Four epithelial 52/55-patient pseudobulk sensitivity analyses using recovered model/GSEA details |
| `05_core60_analysis.R` | ESTROGEN_RESPONSE_EARLY leading-edge recurrence, core60, patient score and gene statistics |
| `06_published_program_comparison.R` | Published epithelial marker-program overlap, correlations, pN associations and residual model |
| `07_cross_celltype_analysis.R` | Myeloid/T/B/fibroblast Hallmark analysis and epithelial-core projection |
| `08_bulk_projection_meta.R` | core57 projection into three bulk cohorts and REML random-effects meta-analysis |
| `09_reproducibility_checks.R` | Locked-result checks and core60 g:Profiler ORA |

`run_all.R` executes the modules sequentially.

## Important recovered execution details

Several details were recovered from the earlier execution conversation and materially improve reproducibility over a Methods-only reconstruction:

- epithelial 52/55 limma fits used `eBayes(..., robust = TRUE)`;
- GSE160269 Hallmark analyses used the project-local `h.all.v2025.1.Hs.symbols.gmt`;
- each epithelial cohort was filtered/TMM-normalized once, then the same DGEList was reused for binary-T and four-level-T models;
- the original fgsea RNG behavior is preserved in script 04;
- myeloid/T/B/fibroblast fits used ordinary `eBayes(fit)` rather than robust eBayes;
- published marker programs were read from MOESM8 sheet `Genes`, rows 4:33, eight columns;
- the core57 standardization/regression/meta-analysis block is directly represented in the retained R history.

## Expected locked results

The final audit is expected to recover/check, among other values:

- bulk candidate Hallmarks: `OXIDATIVE_PHOSPHORYLATION`, `ADIPOGENESIS`, `GLYCOLYSIS`, `ESTROGEN_RESPONSE_EARLY`;
- tumor-vs-normal positive control: Spearman rho about 0.922 and 87.3% directional concordance;
- epithelial main cohort: 52 patients (28 pN0, 24 pN+), 44,122 cells;
- ESTROGEN_RESPONSE_EARLY leading-edge sizes 63, 72, 73 and 65; union 77; recurrence counts 7/4/6/60; core60 size 60;
- epithelial core60 pN coefficient about 0.58 (95% CI 0.35-0.81, P about 7.8e-6);
- myeloid ESTROGEN_RESPONSE_EARLY NES/FDR about 1.51/0.023 and 1.58/0.028;
- core57 cohort coefficients about 0.007, 0.277 and 0.653; REML pooled beta about 0.287, P about 0.148, I2 about 56.6%;
- corrected eight-program BH FDR: Epi1/Epi2 set2 0.004585 and Stress 0.007954.

## Running locally

The scripts do not overwrite the original project outputs. Point them to the original project directory and a separate output directory:

```r
setwd("<path-to-this-repository>")

Sys.setenv(
  ESCC_DATA_ROOT = "<path-to-local-ESCC-project>",
  ESCC_OUT_ROOT  = "<path-to-output-directory>"
)

source("run_all.R")
```

The repository intentionally relies on preserved intermediate objects for stages whose original upstream source was not fully recoverable. These large local intermediates are not committed to GitHub. Public source datasets should be obtained from GEO/GDC using the accession identifiers listed in `DATA_AVAILABILITY.md`.

## Software environment

The original final environment snapshot is included at `provenance/sessionInfo_final.txt`. Key versions were:

- R 4.6.1 on Windows 11
- limma 3.68.5
- edgeR 4.10.5
- fgsea 1.38.0
- msigdbr 26.1.1
- metafor 5.2.1
- gprofiler2 0.2.4

## What is deliberately not claimed

The repository does **not** label reconstructed files as original historical scripts. In particular, the complete historical source for initial GEO download/reannotation, every TCGA clinical-filtering line, the earliest cell-to-patient aggregation pass, and final Figure 1-5 styling was not preserved. Where only partial execution evidence exists, that limitation is stated explicitly.

Exact historical Figure 1-5 plotting source was not consistently recoverable. Final figure files and result tables exist, but no invented plotting code is presented as the original source.


## Citation and license

This release candidate is prepared for a `v1.0.0` public release. Code is licensed under the MIT License; see `LICENSE`.

Citation metadata are provided in `CITATION.cff` for the sole author, Jiafu Qian. The author's ORCID is included in `CITATION.cff`. The Zenodo DOI and release date will be added after the public `v1.0.0` release is archived.

Planned canonical repository URL: `https://github.com/yui577/ESCC-LNM-transcriptomic-program`.
