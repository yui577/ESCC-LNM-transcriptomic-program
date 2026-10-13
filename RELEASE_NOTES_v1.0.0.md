# v1.0.0 — submission release candidate

This is the planned first public software release accompanying the manuscript:

> *Cross-cohort transcriptomic analysis reveals heterogeneous bulk signatures and a recurrent epithelial transcriptional program associated with pathological nodal involvement in esophageal squamous cell carcinoma*

## Validation

- Local validation date: 2026-10-13
- Locked reproducibility checks: **28/28 PASS**
- Environment snapshot: `provenance/sessionInfo_final.txt`
- Validation table: `VALIDATION_RESULTS.tsv`
- Validation interpretation: `VALIDATION_STATUS.md`

## Included

- analysis scripts `01`–`09` and shared helpers;
- `run_all.R` orchestration script;
- data/input manifest and availability documentation;
- code provenance and recovery-status documentation;
- locked audit-reference tables used for numerical validation;
- retained late-session R history and final R session information;
- MIT License;
- completed citation metadata for the sole author, Jiafu Qian, including ORCID; Zenodo DOI/release date remain to be added after archiving.

## Scope

This release is an audited reproducibility reconstruction of analyses originally run interactively in RStudio. It is not presented as byte-for-byte recovery of every historical source script. The complete historical source for initial GEO download/reannotation, all TCGA clinical-filtering lines, the earliest GSE160269 aggregation pass, and exact final figure styling was not fully preserved.

The 28/28 PASS establishes numerical reproducibility of the locked manuscript-result chain from the preserved project inputs/intermediate objects used in the study. It does not imply a fully standalone raw-data-to-publication rebuild from GitHub alone.

## Before tagging the final public v1.0.0 release

1. Author ORCID added to `CITATION.cff`: https://orcid.org/0009-0005-7143-058X.
2. Add the final release date and Zenodo DOI to `CITATION.cff` after archiving the GitHub release.
3. Verify that the manuscript data/code availability statement uses the final GitHub URL and Zenodo DOI.
