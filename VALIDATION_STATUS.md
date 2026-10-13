# Validation status

**Status: PASS — 28/28 locked reproducibility checks passed on 2026-10-13.**

The reconstructed analysis was executed locally against the original ESCC project data root on Windows/R and the final reproducibility checklist returned `TRUE` for all 28 locked checks.

Key exact/near-exact validations included:

- 4 prioritized bulk Hallmarks recovered exactly.
- Tumor-vs-normal positive-control Spearman rho = 0.921761; directional concordance = 0.873147.
- GSE160269 epithelial cohort: 52 patients (28 pN0, 24 pN+), 44,122 epithelial cells.
- ESTROGEN_RESPONSE_EARLY leading-edge union = 77; recurrence distribution = 7/4/6/60; core60 = 60.
- core60 main pN coefficient = 0.580286; P = 7.792559e-06.
- Residual model pN coefficient = 0.519604; P = 0.017248.
- Myeloid ESTROGEN_RESPONSE_EARLY: saved and reconstructed results matched exactly for both models:
  - binary-T: NES = 1.508108; FDR = 0.02298321.
  - T1-T4: NES = 1.580379; FDR = 0.02810000.
- core57 bulk coefficients = 0.007336, 0.277291, 0.652863.
- REML pooled core57 beta = 0.2866; P = 0.1476; I2 = 56.58%.
- Corrected 8-program BH values matched the final manuscript values.

This validation establishes numerical reproducibility of the locked manuscript-result chain from the preserved project inputs/intermediate objects used by this repository. It does not convert partially recovered historical source into original source; provenance limitations remain documented in `CODE_PROVENANCE.tsv` and `RECOVERY_STATUS.md`.
