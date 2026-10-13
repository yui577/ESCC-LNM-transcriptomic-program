# Recovery status by analysis module

| Module | Status | What is actually supported |
|---|---|---|
| GEO adjusted pN limma | High | Final design, robust eBayes and topTable commands recovered; saved gene-level matrices/clinical tables available. |
| GPL18109 -> GENCODE40 | Moderate | Exact-full-length/unique mapping logic, 32,200 probes, 17,745 genes and arithmetic-mean collapse verified; original full script not preserved. |
| Bulk Hallmark GSEA | High | Ranking, msigdbr Hallmark, seed=123 and fgseaMultilevel parameters recovered. |
| GEO tumor-vs-normal positive control | High for result, moderate for source | Saved executed gene-level concordance object permits exact rho/direction recomputation. |
| TCGA adjusted edgeR/voom/limma | High | DGEList/filterByExpr/TMM/voom/robust limma commands recovered; saved 83-patient voom object available. |
| TCGA GDC/clinical restriction | Moderate | GDC STAR-count acquisition recovered; exact clinical complete-case filtering source not fully retained. |
| GSE160269 epithelial 52/55 | High | Command-level limma/GSEA details recovered, including robust eBayes, local GMT, one filtering pass per cohort and RNG behavior. |
| core60 recurrence | High | Actual recurrence/intersection code pattern recovered and executed outputs preserved. |
| core60 patient score | High | preserved execution transcript supports unfiltered 17,012-gene pseudobulk -> fresh TMM/logCPM -> gene-wise z -> mean score and T-adjusted models. |
| Published epithelial programs | Moderate-high | Exact marker parsing, unfiltered-pseudobulk TMM/logCPM score basis and residual model recovered; association-table wrapper is curated and checked against locked outputs. |
| Myeloid cross-cell | High | Main/T4 limma/GSEA pattern directly recovered. |
| T/B/Fibroblast cross-cell | Moderate-high | Parallel pipeline supported by saved objects and Methods; full historical command text not all retained. |
| core57/meta | Very high | Direct late-session R history plus checkpoint and saved meta outputs. |
| ORA | High for parameters/result | gprofiler2 parameters and final output locked; audit script reruns from core60/background. |
| Figure 1-5 styling/source | Low | Final files/object names exist; complete historical plotting source was not preserved and is not fabricated. |
