# Data availability and local inputs

The study uses public transcriptomic datasets available from the NCBI Gene Expression Omnibus under accession numbers **GSE53624**, **GSE53622**, and **GSE160269**, and from the National Cancer Institute Genomic Data Commons under project **TCGA-ESCA**, restricted in the analysis to histologically confirmed ESCC.

Large public raw data files are intentionally not duplicated in this repository. The validated reconstruction starts from preserved intermediate objects for analytical stages whose complete historical upstream source code was not retained. `DATA_MANIFEST.tsv` lists the expected local inputs used for the validated manuscript-result chain.

For local verification against a project directory containing those inputs, set for example:

```r
Sys.setenv(
  ESCC_DATA_ROOT = "<path-to-local-ESCC-project>",
  ESCC_OUT_ROOT  = "<path-to-output-directory>"
)
```

The 28/28 validation result therefore establishes numerical reproducibility of the locked manuscript-result chain from the preserved inputs/intermediate objects used in the study; it is not represented as a complete raw-data-to-publication rebuild from GitHub alone.
