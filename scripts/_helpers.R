# Helper functions for the reconstructed ESCC analysis scripts.
# STATUS: reconstructed for reproducibility from saved analysis objects,
# manuscript methods, and executed R history. This is not an original source file.

get_data_root <- function() {
  root <- Sys.getenv("ESCC_DATA_ROOT", unset = ".")
  normalizePath(root, winslash = "/", mustWork = TRUE)
}

get_out_root <- function() {
  out <- Sys.getenv("ESCC_OUT_ROOT", unset = "reconstructed_outputs")
  if (!dir.exists(out)) dir.create(out, recursive = TRUE)
  normalizePath(out, winslash = "/", mustWork = TRUE)
}

p_data <- function(...) file.path(get_data_root(), ...)
p_out  <- function(...) file.path(get_out_root(), ...)

assert_file <- function(path) {
  if (!file.exists(path)) stop("Required file not found: ", path)
  invisible(path)
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  invisible(path)
}

make_n_group <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  out <- ifelse(x %in% c(1, "1", "N1", "N2", "N3", "Nplus", "pN+"), "Nplus", "N0")
  factor(out, levels = c("N0", "Nplus"))
}

make_t_stage4 <- function(x) {
  x <- as.character(x)
  x <- sub("^pT", "T", x)
  x <- ifelse(grepl("^[1-4]$", x), paste0("T", x), x)
  factor(x, levels = c("T1", "T2", "T3", "T4"))
}

make_t_group <- function(x) {
  t4 <- make_t_stage4(x)
  factor(ifelse(t4 %in% c("T1", "T2"), "Early_T", "Advanced_T"),
         levels = c("Early_T", "Advanced_T"))
}


# Read saved GEO gene-level expression objects. The executed history shows that
# these saved tables contain two annotation columns (gene_id, gene_symbol)
# followed by patient expression columns. Strip annotation columns before
# numeric conversion and use gene_symbol as row names.
read_geo_gene_matrix <- function(path) {
  assert_file(path)
  x <- readRDS(path)

  if (is.data.frame(x)) {
    nms <- names(x)
    ann <- intersect(c("gene_id", "gene_symbol"), nms)

    if ("gene_symbol" %in% nms) {
      genes <- as.character(x$gene_symbol)
    } else {
      genes <- rownames(x)
    }

    sample_cols <- setdiff(nms, ann)
    if (!length(sample_cols)) stop("No expression sample columns found in: ", path)
    m <- as.matrix(x[, sample_cols, drop = FALSE])

    if (!is.null(genes) && length(genes) == nrow(m) &&
        all(!is.na(genes)) && all(nzchar(genes))) {
      if (anyDuplicated(genes)) {
        stop("Duplicated gene symbols found in saved gene-level object: ", path)
      }
      rownames(m) <- genes
    }
  } else {
    m <- as.matrix(x)
    # Defensive handling in case annotation columns survived in a matrix-like object.
    ann <- intersect(c("gene_id", "gene_symbol"), colnames(m))
    if (length(ann)) {
      genes <- if ("gene_symbol" %in% colnames(m)) as.character(m[, "gene_symbol"]) else rownames(m)
      m <- m[, setdiff(colnames(m), ann), drop = FALSE]
      if (!is.null(genes) && length(genes) == nrow(m)) rownames(m) <- genes
    }
  }

  suppressWarnings(storage.mode(m) <- "numeric")
  if (anyNA(m)) {
    stop("NA values appeared during numeric conversion of GEO expression matrix: ", path)
  }
  m
}

align_clinical <- function(expr, clin, id_col = "patient_id") {
  if (is.null(colnames(expr))) stop("Expression matrix has no sample column names.")
  if (!id_col %in% names(clin)) stop("Clinical table lacks column: ", id_col)
  idx <- match(colnames(expr), clin[[id_col]])
  if (anyNA(idx)) {
    stop("Clinical matching failed for: ", paste(colnames(expr)[is.na(idx)], collapse = ", "))
  }
  clin[idx, , drop = FALSE]
}

hallmark_sets <- function() {
  if (!requireNamespace("msigdbr", quietly = TRUE)) stop("Package 'msigdbr' is required.")
  msig <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
  split(msig$gene_symbol, msig$gs_name)
}

hallmark_sets_gse160269 <- function() {
  if (!requireNamespace("fgsea", quietly = TRUE)) stop("Package 'fgsea' is required.")
  gmt <- p_data("05_scRNA_GSE160269", "h.all.v2025.1.Hs.symbols.gmt")
  assert_file(gmt)
  fgsea::gmtPathways(gmt)
}

run_fgsea <- function(stats, pathways, min_size = 15, max_size = 500, eps = NULL, seed = 123) {
  if (!requireNamespace("fgsea", quietly = TRUE)) stop("Package 'fgsea' is required.")
  stats <- stats[is.finite(stats)]
  stats <- sort(stats, decreasing = TRUE)
  if (!is.null(seed)) set.seed(seed)
  args <- list(pathways = pathways, stats = stats, minSize = min_size, maxSize = max_size)
  if (!is.null(eps)) args$eps <- eps
  res <- do.call(fgsea::fgseaMultilevel, args)
  res[order(res$padj, -abs(res$NES)), ]
}

read_pseudobulk_tsv <- function(path) {
  assert_file(path)
  x <- utils::read.delim(path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  m <- as.matrix(x)
  storage.mode(m) <- "numeric"
  m
}

align_pseudobulk_meta <- function(counts, meta) {
  if (is.null(colnames(counts))) stop("Pseudobulk count matrix has no patient column names.")
  idx <- match(colnames(counts), meta$patient_id)
  if (anyNA(idx)) stop("Unable to match pseudobulk columns to metadata.")
  meta <- meta[idx, , drop = FALSE]
  stopifnot(identical(as.character(meta$patient_id), colnames(counts)))
  meta
}

prepare_pseudobulk_dge <- function(counts, meta, filter_formula) {
  if (!requireNamespace("edgeR", quietly = TRUE)) stop("Package 'edgeR' is required.")
  meta <- align_pseudobulk_meta(counts, meta)
  filter_design <- stats::model.matrix(filter_formula, data = meta)
  y <- edgeR::DGEList(counts = counts)
  keep <- edgeR::filterByExpr(y, design = filter_design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  # edgeR >= 4.10 renamed calcNormFactors() to normLibSizes().
  if ("normLibSizes" %in% getNamespaceExports("edgeR")) {
    y <- edgeR::normLibSizes(y, method = "TMM")
  } else {
    y <- edgeR::calcNormFactors(y, method = "TMM")
  }
  list(y = y, meta = meta, filter_design = filter_design)
}

run_limma_from_dge <- function(y, meta, design_formula, robust = FALSE, trend = FALSE, voom_plot = FALSE) {
  if (!requireNamespace("limma", quietly = TRUE)) stop("Package 'limma' is required.")
  if (!identical(colnames(y), as.character(meta$patient_id))) {
    idx <- match(colnames(y), meta$patient_id)
    if (anyNA(idx)) stop("Unable to align DGEList columns to metadata.")
    meta <- meta[idx, , drop = FALSE]
  }
  design <- stats::model.matrix(design_formula, data = meta)
  v <- limma::voom(y, design = design, plot = voom_plot)
  fit <- limma::lmFit(v, design)
  fit <- limma::eBayes(fit, trend = trend, robust = robust)
  tt <- limma::topTable(fit, coef = "N_groupNplus", number = Inf, sort.by = "P")
  tt$gene <- rownames(tt)
  list(y = y, voom = v, fit = fit, table = tt, design = design, meta = meta)
}

run_pseudobulk_limma <- function(counts, meta, design_formula, filter_formula = design_formula) {
  prep <- prepare_pseudobulk_dge(counts, meta, filter_formula)
  run_limma_from_dge(prep$y, prep$meta, design_formula)
}

# Write a data.frame containing fgsea/list-columns to a human-readable CSV.
# The original object should still be saved as RDS so list columns (especially
# leadingEdge) remain lossless. CSV gets semicolon-delimited text versions.
write_flat_csv <- function(x, path, row.names = FALSE) {
  y <- as.data.frame(x, stringsAsFactors = FALSE)
  is_list <- vapply(y, is.list, logical(1))
  if (any(is_list)) {
    for (nm in names(y)[is_list]) {
      y[[nm]] <- vapply(
        y[[nm]],
        function(z) {
          if (is.null(z) || length(z) == 0L) return("")
          paste(as.character(z), collapse = ";")
        },
        character(1)
      )
    }
  }
  utils::write.csv(y, path, row.names = row.names)
  invisible(path)
}

z_by_gene <- function(expr) {
  z <- t(scale(t(expr)))
  z
}

patient_mean_z_score <- function(expr, genes) {
  use <- intersect(genes, rownames(expr))
  if (length(use) == 0) stop("No requested genes were found in expression matrix.")
  z <- z_by_gene(expr[use, , drop = FALSE])
  score <- colMeans(z, na.rm = TRUE)
  list(score = score, genes = use, z = z)
}

extract_coef <- function(fit, term = "N_groupNplus") {
  s <- summary(fit)$coefficients
  if (!term %in% rownames(s)) stop("Coefficient not found: ", term)
  est <- unname(s[term, "Estimate"])
  se <- unname(s[term, "Std. Error"])
  p <- unname(s[term, "Pr(>|t|)"])

  # Use the model-specific 95% CI (t-based for lm) rather than the rough
  # estimate +/- 1.96*SE approximation. This matches the saved manuscript
  # tables, including the corrected published-program association table.
  ci <- stats::confint(fit, parm = term, level = 0.95)
  data.frame(beta = est, se = se, p_value = p,
             ci_low = unname(ci[1]), ci_high = unname(ci[2]))
}
