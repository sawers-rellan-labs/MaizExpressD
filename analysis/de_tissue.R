#!/usr/bin/env Rscript
# M4: differential expression between tissues for every contrast in metadata/contrasts.tsv.
#
#   Rscript analysis/de_tissue.R                  # all contrasts -> results/de.parquet, results/de_summary.tsv
#   Rscript analysis/de_tissue.R <contrast_id>... # selected contrasts -> results/de_parts/<contrast_id>.parquet/.tsv
#   Rscript analysis/de_tissue.R --merge          # merge results/de_parts/ into results/de.parquet and de_summary.tsv
#   Rscript analysis/de_tissue.R --cpm            # recompute the CPM columns of the parts in results/de_parts/, no refit
#
# Adapted from MaizExpressC's analysis/de_per_study.R (approved at the M4 checkpoint, Rubén 30 Oct 2026):
#   panel: the diversity panel, paired by line (pair_key: line, or line and sampling group for day/night):
#     edgeR quasi-likelihood ~ pair + tissue (TMM, glmQLFit robust, QL F test of tissue) on the libraries passing
#     library QC of the pairs with libraries on both sides; repeat libraries of a line enter as extra libraries of
#     its pair. Not DESeq2 as in C's Zhou2022 panel (23-25 inbreds): with 180+ pairs (~186 coefficients) DESeq2's
#     Cox-Reid gene-wise dispersion search failed for 38-55% of genes (dispersion ~118 instead of ~0.01, one inflated
#     lfcSE for all of them); edgeR QL estimates the extra variance from the residual df (Rubén, 30 Oct 2026). lfc_se
#     is |log2FC| / sqrt(F) and stat sign(log2FC) * sqrt(F), so ashr shrinks it as lfcShrink does for DESeq2;
#   B73: the B73 libraries alone, DESeq2 ~ tissue (replicated groups).
# Libraries in metadata/library_exclusions.tsv are left out. Counts without length offsets (3' tags quantified
# with salmon --noLengthCorrection). log2FC = log2(numerator / denominator), as in contrasts.tsv. Unshrunk log2FC
# and lfcSE are kept; ashr-shrunk log2FC is for display. A gene is tested if it has >= MIN_COUNT reads in
# >= MIN_SHARE of the fit's libraries (and at least 2), so sparse genes do not slow the large panel fits.
# single_library (MaizExpressN rule, as in C): the fold change rests on one library; such genes keep their log2FC
# for display but get pvalue and padj NA.
# Outputs: gene x contrast table (gene_id, contrast_id, base_mean, log2fc, lfc_se, log2fc_shrunk, stat, pvalue,
# padj, single_library, cpm_numerator_mean, cpm_denominator_mean) and one summary row per contrast. The CPM columns
# are the mean CPM (reads on genes per million) of the fit's numerator and denominator libraries, not salmon's TPM,
# which divides 3' tags by transcript length; the parts fitted on 30-31 Oct 2026 got them with --cpm (they had TPM).

suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
  library(arrow)
  library(BiocParallel)
  library(edgeR)
})
stopifnot(requireNamespace("ashr", quietly = TRUE))
options(warn = 1)
register(MulticoreParam(as.integer(Sys.getenv("DE_WORKERS", "6"))))

PADJ_DEG <- 0.05; LFC_DEG <- 1      # DEG call for the summary table only
MIN_COUNT <- 10; MIN_SHARE <- 0.1   # gene filter (see header)
SPIKE_SHARE <- 0.9; SPIKE_PADJ <- 0.05
GENE_LIMIT <- as.integer(Sys.getenv("DE_GENE_LIMIT", "0"))   # > 0: first N tested genes only (timing runs)
PARTS <- "results/de_parts"

samples <- fread("metadata/samples.tsv", colClasses = "character", na.strings = NULL)
cons <- fread("metadata/contrasts.tsv", colClasses = "character", na.strings = NULL)
excl <- fread("metadata/library_exclusions.tsv", colClasses = "character")
use <- samples[include == "TRUE" & !(sample_id %in% excl$sample_id)]
split_ids <- function(x) setdiff(strsplit(x, ";", fixed = TRUE)[[1]], "")

is_single_library <- function(N, D, lfc, padj) {
  sn <- rowSums(N); sd_ <- rowSums(D); mn <- apply(N, 1, max); md <- apply(D, 1, max)
  share <- pmax(mn, md) / (sn + sd_)
  up <- lfc > 0
  rest_hi <- ifelse(up, (sn - mn) / (ncol(N) - 1), (sd_ - md) / (ncol(D) - 1))
  lo <- ifelse(up, sd_ / ncol(D), sn / ncol(N))
  flag <- (sn + sd_ > 0 & share > SPIKE_SHARE) | (!is.na(padj) & padj < SPIKE_PADJ & rest_hi <= lo)
  flag %in% TRUE
}

read_all <- function(file, ids = NULL) {    # ids: only these libraries (keeps parallel runs small)
  tabs <- lapply(Sys.glob(file.path("results", "Kremling2018", "*", file)), function(f) {
    if (is.null(ids)) return(fread(f))
    h <- names(fread(f, nrows = 0))
    fread(f, select = c(h[1:2], intersect(h, ids)))
  })
  genes <- tabs[[1]]$gene_id
  stopifnot(all(sapply(tabs, function(d) identical(d$gene_id, genes))))
  tabs <- tabs[sapply(tabs, ncol) > 2]       # batches holding none of the libraries
  m <- do.call(cbind, lapply(tabs, function(d) as.matrix(d[, -(1:2)])))
  rownames(m) <- genes
  m
}

cpm_means <- function(cnt, num, den, genes) {
  cp <- t(t(cnt[, c(num, den), drop = FALSE]) / colSums(cnt[, c(num, den), drop = FALSE])) * 1e6
  list(cpm_numerator_mean = rowMeans(cp[genes, num, drop = FALSE]), cpm_denominator_mean = rowMeans(cp[genes, den, drop = FALSE]))
}

contrast_libraries <- function(c1) {   # the libraries a contrast is fitted on: QC, and for the panel the pairs with both tissues
  num <- intersect(split_ids(c1$numerator_samples), use$sample_id)
  den <- intersect(split_ids(c1$denominator_samples), use$sample_id)
  dropped <- setdiff(c(split_ids(c1$numerator_samples), split_ids(c1$denominator_samples)), c(num, den))
  meta <- as.data.frame(use[match(c(num, den), sample_id)])
  rownames(meta) <- meta$sample_id
  meta$tissue <- factor(ifelse(meta$sample_id %in% num, c1$numerator_tissue, c1$denominator_tissue),
                        levels = c(c1$denominator_tissue, c1$numerator_tissue))
  if (c1$design == "panel") {
    meta$pair <- if (c1$pair_key == "genotype|sampling_group") paste(meta$genotype, meta$sampling_group, sep = "|") else meta$genotype
    both <- intersect(meta$pair[meta$tissue == c1$numerator_tissue], meta$pair[meta$tissue == c1$denominator_tissue])
    dropped <- c(dropped, meta$sample_id[!meta$pair %in% both])
    meta <- meta[meta$pair %in% both, ]
    meta$pair <- factor(meta$pair)
    design <- ~ pair + tissue
  } else {
    design <- ~ tissue
  }
  num <- intersect(num, rownames(meta)); den <- intersect(den, rownames(meta))
  list(meta = meta, num = num, den = den, dropped = dropped, design = design)
}

run_contrast <- function(c1, cnt) {
  L <- contrast_libraries(c1)
  meta <- L$meta; num <- L$num; den <- L$den; dropped <- L$dropped; design <- L$design
  x <- round(cnt[, rownames(meta)])
  keep <- rowSums(x >= MIN_COUNT) >= max(2, ceiling(MIN_SHARE * ncol(x)))
  if (GENE_LIMIT > 0) keep[which(keep)[-seq_len(min(GENE_LIMIT, sum(keep)))]] <- FALSE
  t0 <- Sys.time()
  if (c1$design == "panel") {
    y <- normLibSizes(DGEList(x[keep, ]))
    X <- model.matrix(design, meta)
    stopifnot(colnames(X)[ncol(X)] == paste0("tissue", c1$numerator_tissue))
    fit <- glmQLFit(y, X, robust = TRUE)
    tt <- glmQLFTest(fit, coef = ncol(X))$table
    tt$F <- pmax(tt$F, 0)     # edgeR 4 QL can give F slightly below 0 for genes with no change
    se <- ifelse(tt$F > 0, abs(tt$logFC) / sqrt(tt$F), NA_real_)
    ok <- is.finite(se) & se > 0
    shr <- rep(NA_real_, nrow(tt))
    shr[ok] <- ashr::ash(tt$logFC[ok], se[ok], mixcompdist = "normal", method = "shrink")$result$PosteriorMean
    eff <- y$samples$lib.size * y$samples$norm.factors
    nc <- t(t(x[keep, ]) / (eff / exp(mean(log(eff)))))
    d <- data.table(gene_id = rownames(tt), contrast_id = c1$contrast_id, base_mean = rowMeans(nc),
                    log2fc = tt$logFC, lfc_se = se, log2fc_shrunk = shr, stat = sign(tt$logFC) * sqrt(tt$F),
                    pvalue = tt$PValue, padj = p.adjust(tt$PValue, "BH"))
    note <- sprintf("prior QL df median %.1f", median(fit$df.prior))
  } else {
    dds <- DESeq(DESeqDataSetFromMatrix(x[keep, ], meta, design), quiet = TRUE, parallel = TRUE)
    rn <- paste0("tissue_", c1$numerator_tissue, "_vs_", c1$denominator_tissue)
    res <- results(dds, name = rn, alpha = PADJ_DEG)
    sh <- suppressMessages(lfcShrink(dds, coef = rn, res = res, type = "ashr", quiet = TRUE))
    nc <- counts(dds, normalized = TRUE)
    d <- data.table(gene_id = rownames(res), contrast_id = c1$contrast_id, base_mean = res$baseMean,
                    log2fc = res$log2FoldChange, lfc_se = res$lfcSE, log2fc_shrunk = sh$log2FoldChange,
                    stat = res$stat, pvalue = res$pvalue, padj = res$padj)
    note <- ""
  }
  secs <- as.numeric(Sys.time() - t0, units = "secs")
  d[, single_library := is_single_library(nc[gene_id, num, drop = FALSE], nc[gene_id, den, drop = FALSE], log2fc, padj)]
  d[single_library == TRUE, `:=`(pvalue = NA_real_, padj = NA_real_)]
  d[, c("cpm_numerator_mean", "cpm_denominator_mean") := cpm_means(cnt, num, den, d$gene_id)]
  s <- data.table(contrast_id = c1$contrast_id, contrast_class = c1$contrast_class, design = c1$design,
                  model = if (c1$design == "panel") "edgeR QL ~ pair + tissue" else "DESeq2 ~ tissue",
                  n_pairs = if (c1$design == "panel") nlevels(meta$pair) else NA_integer_,
                  n_numerator = length(num), n_denominator = length(den),
                  excluded_libraries = length(unique(dropped)), genes_tested = nrow(d),
                  n_up = d[padj < PADJ_DEG & log2fc >= LFC_DEG, .N], n_down = d[padj < PADJ_DEG & log2fc <= -LFC_DEG, .N],
                  n_padj = d[padj < PADJ_DEG, .N], n_single_library = d[single_library == TRUE, .N], fit_seconds = round(secs), note = note)
  list(de = d, summary = s)
}

args <- commandArgs(trailingOnly = TRUE)
if (identical(args, "--merge")) {
  ids <- intersect(cons$contrast_id, sub("\\.parquet$", "", list.files(PARTS, "\\.parquet$")))
  stopifnot(length(ids) == nrow(cons))
  write_parquet(rbindlist(lapply(ids, function(i) read_parquet(file.path(PARTS, paste0(i, ".parquet"))))), "results/de.parquet")
  s <- rbindlist(lapply(ids, function(i) fread(file.path(PARTS, paste0(i, ".tsv")))), fill = TRUE)
  s[model == "~ tissue", model := "DESeq2 ~ tissue"]      # B73 parts written before the panel switch to edgeR
  fwrite(s, "results/de_summary.tsv", sep = "\t")
  message(sprintf("merged %d contrasts into results/de.parquet and results/de_summary.tsv", length(ids)))
  quit(save = "no")
}
if (identical(args, "--cpm")) {
  for (f in list.files(PARTS, "\\.parquet$", full.names = TRUE)) {
    d <- as.data.table(read_parquet(f))
    c1 <- cons[contrast_id == d$contrast_id[1]]
    L <- contrast_libraries(c1)
    s <- fread(sub("parquet$", "tsv", f))
    stopifnot(length(L$num) == s$n_numerator, length(L$den) == s$n_denominator)   # the fit's libraries
    old <- intersect(c("tpm_numerator_mean", "tpm_denominator_mean"), names(d))
    if (length(old)) d[, (old) := NULL]
    d[, c("cpm_numerator_mean", "cpm_denominator_mean") := cpm_means(read_all("salmon.merged.gene_counts.tsv", c(L$num, L$den)),
                                                                     L$num, L$den, d$gene_id)]
    write_parquet(d, f)
    message(sprintf("%-52s CPM columns from %d + %d libraries", c1$contrast_id, length(L$num), length(L$den)))
  }
  quit(save = "no")
}
todo <- if (length(args)) cons[contrast_id %in% args] else cons
stopifnot(all(args %in% cons$contrast_id))
ids <- unique(unlist(lapply(c(todo$numerator_samples, todo$denominator_samples), split_ids)))
cnt <- read_all("salmon.merged.gene_counts.tsv", ids)
out <- list(); summ <- list()
for (i in seq_len(nrow(todo))) {
  c1 <- todo[i]
  r <- run_contrast(c1, cnt)
  message(sprintf("%-52s %5d genes  %4d up  %4d down  %6d s", c1$contrast_id, r$summary$genes_tested, r$summary$n_up,
                  r$summary$n_down, r$summary$fit_seconds))
  if (length(args) && GENE_LIMIT == 0) {
    dir.create(PARTS, showWarnings = FALSE)
    write_parquet(r$de, file.path(PARTS, paste0(c1$contrast_id, ".parquet")))
    fwrite(r$summary, file.path(PARTS, paste0(c1$contrast_id, ".tsv")), sep = "\t")
  }
  out[[i]] <- r$de; summ[[i]] <- r$summary
}
if (!length(args) && GENE_LIMIT == 0) {
  write_parquet(rbindlist(out), "results/de.parquet")
  fwrite(rbindlist(summ), "results/de_summary.tsv", sep = "\t")
}
