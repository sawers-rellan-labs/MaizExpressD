# Shared by the M4 expression scripts (expression_tissue.R, tissue_specificity.R, line_variation.R): the libraries
# they use and the gene count matrix.
#
# Libraries: include, passing library QC (metadata/library_exclusions.tsv), without B73-10 (M1: out of line-level
# statistics) and without the pilot L3Base libraries (M3: out of the L3Base line values; flowcell HHV55BGXX).
# CPM = gene reads / the library's reads on genes x 1e6, from salmon's gene counts. Not salmon's TPM: with
# --noLengthCorrection salmon still divides reads by the full transcript length when it reports TPM (EffectiveLength =
# Length in quant.sf), which 3' tags (one tag per transcript whatever its length) do not need.

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(matrixStats)
})

PILOT_FC <- "HHV55BGXX"
TISSUES <- c("GRoot", "GShoot", "L3Base", "L3Mid", "LMid", "L3Tip", "LMAD", "LMAN", "Kern")   # display order
LEI_TISSUES <- c("GRoot", "GShoot", "L3Base", "L3Tip", "LMAD", "LMAN", "Kern")               # Lei et al. 2025's seven

samples <- fread("metadata/samples.tsv", colClasses = "character", na.strings = NULL)
excl <- fread("metadata/library_exclusions.tsv", colClasses = "character")
libs <- samples[include == "TRUE" & !(sample_id %in% excl$sample_id) & genotype != "B73-10" &
                !(tissue == "L3Base" & flowcell == PILOT_FC)]
libs[, tissue := factor(tissue, levels = TISSUES)]
setorder(libs, tissue, genotype, sample_id)

read_counts <- function(ids) {   # gene x library counts for these libraries, from every batch's merged table
  tabs <- lapply(Sys.glob(file.path("results", "Kremling2018", "*", "salmon.merged.gene_counts.tsv")), function(f) {
    h <- names(fread(f, nrows = 0))
    fread(f, select = c(h[1:2], intersect(h, ids)))
  })
  genes <- tabs[[1]]$gene_id
  stopifnot(all(sapply(tabs, function(d) identical(d$gene_id, genes))))
  tabs <- tabs[sapply(tabs, ncol) > 2]
  m <- do.call(cbind, lapply(tabs, function(d) as.matrix(d[, -(1:2)])))
  rownames(m) <- genes
  stopifnot(setequal(colnames(m), ids))
  m[, ids, drop = FALSE]
}

cpm <- function(x) t(t(x) / colSums(x)) * 1e6
