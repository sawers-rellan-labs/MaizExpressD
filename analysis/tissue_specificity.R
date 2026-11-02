#!/usr/bin/env Rscript
# M4: tissue specificity on the seven tissues of Lei et al. 2025 (item M4_tissue_specificity, approved by Rubén
# 30 Oct 2026); L3Mid and LMid are shown in the app but not scored.
#
#   Rscript analysis/tissue_specificity.R      (after analysis/expression_tissue.R)
#
# Input: the line means of results/expression_line.parquet (one value per line and tissue, as Lei's matrices).
# Genes: detected (CPM > 0) in more than half of the lines of at least one of the seven tissues (Lei's expression
# filter, per tissue).
# (1) Lei's method, as in Lei's code (github.com/leilei37/maizeGTEx, Tissue-specificGeneFindingUsingTtest.R; after
#     Finucane et al. 2018), which is narrower than the paper's text: only genes detected in more than half of the
#     lines of all seven tissues (Lei: 15,505); per gene, log(CPM + 0.25) of every line and tissue, scaled, regressed on
#     a tissue indicator (1 for the tested tissue, -1 for the six others), one tissue at a time; t statistic of the
#     tissue effect. A gene is specific to a tissue when its t is in that tissue's top 5% and its median over lines is
#     more than 2-fold the highest median of the other tissues (so only the top tissue can qualify).
# (2) tau (Yanai et al. 2005) on log2(median CPM + 1) of the seven tissues, for every gene above; tissue-specific
#     when tau >= 0.8. NA when no tissue reaches a median of 1 CPM (tau is unstable for barely expressed genes).
# Output: results/tissue_specificity.parquet, one row per gene: gene_id, top_tissue (highest median CPM),
#   max_median_cpm, fold_vs_next (that median / the next highest), lei_scored (detected in all seven tissues),
#   t_<tissue> for the seven tissues (NA unless lei_scored), lei_tissue (NA unless specific by (1)), tau,
#   tau_specific.

source("analysis/common.R")
TOP_SHARE <- 0.05; MIN_FOLD <- 2; TAU_MIN <- 0.8; TAU_MIN_CPM <- 1

summ <- as.data.table(read_parquet("results/expression_tissue.parquet"))[tissue %in% LEI_TISSUES]
genes <- summ[n_detected / n_lines > 0.5, unique(gene_id)]
lei_genes <- summ[n_detected / n_lines > 0.5, .N, by = gene_id][N == length(LEI_TISSUES), gene_id]
lines <- fread("results/expression_lines.tsv")[tissue %in% LEI_TISSUES]
lines[, col := paste(tissue, genotype, sep = "|")]

d <- as.data.table(read_parquet("results/expression_line.parquet"))[tissue %in% LEI_TISSUES & gene_id %in% lei_genes]
y <- matrix(0, length(lei_genes), nrow(lines), dimnames = list(lei_genes, lines$col))   # lines without a row: CPM 0
y[cbind(match(d$gene_id, lei_genes), match(paste(d$tissue, d$genotype, sep = "|"), lines$col))] <- d$cpm
rm(d); invisible(gc())
y <- log(y + 0.25)     # Lei's transform; scaling each gene first leaves t unchanged

n <- ncol(y)
yc <- y - rowMeans(y)
ss_y <- rowSums(yc^2)
tstat <- sapply(LEI_TISSUES, function(tis) {
  x <- ifelse(lines$tissue == tis, 1, -1)
  xc <- x - mean(x); sxx <- sum(xc^2)
  b <- drop(yc %*% xc) / sxx
  s2 <- (ss_y - b^2 * sxx) / (n - 2)
  b / sqrt(s2 / sxx)
})
top5 <- apply(tstat, 2, function(t) rank(-t, ties.method = "first") <= ceiling(TOP_SHARE * length(t)))

med <- dcast(summ[gene_id %in% genes], gene_id ~ tissue, value.var = "median_cpm", fill = 0)
med <- as.matrix(med[match(genes, gene_id), ..LEI_TISSUES])
o <- t(apply(med, 1, order, decreasing = TRUE))
top <- o[, 1]
top_med <- med[cbind(seq_along(genes), top)]
next_med <- med[cbind(seq_along(genes), o[, 2])]
fold <- top_med / next_med                                      # Inf when every other tissue has median 0

lm <- log2(med + 1)
tau <- rowSums(1 - lm / apply(lm, 1, max)) / (ncol(lm) - 1)
tau[top_med < TAU_MIN_CPM] <- NA

res <- data.table(gene_id = genes, top_tissue = factor(LEI_TISSUES[top], levels = LEI_TISSUES),
                  max_median_cpm = top_med, fold_vs_next = round(fold, 3), lei_scored = genes %in% lei_genes)
tcols <- paste0("t_", LEI_TISSUES)
res[, (tcols) := NA_real_]
res[match(lei_genes, gene_id), (tcols) := as.data.table(round(tstat, 3))]
k <- match(lei_genes, genes)
lei <- logical(length(genes))
lei[k] <- fold[k] > MIN_FOLD & top5[cbind(seq_along(lei_genes), top[k])]
res[, `:=`(lei_tissue = factor(ifelse(lei, LEI_TISSUES[top], NA), levels = LEI_TISSUES),
           tau = round(tau, 4), tau_specific = !is.na(tau) & tau >= TAU_MIN)]
write_parquet(res, "results/tissue_specificity.parquet")

message(sprintf("%d genes (detected in > half of the lines of at least one tissue): Lei's t for %d (all seven), tau for %d",
                length(genes), length(lei_genes), sum(!is.na(tau))))
message("tissue  lei  tau>=0.8  both")
for (tis in LEI_TISSUES)
  message(sprintf("%-7s %4d %8d %5d", tis, res[lei_tissue == tis, .N], res[tau_specific & top_tissue == tis, .N],
                  res[lei_tissue == tis & tau_specific, .N]))
