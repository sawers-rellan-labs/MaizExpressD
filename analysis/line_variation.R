#!/usr/bin/env Rscript
# M4: variation among lines, per tissue (item M4_line_variation, approved by Rubén 30 Oct 2026).
#
#   Rscript analysis/line_variation.R      (after analysis/expression_tissue.R)
#
# Libraries and CPM as in analysis/common.R. Per tissue and expressed gene (median CPM over lines >= 5, from
# results/expression_tissue.parquet): v_lines = variance of log2(CPM + 1) among lines, one library per line (the one
# with most reads on genes when a line has repeats; B73 is one of the lines), and v_b73 = the same variance among the
# B73 libraries of that tissue (plant, batch and depth noise within one line; at least 3 libraries).
# share_beyond_b73 = (v_lines - v_b73) / v_lines, 0 when negative; NA without v_b73 (L3Mid, LMid). For LMAN, B73 has
# libraries in both sampling groups, as the lines do; for LMAD, B73 was sampled in group 26 only, so its v_b73 lacks
# the group difference that v_lines holds.
# Output: results/line_variation.parquet: gene_id, tissue, n_lines, median_cpm, v_lines, n_b73, v_b73, share_beyond_b73.
# Also results/b73_unstable.tsv (M5): genes whose B73 libraries disagree on and off in a tissue where B73 has >= 3
# libraries and a mean >= 5 CPM: variance of log2(CPM + 1) among them > 4 (SD > 2 log2 units; the 99th percentile is
# 0.8-2.4 by tissue), pooled within sampling group (B73's night libraries come from both adult-leaf groups, whose dates
# differ; 24 genes that differ only between the dates are not flagged). Mostly nuclear or scaffold copies of
# chloroplast and mitochondrial genes (psbA, rbcL, rps4, petA,
# atpA, ccmFN) that collect organelle reads unevenly, so their CPM and tissue specificity are not biological; the app
# flags them. Columns: gene_id, tissue, n_b73, b73_cpm, v_b73, share_b73_zero (B73 libraries with no reads).

source("analysis/common.R")
MIN_MEDIAN_CPM <- 5; MIN_B73 <- 3; UNSTABLE_V_B73 <- 4

summ <- as.data.table(read_parquet("results/expression_tissue.parquet"))
cnt <- read_counts(libs$sample_id)
depth <- colSums(cnt)
x <- log2(cpm(cnt) + 1)
rm(cnt); invisible(gc())

libs[, depth := depth[sample_id]]
one <- libs[order(-depth)][, .SD[1], by = .(tissue, genotype)]      # one library per line and tissue
res <- list()
for (tis in TISSUES) {
  g <- summ[tissue == tis & median_cpm >= MIN_MEDIAN_CPM, gene_id]
  ids <- one[tissue == tis, sample_id]
  b73 <- libs[tissue == tis & genotype == "B73", sample_id]
  v_lines <- rowVars(x[g, ids, drop = FALSE])
  v_b73 <- if (length(b73) >= MIN_B73) rowVars(x[g, b73, drop = FALSE]) else rep(NA_real_, length(g))
  res[[tis]] <- data.table(gene_id = g, tissue = tis, n_lines = length(ids),
                           median_cpm = summ[tissue == tis][match(g, gene_id), median_cpm],
                           v_lines = v_lines, n_b73 = length(b73), v_b73 = v_b73,
                           share_beyond_b73 = pmax(0, (v_lines - v_b73) / v_lines))
  r <- res[[tis]]
  message(sprintf("%-7s %4d lines %3d B73 libraries %6d genes: median v_lines %.3f, v_b73 %s, share beyond B73 %s (zero for %s)",
                  tis, length(ids), length(b73), length(g), median(v_lines),
                  if (length(b73) >= MIN_B73) sprintf("%.3f", median(v_b73)) else "NA",
                  if (length(b73) >= MIN_B73) sprintf("%.2f", median(r$share_beyond_b73)) else "NA",
                  if (length(b73) >= MIN_B73) sprintf("%.1f%%", 100 * mean(r$share_beyond_b73 == 0)) else "NA"))
}
res <- rbindlist(res)
res[, tissue := factor(tissue, levels = TISSUES)]
cols <- c("v_lines", "v_b73", "share_beyond_b73")
res[, (cols) := lapply(.SD, round, 4), .SDcols = cols]
write_parquet(res, "results/line_variation.parquet")
message(sprintf("line_variation: %d rows", nrow(res)))

uns <- rbindlist(lapply(TISSUES, function(tis) {
  b <- libs[tissue == tis & genotype == "B73"]
  g <- summ[tissue == tis & b73_cpm >= MIN_MEDIAN_CPM, gene_id]
  if (nrow(b) < MIN_B73 || !length(g)) return(NULL)
  grp <- Filter(function(ids) length(ids) >= 2, split(b$sample_id, b$sampling_group))
  ss <- Reduce(`+`, lapply(grp, function(ids) (length(ids) - 1) * rowVars(x[g, ids, drop = FALSE])))
  df <- sum(lengths(grp)) - length(grp)
  xb <- x[g, b$sample_id, drop = FALSE]
  d <- data.table(gene_id = g, tissue = tis, n_b73 = nrow(b), b73_cpm = summ[tissue == tis][match(g, gene_id), b73_cpm],
                  v_b73 = round(ss / df, 3), share_b73_zero = round(rowMeans(xb == 0), 3))
  d[v_b73 > UNSTABLE_V_B73]
}))
fwrite(uns[order(gene_id, tissue)], "results/b73_unstable.tsv", sep = "\t")
message(sprintf("b73_unstable: %d genes (%d gene x tissue rows)", uniqueN(uns$gene_id), nrow(uns)))
