#!/usr/bin/env Rscript
# M4: expression per line and tissue, and per tissue across lines (item M4_expression_per_tissue, approved by Rubén
# 30 Oct 2026).
#
#   Rscript analysis/expression_tissue.R
#
# Libraries and CPM as in analysis/common.R. Per line and tissue: mean CPM of its libraries (repeat libraries
# averaged, M1; B73 = mean of its libraries, both adult-leaf sampling groups for LMAN). Per tissue: quartiles over
# lines, B73 counting as one line; all nine tissues (LMid display only).
# Outputs:
#   results/expression_lines.tsv       one row per line and tissue: tissue, genotype, n_libraries, sample_ids
#   results/expression_line.parquet    gene_id, tissue, genotype, cpm (rows with cpm > 0 only)
#   results/expression_tissue.parquet  gene_id, tissue, n_lines, n_detected (lines with cpm > 0), q25_cpm,
#                                      median_cpm, q75_cpm, b73_cpm (genes detected in the tissue only)

source("analysis/common.R")

cnt <- read_counts(libs$sample_id)
x <- cpm(cnt)
rm(cnt); invisible(gc())

lines <- libs[, .(n_libraries = .N, sample_ids = paste(sample_id, collapse = ";")), by = .(tissue, genotype)]
fwrite(lines, "results/expression_lines.tsv", sep = "\t")

long <- list(); summ <- list()
for (tis in TISSUES) {
  l1 <- libs[tissue == tis]
  g <- factor(l1$genotype)
  W <- sapply(levels(g), function(v) (g == v) / sum(g == v))     # library x line weights: mean of the line's libraries
  m <- x[, l1$sample_id, drop = FALSE] %*% W
  nz <- which(m > 0, arr.ind = TRUE)
  long[[tis]] <- data.table(gene_id = rownames(m)[nz[, 1]], tissue = tis, genotype = colnames(m)[nz[, 2]],
                            cpm = round(m[nz], 3))
  q <- rowQuantiles(m, probs = c(0.25, 0.5, 0.75))
  s <- data.table(gene_id = rownames(m), tissue = tis, n_lines = ncol(m), n_detected = rowSums(m > 0),
                  q25_cpm = q[, 1], median_cpm = q[, 2], q75_cpm = q[, 3],
                  b73_cpm = if ("B73" %in% colnames(m)) m[, "B73"] else NA_real_)
  summ[[tis]] <- s[n_detected > 0]
  message(sprintf("%-7s %4d libraries %4d lines %6d genes detected, %6d with median CPM >= 1", tis, nrow(l1), ncol(m),
                  nrow(summ[[tis]]), sum(s$median_cpm >= 1)))
}
long <- rbindlist(long)
long[, `:=`(tissue = factor(tissue, levels = TISSUES))]
write_parquet(long, "results/expression_line.parquet")
summ <- rbindlist(summ)
summ[, tissue := factor(tissue, levels = TISSUES)]
cols <- c("q25_cpm", "median_cpm", "q75_cpm", "b73_cpm")
summ[, (cols) := lapply(.SD, round, 3), .SDcols = cols]
write_parquet(summ, "results/expression_tissue.parquet")
message(sprintf("expression_line: %d rows; expression_tissue: %d rows, %d genes", nrow(long), nrow(summ), uniqueN(summ$gene_id)))
