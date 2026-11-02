# Extract tables from the Lei et al. 2025 supplement (G3 15:jkaf140, doi:10.1093/g3journal/jkaf140)
# into TSV. Source: metadata/build/source/Lei2025_jkaf140_supplementary_tables.xlsx, from the
# article's supplementary zip (Europe PMC supplementaryFiles API for PMC12405892).
# M1 needs Table S1 (the 298 lines); M5 adds S2-S12 (tissue-specific genes, eQTLs, TWAS, QTLs, their
# colocalization, CUBIC check, candidate genes, trait names and QTL sources), written as stored: every cell as
# text, merged cells left empty below their first row (S8, S10), v3 gene IDs. annotation/build_lei.R lifts them
# to v5 genes. Note rows are dropped: S5 "Category 1: Genes significantly associated in only one tissue. Category 2:
# Genes associated in multiple tissues."; S8 "DTT: Days to Tasseling, LW: Leaf Width, TL: Tassel length, TBN: Tassel
# Branch Number." S10's full references (rows from "Reference:" on) go to Lei2025_TableS10_references.tsv.
#
# Usage: Rscript metadata/build/extract_lei_tables.R

suppressMessages(library(readxl))

src <- "metadata/build/source/Lei2025_jkaf140_supplementary_tables.xlsx"

s1 <- read_excel(src, sheet = "Table S1", skip = 1)
names(s1) <- c("inbred_lei", "origin", "pedigree", "subpopulation")
s1 <- s1[!is.na(s1$inbred_lei), ]
stopifnot(nrow(s1) == 298, !anyDuplicated(s1$inbred_lei))
write.table(s1, "metadata/build/source/Lei2025_TableS1_lines.tsv", sep = "\t", quote = FALSE,
            row.names = FALSE, na = "")
cat("Table S1:", nrow(s1), "lines\n")

tables <- list(
  S2 = list("tissue_specific", c("gene_v3", "tissue_lei", "gene_name")),
  S3 = list("enrichment", c("tissue_lei", "description", "padj", "enrichment_fold", "type")),
  S4 = list("eqtls", c("gene_v3", "eqtl_v3", "eqtl_v5", "group_v3", "group_v5", "tissue_lei", "pvalue")),
  S5 = list("twas", c("gene_v3", "tissue_lei", "traits", "category", "annotation")),
  S6 = list("qtls", c("trait", "qtl_v3", "qtl_v5", "source")),
  S7 = list("coloc_500kb", c("trait", "qtl_v3", "gene_v3", "eqtl_v3", "tissue_lei")),
  S8 = list("coloc_100bp", c("eqtl_v3", "eqtl_v5", "qtl_v3", "qtl_v5", "gene_v3", "trait", "tissue_lei")),
  S9 = list("cubic", c("gene_v3", "trait", "pvalue")),
  S10 = list("candidates", c("eqtl_v3", "eqtl_v5", "qtl_v3", "qtl_v5", "gene_v3", "tissue_lei", "trait", "gene_pos_v3",
                             "gene_pos_v5", "rho", "pvalue", "cubic_pvalue", "role_reference")),
  S11 = list("traits", c("agronomic_trait", "trait_protocol", "tocochromanol_trait", "carotenoid_trait")),
  S12 = list("qtl_sources", c("population", "traits", "qtls", "reference")))
for (s in names(tables)) {
  x <- suppressMessages(read_excel(src, sheet = paste("Table", s), skip = 1, col_types = "text"))
  x <- x[, seq_along(tables[[s]][[2]])]                 # the sheets' extra columns are empty
  names(x) <- tables[[s]][[2]]
  x <- x[rowSums(!is.na(x)) > 0, ]
  if (s == "S5") x <- x[!grepl("^Note:", x$gene_v3), ]
  if (s == "S8") x <- x[!grepl("^Note:", x$eqtl_v3), ]
  if (s == "S10") {
    r <- which(grepl("^Reference:", x$eqtl_v3))
    refs <- data.frame(reference = sub("^Reference: *", "", x$eqtl_v3[r:nrow(x)]))
    write.table(refs, "metadata/build/source/Lei2025_TableS10_references.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
    x <- x[seq_len(r - 1), ]
  }
  write.table(x, sprintf("metadata/build/source/Lei2025_Table%s_%s.tsv", s, tables[[s]][[1]]), sep = "\t", quote = FALSE,
              row.names = FALSE, na = "")
  cat(sprintf("Table %s: %d rows (%s)\n", s, nrow(x), tables[[s]][[1]]))
}
