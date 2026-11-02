#!/usr/bin/env Rscript
# M5: Lei et al. 2025's results (G3 15:jkaf140) on B73 v5 genes (HANDOFF M5: "Lei tables lifted to v5 gene IDs").
#
#   Rscript annotation/build_lei.R      (after metadata/build/extract_lei_tables.R)
#
# Genes: Lei's v3 IDs (GRMZM..., AC..._FG...) to v5 through MaizeGDB's B73v5_to_B73v3 cross-reference (the v3_ids of
# annotation/genes.tsv). Split or merged models and tandem clusters make some matches many-to-many, so every match is
# kept, one row each, labelled
#   map = 1:1 | 1:n (the v3 ID matches n v5 genes) | n:1 (the v5 gene matches several v3 IDs) | n:m | none
# ("none": no v5 match, gene_id empty; the row is kept so counts follow Lei's tables).
# Positions: "Chromosome 2: 232035498bp" -> chrom chr2 and a number, for v3 and v5 separately (v5 positions are Lei's
# own liftover: "NA" where it failed, and some land on another chromosome or a scaffold). S7 gives v3 positions only:
# its v5 positions come from S4 (same gene, eQTL and tissue) and S6 (same trait and QTL).
# Tissues: Lei's codes to ours (metadata/build/tissues.tsv); Lei's L3 (LT, LD and LN together) -> L3Tip+LMAD+LMAN.
# Merged cells (S8; the identity columns of S10) are filled down.
# Writes annotation/lei_<table>.tsv and .parquet:
#   lei_tissue_specific  S2   gene_v3, gene_id, map, tissue_lei, tissue, gene_name
#   lei_eqtls            S4   gene_v3, gene_id, map, tissue_lei, tissue, eqtl_chrom_v3, eqtl_pos_v3, eqtl_chrom_v5,
#                             eqtl_pos_v5, group_chrom_v5, group_pos_v5 (Lei's grouped eQTL), pvalue
#   lei_twas             S5   gene_v3, gene_id, map, tissue_lei, tissue, traits, category, annotation
#   lei_qtls             S6   trait, chrom_v3, pos_v3, chrom_v5, pos_v5, source
#   lei_coloc            S7 (window 500 kb) and S8 (100 bp): gene_v3, gene_id, map, tissue_lei, tissue, trait, window,
#                             eqtl_chrom_v3, eqtl_pos_v3, eqtl_chrom_v5, eqtl_pos_v5, qtl_pos_v3, qtl_pos_v5 (QTLs as
#                             text, on the eQTL's chromosome: some S8 rows list several positions)
#   lei_cubic            S9   gene_v3, gene_id, map, trait, pvalue
#   lei_candidates       S10  gene_v3, symbol_lei, gene_id, map, tissue_lei, tissue, trait, gene_pos_v5, rho, pvalue,
#                             cubic_pvalue, role_reference
#   lei_traits           S11  trait, description, trait_group (agronomic, tocochromanol, carotenoid)

suppressPackageStartupMessages({ library(data.table); library(arrow) })
SRC <- "metadata/build/source/Lei2025_Table"
rd <- function(s) fread(paste0(SRC, s, ".tsv"), colClasses = "character", na.strings = "", sep = "\t", quote = "")

g <- fread("annotation/genes.tsv", select = c("gene_id", "v3_ids"), colClasses = "character", na.strings = NULL)
xw <- g[v3_ids != "", .(gene_v3 = unlist(strsplit(v3_ids, ",", fixed = TRUE))), by = gene_id]
xw[, n_v5 := .N, by = gene_v3][, n_v3 := .N, by = gene_id]
xw[, map := fifelse(n_v5 == 1, fifelse(n_v3 == 1, "1:1", "n:1"), fifelse(n_v3 == 1, "1:n", "n:m"))]
lift <- function(d) {
  x <- merge(d, xw[, .(gene_v3, gene_id, map)], by = "gene_v3", all.x = TRUE, allow.cartesian = TRUE, sort = FALSE)
  x[is.na(map), map := "none"]
  x
}

tis <- fread("metadata/build/tissues.tsv", select = c("tissue", "lei_code"), na.strings = "")[!is.na(lei_code)]
tissue_of <- function(code) {
  m <- c(setNames(tis$tissue, tis$lei_code), L3 = "L3Tip+LMAD+LMAN")
  out <- unname(m[code]); stopifnot(all(!is.na(out[!is.na(code)]))); out
}
pos_chrom <- function(x) ifelse(grepl("^Chromosome \\w+:", x), paste0("chr", sub("^Chromosome (\\w+):.*$", "\\1", x)), NA_character_)
pos_bp <- function(x) {
  ok <- grepl("^Chromosome \\w+: *[0-9]+ *bp$", x)
  out <- rep(NA_real_, length(x)); out[ok] <- as.numeric(sub("^.*: *([0-9]+) *bp$", "\\1", x[ok])); out
}
fill_down <- function(v) { for (i in seq_along(v)[-1]) if (is.na(v[i])) v[i] <- v[i - 1]; v }
out <- function(d, name) {
  fwrite(d, file.path("annotation", paste0(name, ".tsv")), sep = "\t")
  write_parquet(d, file.path("annotation", paste0(name, ".parquet")))
  n_rows <- nrow(d)
  if ("map" %in% names(d)) {
    u <- unique(d[, .(gene_v3, mapped = map != "none")])
    message(sprintf("%-20s %6d rows; %5d v3 genes, %5.1f%% with a v5 match (%d v5 genes)", name, n_rows, nrow(u),
                    100 * mean(u$mapped), uniqueN(d[map != "none", gene_id])))
  } else message(sprintf("%-20s %6d rows", name, n_rows))
}

s2 <- lift(rd("S2_tissue_specific")[, .(gene_v3, tissue_lei, gene_name = fifelse(gene_name == "NA", NA_character_, gene_name))])
s2[, tissue := tissue_of(tissue_lei)]
out(s2[, .(gene_v3, gene_id, map, tissue_lei, tissue, gene_name)], "lei_tissue_specific")

s4 <- rd("S4_eqtls")
stopifnot(!anyNA(pos_bp(s4$eqtl_v3)))
message(sprintf("S4: %d of %d eQTLs without a v5 position, %d on another chromosome or scaffold in v5", sum(is.na(pos_bp(s4$eqtl_v5))),
                nrow(s4), s4[!is.na(pos_bp(eqtl_v5)) & pos_chrom(eqtl_v3) != pos_chrom(eqtl_v5), .N]))
s4 <- s4[, .(gene_v3, tissue_lei, tissue = tissue_of(tissue_lei), eqtl_chrom_v3 = pos_chrom(eqtl_v3), eqtl_pos_v3 = pos_bp(eqtl_v3),
             eqtl_chrom_v5 = pos_chrom(eqtl_v5), eqtl_pos_v5 = pos_bp(eqtl_v5), group_chrom_v5 = pos_chrom(group_v5),
             group_pos_v5 = pos_bp(group_v5), pvalue = as.numeric(pvalue))]
out(lift(s4)[, .(gene_v3, gene_id, map, tissue_lei, tissue, eqtl_chrom_v3, eqtl_pos_v3, eqtl_chrom_v5, eqtl_pos_v5,
                 group_chrom_v5, group_pos_v5, pvalue)], "lei_eqtls")

s5 <- lift(rd("S5_twas"))
s5[, tissue := tissue_of(tissue_lei)]
out(s5[, .(gene_v3, gene_id, map, tissue_lei, tissue, traits, category, annotation)], "lei_twas")

s6 <- rd("S6_qtls")
out(s6[, .(trait, chrom_v3 = pos_chrom(qtl_v3), pos_v3 = pos_bp(qtl_v3), chrom_v5 = pos_chrom(qtl_v5), pos_v5 = pos_bp(qtl_v5),
           source)], "lei_qtls")

s4r <- rd("S4_eqtls")
e5 <- unique(s4r[, .(gene_v3, eqtl_v3, tissue_lei, eqtl_v5)])
q5 <- unique(s6[, .(trait, qtl_v3, qtl_v5)])
stopifnot(!anyDuplicated(e5, by = c("gene_v3", "eqtl_v3", "tissue_lei")), !anyDuplicated(q5, by = c("trait", "qtl_v3")))
s7 <- rd("S7_coloc_500kb")
n7 <- nrow(s7)
s7 <- merge(merge(s7, e5, by = c("gene_v3", "eqtl_v3", "tissue_lei"), all.x = TRUE, sort = FALSE),
            q5, by = c("trait", "qtl_v3"), all.x = TRUE, sort = FALSE)
stopifnot(nrow(s7) == n7)
message(sprintf("S7: eQTL found in S4 for %.1f%% of rows (v5 position %.1f%%), QTL found in S6 for %.1f%% (v5 position %.1f%%)",
                100 * mean(!is.na(s7$eqtl_v5)), 100 * mean(!is.na(pos_bp(s7$eqtl_v5))), 100 * mean(!is.na(s7$qtl_v5)),
                100 * mean(!is.na(pos_bp(s7$qtl_v5)))))
s8 <- rd("S8_coloc_100bp")
for (k in names(s8)) set(s8, j = k, value = fill_down(s8[[k]]))
s8 <- s8[, .(tissue_lei = unlist(strsplit(tissue_lei, " +"))), by = .(gene_v3, trait, eqtl_v3, eqtl_v5, qtl_v3, qtl_v5)]
qtext <- function(x) ifelse(is.na(x) | x == "NA", NA_character_, gsub("Chromosome \\w+: *|bp", "", x))
coloc <- function(d, w) d[, .(gene_v3, tissue_lei, trait, window = w, eqtl_chrom_v3 = pos_chrom(eqtl_v3), eqtl_pos_v3 = pos_bp(eqtl_v3),
                              eqtl_chrom_v5 = pos_chrom(eqtl_v5), eqtl_pos_v5 = pos_bp(eqtl_v5), qtl_pos_v3 = qtext(qtl_v3),
                              qtl_pos_v5 = qtext(qtl_v5))]
co <- lift(rbind(coloc(s7, "500 kb"), coloc(s8, "100 bp")))
co[, tissue := tissue_of(tissue_lei)]
out(co[, .(gene_v3, gene_id, map, tissue_lei, tissue, trait, window, eqtl_chrom_v3, eqtl_pos_v3, eqtl_chrom_v5, eqtl_pos_v5,
           qtl_pos_v3, qtl_pos_v5)], "lei_coloc")

s9 <- lift(rd("S9_cubic")[, .(gene_v3, trait, pvalue = as.numeric(pvalue))])
out(s9[, .(gene_v3, gene_id, map, trait, pvalue)], "lei_cubic")

s10 <- rd("S10_candidates")
for (k in c("eqtl_v3", "eqtl_v5", "qtl_v3", "qtl_v5", "gene_v3", "tissue_lei", "trait", "gene_pos_v3", "gene_pos_v5",
            "role_reference")) set(s10, j = k, value = fill_down(s10[[k]]))
s10[, `:=`(symbol_lei = ifelse(grepl("\\(", gene_v3), sub("^.*\\((.*)\\)$", "\\1", gene_v3), NA_character_),
           gene_v3 = sub("\\(.*$", "", gene_v3))]
s10 <- lift(s10)
s10[, tissue := tissue_of(tissue_lei)]
out(s10[, .(gene_v3, symbol_lei, gene_id, map, tissue_lei, tissue, trait, gene_pos_v5, rho = as.numeric(rho),
            pvalue = as.numeric(pvalue), cubic_pvalue = as.numeric(cubic_pvalue), role_reference)], "lei_candidates")

s11 <- rd("S11_traits")
tr <- rbind(s11[!is.na(agronomic_trait), .(trait = agronomic_trait, description = trait_protocol, trait_group = "agronomic")],
            s11[!is.na(tocochromanol_trait), .(trait = tocochromanol_trait, description = NA_character_, trait_group = "tocochromanol")],
            s11[!is.na(carotenoid_trait), .(trait = carotenoid_trait, description = NA_character_, trait_group = "carotenoid")])
out(tr, "lei_traits")
