#!/usr/bin/env Rscript
# M5: turn annotation/genes.tsv and go_terms.tsv (annotation/build_genes.py) into Parquet, and
# build a search index of every identifier and name.
#   Rscript annotation/build_parquet.R
# Writes annotation/genes.parquet, annotation/go_terms.parquet and annotation/gene_aliases.parquet
# (and, once annotation/fetch_maizemine.py and fetch_literature.py have run, pathways, gwas_hits, variants,
# pangenes, atlas, atlas_samples and literature .parquet)
# (alias, alias_lower, alias_type, gene_id). Search = exact match on alias_lower; an alias that
# maps to several genes returns all of them (v3 and v4 IDs often do).

suppressPackageStartupMessages({ library(data.table); library(arrow) })

g <- fread("annotation/genes.tsv", colClasses = "character", na.strings = NULL)
g[, `:=`(start = as.integer(start), end = as.integer(end), quantified = quantified == "TRUE")]
write_parquet(g, "annotation/genes.parquet")
write_parquet(fread("annotation/go_terms.tsv"), "annotation/go_terms.parquet")

split_col <- function(col, type, strip = identity) {
  d <- g[get(col) != "", .(alias = strip(unlist(strsplit(get(col), ",", fixed = TRUE)))), by = gene_id]
  d[alias != "", .(alias, alias_type = type, gene_id)]
}
ath <- g[ath_orthologs != "", .(lbl = unlist(strsplit(ath_orthologs, ",", fixed = TRUE))), by = gene_id]
# "AGI:symbol" splits at the first colon only (Ensembl writes some symbols with one, e.g. NRT2:1)
ath[, `:=`(agi = sub(":.*$", "", lbl), sym = ifelse(grepl(":", lbl), sub("^[^:]*:", "", lbl), NA_character_))]
aliases <- rbindlist(list(
  g[, .(alias = gene_id, alias_type = "v5", gene_id)],
  split_col("v4_ids", "v4"),
  split_col("v3_ids", "v3"),
  g[symbol != "", .(alias = symbol, alias_type = "symbol", gene_id)],
  split_col("synonyms", "synonym", function(x) sub("^NCBI:", "", x)),
  ath[, .(alias = agi, alias_type = "arabidopsis", gene_id)],
  ath[!is.na(sym) & sym != "", .(alias = sym, alias_type = "arabidopsis", gene_id)],
  ath[!is.na(sym) & grepl(":", sym), .(alias = gsub(":", ".", sym), alias_type = "arabidopsis", gene_id)]  # NRT2:1 -> NRT2.1
))
aliases[, alias_lower := tolower(alias)]
aliases <- unique(aliases, by = c("alias_lower", "gene_id", "alias_type"))
setcolorder(aliases, c("alias", "alias_lower", "alias_type", "gene_id"))
write_parquet(aliases, "annotation/gene_aliases.parquet")

cat(sprintf("genes.parquet: %d genes; gene_aliases.parquet: %d aliases (%d distinct)\n",
            nrow(g), nrow(aliases), uniqueN(aliases$alias_lower)))
print(aliases[, .N, by = alias_type])
for (q in c("pep1", "sh2", "zp27", "CCA1", "Zm00001d046170", "GRMZM2G083841")) {
  hit <- aliases[alias_lower == tolower(q), unique(gene_id)]
  cat(sprintf("%-15s -> %s\n", q, paste(head(hit, 6), collapse = ", ")))
}

# MaizeMine extras (annotation/fetch_maizemine.py): pathways, GWAS hits near genes, protein-changing
# variant counts and pan-gene presence, one Parquet each
if (file.exists("annotation/pangenes.tsv")) {
  for (f in c("pathways", "gwas_hits", "variants", "pangenes")) {
    x <- fread(file.path("annotation", paste0(f, ".tsv")), na.strings = NULL,
               colClasses = if (f == "gwas_hits") list(character = c("p_value", "r2", "pubmed")) else NULL)
    if (f == "gwas_hits") x[, `:=`(p_value = as.numeric(p_value), r2 = as.numeric(r2))]
    write_parquet(x, file.path("annotation", paste0(f, ".parquet")))
    cat(sprintf("%s.parquet: %d rows, %d genes\n", f, nrow(x), uniqueN(x$gene_id)))
  }
}
# B73 expression atlas (fetch_maizemine.py) and gene-specific literature (fetch_literature.py)
if (file.exists("annotation/atlas.tsv.gz")) {
  write_parquet(fread("annotation/atlas_samples.tsv", na.strings = NULL), "annotation/atlas_samples.parquet")
  a <- fread(cmd = "gzip -dc annotation/atlas.tsv.gz", colClasses = c(fpkm = "numeric"))   # no R.utils needed
  write_parquet(a, "annotation/atlas.parquet")
  cat(sprintf("atlas.parquet: %d values, %d genes\n", nrow(a), uniqueN(a$gene_id)))
}
if (file.exists("annotation/literature.tsv")) {
  l <- fread("annotation/literature.tsv", colClasses = c(pmid = "character", year = "character"), na.strings = NULL, quote = "")
  write_parquet(l, "annotation/literature.parquet")
  cat(sprintf("literature.parquet: %d gene-paper links, %d genes\n", nrow(l), uniqueN(l$gene_id)))
}
