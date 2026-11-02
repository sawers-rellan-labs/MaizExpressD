#!/usr/bin/env Rscript
# M7: build the app's data folder from the project outputs (adapted from MaizExpressC's app/prepare_data.R).
#   Rscript app/prepare_data.R            (from the project root; writes app/data/)
#
# Small tables are written whole and loaded at start-up. Per-gene tables of per-line expression and co-expression
# partners are split into ~440 files by gene-ID bucket (the 3 digits after "Zm00001eb"), so a lookup reads one file.
#   genes, gene_aliases, go_terms, pathways, gwas_hits, variants, pangenes, literature .parquet, maizemine_release.json,
#     atlas/<bucket>.parquet, atlas_samples.parquet                annotation (M5, from MaizExpressC)
#   tissues.parquet       the nine tissues: labels, Lei code, stage, lines and libraries in use
#   lines.parquet         lines with libraries in use: Lei name, subpopulation, origin, pedigree, adult-leaf group, tissues
#   line_tissue.parquet   one row per line and tissue: libraries, sampling group and flowcell levels (co-expression views)
#   samples.parquet       libraries in use: tissue, line, sampling group, flowcell, mapped reads and rate
#   contrasts.parquet     the 12 tissue contrasts with their model, line pairs or libraries, and 2-fold DEG counts
#   de.parquet            results/de.parquet (12 contrasts)
#   expr/<bucket>.parquet per-line CPM (results/expression_line.parquet; rows with CPM > 0)
#   expr_tissue.parquet, specificity.parquet, line_variation.parquet, b73_unstable.parquet   (M4, M5)
#   lei_*.parquet         Lei et al. 2025 on v5 genes (annotation/build_lei.R)
#   partners/<bucket>.parquet   co-expression partners (analysis/coexpr_partners.R)
# Libraries left out by QC are not listed, and no reasons are shipped (Methods gives the rules only).

suppressPackageStartupMessages({ library(data.table); library(arrow) })
out <- "app/data"
bucket <- function(g) substr(g, 10, 12)
for (d in c("expr", "partners", "atlas")) {
  dir.create(file.path(out, d), recursive = TRUE, showWarnings = FALSE)
  unlink(Sys.glob(file.path(out, d, "*.parquet")))      # no stale buckets from an earlier build
}
cp <- function(from, to = basename(from)) stopifnot(file.copy(from, file.path(out, to), overwrite = TRUE))

# annotation
for (f in c("genes", "gene_aliases", "go_terms", "pathways", "gwas_hits", "variants", "pangenes", "literature", "atlas_samples"))
  cp(file.path("annotation", paste0(f, ".parquet")))
cp("annotation/source/maizemine/release.json", "maizemine_release.json")
a <- as.data.table(read_parquet("annotation/atlas.parquet"))[, b := bucket(gene_id)]
for (bb in unique(a$b)) write_parquet(a[b == bb][, b := NULL], file.path(out, "atlas", paste0(bb, ".parquet")))
rm(a)

# libraries, lines and tissues in use (as analysis/common.R)
TISSUES <- c("GRoot", "GShoot", "L3Base", "L3Mid", "LMid", "L3Tip", "LMAD", "LMAN", "Kern")
smp <- fread("metadata/samples.tsv", colClasses = "character", na.strings = NULL)
excl <- fread("metadata/library_exclusions.tsv", colClasses = "character")
qc <- fread("metadata/library_qc.tsv", colClasses = "character")
lt <- fread("results/expression_lines.tsv", colClasses = "character")
used <- unlist(strsplit(lt$sample_ids, ";", fixed = TRUE))
s <- smp[sample_id %in% used]
stopifnot(nrow(s) == length(used), !any(s$sample_id %in% excl$sample_id))
s <- qc[, .(sample_id, reads_mapped = as.numeric(reads_mapped), percent_mapped = as.numeric(percent_mapped))][s, on = "sample_id"]
write_parquet(s[, .(sample_id, tissue, genotype, sampling_group, flowcell, plate, replicate, reads_mapped, percent_mapped)],
              file.path(out, "samples.parquet"))

lev <- function(ids, col) { v <- unique(s[match(ids, sample_id), get(col)]); if (length(v) == 1) v else "mixed" }
lt[, n_libraries := as.integer(n_libraries)]
lt[, `:=`(sgroup = lev(strsplit(sample_ids, ";")[[1]], "sampling_group"), flowcell = lev(strsplit(sample_ids, ";")[[1]], "flowcell")),
   by = .(tissue, genotype)]
lt[, g_lines := ifelse(tissue %in% c("LMAD", "LMAN"), sgroup, "all")][, g_flowcell := paste(g_lines, flowcell, sep = "|")]
write_parquet(lt, file.path(out, "line_tissue.parquet"))

SUBPOP <- c(ss = "Stiff stalk", nss = "Non-stiff stalk", ts = "Tropical/subtropical", mixed = "Mixed", popcorn = "Popcorn",
            sweet = "Sweet corn", landrace = "Landrace", NR = "Not reported")
ln <- fread("metadata/lines.tsv", colClasses = "character", na.strings = NULL)[genotype %in% lt$genotype]
ln[, `:=`(subpopulation = fifelse(genotype_class == "landrace", "Landrace", SUBPOP[subpopulation]),   # MR landraces: not in Lei
          n_libraries = NULL, tissues_with_data = NULL)]
ln <- lt[, .(n_libraries = sum(n_libraries), n_tissues = .N, tissues = paste(TISSUES[TISSUES %in% tissue], collapse = ", ")),
         by = genotype][ln, on = "genotype"]
write_parquet(ln[, .(genotype, genotype_lei, genotype_class, subpopulation, origin, pedigree, adult_leaf_group, n_libraries,
                     n_tissues, tissues, note)], file.path(out, "lines.parquet"))

ts <- fread("metadata/build/tissues.tsv", colClasses = "character", na.strings = NULL)
ts <- lt[, .(n_lines = .N, n_libraries = sum(n_libraries)), by = tissue][ts, on = "tissue"]
ts[, order := match(tissue, TISSUES)]
# LMid: confirmed by expression as the middle of leaf 3 (M3); public text without the curation notes
ts[tissue == "LMid", `:=`(tissue_label = "Leaf 3 middle (pilot)", tissue_detail = "middle of leaf 3, pilot libraries (display only)",
                          stage = "seedling, leaf 3", stage_class = "seedling")]
write_parquet(ts[order(order), .(tissue, tissue_label, lei_code, tissue_group, tissue_detail, stage, stage_class, order, n_lines,
                                 n_libraries)], file.path(out, "tissues.parquet"))

# contrasts and DE
cons <- fread("metadata/contrasts.tsv", colClasses = "character", na.strings = NULL)
summ <- fread("results/de_summary.tsv", colClasses = "character", na.strings = NULL)
cm <- summ[, .(contrast_id, model, n_pairs = as.integer(n_pairs), n_num = as.integer(n_numerator), n_den = as.integer(n_denominator),
               genes_tested = as.integer(genes_tested), n_up = as.integer(n_up), n_down = as.integer(n_down))][cons, on = "contrast_id"]
write_parquet(cm[, .(contrast_id, contrast_class, design, numerator_tissue, denominator_tissue, sampling_group, replicated, model,
                     n_pairs, n_num, n_den, genes_tested, n_up, n_down, sign_note, notes, tags)], file.path(out, "contrasts.parquet"))
cp("results/de.parquet")

# expression per line and tissue, by bucket; tissue summaries, specificity, variation
e <- as.data.table(read_parquet("results/expression_line.parquet"))
e[, `:=`(tissue = as.character(tissue), b = bucket(gene_id))]
for (bb in unique(e$b)) write_parquet(e[b == bb, .(gene_id, tissue, genotype, cpm)], file.path(out, "expr", paste0(bb, ".parquet")))
message(sprintf("expr: %d rows in %d buckets", nrow(e), uniqueN(e$b)))
rm(e); invisible(gc())
cp("results/expression_tissue.parquet", "expr_tissue.parquet")
cp("results/tissue_specificity.parquet", "specificity.parquet")
cp("results/line_variation.parquet")
write_parquet(fread("results/b73_unstable.tsv"), file.path(out, "b73_unstable.parquet"))

# Lei et al. 2025 and co-expression partners
for (f in Sys.glob("annotation/lei_*.parquet")) cp(f)
stopifnot(all(file.copy(Sys.glob("results/partners/*.parquet"), file.path(out, "partners"), overwrite = TRUE)))

writeLines(format(Sys.time(), "%d %b %Y %H:%M"), file.path(out, "BUILT"))
message(sprintf("app data written to %s: %d libraries, %d lines, %d contrasts", out, nrow(s), nrow(ln), nrow(cm)))
