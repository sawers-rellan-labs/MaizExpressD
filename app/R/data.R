# Data access for the MaizExpressD app. Shiny sources app/R/*.R before app.R (alphabetically),
# so the packages these helpers need are loaded here.
suppressPackageStartupMessages({ library(data.table); library(arrow); library(plotly) })

# The data folder is app/data by default; set MAIZEXPRESSD_DATA to point elsewhere (built by app/prepare_data.R).
DATA_DIR <- Sys.getenv("MAIZEXPRESSD_DATA", "data")
rp <- function(f) as.data.table(arrow::read_parquet(file.path(DATA_DIR, f)))

GENES       <- rp("genes.parquet")
ALIASES     <- rp("gene_aliases.parquet")
GO_TERMS    <- rp("go_terms.parquet")
TISSUE_TAB  <- rp("tissues.parquet")
LINES       <- rp("lines.parquet")
LINE_TISSUE <- rp("line_tissue.parquet")
SAMPLES     <- rp("samples.parquet")
CONTRASTS   <- rp("contrasts.parquet")
DE          <- rp("de.parquet")
EXPR_T      <- rp("expr_tissue.parquet")[, tissue := as.character(tissue)]
SPEC        <- rp("specificity.parquet")
LVAR        <- rp("line_variation.parquet")[, tissue := as.character(tissue)]
UNSTABLE    <- rp("b73_unstable.parquet")
LEI_TS      <- rp("lei_tissue_specific.parquet")
LEI_EQTL    <- rp("lei_eqtls.parquet")
LEI_TWAS    <- rp("lei_twas.parquet")
LEI_COLOC   <- rp("lei_coloc.parquet")
LEI_CUBIC   <- rp("lei_cubic.parquet")
LEI_CAND    <- rp("lei_candidates.parquet")
LEI_TRAITS  <- rp("lei_traits.parquet")
PATHWAYS    <- rp("pathways.parquet")
GWAS_HITS   <- rp("gwas_hits.parquet")
VARIANTS    <- rp("variants.parquet")
PANGENES    <- rp("pangenes.parquet")
LITERATURE  <- rp("literature.parquet")
ATLAS_SAMPLES <- rp("atlas_samples.parquet")
MAIZEMINE   <- tryCatch(jsonlite::fromJSON(file.path(DATA_DIR, "maizemine_release.json")), error = function(e) list())
BUILT       <- tryCatch(readLines(file.path(DATA_DIR, "BUILT"))[1], error = function(e) "")
setkey(GENES, gene_id); setkey(ALIASES, alias_lower)
for (x in list(DE, EXPR_T, SPEC, LVAR, UNSTABLE, LEI_TS, LEI_EQTL, LEI_TWAS, LEI_COLOC, LEI_CUBIC, LEI_CAND,
               PATHWAYS, GWAS_HITS, VARIANTS, PANGENES, LITERATURE)) setkey(x, gene_id)
go_name <- setNames(GO_TERMS$go_name, GO_TERMS$go_id)

TISSUES <- TISSUE_TAB$tissue                                   # display order
LEI_TISSUES <- c("GRoot", "GShoot", "L3Base", "L3Tip", "LMAD", "LMAN", "Kern")
TISSUE_LABEL <- setNames(TISSUE_TAB$tissue_label, TISSUES)
# root brown, shoot blue, the leaf-3 gradient in greens (base light, tip dark), adult leaf by day amber and by
# night indigo, kernel pink (Okabe-Ito based)
TISSUE_COL <- c(GRoot = "#7A5C3A", GShoot = "#3C8DC4", L3Base = "#9BC97E", L3Mid = "#5DA35A", LMid = "#7FAE9A",
                L3Tip = "#1E6B34", LMAD = "#D99000", LMAN = "#3D3A8C", Kern = "#C46A9E")
SUBPOPS <- c("Stiff stalk", "Non-stiff stalk", "Tropical/subtropical", "Mixed", "Popcorn", "Sweet corn", "Landrace", "Not reported")
SUBPOP_COL <- c(`Stiff stalk` = "#0072B2", `Non-stiff stalk` = "#E69F00", `Tropical/subtropical` = "#009E73", Mixed = "#9AA3AB",
                Popcorn = "#D55E00", `Sweet corn` = "#CC79A7", Landrace = "#56B4E9", `Not reported` = "#C9CED3")
CLASSES <- c("leaf_maturation", "leaf_stage", "day_night", "organ")
CLASS_LABEL <- c(leaf_maturation = "Leaf 3 gradient", leaf_stage = "Adult vs seedling leaf", day_night = "Night vs day",
                 organ = "Shoot vs root")
CLASS_COL <- c(leaf_maturation = "#1E6B34", leaf_stage = "#D99000", day_night = "#3D3A8C", organ = "#7A5C3A")
LEI_TISSUE_CODE <- c(RT = "GRoot", SH = "GShoot", LB = "L3Base", LT = "L3Tip", LD = "LMAD", LN = "LMAN", KN = "Kern")

CONTRASTS[, `:=`(class_label = CLASS_LABEL[contrast_class],
                 label = paste0(numerator_tissue, " vs ", denominator_tissue,
                                fifelse(design == "panel", fifelse(nzchar(sampling_group), paste0(" · group ", sampling_group), ""), " · B73")),
                 n_text = fifelse(design == "panel", sprintf("%d lines", n_pairs), sprintf("%d vs %d B73 libraries", n_num, n_den)))]
CONTRASTS <- CONTRASTS[order(match(contrast_class, CLASSES), design != "panel", sampling_group != "", -n_pairs)]   # main panel fit first

bucket <- function(g) substr(g, 10, 12)
read_bucket <- function(dir, gid) {   # gid: one gene, or several of one bucket
  f <- file.path(DATA_DIR, dir, paste0(bucket(gid[1]), ".parquet"))
  if (!file.exists(f)) return(data.table())
  as.data.table(arrow::read_parquet(f))[gene_id %in% gid]
}

# CPM of every line in every tissue (the per-line tables keep CPM > 0 only), with the line's subpopulation
get_expr <- function(gid) {
  x <- read_bucket("expr", gid)
  d <- LINE_TISSUE[, .(tissue, genotype, n_libraries)]
  d[x, cpm := i.cpm, on = c("tissue", "genotype")]
  d[is.na(cpm), cpm := 0]
  d[LINES, subpopulation := i.subpopulation, on = "genotype"]
  d[]
}
get_atlas <- function(gid) read_bucket("atlas", gid)
MAIZEMINE_URL <- "https://maizemine.rnet.missouri.edu/maizemine/portal.do"
maizemine_link <- function(id, class) paste0(MAIZEMINE_URL, "?externalids=", utils::URLencode(id, reserved = TRUE), "&class=", class)

# strongest co-expression partners (analysis/coexpr_partners.R); partner IDs are stored as their number
get_partners <- function(gid) {
  x <- read_bucket("partners", gid)
  if (!nrow(x)) return(x)
  x[, partner_id := sprintf("Zm00001eb%06d", partner)][, partner := NULL][]
}

resolve_tokens <- function(tokens) {
  tokens <- unique(trimws(tokens[nzchar(trimws(tokens))]))
  hits <- lapply(tokens, function(t) unique(ALIASES[.(tolower(t)), gene_id, nomatch = 0]))
  list(genes = unique(unlist(hits)), missing = tokens[lengths(hits) == 0])
}
gene_label <- function(ids) {
  g <- GENES[.(ids)]
  ifelse(!is.na(g$symbol) & nzchar(g$symbol), paste0(g$symbol, " (", g$gene_id, ")"), ids)
}
short_name <- function(g) { x <- GENES[.(g)]; if (nrow(x) && !is.na(x$symbol) && nzchar(x$symbol)) x$symbol else g }

# v5 chromosome starts for genome-wide plots (from the gene models; scaffolds left out)
CHROMS <- GENES[grepl("^chr[0-9]+$", chrom), .(len = max(end)), by = chrom][order(as.integer(sub("chr", "", chrom)))]
CHROMS[, offset := cumsum(as.numeric(c(0, head(len, -1))))]
