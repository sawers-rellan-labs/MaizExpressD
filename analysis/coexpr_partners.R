#!/usr/bin/env Rscript
# M6: strongest co-expression partners of every gene across lines, within each tissue (HANDOFF M6), adapted from
# MaizExpressC's app/prepare_partners.R. Run after analysis/expression_tissue.R and analysis/line_variation.R:
#   Rscript analysis/coexpr_partners.R
#
# Unit = tissue (the seven of Lei et al. 2025; L3Mid has 26 lines and LMid 10, too few for ranked lists). Pearson r of
# log2(CPM + 1) of the line means (results/expression_line.parquet; zeros counted), one value per line, B73 one of the
# lines. Two views, as C's "across" and "within conditions":
#   lines            across lines; adult leaf (LMAD, LMAN) centred within sampling group first, since the two groups
#                    were sampled on different dates (M1);
#   within_flowcell  also centred within flowcell (partial r given flowcell): flowcell explains up to 61-81% of a top
#                    principal component of the line values in L3Tip, LMAD and LMAN, partly together with
#                    subpopulation, so this view takes the sequencing batch out at the cost of some genetic signal.
# A line whose libraries span groups or flowcells is its own level. df = lines - levels - 1; a unit counts for a pair
# only when both genes have mean CPM >= 1 there, neither is B73-unstable there (results/b73_unstable.tsv), and it
# has >= 3 df; Fisher z with SE 1/sqrt(df - 1). Scopes: each tissue alone, and All (DerSimonian-Laird pooling over
# the seven tissues; a gene gets a list when expressed in at least 4, a pair is ranked when at least 4 have both).
# Pairs are ranked by r (pooled r for All); top K per gene, scope, view and direction.
#   results/partners/<bucket>.parquet   gene_id, scope, view, direction (+/-), rank, partner (number of its ID),
#                                       units, r, lo, hi, i2 (bucket = characters 10-12 of the gene ID, as in C)
#   results/coexpr_strongest_pairs.tsv  the strongest pairs of the All scope, both views

source("analysis/common.R")
suppressPackageStartupMessages(library(parallel))
MIN_CPM <- 1; MIN_DF <- 3; BLOCK <- 3000; K <- 25; CORES <- 5
SCOPES <- strsplit(Sys.getenv("PARTNER_SCOPES", paste(c(LEI_TISSUES, "All"), collapse = ",")), ",")[[1]]   # subset for testing
bucket <- function(g) substr(g, 10, 12)

S <- fread("results/expression_lines.tsv")[tissue %in% LEI_TISSUES]
lev <- function(col) sapply(strsplit(S$sample_ids, ";", fixed = TRUE), function(ids) {
  v <- unique(libs[match(ids, sample_id), get(col)]); if (length(v) == 1) v else "mixed" })
S[, `:=`(unit = tissue, col = paste(tissue, genotype, sep = "|"), sgroup = lev("sampling_group"), flowcell = lev("flowcell"))]
S[, `:=`(g_lines = ifelse(tissue %in% c("LMAD", "LMAN"), sgroup, "all"))]
S[, g_flowcell := paste(g_lines, flowcell, sep = "|")]
units <- S[, .(n = .N, lv_lines = uniqueN(g_lines), lv_flowcell = uniqueN(g_flowcell)), by = unit]

d <- as.data.table(read_parquet("results/expression_line.parquet"))[tissue %in% LEI_TISSUES]
m <- d[, .(s = sum(cpm)), by = .(gene_id, unit = as.character(tissue))][units, on = "unit"][, mean_cpm := s / n]
uns <- fread("results/b73_unstable.tsv")[, .(gene_id, unit = tissue)]
m <- m[!uns, on = c("gene_id", "unit")]
G_all <- sort(unique(m[mean_cpm >= MIN_CPM, gene_id]))
d <- d[gene_id %in% G_all]
CP <- matrix(0, length(G_all), nrow(S))
CP[cbind(match(d$gene_id, G_all), match(paste(d$tissue, d$genotype, sep = "|"), S$col))] <- d$cpm
rm(d); invisible(gc())
L <- log2(CP + 1)
EXPR <- matrix(FALSE, length(G_all), nrow(units), dimnames = list(NULL, units$unit))   # mean CPM >= 1, not unstable
e <- m[mean_cpm >= MIN_CPM & gene_id %in% G_all]
EXPR[cbind(match(e$gene_id, G_all), match(e$unit, units$unit))] <- TRUE

# per unit: rows centred within the view's levels and scaled to unit norm, so r = Z_a . Z_b
prep <- function(G, us, view) {
  gi <- match(G, G_all)
  lapply(us, function(uu) {
    j <- which(S$unit == uu); grp <- S[[paste0("g_", view)]][j]
    df <- length(j) - uniqueN(grp) - 1
    H <- outer(grp, unique(grp), "==") * 1
    X <- L[gi, j, drop = FALSE]; X <- X - X %*% (H %*% (t(H) / colSums(H)))   # minus level means
    ss <- sqrt(rowSums(X^2))
    ok <- EXPR[gi, uu] & ss > 1e-12
    X[!ok, ] <- 0; X[ok, ] <- X[ok, ] / ss[ok]
    list(Z = X, ok = ok, w = df - 1)
  })
}

# pooled r for genes ia x genes ib (MaizExpressC's block_pair), then each gene's top K partners in both directions
block_pair <- function(U, ia, ib, kmin) {
  same <- identical(ia, ib)
  zfun <- function(u) atanh(pmin(pmax(tcrossprod(u$Z[ia, , drop = FALSE], u$Z[ib, , drop = FALSE]), -0.9999), 0.9999))
  mk <- function(u) outer(u$ok[ia], u$ok[ib])
  k <- SW <- SW2 <- SWZ <- SWZ2 <- 0
  for (u in U) { M <- mk(u); z <- zfun(u) * M
    k <- k + M; SW <- SW + u$w * M; SW2 <- SW2 + u$w^2 * M; SWZ <- SWZ + u$w * z; SWZ2 <- SWZ2 + u$w * z^2 }
  Q <- SWZ2 - SWZ^2 / SW
  tau2 <- pmax(0, (Q - (k - 1)) / (SW - SW2 / SW)); tau2[k < 2] <- 0
  A <- B <- 0
  for (u in U) { M <- mk(u); ws <- M / (1 / u$w + tau2); A <- A + ws * zfun(u); B <- B + ws }
  zr <- A / B; s <- sqrt(1 / B)
  r <- tanh(zr); r[k < kmin] <- NA
  if (same) diag(r) <- NA
  lo <- tanh(zr - 1.96 * s); hi <- tanh(zr + 1.96 * s)
  i2 <- ifelse(Q > 0, pmax(0, (Q - (k - 1)) / Q), 0)
  tops <- function(R, gi, pj, tr) {
    rbindlist(lapply(c(TRUE, FALSE), function(dec) {
      idx <- lapply(seq_len(nrow(R)), function(i) head(order(R[i, ], decreasing = dec, na.last = NA), K))
      a <- rep(seq_along(idx), lengths(idx)); b <- unlist(idx)
      ij <- if (tr) cbind(b, a) else cbind(a, b)
      data.table(g = gi[a], p = pj[b], direction = if (dec) "+" else "-", units = as.integer(k[ij]),
                 r = r[ij], lo = lo[ij], hi = hi[ij], i2 = i2[ij])
    }))
  }
  rbind(tops(r, ia, ib, FALSE), if (!same) tops(t(r), ib, ia, TRUE))
}

res <- list()
for (sc in SCOPES) for (view in c("lines", "flowcell")) {
  t0 <- Sys.time()
  uu <- units[(sc == "All" | unit == sc) & n - get(paste0("lv_", view)) - 1 >= MIN_DF, unit]
  G <- m[unit %in% uu & mean_cpm >= MIN_CPM, .N, by = gene_id][N >= length(uu) / 2, sort(gene_id)]
  kmin <- if (sc == "All") as.integer(ceiling(length(uu) / 2)) else 1L
  U <- prep(G, uu, view)
  blocks <- split(seq_along(G), ceiling(seq_along(G) / BLOCK))
  bp <- CJ(a = seq_along(blocks), b = seq_along(blocks))[a <= b]
  x <- mclapply(seq_len(nrow(bp)), function(i) block_pair(U, blocks[[bp$a[i]]], blocks[[bp$b[i]]], kmin),
                mc.cores = CORES, mc.preschedule = FALSE)
  bad <- vapply(x, inherits, TRUE, "try-error"); if (any(bad)) stop(x[[which(bad)[1]]])
  x <- rbindlist(x)
  x[, rank := frank(if (direction[1] == "+") -r else r, ties.method = "first"), by = .(g, direction)]
  x <- x[rank <= K]
  vname <- if (view == "lines") "lines" else "within_flowcell"
  x[, `:=`(gene_id = G[g], partner_id = G[p], scope = sc, view = vname)]
  res[[length(res) + 1]] <- x[, .(gene_id, scope, view, direction, rank = as.integer(rank), partner_id, units,
                                  r = round(r, 3), lo = round(lo, 3), hi = round(hi, 3), i2 = round(i2, 2))]
  message(sprintf("%-6s %-15s %d units, %5d genes, %.0f s", sc, vname, length(uu), length(G),
                  as.numeric(Sys.time() - t0, units = "secs")))
  rm(U, x); invisible(gc())
}
P <- rbindlist(res)

dir.create("results/partners", showWarnings = FALSE)
unlink(Sys.glob("results/partners/*.parquet"))
P[, b := bucket(gene_id)]
stopifnot(all(grepl("^Zm00001eb[0-9]{6}$", P$partner_id)))   # stored as their number, as in C (get_partners rebuilds)
for (bb in unique(P$b))
  write_parquet(P[b == bb, .(gene_id, scope, view, direction, rank, partner = as.integer(substr(partner_id, 10, 15)),
                             units, r, lo, hi, i2)], file.path("results/partners", paste0(bb, ".parquet")), compression = "zstd")
P[, b := NULL]
message(sprintf("partners: %d rows for %d genes in %d buckets", nrow(P), uniqueN(P$gene_id), uniqueN(bucket(P$gene_id))))

if (!"All" %in% SCOPES) quit(save = "no")
A <- as.data.table(read_parquet("annotation/genes.parquet"))[, .(gene_id, symbol, description, chrom, start, end)]
top <- P[scope == "All"][, `:=`(g1 = pmin(gene_id, partner_id), g2 = pmax(gene_id, partner_id))]
top <- unique(top, by = c("view", "g1", "g2"))[, `:=`(gene_id = g1, partner_id = g2)]
top <- rbindlist(lapply(split(top, top$view), function(d) rbind(d[order(-r)][1:1000], d[order(r)][1:200])))
a <- A[match(top$gene_id, gene_id)]; p <- A[match(top$partner_id, gene_id)]
top[, `:=`(symbol_a = a$symbol, symbol_b = p$symbol, description_a = a$description, description_b = p$description,
           distance_kb = ifelse(a$chrom == p$chrom, round(pmax(0, pmax(a$start, p$start) - pmin(a$end, p$end)) / 1e3, 1), NA_real_))]
fwrite(top[, .(view, gene_a = gene_id, symbol_a, gene_b = partner_id, symbol_b, tissues = units, pooled_r = r, lo, hi, i2,
               distance_kb, description_a, description_b)], "results/coexpr_strongest_pairs.tsv", sep = "\t")
message("strongest pairs written to results/coexpr_strongest_pairs.tsv")
