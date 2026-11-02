#!/usr/bin/env Rscript
# Library QC for Kremling2018 (M3). Strict by design: quality over quantity (rules from MaizExpressN and C,
# adapted to one 3' tag library per line and tissue, plus two identity checks new for D).
#
# Reads metadata/library_qc.tsv (pipeline/library_qc.py on Hazel), metadata/samples.tsv and the gene count
# tables in results/Kremling2018/<batch>/. Writes metadata/library_exclusions.tsv (study_id, batch, sample_id;
# no reasons) and, in private/library_qc/, the reasons and every per-library measure.
# A library is excluded when any of these holds:
#   1. depth: fewer than MIN_MAPPED reads mapped by salmon (3' tag floor of 1 M, approved at the M3 checkpoint
#      with the rest of these rules, Rubén 30 Oct 2026; the impact of 1.5 and 2 M is printed alongside);
#   2. mapping: more than 10 points below what its tissue x flowcell and its line predict (additive fit; C
#      compared with the median of the same study, tissue and genotype, which needs replicates). No absolute
#      floor: 3' reads that run past annotated transcript ends count as decoys, so batch medians are 46-64%;
#   3. library type: not SF (the stranded-forward 3' prep);
#   4. replicates (B73 with B73-10 per tissue, the other repeat libraries; adult leaf within sampling group;
#      log2-CPM r over genes with tissue median CPM >= 20, so 1-2 M-read libraries are not judged on depth):
#      - in a group of >= 3, its mean r with the others is > 0.03 below the tissue's median within-group r and
#        > 0.03 below the others' r among themselves;
#      - in a group of 2, the pair's r is > 0.05 below the tissue's median within-group r (both go);
#   5. tissue carry-over (MaizExpressN rule 5, for every pair of organs: root, shoot, leaf, kernel): another
#      organ's markers make up >= 3x the tissue median share of mRNA, >= 0.5 points above it, >= 3x the median
#      of its replicates (if any), and >= 20% of those markers are >= 3x their tissue median. Markers: the 100
#      genes most expressed in the other organ (median CPM >= 20) that are >= 20x their median in the tissue;
#   6. tissue identity (new): its profile is closer to another tissue's than to its own (r of deviations from
#      the mean of the tissue profiles; its own tissue's profile leaves it out). L3Mid and LMid count as one
#      tissue; day vs night in adult leaf is decided by a day/night score instead (mean deviation over genes
#      that change >= 2^1.5-fold between a line's LMAD and LMAN, in the same direction in >= 90% of lines);
#   7. line identity (new): its genotype signal (log deviations from the tissue median at genes with >= 20
#      expected reads) is not clearly shared with its own line's libraries in other tissues: median similarity
#      < 2 MADs above its similarity to other lines' libraries (a library from another line sits near 0), when the line has >= 2 such reference libraries
#      (passing rules 1-6, >= 1 M mapped reads). B73, B73-10 and B73Htrhm count as one line here, as do W22 and
#      W22R-rstd, and NC296 and NC296A.
# The pilot leaf-3 base libraries (flowcell HHV55BGXX, submission 10344825) are a set of their own for rule 4:
# they are not replicates of the main L3Base libraries.
# Usage (project root): Rscript metadata/build/library_exclusions.R [min_mapped_millions]
suppressMessages({library(data.table); library(matrixStats)})

args <- commandArgs(trailingOnly = TRUE)
MIN_MAPPED <- if (length(args)) as.numeric(args[1]) * 1e6 else 1e6
FLOORS <- c(1e6, 1.5e6, 2e6)                 # impact of each floor, for the checkpoint
MAX_BELOW <- 10
REP_MIN_CPM <- 20; R_MARGIN_3 <- 0.03; R_MARGIN_2 <- 0.05
CARRY_FOLD <- 3; CARRY_EXCESS <- 0.5; CARRY_BREADTH <- 0.2
MARKER_CPM <- 20; MARKER_FOLD <- 20; N_MARKERS <- 100
DN_LFC <- 1.5; DN_CONS <- 0.9
ID_EXPECTED <- 20; ID_Z <- 2; ID_MIN_REF <- 2; ID_REF_MAPPED <- 1e6
FAMILY <- c("B73-10" = "B73", B73Htrhm = "B73", "W22R-rstd" = "W22", NC296A = "NC296")
PILOT_FC <- "HHV55BGXX"
ORGAN <- c(GRoot = "root", GShoot = "shoot", L3Base = "leaf", L3Mid = "leaf", LMid = "leaf", L3Tip = "leaf",
           LMAD = "leaf", LMAN = "leaf", Kern = "kernel")
TISSUES <- names(ORGAN)
stopifnot(file.exists("metadata/samples.tsv"))

q <- fread("metadata/library_qc.tsv")
s <- fread("metadata/samples.tsv")[include == TRUE]
tabs <- lapply(Sys.glob("results/Kremling2018/*/salmon.merged.gene_counts.tsv"), fread)
X <- do.call(cbind, lapply(tabs, function(x) as.matrix(x[, -(1:2)])))
rownames(X) <- tabs[[1]]$gene_id
stopifnot(all(sapply(tabs, function(x) identical(x$gene_id, rownames(X)))))
rm(tabs)
s <- s[sample_id %in% colnames(X)]
X <- X[, s$sample_id]
q <- q[match(s$sample_id, sample_id)]
stopifnot(!anyNA(q$sample_id), nrow(s) == ncol(X))
ids <- s$sample_id
N <- colSums(X)
cpm <- t(t(X) / N) * 1e6
tis <- s$tissue
fam <- ifelse(s$genotype %in% names(FAMILY), FAMILY[s$genotype], s$genotype)
set <- ifelse(tis == "L3Base" & s$flowcell == PILOT_FC, "L3Base_pilot", tis)
grp <- paste(set, ifelse(s$genotype == "B73-10", "B73", s$genotype), ifelse(tis %in% c("LMAD", "LMAN"), s$sampling_group, ""))
tmed_of <- function(ok) sapply(TISSUES, function(t) rowMedians(cpm[, tis == t & ok, drop = FALSE]))
Lc <- log2(cpm + 1)

# ---- precomputed, independent of the depth floor ----
# rule 2: additive fit of % mapped
fit <- lm(percent_mapped ~ tissue:flowcell + genotype, data = cbind(q[, .(percent_mapped)], s[, .(tissue, flowcell, genotype)]))
map_expected <- fitted(fit)
# rule 4: r between libraries of the same tissue
rcor <- lapply(setNames(TISSUES, TISSUES), function(t) {
  k <- which(tis == t)
  g <- rowMedians(cpm[, k]) >= REP_MIN_CPM
  r <- cor(log2(cpm[g, k] + 1)); dimnames(r) <- list(ids[k], ids[k]); r
})
# rule 7: genotype signal and similarity between all libraries
tm_all <- tmed_of(rep(TRUE, length(ids)))
gk <- rowMaxs(tm_all) >= 5
E <- tm_all[gk, tis] * rep(N / 1e6, each = sum(gk))
D <- log2((X[gk, ] + 1) / (E + 1))
D[E < ID_EXPECTED] <- 0
D <- pmin(pmax(D, -5), 5)
D <- sweep(D, 2, colSums(D) / pmax(colSums(D != 0), 1))
D[E < ID_EXPECTED] <- 0
nd <- sqrt(colSums(D^2))
SIM <- crossprod(D) / outer(nd, nd)
rm(D, E)

run <- function(floor) {
  why <- setNames(vector("list", length(ids)), ids)
  add <- function(i, rule, txt) for (k in seq_along(i)) why[[i[k]]] <<- c(why[[i[k]]], paste0("[", rule, "] ", txt[k]))
  fails <- function() lengths(why) > 0
  m <- list()
  # rules 1-3
  i <- which(N < floor); add(i, 1, sprintf("%.2f M mapped reads (< %.1f M)", N[i] / 1e6, floor / 1e6))
  res <- q$percent_mapped - map_expected
  i <- which(res < -MAX_BELOW)
  add(i, 2, sprintf("%.1f%% mapped vs %.1f%% expected for its tissue, flowcell and line", q$percent_mapped[i], map_expected[i]))
  i <- which(q$library_type != "SF"); add(i, 3, sprintf("library type %s vs SF", q$library_type[i]))
  pass13 <- !fails()
  # rule 4
  rep_r <- setNames(rep(NA_real_, length(ids)), ids)
  for (t in TISSUES) {
    R <- rcor[[t]]
    k <- which(tis == t & pass13)
    G <- split(ids[k], grp[k]); G <- G[lengths(G) >= 2]
    within <- unlist(lapply(G, function(g) { r <- R[g, g]; r[upper.tri(r)] }))
    if (length(within) < 5) next
    wmed <- median(within)
    for (g in G) {
      if (length(g) == 2) {
        r <- R[g[1], g[2]]; rep_r[g] <- r
        if (r < wmed - R_MARGIN_2)
          add(g, 4, rep(sprintf("its 2 replicates disagree: r = %.3f vs tissue within-group median %.3f", r, wmed), 2))
      } else for (x in g) {
        o <- setdiff(g, x)
        r_x <- mean(R[x, o]); ro <- R[o, o]; r_o <- mean(ro[upper.tri(ro)])
        rep_r[x] <- r_x
        if (r_x < wmed - R_MARGIN_3 && r_o - r_x > R_MARGIN_3)
          add(x, 4, sprintf("disagrees with its replicates: r = %.3f with them vs %.3f among them (tissue within-group median %.3f)", r_x, r_o, wmed))
      }
    }
  }
  m$rep_r <- rep_r
  pass14 <- !fails()
  # rule 5
  tm <- tmed_of(pass14)
  om <- sapply(unique(ORGAN), function(o) rowMedians(cpm[, ORGAN[tis] == o & pass14, drop = FALSE]))
  carry_fold <- setNames(rep(NA_real_, length(ids)), ids); carry_organ <- setNames(rep(NA_character_, length(ids)), ids)
  for (t in TISSUES) for (o in setdiff(colnames(om), ORGAN[t])) {
    g <- which(om[, o] >= MARKER_CPM & (om[, o] + 1) / (tm[, t] + 1) >= MARKER_FOLD)
    g <- g[order(-om[g, o])][seq_len(min(N_MARKERS, length(g)))]
    k <- which(tis == t & pass14)
    sh <- colSums(cpm[g, k, drop = FALSE]) / 1e4
    med <- median(sh)
    fold <- (sh + 1e-4) / (med + 1e-4)
    breadth <- colMeans((cpm[g, k, drop = FALSE] + 1) / (tm[g, t] + 1) >= 3)
    for (j in seq_along(k)) {
      x <- k[j]
      if (is.na(carry_fold[x]) || fold[j] > carry_fold[x]) { carry_fold[x] <- fold[j]; carry_organ[x] <- o }
      reps <- setdiff(k[grp[k] == grp[x]], x)
      fold_rep <- if (length(reps)) (sh[j] + 1e-4) / (median(sh[match(reps, k)]) + 1e-4) else NA
      if (fold[j] >= CARRY_FOLD && sh[j] - med >= CARRY_EXCESS && breadth[j] >= CARRY_BREADTH && (is.na(fold_rep) || fold_rep >= CARRY_FOLD))
        add(x, 5, sprintf("%s carry-over: %s markers %.2f%% of mRNA vs %.2f%% tissue median (%.1fx; %s; %.0f%% of markers >= 3x)",
                          o, o, sh[j], med, fold[j], if (is.na(fold_rep)) "no replicates" else sprintf("%.1fx its replicates", fold_rep),
                          100 * breadth[j]))
    }
  }
  m$carry_fold <- carry_fold; m$carry_organ <- carry_organ
  # rule 6
  tm13 <- tmed_of(pass13)
  L <- Lc[rowMaxs(tm13) >= 5, ]
  S <- sapply(TISSUES, function(t) rowSums(L[, tis == t & pass13, drop = FALSE]))
  n <- sapply(TISSUES, function(t) sum(tis == t & pass13))
  C0 <- sweep(S, 2, n, "/")
  sc <- t(sapply(seq_along(ids), function(i) {
    C <- C0
    if (pass13[i]) C[, tis[i]] <- (S[, tis[i]] - L[, i]) / (n[tis[i]] - 1)
    g <- rowMeans(C)
    cor(L[, i] - g, C - g)[1, ]
  }))
  colnames(sc) <- TISSUES
  nearest <- TISSUES[max.col(sc)]
  same6 <- function(a, b) a == b | (a %in% c("L3Mid", "LMid") & b %in% c("L3Mid", "LMid")) | (a %in% c("LMAD", "LMAN") & b %in% c("LMAD", "LMAN"))
  own <- sc[cbind(seq_along(ids), match(tis, TISSUES))]
  i <- which(!same6(nearest, tis))
  add(i, 6, sprintf("profile closer to %s (r %.3f) than to %s (r %.3f)", nearest[i], sc[cbind(i, match(nearest[i], TISSUES))], tis[i], own[i]))
  ad <- s[tis %in% c("LMAD", "LMAN") & pass13, .(sample_id, genotype, tissue, sampling_group)]
  pr <- merge(ad[tissue == "LMAD"], ad[tissue == "LMAN"], by = c("genotype", "sampling_group"), suffixes = c(".d", ".n"))
  pr <- pr[!duplicated(pr[, .(genotype, sampling_group)])]
  dn <- Lc[, pr$sample_id.n] - Lc[, pr$sample_id.d]
  lfc <- rowMedians(dn); cons <- rowMeans(sign(dn) == sign(lfc))
  dg <- rowMedians(cpm[, c(pr$sample_id.n, pr$sample_id.d)]) >= 10 & abs(lfc) >= DN_LFC & cons >= DN_CONS
  mid <- (rowMedians(Lc[dg, tis == "LMAD" & pass13]) + rowMedians(Lc[dg, tis == "LMAN" & pass13])) / 2
  night <- setNames(rep(NA_real_, length(ids)), ids)
  a <- which(tis %in% c("LMAD", "LMAN"))
  night[a] <- as.vector(crossprod(sign(lfc[dg]), Lc[dg, a] - mid)) / sum(dg)
  md <- median(night[tis == "LMAD" & pass13]); mn <- median(night[tis == "LMAN" & pass13])
  i <- which((tis == "LMAD" & night > 0) | (tis == "LMAN" & night < 0))
  add(i, 6, sprintf("day/night score %.2f (LMAD median %.2f, LMAN median %.2f)", night[i], md, mn))
  m$own_r <- own; m$nearest <- nearest; m$nearest_r <- sc[cbind(seq_along(ids), match(nearest, TISSUES))]
  m$night <- night; m$n_diurnal <- sum(dg); m$n_pairs <- nrow(pr)
  # rule 7
  ref <- !fails() & N >= ID_REF_MAPPED
  idm <- rbindlist(lapply(seq_along(ids), function(i) {
    o <- ref & tis != tis[i]; o[i] <- FALSE
    bg <- SIM[i, o & fam != fam[i]]
    mm <- tapply(SIM[i, o & fam != fam[i]], fam[o & fam != fam[i]], median)
    b <- median(bg); d <- mad(bg)
    own_i <- if (any(o & fam == fam[i])) median(SIM[i, o & fam == fam[i]]) else NA_real_
    list(id_n_ref = sum(o & fam == fam[i]), id_own = own_i, id_z_own = (own_i - b) / d,
         id_alt_line = names(mm)[which.max(mm)], id_z_alt = (max(mm) - b) / d)
  }))
  i <- which(idm$id_n_ref >= ID_MIN_REF & idm$id_z_own < ID_Z)
  add(i, 7, sprintf("not clearly its line: similarity to its line's %d libraries in other tissues %.1f MADs above other lines (closest other line %s, %.1f)",
                    idm$id_n_ref[i], idm$id_z_own[i], idm$id_alt_line[i], idm$id_z_alt[i]))
  m$id <- idm
  list(why = why, m = m)
}

summ <- function(why) {
  ex <- lengths(why) > 0
  d <- data.table(tissue = tis, set = set, genotype = s$genotype, ex = ex)
  per <- d[, .(libraries = .N, excluded = sum(ex), kept = sum(!ex),
               lines_kept = uniqueN(genotype[!ex & !(genotype %in% c("B73-10"))])), .(tissue)]
  list(ex = ex, per = per)
}
rules_of <- function(why) sapply(why, function(w) paste(unique(sub("^\\[(\\d)\\].*", "\\1", w)), collapse = ";"))

cat("Library QC, Kremling2018:", length(ids), "libraries\n")
fl <- sort(unique(c(FLOORS, MIN_MAPPED)))
runs <- setNames(lapply(fl, run), fl)
impact <- lapply(fl, function(f) cbind(floor_M = f / 1e6, summ(runs[[as.character(f)]]$why)$per))
why <- runs[[as.character(MIN_MAPPED)]]$why; m <- runs[[as.character(MIN_MAPPED)]]$m
rl <- rules_of(why)
ex <- lengths(why) > 0
cat(sprintf("\nfloor %.1f M: %d of %d libraries excluded\n", MIN_MAPPED / 1e6, sum(ex), length(ids)))
cat("libraries flagged by each rule, by floor (in brackets: by that rule alone):\n")
print(rbindlist(lapply(names(runs), function(f) {
  r <- rules_of(runs[[f]]$why)
  as.data.table(c(list(floor_M = as.numeric(f) / 1e6), setNames(lapply(1:7, function(k) {
    has <- grepl(as.character(k), r); sprintf("%d (%d)", sum(has), sum(has & r == as.character(k)))
  }), paste0("rule", 1:7)), list(excluded = sum(lengths(runs[[f]]$why) > 0))))
})))
cat(sprintf("day/night score from %d genes over %d LMAD/LMAN pairs\n", m$n_diurnal, m$n_pairs))
cat("\nper tissue at each floor (lines_kept counts genotypes with >= 1 library kept; B73-10 counted with B73):\n")
print(dcast(rbindlist(impact), tissue ~ floor_M, value.var = c("excluded", "lines_kept")))
cat("\nlines with every main tissue kept (GRoot, GShoot, L3Base, L3Tip, LMAD or LMAN, Kern), by floor:\n")
for (f in names(runs)) {
  w <- runs[[f]]$why
  k <- data.table(genotype = s$genotype, tissue = tis, set = set, ok = lengths(w) == 0)[ok == TRUE & set != "L3Base_pilot"]
  a <- k[, .(all = all(c("GRoot", "GShoot", "L3Base", "L3Tip", "Kern") %in% tissue) & any(c("LMAD", "LMAN") %in% tissue)), genotype]
  cat(sprintf("  %.1f M: %d lines\n", as.numeric(f) / 1e6, sum(a$all)))
}
cat("\nlines with both libraries kept for the proposed paired contrasts (adult leaf: same sampling group), by floor:\n")
pairs_of <- function(w) {
  k <- data.table(genotype = s$genotype, tissue = tis, set = set, sg = s$sampling_group, ok = lengths(w) == 0)[ok == TRUE & set != "L3Base_pilot"]
  has <- function(a, b) length(intersect(k[tissue == a, genotype], k[tissue == b, genotype]))
  dn <- merge(unique(k[tissue == "LMAD", .(genotype, sg)]), unique(k[tissue == "LMAN", .(genotype, sg)]))
  data.table(L3Tip_vs_L3Base = has("L3Base", "L3Tip"), LMAN_vs_LMAD = uniqueN(dn$genotype),
             LMAD_vs_L3Tip = has("L3Tip", "LMAD"), GShoot_vs_GRoot = has("GRoot", "GShoot"))
}
print(rbindlist(lapply(names(runs), function(f) cbind(floor_M = as.numeric(f) / 1e6, pairs_of(runs[[f]]$why)))))

# B73 controls: agreement within and across flowcells (libraries kept)
cat("\nB73 controls (B73 and B73-10, libraries kept): median r within / across flowcells; different lines on one flowcell\n")
b73 <- s$genotype %in% c("B73", "B73-10")
bt <- rbindlist(lapply(c("GRoot", "GShoot", "L3Base", "L3Tip", "LMAD", "LMAN", "Kern"), function(t) {
  k <- which(tis == t & !ex & set != "L3Base_pilot"); r <- rcor[[t]][ids[k], ids[k]]
  ut <- which(upper.tri(r), arr.ind = TRUE); a <- k[ut[, 1]]; b <- k[ut[, 2]]
  sf <- s$flowcell[a] == s$flowcell[b]; bb <- b73[a] & b73[b]; dl <- s$genotype[a] != s$genotype[b] & !bb
  sg <- s$sampling_group[a] == s$sampling_group[b] | !(t %in% c("LMAD", "LMAN"))
  data.table(tissue = t, b73 = sum(b73[k]), flowcells = uniqueN(s$flowcell[k][b73[k]]),
             within_fc = round(median(r[ut][bb & sf & sg]), 3), across_fc = round(median(r[ut][bb & !sf & sg]), 3),
             across_group = if (t %in% c("LMAD", "LMAN")) round(median(r[ut][bb & !sg]), 3) else NA_real_,
             lines_same_fc = round(median(r[ut][dl & sf]), 3))
}))
print(bt)

# outputs
dir.create("private/library_qc", recursive = TRUE, showWarnings = FALSE)
fwrite(data.table(study_id = "Kremling2018", batch = q$study_id, sample_id = ids)[ex], "metadata/library_exclusions.tsv", sep = "\t")
fwrite(data.table(batch = q$study_id, sample_id = ids, tissue = tis, genotype = s$genotype, rules = rl,
                  reason = sapply(why, paste, collapse = "; "))[ex], "private/library_qc/library_exclusions_reasons.tsv", sep = "\t")
fwrite(cbind(data.table(batch = q$study_id, sample_id = ids, tissue = tis, set = set, genotype = s$genotype, flowcell = s$flowcell,
                        plate = s$plate, sampling_group = s$sampling_group, reads_mapped = q$reads_mapped, percent_mapped = q$percent_mapped,
                        mapping_expected = round(map_expected, 1), replicate_r = round(m$rep_r, 3),
                        carry_max_fold = round(m$carry_fold, 2), carry_organ = m$carry_organ,
                        own_tissue_r = round(m$own_r, 3), nearest_tissue = m$nearest, nearest_r = round(m$nearest_r, 3),
                        night_score = round(m$night, 3)),
             m$id[, .(id_n_ref, id_own = round(id_own, 3), id_z_own = round(id_z_own, 2), id_alt_line, id_z_alt = round(id_z_alt, 2))],
             data.table(excluded = ex, rules = rl)),
       "private/library_qc/library_qc_measures.tsv", sep = "\t")
cat("\nwrote metadata/library_exclusions.tsv, private/library_qc/library_exclusions_reasons.tsv, private/library_qc/library_qc_measures.tsv\n")
