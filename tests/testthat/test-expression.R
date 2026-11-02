# Expression per line and tissue (analysis/expression_tissue.R), tissue specificity (analysis/tissue_specificity.R)
# and variation among lines (analysis/line_variation.R).

test_that("tissue summaries are taken over lines, zeros included, B73 counted once", {
  e <- app_env()
  x <- e$get_expr(GENE[["pep1"]])
  s <- pq("results", "expression_tissue.parquet")[gene_id == GENE[["pep1"]]]
  for (tis in unique(as.character(s$tissue))) {
    v <- x[tissue == tis, cpm]; s1 <- s[tissue == tis]
    expect_equal(s1$n_lines, length(v), label = tis)
    expect_equal(s1$n_detected, sum(v > 0), label = tis)
    expect_equal(s1$median_cpm, median(v), tolerance = 1e-3, label = tis)
    expect_equal(s1$q75_cpm, unname(quantile(v, 0.75)), tolerance = 1e-3, label = tis)
    if (!is.na(s1$b73_cpm)) expect_equal(s1$b73_cpm, x[tissue == tis & genotype == "B73", cpm], tolerance = 1e-3)
  }
})

test_that("the line counts per tissue match the libraries in use", {
  lt <- tsv("results", "expression_lines.tsv")
  expect_equal(nrow(lt), uniqueN(lt[, .(tissue, genotype)]))
  libs <- unlist(lapply(lt$sample_ids, split_ids))
  expect_length(libs, 1836)
  expect_equal(uniqueN(lt$genotype), 296)
  expect_false("B73-10" %in% lt$genotype)
  smp <- tsv("metadata", "samples.tsv")[sample_id %in% libs]
  expect_length(intersect(libs, tsv("metadata", "library_exclusions.tsv")$sample_id), 0)
  expect_false(any(smp$tissue == "L3Base" & smp$flowcell == PILOT_FC))
})

test_that("tau is computed on log2(median CPM + 1) of the seven tissues", {
  sp <- pq("results", "tissue_specificity.parquet")
  et <- pq("results", "expression_tissue.parquet")[tissue %in% LEI_TISSUES]
  set.seed(1); g <- c(GENE[["pep1"]], GENE[["zp27"]], sample(sp[!is.na(tau), gene_id], 20))
  for (id in g) {
    med <- setNames(numeric(7), LEI_TISSUES)
    m <- et[gene_id == id]; med[as.character(m$tissue)] <- m$median_cpm
    x <- log2(med + 1)
    expect_equal(sp[gene_id == id, tau], sum(1 - x / max(x)) / 6, tolerance = 1e-3, label = id)
    expect_equal(as.character(sp[gene_id == id, top_tissue]), names(which.max(med)), label = id)
  }
  expect_true(all(is.na(sp[max_median_cpm < 1, tau])))
  expect_true(all(sp[tau_specific == TRUE, tau] >= 0.8))
})

test_that("Lei's method calls a gene only for its top tissue, in that tissue's top 5% of t and > 2-fold the next", {
  sp <- pq("results", "tissue_specificity.parquet")
  expect_true(all(sp[!is.na(lei_tissue), lei_scored]))
  expect_true(all(sp[!is.na(lei_tissue), as.character(lei_tissue) == as.character(top_tissue)]))
  expect_true(all(sp[!is.na(lei_tissue), fold_vs_next] > 2))
  for (tis in LEI_TISSUES) {
    t <- sp[lei_scored == TRUE, get(paste0("t_", tis))]
    cut <- sort(t, decreasing = TRUE)[ceiling(0.05 * length(t))]
    expect_true(all(sp[lei_tissue == tis, get(paste0("t_", tis))] >= cut - 1e-3), label = tis)
  }
  expect_true(all(is.na(sp[lei_scored == FALSE, t_GRoot])))
})

test_that("known tissue-specific genes are called, and the calls agree with Lei et al.'s list", {
  sp <- pq("results", "tissue_specificity.parquet")
  expect_true(sp[gene_id == GENE[["zp27"]], tau_specific & top_tissue == "Kern"])     # a zein
  expect_equal(as.character(sp[gene_id == GENE[["sh2"]], lei_tissue]), "Kern")       # starch synthesis in kernels
  lei <- pq("annotation", "lei_tissue_specific.parquet")[map == "1:1" & tissue %in% LEI_TISSUES]
  x <- merge(lei[, .(gene_id, tissue)], sp[, .(gene_id, top_tissue = as.character(top_tissue))], by = "gene_id")
  expect_gt(nrow(x), 1000)
  expect_gt(x[, mean(tissue == top_tissue)], 0.95)                                  # 96-100% by tissue (M5)
})

test_that("variation among lines: the share beyond B73 follows its definition", {
  lv <- pq("results", "line_variation.parquet")
  expect_true(all(lv$median_cpm >= 5))
  expect_true(all(is.na(lv[tissue %in% c("L3Mid", "LMid"), v_b73])))
  x <- lv[!is.na(v_b73)]
  expect_true(all(x$n_b73 >= 3))
  expect_true(all(x$share_beyond_b73 >= 0 & x$share_beyond_b73 <= 1))
  expect_equal(x$share_beyond_b73, pmax(0, (x$v_lines - x$v_b73) / x$v_lines), tolerance = 1e-3)
  expect_true(all(is.na(lv[is.na(v_b73), share_beyond_b73])))
})

test_that("genes unstable in B73 are flagged by their definition and kept out of that tissue's partners", {
  u <- tsv("results", "b73_unstable.tsv")
  expect_true(all(as.numeric(u$v_b73) > 4 & as.integer(u$n_b73) >= 3 & as.numeric(u$b73_cpm) >= 5))
  al <- pq("annotation", "gene_aliases.parquet")
  expect_true(any(u$gene_id %in% al[alias_lower == "psba", gene_id]))                # organelle gene copies
  skip_if_not(dir.exists(file.path(ROOT, "results", "partners")), "partner lists not built")
  p <- rbindlist(lapply(unique(substr(u$gene_id, 10, 12)), function(b)
    pq("results", "partners", paste0(b, ".parquet"))[, .(gene_id, scope)]))
  expect_equal(nrow(merge(unique(p), u[, .(gene_id, scope = tissue)], by = c("gene_id", "scope"))), 0)
})
