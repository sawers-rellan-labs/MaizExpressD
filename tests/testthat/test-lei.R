# Lei et al. 2025 results on v5 genes (annotation/build_lei.R): every v3-v5 match kept and labelled, tissues in
# our codes.

xw <- local({
  g <- tsv("annotation", "genes.tsv")[v3_ids != "", .(gene_id, v3_ids)]
  g[, .(gene_v3 = unlist(strsplit(v3_ids, ",", fixed = TRUE))), by = gene_id]
})

test_that("each lifted gene carries every v5 match of its v3 ID and the right label", {
  for (t in c("tissue_specific", "eqtls", "twas", "coloc", "cubic", "candidates")) {
    d <- unique(pq("annotation", paste0("lei_", t, ".parquet"))[, .(gene_v3, gene_id, map)])
    expect_true(all(d$map %in% c("1:1", "1:n", "n:1", "n:m", "none")), label = t)
    expect_true(all(is.na(d[map == "none", gene_id]) | d[map == "none", gene_id] == ""), label = t)
    m <- d[map != "none"]
    expect_true(all(paste(m$gene_v3, m$gene_id) %in% paste(xw$gene_v3, xw$gene_id)), label = t)
    want <- xw[gene_v3 %in% m$gene_v3]                                      # every match of those v3 IDs
    expect_setequal(paste(m$gene_v3, m$gene_id), paste(want$gene_v3, want$gene_id))
    n_v5 <- xw[, .N, by = gene_v3]; n_v3 <- xw[, .N, by = gene_id]
    one <- m[map == "1:1"]
    expect_true(all(n_v5[match(one$gene_v3, gene_v3), N] == 1 & n_v3[match(one$gene_id, gene_id), N] == 1), label = t)
  }
})

test_that("most of Lei's genes have a v5 match", {
  for (t in c("tissue_specific", "eqtls", "twas")) {
    d <- pq("annotation", paste0("lei_", t, ".parquet"))
    expect_gt(d[, .(mapped = any(map != "none")), by = gene_v3][, mean(mapped)], 0.9, label = t)
  }
})

test_that("Lei's tissue codes become ours", {
  ok <- c(LEI_TISSUES, "L3Tip+LMAD+LMAN")
  for (t in c("tissue_specific", "eqtls", "twas", "coloc", "candidates")) {
    d <- pq("annotation", paste0("lei_", t, ".parquet"))
    expect_true(all(d$tissue %in% ok | is.na(d$tissue)), label = t)
  }
  expect_setequal(unique(pq("annotation", "lei_tissue_specific.parquet")$tissue), ok)       # L3: Lei's leaf-group list
})
