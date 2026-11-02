# Tissue contrasts (metadata/contrasts.tsv, analysis/de_tissue.R): libraries, pairing and the sign convention,
# log2FC = log2(numerator tissue / denominator tissue).

cons <- tsv("metadata", "contrasts.tsv")
smp <- tsv("metadata", "samples.tsv")
setkey(smp, sample_id)
excl <- tsv("metadata", "library_exclusions.tsv")$sample_id

test_that("metadata/build/validate_contrasts.py finds no problem", {
  py <- Sys.which("python3"); skip_if(!nzchar(py), "python3 not found")
  old <- setwd(ROOT); on.exit(setwd(old))
  out <- system2(py, c("-I", file.path("metadata", "build", "validate_contrasts.py")), stdout = TRUE, stderr = TRUE)
  expect_match(out[1], "^12 contrasts; 0 problems$")
})

test_that("each side holds included libraries of its own tissue, and the sides do not overlap", {
  for (i in seq_len(nrow(cons))) {
    c1 <- cons[i]; num <- split_ids(c1$numerator_samples); den <- split_ids(c1$denominator_samples)
    expect_true(length(num) > 0 && length(den) > 0, info = c1$contrast_id)
    expect_true(all(c(num, den) %in% smp$sample_id), info = c1$contrast_id)
    expect_true(all(smp[c(num, den), include] == "TRUE"), info = c1$contrast_id)
    expect_length(intersect(num, den), 0)
    expect_true(all(smp[num, tissue] == c1$numerator_tissue), info = c1$contrast_id)   # the tissue column (CML505 swap)
    expect_true(all(smp[den, tissue] == c1$denominator_tissue), info = c1$contrast_id)
  }
})

test_that("panel contrasts leave out B73, B73 contrasts hold B73 only; B73-10 and the pilot L3Base are in none", {
  for (i in seq_len(nrow(cons))) {
    c1 <- cons[i]; s <- smp[c(split_ids(c1$numerator_samples), split_ids(c1$denominator_samples))]
    if (c1$design == "panel") expect_false(any(s$genotype %in% c("B73", "B73-10")), info = c1$contrast_id)
    else expect_true(all(s$genotype == "B73"), info = c1$contrast_id)
    expect_false(any(s$tissue == "L3Base" & s$flowcell == PILOT_FC), info = c1$contrast_id)
  }
  expect_setequal(unique(cons$design), c("panel", "B73"))
})

test_that("day and night are paired within line and sampling group; the group contrasts hold one group", {
  dn <- cons[contrast_class == "day_night" & design == "panel"]
  expect_true(all(dn$pair_key == "genotype|sampling_group"))
  for (i in which(nzchar(dn$sampling_group))) {
    s <- smp[c(split_ids(dn$numerator_samples[i]), split_ids(dn$denominator_samples[i]))]
    expect_true(all(s$sampling_group == dn$sampling_group[i]), info = dn$contrast_id[i])
  }
  expect_setequal(dn[nzchar(sampling_group), sampling_group], c("8", "26"))
})

test_that("panel fits use every line with libraries passing QC on both sides", {
  sm <- tsv("results", "de_summary.tsv")
  for (i in which(cons$design == "panel")) {
    c1 <- cons[i]
    key <- function(ids) {
      s <- smp[setdiff(ids, excl)]
      if (c1$pair_key == "genotype|sampling_group") paste(s$genotype, s$sampling_group) else s$genotype
    }
    n <- length(intersect(key(split_ids(c1$numerator_samples)), key(split_ids(c1$denominator_samples))))
    expect_equal(as.integer(sm[contrast_id == c1$contrast_id, n_pairs]), n, label = c1$contrast_id)
  }
})

test_that("log2FC follows the numerator / denominator orientation in every contrast, panel and B73", {
  de <- pq("results", "de.parquet")
  expect_setequal(unique(de$contrast_id), cons$contrast_id)
  for (cid in cons$contrast_id) {
    d <- de[contrast_id == cid & base_mean > 50]
    r <- cor(d$log2fc, log2((d$cpm_numerator_mean + 1) / (d$cpm_denominator_mean + 1)))
    expect_gt(r, 0.9, label = cid)
  }
})

test_that("marker genes change in the expected direction", {
  de <- pq("results", "de.parquet")
  lfc <- function(g, cid) de[gene_id == GENE[[g]] & contrast_id == paste0("Kremling2018_", cid), log2fc]
  for (d in c("panel", "B73")) {
    expect_gt(lfc("pep1", paste0("organ_GShoot_vs_GRoot_", d)), 3)            # C4 photosynthesis: shoot, not root
    expect_gt(lfc("ssu1", paste0("organ_GShoot_vs_GRoot_", d)), 3)
    expect_gt(lfc("pep1", paste0("leaf_maturation_L3Tip_vs_L3Base_", d)), 3)  # mature tip vs dividing base of leaf 3
    expect_lt(lfc("pcna2", paste0("leaf_maturation_L3Tip_vs_L3Base_", d)), -3)
  }
  for (cid in sub("^Kremling2018_", "", cons[contrast_class == "day_night", contrast_id]))
    for (g in c("cca1", "lhy1")) expect_lt(lfc(g, cid), -3, label = paste(g, cid))  # clock genes: mid-night vs mid-day
})
