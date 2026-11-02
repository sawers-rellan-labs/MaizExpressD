# Co-expression across lines (app/R/coexpr.R and analysis/coexpr_partners.R): per-tissue r, the partial correlations
# behind the two views, DerSimonian-Laird pooling, and the precomputed partner lists.

line_values <- function(e, ga, gb, tis, view) {
  x <- merge(e$get_expr(ga)[tissue == tis, .(genotype, a = cpm)], e$get_expr(gb)[tissue == tis, .(genotype, b = cpm)], by = "genotype")
  x[e$LINE_TISSUE[tissue == tis], lv := if (view == "lines") i.g_lines else i.g_flowcell, on = "genotype"]
  x
}

test_that("per-tissue r is the correlation of log2(CPM + 1) across lines, zeros included", {
  e <- app_env()
  u <- e$coexpr_tissues(e$get_expr(GENE[["pep1"]]), e$get_expr(GENE[["pdk1"]]), within_flowcell = FALSE, min_cpm = 1)
  tt <- u[tissue %in% c("GRoot", "GShoot", "L3Base", "L3Tip", "Kern") & is.finite(r), tissue]   # one level: no centring
  expect_gte(length(tt), 3)
  expect_true(is.na(u[tissue == "Kern", r]))                                            # pdk1 below 1 CPM in kernels
  for (tis in tt) {
    x <- line_values(e, GENE[["pep1"]], GENE[["pdk1"]], tis, "lines")
    expect_equal(u[tissue == tis, n], nrow(x))
    expect_equal(u[tissue == tis, r], cor(log2(x$a + 1), log2(x$b + 1)), label = tis)
    expect_equal(u[tissue == tis, se], 1 / sqrt(nrow(x) - 3), label = tis)              # df = n - 1 level - 1
  }
})

test_that("adult leaf is centred within sampling group, and within flowcell gives the partial r given flowcell", {
  e <- app_env()
  xa <- e$get_expr(GENE[["pep1"]]); xb <- e$get_expr(GENE[["pdk1"]])
  for (view in c("lines", "within_flowcell")) {
    u <- e$coexpr_tissues(xa, xb, within_flowcell = view == "within_flowcell", min_cpm = 1)
    for (tis in c("LMAD", "L3Base")) {
      x <- line_values(e, GENE[["pep1"]], GENE[["pdk1"]], tis, view)
      g <- uniqueN(x$lv)
      pr <- if (g > 1) cor(resid(lm(log2(a + 1) ~ lv, x)), resid(lm(log2(b + 1) ~ lv, x))) else cor(log2(x$a + 1), log2(x$b + 1))
      expect_equal(u[tissue == tis, r], pr, label = paste(view, tis))
      expect_equal(u[tissue == tis, se], 1 / sqrt(nrow(x) - g - 2), label = paste(view, tis))
    }
  }
  expect_gt(e$LINE_TISSUE[tissue == "LMAD", uniqueN(g_lines)], 1)                       # groups 8 and 26
  expect_equal(e$LINE_TISSUE[tissue == "L3Base", uniqueN(g_lines)], 1)
})

test_that("pooled correlation matches metafor DerSimonian-Laird on Fisher z over the seven tissues", {
  skip_if_not_installed("metafor")
  e <- app_env()
  u <- e$coexpr_tissues(e$get_expr(GENE[["pep1"]]), e$get_expr(GENE[["pdk1"]]), FALSE, 1)
  P <- e$pool_r(u); ok <- u[tissue %in% LEI_TISSUES & is.finite(r)]
  m <- metafor::rma(yi = z, sei = se, data = ok, method = "DL")
  expect_equal(P$k, nrow(ok))
  expect_equal(P$r, tanh(as.numeric(m$b)))
  expect_equal(P$lo, tanh(m$ci.lb), tolerance = 1e-4)                                   # app uses 1.96, metafor qnorm(0.975)
  expect_equal(P$i2, m$I2 / 100, tolerance = 1e-6)
})

test_that("a gene with itself correlates perfectly", {
  e <- app_env()
  x <- e$get_expr(GENE[["pep1"]])
  u <- e$coexpr_tissues(x, x, FALSE, 1)
  expect_true(all(abs(u[is.finite(r), r] - 1) < 1e-9))
})

test_that("precomputed partner lists are well formed and match the two-gene view", {
  e <- app_env()
  p <- e$get_partners(GENE[["pep1"]])
  expect_true(all(p$scope %in% c(LEI_TISSUES, "All")))
  expect_setequal(unique(p$view), c("lines", "within_flowcell"))
  expect_false(any(p$partner_id == GENE[["pep1"]]))
  for (s in split(p, by = c("scope", "view", "direction"))) {
    expect_equal(sort(s$rank), seq_len(nrow(s)))
    s <- s[order(rank)]
    expect_true(all(diff(if (s$direction[1] == "+") -s$r else s$r) >= -1e-9))
    expect_lte(nrow(s), 25)
  }
  chk <- p[scope == "All" & rank <= 2]
  for (i in seq_len(nrow(chk))) {
    u <- e$coexpr_tissues(e$get_expr(GENE[["pep1"]]), e$get_expr(chk$partner_id[i]), chk$view[i] == "within_flowcell", 1)
    P <- e$pool_r(u)
    expect_equal(P$k, chk$units[i])
    expect_equal(P$r, chk$r[i], tolerance = 1e-3)
    expect_equal(P$lo, chk$lo[i], tolerance = 1e-3)
  }
  one <- p[scope == "GShoot" & view == "lines" & direction == "+" & rank == 1]
  u <- e$coexpr_tissues(e$get_expr(GENE[["pep1"]]), e$get_expr(one$partner_id), FALSE, 1)
  expect_equal(u[tissue == "GShoot", r], one$r, tolerance = 1e-3)
})

test_that("known partners: pdk1 with pep1 (C4 cycle), cca1 with lhy1 (clock), o2 with zeins", {
  e <- app_env()
  pos <- function(g, sc) e$get_partners(GENE[[g]])[scope == sc & view == "lines" & direction == "+"]
  expect_true(GENE[["pdk1"]] %in% pos("pep1", "All")$partner_id)
  expect_true(GENE[["lhy1"]] %in% pos("cca1", "All")$partner_id)
  zein <- e$GENES[grepl("zein", description, ignore.case = TRUE) | grepl("^z1|^zp|^zein", symbol), gene_id]
  expect_gte(sum(pos("o2", "Kern")$partner_id %in% zein), 3)
})
