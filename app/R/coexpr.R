# Co-expression of two genes across lines (Co-expression page), as analysis/coexpr_partners.R: within each tissue,
# the Pearson correlation of log2(CPM + 1) over the lines in view (one value per line; zeros counted). Adult leaf
# (LMAD, LMAN) is centred within sampling group first; "within flowcell" also centres within flowcell (partial r
# given flowcell). Fisher z with SE 1/sqrt(df - 1), df = lines - levels - 1. Tissues where either gene has mean CPM
# below the threshold, or with fewer than 3 df, are left out. The seven tissues of Lei et al. 2025 are pooled with
# DerSimonian-Laird random effects on Fisher z; L3Mid and LMid are shown but not pooled (few lines).

MIN_DF <- 3

coexpr_tissues <- function(xa, xb, within_flowcell = FALSE, min_cpm = 1) {
  x <- merge(xa[, .(tissue, genotype, subpopulation, a = cpm)], xb[, .(tissue, genotype, b = cpm)], by = c("tissue", "genotype"))
  if (!nrow(x)) return(data.table())
  x <- LINE_TISSUE[, .(tissue, genotype, lv = if (within_flowcell) g_flowcell else g_lines)][x, on = c("tissue", "genotype")]
  x[, `:=`(la = log2(a + 1), lb = log2(b + 1))]
  u <- x[, {
    g <- uniqueN(lv)
    ra <- la - ave(la, lv); rb <- lb - ave(lb, lv)
    df <- .N - g - 1
    expressed <- mean(a) >= min_cpm && mean(b) >= min_cpm
    ok <- expressed && df >= MIN_DF && sd(ra) > 0 && sd(rb) > 0
    list(n = .N, levels = g, df = df, mean_a = mean(a), mean_b = mean(b),
         r = if (ok) sum(ra * rb) / sqrt(sum(ra^2) * sum(rb^2)) else NA_real_,
         reason = if (ok) "" else if (!expressed) sprintf("mean CPM < %g", min_cpm) else if (df < MIN_DF) "too few lines" else "no variation")
  }, by = tissue]
  u[, z := atanh(pmin(pmax(r, -0.9999), 0.9999))]
  u[, se := fifelse(is.na(r), NA_real_, 1 / sqrt(df - 1))]
  u[, `:=`(lo = tanh(z - 1.96 * se), hi = tanh(z + 1.96 * se), p = 2 * pnorm(-abs(z / se)))]
  u[order(match(tissue, TISSUES))]
}

# DerSimonian-Laird on Fisher z over the seven main tissues
pool_r <- function(u) {
  u <- u[tissue %in% LEI_TISSUES & is.finite(z) & is.finite(se)]
  if (nrow(u) < 2) return(data.table())
  w <- 1 / u$se^2; k <- nrow(u)
  zf <- sum(w * u$z) / sum(w); Q <- sum(w * (u$z - zf)^2)
  tau2 <- max(0, (Q - (k - 1)) / (sum(w) - sum(w^2) / sum(w)))
  ws <- 1 / (u$se^2 + tau2); zr <- sum(ws * u$z) / sum(ws); s <- sqrt(1 / sum(ws))
  data.table(k = k, n = sum(u$n), r = tanh(zr), lo = tanh(zr - 1.96 * s), hi = tanh(zr + 1.96 * s),
             p = 2 * pnorm(-abs(zr / s)), i2 = if (Q > 0) max(0, (Q - (k - 1)) / Q) else 0)
}

coexpr_forest <- function(u, P, la, lb) {
  d <- u[is.finite(r)]
  if (!nrow(d)) return(empty_plot("No tissue has both genes expressed"))
  rows <- rbind(d[, .(tissue, lab = TISSUE_LABEL[tissue], r, lo, hi, pooled = FALSE,
                      hover = paste0("<b>", TISSUE_LABEL[tissue], "</b><br>r = ", sprintf("%.2f", r), " [", sprintf("%.2f", lo), ", ",
                                     sprintf("%.2f", hi), "]<br>", n, " lines · mean CPM ", signif(mean_a, 3), " (", la, ") / ",
                                     signif(mean_b, 3), " (", lb, ")"))],
                if (nrow(P)) P[, .(tissue = "All", lab = sprintf("Pooled, %d tissues", k), r, lo, hi, pooled = TRUE,
                                   hover = paste0("<b>Pooled</b> over ", k, " tissues<br>r = ", sprintf("%.2f", r), " [", sprintf("%.2f", lo),
                                                  ", ", sprintf("%.2f", hi), "] · I² = ", round(100 * i2), "%<br>p = ", fmtp(p)))])
  rows[, row := .I]
  p <- plot_ly(height = max(240, 26 * nrow(rows) + 70))
  for (i in seq_len(nrow(rows))) {
    q <- rows[i]; col <- if (q$pooled) INK else TISSUE_COL[[q$tissue]]
    p <- add_trace(p, data = q, x = ~r, y = ~row, type = "scatter", mode = "markers", showlegend = FALSE,
                   error_x = list(type = "data", symmetric = FALSE, array = q$hi - q$r, arrayminus = q$r - q$lo, color = col,
                                  thickness = if (q$pooled) 2.2 else 1.4, width = 3),
                   marker = list(color = col, size = if (q$pooled) 13 else 9, symbol = if (q$pooled) "diamond" else "circle",
                                 line = list(color = "white", width = 1)), text = ~hover, hoverinfo = "text")
  }
  base_layout(p,
    xaxis = list(title = "correlation across lines of log2(CPM + 1)", range = c(-1.05, 1.05), zeroline = TRUE, zerolinecolor = MUTED, gridcolor = GRID),
    yaxis = list(title = "", tickvals = rows$row, ticktext = rows$lab, autorange = "reversed", gridcolor = "rgba(0,0,0,0)"),
    margin = list(l = 10, r = 10, t = 10, b = 40))
}

coexpr_scatter <- function(xa, xb, tis, la, lb, within_flowcell) {
  x <- merge(xa[tissue == tis, .(genotype, subpopulation, a = cpm)], xb[tissue == tis, .(genotype, b = cpm)], by = "genotype")
  if (!nrow(x)) return(empty_plot("Pick a tissue"))
  x <- LINE_TISSUE[tissue == tis, .(genotype, lv = if (within_flowcell) g_flowcell else g_lines)][x, on = "genotype"]
  x[, `:=`(la_ = log2(a + 1), lb_ = log2(b + 1))]
  centred <- uniqueN(x$lv) > 1
  if (centred) x[, `:=`(la_ = la_ - ave(la_, lv), lb_ = lb_ - ave(lb_, lv))]
  p <- plot_ly(height = 380)
  for (sp in intersect(SUBPOPS, unique(x$subpopulation))) {
    d <- x[subpopulation == sp]
    p <- add_trace(p, data = d, x = ~la_, y = ~lb_, type = "scatter", mode = "markers", name = sp,
                   marker = list(color = SUBPOP_COL[[sp]], size = if (sp == "Not reported") 6 else 8, line = list(color = "white", width = 0.6)),
                   hovertext = ~paste0("<b>", genotype, "</b> · ", subpopulation, "<br>", la, " ", cpm_txt(a), " CPM · ", lb, " ",
                                       cpm_txt(b), " CPM"), hoverinfo = "text")
  }
  ax <- function(g) paste0(g, if (centred) " (log2 CPM+1, centred)" else " log2(CPM + 1)")
  base_layout(p, xaxis = list(title = ax(la), gridcolor = GRID, zeroline = centred), yaxis = list(title = ax(lb), gridcolor = GRID, zeroline = centred),
              legend = list(orientation = "h", y = -0.22, font = list(size = 11)), margin = list(l = 55, r = 10, t = 10, b = 60))
}
