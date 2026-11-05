# Plotly figures for the MaizExpressD app.

INK <- "#22303A"; MUTED <- "#6B7884"; GRID <- "#E6E9EC"
FONT <- list(family = "IBM Plex Sans, system-ui, sans-serif", size = 12, color = INK)
base_layout <- function(p, ...) {
  layout(p, font = FONT, paper_bgcolor = "rgba(0,0,0,0)", plot_bgcolor = "rgba(0,0,0,0)",
         hoverlabel = list(bgcolor = "white", bordercolor = "#C9D1D8", font = list(color = INK, size = 12)),
         ...) |> config(displaylogo = FALSE, modeBarButtonsToRemove = c("lasso2d", "select2d", "autoScale2d"))
}
fmt <- function(x, d = 2) ifelse(is.na(x), "–", formatC(x, format = "f", digits = d, flag = "+"))
# p-values below 1e-300 have underflowed (the strongest ones to exactly 0): shown as "< 1e-300" in hover text and tables
P_MIN <- 1e-300; P_MIN_TXT <- "< 1e-300"
fmtp <- function(p) ifelse(is.na(p), "–", ifelse(p < P_MIN, P_MIN_TXT,
                                                 ifelse(p < 1e-4, formatC(p, format = "e", digits = 1), formatC(p, format = "f", digits = 4))))
p_tab <- function(p) fifelse(p < P_MIN, 0, signif(p, 3))   # table column: 0 prints as P_MIN_TXT (formatSignif zero.print), sorts first
cpm_txt <- function(x) ifelse(x >= 100, formatC(x, format = "f", digits = 0, big.mark = ","), as.character(signif(x, 3)))
empty_plot <- function(msg) {
  plot_ly(type = "scatter", mode = "markers") |>
    base_layout(xaxis = list(visible = FALSE), yaxis = list(visible = FALSE),
                annotations = list(list(text = msg, showarrow = FALSE, font = list(color = MUTED, size = 13))))
}
CPM_FLOOR <- 0.1   # 0 CPM is drawn at this value on the log scale

# a fixed horizontal offset per line, so one line sits at the same place in every tissue
line_jitter <- function(g) { i <- match(g, sort(LINES$genotype)); ((i * 0.6180339887) %% 1 - 0.5) * 0.62 }

# per-line CPM by tissue: a box over lines (median, quartiles) and one point per line, coloured by subpopulation;
# B73 is a black diamond and highlighted lines are ringed and labelled
expr_plot <- function(x, logy = TRUE, highlight = character()) {
  if (!nrow(x)) return(empty_plot("No lines in view"))
  x <- copy(x)[, tis := factor(tissue, levels = TISSUES)]
  x <- x[!is.na(tis)][, ix := as.integer(droplevels(tis))]
  shown <- levels(droplevels(x$tis))
  x[, `:=`(y = if (logy) pmax(cpm, CPM_FLOOR) else cpm, px = ix + line_jitter(genotype))]
  x[, hover := paste0("<b>", genotype, "</b> · ", subpopulation, "<br>", TISSUE_LABEL[tissue], ": ", cpm_txt(cpm), " CPM",
                      ifelse(n_libraries > 1, paste0(" (mean of ", n_libraries, " libraries)"), ""))]
  p <- plot_ly(height = 430)
  for (t in shown) {
    d <- x[tissue == t]
    p <- add_trace(p, type = "box", x = d$ix, y = d$y, name = t, boxpoints = FALSE, width = 0.72, showlegend = FALSE,
                   fillcolor = paste0(TISSUE_COL[[t]], "22"), line = list(color = TISSUE_COL[[t]], width = 1.4), hoverinfo = "skip")
  }
  for (sp in intersect(SUBPOPS, unique(x$subpopulation))) {
    d <- x[subpopulation == sp & genotype != "B73"]
    p <- add_trace(p, data = d, x = ~px, y = ~y, type = "scatter", mode = "markers", name = sp, legendgroup = sp,
                   marker = list(color = SUBPOP_COL[[sp]], size = 5.5, opacity = 0.8, line = list(color = "white", width = 0.4)),
                   text = ~hover, hoverinfo = "text")
  }
  b <- x[genotype == "B73"]
  if (nrow(b)) p <- add_trace(p, data = b, x = ~px, y = ~y, type = "scatter", mode = "markers", name = "B73",
                              marker = list(color = INK, size = 10, symbol = "diamond", line = list(color = "white", width = 1)),
                              text = ~hover, hoverinfo = "text")
  h <- x[genotype %in% highlight]
  if (nrow(h)) p <- add_trace(p, data = h, x = ~px, y = ~y, type = "scatter", mode = "markers+text", name = "highlighted",
                              text = ~genotype, textposition = "middle right", textfont = list(size = 10, color = INK),
                              hovertext = ~hover, hoverinfo = "text",
                              marker = list(color = "rgba(0,0,0,0)", size = 13, line = list(color = "#C0392B", width = 2)))
  base_layout(p,
    xaxis = list(title = "", tickvals = seq_along(shown), ticktext = paste0("<b>", shown, "</b><br>", TISSUE_LABEL[shown]),
                 range = c(0.4, length(shown) + 0.6), zeroline = FALSE, showgrid = FALSE, tickfont = list(size = 11)),
    yaxis = list(title = "CPM (one point per line)", type = if (logy) "log" else "linear", gridcolor = GRID,
                 dtick = if (logy) 1 else NULL, exponentformat = "none"),
    legend = list(orientation = "h", y = 1.08, x = 0, font = list(size = 11)), margin = list(l = 60, r = 10, t = 30, b = 60))
}

# the 12 tissue contrasts: log2FC with 95% CI; panel (paired by line) filled, B73 open
contrast_plot <- function(d) {
  if (!nrow(d)) return(empty_plot("Not tested in any contrast (too few reads)"))
  d <- CONTRASTS[, .(contrast_id, contrast_class, design, label, n_text, sign_note)][d, on = "contrast_id", nomatch = 0]
  d <- d[order(match(contrast_id, CONTRASTS$contrast_id))][, row := .I]
  d[, ci := fifelse(is.finite(lfc_se), 1.96 * lfc_se, NA_real_)]
  d[, hover := paste0("<b>", label, "</b> · ", CLASS_LABEL[contrast_class], "<br>", sign_note, "<br>log2FC ", fmt(log2fc),
                      ifelse(is.na(ci), "", paste0(" [", fmt(log2fc - ci), ", ", fmt(log2fc + ci), "]")), " · padj ", fmtp(padj),
                      "<br>mean CPM ", cpm_txt(cpm_numerator_mean), " vs ", cpm_txt(cpm_denominator_mean), " · ", n_text,
                      ifelse(single_library %in% TRUE, "<br><i>carried by one library: no p-value</i>", ""))]
  lim <- min(14, max(2, ceiling(max(abs(d$log2fc) + ifelse(is.na(d$ci), 0, d$ci), na.rm = TRUE) + 0.3)))
  p <- plot_ly(height = max(220, 26 * nrow(d) + 60))
  for (cl in intersect(CLASSES, unique(d$contrast_class))) for (ds in c("panel", "B73")) {
    q <- d[contrast_class == cl & design == ds]
    if (!nrow(q)) next
    p <- add_trace(p, data = q, x = ~log2fc, y = ~row, type = "scatter", mode = "markers", showlegend = FALSE,
                   error_x = list(type = "data", array = q$ci, color = CLASS_COL[[cl]], thickness = 1.6, width = 3),
                   marker = list(color = if (ds == "panel") CLASS_COL[[cl]] else "white", size = 10,
                                 symbol = ifelse(!is.na(q$padj) & q$padj < 0.05, if (ds == "panel") "circle" else "diamond",
                                                 if (ds == "panel") "circle-dot" else "diamond-dot"),
                                 line = list(color = CLASS_COL[[cl]], width = 1.8)),
                   text = ~hover, hoverinfo = "text")
  }
  base_layout(p,
    xaxis = list(title = "log2(first / second tissue)", range = c(-lim, lim), zeroline = TRUE, zerolinecolor = MUTED, gridcolor = GRID),
    yaxis = list(title = "", tickvals = d$row, ticktext = paste0(d$label, "  <span style='color:#6B7884'>", d$n_text, "</span>"),
                 autorange = "reversed", gridcolor = "rgba(0,0,0,0)", tickfont = list(size = 11)),
    margin = list(l = 10, r = 10, t = 10, b = 40))
}

contrast_legend <- function() {
  htmltools::span(class = "legend", lapply(CLASSES, function(cl)
    htmltools::span(class = "lg", htmltools::span(class = "sw", style = paste0("background:", CLASS_COL[[cl]])), CLASS_LABEL[[cl]])),
    htmltools::span(class = "lg", htmltools::span(class = "sw open"), "B73 libraries (open)"))
}

volcano_plot <- function(m) {
  if (!nrow(m)) return(empty_plot("No genes"))
  m <- copy(m)[, y := -log10(pmax(padj, 1e-300))]
  m[, col := fifelse(!is.na(padj) & padj < 0.05 & abs(log2fc) >= 1, fifelse(log2fc > 0, "#B2182B", "#2166AC"), "#B8C0C7")]
  plot_ly(m, x = ~log2fc, y = ~y, type = "scattergl", mode = "markers", key = ~gene_id,
          marker = list(color = ~col, size = 4.5, opacity = 0.7),
          hovertext = ~paste0(label, "<br>log2FC ", fmt(log2fc), " · padj ", fmtp(padj), "<br>mean CPM ", cpm_txt(cpm_numerator_mean),
                              " vs ", cpm_txt(cpm_denominator_mean)),
          hoverinfo = "text", source = "volcano", height = 420) |>
    event_register("plotly_click") |>
    base_layout(xaxis = list(title = "log2(first / second tissue)", zeroline = TRUE, zerolinecolor = MUTED, gridcolor = GRID),
                yaxis = list(title = "−log10 adjusted p", gridcolor = GRID), margin = list(l = 55, r = 10, t = 10, b = 45))
}

heat_plot <- function(Z, hov, zscore) {
  if (!nrow(Z)) return(empty_plot("No matched genes"))
  z <- as.matrix(Z[, -1, with = FALSE])
  plot_ly(x = colnames(Z)[-1], y = Z$gene, z = z, type = "heatmap", xgap = 2, ygap = 2, text = hov, hoverinfo = "text",
          zmin = if (zscore) -2.5 else 0, zmax = if (zscore) 2.5 else NULL, zmid = if (zscore) 0 else NULL,
          colorscale = if (zscore) list(c(0, "#2166AC"), c(0.5, "#F2F2F0"), c(1, "#B2182B"))
                       else list(c(0, "#F4F6F2"), c(0.35, "#A9CF8F"), c(0.7, "#2E7D45"), c(1, "#123D22")),
          colorbar = list(title = if (zscore) "z (row)" else "log2(median CPM + 1)", thickness = 12),
          height = max(280, 22 * nrow(Z) + 170)) |>
    base_layout(xaxis = list(tickangle = -30), yaxis = list(autorange = "reversed", tickfont = list(family = "IBM Plex Mono, monospace", size = 11)),
                margin = list(l = 10, r = 10, t = 10, b = 90))
}

# Lei's eQTLs along the v5 genome (Lei's liftover), coloured by tissue; the gene itself as a vertical line
eqtl_plot <- function(e, g) {
  e <- e[!is.na(eqtl_pos_v5) & eqtl_chrom_v5 %in% CHROMS$chrom]
  if (!nrow(e)) return(empty_plot("No eQTL with a v5 position"))
  e <- CHROMS[, .(eqtl_chrom_v5 = chrom, offset)][e, on = "eqtl_chrom_v5"]
  e[, `:=`(x = offset + eqtl_pos_v5, y = -log10(pvalue))]
  gx <- if (nrow(g) && g$chrom %in% CHROMS$chrom) CHROMS[chrom == g$chrom, offset] + (g$start + g$end) / 2 else NA
  e[, kind := fifelse(!is.na(gx) & eqtl_chrom_v5 == g$chrom & abs(eqtl_pos_v5 - (g$start + g$end) / 2) < 1e6, "cis (< 1 Mb)", "trans")]
  e[, hover := paste0("<b>", TISSUE_LABEL[tissue], "</b> eQTL<br>", eqtl_chrom_v5, ":", format(eqtl_pos_v5, big.mark = ","), " (v5) · ", kind,
                      "<br>p = ", fmtp(pvalue), "<br>Lei gene ", gene_v3, " (", map, ")")]
  p <- plot_ly(height = 300)
  for (t in intersect(TISSUES, unique(e$tissue))) {
    d <- e[tissue == t]
    p <- add_trace(p, data = d, x = ~x, y = ~y, type = "scatter", mode = "markers", name = TISSUE_LABEL[[t]],
                   marker = list(color = TISSUE_COL[[t]], size = 10, line = list(color = "white", width = 1)), text = ~hover, hoverinfo = "text")
  }
  shapes <- lapply(CHROMS$offset[-1], function(o) list(type = "line", x0 = o, x1 = o, yref = "paper", y0 = 0, y1 = 1,
                                                       line = list(color = GRID, width = 1)))
  if (!is.na(gx)) shapes <- c(shapes, list(list(type = "line", x0 = gx, x1 = gx, yref = "paper", y0 = 0, y1 = 1,
                                                line = list(color = "#C0392B", width = 1.5, dash = "dot"))))
  base_layout(p, shapes = shapes,
    xaxis = list(title = "", tickvals = CHROMS$offset + CHROMS$len / 2, ticktext = sub("chr", "", CHROMS$chrom), range = c(0, sum(CHROMS$len)),
                 showgrid = FALSE, zeroline = FALSE),
    yaxis = list(title = "−log10 p", gridcolor = GRID, rangemode = "tozero"),
    legend = list(orientation = "h", y = 1.15, x = 0), margin = list(l = 55, r = 10, t = 30, b = 30))
}

subpop_legend <- function() {
  htmltools::span(class = "legend", lapply(SUBPOPS, function(s)
    htmltools::span(class = "lg", htmltools::span(class = "sw", style = paste0("background:", SUBPOP_COL[[s]])), s)),
    htmltools::span(class = "lg", htmltools::span(class = "sw dia"), "B73"))
}
