# Gene-card rows from MaizeMine (MaizeGDB's data warehouse), built by annotation/fetch_maizemine.py:
# metabolic pathways (MaizeGDB E2P2), GWAS hits within 10 kb (GWAS Atlas, Wallace 2014 NAM GWAS),
# protein-changing variants (Ensembl variant set, with SIFT) and the gene's pan-gene group across
# the B73 and NAM founder genomes. Links open the MaizeMine or Ensembl record.

pathway_list <- function(pw) {
  if (!nrow(pw)) return(span(class = "muted", "none in the MaizeGDB pathway set"))
  tagList(lapply(seq_len(nrow(pw)), function(i) tagList(
    if (i > 1) "; ",
    a(href = maizemine_link(pw$pathway_id[i], "Pathway"), target = "_blank", pw$pathway_name[i]),
    " ", tags$a(href = "#", class = "pw-go small", title = "Show this pathway's genes in Gene list",
                onclick = sprintf("Shiny.setInputValue('pw_pick', '%s', {priority: 'event'}); return false;", pw$pathway_id[i]),
                "genes \u203a"))))
}

gwas_list <- function(h, show = 5) {
  if (!nrow(h)) return(span(class = "muted", "no GWAS hit within 10 kb"))
  t <- h[, .(n = .N, d = min(distance_bp), src = paste(sort(unique(sub(" GWAS$", "", source))), collapse = ", ")), by = trait][order(d, -n)]
  one <- function(i) tagList(if (i > 1) "; ", span(title = sprintf("%d hit%s, nearest %s (%s)", t$n[i], if (t$n[i] > 1) "s" else "",
                                                                  if (t$d[i] == 0) "inside the gene" else paste(format(t$d[i], big.mark = ","), "bp away"), t$src[i]),
                                                   gsub("_", " ", t$trait[i]), if (t$d[i] == 0) span(class = "tag", "in gene")))
  tagList(sprintf("%d trait%s: ", nrow(t), if (nrow(t) > 1) "s" else ""),
          lapply(seq_len(min(show, nrow(t))), one),
          if (nrow(t) > show) tags$details(class = "more", tags$summary(sprintf("all %d", nrow(t))),
                                           div(class = "small", lapply(seq_len(nrow(t)), one))))
}

variant_summary <- function(v, gene_id) {
  link <- a(href = paste0("https://plants.ensembl.org/Zea_mays/Gene/Variation_Gene/Table?db=core;g=", gene_id),
            target = "_blank", "Ensembl variant table ↗")
  if (!nrow(v)) return(tagList(span(class = "muted", "no protein-changing variant listed"), " · ", link))
  v <- v[1]
  lof <- c(`stop gained` = v$stop_gained, frameshift = v$frameshift, `splice site` = v$splice_site,
           `start lost` = v$start_lost, `stop lost` = v$stop_lost)
  lof <- lof[lof > 0]
  tagList(
    if (length(lof)) span(class = "tag lof", paste(sprintf("%d %s", lof, names(lof)), collapse = ", ")),
    if (length(lof)) " ",
    sprintf("%d missense (%d predicted deleterious by SIFT)", v$missense, v$missense_deleterious),
    " · ", link)
}

pangene_presence <- function(pg, show = 6) {
  if (!nrow(pg)) return(span(class = "muted", "not in a MaizeGDB pan-gene group"))
  pg <- pg[1]
  miss <- if (nzchar(pg$missing)) strsplit(pg$missing, ",", fixed = TRUE)[[1]] else character()
  tagList(
    a(href = maizemine_link(pg$pangene_id, "PangeneGroup"), target = "_blank", pg$pangene_id), ": ",
    if (!length(miss)) sprintf("in all %d genomes (B73 and %d NAM founders)", pg$n_genomes, pg$n_genomes - 1L)
    else tagList(sprintf("in %d of %d genomes; no copy in ", pg$genomes_present, pg$n_genomes),
                 paste(head(miss, show), collapse = ", "),
                 if (length(miss) > show) span(class = "muted", sprintf(" and %d more", length(miss) - show))),
    if (pg$b73_genes_in_group > 1) span(class = "muted", sprintf(" · %d B73 genes in this group", pg$b73_genes_in_group)))
}

# presence of the gene's pan-gene group in each of the 26 genomes (B73 first), as a strip of cells
pangene_strip <- function(pg) {
  g <- MAIZEMINE$genomes
  if (!nrow(pg) || is.null(g)) return(NULL)
  miss <- strsplit(pg$missing[1], ",", fixed = TRUE)[[1]]
  span(class = "pg-strip", lapply(g, function(x) span(class = paste("pg-cell", if (x %in% miss) "off" else "on"),
                                                       title = paste0(x, if (x %in% miss) ": no copy" else ": has a copy"))))
}

literature_list <- function(lit, gene_id, show = 3) {
  if (!nrow(lit)) return(span(class = "muted", "no gene-specific paper in NCBI Gene"))
  lit <- lit[order(-suppressWarnings(as.integer(year)), pmid)]
  one <- function(i) tagList(if (i > 1) "; ",
    a(href = paste0("https://pubmed.ncbi.nlm.nih.gov/", lit$pmid[i], "/"), target = "_blank",
      title = sprintf("%s (%s)%s", lit$title[i], lit$journal[i], if (lit$n_genes[i] > 1) sprintf(", %d genes", lit$n_genes[i]) else ""),
      sprintf("%s %s", sub(" .*$", "", lit$first_author[i]), lit$year[i])))
  tagList(sprintf("%d gene-specific paper%s: ", nrow(lit), if (nrow(lit) > 1) "s" else ""),
          lapply(seq_len(min(show, nrow(lit))), one),
          if (nrow(lit) > show) tagList(" ", a(href = paste0("https://pubmed.ncbi.nlm.nih.gov/?term=", paste(lit$pmid, collapse = "+OR+"), "[uid]"),
                                               target = "_blank", sprintf("all %d in PubMed ↗", nrow(lit)))))
}

ATLAS_COL <- c(Root = "#7A5C3A", Leaf = "#009E73", Shoot = "#0072B2", Reproductive = "#E69F00", Seed = "#CC79A7")

# bars of mean FPKM per B73 atlas sample, ordered and coloured by organ group
atlas_plot <- function(x, log_scale = TRUE) {
  if (!nrow(x) || !nrow(ATLAS_SAMPLES)) return(plotly_empty() |> layout(title = list(text = "not in the B73 expression atlas", font = list(size = 13))))
  d <- ATLAS_SAMPLES[x, on = "sample", nomatch = 0][order(order)]
  d[is.na(fpkm), fpkm := 0]
  d[, lab := sprintf("%s%s<br>%s: %s FPKM", tissue, fifelse(nzchar(stage), paste0(", ", stage), ""), organ_group, signif(fpkm, 3))]
  d[, x := factor(sample, levels = sample)]
  p <- plot_ly(height = 240)
  for (og in intersect(names(ATLAS_COL), unique(d$organ_group))) {
    s <- d[organ_group == og]
    p <- add_bars(p, data = s, x = ~x, y = ~fpkm, name = og, marker = list(color = ATLAS_COL[[og]]),
                  text = ~lab, hoverinfo = "text", textposition = "none")
  }
  p |> layout(barmode = "overlay", bargap = 0.15, margin = list(l = 50, r = 10, t = 10, b = 10),
              xaxis = list(title = "", showticklabels = FALSE, categoryorder = "array", categoryarray = levels(d$x)),
              yaxis = list(title = "FPKM", type = if (log_scale) "log" else "linear", rangemode = "tozero"),
              legend = list(orientation = "h", y = 1.15, x = 0)) |>
    config(displaylogo = FALSE, modeBarButtonsToRemove = c("lasso2d", "select2d"))
}

atlas_summary <- function(x, n = 3) {
  if (!nrow(x) || !nrow(ATLAS_SAMPLES)) return("")
  d <- ATLAS_SAMPLES[x, on = "sample", nomatch = 0][order(-fpkm)][seq_len(min(n, .N))]
  paste("highest:", paste(sprintf("%s%s (%s)", d$tissue, fifelse(nzchar(d$stage), paste0(", ", d$stage), ""), signif(d$fpkm, 3)), collapse = "; "))
}
