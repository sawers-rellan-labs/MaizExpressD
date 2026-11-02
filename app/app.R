# MaizExpressD: maize gene expression across tissues and stages in the 282 diversity panel (Kremling et al. 2018),
# with Lei et al. 2025's eQTL, TWAS and QTL results. Adapted from MaizExpressC.
#   shiny::runApp("app")        (from the project root; data in app/data, built by app/prepare_data.R)
# Set MAIZEXPRESSD_DATA to use a data folder elsewhere.

suppressPackageStartupMessages({
  library(shiny); library(bslib); library(plotly); library(DT); library(data.table); library(arrow)
  library(markdown)   # used by includeMarkdown(); loaded explicitly so deployment tools see the dependency
})
# helpers in R/ (coexpr.R, data.R, maizemine.R, plots.R) are sourced by Shiny before this file

theme <- bs_theme(
  version = 5, bg = "#FFFFFF", fg = "#22303A", primary = "#2F5D3A", secondary = "#6B7884",
  base_font = font_google("IBM Plex Sans", local = FALSE), heading_font = font_google("Spectral", wght = c(500, 600), local = FALSE),
  code_font = font_google("IBM Plex Mono", local = FALSE), "border-radius" = "8px", "font-size-base" = "0.93rem"
)

START_GENE <- "Zm00001eb219040"                                       # ga20ox4: kernel-high, Lei eQTL and TWAS hit
CX_DEFAULT <- c(a = "Zm00001eb383680", b = "Zm00001eb287770")         # pep1, pdk1 (C4 photosynthesis)
GENE_LIST <- "pep1 pdk1 ppck1 cah1 ssu1 cca1 lhy1 gi1 pcna2 rtcs1 rth3 o2 zp27 z1A-1 sh2 bt2 ae1 ga20ox4"
LINE_CHOICES <- sort(LINES$genotype)

filters <- sidebar(
  width = 270, open = "desktop",
  checkboxGroupInput("f_tissue", "Tissues", choices = setNames(TISSUES, paste0(TISSUES, " · ", TISSUE_LABEL[TISSUES])), selected = TISSUES),
  checkboxGroupInput("f_subpop", "Subpopulations (lines)", choices = SUBPOPS, selected = SUBPOPS),
  tags$small(class = "text-muted", "Lei et al. 2025's groups; the MR landraces are not in Lei's set."),
  selectizeInput("f_highlight", "Highlight lines", choices = LINE_CHOICES, multiple = TRUE, options = list(placeholder = "e.g. Mo17, B97")),
  input_switch("f_log", "Log scale for CPM", TRUE),
  actionButton("reset", "Reset", class = "btn-sm btn-outline-secondary"),
  textOutput("n_lines", inline = TRUE)
)

search_input <- function(id, placeholder, label = NULL)
  selectizeInput(id, label, choices = NULL, width = "100%",
                 options = list(placeholder = placeholder, maxOptions = 40, searchField = c("label", "search"), loadThrottle = 250))

lookup_page <- nav_panel(
  "Gene lookup", value = "lookup",
  card(class = "searchcard",
       search_input("gene", "Search a gene: symbol (pep1, ga20ox4, o2), v5/v4/v3 ID, Arabidopsis name (CCA1, GA20OX2) or description"),
       uiOutput("genecard")),
  card(full_screen = TRUE, card_header("Expression across tissues and lines", span(class = "hint", "one point per line; box = median and quartiles over lines")),
       plotlyOutput("expr", height = "auto"),
       p(class = "hint", "CPM: reads on the gene per million reads on genes (3′ tags, no length correction); a line with several libraries in a tissue shows their mean. ",
         "On the log scale, 0 CPM is drawn at 0.1. Each line keeps the same horizontal place in every tissue. Lines and tissues follow the sidebar.")),
  layout_columns(col_widths = c(7, 5),
    card(full_screen = TRUE, card_header("Tissue contrasts", span(class = "hint", "log2FC with 95% CI; dotted marker: padj ≥ 0.05")),
         div(class = "legendrow", contrast_legend()), plotlyOutput("contrasts", height = "auto"),
         p(class = "hint", "Panel contrasts pair each line with itself (edgeR quasi-likelihood, ~ line + tissue); B73 contrasts use B73's replicate libraries (DESeq2). ",
           "Adult leaf by night vs day is fitted within each sampling date as well as across both.")),
    card(full_screen = TRUE, card_header("Variation among lines", span(class = "hint", "per tissue")),
         DTOutput("variation"),
         p(class = "hint", "SD: standard deviation of log2(CPM + 1) among lines (one library each). Beyond B73: the share of that variance above the variance among B73's own replicate libraries ",
           "(plant, batch and depth noise), for genes with median CPM ≥ 5. Shown when B73 has 3 or more libraries in the tissue."))),
  card(full_screen = TRUE, card_header("Baseline expression in B73", span(class = "hint", textOutput("atlas_hint", inline = TRUE))),
       input_switch("atlas_log", "log scale", TRUE), plotlyOutput("atlas", height = "auto"),
       p(class = "hint", "MaizeGDB expression atlas via MaizeMine: mean FPKM per tissue and stage in B73 (full-length RNA-seq, other plants and conditions). ",
         "Hover a bar for the tissue; 0 FPKM has no bar on the log scale.")),
  card(card_header("Contrast results", downloadButton("dl_gene", "Download CSV", class = "btn-sm")), DTOutput("restable"))
)

traits_page <- nav_panel(
  "Traits & eQTL", value = "traits",
  card(card_header(textOutput("tr_title", inline = TRUE), span(class = "hint", "Lei et al. 2025, G3 15:jkaf140 (B73 v3 genes, matched to v5)")),
       uiOutput("tr_summary"),
       plotlyOutput("eqtl_plot", height = "auto"),
       p(class = "hint", "eQTLs of the gene's expression in six tissues (kernel had too few), at Lei's v5 positions. The dotted red line is the gene. ",
         "Lei mapped v3 gene models; where a v3 gene matches several v5 genes (or the reverse), its results are shown on each, with the match type.")),
  layout_columns(col_widths = c(6, 6),
    card(card_header("eQTLs"), DTOutput("tr_eqtl")),
    card(card_header("Trait associations (TWAS) and colocalized QTLs"), tags$h6("TWAS"), DTOutput("tr_twas"),
         tags$h6("QTLs colocalized with the gene's eQTLs"), DTOutput("tr_coloc"))),
  card(card_header("All TWAS genes", span(class = "hint", "click a row to open the gene")), DTOutput("tr_all"),
       p(class = "hint", "Category 1: associated in one tissue only; category 2: in several. Trait codes are listed in Methods."))
)

lines_page <- nav_panel(
  "Lines", value = "lines",
  card(card_header("Lines", span(class = "hint", "click a line to see its libraries and its expression of the current gene")),
       DTOutput("line_table")),
  layout_columns(col_widths = c(7, 5),
    card(full_screen = TRUE, card_header(textOutput("line_title", inline = TRUE)), plotlyOutput("line_expr", height = "auto"),
         p(class = "hint", "The selected line is ringed in red among all lines in view.")),
    card(card_header("Libraries and rank among lines"), DTOutput("line_rank"), tags$h6("Libraries"), DTOutput("line_libs")))
)

coexpr_page <- nav_panel(
  "Co-expression", value = "coexpr",
  card(class = "searchcard",
    layout_columns(col_widths = c(5, 5, 2),
      search_input("cg_a", "Search gene A", "Gene A"), search_input("cg_b", "Search gene B", "Gene B"),
      div(input_switch("cg_fc", "Within flowcell", FALSE), numericInput("cg_min", "Min mean CPM", 1, min = 0, step = 0.5, width = "100%"))),
    p(class = "hint", "Within each tissue, the correlation of log2(CPM + 1) across lines (one value per line); adult leaf is centred within its two sampling dates. ",
      "“Within flowcell” also removes each sequencing flowcell’s mean, which takes out batch effects but also some genetic differences between groups of lines.")),
  card(card_header(textOutput("pt_title", inline = TRUE), span(class = "hint", "precomputed · click a row to make it gene B"),
                   downloadButton("dl_pt", "Download CSV", class = "btn-sm")),
    radioButtons("pt_scope", NULL, c("All tissues" = "All", setNames(LEI_TISSUES, LEI_TISSUES)), inline = TRUE),
    navset_underline(
      nav_panel("Positively co-expressed", div(class = "pt-tables", DTOutput("pt_pos"))),
      nav_panel("Negatively co-expressed", div(class = "pt-tables", DTOutput("pt_neg")))),
    p(class = "hint", "Top 25 partners by r in one tissue, or by pooled r over the seven main tissues (k: tissues where both are expressed); all lines, min mean CPM 1, ",
      "following the “Within flowcell” switch. Genes whose B73 libraries disagree on and off are left out of that tissue. ",
      "Near-identical paralogs share reads and can co-vary for that reason: “nearby” marks partners within 100 kb.")),
  layout_columns(col_widths = c(6, 6),
    card(full_screen = TRUE, card_header("Co-expression by tissue", span(class = "hint", "r with 95% CI · diamond: pooled (random effects)")),
         plotlyOutput("cx_forest", height = "auto"), uiOutput("cx_note")),
    card(full_screen = TRUE, card_header("Lines in one tissue"), selectInput("cx_tissue", NULL, choices = NULL),
         plotlyOutput("cx_scatter", height = "auto")))
)

rank_page <- nav_panel(
  "Rank genes", value = "rank",
  card(card_header("By tissue contrast", span(class = "hint", "click a point or a row to open the gene")),
       layout_columns(col_widths = c(6, 3, 3),
         selectInput("rk_contrast", "Contrast", choices = setNames(CONTRASTS$contrast_id, paste(CONTRASTS$label, "·", CONTRASTS$n_text)),
                     selected = "Kremling2018_leaf_maturation_L3Tip_vs_L3Base_panel"),
         numericInput("rk_lfc", "|log2FC| ≥", 1, min = 0, step = 0.5),
         selectInput("rk_sort", "Rank by", c("adjusted p" = "padj", "log2FC, up" = "up", "log2FC, down" = "down"))),
       plotlyOutput("volcano", height = "auto"), DTOutput("rk_table")),
  card(card_header("Tissue-specific genes"),
       layout_columns(col_widths = c(3, 9),
         selectInput("sp_tissue", "Tissue", choices = setNames(LEI_TISSUES, paste(LEI_TISSUES, "·", TISSUE_LABEL[LEI_TISSUES]))),
         radioButtons("sp_method", "Call", inline = TRUE,
                      c("tau ≥ 0.8 (restricted to the tissue)" = "tau", "Lei's method, this data" = "lei_here",
                        "Lei et al. 2025, as published" = "lei_paper"))),
       DTOutput("sp_table"),
       p(class = "hint", "tau (Yanai et al. 2005) on the per-tissue medians of the line values, log2(CPM + 1), seven tissues. ",
         "Lei's method: genes expressed in all seven tissues whose t statistic for the tissue is in its top 5% and whose median is over 2-fold every other tissue's; ",
         "so it finds tissue-enriched genes that are expressed everywhere, and tau finds genes restricted to one tissue. ",
         "Lei's published lists (B73 v3) are shown on v5 genes with a one-to-one match."))
)

genelist_page <- nav_panel(
  "Gene list", value = "genelist",
  card(card_header("Gene list"),
       p(class = "lede", "Paste symbols or IDs (v5, v4, v3 or Arabidopsis), separated by spaces, commas or new lines. An identifier that maps to several genes brings in all of them. ",
         "The starting list holds C4 and photosynthesis genes, clock genes, a cell-division gene, root genes and kernel genes."),
       textAreaInput("gl_text", NULL, rows = 3, width = "100%", value = GENE_LIST),
       layout_columns(col_widths = c(3, 3, 2, 4), actionButton("gl_go", "Show heatmap", class = "btn-primary"),
                      selectInput("gl_order", NULL, c("as entered" = "input", "by top tissue" = "top", "by similarity" = "cluster")),
                      input_switch("gl_z", "z per gene", FALSE), uiOutput("gl_matched"))),
  card(full_screen = TRUE, card_header("Genes × tissues", span(class = "hint", "median CPM over the lines in view"),
                                       downloadButton("dl_list", "Download CSV", class = "btn-sm")),
       plotlyOutput("heat", height = "auto"))
)

data_page <- nav_panel(
  "Data", value = "data",
  card(card_header("Kremling et al. 2018: 282 diversity panel, seven tissues"), uiOutput("kpis"),
       p(class = "lede", "3′ mRNA-seq of the maize 282 diversity panel (Kremling et al. 2018, Nature 555:520; SRA PRJNA383416), requantified on B73 RefGen_v5. ",
         "Tissue definitions follow Lei et al. 2025 and Giri et al. 2021."), DTOutput("tissue_table")),
  card(card_header("Tissue contrasts"), DTOutput("contrast_table"),
       p(class = "hint", "Up and down: genes at padj < 0.05 and at least 2-fold. With 180 or more line pairs most genes reach padj < 0.05, so rank by effect size."))
)

methods_page <- nav_panel("Methods", value = "methods", card(class = "methods", includeMarkdown("methods.md")))

ui <- page_navbar(
  id = "nav", title = span(class = "brand", "MaizExpressD"), theme = theme, sidebar = filters, fillable = FALSE,
  header = tags$head(tags$link(rel = "stylesheet", href = "styles.css")),
  lookup_page, traits_page, lines_page, coexpr_page, rank_page, genelist_page, data_page, methods_page,
  nav_spacer(), nav_item(span(class = "built", paste("data", BUILT)))
)

server <- function(input, output, session) {
  dt_opts <- function(...) list(dom = "tip", scrollX = TRUE, ...)
  note_dt <- function(msg) datatable(data.table(note = msg), rownames = FALSE, selection = "none", options = list(dom = "t"))

  # ---- gene search choices (server side)
  alias_txt <- ALIASES[, .(search = paste(unique(alias), collapse = " ")), by = gene_id]
  choices <- GENES[alias_txt, on = "gene_id"][, .(value = gene_id,
    label = paste0(ifelse(nzchar(symbol), paste0(symbol, " · "), ""), gene_id,
                   ifelse(nzchar(description), paste0(" · ", substr(description, 1, 60)), "")),
    search = paste(search, description, full_name))]
  updateSelectizeInput(session, "gene", choices = choices, server = TRUE, selected = START_GENE)
  gene <- reactive({ g <- input$gene; if (is.null(g) || !nzchar(g)) START_GENE else g })
  go_gene <- function(g) { updateSelectizeInput(session, "gene", choices = choices, server = TRUE, selected = g); nav_select("nav", "lookup") }

  observeEvent(input$reset, {
    updateCheckboxGroupInput(session, "f_tissue", selected = TISSUES)
    updateCheckboxGroupInput(session, "f_subpop", selected = SUBPOPS)
    updateSelectizeInput(session, "f_highlight", selected = character())
    update_switch("f_log", value = TRUE)
  })
  lines_in_view <- reactive(LINES[subpopulation %in% input$f_subpop, genotype])
  output$n_lines <- renderText(sprintf("%d of %d lines", length(lines_in_view()), nrow(LINES)))
  in_view <- function(x) x[tissue %in% input$f_tissue & genotype %in% lines_in_view()]

  ex_all <- reactive(get_expr(gene()))
  ex <- reactive(in_view(ex_all()))

  # ---- gene card
  output$genecard <- renderUI({
    gid <- gene(); g <- GENES[.(gid)]
    t <- EXPR_T[.(gid), nomatch = 0]
    top <- if (nrow(t)) t[which.max(median_cpm)] else NULL
    s <- SPEC[.(gid), nomatch = 0]
    de <- DE[.(gid), nomatch = 0]
    eq <- LEI_EQTL[.(gid), nomatch = 0]; tw <- LEI_TWAS[.(gid), nomatch = 0]; lts <- LEI_TS[.(gid), nomatch = 0]
    co <- LEI_COLOC[.(gid), nomatch = 0]; un <- UNSTABLE[.(gid), nomatch = 0]
    gos <- if (nzchar(g$go_terms)) strsplit(g$go_terms, ",")[[1]] else character()
    pred <- grepl("^EnTAP", g$description_source)
    stat <- function(v, l) div(class = "stat", div(class = "v", v), div(class = "l", l))
    lei_txt <- if (!nrow(eq) && !nrow(tw) && !nrow(lts) && !nrow(co)) span(class = "muted", "no eQTL, TWAS or tissue-specific call") else tagList(
      if (nrow(lts)) paste0("tissue-specific in ", paste(unique(lts$tissue_lei), collapse = ", "), "; "),
      if (nrow(eq)) sprintf("eQTLs in %s; ", paste(intersect(TISSUES, unique(eq$tissue)), collapse = ", ")),
      if (nrow(tw)) sprintf("TWAS: %s; ", paste(unique(paste0(tw$traits, " (", tw$tissue, ")")), collapse = "; ")),
      if (nrow(co)) sprintf("eQTLs colocalize with QTLs for %d traits; ", uniqueN(co$trait)),
      tags$a(href = "#", onclick = "Shiny.setInputValue('go_traits', Math.random()); return false;", "details ›"),
      if (any(c(eq$map, tw$map, lts$map) != "1:1")) span(class = "muted", " (v3 gene matches several v5 genes)"))
    div(class = "genecard",
      div(class = "gc-main",
        div(class = "sym", if (nzchar(g$symbol)) g$symbol else g$gene_id,
            if (nrow(un)) span(class = "tag warn", title = "Its B73 replicate libraries disagree on and off: often a nuclear copy of an organelle gene collecting organelle reads. Treat its CPM with care.",
                               paste("unstable in B73:", paste(intersect(TISSUES, un$tissue), collapse = ", ")))),
        div(class = "ids", span(class = "mono", g$gene_id),
            a(href = paste0("https://www.maizegdb.org/gene_center/gene/", g$gene_id), target = "_blank", "MaizeGDB gene page ↗")),
        tags$dl(
          tags$dt("Description"), tags$dd(if (nzchar(g$description)) g$description else span(class = "muted", "none"),
                                          if (pred) span(class = "tag", "predicted (EnTAP)")),
          if (nzchar(g$full_name)) tagList(tags$dt("MaizeGDB name"), tags$dd(g$full_name)),
          tags$dt("Location (v5)"), tags$dd(class = "mono", sprintf("%s:%s–%s (%s)", g$chrom, format(g$start, big.mark = ","), format(g$end, big.mark = ","), g$strand)),
          tags$dt("v4 / v3"), tags$dd(class = "mono", paste(c(if (nzchar(g$v4_ids)) g$v4_ids else "–", if (nzchar(g$v3_ids)) g$v3_ids else "–"), collapse = " / ")),
          tags$dt("Arabidopsis"), tags$dd(if (nzchar(g$ath_best_hit)) tagList(span(class = "mono", g$ath_best_hit), if (nzchar(g$ath_symbol)) paste0(" ", g$ath_symbol),
                                          if (lengths(strsplit(g$ath_orthologs, ",")) > 1) span(class = "muted", paste0("  (orthologs: ", gsub(",", ", ", g$ath_orthologs), ")"))) else span(class = "muted", "no Ensembl ortholog")),
          tags$dt("Lei et al. 2025"), tags$dd(lei_txt),
          tags$dt("Pathways"), tags$dd(pathway_list(PATHWAYS[.(gid), nomatch = 0])),
          tags$dt("GWAS hits"), tags$dd(gwas_list(GWAS_HITS[.(gid), nomatch = 0])),
          tags$dt("Variants"), tags$dd(variant_summary(VARIANTS[.(gid), nomatch = 0], gid)),
          tags$dt("Pan-genome"), tags$dd(pangene_presence(PANGENES[.(gid), nomatch = 0]), " ", pangene_strip(PANGENES[.(gid), nomatch = 0])),
          tags$dt("Literature"), tags$dd(literature_list(LITERATURE[.(gid), nomatch = 0], gid)),
          if (length(gos)) tagList(tags$dt("GO"), tags$dd(class = "small", paste(head(paste0(go_name[gos], " (", gos, ")"), 6), collapse = "; "),
                                                         if (length(gos) > 6) span(class = "muted", sprintf(" … %d more", length(gos) - 6)))))),
      div(class = "stats",
        stat(if (!is.null(top)) cpm_txt(top$median_cpm) else "–",
             if (!is.null(top)) sprintf("highest median CPM: %s (%d lines)", TISSUE_LABEL[top$tissue], top$n_lines) else "not detected"),
        stat(if (nrow(s) && !is.na(s$tau)) sprintf("%.2f", s$tau) else "–",
             if (!nrow(s) || is.na(s$tau)) "tau: no tissue median reaches 1 CPM"
             else paste0(if (s$tau_specific) sprintf("tau ≥ 0.8: specific to %s", s$top_tissue) else "tau (0 = even, 1 = one tissue)",
                         if (!is.na(s$lei_tissue)) sprintf(" · Lei's method: %s", s$lei_tissue) else "")),
        stat(sprintf("%d / %d", de[!is.na(padj) & padj < 0.05 & abs(log2fc) >= 1, .N], nrow(de)),
             "contrasts with a 2-fold change at padj < 0.05"),
        stat(if (nrow(eq)) uniqueN(eq$tissue) else "0", if (nrow(eq)) sprintf("tissues with an eQTL (Lei; %d eQTLs)", nrow(unique(eq[, .(tissue, eqtl_pos_v3)]))) else "tissues with an eQTL (Lei)")))
  })
  observeEvent(input$go_traits, nav_select("nav", "traits"))

  # ---- expression, contrasts, variation
  output$expr <- renderPlotly(expr_plot(ex(), isTRUE(input$f_log), input$f_highlight))
  de_gene <- reactive(DE[.(gene()), nomatch = 0])
  output$contrasts <- renderPlotly(contrast_plot(de_gene()))
  output$variation <- renderDT({
    t <- EXPR_T[.(gene()), nomatch = 0][tissue %in% input$f_tissue]
    if (!nrow(t)) return(note_dt("not detected in the tissues in view"))
    v <- LVAR[.(gene()), nomatch = 0]
    t <- v[, .(tissue, v_lines, v_b73, share_beyond_b73)][t, on = "tissue"][order(match(tissue, TISSUES))]
    datatable(t[, .(tissue, lines = n_lines, detected = paste0(round(100 * n_detected / n_lines), "%"),
                    `median CPM [IQR]` = sprintf("%s [%s–%s]", cpm_txt(median_cpm), cpm_txt(q25_cpm), cpm_txt(q75_cpm)),
                    B73 = ifelse(is.na(b73_cpm), "–", cpm_txt(b73_cpm)), SD = ifelse(is.na(v_lines), "–", sprintf("%.2f", sqrt(v_lines))),
                    `beyond B73` = ifelse(is.na(share_beyond_b73), "–", paste0(round(100 * share_beyond_b73), "%")))],
              rownames = FALSE, selection = "none", class = "compact stripe", options = list(dom = "t", ordering = FALSE, scrollX = TRUE))
  })
  atlas_x <- reactive(get_atlas(gene()))
  output$atlas <- renderPlotly(atlas_plot(atlas_x(), isTRUE(input$atlas_log)))
  output$atlas_hint <- renderText(atlas_summary(atlas_x()))
  res_tab <- reactive({
    d <- CONTRASTS[de_gene(), on = "contrast_id", nomatch = 0][order(match(contrast_id, CONTRASTS$contrast_id))]
    d[, .(contrast = label, class = class_label, lines_or_libraries = n_text, model, log2FC = round(log2fc, 3), SE = round(lfc_se, 3),
          `log2FC (shrunk)` = round(log2fc_shrunk, 3), padj = signif(padj, 3), `CPM first` = round(cpm_numerator_mean, 2),
          `CPM second` = round(cpm_denominator_mean, 2), `one library` = single_library %in% TRUE, contrast_id)]
  })
  output$restable <- renderDT(datatable(res_tab()[, !"contrast_id"], rownames = FALSE, class = "compact stripe", options = dt_opts(pageLength = 12)) |>
                                formatSignif("padj", 3))
  output$dl_gene <- downloadHandler(filename = function() paste0("MaizExpressD_", gene(), "_contrasts.csv"),
                                    content = function(f) fwrite(res_tab(), f))

  # ---- traits & eQTL (Lei et al. 2025)
  output$tr_title <- renderText(paste("Lei et al. 2025 results for", short_name(gene())))
  output$tr_summary <- renderUI({
    gid <- gene(); eq <- LEI_EQTL[.(gid), nomatch = 0]; tw <- LEI_TWAS[.(gid), nomatch = 0]
    co <- LEI_COLOC[.(gid), nomatch = 0]; lts <- LEI_TS[.(gid), nomatch = 0]; cu <- LEI_CUBIC[.(gid), nomatch = 0]
    ca <- LEI_CAND[.(gid), nomatch = 0]
    k <- function(v, l) div(class = "stat", div(class = "v", v), div(class = "l", l))
    v3 <- unique(c(eq$gene_v3, tw$gene_v3, co$gene_v3, lts$gene_v3))
    tagList(div(class = "kpis",
      k(nrow(unique(eq[, .(tissue, eqtl_pos_v3)])), sprintf("eQTLs in %d tissues", uniqueN(eq$tissue))),
      k(nrow(tw), if (nrow(tw)) paste("TWAS hits:", paste(unique(tw$tissue), collapse = ", ")) else "TWAS hits"),
      k(uniqueN(co$trait), "traits with a colocalized QTL"),
      k(if (nrow(lts)) paste(unique(lts$tissue_lei), collapse = ", ") else "–", "Lei's tissue-specific call")),
      if (length(v3)) p(class = "hint", "Lei's v3 gene", if (length(v3) > 1) "s", ": ", span(class = "mono", paste(v3, collapse = ", ")),
                        " · match ", paste(unique(c(eq$map, tw$map, co$map, lts$map)), collapse = ", ")),
      if (nrow(cu)) p(class = "hint", "Checked in the CUBIC population (Lei Table S9): ", paste(sprintf("%s p = %s", cu$trait, signif(cu$pvalue, 2)), collapse = "; ")),
      if (nrow(ca)) p(class = "hint", "Lei's candidate (Table S10): ", paste(unique(na.omit(c(ca$symbol_lei, ca$role_reference))), collapse = "; ")))
  })
  output$eqtl_plot <- renderPlotly(eqtl_plot(LEI_EQTL[.(gene()), nomatch = 0], GENES[.(gene())]))
  output$tr_eqtl <- renderDT({
    eq <- LEI_EQTL[.(gene()), nomatch = 0]
    if (!nrow(eq)) return(note_dt("no eQTL in Lei et al. 2025"))
    g <- GENES[.(gene())]
    eq[, cis := fifelse(eqtl_chrom_v5 == g$chrom & abs(eqtl_pos_v5 - (g$start + g$end) / 2) < 1e6, "cis", "trans", na = "–")]
    datatable(eq[order(match(tissue, TISSUES), pvalue), .(tissue, `eQTL (v5)` = ifelse(is.na(eqtl_pos_v5), "–", sprintf("%s:%s", eqtl_chrom_v5, format(eqtl_pos_v5, big.mark = ","))),
                     `cis/trans` = cis, p = signif(pvalue, 3), `Lei gene` = gene_v3, match = map)],
              rownames = FALSE, selection = "none", class = "compact stripe", options = dt_opts(pageLength = 8)) |> formatSignif("p", 3)
  })
  output$tr_twas <- renderDT({
    tw <- LEI_TWAS[.(gene()), nomatch = 0]
    if (!nrow(tw)) return(note_dt("no TWAS hit"))
    datatable(tw[, .(tissue, traits, category, annotation, `Lei gene` = gene_v3, match = map)], rownames = FALSE, selection = "none",
              class = "compact stripe", options = list(dom = "t", scrollX = TRUE))
  })
  output$tr_coloc <- renderDT({
    co <- LEI_COLOC[.(gene()), nomatch = 0]
    if (!nrow(co)) return(note_dt("no colocalized QTL"))
    datatable(co[order(window != "100 bp", trait), .(trait, tissue, window, `QTL (v5)` = qtl_pos_v5, `eQTL (v5)` = ifelse(is.na(eqtl_pos_v5), "–",
                     paste0(eqtl_chrom_v5, ":", format(eqtl_pos_v5, big.mark = ","))))],
              rownames = FALSE, selection = "none", class = "compact stripe", options = dt_opts(pageLength = 8))
  })
  tw_all <- LEI_TWAS[map != "none"][GENES[, .(gene_id, symbol)], on = "gene_id", nomatch = 0][order(tissue, traits)]
  output$tr_all <- renderDT(datatable(tw_all[, .(gene = fifelse(nzchar(symbol), symbol, gene_id), gene_id, tissue, traits, category, annotation)],
                                      rownames = FALSE, selection = "single", class = "compact stripe hover", options = dt_opts(pageLength = 10)))
  observeEvent(input$tr_all_rows_selected, go_gene(tw_all$gene_id[input$tr_all_rows_selected]))

  # ---- lines
  line_tab <- reactive(LINES[subpopulation %in% input$f_subpop][order(genotype)])
  output$line_table <- renderDT(datatable(line_tab()[, .(line = genotype, `Lei name` = genotype_lei, subpopulation, origin, pedigree,
                                                         `adult-leaf group` = adult_leaf_group, tissues, libraries = n_libraries)],
                                          rownames = FALSE, selection = list(mode = "single", selected = which(line_tab()$genotype == "B73")),
                                          class = "compact stripe hover", options = dt_opts(pageLength = 10)))
  sel_line <- reactive({ i <- input$line_table_rows_selected; if (length(i)) line_tab()$genotype[i] else "B73" })
  output$line_title <- renderText(sprintf("%s in %s", short_name(gene()), sel_line()))
  output$line_expr <- renderPlotly(expr_plot(ex(), isTRUE(input$f_log), unique(c(sel_line(), input$f_highlight))))
  output$line_rank <- renderDT({
    x <- ex()[, .(genotype, cpm, rank = frank(-cpm, ties.method = "min"), n = .N), by = tissue][genotype == sel_line()]
    if (!nrow(x)) return(note_dt("this line has no library in the tissues in view"))
    datatable(x[order(match(tissue, TISSUES)), .(tissue, CPM = cpm_txt(cpm), `rank among lines` = sprintf("%d of %d", rank, n))],
              rownames = FALSE, selection = "none", class = "compact stripe", options = list(dom = "t", ordering = FALSE))
  })
  output$line_libs <- renderDT(datatable(SAMPLES[genotype == sel_line()][order(match(tissue, TISSUES)),
    .(library = sample_id, tissue, `sampling group` = sampling_group, flowcell, `mapped reads (M)` = round(reads_mapped / 1e6, 2), `% mapped` = round(percent_mapped, 1))],
    rownames = FALSE, selection = "none", class = "compact stripe", options = dt_opts(pageLength = 10)))

  # ---- co-expression
  updateSelectizeInput(session, "cg_a", choices = choices, server = TRUE, selected = CX_DEFAULT[["a"]])
  updateSelectizeInput(session, "cg_b", choices = choices, server = TRUE, selected = CX_DEFAULT[["b"]])
  ga <- reactive({ g <- input$cg_a; if (is.null(g) || !nzchar(g)) CX_DEFAULT[["a"]] else g })
  gb <- reactive({ g <- input$cg_b; if (is.null(g) || !nzchar(g)) CX_DEFAULT[["b"]] else g })
  pt_all <- reactive(get_partners(ga()))
  pt <- reactive({
    p <- pt_all(); if (!nrow(p)) return(p)
    p <- p[scope == input$pt_scope & view == if (isTRUE(input$cg_fc)) "within_flowcell" else "lines"][order(direction, rank)]
    g <- GENES[.(p$partner_id)]; a <- GENES[.(ga())]
    d <- ifelse(g$chrom == a$chrom, pmax(0, pmax(g$start, a$start) - pmin(g$end, a$end)) / 1e3, Inf)
    p[, `:=`(symbol = fifelse(is.na(g$symbol), "", g$symbol), description = fifelse(is.na(g$description), "", g$description),
             nearby = fifelse(!is.na(d) & d < 100, sprintf("%.0f kb", d), ""))]
  })
  output$pt_title <- renderText(paste("Strongest partners of", short_name(ga())))
  pt_dt <- function(d) {
    if (!nrow(d)) return(note_dt(sprintf("No list for %s here: it needs mean CPM ≥ 1 in this tissue (in at least 4 of the 7 for All tissues).", short_name(ga()))))
    esc <- htmltools::htmlEscape
    short <- fifelse(nchar(d$description) > 60, paste0(substr(d$description, 1, 59), "…"), d$description)
    datatable(d[, .(`#` = rank, gene = esc(fifelse(nzchar(symbol), symbol, "–")), `gene ID` = partner_id,
                    description = sprintf('<span title="%s">%s</span>', esc(description), esc(short)),
                    `r [95% CI]` = sprintf("%.2f [%.2f, %.2f]", r, lo, hi), k = units, nearby)],
              rownames = FALSE, escape = FALSE, selection = "single", class = "compact stripe hover nowrap",
              options = list(dom = "tp", pageLength = 10, ordering = FALSE, scrollX = TRUE))
  }
  output$pt_pos <- renderDT(pt_dt(pt()[direction == "+"]))
  output$pt_neg <- renderDT(pt_dt(pt()[direction == "-"]))
  pick_partner <- function(dir, i) {
    d <- pt()[direction == dir]
    if (length(i) && i <= nrow(d)) updateSelectizeInput(session, "cg_b", choices = choices, server = TRUE, selected = d$partner_id[i])
  }
  observeEvent(input$pt_pos_rows_selected, pick_partner("+", input$pt_pos_rows_selected))
  observeEvent(input$pt_neg_rows_selected, pick_partner("-", input$pt_neg_rows_selected))
  output$dl_pt <- downloadHandler(
    filename = function() sprintf("MaizExpressD_partners_%s_%s_%s.csv", short_name(ga()), input$pt_scope, if (isTRUE(input$cg_fc)) "within_flowcell" else "lines"),
    content = function(f) fwrite(pt()[, .(gene_id, scope, view, direction = fifelse(direction == "+", "positive", "negative"), rank, partner_id,
                                          partner_symbol = symbol, tissues = units, r, ci_low = lo, ci_high = hi, i2, nearby, partner_description = description)], f))
  xa <- reactive(in_view(get_expr(ga()))); xb <- reactive(in_view(get_expr(gb())))
  cx_u <- reactive(coexpr_tissues(xa(), xb(), isTRUE(input$cg_fc), max(0, input$cg_min %||% 1)))
  output$cx_forest <- renderPlotly(coexpr_forest(cx_u(), pool_r(cx_u()), short_name(ga()), short_name(gb())))
  output$cx_note <- renderUI({
    u <- cx_u(); if (!nrow(u)) return(NULL)
    out <- u[!is.finite(r)]
    p(class = "hint", sprintf("%d tissues shown; the pooled estimate uses the seven main ones in view.", u[is.finite(r), .N]),
      if (nrow(out)) paste0(" Left out: ", paste(sprintf("%s (%s)", out$tissue, out$reason), collapse = "; "), "."))
  })
  observe({
    u <- cx_u()[is.finite(r)]; ch <- u$tissue; cur <- isolate(input$cx_tissue)
    updateSelectInput(session, "cx_tissue", choices = setNames(ch, TISSUE_LABEL[ch]),
                      selected = if (!is.null(cur) && cur %in% ch) cur else if (nrow(u)) u[which.max(abs(r)), tissue] else character())
  })
  output$cx_scatter <- renderPlotly({ req(input$cx_tissue)
    coexpr_scatter(xa(), xb(), input$cx_tissue, short_name(ga()), short_name(gb()), isTRUE(input$cg_fc)) })

  # ---- rank genes
  rank_data <- reactive({
    m <- DE[contrast_id == input$rk_contrast]
    m <- GENES[, .(gene_id, symbol, description)][m, on = "gene_id"]
    m[, label := ifelse(nzchar(symbol), paste0(symbol, " · ", gene_id), gene_id)]
  })
  output$volcano <- renderPlotly(volcano_plot(rank_data()))
  rank_tab <- reactive({
    m <- rank_data()[!is.na(padj) & padj < 0.05 & abs(log2fc) >= input$rk_lfc]
    m <- switch(input$rk_sort, padj = m[order(padj, -abs(log2fc))], up = m[order(-log2fc)], down = m[order(log2fc)])
    m[, .(gene_id, symbol, description = substr(description, 1, 60), log2FC = round(log2fc, 2), SE = round(lfc_se, 3), padj = signif(padj, 3),
          `CPM first` = round(cpm_numerator_mean, 1), `CPM second` = round(cpm_denominator_mean, 1))]
  })
  output$rk_table <- renderDT(datatable(rank_tab(), rownames = FALSE, selection = "single", class = "compact stripe hover", options = dt_opts(pageLength = 15)) |>
                                formatSignif("padj", 3))
  observeEvent(input$rk_table_rows_selected, go_gene(rank_tab()$gene_id[input$rk_table_rows_selected]))
  observeEvent(event_data("plotly_click", source = "volcano"), { k <- event_data("plotly_click", source = "volcano")$key; if (length(k)) go_gene(k[[1]]) })
  sp_tab <- reactive({
    tis <- input$sp_tissue
    ids <- switch(input$sp_method,
                  tau = SPEC[tau_specific & top_tissue == tis, gene_id],
                  lei_here = SPEC[lei_tissue == tis, gene_id],
                  lei_paper = LEI_TS[map == "1:1" & tissue == tis, unique(gene_id)])
    s <- SPEC[.(ids), nomatch = 0]
    m <- EXPR_T[tissue == tis, .(gene_id, median_cpm)][s, on = "gene_id"]
    m <- GENES[, .(gene_id, symbol, description)][m, on = "gene_id"]
    tcol <- paste0("t_", tis)
    m[, t := get(tcol)]
    m[order(-tau, -median_cpm), .(gene_id, symbol, description = substr(description, 1, 60), `median CPM` = round(median_cpm, 1),
                                  `top tissue` = top_tissue, `fold over next` = ifelse(is.infinite(fold_vs_next), "∞", sprintf("%.1f", fold_vs_next)),
                                  tau = round(tau, 2), `Lei t` = round(t, 1))]
  })
  output$sp_table <- renderDT(datatable(sp_tab(), rownames = FALSE, selection = "single", class = "compact stripe hover",
                                        caption = sprintf("%d genes", nrow(sp_tab())), options = dt_opts(pageLength = 10)))
  observeEvent(input$sp_table_rows_selected, go_gene(sp_tab()$gene_id[input$sp_table_rows_selected]))

  # ---- gene list
  gl_req <- reactiveVal(NULL)
  observeEvent(input$gl_go, gl_req(list(tokens = strsplit(input$gl_text, "[\\s,;]+", perl = TRUE)[[1]])), ignoreNULL = FALSE)
  observeEvent(input$pw_pick, {
    ids <- PATHWAYS[pathway_id == input$pw_pick, unique(gene_id)]
    if (!length(ids)) return()
    updateTextAreaInput(session, "gl_text", value = paste(ids, collapse = "\n"))
    gl_req(list(tokens = ids, exact = TRUE))
    nav_select("nav", "genelist")
  })
  gl <- reactive({
    q <- gl_req(); req(q)
    r <- if (isTRUE(q$exact)) list(genes = q$tokens, missing = character()) else resolve_tokens(q$tokens)
    ids <- head(r$genes, 300)
    if (!length(ids)) return(list(r = r, Z = data.table(), hov = NULL))
    tis <- intersect(TISSUES, input$f_tissue)
    # medians over the lines in view (the subpopulation filter applies), from the per-line tables
    e <- rbindlist(lapply(split(ids, bucket(ids)), function(g) read_bucket("expr", g)))
    lt <- LINE_TISSUE[tissue %in% tis & genotype %in% lines_in_view(), .(tissue, genotype)]
    full <- lt[, .(gene_id = ids), by = .(tissue, genotype)]
    full[e, cpm := i.cpm, on = c("gene_id", "tissue", "genotype")][is.na(cpm), cpm := 0]
    d <- full[, .(m = median(cpm)), by = .(gene_id, tissue)]
    Z <- dcast(d, gene_id ~ tissue, value.var = "m")[, c("gene_id", tis), with = FALSE]
    Z <- Z[match(intersect(ids, Z$gene_id), gene_id)]
    M <- as.matrix(Z[, -1]); raw <- M; L <- log2(M + 1)
    if (isTRUE(input$gl_z)) { L <- t(scale(t(L))); L[!is.finite(L)] <- 0 }
    o <- switch(input$gl_order, input = seq_len(nrow(Z)), top = order(match(tis[max.col(raw, ties.method = "first")], tis), -apply(raw, 1, max)),
                cluster = if (nrow(Z) > 2) hclust(dist(L))$order else seq_len(nrow(Z)))
    Z <- Z[o]; L <- L[o, , drop = FALSE]; raw <- raw[o, , drop = FALSE]
    hov <- matrix(sprintf("%s<br>%s<br>median %s CPM", rep(gene_label(Z$gene_id), ncol(raw)), rep(TISSUE_LABEL[tis], each = nrow(raw)), cpm_txt(raw)), nrow = nrow(raw))
    out <- data.table(gene = gene_label(Z$gene_id), L); setnames(out, c("gene", tis))
    list(r = r, Z = out, hov = hov, raw = cbind(Z[, .(gene_id)], raw))
  })
  output$gl_matched <- renderUI({ r <- gl()$r
    p(class = "hint", sprintf("%d genes matched", length(r$genes)), if (length(r$missing)) span(class = "warn", paste0(" · not found: ", paste(r$missing, collapse = ", ")))) })
  output$heat <- renderPlotly({ x <- gl(); if (!nrow(x$Z)) return(empty_plot("No matched genes")); heat_plot(x$Z, x$hov, isTRUE(input$gl_z)) })
  output$dl_list <- downloadHandler(filename = function() "MaizExpressD_gene_list_median_cpm.csv", content = function(f) fwrite(gl()$raw, f))

  # ---- data
  output$kpis <- renderUI({
    k <- function(v, l) div(class = "stat", div(class = "v", v), div(class = "l", l))
    div(class = "kpis", k(format(nrow(SAMPLES), big.mark = ","), "libraries in use"), k(nrow(LINES), "lines"), k(length(TISSUES), "tissues"),
        k(nrow(CONTRASTS), "tissue contrasts"), k(format(GENES[quantified == TRUE, .N], big.mark = ","), "genes quantified"))
  })
  output$tissue_table <- renderDT(datatable(TISSUE_TAB[, .(tissue, name = TISSUE_LABEL[tissue], `Lei code` = lei_code, stage, lines = n_lines,
                                                           libraries = n_libraries, detail = tissue_detail)],
                                            rownames = FALSE, selection = "none", class = "compact stripe", options = list(dom = "t", scrollX = TRUE, ordering = FALSE)))
  output$contrast_table <- renderDT(datatable(CONTRASTS[, .(contrast = label, class = class_label, `lines or libraries` = n_text, model,
                                                            `genes tested` = genes_tested, up = n_up, down = n_down, notes)],
                                              rownames = FALSE, selection = "none", class = "compact stripe", options = list(dom = "t", scrollX = TRUE, ordering = FALSE)))
}

shinyApp(ui, server)
