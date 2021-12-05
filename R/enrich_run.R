#' Enrichment analysis of go
#'
#' run enrich_go and plot results
#'
#' @inheritParams enrich_go
#' @param showCategory Category numbers to show
#' @param dir where to save results files
#' @param prefix a prefix of file names in this step
#'
#' @importFrom utils write.table
#' @importFrom graphics strwidth
#'
#' @return a list result files of GO
#' @export
#'
#' @examples
#' \dontrun{
#' enrichGO_run(DEG_df, x = "log2FoldChange", y = "pvalue", dir = tempdir())
#' }
enrichGO_run <- function(deg_data, x, y, cut_FC = 1, cut_P = 0.05, showCategory = 10, dir = ".", prefix = "3-EnrichGO",
                       OrgDb = 'org.Hs.eg.db', keyType = "SYMBOL", ont = "ALL", simplify = TRUE,
                       pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2,minGSSize = 10,
                       maxGSSize = 500, readable = FALSE, pool = FALSE,
                       label = c("Down", "Stable", "Up"),
                       label_ns = "Stable",
                       mc.cores = 1L) {

  if (!fs::dir_exists(dir)) {
    fs::dir_create(dir)
  }

  go_resl <- enrich_go(deg_data, x, y, cut_FC = cut_FC, cut_P = cut_P,
            OrgDb = OrgDb, keyType = keyType, ont = ont, simplify = simplify,
            pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod, qvalueCutoff = qvalueCutoff,minGSSize = minGSSize,
            maxGSSize = maxGSSize, readable = readable, pool = pool,
            label = label,
            label_ns = label_ns,
            mc.cores = 1L)

  go_res_file <- glue("{dir}/{prefix}_go_result.Rdata")
  save(go_resl,file = go_res_file)
  message(glue::glue("Result of enrichGO stored in {go_res_file}"))

  tmp <- lapply(seq_along(go_resl), function(x)
    write.table(go_resl[[x]]@result,file = glue::glue('{dir}/{prefix}_gene_{names(go_resl)[[x]]}_GO_enrichment.csv'))
  )

  if (ont == "ALL") {
    plots <- lapply(go_resl, function(x)
      enhance_barplot(x@result,showCategory=showCategory,split = 'ONTOLOGY')
    )
    text_w_l <- lapply(go_resl, function(x)
      max(strwidth(x@result$Description,units = "inch"))
    )
    lapply(seq_along(plots), function(x) {

      text_w <- text_w_l[[x]]
      ggsave(plot = plots[[x]],
             filename = glue::glue("{dir}/{prefix}_barplot-gene_{names(plots)[[x]]}_GO_enrichment.pdf"),
             height = 880*2.5/300, width = text_w+1, dpi = 300)

    }
    )
  } else {
    plots <- lapply(go_resl, function(x)
      enhance_barplot(x@result,showCategory=showCategory)
    )
    text_w_l <- lapply(go_resl, function(x)
      max(strwidth(x@result$Description,units = "inch"))
    )
    lapply(seq_along(plots), function(x) {

      text_w <- text_w_l[[x]]
      ggsave(plot = plots[[x]], filename = glue::glue("{dir}/{prefix}_barplot-gene_{names(plots)[[x]]}_GO_enrichment.pdf"),
             height = 880*2.5/300/3, width = text_w+1, dpi = 300)

    }

    )
  }

  message(glue::glue("Result of enrichGO ploted in {dir}"))

}

#' Enrichment analysis of kegg
#'
#' run enrich_kegg and plot results
#'
#' @inheritParams enrich_kegg
#' @inheritParams kegg_barplot
#' @param top top rows for up and down
#' @param dir where to save results files
#' @param prefix a prefix of file names in this step
#'
#' @return a list result files of KEGG
#' @export
#'
#' @examples
#' \dontrun{
#' enrichKEGG_run(deg_data = DEG_df, x = "log2FoldChange", y = "pvalue", cut_FC = 1,
#' cut_P = 0.05, top = -10)
#' }
enrichKEGG_run <- function(deg_data, x, y, cut_FC = 1, cut_P = 0.05, top = -10,
                           dir= ".", prefix = "4-EnrichKEGG",
                           organism = "hsa",
                           keyType = "kegg",
                           OrgDb = 'org.Hs.eg.db',
                           pvalueCutoff = 0.05,
                           pAdjustMethod = "BH",
                           minGSSize = 10,
                           maxGSSize = 500,
                           qvalueCutoff = 0.2,
                           use_internal_data = FALSE,
                           label = c("Down", "Stable", "Up"),
                           label_ns = "Stable",
                           mc.cores = 1L,down_label = "Down"
                           ) {

  if (!fs::dir_exists(dir)) {
    fs::dir_create(dir)
  }
  # enrich_kegg
  kegg_resl <- enrich_kegg(deg_data = deg_data, x = x, y = y, cut_FC = cut_FC,
                           cut_P = cut_P,
                           organism = organism,
                           keyType = keyType,
                           OrgDb = OrgDb,
                           pvalueCutoff = pvalueCutoff,
                           pAdjustMethod = pAdjustMethod,
                           minGSSize = minGSSize,
                           maxGSSize = maxGSSize,
                           qvalueCutoff = qvalueCutoff,
                           use_internal_data = use_internal_data,
                           label = label,
                           label_ns = label_ns,
                           mc.cores = mc.cores)

  kegg_res_file <- glue("{dir}/{prefix}_kegg_result.Rdata")
  save(kegg_resl,file = kegg_res_file)
  message(glue::glue("Result of enrich_kegg stored in {dir}/{prefix}_{kegg_res_file}"))

  text_w <- max(strwidth(kegg_resl$Down@result$Description,units = "inch"),
      strwidth(kegg_resl$Up@result$Description,units = "inch"))

  p <- kegg_barplot(kegg_resl, top = top, down_label = down_label)
  ggsave(plot = p, filename = glue::glue("{dir}/{prefix}_up_and_down_KEGG.pdf"),height = 5, width = text_w + 3)
  message(glue::glue("Result of enrichKEGG ploted in {dir}"))

}

#' run enrichgesKEGG
#'
#' run enrichgesKEGG and make output files
#'
#' @param dir where to save results files
#' @param prefix a prefix of file names in this step
#' @inheritParams enrich_gsekegg
#' @inheritParams geskegg_barplot
#'
#' @importFrom glue glue
#' @importFrom plyr l_ply
#'
#' @return a list of file
#' @export
#'
#' @examples
#' \dontrun{
#' enrichgesKEGG_run(deg_data = DEG_df, x = "log2FoldChange", dir = tempdir(),eps = 0)
#' }
enrichgesKEGG_run <- function(deg_data,x,dir= ".", prefix = "5-EnrichgseKEGG",
                              pvalue_cut = 0.1, enrichmentScore_cut = 0.5, top = -10,
                              organism = "hsa",
                              keyType = "kegg",
                              OrgDb = 'org.Hs.eg.db',
                              exponent = 1,
                              minGSSize = 10,
                              maxGSSize = 500,
                              eps = 1e-10,
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "BH",
                              verbose = TRUE,
                              use_internal_data = FALSE,
                              seed = FALSE,
                              by = "fgsea") {

  if (!fs::dir_exists(dir)) {
    fs::dir_create(dir)
  }

  gsekegg_res <- enrich_gsekegg(deg_data = deg_data, x = x,
                                organism = organism,
                                keyType = keyType,
                                OrgDb = OrgDb,
                                exponent = exponent,
                                minGSSize = minGSSize,
                                maxGSSize = maxGSSize,
                                eps = eps,
                                pvalueCutoff = pvalueCutoff,
                                pAdjustMethod = pAdjustMethod,
                                verbose = verbose,
                                use_internal_data = use_internal_data,
                                seed = seed,
                                by = by)

  gsekegg_res_file <- glue("{dir}/{prefix}_gsekegg_result.Rdata")
  save(gsekegg_res,file = gsekegg_res_file)
  message(glue::glue("Result of enrich_gsekegg stored in {dir}/{prefix}_{gsekegg_res_file}"))

  write.csv(gsekegg_res@result,glue('{dir}/{prefix}_kegg.gsea.csv'))

  text_w <- max(strwidth(gsekegg_res@result$Description,units = "inch"))

  p = geskegg_barplot(gsekegg_res,pvalue_cut = pvalue_cut, enrichmentScore_cut = enrichmentScore_cut, top = top)
  ggsave(plot = p,filename = glue("{dir}/{prefix}_gseKEGG_barplot.pdf"),height = 5, width = text_w+3)

  message(glue::glue("Barplot of enrichgseKEGG ploted in {dir}"))

  plots_l <- enhance_gseplot(gsekegg_res, top = top,
              pvalue_cut = 0.1, enrichmentScore_cut = 0.5)

  # li_down = structure(plots_l[["down_plots"]], class = c("gglist", "ggplot"))
  # print.gglist = function(x, ...) {plyr::l_ply(x, print, ...)}
  # ggsave(li_down, file = glue("{dir}/{prefix}_gseKEGG_down_gseplot.pdf"))
  #
  # li_up = structure(plots_l[["up_plots"]], class = c("gglist", "ggplot"))
  # ggsave(li_up, file = glue("{dir}/{prefix}_gseKEGG_up_gseplot.pdf"))

  pdf(glue("{dir}/{prefix}_gseKEGG_down_gseplot.pdf"),height = 3,width = 4)
  invisible(lapply(plots_l[["down_plots"]], print))
  dev.off()

  pdf(glue("{dir}/{prefix}_gseKEGG_up_gseplot.pdf"),height = 3,width = 4)
  invisible(lapply(plots_l[["up_plots"]], print))
  dev.off()

  message(glue::glue("Gseplot of enrichgseKEGG ploted in {dir}"))

}


