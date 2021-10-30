#' Enrichment analysis of go
#'
#' run enrich_go and plot results
#'
#' @inheritParams enrich_go
#' @param top top rows
#' @param dir where to save results files
#'
#' @importFrom utils write.table
#'
#' @return a list result files of GO
#' @export
#'
#' @examples
#' \dontrun{
#' enrichGO_run(DEG_df, x = "log2FoldChange", y = "pvalue", dir = tempdir())
#' }
enrichGO_run <- function(deg_data, x, y, cut_FC = 1, cut_P = 0.05, top = 10, dir = ".",
                       OrgDb = 'org.Hs.eg.db', keyType = "SYMBOL", ont = "ALL", simplify = TRUE,
                       pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2,minGSSize = 10,
                       maxGSSize = 500, readable = FALSE, pool = FALSE,
                       label = c("Down", "Stable", "Up"),
                       label_ns = "Stable",
                       mc.cores = 1L) {

  go_resl <- enrich_go(deg_data, x, y, cut_FC = cut_FC, cut_P = cut_P,
            OrgDb = OrgDb, keyType = keyType, ont = ont, simplify = simplify,
            pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod, qvalueCutoff = qvalueCutoff,minGSSize = minGSSize,
            maxGSSize = maxGSSize, readable = readable, pool = pool,
            label = label,
            label_ns = label_ns,
            mc.cores = 1L)

  go_res_file <- glue("{dir}/go_result.Rdata")
  save(go_resl,file = go_res_file)
  message(glue::glue("Result of enrichGO stored in {dir}/{go_res_file}"))

  tmp <- lapply(seq_along(go_resl), function(x)
    write.table(go_resl[[x]]@result,file = glue::glue('{dir}/gene_{names(go_resl)[[x]]}_GO_enrichment.csv'))
  )

  if (ont == "ALL") {
    plots <- lapply(go_resl, function(x)
      enhance_barplot(x@result,top=top,split = ONTOLOGY)
    )
    lapply(seq_along(plots), function(x)
      ggsave(plot = plots[[x]], filename = glue::glue("{dir}/barplot-gene_{names(plots)[[x]]}_GO_enrichment.pdf"),
             height = 880*2.5/300, width = 780*2.5/300, dpi = 300)
    )
  } else {
    plots <- lapply(go_resl, function(x)
      enhance_barplot(x@result,top=top)
    )
    lapply(seq_along(plots), function(x)
      ggsave(plot = plots[[x]], filename = glue::glue("{dir}/barplot-gene_{names(plots)[[x]]}_GO_enrichment.pdf"),
             height = 880*2.5/300/3, width = 780*2.5/300, dpi = 300)
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
#'
#' @return a list result files of KEGG
#' @export
#'
#' @examples
#' \dontrun{
#' enrichKEGG_run(deg_data = DEG_df, x = "log2FoldChange", y = "pvalue", cut_FC = 1,
#' cut_P = 0.05, top = 10)
#' }
enrichKEGG_run <- function(deg_data, x, y, cut_FC = 1, cut_P = 0.05, top = 10, dir= ".",
                           organism = "hsa",
                           keyType = "kegg",
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
  # enrich_kegg
  kegg_resl <- enrich_kegg(deg_data = deg_data, x = x, y = y, cut_FC = 1,
                           cut_P = cut_P,
                           organism = organism,
                           keyType = keyType,
                           pvalueCutoff = pvalueCutoff,
                           pAdjustMethod = pAdjustMethod,
                           minGSSize = minGSSize,
                           maxGSSize = maxGSSize,
                           qvalueCutoff = qvalueCutoff,
                           use_internal_data = use_internal_data,
                           label = label,
                           label_ns = label_ns,
                           mc.cores = mc.cores)
  p <- kegg_barplot(kegg_resl, top = top, down_label = down_label)
  ggsave(plot = p, filename = glue::glue("{dir}/up_and_down_KEGG.pdf"),height = 7.01,width = 6.7)

}
