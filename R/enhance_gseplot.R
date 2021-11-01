#' plot gseKEGG plot
#'
#' plot gseKEGG plots by up and down
#'
#' @param data output from enrich_gseKEGG
#' @param top filter top by pvalue
#' @param pvalue_cut filter cut of pvalue
#' @param enrichmentScore_cut filter cut of enrichmentScore
#'
#' @importFrom enrichplot gseaplot2
#'
#' @return a list contains up and down
#' @export
#'
#' @examples
#' \dontrun{
#' gsekegg_res <- enrich_gsekegg(DEG_df,x = "log2FoldChange")
#' plots_l <- enhance_gseplot(gsekegg_res)
#' }
enhance_gseplot <- function(data, top = 10,
                            pvalue_cut = 0.1, enrichmentScore_cut = 0.5) {

  down_kegg<-data[data$pvalue< pvalue_cut & data$enrichmentScore < -enrichmentScore_cut,]
  up_kegg<-data[data$pvalue< pvalue_cut & data$enrichmentScore > enrichmentScore_cut,]

  # down_kegg <- down_kegg[order(down_kegg[,"pvalue"])[1:top],]
  # up_kegg <- up_kegg[order(up_kegg[,"pvalue"])[1:top],]
  down_kegg <- down_kegg %>% top_n(top,wt = pvalue)
  up_kegg <- up_kegg %>% top_n(top,wt = pvalue)

  plot_gseplot <- function(data,data_ud,x) {
    gseaplot2(data,data_ud$ID[x],
              title=data_ud$Description[x],pvalue_table = FALSE,subplots = 1:3,color = "blue") +
      labs(caption = paste(sprintf('Pvalue:  %.3f _ P.adjust: %.3f', data_ud$pvalue[x], data_ud$p.adjust[x]))) +
      theme(plot.caption = element_text(family = "Times", size = 8, face = "italic", colour = "dodgerblue"))
  }

  down_plots <- lapply(1:nrow(down_kegg), function(x)
    plot_gseplot(data,down_kegg,x)
      )

  up_plots <- lapply(1:nrow(up_kegg), function(x)
    plot_gseplot(data,up_kegg,x)
  )

  resl <- list(up_plots, down_plots)

  names(resl) <- c("up_plots", "down_plots")

  return(resl)

}
