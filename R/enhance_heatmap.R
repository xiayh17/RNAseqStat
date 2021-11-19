#' Heatmap for DEG data frame
#'
#' default will return a top 100 deg heatmap in p value = 0.05
#'
#' @param counts_data a counts data frame of rows in genes and columns in samples
#' @param x which column is log FC
#' @param y which column is P value
#' @param cut_P threshold value of P value, can set for every cut_FC number in numeric vector format
#' @param top a single number or a length of 2 numeric vector, if 2 numeric vector, first one is top max logFC.
#' @param deg_data a DEG data frame contains logFC and p value
#' @param group_list a character vector ordered by samples in counts_data
#' @param dir a directory to store results
#' @param prefix a prefix of file names in this step
#' @param palette a color palette for plots
#'
#'
#' @importFrom edgeR cpm
#' @importFrom pheatmap pheatmap
#' @importFrom glue glue
#'
#' @return a heatmap plot file
#' @export
#'
#' @examples
#' enhance_heatmap(counts_input, DEG_df, group_list,
#'     x = "log2FoldChange", y = "pvalue", dir = tempdir())
enhance_heatmap <- function(counts_data, deg_data, group_list, x, y, top = 50, cut_P = 0.05, dir = ".", prefix = "2-DEG", palette = RColorBrewer::brewer.pal(3,"Set2")[1:2]) {

  choose_gene <- top_deg(deg_data,x = x, y = y, top = top, cut_P = cut_P)

  if (is.na(dir)) {
    filename = NA
  } else {
    filename = glue('{dir}/{prefix}_top{top*2}_heatmap.pdf')
  }

  exprSet=log(edgeR::cpm(counts_data)+1)
  choose_matrix=exprSet[choose_gene,]
  choose_matrix=t(scale(t(choose_matrix)))
  choose_matrix[choose_matrix>2]=2
  choose_matrix[choose_matrix< -2]= -2
  colD=data.frame(Groups=group_list)
  rownames(colD)=colnames(exprSet)
  names(palette) <- unique(group_list)
  pheatmap(choose_matrix,
           annotation_col = colD,fontsize = 12,
           width = (ncol(choose_matrix)*0.3+2.2) *2,
           height = 550/100*3*nrow(choose_matrix)/100,
           annotation_colors = list(Groups = palette),
           filename = filename)
}


