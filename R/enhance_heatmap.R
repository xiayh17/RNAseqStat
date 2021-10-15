#' Heatmap for DEG data frame
#'
#' default will return a top 100 deg heatmap in p value = 0.05
#'
#' @param counts_data a counts data frame of rows in genes and columns in samples
#' @param x which column is log FC
#' @param y which column is P value
#' @param deg_data a DEG data frame contains logFC and p value
#' @param group_list a character vector ordered by samples in counts_data
#'
#' @importFrom edgeR cpm
#' @importFrom pheatmap pheatmap
#' @importFrom glue glue
#'
#' @return a heatmap plot file
#' @export
#'
#' @examples
#' enhance_heatmap(counts_input, deg_data, group_list, x = "log2FoldChange", y = "pvalue")
enhance_heatmap <- function(counts_data, deg_data, group_list, x, y, top = 50, cut_P = 0.05, prefix = "2-DEG") {
  choose_gene <- top_deg(deg_data,x = x, y = y, top = top, cut_P = cut_P)
  exprSet=log(edgeR::cpm(counts_data)+1)
  choose_matrix=exprSet[choose_gene,]
  choose_matrix=t(scale(t(choose_matrix)))
  choose_matrix[choose_matrix>2]=2
  choose_matrix[choose_matrix< -2]= -2
  colD=data.frame(Groups=group_list)
  rownames(colD)=colnames(exprSet)
  pheatmap(choose_matrix,
           annotation_col = colD,fontsize = 12,
           width = 400/100*2,
           height = 550/100*3,
           filename = glue('{prefix}_top100_heatmap.pdf'))
}


