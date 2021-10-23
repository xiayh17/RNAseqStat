#' Choose top genes from DEG data frame
#'
#' For example, if top = 50, head 50 and tail 50 will return from a
#' ordered DEG data frame after filtering by threshold value of P value
#'
#' @param deg_data a DEG data frame contains logFC and p value
#' @param x which column is log FC
#' @param y which column is P value
#' @param top a single number or a length of 2 numeric vector, if 2 numeric vector, first one is top max logFC.
#' @param cut_P a single number for threshold value of P value
#'
#' @importFrom utils head tail
#' @importFrom data.table data.table .SD
#'
#' @return a character vector of top genes
#' @export
#'
#' @examples
#' top_deg(DEG_df, x = "log2FoldChange", y = "pvalue", 50, 0.05)
top_deg <- function(deg_data, x, y, top, cut_P = 0.05) {

  data_f <- deg_data[which(deg_data[,y] <= cut_P),]

  d <- data.table(data_f, key=x, keep.rownames = TRUE)

  if (length(top) == 1) {
    td <- d[, head(.SD, top)]
    hd <-  d[, tail(.SD, top)]
  } else if (length(top) == 2) {
    td <- d[, head(.SD, top[2])]
    hd <-  d[, tail(.SD, top[1])] # the descending order
  } else {
    stop("top should be a single number or a length of 2 numeric vector")
  }

  hd_names <- hd$rn
  td_names <- td$rn

  choose_names <- c(hd_names, td_names)

  return(choose_names)
}

#' cut data
#'
#' group data for plot DEG volcano by FDR and logFC
#'
#' @param deg_data a DEG data frame contains logFC and p value
#' @param x which column is log FC
#' @param y which column is P value
#' @param cut_FC a single number character or numeric vector in threshold value of log FC
#' @param cut_P threshold value of P value, can set for every cut_FC number in numeric vector format
#' @param label symbol word for groups
#' @param label_ns which group is the stable group
#'
#' @return add group column in data frame
#' @export
#'
#' @examples
#' cut_much(DEG_df, x = "log2FoldChange", y = "pvalue", cut_FC = 1,cut_P = 0.05,
#' label = c("Down","Stable","Up"), label_ns = "Stable")
cut_much <- function(deg_data,x, y,cut_FC = 1,cut_P = 0.05,label = c("Down","Stable","Up"), label_ns = "Stable") {

  if (length(cut_FC) == 1) {
    cut_FC <- c(-cut_FC, cut_FC)
  }
  if (length(cut_P) == 1) {
    cut_P <- rep(cut_P,length(cut_FC))
  }

  label_cg <- base::setdiff(label,label_ns)
  names(cut_P) <- label_cg
  deg_data$group <- cut(deg_data[,x],
                        breaks = c(-Inf,cut_FC,Inf),
                        labels = label)
  index = list()
  for (i in label_cg) {
    index[[i]] <- setdiff(which(deg_data$group == i), which(deg_data$group == i & deg_data[,y] < cut_P[i]))
    deg_data$group[index[[i]]] <- label_ns
  }
  message(paste(label, as.character(table(deg_data$group)),"\n"))
  return(deg_data)
}
