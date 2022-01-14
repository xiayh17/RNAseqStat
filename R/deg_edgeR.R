#' Run degeR
#'
#' A integrated function for run edgeR in a counts data and return results files.
#'
#' @param counts_data a counts data frame of rows in genes and columns in samples
#' @param group_list a list ordered by samples in counts_data
#' @param dir a directory to store results
#' @param control_group the name of the denominator level for the fold change (Control group)
#' @param x which column is log FC
#' @param y which column is P value
#' @param prefix a prefix of file names in this step
#'
#' @importFrom glue glue
#' @importFrom ggplot2 ggsave
#' @importFrom fs dir_exists dir_create
#'
#' @return a directory contains figures and csv files  and a deg data frame
#' @export
#'
#' @examples
#' deg_edgeR(counts_input,group_list,
#'           control_group = "C",
#'            x = "logFC", y = "PValue",
#'            dir = tempdir(), prefix = "2-DEG_edgeR")
deg_edgeR <- function(counts_data,group_list,
                       control_group,x = "logFC", y = "PValue",
                       dir = ".",prefix = "2-DEG_edgeR") {

  if (!fs::dir_exists(dir)) {
    fs::dir_create(dir)
  }

  deg_data <- run_edgeR(counts_data, group_list, control_group= control_group)

  enhance_heatmap(counts_data, deg_data, group_list, x = x, y = y, dir = dir, prefix = prefix)
  message(glue("{emoji('deciduous_tree')} edgeR heatmap results were store in {dir}."))

  res <- enhance_volcano(deg_data,x = x, y = y,
                         label = c("Down","Stable","Up"), label_ns = "Stable",
                         palette =  c("#2874C5", "grey", "#f87669"),
                         cut_FC = "auto",cut_P = 0.05,top = 10, size = 2.0,expand = c(0.25,0.25),
                         genes_list = "top", highlight = NULL)
  ggsave(res,filename = glue("{dir}/{prefix}_volcano.pdf"), width = 1600,height = 1600,units = "px",limitsize = FALSE)
  message(glue("{emoji('volcano')} edgeR volcano results were store in {dir}."))

  return(deg_data)

}

#' Basic produce of edgeR
#'
#' A basic function to get data and produce results of edgeR
#'
#' @param counts_data a counts data frame of rows in genes and columns in samples
#' @param group_list a character vector ordered by samples in counts_data
#' @param control_group the name of the denominator level for the fold change (Control group)
#'
#' @importFrom stats relevel model.matrix
#' @importFrom edgeR DGEList cpm calcNormFactors estimateGLMCommonDisp estimateGLMTrendedDisp estimateGLMTagwiseDisp glmFit glmLRT topTags
#'
#' @return a DEG data frame
#' @export
#'
#' @examples
#' run_edgeR(counts_input, group_list, control_group= "C")
run_edgeR <- function(counts_data, group_list, control_group) {
  g=factor(group_list)
  g=relevel(g,control_group)

  d <- DGEList(counts=counts_data,group=g)
  keep <- rowSums(cpm(d)>1) >= 2

  d <- d[keep, , keep.lib.sizes=FALSE]
  d$samples$lib.size <- colSums(d$counts)
  d <- calcNormFactors(d)

  dge=d
  design <- model.matrix(~0+g)
  # rownames(design)<-colnames(dge)
  # colnames(design)<-levels(g)

  dge <- estimateGLMCommonDisp(dge,design)
  dge <- estimateGLMTrendedDisp(dge, design)
  dge <- estimateGLMTagwiseDisp(dge, design)

  fit <- glmFit(dge, design)
  # https://www.biostars.org/p/110861/
  lrt <- glmLRT(fit, contrast=c(-1,1))
  nrDEG=topTags(lrt, n=nrow(dge))
  nrDEG=as.data.frame(nrDEG)

  return(nrDEG)

}
