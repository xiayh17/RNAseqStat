#' Run limma
#'
#' A integrated function for run limma in a counts data and return results files.
#'
#' @param counts_data a counts data frame of rows in genes and columns in samples
#' @param group_list a list ordered by samples in counts_data
#' @param dir a directory to store results
#' @param test_group the name of the numerator level for the fold change (Test group)
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
#' deg_limma(counts_input,group_list,
#'           test_group = "T", control_group = "C",
#'            x = "logFC", y = "P.Value",
#'            dir = tempdir(), prefix = "2-DEG_limma")
deg_limma <- function(counts_data,group_list,
                       test_group,control_group,x,y,
                       dir = ".",prefix = "2-DEG_limma") {

  if (!fs::dir_exists(dir)) {
    fs::dir_create(dir)
  }

  deg_data <- run_limma(counts_data = counts_data,group_list = group_list,
                         test_group = test_group,control_group=control_group)

  enhance_heatmap(counts_data, deg_data, group_list, x = x, y = y, dir = dir, prefix = prefix)
  message(glue("{emoji('deciduous_tree')} limma heatmap results were store in {dir}."))

  res <- enhance_volcano(deg_data,x = x, y = y,
                         label = c("Down","Stable","Up"), label_ns = "Stable",
                         palette =  c("#2874C5", "grey", "#f87669"),
                         cut_FC = "auto",cut_P = 0.05,top = 10, size = 2.0,expand = c(0.25,0.25),
                         genes_list = "top", highlight = NULL)
  ggsave(res,filename = glue("{dir}/{prefix}_volcano.pdf"), width = 1600,height = 1600,units = "px",limitsize = FALSE)
  message(glue("{emoji('volcano')} limma volcano results were store in {dir}."))

  return(deg_data)

}

#' Basic produce of limma
#'
#' A basic function to get data and produce results of limma
#'
#' @param counts_data a counts data frame of rows in genes and columns in samples
#' @param group_list a character vector ordered by samples in counts_data
#' @param test_group the name of the numerator level for the fold change (Test group)
#' @param control_group the name of the denominator level for the fold change (Control group)
#'
#' @importFrom stats model.matrix na.omit
#' @importFrom edgeR DGEList cpm calcNormFactors
#' @importFrom limma voom lmFit makeContrasts contrasts.fit eBayes topTable
#'
#' @return a DEG data frame
#' @export
#'
#' @examples
#' run_limma(counts_input, group_list, control_group= "C", test_group = "T")
run_limma <- function(counts_data, group_list, control_group, test_group) {

  design <- model.matrix(~0+factor(group_list))
  colnames(design)=levels(factor(group_list))
  rownames(design)=colnames(counts_data)

  dge <- DGEList(counts=counts_data)
  dge <- calcNormFactors(dge)
  logCPM <- cpm(dge, log=TRUE, prior.count=3)

  v <- voom(dge,design,plot=TRUE, normalize.method="quantile")
  fit <- lmFit(v, design)

  con=paste0(test_group,'-',control_group)

  cont.matrix=makeContrasts(contrasts=c(con),levels = design)
  fit2=contrasts.fit(fit,cont.matrix)
  fit2=eBayes(fit2)

  tempOutput = topTable(fit2, coef=con, number=Inf)
  DEG_limma_voom = na.omit(tempOutput)

  return(DEG_limma_voom)

}
