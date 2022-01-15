#' Run DEseq2
#'
#' A integrated function for run DEseq2 in a counts data and return results files.
#'
#' @param counts_data a counts data frame of rows in genes and columns in samples
#' @param group_list a list ordered by samples in counts_data
#' @param dir a directory to store results
#' @param test_group the name of the numerator level for the fold change (Test group)
#' @param control_group the name of the denominator level for the fold change (Control group)
#' @param qc qc plots
#' @param x which column is log FC
#' @param y which column is P value
#' @param prefix a prefix of file names in this step
#' @param parallel if FALSE, no parallelization. if TRUE, parallel execution using BiocParallel
#'
#' @importFrom glue glue
#' @importFrom ggplot2 ggsave
#' @importFrom fs dir_exists dir_create
#'
#' @return a directory contains figures and csv files and a deg data frame
#' @export
#'
#' @examples
#' \dontrun{
#' deg_DESeq2(counts_input,group_list,
#'           test_group = "T", control_group = "C", qc = TRUE,
#'            x = "log2FoldChange", y = "pvalue",
#'            dir = tempdir(), prefix = "2-DEG_DEseq2")
#' }
deg_DESeq2 <- function(counts_data,group_list, parallel = F,
                       test_group,control_group,qc = TRUE,x,y,
                       dir = ".",prefix = "2-DEG_DEseq2") {

  if (!fs::dir_exists(dir)) {
    fs::dir_create(dir)
  }

    deg_data <- run_DESeq2(counts_data = counts_data,group_list = group_list, parallel = parallel,
                           test_group = test_group,control_group=control_group,qc = qc,dir = dir,prefix = prefix)

    enhance_heatmap(counts_data, deg_data, group_list, x = x, y = y, dir = dir, prefix = prefix)
    message(glue("{emoji('deciduous_tree')} DESeq2 heatmap results were store in {dir}."))

    res <- enhance_volcano(deg_data,x = x, y = y,
                    label = c("Down","Stable","Up"), label_ns = "Stable",
                    palette =  c("#2874C5", "grey", "#f87669"),
                    cut_FC = "auto",cut_P = 0.05,top = 10, size = 2.0,expand = c(0.25,0.25),
                    genes_list = "top", highlight = NULL)
    ggsave(res,filename = glue("{dir}/{prefix}_volcano.pdf"), width = 1600,height = 1600,units = "px",limitsize = FALSE)
    message(glue("{emoji('volcano')} DESeq2 volcano results were store in {dir}."))

    return(deg_data)

}



#' Basic produce of DESeq2
#'
#' A basic function to get data and produce results of DESeq2
#'
#' @param counts_data a counts data frame of rows in genes and columns in samples
#' @param group_list a character vector ordered by samples in counts_data
#' @param test_group the name of the numerator level for the fold change (Test group)
#' @param control_group the name of the denominator level for the fold change (Control group)
#' @param qc qc plots
#' @param dir a directory to store results
#' @param prefix a prefix of file names in this step
#' @param parallel if FALSE, no parallelization. if TRUE, parallel execution using BiocParallel
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#' @importFrom stats na.omit
#'
#' @return a DEG data frame
#' @export
#'
#' @examples
#' \dontrun{
#' run_DESeq2(counts_input, group_list,test_group = "T", control_group = "C", dir = tempdir())
#' }
run_DESeq2 <- function(counts_data,group_list,test_group,control_group,qc = TRUE,dir = ".",prefix = "2-DEG_DEseq2",parallel = FALSE) {

  colData <- data.frame(row.names=colnames(counts_data),
                         group_list=group_list)

  colData$group_list <- factor(colData$group_list)

  dds <- DESeqDataSetFromMatrix(countData = counts_data,
                                colData = colData,
                                design = ~ group_list)
  dds <- DESeq(dds)

  if (qc) {
    DESeq2_qc(counts_data,dds,dir = dir,prefix = prefix)
  }

  res <- results(dds,parallel = parallel,
                 contrast=c("group_list",test_group,control_group))
  resOrdered <- res[order(res$padj),]

  DEG =as.data.frame(resOrdered)
  DEG = na.omit(DEG)
  return(DEG)
}

#' QC for DESeq2
#'
#' plot dispersions and RAWvsNORM
#'
#' @param counts_data a counts data frame of rows in genes and columns in samples
#' @param dds a DESeqDataSet class data set
#'
#' @importFrom glue glue
#' @importFrom grDevices pdf dev.off rainbow
#' @importFrom graphics par boxplot hist
#' @importFrom DESeq2 plotDispEsts rlogTransformation varianceStabilizingTransformation
#' @importMethodsFrom SummarizedExperiment assay
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#'
#' @return dispersions and RAWvsNORM firures
#'
#' @noRd
#' @examples
#' \dontrun{
#' DESeq2_qc(counts_input,dds, dir = tempdir())
#' }
DESeq2_qc <- function(counts_data, dds, dir = ".", prefix = "2-DEG") {

  pdf(glue("{dir}/{prefix}_qc_dispersions.pdf"), 18, 18, pointsize=35)
  plotDispEsts(dds, main="Dispersion plot",cex = 0.45)
  dev.off()

  if (length(rownames(dds@colData)) >=50 ) {
    rld <- varianceStabilizingTransformation(dds)
  } else {

    rld <- rlogTransformation(dds)

  }

  exprMatrix_rlog=assay(rld)

  pdf(glue("{dir}/{prefix}_RAWvsNORM.pdf"),height = 8,width = 8)
  par(cex = 0.7)
  n.sample=ncol(counts_data)
  if(n.sample>40) par(cex = 0.5)
  cols <- rainbow(n.sample*1.2)
  par(mfrow=c(2,2))
  boxplot(counts_data, col = cols,main="expression value",las=2)
  boxplot(exprMatrix_rlog, col = cols,main="expression value",las=2)
  hist(as.matrix(counts_data))
  hist(exprMatrix_rlog)
  dev.off()
}

