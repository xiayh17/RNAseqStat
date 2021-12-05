#' run all step in default setting
#'
#' Only need count data, group, OrgDb, dir
#'
#' @param count_data a counts data frame of rows in genes and columns in samples
#' @param group_list a list ordered by samples in counts_data
#' @param OrgDb OrgDb
#' @param dir a directory to store results
#' @param test_group which one is test group in your group_list
#' @param control_group which one is control group in your group_list
#' @param parallel if FALSE, no parallelization. if TRUE, parallel execution using BiocParallel
#'
#' @return a dir contains all results
#' @export
#'
#' @examples
#' \dontrun{
#' runAll(count_data = counts_input, group_list = group_list, OrgDb = 'org.Hs.eg.db', dir = tempdir())
#' }
runAll <- function(count_data, group_list, OrgDb = 'org.Hs.eg.db', dir = ".",test_group = "T", control_group = "C",parallel = FALSE) {

  message(glue("Step1: Check you data"))
  pre_check(counts_data = count_data, group_list = group_list, dir = dir)

  message(glue("Step2: DEG analysis"))
  deg_results <- deg_run(count_data, group_list, test_group = test_group, control_group = control_group,dir = dir,parallel = parallel)

  message(glue("Step3: EnrichGO analysis"))
  enrichGO_run(deg_results@deg_df_limma, x = "logFC", y = "P.Value",
               OrgDb = OrgDb, dir = dir,prefix = "3-EnrichGO-limma")
  enrichGO_run(deg_results@deg_df_edgeR, x = "logFC", y = "PValue",
               OrgDb = OrgDb, dir = dir,prefix = "3-EnrichGO-edgeR")
  enrichGO_run(deg_results@deg_df_DESeq2, x = "log2FoldChange", y = "pvalue",
               OrgDb = OrgDb, dir = dir,prefix = "3-EnrichGO-DESeq2")

  message(glue("Step4: EnrichKEGG analysis"))
  enrichKEGG_run(deg_results@deg_df_limma, x = "logFC", y = "P.Value",
                 OrgDb = OrgDb, dir = dir,prefix = "4-EnrichKEGG-limma")
  enrichKEGG_run(deg_results@deg_df_edgeR, x = "logFC", y = "PValue",
                 OrgDb = OrgDb, dir = dir,prefix = "4-EnrichKEGG-edgeR")
  enrichKEGG_run(deg_results@deg_df_DESeq2, x = "log2FoldChange", y = "pvalue",
                 OrgDb = OrgDb, dir = dir,prefix = "4-EnrichKEGG-DESeq2")

  message(glue("Step5: EnrichgseKEGG analysis"))
  enrichgesKEGG_run(deg_results@deg_df_limma, x = "logFC", OrgDb = OrgDb,
                    dir = dir,prefix = "5-EnrichgseKEGG-limma")
  enrichgesKEGG_run(deg_results@deg_df_edgeR, x = "logFC", OrgDb = OrgDb,
                    dir = dir,prefix = "5-EnrichgseKEGG-edgeR")
  enrichgesKEGG_run(deg_results@deg_df_DESeq2, x = "log2FoldChange", OrgDb = OrgDb,
                    dir = dir,prefix = "5-EnrichgseKEGG-DESeq2")

}
