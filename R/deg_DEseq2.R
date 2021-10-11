#' Pre-check you data by PCA and Heat maps
#'
#' A previous check function for overview the data sets.
#'
#' @param counts_data a counts data frame of rows in genes and columns in samples
#' @param group_list a list ordered by samples in counts_data
#' @param dir a directory to store results
#'
#' @importFrom glue glue
#' @importFrom edgeR cpm
#' @importFrom fs dir_exists dir_create
#'
#' @return a directory contains figures for data QC check
#' @export
#'
#' @examples
#' pre_check(counts_input, group_list)
# deg_DESeq2 <- function(counts_data,group_list,dir = ".") {
#
# }
#
# run_DESeq2 <- function(counts_data,group_list) {
#
# }



