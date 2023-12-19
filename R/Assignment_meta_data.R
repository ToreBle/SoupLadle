#' Bulk Assignments to meta data
#'
#' @param assignment_df Assignemnt data.frame from the SNP_based_assignment
#' @param clusters Demultiplexing file of the scRNA-seq data (e.g. cluster.tsv from SouporCell)
#' @param output_csv Logical of wheter an output file will be generated
#' @param output_file Path for the output csv file that can be used as meta data in any scRNA-seq processing object
#'
#' @return meta data for Seurat
#' @export
#'
#' @examples

Assignment_meta_data <- function(assignment_df, clusters, output_csv = TRUE, output_file = "SoupLadle_Assignments.csv") {
  colnames(clusters)[colnames(clusters) == "assignment"] <- "SouporCell_Assignment"
  bulk_assignments <- assignment_df$Reference[match(clusters$SouporCell_Assignment, assignment_df$Cluster)]
  clusters$SoupLadle <- ifelse(is.na(bulk_assignments), clusters$status, bulk_assignments)

  if (output_csv) {
    write.csv(clusters, file = output_file)
  }

  return(clusters)
}
