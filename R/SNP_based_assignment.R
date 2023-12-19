# Calculate a dissimilarity matrix for Hamming distances between the SNPs of the reference
# and the deconvoluted single-cell clusters
#' Dissimilarity Matrix of binarized SNP matrix
#'
#' @param matrix1 Binarized SNP matrix from demultiplexed single-cell RNA-seq experiments
#' @param matrix2 Binarized SNP matrix from called SNPs of bulkRNA
#'
#' @return Dissimilarity Matrix
#' @export
#'
#' @examples
SNP_based_assignment <- function(matrix1, matrix2) {
  if (nrow(matrix1) != nrow(matrix2)) {
    stop("SNP Matrices must have the same number of patients and demultiplexed clusters")
  }

  num_rows <- nrow(matrix1)
  dissimilarity_matrix <- matrix(0, nrow = num_rows, ncol = num_rows)

  # Get row names from the original matrices
  rownames_matrix1 <- rownames(matrix1)
  rownames_matrix2 <- rownames(matrix2)

  for (i in 1:num_rows) {
    row1 <- matrix1[i, , drop = FALSE]

    for (j in 1:num_rows) {
      row2 <- matrix2[j, , drop = FALSE]

      # Ensure both rows have the same length
      if (length(row1) != length(row2)) {
        stop("Rows must have the same length")
      }

      # Calculate Hamming distance for the current pair of rows
      dissimilarity_matrix[i, j] <- sum(row1 != row2)
    }
  }
  # Set row and column names for the dissimilarity matrix
  rownames(dissimilarity_matrix) <- rownames_matrix1
  colnames(dissimilarity_matrix) <- rownames_matrix2

  # Hungarian algorithm for assignment
  lsap_result <- clue::solve_LSAP(dissimilarity_matrix)

  assignment_df <- data.frame(
    "Reference" = colnames(dissimilarity_matrix)[lsap_result],
    "Cluster" = rownames(dissimilarity_matrix)
  )

  result <- list(Dissimilarity_Matrix = dissimilarity_matrix, Assignments = assignment_df)
  return(result)
}
