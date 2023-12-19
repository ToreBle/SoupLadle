#Replace the variant with a binary string
#' Binarization of the SNP profile
#'
#' @param matrix SNP matrix that will be binarized
#'
#' @return binarized SNP matrix
#' @export
#'
#' @examples
Binarize_genotype_values <- function(matrix) {
  matrix[matrix == "0/0"] <- "100"
  matrix[matrix == "0/1"] <- "010"
  matrix[matrix == "1/1"] <- "001"
  matrix[matrix == "0/2"] <- "011"
  matrix[matrix == "0/3"] <- "011"
  matrix[is.na(matrix)] <- "000"

  return(matrix)
}
