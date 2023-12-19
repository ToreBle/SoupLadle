#' Read vcf files
#'
#' @param vcf_file Path to vcf file to load (e.g. bulkRNA-seq or scRNA-seq variants)
#' @param gt element to extract from vcf genotype data. Common options include 'DP', 'GT' and 'GQ'.
#' @param common_SNP_filter Filter out all SNPs that have a common genotype in all samples or clusters
#'
#' @return A tibble of the loaded vcf file
#' @export
#'
#' @examples
#'
read_extract_vcf <- function(vcf_file, gt = NULL, common_SNP_filter = TRUE) {
  cat(paste("Processing the file: ", vcf_file, "\n", sep=""))
  VCF <- vcfR::read.vcfR(vcf_file)

  if (is.null(gt)) {
    gt <- vcfR::extract.gt(VCF, element = 'GT')
  }

  if (common_SNP_filter) {
    gt <- filter_rows(gt)
    cat(paste("Filtered variant: ", nrow(gt), sep=""))
  }

  return(gt)
}

# Filter out SNPs that are common in all samples
filter_rows <- function(dataframe) {
  # Check for rows where more than one non-NA value exists
  discard_rows <- apply(dataframe, 1, function(row) {
    # Filter out NA values
    non_na_values <- na.omit(row)
    # Check if more than one non-NA value exists
    length(unique(non_na_values)) <= 1
  })

  # Subset the dataframe to keep rows that are not common for all samples
  filtered_df <- dataframe[!discard_rows | rowSums(!is.na(dataframe)) <= 1, ]

  return(filtered_df)
}


