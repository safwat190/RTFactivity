#' Perform FGSEA on a ranked gene list using GRN-based gene sets
#'
#' This function performs fast Gene Set Enrichment Analysis (FGSEA) using
#' a set of transcription factor target gene sets derived from a gene regulatory network (GRN)
#' and a pre-ranked gene list (e.g., logFC values from DEGs).
#'
#' @param geneset_GRN A list of gene sets (named list). Each element corresponds to a TF and contains a character vector of target genes.
#' @param ranked_list A named numeric vector of gene-level statistics (e.g., logFC). Names must be gene symbols.
#' @param minSize Minimum size of gene sets to test (default: 10).
#' @param maxSize Maximum size of gene sets to test (default: 1000).
#'
#' @return A data frame containing the FGSEA results with columns such as pathway (TF), NES, pval, padj, and leadingEdge.
#' @export
#'
calculate_relative_score <- function(geneset_GRN, ranked_list, minSize = 10, maxSize = 1000) {
  # Sort ranked list in decreasing order
  ranked_list <- sort(ranked_list, decreasing = TRUE)

  # Run fgsea
  fgsea_res <- fgsea::fgsea(
    pathways = geneset_GRN,
    stats = ranked_list,
    minSize = minSize,
    maxSize = maxSize
  )

  # Return result as a data frame
  return(as.data.frame(fgsea_res))
}
