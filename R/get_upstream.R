#' Recursively identify upstream transcription factors of a target gene from a gene regulatory network,
#' optionally restricting to targets among DEGs.
#'
#' @param target_gene Character. The name of the gene for which upstream TFs are to be identified.
#' @param geneset_GRN Named list where each element is a character vector of target genes regulated by a TF.
#' @param max_depth Integer. Maximum recursion depth (levels of upstream regulation to explore). Default is 2.
#' @param visited Character vector. Keeps track of visited genes to avoid cycles. Default is empty (should not be set by user).
#' @param exclude_self Logical. Whether to exclude self-regulatory edges (TF regulates itself). Default is FALSE.
#' @param deg_genes Optional character vector of DEG gene names to restrict target genes considered. Default is NULL (no restriction).
#'
#' @return A list with:
#'   - upstream_TFs: Character vector of TFs that regulate the target gene directly or indirectly (up to max_depth).
#'   - edges: Data frame with two columns `TF` and `Target` representing the regulatory edges leading to the target gene.
#'
get_upstream <- function(target_gene, geneset_GRN, max_depth = 2, visited = character(), exclude_self = FALSE, deg_genes = NULL) {
  # Stop if max depth reached or gene visited (to avoid cycles)
  if (max_depth == 0 || target_gene %in% visited) {
    return(list(
      upstream_TFs = character(0),
      edges = data.frame(TF = character(0), Target = character(0), stringsAsFactors = FALSE)
    ))
  }

  visited <- c(visited, target_gene)

  # Find TFs that regulate the target_gene
  direct_TFs <- names(Filter(function(targets) target_gene %in% targets, geneset_GRN))

  # If deg_genes specified, filter TFs to only those present in deg_genes
  if (!is.null(deg_genes)) {
    direct_TFs <- intersect(direct_TFs, deg_genes)
  }

  if (length(direct_TFs) == 0) {
    return(list(
      upstream_TFs = character(0),
      edges = data.frame(TF = character(0), Target = character(0), stringsAsFactors = FALSE)
    ))
  }

  edges <- data.frame(TF = direct_TFs, Target = target_gene, stringsAsFactors = FALSE)

  if (exclude_self) {
    edges <- edges[edges$TF != edges$Target, , drop = FALSE]
  }

  upstream_TFs <- direct_TFs

  for (tf in direct_TFs) {
    res <- get_upstream(tf, geneset_GRN, max_depth - 1, visited, exclude_self, deg_genes)
    upstream_TFs <- c(upstream_TFs, res$upstream_TFs)
    edges <- rbind(edges, res$edges)
  }

  list(
    upstream_TFs = unique(upstream_TFs),
    edges = unique(edges)
  )
}
