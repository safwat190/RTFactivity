#' Recursively identify downstream target genes of a transcription factor from a gene regulatory network,
#' optionally restricting targets to DEGs.
#'
#' @param tf_gene Character. The name of the transcription factor (or gene) for which downstream targets are to be identified.
#' @param geneset_GRN Named list where each element is a character vector of target genes regulated by a TF.
#' @param max_depth Integer. Maximum recursion depth (levels of downstream targets to explore). Default is 2.
#' @param visited Character vector. Keeps track of visited genes to avoid cycles. Default is empty (should not be set by user).
#' @param exclude_self Logical. Whether to exclude self-regulatory edges (TF regulates itself). Default is FALSE.
#' @param deg_genes Optional character vector of DEG gene names to restrict downstream targets considered. Default is NULL (no restriction).
#'
#' @return A list with:
#'   - downstream_targets: Character vector of target genes regulated directly or indirectly by the TF (up to max_depth).
#'   - edges: Data frame with two columns `TF` and `Target` representing regulatory edges downstream from the TF.
#'
get_downstream <- function(tf_gene, geneset_GRN, max_depth = 2, visited = character(), exclude_self = FALSE, deg_genes = NULL) {
  # Stop recursion if max depth reached, gene already visited, or no targets for this gene
  if (max_depth == 0 || tf_gene %in% visited || is.null(geneset_GRN[[tf_gene]])) {
    return(list(
      downstream_targets = character(0),
      edges = data.frame(TF = character(0), Target = character(0), stringsAsFactors = FALSE)
    ))
  }

  visited <- c(visited, tf_gene)

  # Get direct targets regulated by this TF
  direct_targets <- geneset_GRN[[tf_gene]]

  # If deg_genes specified, filter direct_targets to only those in deg_genes
  if (!is.null(deg_genes)) {
    direct_targets <- intersect(direct_targets, deg_genes)
  }

  if (length(direct_targets) == 0) {
    return(list(
      downstream_targets = character(0),
      edges = data.frame(TF = character(0), Target = character(0), stringsAsFactors = FALSE)
    ))
  }

  edges <- data.frame(TF = tf_gene, Target = direct_targets, stringsAsFactors = FALSE)

  if (exclude_self) {
    edges <- edges[edges$TF != edges$Target, , drop = FALSE]
  }

  downstream_targets <- direct_targets

  for (target in direct_targets) {
    res <- get_downstream(target, geneset_GRN, max_depth - 1, visited, exclude_self, deg_genes)
    downstream_targets <- c(downstream_targets, res$downstream_targets)
    edges <- rbind(edges, res$edges)
  }

  list(
    downstream_targets = unique(downstream_targets),
    edges = unique(edges)
  )
}
