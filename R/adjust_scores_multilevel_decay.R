#' Adjust transcription factor scores by propagating influence through multilevel targets with exponential decay
#'
#' @param gsea_res Named numeric vector of original TF scores (e.g., from GSEA).
#' @param geneset_GRN A named list where each element is a character vector of target genes regulated by a TF. The names of the list are TFs.
#' @param max_depth Integer. Maximum depth of target levels to propagate scores through. Default is 3.
#' @param decay_factor Numeric between 0 and 1. Decay factor applied exponentially by depth level to weight downstream TF contributions. Default is 0.5.
#'
#' @return Data frame with columns:
#'   - original.score: original TF scores,
#'   - adjusted.score: TF scores adjusted by incorporating multilevel target influences weighted by decay.
#'
#' @examples
#' gsea_scores <- c(TF1 = 2, TF2 = 1.5, TF3 = 3)
#' grn <- list(TF1 = c("TF2", "GeneA"), TF2 = c("TF3"), TF3 = c("GeneB"))
#' adjust_scores_multilevel_decay(gsea_scores, grn, max_depth = 2, decay_factor = 0.6)
adjust_scores_multilevel_decay <- function(gsea_res, geneset_GRN, max_depth = 3, decay_factor = 0.5) {

  # Initialize adjusted scores with the original TF scores
  adjusted_scores <- gsea_res

  # Recursive helper function to get all targets of a TF at a given depth level
  get_targets_at_depth <- function(tf, depth, visited = character()) {
    # Stop recursion if max depth reached, TF already visited (to avoid cycles), or no targets for TF
    if (depth == 0 || tf %in% visited || is.null(geneset_GRN[[tf]])) return(character(0))

    # Mark current TF as visited
    visited <- c(visited, tf)

    # Get direct targets of the current TF
    direct_targets <- geneset_GRN[[tf]]

    # Recursively get targets of targets at depth-1 for those targets that are TFs themselves
    deeper_targets <- unlist(lapply(direct_targets, function(tgt) {
      if (tgt %in% names(geneset_GRN)) {
        get_targets_at_depth(tgt, depth - 1, visited)
      } else character(0)
    }))

    # Return unique combination of direct targets and recursively found deeper targets
    unique(c(direct_targets, deeper_targets))
  }

  # Iterate over each TF in the original scores vector
  for (tf in names(gsea_res)) {
    additive_score <- 0

    # Only proceed if max_depth is at least 1 (to consider targets)
    if (max_depth >= 1) {
      # For each depth level from 1 to max_depth
      for (depth in 1:max_depth) {
        # Retrieve all targets of the TF at this depth
        targets <- get_targets_at_depth(tf, depth)

        # Filter targets to keep only those that are TFs with available scores
        tf_targets <- intersect(targets, names(gsea_res))

        if (length(tf_targets) > 0) {
          # Extract the scores of these target TFs
          target_scores <- gsea_res[tf_targets]

          # Normalize the sum of scores by the number of target TFs at this depth
          normalized_score <- sum(target_scores) / length(tf_targets)

          # Apply exponential decay weighting by depth and add to the additive score
          additive_score <- additive_score + normalized_score * (decay_factor ^ depth)
        }
      }
    }

    # Update the adjusted score by adding the weighted sum of downstream TF scores
    adjusted_scores[tf] <- adjusted_scores[tf] + additive_score
  }

  # Combine original and adjusted scores into a data frame
  df <- data.frame(
    TF = names(gsea_res),
    original.score = gsea_res,
    adjusted.score = adjusted_scores,
    stringsAsFactors = FALSE
  )

  # Optionally reset row names to sequential numbers
  rownames(df) <- seq_len(nrow(df))

  return(df)
}
