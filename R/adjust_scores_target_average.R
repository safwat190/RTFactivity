#' Adjust TF scores and return direct, indirect, and combined scores as a data frame
#'
#' This function computes three versions of transcription factor (TF) scores:
#' - Original TF scores (alpha = 1),
#' - Mean scores of direct target TFs (alpha = 0),
#' - Weighted combination of the two based on input alpha.
#'
#' The resulting data frame includes TF names as a column, with row names set to simple sequential numbers.
#'
#' @param gsea_res Named numeric vector of original TF scores (e.g., NES scores from GSEA).
#'                  The names should be TF names.
#' @param geneset_GRN Named list where each element contains a character vector of direct target genes for each TF.
#' @param alpha Numeric value between 0 and 1 specifying the weight of the original TF score in the combined score.
#'              - alpha = 1 means no adjustment (original scores only),
#'              - alpha = 0 means score is fully replaced by mean target score.
#'
#' @return Data frame with columns:
#'   - TF: transcription factor names,
#'   - direct.alpha.1: original TF scores,
#'   - indirect.alpha.0: mean scores of direct target TFs,
#'   - combined.alpha.of.interest: weighted combination of original and target scores.
#'   Row names are numeric sequential indices (1, 2, 3, ...) instead of TF names.
adjust_scores_target_average <- function(gsea_res, geneset_GRN, alpha = 0.5) {
  # Initialize numeric vectors to store scores, preserving TF names
  direct_scores <- numeric(length(gsea_res))
  indirect_scores <- numeric(length(gsea_res))
  combined_scores <- numeric(length(gsea_res))

  names(direct_scores) <- names(gsea_res)
  names(indirect_scores) <- names(gsea_res)
  names(combined_scores) <- names(gsea_res)

  # Loop over each TF to compute scores
  for (tf in names(gsea_res)) {
    # Original TF score (alpha = 1)
    direct_scores[tf] <- gsea_res[tf]

    # Get direct targets of the TF
    direct_targets <- geneset_GRN[[tf]]
    # Keep only targets that are TFs with available scores
    tf_targets <- intersect(direct_targets, names(gsea_res))

    # Mean score of direct target TFs, or 0 if none
    if (length(tf_targets) > 0) {
      indirect_scores[tf] <- mean(gsea_res[tf_targets])
    } else {
      indirect_scores[tf] <- 0
    }

    # Weighted combination of original and mean target scores
    combined_scores[tf] <- alpha * direct_scores[tf] + (1 - alpha) * indirect_scores[tf]
  }

  # Combine all scores into a data frame
  df <- data.frame(
    TF = names(gsea_res),
    direct.alpha.1 = direct_scores,
    indirect.alpha.0 = indirect_scores,
    combined.alpha.of.interest = combined_scores,
    stringsAsFactors = FALSE
  )

  # Set row names to sequential numbers (1, 2, 3, ...)
  rownames(df) <- seq_len(nrow(df))

  return(df)
}
