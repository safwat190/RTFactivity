# Function to summarize number of targets of TFs
#' Count number of targets per transcription factor and overlap with DEGs
#'
#' @param targets_list A named list where each element is a character vector of gene targets for a transcription factor (TF). The names of the list should be TF names.
#' @param deg_genes Optional character vector of differentially expressed genes (DEGs). If provided, the function will count how many targets of each TF are in this DEG list.
#'
#' @return A data frame with:
#'   - TF: the name of the transcription factor
#'   - n_targets: total number of targets for that TF
#'   - n_targets_in_DEGs (optional): number of targets that are also DEGs
#'   The data frame is sorted in descending order by `n_targets_in_DEGs` if DEGs are provided, or by `n_targets` otherwise. Row numbers are reset.
#'
count_targets <- function(targets_list, deg_genes = NULL) {
  # Create a data frame with TF names and number of targets
  tf_target_counts <- data.frame(
    TF = names(targets_list),
    n_targets = sapply(targets_list, length),
    stringsAsFactors = FALSE
  )

  # If DEGs are provided, count how many targets are also DEGs
  if (!is.null(deg_genes)) {
    tf_target_counts$n_targets_in_DEGs <- sapply(targets_list, function(targets) {
      sum(targets %in% deg_genes)
    })

    # Sort by number of DEG-overlapping targets
    tf_target_counts <- tf_target_counts[order(-tf_target_counts$n_targets_in_DEGs), ]
  } else {
    # Sort by total number of targets if no DEG info
    tf_target_counts <- tf_target_counts[order(-tf_target_counts$n_targets), ]
  }

  # Reset row names to avoid mismatched indices after sorting
  row.names(tf_target_counts) <- NULL

  return(tf_target_counts)
}
