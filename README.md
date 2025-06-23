# RTFactivity

**RTFactivity** is an R package to estimate the **Relative Transcription
Factor (TF) Activity** between two biological conditions using gene
regulatory networks (GRNs) and differentially expressed genes (DEGs). It
computes both **direct** and **cascade (indirect)** effects of
transcription factors, enabling discovery of key regulators that drive
**cell state transitions**, such as those involved in wound healing,
cell reprogramming, or regeneration.

------------------------------------------------------------------------

## ✨ Key Features

-   Calculates **direct TF scores** based on immediate targets  
-   Calculates **cascade TF scores** through second-order targets
    (targets-of-targets)  
-   Compatible with any user-provided GRN (global or cell
    type–specific)  
-   Accepts standard DE results (Seurat, DESeq2, edgeR, etc.)  
-   Designed for applications in **regeneration**, **differentiation**,
    and **disease modeling**

------------------------------------------------------------------------

## 📦 Installation

    install.packages("devtools")
    devtools::install_github("safwat190/RTFactivity")

## 🚀 Quick Start

    library(RTFactivity)

    # Load example data:
    # example_DEGs: differentially expressed genes in fibroblasts between healing and nonhealing ulcer)
    # example_GRN: fibroblast gene regulatory network (GRN)
    data("example_DEGs")
    data("example_GRN")

    # convert GRN to gene sets
    geneset_GRN <- split(example_GRN$Target, example_GRN$TF)

    # order DEGs in descending order
    example_DEGs <- example_DEGs[order(-example_DEGs$avg_log2FC), ]

    # prepare ranked gene list from DEGs
    ranked_list = example_DEGs$avg_log2FC
    names(ranked_list) = example_DEGs$gene
    ranked_list <- sort(ranked_list, decreasing = TRUE)

    # check number of targets for each TF which are among DEGs or not
    tf_summary <- count_targets(geneset_GRN, deg_genes = example_DEGs$gene)
    head(tf_summary)
    hist(tf_summary$n_targets, breaks = 30)
    hist(tf_summary$n_targets_in_DEGs, breaks = 30)
    plot(tf_summary$n_targets, tf_summary$n_targets_in_DEGs)
    summary(tf_summary$n_targets_in_DEGs)

    # calculate the relative activity of each TF in fibroblasts between healing and non-healing DFU 
    gsea_res_GRN = calculate_relative_score(geneset_GRN = geneset_GRN,
                                            ranked_list = ranked_list,
                                            minSize = 10,
                                            maxSize = 1000)
    View(gsea_res_GRN)

    # Get TFs with adjusted p value < 0.05 
    gsea_res_GRN_sig = subset(gsea_res_GRN, subset = padj < 0.05)

    sig_TFs <- gsea_res_GRN %>%
      filter(padj < 0.05) %>%
      arrange(desc(NES)) %>%
      pull(pathway)

    # Plot heatmap for top 5 with positive and negative scores
    # get TF names
    top_TF = 5
    sig_TFs = c(head(sig_TFs,top_TF),tail(sig_TFs,top_TF))

    # Create a dataframe to store their targets among DEGs as binaries
    targets_df <- data.frame(row.names = example_DEGs$gene)

    for (tf in sig_TFs) {
      targets_df[[tf]] <- ifelse(rownames(targets_df) %in% geneset_GRN[[tf]], 1, 0)
    }
    head(targets_df)
    tail(targets_df)

    # annotate upregulated and downregulated genes in the heatmap
    annotation_row <- data.frame(
      direction = ifelse(example_DEGs$avg_log2FC > 0, "positive", "negative")
    )
    rownames(annotation_row) <- example_DEGs$gene

    pheatmap(targets_df,
             cluster_rows = F,
             cluster_cols = F,
             color = colorRampPalette(c("white", "firebrick3"))(50),
             angle_col = 90,
             annotation_row = annotation_row,
             na_col = "grey90",
             border_color = 'white',
             show_rownames = FALSE
    )

    # calculate indirect scores based on targets of targets
    adjusted_scores <- adjust_scores_target_average(
      gsea_res = setNames(gsea_res_GRN_sig$NES, gsea_res_GRN_sig$pathway),
      geneset_GRN = geneset_GRN,
      alpha = 0.5
    )

    # plot direct against indirect scores
    adjusted_scores$sign_category <- with(adjusted_scores, ifelse(
      direct.alpha.1 > 0 & indirect.alpha.0 > 0, "both_positive",
      ifelse(direct.alpha.1 < 0 & indirect.alpha.0 < 0, "both_negative", "non_consistent")
    ))

    ggplot(adjusted_scores, aes(x = direct.alpha.1, y = indirect.alpha.0, color = sign_category)) +
      geom_point(size = 3, alpha = 0.8) +          
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
      scale_color_manual(
        values = c(
          both_positive = "red",
          both_negative = "darkgreen",
          non_consistent = "grey"
        )
      ) +
      labs(
        title = "Scatter plot of Direct vs Indirect TF Scores",
        x = "Direct TF Score (alpha = 1)",
        y = "Indirect TF Score (alpha = 0)",
        color = "Sign Category",
        caption = "Source: Your adjusted_scores dataframe"
      ) +
      theme_minimal(base_size = 15) +                                
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title = element_text(),
        panel.grid.minor = element_blank()
      )

    # View results
    head(adjusted_scores)

    # Get direct downstream targets for the TF NR3C1
    res <- get_downstream("NR3C1", geneset_GRN, max_depth = 1, exclude_self = T)
    res$edges
    res$downstream_targets

    # Get indirect downstream targets for the TF NR3C1 at 2 levels
    res <- get_downstream("NR3C1", geneset_GRN, max_depth = 2, exclude_self = T)
    res$edges
    res$downstream_targets

    # Get direct upstream targets for the TF NR3C1
    res <- get_upstream("NR3C1", geneset_GRN, max_depth = 1, exclude_self = T)
    res$edges
    res$upstream_TFs

    # Get indirect upstream targets for the TF NR3C1 at 2 levels
    res <- get_upstream("NR3C1", geneset_GRN, max_depth = 2, exclude_self = T)
    res$edges
    res$upstream_TFs

## 📂 Input Format

Your input should include:

    DEGs table (data frame):

        gene: gene symbol

        logFC: log fold-change between condition A vs B

        adj.pval: adjusted p-value

    Gene Regulatory Network (GRN) (data frame):

        TF: transcription factor name

        Target: gene regulated by TF
        

## 📈 Output

The function calculate\_relative\_score() returns a data frame with the
following columns:

    TF: transcription factor name

    direct.alpha.1: TF score based on direct targets only

    indirect.alpha.0: TF score based on indirect targets only

    combined.alpha.of.interest: weighted average TF score based on alpha (between 0 and 1)

## 🌍 Applications

RTFactivity was developed to identify TFs driving beneficial cell state
transitions such as in:

    ✅ Wound healing (e.g., reprogramming fibroblasts to promote angiogenesis)

    ✅ Beta cell differentiation

    ✅ Spinal cord regeneration

## 📚 Citation

Ahmed S. Abouhashem, et al. (2025). Identification of Skin Multicellular
Reprogramming Factors as Potential Treatment for Non-Healing Diabetic
Foot Ulcers. \[Under Revision\] RTFactivity GitHub:
<https://github.com/safwat190/RTFactivity>

Korotkevich G, Sukhov V, Sergushichev A (2019). “Fast gene set
enrichment analysis.” bioRxiv. <doi:10.1101/060012>

Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L.,
Gillette, M. A., Paulovich, A., Pomeroy, S. L., Golub, T. R., Lander, E.
S., & Mesirov, J. P. (2005). Gene set enrichment analysis: a
knowledge-based approach for interpreting genome-wide expression
profiles. Proceedings of the National Academy of Sciences of the United
States of America, 102(43), 15545–15550.
<https://doi.org/10.1073/pnas.0506580102>

## 🤝 Contributing

Pull requests and feedback are welcome. Please open an issue or submit
changes via PR.

## 🧠 Maintainer

Dr. Ahmed S. Abouhashem  
Research Assistant Professor  
University of Pittsburgh  
📧 <A.S.A@pitt.edu>
