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

    # Load example data
    data("example_degs")
    data("example_grn")

    # Run TF activity estimation
    results <- compute_TF_activity(degs = example_degs, grn = example_grn)

    # View results
    head(results)

## 📂 Input Format

Your input should include:

    DEGs table (data frame):

        gene: gene symbol

        logFC: log fold-change between condition A vs B

        adj.pval: adjusted p-value

    Gene Regulatory Network (GRN) (data frame):

        TF: transcription factor name

        Target: gene regulated by TF

        (optional) score: edge weight or confidence
        

## 📈 Output

The function compute\_TF\_activity() returns a data frame with the
following columns:

    TF: transcription factor name

    direct_score: effect from immediate (first-order) target DEGs

    cascade_score: effect from downstream (second-order) targets

    n_direct_targets: number of direct targets used

    n_cascade_targets: number of cascade targets used

## 🌍 Applications

RTFactivity was developed to identify TFs driving beneficial cell state
transitions across:

    ✅ Wound healing (e.g., reprogramming fibroblasts to promote angiogenesis)

    ✅ Diabetic foot ulcers (healing vs non-healing)

    ✅ Beta cell differentiation

    ✅ Spinal cord regeneration

## 📚 Citation

Ahmed S. Abouhashem, et al. (2024). Identification of Skin Multicellular
Reprogramming Factors as Potential Treatment for Non-Healing Diabetic
Foot Ulcers. \[Under Revision\] RTFactivity GitHub:
<https://github.com/safwat190/RTFactivity>

## 🧪 Vignette

For a full analysis pipeline, see the vignette:

    vignette("RTFactivity-intro")

## 🤝 Contributing

Pull requests and feedback are welcome. Please open an issue or submit
changes via PR.

## 🧠 Maintainer

Dr. Ahmed S. Abouhashem, Assistant Professor, University of Pittsburgh
📧 <A.S.A@pitt.edu>
