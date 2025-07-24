# AE-IPF RNA-Seq Analysis in R: DEG Visualization and Functional Enrichment

This repository provides a comprehensive R-based pipeline for differential gene expression (DEG) visualization and pathway enrichment analysis in murine models of Acute Exacerbation of Idiopathic Pulmonary Fibrosis (AE-IPF). It includes scripts for volcano plots, gene set enrichment analysis (GSEA), and GO/KEGG functional enrichment.

## Project Structure

| File                   | Description                                                                                         |
| ---------------------- | --------------------------------------------------------------------------------------------------- |
| `volcano_plot.R`       | Script for generating volcano plots with gradient coloring and gene labeling.                       |
| `gsea_analysis.R`      | Pipeline for GSEA using the `fgsea` package.                                                        |
| `go_kegg_enrichment.R` | Functional enrichment and visualization based on `clusterProfiler` and `simplifyEnrichment`.      |
| `README.md`            | Project documentation and usage instructions.                                                       |

## Key Dependencies

This project uses R (≥ 4.0) and the following packages:
* `ggplot2`
* `ggrepel`
* `dplyr`
* `data.table`
* `fgsea`
* `clusterProfiler`
* `org.Mm.eg.db`
* `enrichplot`
* `KEGGREST`
* `simplifyEnrichment`
* `GO.db`
* `scales`

To install core Bioconductor packages:
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "clusterProfiler", "org.Mm.eg.db", "enrichplot",
  "KEGGREST", "fgsea", "simplifyEnrichment", "GO.db"
))
```

## Scripts Overview

### 1. `volcano_plot.R`
*   **Input**: DEGs table with columns such as `symbol`, `log2FoldChange`, `padj`, `regulate`.
*   **Features**:
    *   Draws a volcano plot with a color gradient based on –log10(padj).
    *   Adds labels for pre-specified genes and the top significant up/downregulated genes.
*   **Output**: Saves the plot as a high-resolution PDF file.

### 2. `gsea_analysis.R`
*   **Input**: DEGs table with log fold changes and adjusted p-values.
*   **Process**:
    1.  Generates a `.rnk` file by combining `log2FC` and `–log10(padj)`.
    2.  Performs GSEA using `fgsea` with local `.gmt` files (e.g., MSigDB gene sets).
*   **Output**:
    *   Table of significant pathways.
    *   Bubble plot (NES vs adjusted p-value).
    *   Barplot of the top 10 pathways.
    *   Enrichment plots for leading pathways (saved to PDF).

### 3. `go_kegg_enrichment.R`
*   **Input**: A list of genes in SYMBOL format.
*   **Process**:
    1.  Converts gene symbols to Entrez IDs.
    2.  Runs `enrichGO()` for Gene Ontology analysis (BP, MF, CC, or ALL).
    3.  Runs `enrichKEGG()` for pathway analysis based on the KEGG REST API.
*   **Visualization**:
    *   Dotplots for GO and KEGG results.
    *   Custom barplots to categorize GO terms based on keywords (e.g., ECM, immune, PI3K-AKT).
    *   Semantic clustering of GO terms using `simplifyEnrichment` based on a similarity matrix.

## Input File Requirements

| Script                 | Required File | Format                                             |
| ---------------------- | ------------- | -------------------------------------------------- |
| `volcano_plot.R`       | DEGs table    | CSV with columns `symbol`, `log2FoldChange`, `padj` |
| `gsea_analysis.R`      | Same as above | `.rnk` file is auto-generated from the input table |
| `go_kegg_enrichment.R` | Gene list     | Single-column CSV with gene SYMBOLs                |

## Expected Output Files

*   `volcano_plot_final.pdf`: Volcano plot with significant gene labels.
*   `GSEA_significant_pathways.csv`: Top GSEA results with NES and padj.
*   `GO_ALL_enrichment_result.csv`, `KEGG_enrichment_result.csv`: Full enrichment outputs.
*   `GO_ALL_dotplot.pdf`, `KEGG_dotplot.pdf`: Dotplot figures.
*   `GO_Category_Barplot.pdf`: Functional group summary plot.
*   `simplifyGO_cluster_BP.pdf` (optional): GO term similarity clustering.

## Notes

*   The pipeline is configured for mouse annotation (`org.Mm.eg.db`). For human studies, replace it with `org.Hs.eg.db`.
*   Make sure to update all `setwd()` paths in the scripts to point to your working directory.
*   Custom keyword matching for GO functional categories can be modified in `go_kegg_enrichment.R`.

## Suggested Workflow

1.  Run `volcano_plot.R` to visualize differential expression results and highlight key genes.
2.  Run `gsea_analysis.R` to identify enriched regulatory pathways from the full ranked gene list.
3.  Use `go_kegg_enrichment.R` with a list of significant DEGs to explore biological processes and visualize GO categories/clustering.

## Contact

Project maintained by Xinmei Huang. For inquiries or contributions, please open an issue on GitHub.
