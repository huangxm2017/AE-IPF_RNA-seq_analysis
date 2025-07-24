#=======================================
# Gene Set Enrichment Analysis (GSEA) Pipeline using fgsea
# Author: YOUR NAME
#=======================================

# Set working directory (modify this path accordingly)
setwd("your/analysis/path")

# Load required packages
library(fgsea)
library(dplyr)
library(ggplot2)
library(data.table)
library(scales)

#---------------------------------------
# Step 1: Load DE result
#---------------------------------------
deg <- read.csv("your_DE_result.csv", check.names = FALSE)
colnames(deg)[1] <- "Gene"
rownames(deg) <- make.unique(deg$Gene)
deg$Gene <- NULL

# Clean data: remove NA
deg <- deg[!is.na(deg$padj) & deg$padj > 0 & !is.na(deg$log2FoldChange), ]

#---------------------------------------
# Step 2: Create weighted score for GSEA
#---------------------------------------
deg$score <- sign(deg$log2FoldChange) * -log10(deg$padj)
rnk_df <- deg[order(deg$score, decreasing = TRUE), ]
rnk_out <- data.frame(Gene = rownames(rnk_df), Score = rnk_df$score)

# Write .rnk style file
write.table(rnk_out, "GSEA_input_weighted.rnk", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

#---------------------------------------
# Step 3: Load gene sets (GMT file)
#---------------------------------------
gmt_file <- "your_pathways_file.gmt"
if (!file.exists(gmt_file)) stop("GMT file not found. Please check the path.")

pathways <- gmtPathways(gmt_file)

#---------------------------------------
# Step 4: Run GSEA
#---------------------------------------
gene_list <- setNames(rnk_out$Score, rnk_out$Gene)
gene_list <- gene_list[!duplicated(names(gene_list))]

fgseaRes <- fgsea(pathways = pathways, stats = gene_list, nperm = 1000)

# Filter significantly enriched pathways
fgsea_sig <- fgseaRes %>%
  filter(padj < 0.05) %>%
  arrange(desc(abs(NES)))

# Save result table
fgsea_sig$leadingEdge <- sapply(fgsea_sig$leadingEdge, paste, collapse = ",")
write.csv(fgsea_sig, "GSEA_significant_pathways.csv", row.names = FALSE)

#---------------------------------------
# Step 5: Bubble plot - NES vs P-adj
#---------------------------------------
ggplot(fgsea_sig, aes(x = NES, y = padj, size = size, color = NES)) +
  geom_point(alpha = 0.85) +
  scale_color_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
  scale_y_reverse(trans = 'log10', labels = label_number(accuracy = 0.001)) +
  labs(title = "GSEA Bubble Plot",
       x = "Normalized Enrichment Score (NES)",
       y = "Adjusted P-value (FDR)",
       size = "Gene Set Size") +
  theme_minimal(base_size = 13)

#---------------------------------------
# Step 6: Bar plot of Top 10 enriched pathways
#---------------------------------------
top10 <- fgsea_sig %>%
  arrange(desc(abs(NES))) %>%
  head(10)

ggplot(top10, aes(x = reorder(pathway, NES), y = NES, fill = NES > 0)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("red", "steelblue")) +
  labs(title = "Top 10 Enriched Pathways",
       x = NULL,
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

#---------------------------------------
# Step 7: Enrichment curve for one pathway
#---------------------------------------
plotEnrichment(pathways[[fgsea_sig$pathway[1]]], gene_list) +
  labs(title = paste("Enrichment Plot:", fgsea_sig$pathway[1]))

#---------------------------------------
# Step 8: Save enrichment plots for top 30 pathways
#---------------------------------------
pdf("GSEA_Top30_Enrichment_Plots.pdf", width = 7, height = 5)

for (pw in fgsea_sig$pathway[1:30]) {
  try({
    print(
      plotEnrichment(pathways[[pw]], gene_list) +
        labs(title = paste("Enrichment:", pw))
    )
  }, silent = TRUE)
}

dev.off()