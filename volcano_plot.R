## Volcano Plot with Gradient Color and Annotated Genes
## (c) Author, Year
## This script generates a volcano plot with gradient coloring, significance thresholds, and label annotations for specific genes.

# Load required packages
library(ggplot2)
library(ggrepel)
library(tidyverse)

# Set working directory (TO BE MODIFIED BY USER)
setwd("your/working/path/here")

# Load DESeq2-like results
data <- read.csv("your_data_file.csv", header = TRUE, check.names = FALSE)

# Auto-fill rownames from first column and remove it from data frame
row.names(data) <- make.names(data[[1]], unique = TRUE)
data <- data[, -1]

# Check variable format (optional debugging)
class(data$Log2FoldChange)
class(data$padj)

# Add gene symbol and empty label column
data$symbol <- rownames(data)
data$label <- NA

# Example: manually assign labels to specific genes you care about
data$label[which(data$symbol == "YourGene1")] <- "YourGene1"

# Automatically label top genes (adjust N as needed)
# Top upregulated genes
up_data <- data %>%
  filter(regulate == 'up') %>%
  distinct(symbol, .keep_all = TRUE) %>%
  top_n(10, -log10(padj))

# Top downregulated genes
down_data <- data %>%
  filter(regulate == 'down') %>%
  distinct(symbol, .keep_all = TRUE) %>%
  top_n(10, -log10(padj))

# Basic volcano plot with gradient color
g <- ggplot(data, aes(x = Log2FoldChange, y = -log10(padj))) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey60") +
  geom_point(aes(size = -log10(padj), color = -log10(padj))) +
  scale_color_gradientn(values = c(0, 0.1, 0.2, 0.5, 1),
                        colors = c("#39489f", "#39bbec", "#f9ed36", "#f38466", "#b81f25")) +
  scale_size_continuous(range = c(0.5, 3)) +
  theme_bw() +
  theme(panel.grid = element_blank())

# Add manually specified labels using geom_label_repel
g_label <- g +
  geom_label_repel(data = data, aes(label = label),
                   size = 4,
                   box.padding = unit(0.5, "lines"),
                   point.padding = unit(0.8, "lines"),
                   segment.color = "red",
                   show.legend = FALSE,
                   max.overlaps = 10000)

# Add top up/downregulated genes using symbol annotation
g_label_full <- g +
  geom_label_repel(data = up_data,
                   aes(x = Log2FoldChange, y = -log10(padj), label = symbol),
                   size = 3,
                   segment.color = "grey50",
                   force = 5,
                   force_pull = 1,
                   seed = 123) +
  geom_label_repel(data = down_data,
                   aes(x = Log2FoldChange, y = -log10(padj), label = symbol),
                   size = 3,
                   segment.color = "grey50",
                   force = 5,
                   force_pull = 1,
                   seed = 123) +
  coord_cartesian(xlim = c(-10, 10), ylim = c(0, 50)) +
  scale_x_continuous(breaks = seq(-10, 10, by = 2)) +
  scale_y_continuous(breaks = seq(0, 50, by = 10), limits = c(0, 50)) +
  labs(title = "Volcano Plot with Gene Annotations",
       x = "Log2 Fold Change",
       y = "-log10(Adjusted P-value)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(3, "pt")
  )

# Export as PDF
ggsave("volcano_plot_final.pdf", plot = g_label_full,
       width = 8, height = 6, dpi = 300)