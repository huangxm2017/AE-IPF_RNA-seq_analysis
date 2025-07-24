#==========================================
# GO / KEGG Functional Enrichment Analysis
#==========================================
# Author: Your Name
# Description: This script performs GO and KEGG enrichment 
# analysis and visualization using clusterProfiler and simplifyEnrichment.
# For Mus musculus (mouse), adjust organism/human databases if needed.

#==========================================
# Load or install required packages
#==========================================
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install necessary Bioconductor packages
BiocManager::install(c("clusterProfiler", "org.Mm.eg.db", "enrichplot", "KEGGREST", "simplifyEnrichment", "GO.db"), update = TRUE, ask = FALSE)

# Load libraries
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(KEGGREST)
library(ggplot2)
library(dplyr)
library(simplifyEnrichment)
library(GO.db)

#==========================================
# STEP 1: Load gene list
#==========================================
setwd("your/working/directory")
gene_list_raw <- read.csv("your_gene_list.csv", header = FALSE, stringsAsFactors = FALSE)
gene_symbols <- unique(gene_list_raw$V1)

# Convert SYMBOL to ENTREZID
gene_df <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
cat("Number of successfully mapped genes:", nrow(gene_df), "\n")

#==========================================
# STEP 2: KEGG pathway enrichment
#==========================================
ekegg <- enrichKEGG(
  gene = gene_df$ENTREZID,
  organism = "mmu",
  pvalueCutoff = 0.05,
  keyType = "kegg"
)

#==========================================
# STEP 3: GO enrichment (ALL → can be "BP", "CC", "MF")
#==========================================
ego <- enrichGO(
  gene = gene_df$ENTREZID,
  OrgDb = org.Mm.eg.db,
  ont = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

# Save results
write.csv(as.data.frame(ekegg), "KEGG_enrichment_result.csv", row.names = FALSE)
write.csv(as.data.frame(ego), "GO_ALL_enrichment_result.csv", row.names = FALSE)

#==========================================
# STEP 4: Dotplot visualization (PDF output)
#==========================================
pdf("GO_ALL_dotplot.pdf", width = 9, height = 6)
dotplot(ego, showCategory = 20, title = "GO Enrichment") + theme_minimal()
dev.off()

pdf("KEGG_dotplot.pdf", width = 9, height = 6)
dotplot(ekegg, showCategory = 20, title = "KEGG Pathway Enrichment") + theme_minimal()
dev.off()

#==========================================
# STEP 5: Functional Classification (GO Term Categorization & Barplot)
#==========================================
go_result <- read.csv("GO_BP_enrichment_result.csv")

# Define keyword-based annotation
category_patterns <- list(
  "Cell adhesion" = c("adhesion"),
  "ECM organization" = c("extracellular matrix", "ECM"),
  "Inflammation" = c("inflammatory", "immune", "cytokine"),
  "Fibroblast" = c("fibroblast"),
  "PI3K-AKT" = c("PI3K", "AKT", "mTOR", "phosphatidylinositol")
)

go_result$Category <- "Other"

# Assign category based on matches in description
for (cat in names(category_patterns)) {
  patterns <- category_patterns[[cat]]
  for (p in patterns) {
    matches <- grepl(p, go_result$Description, ignore.case = TRUE)
    go_result$Category[matches] <- cat
  }
}

# Calculate summary
term_summary <- go_result %>%
  group_by(Category) %>%
  summarise(Term_Count = n()) %>%
  mutate(Percentage = round(100 * Term_Count / sum(Term_Count), 1))

# Filter and plot
filtered_summary <- term_summary %>%
  filter(Category != "Other") %>%
  arrange(desc(Term_Count))

ggplot(filtered_summary, aes(x = reorder(Category, -Term_Count), y = Term_Count, fill = Category)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_text(aes(label = paste0(Percentage, "%")), vjust = -0.3, size = 3.5) +
  scale_y_continuous(limits = c(0, 50), expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = "Selected GO Functional Categories",
    x = NULL,
    y = "Number of GO Terms"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10, angle = 30, hjust = 1),
    axis.title.y = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "none"
  )

#==========================================
# STEP 6: Term Clustering with simplifyEnrichment (GO BP, MF, CC)
#==========================================

#---- GO = BP
go_bp <- read.csv("GO_BP_enrichment_result.csv")
go_bp_f <- subset(go_bp, p.adjust < 0.05)
go_ids_bp <- as.character(go_bp_f$ID)
mat_bp <- GO_similarity(go_ids_bp, ont = "BP", db = 'org.Mm.eg.db', measure = "Rel")
df_bp <- simplifyGO(mat_bp)

#---- GO = MF
go_mf <- read.csv("GO_MF_enrichment_result.csv")
go_mf_f <- subset(go_mf, p.adjust < 0.05)
mat_mf <- GO_similarity(as.character(go_mf_f$ID), ont = "MF", db = 'org.Mm.eg.db', measure = "Rel")
df_mf <- simplifyGO(mat_mf)

#---- GO = CC
go_cc <- read.csv("GO_CC_enrichment_result.csv")
go_cc_f <- subset(go_cc, p.adjust < 0.05)
mat_cc <- GO_similarity(as.character(go_cc_f$ID), ont = "CC", db = 'org.Mm.eg.db', measure = "Rel")
df_cc <- simplifyGO(mat_cc)