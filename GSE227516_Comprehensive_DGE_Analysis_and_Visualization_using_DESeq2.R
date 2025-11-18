################################################################################
#
# GSE227516: Exercise activates AMPK in mouse pancreatic islet leading 
#            to decreased senescence
#
# Comprehensive DGE Analysis and Visualization using DESeq2
# (Version with automatic plot saving)
#
################################################################################


## 1. Define Output Path ##
# Define the path where all plots and results will be saved
OUTPUT_PATH <- "D:/DOWNLOADS/"

## 2. Create Output Directory ##
# Create the directory if it doesn't already exist
if (!dir.exists(OUTPUT_PATH)) {
  dir.create(OUTPUT_PATH, recursive = TRUE)
  print(paste("Created directory:", OUTPUT_PATH))
} else {
  print(paste("Output directory already exists:", OUTPUT_PATH))
}

## 3. Install BiocManager ##
# (if not already installed)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

## 4. Install DESeq2 ##
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2")
}

## 5. Install RUVSeq ##
if (!requireNamespace("RUVSeq", quietly = TRUE)) {
  BiocManager::install("RUVSeq")
}

## 6. Install pheatmap ##
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}

## 7. Install RColorBrewer ##
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  install.packages("RColorBrewer")
}

## 8. Install ggplot2 ##
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

## 9. Install ggrepel ##
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel")
}

## 10. Load Packages ##
library(DESeq2)
library(RUVSeq)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)

################################################################################

## 11. Import Gene Counts ##
# Note: Update this file path to your local directory
COUNTS <- read.csv(file="D:/DOWNLOADS/GSE227516_counts.csv", header=TRUE, row.names=1)
head(COUNTS)

## 12. Import Sample Information (Metadata) ##
# Note: Update this file path to your local directory
META <- read.csv(file="D:/DOWNLOADS/sample_information.csv", header=TRUE)
head(META)

################################################################################

## 13. Reorder Count Matrix Columns ##
# This ensures the order of columns in COUNTS matches the order of rows in META
COUNTS <- COUNTS[, c("P1", paste0("P", 2:9), "P10")]
head(COUNTS)

## 14. Format Counts to Integer Matrix ##
# DESeq2 requires non-negative integers
COUNTS <- round(COUNTS)
COUNTS <- as.matrix(COUNTS)

## 15. Check Unique Conditions ##
unique(META$condition)

################################################################################

## 16. Create DESeqDataSet Object ##
dds <- DESeqDataSetFromMatrix(countData = COUNTS, 
                              colData=META, 
                              design=~condition)
dim(dds)

## 17. Pre-filter Low Count Genes ##
# Removes genes with a row mean count less than 10
threshold <- 10
dds <- dds[ rowMeans(counts(dds)) >= threshold,]
dim(dds) # Check new dimensions after filtering

## 18. Run DESeq Analysis ##
# This single function runs:
# 1. estimateSizeFactors()
# 2. estimateDispersions()
# 3. nbinomWaldTest()
prdds <- DESeq(dds)
prdds

################################################################################

## 19. Get Normalized Counts ##
norm_counts <- counts(prdds, normalized = TRUE)
norm_counts <- as.data.frame(norm_counts)
head(norm_counts)

## 20. Apply Transformations (rld/vsd) ##
# These transformed values are for QC plots (PCA, heatmaps), 
# NOT for the differential expression testing itself.
mks <- estimateSizeFactors(dds) # Needed for log2(x+1) plot
rld <- rlogTransformation(prdds, blind = FALSE)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

################################################################################

## 21. Save QC Scatter Plots (log2 vs rlog vs vst) ##
png(file.path(OUTPUT_PATH, "GSE227516_QC_Scatter_Plots.png"), width=1800, height=600, res=150)
par(mfrow=c(1, 3))
lims <- c(-2, 20)
plot(log2(counts(mks, normalized=TRUE)[,1:2] + 1),pch=16, cex=0.3, main="log2(x + 1)", xlim=lims, ylim=lims)
plot(assay(rld)[,1:2], pch=16, cex=0.3, main="R log (rld)", xlim=lims, ylim=lims)
plot(assay(vsd)[,1:2], pch=16, cex=0.3, main="VST (vsd)", xlim=lims, ylim=lims)
dev.off() # Close the PNG device
par(mfrow=c(1, 1)) # Reset plot window

## 22. Save QC Histograms (Transformed Counts) ##
png(file.path(OUTPUT_PATH, "GSE227516_QC_Histograms.png"), width=1800, height=600, res=150)
par(mfrow=c(1, 3))
hist(counts(mks), main="Raw Counts")
hist(assay(rld), main="R log (rld)")
hist(assay(vsd), main="VST (vsd)")
dev.off() # Close the PNG device
par(mfrow=c(1, 1)) # Reset plot window

## 23. Save Sample-to-Sample Distance Heatmap ##
sample_dist <- dist(t(assay(rld)))
sample_dist_matrix <- as.matrix(sample_dist)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pheatmap(sample_dist_matrix,
         clustering_distance_rows=sample_dist,
         clustering_distance_cols=sample_dist,
         col=colors,
         main = "Sample-to-Sample Distance (rld)",
         filename = file.path(OUTPUT_PATH, "GSE227516_QC_Sample_Distance_Heatmap.png"))

## 24. Save Principal Component Analysis (PCA) Plot ##
pca_data <- plotPCA(rld, intgroup = c("condition"), returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
  geom_point(size = 2, aes(color = condition)) +
  geom_text_repel(aes(label = rownames(pca_data)), nudge_x = 0, nudge_y = 0) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("GSE227516: PCA of Samples (rld)")

print(pca_plot) # Show plot in R
ggsave(file.path(OUTPUT_PATH, "GSE227516_QC_PCA_Plot.png"), 
       plot = pca_plot, 
       width = 8, 
       height = 6, 
       dpi = 300)

## 25. Save Dispersion Estimates Plot ##
png(file.path(OUTPUT_PATH, "GSE227516_QC_Dispersion_Plot.png"), width=800, height=600, res=100)
plotDispEsts(prdds, main = "Dispersion plot", 
             genecol="gray20", fitcol="red", 
             finalcol="dodgerblue3" 
) 
dev.off() # Close the PNG device

################################################################################

## 26. Extract and Order Results ##
# alpha = 0.05 sets the adjusted p-value (padj) cutoff for significance
res05 <- results(prdds, alpha = 0.05)
res05 <- na.omit(res05) # Remove genes with NA values

res05ordered <- res05[order(res05$padj),]
head(as.data.frame(res05ordered))

## 27. Save MA Plot ##
png(file.path(OUTPUT_PATH, "GSE227516_Results_MA_Plot.png"), width=800, height=600, res=100)
DESeq2::plotMA(
  res05, 
  main=paste("GSE227516: Sedentary vs Exercise (alpha=0.05)"), 
  ylim=c(-5,10),
  cex=0.5, 
  colNonSig=adjustcolor("gray20", alpha.f=0.5), 
  colSig=adjustcolor("dodgerblue3", alpha.f=0.5) 
)
abline(h = 1, col = '#ff0000' , lwd = 1)
abline(h = -1, col= '#ff0000', lwd = 1)
dev.off() # Close the PNG device

## 28. Add Gene Status for Volcano Plot ##
res05$gene_status <- ifelse(
  res05$padj < 0.05, 
  ifelse(
    res05$log2FoldChange > 1, 
    "Up-Regulated",
    ifelse(
      res05$log2FoldChange < -1, 
      "Down-Regulated", 
      "Non-significant"
    )
  ), 
  "Non-significant"
)

## 29. Save Volcano Plot ##
volcano_plot <- ggplot(
  as.data.frame(res05), # Convert to data.frame for ggplot
  aes(x = log2FoldChange, y = -log10(padj), color = factor(gene_status))
) +
  geom_point(size = 1, alpha = 0.7) +
  scale_color_manual(values = c("Up-Regulated" = brewer.pal(3, "Set1")[1],
                                "Down-Regulated" = brewer.pal(3, "Set1")[2],
                                "Non-significant" = "grey")) +
  theme_minimal() +
  ggtitle("GSE227516: Volcano Plot of Differentially Expressed Genes") +
  xlab("log2 Fold Change") +
  ylab("-log10(Adjusted p-value)") +
  theme(legend.title = element_blank()) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed")

print(volcano_plot) # Show plot in R
ggsave(file.path(OUTPUT_PATH, "GSE227516_Results_Volcano_Plot.png"),
       plot = volcano_plot,
       width = 8,
       height = 6,
       dpi = 300)

################################################################################

## 30. Filter Significant Genes ##
# (padj < 0.05 & |LFC| > 1)
sig_genes <- as.data.frame(res05[res05$padj < 0.05 & abs(res05$log2FoldChange) > 1, ])
head(sig_genes)

## 31. Create Subsets for Up/Down Genes ##
up_genes <- subset(sig_genes, log2FoldChange > 0)
down_genes <- subset(sig_genes,log2FoldChange < 0)
head(up_genes)
head(down_genes)

## 32. Save Heatmap of Top 10 Up-Regulated Genes ##
top_up <- head(up_genes[order(up_genes$log2FoldChange, decreasing = TRUE), ], 10)
top_up_exp <- assay(rld)[rownames(top_up), ] # Get rlog-transformed counts

pheatmap(top_up_exp,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         scale = "row", # Scale expression by row (gene)
         show_colnames = TRUE,
         col=brewer.pal(name="RdBu", n=11),
         main = "GSE227516: Top 10 Up-Regulated Genes (rld counts)",
         filename = file.path(OUTPUT_PATH, "GSE227516_Results_Heatmap_Top10_Up.png"))

## 33. Save Heatmap of Top 10 Down-Regulated Genes ##
top_down <- head(down_genes[order(down_genes$log2FoldChange, decreasing = FALSE), ], 10) # 'decreasing = FALSE' for most negative
top_down_exp <- assay(rld)[rownames(top_down), ]

pheatmap(top_down_exp,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         scale = "row",
         show_colnames = TRUE,
         col=brewer.pal(name="RdBu", n=11),
         main = "GSE227516: Top 10 Down-Regulated Genes (rld counts)",
         filename = file.path(OUTPUT_PATH, "GSE227516_Results_Heatmap_Top10_Down.png"))

################################################################################

## 34. Summarize DE Gene Counts ##
n_total <- nrow(COUNTS)
n_de_genes <- nrow(sig_genes)
n_up_genes <- nrow(up_genes)
n_down_genes <- nrow(down_genes)

print(paste("Total genes analyzed:", n_total))
print(paste("Total significant DE genes:", n_de_genes))
print(paste("Up-regulated genes:", n_up_genes))
print(paste("Down-regulated genes:", n_down_genes))

## 35. Save Significant DE Genes to CSV (Downloads Path) ##
# Add gene status and gene ID as columns for the final file
sig_genes$gene_status <- ifelse(sig_genes$log2FoldChange > 0, "Up-Regulated", "Down-Regulated")
sig_genes$gene_id <- rownames(sig_genes)

# This line saves the CSV to your D:\DOWNLOADS path
write.table(sig_genes, 
            file = file.path(OUTPUT_PATH, "GSE227516_significant_DE_genes.csv"), 
            sep = ",", 
            row.names = FALSE)
print(paste("Successfully saved significant genes list to:", file.path(OUTPUT_PATH, "GSE227516_significant_DE_genes.csv")))

################################################################################

## 36. Record Session Information ##
# Records R version and package versions for reproducibility
sessionInfo()

################################################################################
# End of Script
################################################################################