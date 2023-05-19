
library(DESeq2)
library(ggplot2)

# Mock TPM counts data
counts <- data.frame(
    Sample1 = c(500, 1000, 2000, 1500, 3000),
    Sample2 = c(1000, 800, 1500, 1200, 2500),
    Sample3 = c(600, 1200, 1800, 1600, 2800),
    Sample4 = c(800, 900, 1400, 1300, 2600),
    Sample5 = c(400, 700, 2200, 1700, 3200),
    row.names = c("Gene1", "Gene2", "Gene3", "Gene4", "Gene5")
)

# Mock sample metadata table with time points
metadata <- data.frame(
    sample = colnames(counts),
    time_point = c(0, 1, 2, 3, 4)
)

# Transpose counts data
counts <- t(counts)

# Create DESeqDataSet object with time point as a covariate
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ time_point)

# Perform differential expression analysis
dds <- DESeq(dds)

# Extract differential expression results
results <- results(dds)

# Create a data frame for plotting
plot_data <- as.data.frame(results)

# Add log fold change and negative log10 p-value columns
plot_data$log2FoldChange <- results$log2FoldChange
plot_data$neg_log10_padj <- -log10(results$padj)

# Volcano plot
volcano_plot <- ggplot(plot_data, aes(x = log2FoldChange, y = neg_log10_padj)) +
    geom_point(size = 3, alpha = 0.8) +
    xlim(c(-5, 5)) +
    ylim(c(0, 10)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    labs(x = "log2 Fold Change", y = "-log10(adjusted p-value)")

print(volcano_plot)


# Visualize the results using a volcano plot
volcano_plot <- ggplot(results, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(size = 3, alpha = 0.8) +
    xlim(c(-5, 5)) +
    ylim(c(0, 10)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    labs(x = "log2 Fold Change", y = "-log10(adjusted p-value)")

print(volcano_plot)
