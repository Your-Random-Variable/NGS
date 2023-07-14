library(DESeq2)
library(tidyverse)

# Read in count data and metadata
count_data <- read.table("counts_all.tsv", header=TRUE, row.names=1)
metadata <- read.table("meta.txt", header=TRUE)

# Convert condition column to a factor
metadata$condition <- factor(metadata$condition)

# Create DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = metadata, design = ~condition)

# Run DESeq2
dds <- DESeq(dds)

# PCA plot of samples
plotPCA(rlog(dds), intgroup="condition")

# Plot of dispersion estimates
plotDispEsts(dds)

# Filter out low count genes
dds <- dds[ rowSums(counts(dds)) > 10, ]

# Run DESeq2 on filtered dataset
dds <- DESeq(dds)

# Get differential expression results
res <- results(dds)

# Annotate results with gene names
gene_names <- rownames(count_data)
gene_names_filt <- rownames(counts(dds))
res <- data.frame(gene_name = gene_names_filt, res)
# Subset rows with complete cases
res <- res[complete.cases(res), ]

# Write annotated results to file
write.table(res, file = "DE_genes_annotated.tsv", sep = "\t", quote = TRUE, row.names = FALSE)

# Plot of log2 fold changes vs. mean expression
ggplot(res, aes(x = baseMean, y = log2FoldChange)) +
  geom_point(alpha = 0.5) +
  scale_x_log10() +
  ggtitle("DE Genes") +
  xlab("Mean Expression") +
  ylab("Log2 Fold Change")

# Plot of p-value vs. log2 fold change
ggplot(res, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(alpha = 0.5) +
  ggtitle("DE Genes") +
  xlab("Log2 Fold Change") +
  ylab("-Log10 P-value")

# Subset results with p-value < 0.05
sig_res <- res[res$pvalue < 0.05, ]
# Subset rows with complete cases
sig_res <- sig_res[complete.cases(sig_res), ]
write.table(sig_res, file = "SIG.tsv", sep = "\t", quote = TRUE, row.names = FALSE)


# Subset up-regulated genes
sig_up <- sig_res[sig_res$log2FoldChange > 0,]

# Write up-regulated genes to file
write.table(sig_up, file = "up.tsv", sep = "\t", quote = TRUE, row.names = FALSE)

# Subset down-regulated genes
sig_down <- sig_res[sig_res$log2FoldChange < 0,]

# Write down-regulated genes to file
write.table(sig_down, file = "down.tsv", sep = "\t", quote = TRUE, row.names = FALSE)
