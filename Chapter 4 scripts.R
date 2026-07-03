## Load featureCounts matrix

# Input file: RNA-seq featureCounts output
# Note: raw count data is stored locally and is not included in this repository.

counts_file <- "data/featurecounts.txt"

# Import RNA-seq featureCounts output
counts_raw <- readr::read_tsv(
  counts_file,
  comment = "#",
  show_col_types = FALSE
)

# Import sample metadata
metadata_file <- "data/metadata.tsv"

meta <- read.delim(
  metadata_file,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)



## Prepare featureCounts matrix

# featureCounts output contains gene annotation columns followed by sample count columns.
# This section extracts the sample counts and creates a gene-by-sample count matrix.

gene_col <- "Geneid"
stopifnot(gene_col %in% names(counts_raw))

annot_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
sample_cols <- setdiff(names(counts_raw), intersect(names(counts_raw), annot_cols))

counts_mat <- counts_raw %>%
  dplyr::select(dplyr::all_of(c(gene_col, sample_cols))) %>%
  tibble::column_to_rownames(gene_col) %>%
  as.matrix()

storage.mode(counts_mat) <- "integer"

# Check dimensions and sample names if needed
# dim(counts_mat)
# colnames(counts_mat)



## Library size quality check

# This plot was used as an initial quality check to compare sequencing depth across samples.

# Calculate total counts per sample
lib_sizes <- colSums(counts_mat)
summary(lib_sizes)

# Plot library sizes in millions
lib_sizes_millions <- lib_sizes / 1e6

plot(
  lib_sizes_millions,
  las = 2,
  ylab = "Total counts per sample (millions)",
  xlab = "Sample",
  main = "Library sizes from featureCounts"
)



## Raw count distribution quality check

# Log-transform raw counts for visualisation only
log_counts <- log2(counts_mat + 1)

# Boxplot of count distributions across samples
boxplot(
  log_counts,
  las = 2,
  outline = FALSE,
  ylab = "log2(count + 1)",
  main = "Raw count distributions per sample"
)



## Count matrix sanity checks

# Check how sparse the count matrix is
zero_frac_by_gene <- rowMeans(counts_mat == 0)
summary(zero_frac_by_gene)

# Check for duplicated gene IDs
sum(duplicated(rownames(counts_mat)))

# Check that counts are non-negative integers
any(counts_mat < 0)
any(counts_mat != round(counts_mat))

# Identify the ten most highly expressed genes
top_genes <- sort(rowSums(counts_mat), decreasing = TRUE)[1:10]
top_genes



## Prepare sample metadata for DESeq2

# Identify sample columns in the featureCounts output
sample_cols <- grep("^06\\.Align/.+\\.sorted\\.bam$", names(counts_raw), value = TRUE)

# Import metadata
metadata_file <- "data/metadata.tsv"

meta <- read.delim(
  metadata_file,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)

colnames(meta) <- c("Well", "Group")

# Create sample names that match the featureCounts column names
meta$Sample <- paste0("06.Align/", meta$Well, ".sorted.bam")

# Check that all count matrix samples are present in metadata
missing_in_meta <- setdiff(sample_cols, meta$Sample)
missing_in_counts <- setdiff(meta$Sample, sample_cols)

if (length(missing_in_meta) > 0) {
  stop("These count matrix samples are missing in metadata:\n", paste(missing_in_meta, collapse = "\n"))
}

if (length(missing_in_counts) > 0) {
  warning("These metadata samples are not present in count matrix:\n", paste(missing_in_counts, collapse = "\n"))
}

# Reorder metadata to match count matrix sample order
meta_aligned <- meta[match(sample_cols, meta$Sample), ]

# Check aligned metadata
meta_aligned



## Prepare count matrix and metadata for DESeq2

# Keep only BAM/sample columns
sample_cols <- grep("^06\\.Align/.+\\.sorted\\.bam$", names(counts_raw), value = TRUE)

count_matrix <- counts_raw[, sample_cols]

# Use Geneid as row names
rownames(count_matrix) <- counts_raw$Geneid

# Convert to matrix
count_matrix <- as.matrix(count_matrix)
mode(count_matrix) <- "numeric"

# Use metadata aligned earlier
meta <- meta_aligned

# Set row names to sample names
rownames(meta) <- meta$Sample

# Reorder metadata to match count matrix columns
meta <- meta[colnames(count_matrix), ]

# Check that metadata and count matrix samples are in the same order
all(rownames(meta) == colnames(count_matrix))



## Prepare count matrix and sample metadata

# Extract sample columns and create count matrix
sample_cols <- grep("^06\\.Align/.+\\.sorted\\.bam$", names(counts_raw), value = TRUE)

count_matrix <- counts_raw[, sample_cols]
rownames(count_matrix) <- counts_raw$Geneid
count_matrix <- as.matrix(count_matrix)
mode(count_matrix) <- "numeric"

# Align metadata to the count matrix
meta <- meta_aligned
rownames(meta) <- meta$Sample
meta <- meta[colnames(count_matrix), , drop = FALSE]

stopifnot(all(rownames(meta) == colnames(count_matrix)))

# Parse condition and replicate information from Group
meta$Replicate <- sub("^.*_(\\d+)$", "\\1", meta$Group)
meta$Condition <- sub("_(\\d+)$", "", meta$Group)
meta$Replicate <- factor(meta$Replicate)
meta$Condition <- factor(meta$Condition)

# Create well ID for plotting labels
meta$Well <- sub("^06\\.Align/|\\.sorted\\.bam$", "", rownames(meta))



## DESeq2 object creation and variance stabilising transformation

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData   = meta,
  design    = ~ Condition
)

# Filter low-count genes
dds <- dds[rowSums(counts(dds)) > 10, ]

# Variance stabilising transformation
vsd <- vst(dds, blind = TRUE)



## PCA for all samples

# Generate PCA data
pca_data <- plotPCA(vsd, intgroup = "Condition", returnData = TRUE)

pca_data$Condition <- as.character(pca_data$Condition)

unique(pca_data$Condition)

# Rename conditions for plotting
pca_data$Condition <- recode(
  pca_data$Condition,
  "1:1000_EPR" = "1:1,000 EPR11334",
  "1:100000_EPR" = "1:100,000 EPR11334",
  "NgR3_siRNA" = "NgR3 siRNA",
  "Scrambled_siRNA" = "Scrambled siRNA",
  "IGG" = "IgG control"
)

pca_data$Condition <- factor(
  pca_data$Condition,
  levels = c(
    "IgG control",
    "1:1,000 EPR11334",
    "1:100,000 EPR11334",
    "NgR3 siRNA",
    "Scrambled siRNA"
  )
)

percentVar <- round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(PC1, PC2, color = Condition, label = name)) +
  geom_point(size = 4) +
  geom_text(vjust = -1) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  coord_fixed()

ggplot(pca_data, aes(PC1, PC2, color = Condition)) +
  geom_point(size = 6) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw(base_size = 18) +
  coord_fixed()

# Save PCA plot
p <- ggplot(pca_data, aes(PC1, PC2, color = Condition)) +
  geom_point(size = 6) +
  xlab(paste0("PC1 (", percentVar[1], "%)")) +
  ylab(paste0("PC2 (", percentVar[2], "%)")) +
  theme_bw(base_size = 18) +
  coord_fixed()

ggsave("/home/yougr345/PCA_plot.png", plot = p, width = 8, height = 6, dpi = 300)



## RNA-seq DEG counts barplot

# Count significant upregulated and downregulated genes

count_deg <- function(res, padj_cutoff = 0.05) {
  res_df <- as.data.frame(res)
  res_df <- res_df[!is.na(res_df$padj), ]
  
  up_count <- sum(res_df$padj < padj_cutoff & res_df$log2FoldChange > 0)
  down_count <- sum(res_df$padj < padj_cutoff & res_df$log2FoldChange < 0)
  
  data.frame(
    Upregulated = up_count,
    Downregulated = down_count
  )
}

# Build DEG count table

deg_counts <- bind_rows(
  count_deg(res_IGG_vs_1000),
  count_deg(res_IGG_vs_100000),
  count_deg(res_1000_vs_100000),
  count_deg(res_NgR3_vs_scrambled)
)

# Apply comparison labels

deg_counts$Comparison <- c(
  "IgG control vs. 1:1,000 EPR11334",
  "IgG control vs. 1:100,000 EPR11334",
  "1:1,000 EPR11334 vs. 1:100,000 EPR11334",
  "NgR3 siRNA vs. scrambled siRNA"
)

# Set comparison order

deg_counts$Comparison <- factor(
  deg_counts$Comparison,
  levels = deg_counts$Comparison
)

# Convert to long format

deg_long <- deg_counts %>%
  pivot_longer(
    cols = c(Upregulated, Downregulated),
    names_to = "Direction",
    values_to = "Count"
  )

# Plot DEG counts

p <- ggplot(deg_long, aes(x = Comparison, y = Count, fill = Direction)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  geom_text(
    aes(label = Count),
    position = position_dodge(width = 0.8),
    vjust = -0.4,
    size = 4
  ) +
  scale_fill_manual(
    values = c(
      "Upregulated" = "firebrick",
      "Downregulated" = "steelblue"
    ),
    labels = c(
      "Upregulated" = "log2FC > 0",
      "Downregulated" = "log2FC < 0"
    )
  ) +
  scale_y_log10() +
  labs(
    x = "Comparison",
    y = "Number of Genes (log scale)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 12),
    axis.title.x = element_text(margin = margin(t = 20)),
    axis.title.y = element_text(margin = margin(r = 15)),
    legend.title = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0.5, margin = margin(b = 15)),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 20, r = 30, b = 80, l = 30)
  )

p

# Save DEG counts barplot

ggsave(
  filename = "/home/yougr345/DEG_counts_barplot.png",
  plot = p,
  width = 14,
  height = 8,
  dpi = 300
)


## Volcano plots for DEG comparisons

## 1) IgG control vs 1:1,000 EPR11334

# Create comparison factor
meta$Condition <- sub("_(\\d+)$", "", meta$Group)
meta$comparison_group <- NA_character_

# IgG samples
meta$comparison_group[grepl("IGG", meta$Condition, ignore.case = TRUE)] <- "IGG"

# 1:1,000 EPR11334 samples
meta$comparison_group[
  grepl("EPR", meta$Condition, ignore.case = TRUE) &
    grepl("1000", meta$Condition) &
    !grepl("100000", meta$Condition)
] <- "EPR_1to1000"

meta$comparison_group <- factor(
  meta$comparison_group,
  levels = c("EPR_1to1000", "IGG")
)

# Keep only samples for this comparison
keep <- !is.na(meta$comparison_group)

meta_sub  <- meta[keep, , drop = FALSE]
count_sub <- count_matrix[, rownames(meta_sub), drop = FALSE]

# Check sample numbers
print(table(meta_sub$comparison_group))

# Run DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = count_sub,
  colData   = meta_sub,
  design    = ~ comparison_group
)

dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)

res_IGG_vs_1000 <- results(
  dds,
  contrast = c("comparison_group", "IGG", "EPR_1to1000")
)

# Save DESeq2 results object
saveRDS(
  res_IGG_vs_1000,
  file = "~/res_IGG_vs_1000.rds"
)

# Prepare volcano plot data
res_df <- as.data.frame(res_IGG_vs_1000)
res_df$Geneid <- rownames(res_df)

# Remove NA rows
res_df <- res_df[!is.na(res_df$padj) & !is.na(res_df$log2FoldChange), ]

res_df$neglog10padj <- -log10(res_df$padj)
res_df$neglog10padj[is.infinite(res_df$neglog10padj)] <- NA_real_

# Classify significance
res_df$Sig <- "Not sig"
res_df$Sig[
  res_df$padj < 0.05 &
    abs(res_df$log2FoldChange) >= 1
] <- "Significant"

# Volcano plot
ggplot(res_df, aes(x = log2FoldChange, y = neglog10padj)) +
  geom_point(aes(color = Sig), size = 2) +
  scale_color_manual(values = c(
    "Not sig" = "grey70",
    "Significant" = "red"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  labs(
    x = "log2 fold change (IgG control vs 1:1,000 EPR11334)",
    y = "-log10 adjusted p-value"
  ) +
  coord_cartesian(xlim = c(-3, 3)) +
  theme_bw(base_size = 16) +
  theme(legend.title = element_blank())

png(
  "volcano_IGG_vs_1000.png",
  width = 2400,
  height = 1800,
  res = 300
)

grid.newpage()

pushViewport(viewport(
  width = 0.88,
  height = 0.88,
  x = 0.50,
  y = 0.53
))

print(last_plot(), newpage = FALSE)

dev.off()

# Export significant upregulated and downregulated genes
sig_genes_IGGv1_1000EPR <- subset(res_df, Sig == "Significant")

sig_up_genes_IGGv1_1000EPR <- subset(
  sig_genes_IGGv1_1000EPR,
  log2FoldChange > 0
)

sig_down_genes_IGGv1_1000EPR <- subset(
  sig_genes_IGGv1_1000EPR,
  log2FoldChange < 0
)

write.csv(
  sig_up_genes_IGGv1_1000EPR,
  "sig_up_genes_IGGv1_1000EPR.csv",
  row.names = FALSE
)

write.csv(
  sig_down_genes_IGGv1_1000EPR,
  "sig_down_genes_IGGv1_1000EPR.csv",
  row.names = FALSE
)


## 2) IgG control vs 1:100,000 EPR11334

# Create comparison factor
meta$Condition <- sub("_(\\d+)$", "", meta$Group)
meta$comparison_group <- NA_character_

# IgG samples
meta$comparison_group[
  grepl("IGG", meta$Condition, ignore.case = TRUE)
] <- "IGG"

# 1:100,000 EPR11334 samples
meta$comparison_group[
  grepl("EPR", meta$Condition, ignore.case = TRUE) &
    grepl("100000", meta$Condition)
] <- "EPR_1to100000"

meta$comparison_group <- factor(
  meta$comparison_group,
  levels = c("EPR_1to100000", "IGG")
)

# Keep only samples for this comparison
keep <- !is.na(meta$comparison_group)

meta_sub  <- meta[keep, , drop = FALSE]
count_sub <- count_matrix[, rownames(meta_sub), drop = FALSE]

print(table(meta_sub$comparison_group))

# Run DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = count_sub,
  colData   = meta_sub,
  design    = ~ comparison_group
)

dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)

res_IGG_vs_100000 <- results(
  dds,
  contrast = c("comparison_group", "IGG", "EPR_1to100000")
)

# Save DESeq2 results object
saveRDS(
  res_IGG_vs_100000,
  file = "~/res_IGG_vs_100000.rds"
)

# Save CSV version of DESeq2 results
res_IGG_vs_100000_df <- as.data.frame(res_IGG_vs_100000)
res_IGG_vs_100000_df$Geneid <- rownames(res_IGG_vs_100000_df)

write.csv(
  res_IGG_vs_100000_df,
  file = "res_IGG_vs_100000.csv",
  row.names = FALSE
)

# Prepare volcano plot data
res_df <- as.data.frame(res_IGG_vs_100000)
res_df$Geneid <- rownames(res_df)

res_df$neglog10padj <- -log10(res_df$padj)
res_df$neglog10padj[is.infinite(res_df$neglog10padj)] <- NA_real_

# Classify significance
res_df$Sig <- "Not sig"
res_df$Sig[
  !is.na(res_df$padj) &
    res_df$padj < 0.05 &
    abs(res_df$log2FoldChange) >= 1
] <- "Significant"

# Volcano plot
ggplot(res_df, aes(x = log2FoldChange, y = neglog10padj)) +
  geom_point(aes(color = Sig), size = 2) +
  scale_color_manual(values = c(
    "Not sig" = "grey70",
    "Significant" = "red"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  labs(
    x = "log2 fold change (IgG control vs 1:100,000 EPR11334)",
    y = "-log10 adjusted p-value"
  ) +
  theme_bw(base_size = 16) +
  theme(legend.title = element_blank())

png(
  "volcano_IGG_vs_100000.png",
  width = 2400,
  height = 1800,
  res = 300
)

grid.newpage()

pushViewport(viewport(
  width = 0.88,
  height = 0.88,
  x = 0.50,
  y = 0.53
))

print(last_plot(), newpage = FALSE)

dev.off()

# Export significant upregulated and downregulated genes
sig_genes_IGGv1_100000EPR <- subset(res_df, Sig == "Significant")

sig_up_genes_IGGv1_100000EPR <- subset(
  sig_genes_IGGv1_100000EPR,
  log2FoldChange > 0
)

sig_down_genes_IGGv1_100000EPR <- subset(
  sig_genes_IGGv1_100000EPR,
  log2FoldChange < 0
)

write.csv(
  sig_up_genes_IGGv1_100000EPR,
  "sig_up_genes_IGGv1_100000EPR.csv",
  row.names = FALSE
)

write.csv(
  sig_down_genes_IGGv1_100000EPR,
  "sig_down_genes_IGGv1_100000EPR.csv",
  row.names = FALSE
)


## 3) EPR11334 dose comparison: 1:1,000 vs 1:100,000

# Create comparison factor
meta$Condition <- sub("_(\\d+)$", "", meta$Group)
meta$comparison_group <- NA_character_

# 1:1,000 EPR11334 samples
meta$comparison_group[
  grepl("EPR", meta$Condition, ignore.case = TRUE) &
    grepl("1000", meta$Condition) &
    !grepl("100000", meta$Condition)
] <- "EPR_1to1000"

# 1:100,000 EPR11334 samples
meta$comparison_group[
  grepl("EPR", meta$Condition, ignore.case = TRUE) &
    grepl("100000", meta$Condition)
] <- "EPR_1to100000"

meta$comparison_group <- factor(
  meta$comparison_group,
  levels = c("EPR_1to100000", "EPR_1to1000")
)

# Keep only samples for this comparison
keep <- !is.na(meta$comparison_group)

meta_sub  <- meta[keep, , drop = FALSE]
count_sub <- count_matrix[, rownames(meta_sub), drop = FALSE]

print(table(meta_sub$comparison_group))

# Run DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = count_sub,
  colData   = meta_sub,
  design    = ~ comparison_group
)

dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)

res_1000_vs_100000 <- results(
  dds,
  contrast = c("comparison_group", "EPR_1to1000", "EPR_1to100000")
)

# Save DESeq2 results object
saveRDS(
  res_1000_vs_100000,
  file = "~/res_1000_vs_100000.rds"
)

list.files("~", pattern = "res_1000_vs_100000")

# Save CSV version of DESeq2 results
res_1000_vs_100000_df <- as.data.frame(res_1000_vs_100000)
res_1000_vs_100000_df$Geneid <- rownames(res_1000_vs_100000_df)

write.csv(
  res_1000_vs_100000_df,
  file = "res_1000_vs_100000.csv",
  row.names = FALSE
)

# Prepare volcano plot data
res_df <- as.data.frame(res_1000_vs_100000)
res_df$Geneid <- rownames(res_df)

res_df$neglog10padj <- -log10(res_df$padj)
res_df$neglog10padj[is.infinite(res_df$neglog10padj)] <- NA_real_

# Classify significance
res_df$Sig <- "Not sig"
res_df$Sig[
  !is.na(res_df$padj) &
    res_df$padj < 0.05 &
    abs(res_df$log2FoldChange) >= 1
] <- "Significant"

# Volcano plot
ggplot(res_df, aes(x = log2FoldChange, y = neglog10padj)) +
  geom_point(aes(color = Sig), size = 2) +
  scale_color_manual(values = c(
    "Not sig" = "grey70",
    "Significant" = "red"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  labs(
    x = "log2 fold change (1:1,000 EPR11334 vs 1:100,000 EPR11334)",
    y = "-log10 adjusted p-value"
  ) +
  theme_bw(base_size = 16) +
  theme(
    legend.title = element_blank(),
    plot.margin = margin(t = 20, r = 40, b = 20, l = 20)
  )

# Export significant upregulated and downregulated genes
sig_genes_1_1000EPRv1_100000EPR <- subset(res_df, Sig == "Significant")

sig_up_genes_1_1000EPRv1_100000EPR <- subset(
  sig_genes_1_1000EPRv1_100000EPR,
  log2FoldChange > 0
)

sig_down_genes_1_1000EPRv1_100000EPR <- subset(
  sig_genes_1_1000EPRv1_100000EPR,
  log2FoldChange < 0
)

write.csv(
  sig_up_genes_1_1000EPRv1_100000EPR,
  "sig_up_genes_1_1000EPRv1_100000EPR.csv",
  row.names = FALSE
)

write.csv(
  sig_down_genes_1_1000EPRv1_100000EPR,
  "sig_down_genes_1_1000EPRv1_100000EPR.csv",
  row.names = FALSE
)

png(
  "volcano_1000_vs_100000.png",
  width = 2400,
  height = 1800,
  res = 300
)

grid.newpage()

pushViewport(viewport(
  width = 0.88,
  height = 0.88,
  x = 0.50,
  y = 0.53
))

print(last_plot(), newpage = FALSE)

dev.off()


## 4) NgR3 siRNA vs scrambled siRNA

library(DESeq2)
library(ggplot2)

# Create comparison factor
meta$Condition <- sub("_(\\d+)$", "", meta$Group)
meta$comparison_group <- NA_character_

# Scrambled siRNA samples
meta$comparison_group[
  grepl("scram|scrambled", meta$Condition, ignore.case = TRUE)
] <- "Scrambled_siRNA"

# NgR3 siRNA samples
meta$comparison_group[
  grepl("NgR3", meta$Condition, ignore.case = TRUE)
] <- "NgR3_siRNA"

meta$comparison_group <- factor(
  meta$comparison_group,
  levels = c("Scrambled_siRNA", "NgR3_siRNA")
)

# Keep only samples for this comparison
keep <- !is.na(meta$comparison_group)

meta_sub  <- meta[keep, , drop = FALSE]
count_sub <- count_matrix[, rownames(meta_sub), drop = FALSE]

# Check sample numbers
print(table(meta_sub$comparison_group))

# Run DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = count_sub,
  colData   = meta_sub,
  design    = ~ comparison_group
)

dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)

# NgR3 vs scrambled siRNA
res_NgR3_vs_scrambled <- results(
  dds,
  contrast = c("comparison_group", "NgR3_siRNA", "Scrambled_siRNA")
)

# Save DESeq2 results object
saveRDS(
  res_NgR3_vs_scrambled,
  file = "~/res_NgR3_vs_scrambled.rds"
)

# Prepare volcano plot data
res_df <- as.data.frame(res_NgR3_vs_scrambled)
res_df$Geneid <- rownames(res_df)

# Remove NA rows
res_df <- res_df[!is.na(res_df$padj) & !is.na(res_df$log2FoldChange), ]

res_df$neglog10padj <- -log10(res_df$padj)
res_df$neglog10padj[is.infinite(res_df$neglog10padj)] <- NA_real_

# Classify significance
res_df$Sig <- "Not sig"
res_df$Sig[
  res_df$padj < 0.05 &
    abs(res_df$log2FoldChange) >= 1
] <- "Significant"

# Volcano plot
ggplot(res_df, aes(x = log2FoldChange, y = neglog10padj)) +
  geom_point(aes(color = Sig), size = 2) +
  scale_color_manual(values = c(
    "Not sig" = "grey70",
    "Significant" = "red"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  labs(
    x = "log2 fold change (NgR3 siRNA vs Scrambled siRNA)",
    y = "-log10 adjusted p-value"
  ) +
  theme_bw(base_size = 16) +
  theme(legend.title = element_blank())

png("volcano_NgR3_vs_scrambled.png", width = 2400, height = 1800, res = 300)

grid.newpage()

pushViewport(viewport(
  width = 0.88,
  height = 0.88,
  x = 0.50,
  y = 0.53
))

print(last_plot(), newpage = FALSE)

dev.off()

# Export significant upregulated and downregulated genes
sig_genes_NgR3_siRNAvScrambled_siRNA <- subset(res_df, Sig == "Significant")

sig_up_genes_NgR3_siRNAvScrambled_siRNA <- subset(
  sig_genes_NgR3_siRNAvScrambled_siRNA,
  log2FoldChange > 0
)

sig_down_genes_NgR3_siRNAvScrambled_siRNA <- subset(
  sig_genes_NgR3_siRNAvScrambled_siRNA,
  log2FoldChange < 0
)

write.csv(
  sig_up_genes_NgR3_siRNAvScrambled_siRNA,
  "sig_up_genes_NgR3_siRNAvScrambled_siRNA.csv",
  row.names = FALSE
)

write.csv(
  sig_down_genes_NgR3_siRNAvScrambled_siRNA,
  "sig_down_genes_NgR3_siRNAvScrambled_siRNA.csv",
  row.names = FALSE
)



## Heatmap 1: Top 100 DESeq2 genes from IgG vs 1:100,000 EPR across all samples

## 1) Use existing VSD object

stopifnot(exists("vsd"))

mat <- assay(vsd)
meta <- as.data.frame(colData(vsd))

# Rename groups
meta$group <- as.character(meta$Condition)
meta$group[meta$group == "IGG"] <- "IgG control"
meta$group[meta$group == "1:1000_EPR"] <- "1:1,000 EPR11334"
meta$group[meta$group == "1:100000_EPR"] <- "1:100,000 EPR11334"
meta$group[meta$group == "NgR3_siRNA"] <- "NgR3 siRNA"
meta$group[meta$group == "Scrambled_siRNA"] <- "Scrambled siRNA"

rownames(meta) <- colnames(mat)
meta <- meta[colnames(mat), , drop = FALSE]

## 2) Get top 100 genes from DESeq2 results

res_df <- as.data.frame(res_IGG_vs_100000)
res_df$EnsemblID <- rownames(res_df)

res_df <- res_df %>%
  filter(!is.na(padj)) %>%
  arrange(padj)

top100 <- head(res_df, 100)

## 3) Map Ensembl IDs to gene symbols

ensembl_ids <- rownames(mat)
ensembl_ids_clean <- sub("\\..*$", "", ensembl_ids)

gene_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(ensembl_ids_clean),
  keytype = "ENSEMBL",
  columns = c("SYMBOL")
)

gene_map <- gene_map[!is.na(gene_map$SYMBOL) & gene_map$SYMBOL != "", ]
gene_map <- gene_map[!duplicated(gene_map$ENSEMBL), ]

matched_symbols <- gene_map$SYMBOL[match(ensembl_ids_clean, gene_map$ENSEMBL)]

keep <- !is.na(matched_symbols) & matched_symbols != ""
mat <- mat[keep, , drop = FALSE]
matched_symbols <- matched_symbols[keep]

rownames(mat) <- matched_symbols

# Remove duplicate gene symbols
row_var <- apply(mat, 1, var)
mat <- mat[order(row_var, decreasing = TRUE), ]
mat <- mat[!duplicated(rownames(mat)), ]

# Match top 100 genes to matrix
top_genes <- intersect(top100$EnsemblID, ensembl_ids)
top_genes_clean <- sub("\\..*$", "", top_genes)

top_symbols <- gene_map$SYMBOL[match(top_genes_clean, gene_map$ENSEMBL)]
top_symbols <- top_symbols[!is.na(top_symbols)]

top_symbols <- intersect(top_symbols, rownames(mat))

## 4) Build heatmap matrix

heatmap_mat <- mat[top_symbols, , drop = FALSE]

# Z-score
heatmap_mat_z <- t(scale(t(heatmap_mat)))
heatmap_mat_z <- heatmap_mat_z[complete.cases(heatmap_mat_z), ]

## 5) Create row annotation

top100$EnsemblID_clean <- sub("\\..*$", "", top100$EnsemblID)

top100$Symbol <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = top100$EnsemblID_clean,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

top100$Direction <- ifelse(
  top100$log2FoldChange > 0,
  "Upregulated in IgG control",
  "Upregulated in 1:100,000 EPR11334"
)

# Keep rows with valid symbols
top100_annot <- top100[, c("Symbol", "Direction", "padj")]
top100_annot <- top100_annot[!is.na(top100_annot$Symbol) & top100_annot$Symbol != "", ]

# Remove duplicate symbols, keeping the most significant one
top100_annot <- top100_annot[order(top100_annot$padj), ]
top100_annot <- top100_annot[!duplicated(top100_annot$Symbol), ]

# Keep only genes present in the heatmap matrix
top100_annot <- top100_annot[top100_annot$Symbol %in% rownames(heatmap_mat_z), ]

row_annotation <- data.frame(
  Direction = top100_annot$Direction,
  row.names = top100_annot$Symbol,
  stringsAsFactors = FALSE
)

# Subset heatmap matrix to match row annotation
heatmap_mat_z <- heatmap_mat_z[rownames(row_annotation), , drop = FALSE]

## 6) Order genes

genes_igg <- rownames(row_annotation)[row_annotation$Direction == "Higher in IGG"]
genes_epr <- rownames(row_annotation)[row_annotation$Direction == "Higher in 1:100,000 EPR"]

ordered_genes <- c()

for (dir in c("Upregulated in IgG control", "Upregulated in 1:100,000 EPR11334")) {
  
  genes <- rownames(row_annotation)[row_annotation$Direction == dir]
  
  if (length(genes) > 1) {
    block <- heatmap_mat_z[genes, , drop = FALSE]
    hc <- hclust(dist(block))
    ordered_genes <- c(ordered_genes, rownames(block)[hc$order])
  } else {
    ordered_genes <- c(ordered_genes, genes)
  }
}

heatmap_mat_z <- heatmap_mat_z[ordered_genes, , drop = FALSE]
row_annotation <- row_annotation[ordered_genes, , drop = FALSE]

gap_rows <- sum(row_annotation$Direction == "Upregulated in IgG control")

## 7) Order columns

desired_group_order <- c(
  "IgG control",
  "1:1,000 EPR11334",
  "1:100,000 EPR11334",
  "Scrambled siRNA",
  "NgR3 siRNA"
)

col_annotation <- data.frame(
  Group = meta$group,
  row.names = rownames(meta)
)

col_annotation$Group <- factor(col_annotation$Group, levels = desired_group_order)

sample_order <- order(col_annotation$Group)

heatmap_mat_z <- heatmap_mat_z[, sample_order]
col_annotation <- col_annotation[colnames(heatmap_mat_z), , drop = FALSE]

## 8) Set colours

annotation_colors <- list(
  Group = c(
    "IgG control" = "#BDBDBD",
    "1:1,000 EPR11334" = "#66C2A5",
    "1:100,000 EPR11334" = "#1B9E77",
    "Scrambled siRNA" = "#E6AB02",
    "NgR3 siRNA" = "#7570B3"
  ),
  Direction = c(
    "Upregulated in IgG control" = "#4DAF4A",
    "Upregulated in 1:100,000 EPR11334" = "#984EA3"
  )
)

breaks_use <- seq(-3, 3, length.out = 101)

## 9) Plot and save heatmap

png(
  "Top100_DESEQ2_genes_IGG_vs_100000_across_all_samples.png",
  width = 4200,
  height = 5200,
  res = 300
)

grid.newpage()

p <- pheatmap(
  heatmap_mat_z,
  annotation_row = row_annotation,
  annotation_col = col_annotation,
  annotation_colors = annotation_colors,
  show_rownames = TRUE,
  show_colnames = FALSE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  gaps_row = gap_rows,
  fontsize = 10,
  fontsize_row = 6.5,
  fontsize_col = 9,
  border_color = NA,
  breaks = breaks_use,
  color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
  silent = TRUE
)

pushViewport(
  viewport(
    x = 0.5,
    y = 0.5,
    width = 0.90,
    height = 0.92
  )
)

grid.draw(p$gtable)
popViewport()

dev.off()


## Heatmap 2: Top 100 DESeq2 genes from NgR3 siRNA vs scrambled siRNA across all samples

## 1) Use existing VSD object

stopifnot(exists("vsd"))

mat <- assay(vsd)
meta <- as.data.frame(colData(vsd))

# Rename groups
meta$group <- as.character(meta$Condition)
meta$group[meta$group == "IGG"] <- "IgG control"
meta$group[meta$group == "1:1000_EPR"] <- "1:1,000 EPR11334"
meta$group[meta$group == "1:100000_EPR"] <- "1:100,000 EPR11334"
meta$group[meta$group == "NgR3_siRNA"] <- "NgR3 siRNA"
meta$group[meta$group == "Scrambled_siRNA"] <- "Scrambled siRNA"

rownames(meta) <- colnames(mat)
meta <- meta[colnames(mat), , drop = FALSE]

## 2) Get top 100 genes from DESeq2 results

res_df <- as.data.frame(res_NgR3_vs_scrambled)
res_df$EnsemblID <- rownames(res_df)

res_df <- res_df %>%
  filter(!is.na(padj)) %>%
  arrange(padj)

top100 <- head(res_df, 100)

## 3) Map Ensembl IDs to gene symbols

ensembl_ids <- rownames(mat)
ensembl_ids_clean <- sub("\\..*$", "", ensembl_ids)

gene_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(ensembl_ids_clean),
  keytype = "ENSEMBL",
  columns = c("SYMBOL")
)

gene_map <- gene_map[!is.na(gene_map$SYMBOL) & gene_map$SYMBOL != "", ]
gene_map <- gene_map[!duplicated(gene_map$ENSEMBL), ]

matched_symbols <- gene_map$SYMBOL[match(ensembl_ids_clean, gene_map$ENSEMBL)]

keep <- !is.na(matched_symbols) & matched_symbols != ""
mat <- mat[keep, , drop = FALSE]
matched_symbols <- matched_symbols[keep]

rownames(mat) <- matched_symbols

# Remove duplicate gene symbols
row_var <- apply(mat, 1, var)
mat <- mat[order(row_var, decreasing = TRUE), ]
mat <- mat[!duplicated(rownames(mat)), ]

# Match top 100 genes to matrix
top_genes <- intersect(top100$EnsemblID, ensembl_ids)
top_genes_clean <- sub("\\..*$", "", top_genes)

top_symbols <- gene_map$SYMBOL[match(top_genes_clean, gene_map$ENSEMBL)]
top_symbols <- top_symbols[!is.na(top_symbols)]

top_symbols <- intersect(top_symbols, rownames(mat))

## 4) Build heatmap matrix

heatmap_mat <- mat[top_symbols, , drop = FALSE]

# Gene-wise z-score
heatmap_mat_z <- t(scale(t(heatmap_mat)))
heatmap_mat_z <- heatmap_mat_z[complete.cases(heatmap_mat_z), ]

## 5) Create row annotation

top100$EnsemblID_clean <- sub("\\..*$", "", top100$EnsemblID)

top100$Symbol <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = top100$EnsemblID_clean,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

top100$Direction <- ifelse(
  top100$log2FoldChange > 0,
  "Upregulated in NgR3 siRNA",
  "Upregulated in Scrambled siRNA"
)

# Keep rows with valid symbols
top100_annot <- top100[, c("Symbol", "Direction", "padj")]
top100_annot <- top100_annot[!is.na(top100_annot$Symbol) & top100_annot$Symbol != "", ]

# Remove duplicate symbols, keeping the most significant one
top100_annot <- top100_annot[order(top100_annot$padj), ]
top100_annot <- top100_annot[!duplicated(top100_annot$Symbol), ]

# Keep only genes present in the heatmap matrix
top100_annot <- top100_annot[top100_annot$Symbol %in% rownames(heatmap_mat_z), ]

row_annotation <- data.frame(
  Direction = top100_annot$Direction,
  row.names = top100_annot$Symbol,
  stringsAsFactors = FALSE
)

# Subset heatmap matrix to match row annotation
heatmap_mat_z <- heatmap_mat_z[rownames(row_annotation), , drop = FALSE]

## 6) Order genes

ordered_genes <- c()

for (dir in c("Upregulated in NgR3 siRNA", "Upregulated in Scrambled siRNA")) {
  
  genes <- rownames(row_annotation)[row_annotation$Direction == dir]
  
  if (length(genes) > 1) {
    block <- heatmap_mat_z[genes, , drop = FALSE]
    hc <- hclust(dist(block))
    ordered_genes <- c(ordered_genes, rownames(block)[hc$order])
  } else {
    ordered_genes <- c(ordered_genes, genes)
  }
}

heatmap_mat_z <- heatmap_mat_z[ordered_genes, , drop = FALSE]
row_annotation <- row_annotation[ordered_genes, , drop = FALSE]

gap_rows <- sum(row_annotation$Direction == "Upregulated in NgR3 siRNA")

## 7) Order columns

desired_group_order <- c(
  "IgG control",
  "1:1,000 EPR11334",
  "1:100,000 EPR11334",
  "Scrambled siRNA",
  "NgR3 siRNA"
)

col_annotation <- data.frame(
  Group = meta$group,
  row.names = rownames(meta)
)

col_annotation$Group <- factor(col_annotation$Group, levels = desired_group_order)

sample_order <- order(col_annotation$Group)

heatmap_mat_z <- heatmap_mat_z[, sample_order]
col_annotation <- col_annotation[colnames(heatmap_mat_z), , drop = FALSE]

## 8) Set colours

annotation_colors <- list(
  Group = c(
    "IgG control" = "#BDBDBD",
    "1:1,000 EPR11334" = "#66C2A5",
    "1:100,000 EPR11334" = "#1B9E77",
    "Scrambled siRNA" = "#E6AB02",
    "NgR3 siRNA" = "#7570B3"
  ),
  Direction = c(
    "Upregulated in NgR3 siRNA" = "#4DAF4A",
    "Upregulated in Scrambled siRNA" = "#984EA3"
  )
)

breaks_use <- seq(-3, 3, length.out = 101)

## 9) Plot and save heatmap

png(
  "Top100_DESEQ2_genes_NgR3_vs_Scrambled_across_all_samples.png",
  width = 4200,
  height = 5200,
  res = 300
)

grid.newpage()

p <- pheatmap(
  heatmap_mat_z,
  annotation_row = row_annotation,
  annotation_col = col_annotation,
  annotation_colors = annotation_colors,
  show_rownames = TRUE,
  show_colnames = FALSE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  gaps_row = gap_rows,
  fontsize = 10,
  fontsize_row = 6.5,
  fontsize_col = 9,
  border_color = NA,
  breaks = breaks_use,
  color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
  silent = TRUE
)

pushViewport(
  viewport(
    x = 0.5,
    y = 0.5,
    width = 0.90,
    height = 0.92
  )
)

grid.draw(p$gtable)
popViewport()

dev.off()



## GSEA hallmark analysis: NgR3 siRNA vs scrambled siRNA

## 1) Create comparison factor

meta$comparison_group <- NA_character_

# Scrambled siRNA samples
meta$comparison_group[
  grepl("Scrambled", meta$Group, ignore.case = TRUE)
] <- "Scrambled_siRNA"

# NgR3 siRNA samples
meta$comparison_group[
  grepl("NgR3", meta$Group, ignore.case = TRUE)
] <- "NgR3_siRNA"

meta$comparison_group <- factor(
  meta$comparison_group,
  levels = c("Scrambled_siRNA", "NgR3_siRNA")
)

## 2) Keep only samples for this comparison

keep <- !is.na(meta$comparison_group)

meta_sub  <- meta[keep, , drop = FALSE]
count_sub <- count_matrix[, rownames(meta_sub), drop = FALSE]

print(table(meta_sub$comparison_group))

## 3) Build DESeq2 dataset

dds_NgR3_vs_Scrambled <- DESeqDataSetFromMatrix(
  countData = count_sub,
  colData   = meta_sub,
  design    = ~ comparison_group
)

dds_NgR3_vs_Scrambled <- dds_NgR3_vs_Scrambled[
  rowSums(counts(dds_NgR3_vs_Scrambled)) > 10,
]

dds_NgR3_vs_Scrambled <- DESeq(dds_NgR3_vs_Scrambled)

resultsNames(dds_NgR3_vs_Scrambled)

## 4) Get DESeq2 results

res_NgR3_vs_Scrambled <- results(
  dds_NgR3_vs_Scrambled,
  contrast = c("comparison_group", "NgR3_siRNA", "Scrambled_siRNA")
)

## 5) Save DESeq2 results

saveRDS(
  res_NgR3_vs_Scrambled,
  file = "~/res_NgR3_vs_Scrambled.rds"
)

res_NgR3_vs_Scrambled_export <- as.data.frame(res_NgR3_vs_Scrambled)
res_NgR3_vs_Scrambled_export$Geneid <- rownames(res_NgR3_vs_Scrambled_export)

write.csv(
  res_NgR3_vs_Scrambled_export,
  file = "res_NgR3_vs_Scrambled.csv",
  row.names = FALSE
)

if (file.exists("~/res_NgR3_vs_Scrambled.rds")) {
  message("DESeq2 result saved successfully")
} else {
  stop("DESeq2 result save failed")
}

## 6) Create ranked gene list for fgsea

gene_list_NgR3_vs_Scrambled <- res_NgR3_vs_Scrambled$log2FoldChange
names(gene_list_NgR3_vs_Scrambled) <- rownames(res_NgR3_vs_Scrambled)

# Remove NA log2 fold changes
gene_list_NgR3_vs_Scrambled <- gene_list_NgR3_vs_Scrambled[
  !is.na(gene_list_NgR3_vs_Scrambled)
]

gene_symbols_NgR3_vs_Scrambled <- mapIds(
  org.Hs.eg.db,
  keys = names(gene_list_NgR3_vs_Scrambled),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

gene_df_NgR3_vs_Scrambled <- data.frame(
  symbol = gene_symbols_NgR3_vs_Scrambled,
  log2FoldChange = gene_list_NgR3_vs_Scrambled,
  stringsAsFactors = FALSE
)

gene_df_NgR3_vs_Scrambled <- gene_df_NgR3_vs_Scrambled[
  !is.na(gene_df_NgR3_vs_Scrambled$symbol),
]

gene_df_NgR3_vs_Scrambled <- gene_df_NgR3_vs_Scrambled %>%
  group_by(symbol) %>%
  slice_max(order_by = abs(log2FoldChange), n = 1, with_ties = FALSE) %>%
  ungroup()

gene_list_NgR3_vs_Scrambled <- gene_df_NgR3_vs_Scrambled$log2FoldChange
names(gene_list_NgR3_vs_Scrambled) <- gene_df_NgR3_vs_Scrambled$symbol

gene_list_NgR3_vs_Scrambled <- sort(gene_list_NgR3_vs_Scrambled, decreasing = TRUE)

## 7) Load hallmark gene sets

msig <- msigdbr(species = "Homo sapiens", category = "H")
pathways <- split(msig$gene_symbol, msig$gs_name)

## 8) Run fgsea

fgseaRes_NgR3_vs_Scrambled <- fgsea(
  pathways = pathways,
  stats = gene_list_NgR3_vs_Scrambled,
  minSize = 15,
  maxSize = 500,
  nperm = 10000
)

fgseaRes_NgR3_vs_Scrambled <- fgseaRes_NgR3_vs_Scrambled %>%
  arrange(padj)

head(fgseaRes_NgR3_vs_Scrambled, 20)

## 9) Save fgsea results

saveRDS(
  fgseaRes_NgR3_vs_Scrambled,
  file = "~/fgseaRes_NgR3_vs_Scrambled.rds"
)

fgsea_export_NgR3_vs_Scrambled <- fgseaRes_NgR3_vs_Scrambled %>%
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ", "))

write.csv(
  fgsea_export_NgR3_vs_Scrambled,
  file = "fgseaRes_NgR3_vs_Scrambled.csv",
  row.names = FALSE
)

write_xlsx(
  fgsea_export_NgR3_vs_Scrambled,
  "gsea_results_NgR3_vs_Scrambled.xlsx"
)


## Hallmark GSEA dot plot: NgR3 siRNA vs scrambled siRNA

## 1) Prepare significant Hallmark pathways

plot_df <- fgseaRes_NgR3_vs_Scrambled %>%
  as.data.frame() %>%
  filter(!is.na(padj)) %>%
  filter(padj < 0.05)

## 2) Clean pathway labels

plot_df$pathway_clean <- NA_character_

plot_df$pathway_clean[plot_df$pathway == "HALLMARK_IL2_STAT5_SIGNALING"] <- "IL2/STAT5 signalling"
plot_df$pathway_clean[plot_df$pathway == "HALLMARK_E2F_TARGETS"] <- "E2F targets"
plot_df$pathway_clean[plot_df$pathway == "HALLMARK_IL6_JAK_STAT3_SIGNALING"] <- "IL6/JAK/STAT3 signalling"
plot_df$pathway_clean[plot_df$pathway == "HALLMARK_INFLAMMATORY_RESPONSE"] <- "Inflammatory response"
plot_df$pathway_clean[plot_df$pathway == "HALLMARK_ESTROGEN_RESPONSE_EARLY"] <- "Estrogen response early"
plot_df$pathway_clean[plot_df$pathway == "HALLMARK_TNFA_SIGNALING_VIA_NFKB"] <- "TNFα/NFκB signalling"
plot_df$pathway_clean[plot_df$pathway == "HALLMARK_UV_RESPONSE_UP"] <- "UV response up"

## 3) Add enrichment direction

plot_df$direction <- ifelse(
  plot_df$NES > 0,
  "Gained with NgR3 knockdown",
  "Lost with NgR3 knockdown"
)

## 4) Set pathway order

plot_df$pathway_clean <- base::factor(
  plot_df$pathway_clean,
  levels = c(
    "IL2/STAT5 signalling",
    "IL6/JAK/STAT3 signalling",
    "Inflammatory response",
    "Estrogen response early",
    "TNFα/NFκB signalling",
    "UV response up",
    "E2F targets"
  )
)

## 5) Split pathways by direction

plot_df_gained <- plot_df %>%
  filter(direction == "Gained with NgR3 knockdown")

plot_df_lost <- plot_df %>%
  filter(direction == "Lost with NgR3 knockdown")

## 6) Set shared plot scales

x_limits <- c(-2, 2)
x_breaks <- c(-2, -1, 0, 1, 2)

colour_limits <- range(-log10(plot_df$padj), na.rm = TRUE)
colour_breaks <- pretty(colour_limits, n = 4)

colour_scale <- scale_colour_gradientn(
  colours = c("deepskyblue", "yellow", "red"),
  limits = colour_limits,
  breaks = colour_breaks,
  guide = guide_colourbar(order = 1)
)

size_scale <- scale_size_continuous(
  range = c(4, 12),
  limits = c(75, 199),
  breaks = c(75, 125, 199),
  labels = c("75", "125", "199"),
  guide = guide_legend(order = 2)
)

x_scale <- scale_x_continuous(
  limits = x_limits,
  breaks = x_breaks
)

## 7) Plot pathways gained with NgR3 knockdown

p_top <- ggplot(plot_df_gained, aes(x = NES, y = pathway_clean)) +
  geom_point(aes(size = size, colour = -log10(padj))) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  colour_scale +
  size_scale +
  x_scale +
  labs(
    title = "Gained with NgR3 knockdown",
    x = NULL,
    y = NULL,
    colour = "-log10(padj)",
    size = "Gene set size"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    legend.position = "none"
  )

## 8) Plot pathways lost with NgR3 knockdown

p_bottom <- ggplot(plot_df_lost, aes(x = NES, y = pathway_clean)) +
  geom_point(aes(size = size, colour = -log10(padj))) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  colour_scale +
  size_scale +
  x_scale +
  labs(
    title = "Lost with NgR3 knockdown",
    x = "Normalized Enrichment Score (NES)",
    y = "Pathway",
    colour = "-log10(padj)",
    size = "Gene set size"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 11),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    legend.position = "right"
  )

## 9) Combine dot plots

hallmark_dotplot <- p_top / p_bottom +
  plot_layout(heights = c(1, 3)) +
  plot_annotation(
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
    )
  )

hallmark_dotplot

## 10) Save Hallmark dot plot

ggsave(
  filename = "GSEA_Hallmark_dotplot_NgR3_siRNA_vs_Scrambled_siRNA.png",
  plot = hallmark_dotplot,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)


## Hallmark GSEA heatmap: NgR3 siRNA vs scrambled siRNA upregulated
## Leading-edge genes from E2F targets

## 1) Create VST expression matrix from full DESeq2 object

vsd <- vst(dds_full, blind = FALSE)
expr_mat <- assay(vsd)

## 2) Keep only NgR3 siRNA and scrambled siRNA samples

samples_keep <- colData(dds_full)$Condition %in% c("NgR3_siRNA", "Scrambled_siRNA")

expr_sub <- expr_mat[, samples_keep, drop = FALSE]
meta_sub <- as.data.frame(colData(dds_full)[samples_keep, ])

table(meta_sub$Condition)

## 3) Select pathway of interest

pathways_of_interest <- c(
  "HALLMARK_E2F_TARGETS"
)

fgsea_top <- fgseaRes_NgR3_vs_Scrambled %>%
  as.data.frame() %>%
  filter(pathway %in% pathways_of_interest)

## 4) Extract leading-edge genes

leading_edge_list <- setNames(fgsea_top$leadingEdge, fgsea_top$pathway)

top_n <- 15
leading_edge_top_list <- lapply(leading_edge_list, function(x) head(x, top_n))
heatmap_genes <- unique(unlist(leading_edge_top_list))

## 5) Map Ensembl IDs to gene symbols

ensembl_ids <- rownames(expr_sub)
ensembl_ids_clean <- sub("\\..*$", "", ensembl_ids)

gene_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(ensembl_ids_clean),
  keytype = "ENSEMBL",
  columns = c("ENSEMBL", "SYMBOL")
)

gene_map <- gene_map[!is.na(gene_map$SYMBOL), ]
gene_map <- gene_map[!duplicated(gene_map$ENSEMBL), ]

expr_df <- as.data.frame(expr_sub)
expr_df$ENSEMBL <- ensembl_ids_clean

expr_df_annot <- merge(expr_df, gene_map, by = "ENSEMBL")

# If multiple Ensembl IDs map to the same gene symbol, keep the row with highest mean expression
expr_df_annot$mean_expr <- rowMeans(expr_df_annot[, colnames(expr_sub), drop = FALSE])

expr_df_annot <- expr_df_annot %>%
  group_by(SYMBOL) %>%
  slice_max(order_by = mean_expr, n = 1, with_ties = FALSE) %>%
  ungroup()

expr_symbol_mat <- as.matrix(expr_df_annot[, colnames(expr_sub), drop = FALSE])
rownames(expr_symbol_mat) <- expr_df_annot$SYMBOL

## 6) Keep leading-edge genes present after mapping

heatmap_genes_present <- intersect(heatmap_genes, rownames(expr_symbol_mat))

expr_heatmap <- expr_symbol_mat[heatmap_genes_present, , drop = FALSE]

## 7) Scale expression by gene

expr_heatmap_scaled <- t(scale(t(expr_heatmap)))
expr_heatmap_scaled <- expr_heatmap_scaled[complete.cases(expr_heatmap_scaled), ]

## 8) Order samples

meta_sub$Replicate <- as.character(meta_sub$Replicate)

order_idx <- order(meta_sub$Condition == "Scrambled_siRNA", as.numeric(meta_sub$Replicate))

expr_heatmap_scaled <- expr_heatmap_scaled[, order_idx, drop = FALSE]
meta_sub <- meta_sub[order_idx, , drop = FALSE]

## 9) Create display labels

display_condition <- ifelse(
  meta_sub$Condition == "NgR3_siRNA",
  "NgR3 siRNA",
  ifelse(
    meta_sub$Condition == "Scrambled_siRNA",
    "Scrambled siRNA",
    as.character(meta_sub$Condition)
  )
)

colnames(expr_heatmap_scaled) <- paste(display_condition, meta_sub$Replicate)

## 10) Create column annotation

col_annot <- data.frame(
  Condition = factor(
    display_condition,
    levels = c("Scrambled siRNA", "NgR3 siRNA")
  )
)

rownames(col_annot) <- colnames(expr_heatmap_scaled)

## 11) Create row annotation

row_annot <- data.frame(
  Pathway = factor(
    rep("E2F targets", nrow(expr_heatmap_scaled)),
    levels = c("E2F targets")
  )
)

rownames(row_annot) <- rownames(expr_heatmap_scaled)

## 12) Set annotation colours

annotation_colors <- list(
  Condition = c(
    "Scrambled siRNA" = "grey70",
    "NgR3 siRNA" = "black"
  ),
  Pathway = c(
    "E2F targets" = "#8A2BE2"
  )
)

## 13) Set heatmap colour scale

breaks <- seq(-2.5, 2.5, length.out = 101)

## 14) Draw heatmap

p <- pheatmap(
  expr_heatmap_scaled,
  annotation_col = col_annot,
  annotation_row = row_annot,
  annotation_colors = annotation_colors,
  cluster_rows = TRUE,
  treeheight_row = 0,
  cluster_cols = FALSE,
  treeheight_col = 0,
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize_row = 8,
  fontsize_col = 10,
  scale = "none",
  breaks = breaks,
  color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
  border_color = "grey60",
  silent = TRUE
)

## 15) Save heatmap

png(
  "NgR3_vs_Scrambled_upregulated_hallmark_pathways_heatmap.png",
  width = 2800,
  height = 1700,
  res = 300
)

grid.newpage()

pushViewport(viewport(
  x = 0.5,
  y = 0.48,
  width = 0.96,
  height = 0.90
))

grid.draw(p$gtable)
popViewport()

dev.off()


## Hallmark GSEA heatmap: NgR3 siRNA vs scrambled siRNA downregulated
## Leading-edge genes from IL2 STAT5 signaling, IL6 JAK STAT3 signaling, inflammatory response, and UV response

## 1) Select pathways

selected_pathways <- c(
  "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_UV_RESPONSE_UP",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_ESTROGEN_RESPONSE_EARLY"
)

## 2) Extract leading-edge genes

fgsea_sub <- as.data.frame(fgseaRes_NgR3_vs_Scrambled) %>%
  dplyr::filter(pathway %in% selected_pathways)

leading_edge_df <- fgsea_sub %>%
  dplyr::select(pathway, leadingEdge) %>%
  tidyr::unnest_longer(leadingEdge) %>%
  dplyr::rename(symbol = leadingEdge)

## 3) Add DESeq2 statistics

res_df <- as.data.frame(res_NgR3_vs_Scrambled) %>%
  rownames_to_column("ensembl")

res_df$ensembl_clean <- sub("\\..*$", "", res_df$ensembl)

res_df$symbol <- mapIds(
  org.Hs.eg.db,
  keys = res_df$ensembl_clean,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

res_df <- res_df %>%
  dplyr::filter(!is.na(symbol), !is.na(log2FoldChange)) %>%
  dplyr::group_by(symbol) %>%
  dplyr::slice_max(order_by = abs(log2FoldChange), n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(symbol, log2FoldChange, padj)

leading_edge_df <- leading_edge_df %>%
  dplyr::left_join(res_df, by = "symbol") %>%
  dplyr::filter(!is.na(log2FoldChange))

## 4) Keep top 15 unique genes per pathway

leading_edge_df <- leading_edge_df %>%
  dplyr::mutate(
    pathway = factor(pathway, levels = selected_pathways)
  )

top_n <- 15
used_symbols <- character(0)
pathway_gene_list <- list()

for (pw in selected_pathways) {
  
  pw_df <- leading_edge_df %>%
    dplyr::filter(pathway == pw, !symbol %in% used_symbols) %>%
    dplyr::arrange(dplyr::desc(abs(log2FoldChange))) %>%
    dplyr::slice_head(n = top_n)
  
  pathway_gene_list[[pw]] <- pw_df
  used_symbols <- c(used_symbols, pw_df$symbol)
}

top_genes_df <- dplyr::bind_rows(pathway_gene_list)

## 5) Add pathway labels

top_genes_df <- top_genes_df %>%
  dplyr::mutate(
    pathway_label = dplyr::case_when(
      pathway == "HALLMARK_IL2_STAT5_SIGNALING" ~ "IL2 STAT5 signalling",
      pathway == "HALLMARK_IL6_JAK_STAT3_SIGNALING" ~ "IL6 JAK STAT3 signalling",
      pathway == "HALLMARK_INFLAMMATORY_RESPONSE" ~ "Inflammatory response",
      pathway == "HALLMARK_UV_RESPONSE_UP" ~ "UV response",
      pathway == "HALLMARK_TNFA_SIGNALING_VIA_NFKB" ~ "TNFα/NFκB signalling",
      pathway == "HALLMARK_ESTROGEN_RESPONSE_EARLY" ~ "Estrogen response early",
      TRUE ~ as.character(pathway)
    )
  )

## 6) Build expression matrix from VST object

vsd <- vst(dds_full, blind = FALSE)
expr_mat <- assay(vsd)

samples_keep <- colData(dds_full)$Condition %in% c("NgR3_siRNA", "Scrambled_siRNA")

expr_sub <- expr_mat[, samples_keep, drop = FALSE]
meta_sub <- as.data.frame(colData(dds_full)[samples_keep, ])

## 7) Map Ensembl IDs to gene symbols

ensembl_ids <- rownames(expr_sub)
ensembl_ids_clean <- sub("\\..*$", "", ensembl_ids)

gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = ensembl_ids_clean,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

expr_df <- as.data.frame(expr_sub) %>%
  rownames_to_column("ensembl")

expr_df$ensembl_clean <- sub("\\..*$", "", expr_df$ensembl)
expr_df$symbol <- gene_symbols

expr_df <- expr_df %>%
  dplyr::filter(!is.na(symbol))

# If multiple Ensembl IDs map to the same symbol, keep strongest DESeq2 result
expr_df <- expr_df %>%
  dplyr::left_join(
    res_df %>% dplyr::select(symbol, log2FoldChange),
    by = "symbol"
  ) %>%
  dplyr::group_by(symbol) %>%
  dplyr::slice_max(order_by = abs(log2FoldChange), n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

selected_symbols <- top_genes_df$symbol

expr_mat2 <- expr_df %>%
  dplyr::filter(symbol %in% selected_symbols) %>%
  dplyr::select(symbol, colnames(expr_sub)) %>%
  tibble::column_to_rownames("symbol") %>%
  as.matrix()

## 8) Scale expression by gene

expr_mat2_scaled <- t(scale(t(expr_mat2)))
expr_mat2_scaled[is.na(expr_mat2_scaled)] <- 0

## 9) Order samples

meta_sub$Replicate <- as.character(meta_sub$Replicate)

order_idx <- order(
  meta_sub$Condition != "Scrambled_siRNA",
  as.numeric(meta_sub$Replicate)
)

expr_mat2_scaled <- expr_mat2_scaled[, order_idx, drop = FALSE]
meta_sub <- meta_sub[order_idx, , drop = FALSE]

display_condition <- ifelse(
  meta_sub$Condition == "Scrambled_siRNA",
  "Scrambled siRNA",
  ifelse(
    meta_sub$Condition == "NgR3_siRNA",
    "NgR3 siRNA",
    as.character(meta_sub$Condition)
  )
)

colnames(expr_mat2_scaled) <- paste(display_condition, meta_sub$Replicate)

## 10) Cluster genes within each pathway

row_order_list <- list()

for (pw_label in c(
  "IL2 STAT5 signalling",
  "IL6 JAK STAT3 signalling",
  "Inflammatory response",
  "UV response",
  "TNFα/NFκB signalling",
  "Estrogen response early"
)) {
  
  genes_this_pathway <- top_genes_df %>%
    dplyr::filter(pathway_label == pw_label) %>%
    dplyr::pull(symbol)
  
  genes_this_pathway <- intersect(genes_this_pathway, rownames(expr_mat2_scaled))
  
  if (length(genes_this_pathway) > 1) {
    
    mat_sub <- expr_mat2_scaled[genes_this_pathway, , drop = FALSE]
    hc <- hclust(dist(mat_sub), method = "complete")
    ordered_genes <- genes_this_pathway[hc$order]
    
  } else {
    
    ordered_genes <- genes_this_pathway
  }
  
  row_order_list[[pw_label]] <- ordered_genes
}

row_order <- unlist(row_order_list, use.names = FALSE)

expr_mat2_scaled <- expr_mat2_scaled[row_order, , drop = FALSE]

## 11) Create row annotation

annotation_row <- top_genes_df %>%
  dplyr::select(symbol, Pathway = pathway_label) %>%
  dplyr::filter(symbol %in% row_order) %>%
  dplyr::distinct(symbol, .keep_all = TRUE)

annotation_row <- annotation_row[match(row_order, annotation_row$symbol), ]
annotation_row <- as.data.frame(annotation_row)
rownames(annotation_row) <- annotation_row$symbol
annotation_row$symbol <- NULL

annotation_row$Pathway <- factor(
  annotation_row$Pathway,
  levels = c(
    "IL2 STAT5 signalling",
    "IL6 JAK STAT3 signalling",
    "Inflammatory response",
    "UV response",
    "TNFα/NFκB signalling",
    "Estrogen response early"
  )
)

## 12) Create column annotation

annotation_col <- data.frame(
  Condition = factor(
    display_condition,
    levels = c("Scrambled siRNA", "NgR3 siRNA")
  )
)

rownames(annotation_col) <- colnames(expr_mat2_scaled)

## 13) Set annotation colours

annotation_colors <- list(
  Condition = c(
    "Scrambled siRNA" = "black",
    "NgR3 siRNA" = "grey70"
  ),
  Pathway = c(
    "IL2 STAT5 signalling" = "#7B2CBF",
    "IL6 JAK STAT3 signalling" = "#B8A1E3",
    "Inflammatory response" = "#2E7D32",
    "UV response" = "#2A9D8F",
    "TNFα/NFκB signalling" = "#DB7093",
    "Estrogen response early" = "#F4A261"
  )
)

## 14) Set heatmap colour scale

breaks <- seq(-2.5, 2.5, length.out = 101)

heatmap_colours <- colorRampPalette(
  rev(brewer.pal(n = 11, name = "RdBu"))
)(100)

## 15) Draw heatmap

ph <- pheatmap(
  expr_mat2_scaled,
  annotation_col = annotation_col,
  annotation_row = annotation_row,
  annotation_colors = annotation_colors,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize_row = 5,
  fontsize_col = 10,
  border_color = "grey60",
  scale = "none",
  breaks = breaks,
  color = heatmap_colours,
  silent = TRUE
)

grid::grid.newpage()
grid::grid.draw(ph$gtable)

## 16) Save heatmap

png(
  "Scrambled_vs_NgR3_downregulated_hallmark_pathways_heatmap.png",
  width = 3200,
  height = 2200,
  res = 300
)

grid::grid.newpage()

grid::pushViewport(
  grid::viewport(
    x = 0.53,
    y = 0.50,
    width = 0.90,
    height = 0.88
  )
)

grid::grid.draw(ph$gtable)
grid::popViewport()

dev.off()



## C2 GSEA analysis: NgR3 siRNA vs scrambled siRNA

## 1) Create comparison factor

meta$Condition <- sub("_(\\d+)$", "", meta$Group)
meta$comparison_group <- NA_character_

# NgR3 siRNA samples
meta$comparison_group[
  grepl("NgR3", meta$Condition, ignore.case = TRUE)
] <- "NgR3_siRNA"

# Scrambled siRNA samples
meta$comparison_group[
  grepl("Scrambled", meta$Condition, ignore.case = TRUE)
] <- "Scrambled_siRNA"

meta$comparison_group <- factor(
  meta$comparison_group,
  levels = c("Scrambled_siRNA", "NgR3_siRNA")
)

## 2) Keep only samples for this comparison

keep <- !is.na(meta$comparison_group)

meta_sub  <- meta[keep, , drop = FALSE]
count_sub <- count_matrix[, rownames(meta_sub), drop = FALSE]

print(table(meta_sub$comparison_group))

## 3) Build DESeq2 dataset

dds_NgR3_vs_Scrambled <- DESeqDataSetFromMatrix(
  countData = count_sub,
  colData   = meta_sub,
  design    = ~ comparison_group
)

dds_NgR3_vs_Scrambled <- dds_NgR3_vs_Scrambled[
  rowSums(counts(dds_NgR3_vs_Scrambled)) > 10,
]

dds_NgR3_vs_Scrambled <- DESeq(dds_NgR3_vs_Scrambled)

resultsNames(dds_NgR3_vs_Scrambled)

## 4) Get DESeq2 results

res_NgR3_vs_Scrambled <- results(
  dds_NgR3_vs_Scrambled,
  contrast = c("comparison_group", "NgR3_siRNA", "Scrambled_siRNA")
)

## 5) Create ranked gene list

gene_list <- res_NgR3_vs_Scrambled$log2FoldChange
names(gene_list) <- rownames(res_NgR3_vs_Scrambled)

gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = names(gene_list),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

gene_df <- data.frame(
  symbol = gene_symbols,
  log2FoldChange = gene_list,
  stringsAsFactors = FALSE
)

# Remove missing symbols
gene_df <- gene_df[!is.na(gene_df$symbol), ]

# Remove duplicate symbols
gene_df <- gene_df %>%
  group_by(symbol) %>%
  slice_max(order_by = abs(log2FoldChange), n = 1, with_ties = FALSE) %>%
  ungroup()

# Convert back to named vector
gene_list <- gene_df$log2FoldChange
names(gene_list) <- gene_df$symbol
gene_list <- sort(gene_list, decreasing = TRUE)

## 6) Load C2 gene sets

msig_c2 <- msigdbr(
  species = "Homo sapiens",
  category = "C2"
)

pathways_c2 <- msig_c2 %>%
  split(x = .$gene_symbol, f = .$gs_name)

## 7) Run fgsea

fgseaRes_c2_NgR3_vs_Scrambled <- fgsea(
  pathways = pathways_c2,
  stats = gene_list,
  minSize = 15,
  maxSize = 500
)

fgseaRes_c2_NgR3_vs_Scrambled <- fgseaRes_c2_NgR3_vs_Scrambled %>%
  arrange(padj)

fgseaRes_c2_NgR3_vs_Scrambled[1:20, ]

## 8) Save C2 GSEA results

saveRDS(
  fgseaRes_c2_NgR3_vs_Scrambled,
  "~/fgseaRes_c2_NgR3_vs_Scrambled.rds"
)

fgsea_export_c2_NgR3_vs_Scrambled <- fgseaRes_c2_NgR3_vs_Scrambled %>%
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ", "))

saveRDS(
  fgsea_export_c2_NgR3_vs_Scrambled,
  file = "~/fgsea_export_c2_NgR3_vs_Scrambled.rds"
)

write.csv(
  fgsea_export_c2_NgR3_vs_Scrambled,
  "gsea_results_C2_NgR3_siRNA_vs_Scrambled_siRNA.csv",
  row.names = FALSE
)


## C2 GSEA dot plot: NgR3 siRNA vs scrambled siRNA

## 1) Select curated C2 pathways

pathways_selected <- c(
  "REACTOME_INTERLEUKIN_10_SIGNALING",
  "TIAN_TNF_SIGNALING_VIA_NFKB",
  "HECKER_IFNB1_TARGETS",
  "KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "KOHN_EMT_EPITHELIAL",
  "ONDER_CDH1_TARGETS_2_DN",
  "WP_DNA_REPLICATION",
  "REACTOME_DNA_STRAND_ELONGATION",
  "KEGG_MEDICUS_REFERENCE_ORIGIN_UNWINDING_AND_ELONGATION",
  "REACTOME_HDR_THROUGH_HOMOLOGOUS_RECOMBINATION_HRR",
  "KEGG_RIBOSOME",
  "REACTOME_NONSENSE_MEDIATED_DECAY_NMD"
)

## 2) Subset fgsea results

plot_df <- fgseaRes_c2_NgR3_vs_Scrambled %>%
  as.data.frame() %>%
  filter(pathway %in% pathways_selected)

## 3) Clean pathway labels

plot_df$pathway_clean <- NA_character_

plot_df$pathway_clean[plot_df$pathway == "REACTOME_INTERLEUKIN_10_SIGNALING"] <- "Interleukin-10 signalling"
plot_df$pathway_clean[plot_df$pathway == "TIAN_TNF_SIGNALING_VIA_NFKB"] <- "TNF/NFκB signalling"
plot_df$pathway_clean[plot_df$pathway == "HECKER_IFNB1_TARGETS"] <- "IFNβ targets"
plot_df$pathway_clean[plot_df$pathway == "KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY"] <- "NOD-like receptor signalling"
plot_df$pathway_clean[plot_df$pathway == "KOHN_EMT_EPITHELIAL"] <- "EMT (epithelial)"
plot_df$pathway_clean[plot_df$pathway == "ONDER_CDH1_TARGETS_2_DN"] <- "CDH1 targets"
plot_df$pathway_clean[plot_df$pathway == "WP_DNA_REPLICATION"] <- "DNA replication"
plot_df$pathway_clean[plot_df$pathway == "REACTOME_DNA_STRAND_ELONGATION"] <- "DNA strand elongation"
plot_df$pathway_clean[plot_df$pathway == "KEGG_MEDICUS_REFERENCE_ORIGIN_UNWINDING_AND_ELONGATION"] <- "Origin unwinding"
plot_df$pathway_clean[plot_df$pathway == "REACTOME_HDR_THROUGH_HOMOLOGOUS_RECOMBINATION_HRR"] <- "Homologous recombination"
plot_df$pathway_clean[plot_df$pathway == "KEGG_RIBOSOME"] <- "Ribosome"
plot_df$pathway_clean[plot_df$pathway == "REACTOME_NONSENSE_MEDIATED_DECAY_NMD"] <- "Nonsense-mediated decay"

## 4) Add enrichment direction

plot_df$direction <- ifelse(
  plot_df$NES > 0,
  "Gained with NgR3 knockdown",
  "Lost with NgR3 knockdown"
)

## 5) Set pathway order

plot_df$pathway_clean <- base::factor(
  plot_df$pathway_clean,
  levels = c(
    "Interleukin-10 signalling",
    "TNF/NFκB signalling",
    "IFNβ targets",
    "NOD-like receptor signalling",
    "EMT (epithelial)",
    "CDH1 targets",
    "DNA replication",
    "DNA strand elongation",
    "Origin unwinding",
    "Homologous recombination",
    "Nonsense-mediated decay",
    "Ribosome"
  )
)

## 6) Create dot plot

dotplot <- ggplot(plot_df, aes(x = NES, y = pathway_clean)) +
  geom_point(aes(size = size, colour = -log10(padj))) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~direction, scales = "free_y", ncol = 1) +
  scale_colour_gradientn(
    colours = c("deepskyblue", "yellow", "red")
  ) +
  scale_size(range = c(4, 12)) +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = "Pathway",
    colour = "-log10(padj)",
    size = "Gene set size"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 11),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    strip.text = element_text(face = "bold", size = 12)
  )

dotplot

## 7) Save C2 dot plot

ggsave(
  filename = "GSEA_C2_dotplot_NgR3_siRNA_vs_Scrambled_siRNA.png",
  plot = dotplot,
  width = 10,
  height = 10,
  dpi = 300,
  bg = "white"
)

## C2 heatmap: NgR3 siRNA vs scrambled siRNA apoptosis, survival, stress response, migration and invasion

## 1) Define genes of interest

genes_heatmap1 <- c(
  # Apoptosis / survival
  "TNF", "TNFAIP3", "NFKBIA", "IRF1",
  
  # Inflammatory / stress response
  "IL6", "CCL2", "CXCL1", "CXCL2", "CXCL8", "ICAM1",
  
  # Migration / invasion
  "EPCAM", "ESRP1", "ESRP2", "ITGB6", "LAMC2", "CLDN1", "ST14", "KRT16", "TJP3"
)

## 2) Check metadata groups

cat("\nAll groups in meta$Group:\n")
print(table(meta$Group))

if ("Condition" %in% colnames(meta)) {
  cat("\nAll treatments in meta$Condition:\n")
  print(table(meta$Condition))
}

meta$comparison_group <- meta$Condition

cat("\nAll groups in meta$comparison_group:\n")
print(table(meta$comparison_group))

## 3) Select NgR3 siRNA and scrambled siRNA samples

meta_sub <- meta[meta$comparison_group %in% c("Scrambled_siRNA", "NgR3_siRNA"), , drop = FALSE]

meta_sub$comparison_group <- droplevels(factor(meta_sub$comparison_group))

cat("\nSelected samples:\n")
print(data.frame(
  Sample = rownames(meta_sub),
  Group = meta_sub$Group,
  comparison_group = meta_sub$comparison_group
))

cat("\nSelected sample counts by treatment:\n")
print(table(meta_sub$comparison_group))

if (!all(c("Scrambled_siRNA", "NgR3_siRNA") %in% meta_sub$comparison_group)) {
  stop("Could not find both Scrambled_siRNA and NgR3_siRNA in meta$comparison_group.")
}

if (!all(table(meta_sub$comparison_group) == 3)) {
  stop("Expected exactly 3 biological replicates for each of Scrambled_siRNA and NgR3_siRNA.")
}

## 4) Clean condition labels and order samples

meta_sub$Condition_label <- ifelse(
  meta_sub$comparison_group == "NgR3_siRNA",
  "NgR3 siRNA",
  "Scrambled siRNA"
)

meta_sub$Condition_label <- factor(
  meta_sub$Condition_label,
  levels = c("NgR3 siRNA", "Scrambled siRNA")
)

meta_sub <- meta_sub[order(meta_sub$Condition_label), , drop = FALSE]

cat("\nFinal ordered samples for heatmap:\n")
print(data.frame(
  Sample = rownames(meta_sub),
  Condition = meta_sub$Condition_label
))

## 5) Get VST expression matrix

expr_mat <- assay(vsd)

cat("\nFull vst matrix dimensions:\n")
print(dim(expr_mat))

expr_mat_sub <- expr_mat[, rownames(meta_sub), drop = FALSE]

cat("\nSubset vst matrix dimensions:\n")
print(dim(expr_mat_sub))

cat("\nFinal heatmap sample names:\n")
print(colnames(expr_mat_sub))

if (ncol(expr_mat_sub) != 6) {
  stop("Expected 6 samples in expr_mat_sub after subsetting.")
}

## 6) Map Ensembl IDs to gene symbols

gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = rownames(expr_mat_sub),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

## 7) Build expression data frame

expr_df <- as.data.frame(expr_mat_sub, stringsAsFactors = FALSE)
expr_df$symbol <- gene_symbols
expr_df$ENSEMBL <- rownames(expr_df)

expr_df1 <- expr_df[expr_df$symbol %in% genes_heatmap1, , drop = FALSE]

cat("\nGenes found before deduplication:\n")
print(expr_df1$symbol)

## 8) Remove duplicate gene symbols

numeric_cols <- colnames(expr_mat_sub)

expr_df1$avg_expr <- rowMeans(expr_df1[, numeric_cols, drop = FALSE], na.rm = TRUE)
expr_df1 <- expr_df1[order(expr_df1$symbol, -expr_df1$avg_expr), ]
expr_df1 <- expr_df1[!duplicated(expr_df1$symbol), ]

cat("\nGenes retained after deduplication:\n")
print(expr_df1$symbol)

cat("\nGenes requested but not found:\n")
print(setdiff(genes_heatmap1, expr_df1$symbol))

## 9) Convert back to matrix

rownames(expr_df1) <- expr_df1$symbol
expr_mat1 <- as.matrix(expr_df1[, numeric_cols, drop = FALSE])

cat("\nexpr_mat1 dimensions before reordering:\n")
print(dim(expr_mat1))

cat("\nexpr_mat1 row names:\n")
print(rownames(expr_mat1))

## 10) Order genes

genes_heatmap1_ordered <- c(
  # Apoptosis / survival
  "TNF", "TNFAIP3", "NFKBIA", "IRF1",
  
  # Inflammatory / stress response
  "IL6", "CCL2", "CXCL1", "CXCL2", "CXCL8", "ICAM1",
  
  # Migration / invasion
  "EPCAM", "ESRP1", "ESRP2", "ITGB6", "LAMC2", "CLDN1", "ST14", "KRT16", "TJP3"
)

genes_heatmap1_present <- genes_heatmap1_ordered[genes_heatmap1_ordered %in% rownames(expr_mat1)]

cat("\nGenes retained for final heatmap in requested order:\n")
print(genes_heatmap1_present)

if (length(genes_heatmap1_present) < 2) {
  stop("Fewer than 2 genes remain after matching. Cannot cluster rows.")
}

expr_mat1 <- expr_mat1[genes_heatmap1_present, , drop = FALSE]

## 11) Reorder columns to match metadata

expr_mat1 <- expr_mat1[, rownames(meta_sub), drop = FALSE]

cat("\nexpr_mat1 dimensions after row/column ordering:\n")
print(dim(expr_mat1))

## 12) Scale expression by gene

expr_mat1_scaled <- t(scale(t(expr_mat1)))

keep_rows <- apply(expr_mat1_scaled, 1, function(x) all(!is.na(x)))
expr_mat1_scaled <- expr_mat1_scaled[keep_rows, , drop = FALSE]

genes_heatmap1_present <- rownames(expr_mat1_scaled)

cat("\nexpr_mat1_scaled dimensions:\n")
print(dim(expr_mat1_scaled))

cat("\nFinal row names in heatmap:\n")
print(rownames(expr_mat1_scaled))

cat("\nFinal column names in heatmap:\n")
print(colnames(expr_mat1_scaled))

if (nrow(expr_mat1_scaled) < 2) {
  stop("Fewer than 2 genes remain after scaling. Cannot cluster rows.")
}

if (ncol(expr_mat1_scaled) < 2) {
  stop("Fewer than 2 samples remain after subsetting. Cannot plot heatmap.")
}

## 13) Create column annotation

annotation_col <- data.frame(
  Condition = meta_sub$Condition_label,
  stringsAsFactors = FALSE
)

rownames(annotation_col) <- rownames(meta_sub)

annotation_col$Condition <- factor(
  annotation_col$Condition,
  levels = c("NgR3 siRNA", "Scrambled siRNA")
)

ann_colors <- list(
  Condition = c(
    "NgR3 siRNA" = "black",
    "Scrambled siRNA" = "grey70"
  )
)

## 14) Create row annotation

annotation_row <- data.frame(
  Phenotype = ifelse(
    genes_heatmap1_present %in% c("TNF", "TNFAIP3", "NFKBIA", "IRF1"),
    "Apoptosis / survival",
    ifelse(
      genes_heatmap1_present %in% c("IL6", "CCL2", "CXCL1", "CXCL2", "CXCL8", "ICAM1"),
      "Inflammatory / stress response",
      ifelse(
        genes_heatmap1_present %in% c("EPCAM", "ESRP1", "ESRP2", "ITGB6", "LAMC2", "CLDN1", "ST14", "KRT14", "KRT16", "TJP3"),
        "Migration / invasion",
        NA
      )
    )
  )
)

rownames(annotation_row) <- genes_heatmap1_present

row_ann_colors <- list(
  Phenotype = c(
    "Apoptosis / survival" = "#800080",
    "Inflammatory / stress response" = "#DB7093",
    "Migration / invasion" = "#2E8B57"
  )
)

all_ann_colors <- c(ann_colors, row_ann_colors)

## 15) Plot heatmap

pheatmap(
  expr_mat1_scaled,
  annotation_col = annotation_col,
  annotation_row = annotation_row,
  annotation_colors = all_ann_colors,
  show_rownames = TRUE,
  labels_row = rownames(expr_mat1_scaled),
  show_colnames = FALSE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_row = 12
)

## 16) Save heatmap

png(
  "NgR3_siRNA_vs_Scrambled_siRNA_Heatmap_apoptosis_stress_migration.png",
  width = 2400,
  height = 1800,
  res = 220
)

pheatmap(
  expr_mat1_scaled,
  annotation_col = annotation_col,
  annotation_row = annotation_row,
  annotation_colors = all_ann_colors,
  show_rownames = TRUE,
  labels_row = rownames(expr_mat1_scaled),
  show_colnames = FALSE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_row = 12
)

dev.off()


## C2 heatmap: NgR3 siRNA vs scrambled siRNA proliferation and replication

## 1) Define genes of interest

genes_heatmap2_ordered <- c(
  # Licensing / origin firing
  "ORC1",
  "CDT1",
  "CDC7",
  
  # Helicase / origin unwinding
  "MCM2",
  "MCM5",
  "MCM6",
  "MCM7",
  
  # DNA synthesis
  "POLA1",
  "PRIM1",
  "PRIM2",
  "PCNA",
  "GINS2",
  
  # Fork stability / elongation
  "RPA1",
  "RPA2",
  "RFC4",
  
  # Replication-associated repair
  "DNA2",
  "RAD51",
  "BRCA1",
  "BRCA2",
  "EXO1"
)

## 2) Check metadata groups

cat("\nAll groups in meta$Group:\n")
print(table(meta$Group))

if ("Condition" %in% colnames(meta)) {
  cat("\nAll treatments in meta$Condition:\n")
  print(table(meta$Condition))
}

if (!"comparison_group" %in% colnames(meta)) {
  meta$comparison_group <- meta$Condition
}

cat("\nAll groups in meta$comparison_group:\n")
print(table(meta$comparison_group))

## 3) Select NgR3 siRNA and scrambled siRNA samples

meta_sub <- meta[meta$comparison_group %in% c("Scrambled_siRNA", "NgR3_siRNA"), , drop = FALSE]

meta_sub$comparison_group <- droplevels(factor(meta_sub$comparison_group))

cat("\nSelected samples:\n")
print(data.frame(
  Sample = rownames(meta_sub),
  Group = meta_sub$Group,
  comparison_group = meta_sub$comparison_group
))

cat("\nSelected sample counts by treatment:\n")
print(table(meta_sub$comparison_group))

if (!all(c("Scrambled_siRNA", "NgR3_siRNA") %in% meta_sub$comparison_group)) {
  stop("Could not find both Scrambled_siRNA and NgR3_siRNA in meta$comparison_group.")
}

if (!all(table(meta_sub$comparison_group) == 3)) {
  stop("Expected exactly 3 biological replicates for each of Scrambled_siRNA and NgR3_siRNA.")
}

## 4) Clean condition labels and order samples

meta_sub$Condition_label <- ifelse(
  meta_sub$comparison_group == "NgR3_siRNA",
  "NgR3 siRNA",
  "Scrambled siRNA"
)

meta_sub$Condition_label <- factor(
  meta_sub$Condition_label,
  levels = c("NgR3 siRNA", "Scrambled siRNA")
)

meta_sub <- meta_sub[order(meta_sub$Condition_label), , drop = FALSE]

cat("\nFinal ordered samples for heatmap:\n")
print(data.frame(
  Sample = rownames(meta_sub),
  Condition = meta_sub$Condition_label
))

## 5) Get VST expression matrix

expr_mat <- assay(vsd)

expr_mat_sub <- expr_mat[, rownames(meta_sub), drop = FALSE]

cat("\nSubset vst matrix dimensions:\n")
print(dim(expr_mat_sub))

cat("\nFinal heatmap sample names:\n")
print(colnames(expr_mat_sub))

if (ncol(expr_mat_sub) != 6) {
  stop("Expected 6 samples in expr_mat_sub after subsetting.")
}

## 6) Map Ensembl IDs to gene symbols

gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = rownames(expr_mat_sub),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

## 7) Build expression data frame

expr_df <- as.data.frame(expr_mat_sub, stringsAsFactors = FALSE)
expr_df$symbol <- gene_symbols
expr_df$ENSEMBL <- rownames(expr_df)

expr_df2 <- expr_df[expr_df$symbol %in% genes_heatmap2_ordered, , drop = FALSE]

cat("\nGenes found before deduplication:\n")
print(expr_df2$symbol)

## 8) Remove duplicate gene symbols

numeric_cols <- colnames(expr_mat_sub)

expr_df2$avg_expr <- rowMeans(expr_df2[, numeric_cols, drop = FALSE], na.rm = TRUE)
expr_df2 <- expr_df2[order(expr_df2$symbol, -expr_df2$avg_expr), ]
expr_df2 <- expr_df2[!duplicated(expr_df2$symbol), ]

cat("\nGenes retained after deduplication:\n")
print(expr_df2$symbol)

cat("\nGenes requested but not found:\n")
print(setdiff(genes_heatmap2_ordered, expr_df2$symbol))

## 9) Convert back to matrix

rownames(expr_df2) <- expr_df2$symbol
expr_mat2 <- as.matrix(expr_df2[, numeric_cols, drop = FALSE])

cat("\nexpr_mat2 dimensions before row ordering:\n")
print(dim(expr_mat2))

cat("\nexpr_mat2 row names:\n")
print(rownames(expr_mat2))

genes_heatmap2_present <- genes_heatmap2_ordered[genes_heatmap2_ordered %in% rownames(expr_mat2)]

cat("\nGenes retained for final heatmap in requested order:\n")
print(genes_heatmap2_present)

if (length(genes_heatmap2_present) < 2) {
  stop("Fewer than 2 genes remain after matching. Cannot cluster rows.")
}

expr_mat2 <- expr_mat2[genes_heatmap2_present, , drop = FALSE]

## 10) Reorder columns to match metadata

expr_mat2 <- expr_mat2[, rownames(meta_sub), drop = FALSE]

cat("\nexpr_mat2 dimensions after row/column ordering:\n")
print(dim(expr_mat2))

## 11) Scale expression by gene

expr_mat2_scaled <- t(scale(t(expr_mat2)))

keep_rows <- apply(expr_mat2_scaled, 1, function(x) all(!is.na(x)))
expr_mat2_scaled <- expr_mat2_scaled[keep_rows, , drop = FALSE]

genes_heatmap2_present <- rownames(expr_mat2_scaled)

cat("\nexpr_mat2_scaled dimensions:\n")
print(dim(expr_mat2_scaled))

cat("\nFinal row names in heatmap:\n")
print(rownames(expr_mat2_scaled))

cat("\nFinal column names in heatmap:\n")
print(colnames(expr_mat2_scaled))

if (nrow(expr_mat2_scaled) < 2) {
  stop("Fewer than 2 genes remain after scaling. Cannot cluster rows.")
}

if (ncol(expr_mat2_scaled) < 2) {
  stop("Fewer than 2 samples remain after subsetting. Cannot plot heatmap.")
}

## 12) Create column annotation

annotation_col <- data.frame(
  Condition = meta_sub$Condition_label,
  stringsAsFactors = FALSE
)

rownames(annotation_col) <- rownames(meta_sub)

annotation_col$Condition <- factor(
  annotation_col$Condition,
  levels = c("NgR3 siRNA", "Scrambled siRNA")
)

ann_colors <- list(
  Condition = c(
    "NgR3 siRNA" = "black",
    "Scrambled siRNA" = "grey70"
  )
)

## 13) Create row annotation

annotation_row <- data.frame(
  `Functional category` = ifelse(
    genes_heatmap2_present %in% c("ORC1", "CDT1", "CDC7"),
    "Licensing / origin firing",
    ifelse(
      genes_heatmap2_present %in% c("MCM2", "MCM5", "MCM6", "MCM7"),
      "Helicase / origin unwinding",
      ifelse(
        genes_heatmap2_present %in% c("POLA1", "PRIM1", "PRIM2", "PCNA", "GINS2"),
        "DNA synthesis",
        ifelse(
          genes_heatmap2_present %in% c("RPA1", "RPA2", "RFC4"),
          "Fork stability / elongation",
          "Replication-associated repair"
        )
      )
    )
  ),
  stringsAsFactors = FALSE
)

rownames(annotation_row) <- genes_heatmap2_present
colnames(annotation_row) <- "Functional category"

annotation_row$`Functional category` <- trimws(annotation_row$`Functional category`)

annotation_row$`Functional category` <- factor(
  annotation_row$`Functional category`,
  levels = c(
    "Licensing / origin firing",
    "Helicase / origin unwinding",
    "DNA synthesis",
    "Fork stability / elongation",
    "Replication-associated repair"
  )
)

colnames(annotation_row) <- "Functional category"

row_ann_colors <- list(
  `Functional category` = c(
    "Licensing / origin firing" = "#66FF33",
    "Helicase / origin unwinding" = "#FF2D95",
    "DNA synthesis" = "#C000FF",
    "Fork stability / elongation" = "#1B9E77",
    "Replication-associated repair" = "#00FFFF"
  )
)

all_ann_colors <- c(ann_colors, row_ann_colors)

## 14) Plot heatmap

pheatmap(
  expr_mat2_scaled,
  annotation_col = annotation_col,
  annotation_row = annotation_row,
  annotation_colors = all_ann_colors,
  show_rownames = TRUE,
  labels_row = rownames(expr_mat2_scaled),
  show_colnames = FALSE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_row = 12
)

## 15) Save heatmap

ph <- pheatmap(
  expr_mat2_scaled,
  annotation_col = annotation_col,
  annotation_row = annotation_row,
  annotation_colors = all_ann_colors,
  show_rownames = TRUE,
  labels_row = rownames(expr_mat2_scaled),
  show_colnames = FALSE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_row = 12,
  silent = TRUE
)

png(
  "NgR3_siRNA_vs_Scrambled_siRNA_Heatmap_proliferation_replication.png",
  width = 2400,
  height = 2000,
  res = 220
)

grid.newpage()
pushViewport(viewport(y = 0.52, height = 0.95))
grid.draw(ph$gtable)

dev.off()



## DESeq2 interaction model: 1:1,000 EPR11334 vs NgR3 siRNA
## Comparison: (1:1,000 EPR11334 vs IgG) vs (NgR3 siRNA vs scrambled siRNA)

## 1) Prepare count matrix

counts_mat_1000 <- as.matrix(counts_raw[, -(1:6)])
rownames(counts_mat_1000) <- counts_raw$Geneid

head(rownames(counts_mat_1000))

## 2) Clean metadata

meta_clean_1000 <- meta[!(rownames(meta) %in% c("NA", "NA.1", "NA.2", "NA.3", "NA.4", "NA.5")), , drop = FALSE]

meta_clean_1000 <- meta_clean_1000[colnames(counts_mat_1000), , drop = FALSE]

stopifnot(identical(colnames(counts_mat_1000), rownames(meta_clean_1000)))

## 3) Keep groups required for the interaction model

keep_conditions_1000 <- c(
  "IGG",
  "1:1000_EPR",
  "Scrambled_siRNA",
  "NgR3_siRNA"
)

meta_sub_1000 <- meta_clean_1000 %>%
  dplyr::filter(Condition %in% keep_conditions_1000)

counts_sub_1000 <- counts_mat_1000[, rownames(meta_sub_1000), drop = FALSE]

meta_sub_1000 <- meta_sub_1000[colnames(counts_sub_1000), , drop = FALSE]

stopifnot(identical(colnames(counts_sub_1000), rownames(meta_sub_1000)))

head(rownames(counts_sub_1000))

## 4) Create interaction model variables

meta_sub_1000$Perturbation_1000 <- dplyr::case_when(
  meta_sub_1000$Condition %in% c("IGG", "1:1000_EPR") ~ "EPR",
  meta_sub_1000$Condition %in% c("Scrambled_siRNA", "NgR3_siRNA") ~ "siRNA"
)

meta_sub_1000$Treatment_1000 <- dplyr::case_when(
  meta_sub_1000$Condition %in% c("IGG", "Scrambled_siRNA") ~ "Control",
  meta_sub_1000$Condition %in% c("1:1000_EPR", "NgR3_siRNA") ~ "Treated"
)

meta_sub_1000$Perturbation_1000 <- factor(
  meta_sub_1000$Perturbation_1000,
  levels = c("siRNA", "EPR")
)

meta_sub_1000$Treatment_1000 <- factor(
  meta_sub_1000$Treatment_1000,
  levels = c("Control", "Treated")
)

table(
  meta_sub_1000$Condition,
  meta_sub_1000$Perturbation_1000,
  meta_sub_1000$Treatment_1000
)

## 5) Build DESeq2 interaction model

dds_int_1000 <- DESeqDataSetFromMatrix(
  countData = counts_sub_1000,
  colData   = meta_sub_1000,
  design    = ~ Perturbation_1000 + Treatment_1000 + Perturbation_1000:Treatment_1000
)

head(rownames(dds_int_1000))

## 6) Filter low-count genes

keep_genes_1000 <- rowSums(counts(dds_int_1000) >= 10) >= 2
dds_int_1000 <- dds_int_1000[keep_genes_1000, ]

## 7) Run DESeq2

dds_int_1000 <- DESeq(dds_int_1000)

## 8) Inspect coefficient names

resultsNames_1000 <- resultsNames(dds_int_1000)
resultsNames_1000

saveRDS(
  dds_int_1000,
  "/home/yougr345/dds_interaction_model_1000EPR_vs_NgR3.rds"
)

## 9) Extract interaction term

interaction_name_1000 <- grep(
  "Perturbation.*Treatment|Treatment.*Perturbation",
  resultsNames_1000,
  value = TRUE
)

interaction_name_1000

res_interaction_1000 <- results(
  dds_int_1000,
  name = interaction_name_1000
)

summary(res_interaction_1000)

res_interaction_1000_df <- as.data.frame(res_interaction_1000)
res_interaction_1000_df$gene_id <- rownames(res_interaction_1000_df)

head(res_interaction_1000_df)

saveRDS(
  res_interaction_1000,
  "/home/yougr345/res_interaction_1000EPR_vs_NgR3.rds"
)

## 10) Rank genes by interaction Wald statistic

rank_df_1000 <- res_interaction_1000_df %>%
  dplyr::filter(!is.na(stat)) %>%
  dplyr::distinct(gene_id, .keep_all = TRUE)

ranked_genes_1000 <- rank_df_1000$stat
names(ranked_genes_1000) <- rank_df_1000$gene_id

# Positive values are relatively enriched in 1:1,000 EPR11334
# Negative values are relatively enriched in NgR3 siRNA
ranked_genes_1000 <- sort(ranked_genes_1000, decreasing = TRUE)

head(ranked_genes_1000)
tail(ranked_genes_1000)

saveRDS(
  ranked_genes_1000,
  "/home/yougr345/ranked_genes_INTERACTION_1000EPR_vs_NgR3_HALLMARK.rds"
)

## 11) Map Ensembl IDs to gene symbols

if (grepl("^ENSG", names(ranked_genes_1000)[1])) {
  
  ens_ids_1000 <- sub("\\..*$", "", names(ranked_genes_1000))
  
  gene_map_1000 <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys    = unique(ens_ids_1000),
    keytype = "ENSEMBL",
    columns = c("ENSEMBL", "SYMBOL")
  ) %>%
    dplyr::filter(!is.na(SYMBOL), SYMBOL != "") %>%
    dplyr::distinct(ENSEMBL, .keep_all = TRUE)
  
  ranked_df_symbols_1000 <- data.frame(
    ENSEMBL = ens_ids_1000,
    stat    = as.numeric(ranked_genes_1000),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::left_join(gene_map_1000, by = "ENSEMBL") %>%
    dplyr::filter(!is.na(SYMBOL), SYMBOL != "") %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::summarise(stat = stat[which.max(abs(stat))], .groups = "drop")
  
  ranked_genes_h_1000 <- ranked_df_symbols_1000$stat
  names(ranked_genes_h_1000) <- ranked_df_symbols_1000$SYMBOL
  ranked_genes_h_1000 <- sort(ranked_genes_h_1000, decreasing = TRUE)
  
} else {
  ranked_genes_h_1000 <- ranked_genes_1000
}


## Hallmark GSEA: 1:1,000 EPR11334 vs NgR3 siRNA interaction model
## Comparison: (1:1,000 EPR11334 vs IgG) vs (NgR3 siRNA vs scrambled siRNA)

## 1) Load Hallmark gene sets

msig_h_1000 <- msigdbr(
  species  = "Homo sapiens",
  category = "H"
)

pathways_h_1000 <- msig_h_1000 %>%
  split(x = .$gene_symbol, f = .$gs_name)

## 2) Run Hallmark fgsea

set.seed(123)

fgseaRes_interaction_H_1000 <- fgsea(
  pathways = pathways_h_1000,
  stats    = ranked_genes_h_1000,
  minSize  = 15,
  maxSize  = 500
)

fgseaRes_interaction_H_1000 <- fgseaRes_interaction_H_1000 %>%
  dplyr::arrange(padj)

print(head(fgseaRes_interaction_H_1000, 20))

sig_fgseaRes_interaction_H_1000 <- fgseaRes_interaction_H_1000 %>%
  dplyr::filter(padj < 0.05)

print(sig_fgseaRes_interaction_H_1000)

## 3) Save Hallmark GSEA results

fgseaRes_interaction_H_1000_clean <- fgseaRes_interaction_H_1000

fgseaRes_interaction_H_1000_clean$leadingEdge <- sapply(
  fgseaRes_interaction_H_1000_clean$leadingEdge,
  function(x) paste(x, collapse = ", ")
)

write.csv(
  fgseaRes_interaction_H_1000_clean,
  file = "~/fgseaRes_INTERACTION_1000EPR_vs_NgR3_HALLMARK_ALL.csv",
  row.names = FALSE
)

sig_fgseaRes_interaction_H_1000_clean <- fgseaRes_interaction_H_1000_clean %>%
  dplyr::filter(padj < 0.05)

write.csv(
  sig_fgseaRes_interaction_H_1000_clean,
  file = "~/fgseaRes_INTERACTION_1000EPR_vs_NgR3_HALLMARK_SIG.csv",
  row.names = FALSE
)


## Hallmark GSEA dot plot: 1:100,000 EPR11334 vs NgR3 siRNA interaction model
## Comparison: (1:100,000 EPR11334 vs IgG) vs (NgR3 siRNA vs scrambled siRNA)

## 1) Prepare significant Hallmark pathways

plot_df <- fgseaRes_interaction_H_1000_clean %>%
  filter(padj < 0.05) %>%
  mutate(
    pathway_label = pathway %>%
      gsub("^HALLMARK_", "", .) %>%
      gsub("_", " ", .) %>%
      tolower() %>%
      tools::toTitleCase(),
    
    pathway_label = case_when(
      pathway_label == "Tnfa Signaling Via Nfkb" ~ "TNFα signalling via NFκB",
      pathway_label == "Il6 Jak Stat3 Signaling" ~ "IL6–JAK–STAT3 signalling",
      pathway_label == "Il2 Stat5 Signaling" ~ "IL2–STAT5 signalling",
      pathway_label == "Tgf Beta Signaling" ~ "TGF-β signalling",
      pathway_label == "Myc Targets V1" ~ "Myc targets V1",
      pathway_label == "Myc Targets V2" ~ "Myc targets V2",
      pathway_label == "Uv Response Up" ~ "UV response up",
      pathway_label == "Uv Response Dn" ~ "UV response down",
      pathway_label == "Dna Repair" ~ "DNA repair",
      pathway_label == "E2f Targets" ~ "E2F targets",
      pathway_label == "Estrogen Response Early" ~ "Estrogen response early",
      pathway_label == "Estrogen Response Late" ~ "Estrogen response late",
      pathway_label == "Inflammatory Response" ~ "Inflammatory response",
      pathway_label == "Oxidative Phosphorylation" ~ "Oxidative phosphorylation",
      pathway_label == "Mitotic Spindle" ~ "Mitotic spindle",
      pathway_label == "Heme Metabolism" ~ "Heme metabolism",
      pathway_label == "Cholesterol Homeostasis" ~ "Cholesterol homeostasis",
      pathway_label == "P53 Pathway" ~ "p53 pathway",
      pathway_label == "Mtorc1 Signaling" ~ "mTORC1 signalling",
      TRUE ~ pathway_label
    ),
    
    category = case_when(
      pathway_label %in% c(
        "IL6–JAK–STAT3 signalling",
        "IL2–STAT5 signalling",
        "Inflammatory response",
        "TGF-β signalling"
      ) ~ "Inflammation / signalling",
      
      pathway_label %in% c(
        "Estrogen response early",
        "UV response down",
        "Heme metabolism"
      ) ~ "Stress / response",
      
      pathway_label %in% c(
        "mTORC1 signalling"
      ) ~ "Growth / signalling",
      
      pathway_label %in% c(
        "DNA repair",
        "p53 pathway"
      ) ~ "Proliferation / repair",
      
      pathway_label %in% c(
        "Cholesterol homeostasis"
      ) ~ "Metabolism",
      
      TRUE ~ "Other"
    ),
    
    neg_log10_padj = -log10(padj)
  )

## 2) Set pathway order

pathway_order <- c(
  "Inflammatory response",
  "IL6–JAK–STAT3 signalling",
  "UV response down",
  "Estrogen response early",
  "mTORC1 signalling",
  "DNA repair",
  "Cholesterol homeostasis"
)

plot_df <- plot_df %>%
  filter(pathway_label %in% pathway_order)

plot_df$pathway_label <- factor(
  plot_df$pathway_label,
  levels = rev(pathway_order)
)

n_pathways <- nrow(plot_df)
label_y <- n_pathways + 0.75

## 3) Create dot plot

p <- ggplot(plot_df, aes(x = NES, y = pathway_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.7) +
  geom_point(aes(size = size, colour = neg_log10_padj), alpha = 0.95) +
  scale_colour_gradientn(
    colours = c("#2C7BB6", "#FEE08B", "#D7191C"),
    name = "-log10(padj)"
  ) +
  scale_size_continuous(
    name = "Gene set size",
    range = c(4, 13)
  ) +
  scale_x_continuous(
    limits = c(-2.2, 2.2),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_discrete(
    expand = expansion(mult = c(0.03, 0.02))
  ) +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = "Pathway"
  ) +
  annotate(
    "text",
    x = -1.35,
    y = label_y,
    label = "Enriched in NgR3 siRNA",
    fontface = "bold",
    size = 4.8,
    hjust = 0.5
  ) +
  annotate(
    "text",
    x = 1.35,
    y = label_y,
    label = "Enriched in 1,000 EPR11334",
    fontface = "bold",
    size = 4.8,
    hjust = 0.5
  ) +
  coord_cartesian(
    ylim = c(1, n_pathways + 1.0),
    clip = "off"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 20, hjust = 0.5, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 13, hjust = 0.5, margin = margin(b = 18)),
    axis.title.x = element_text(face = "bold", size = 16),
    axis.title.y = element_text(face = "bold", size = 16),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    panel.grid.major.y = element_line(colour = "grey85", linewidth = 0.8),
    panel.grid.major.x = element_line(colour = "grey85", linewidth = 0.8),
    panel.grid.minor = element_blank(),
    legend.title = element_text(face = "bold", size = 13),
    legend.text = element_text(size = 11),
    plot.margin = margin(t = 45, r = 70, b = 35, l = 50)
  )

p

## 4) Save dot plot

png(
  "Hallmark_normalised_dotplot_1000EPR_vs_NgR3siRNA_significant.png",
  width = 3200,
  height = 2400,
  res = 250
)

grid.newpage()

p_grob <- ggplotGrob(p)

pushViewport(
  viewport(
    x = 0.48,
    y = 0.5,
    width = 0.95,
    height = 0.90,
    just = c("center", "center")
  )
)

grid.draw(p_grob)
popViewport()

dev.off()




## DESeq2 interaction model: 1:100,000 EPR11334 vs NgR3 siRNA
## Comparison: (1:100,000 EPR11334 vs IgG) vs (NgR3 siRNA vs scrambled siRNA)

## 1) Prepare count matrix

counts_mat_100000 <- as.matrix(counts_raw[, -(1:6)])
rownames(counts_mat_100000) <- counts_raw$Geneid

head(rownames(counts_mat_100000))

## 2) Clean metadata

meta_clean_100000 <- meta[
  !(rownames(meta) %in% c("NA", "NA.1", "NA.2", "NA.3", "NA.4", "NA.5")),
  ,
  drop = FALSE
]

meta_clean_100000 <- meta_clean_100000[colnames(counts_mat_100000), , drop = FALSE]

stopifnot(identical(colnames(counts_mat_100000), rownames(meta_clean_100000)))

## 3) Keep groups required for the interaction model

keep_conditions_100000 <- c(
  "IGG",
  "1:100000_EPR",
  "Scrambled_siRNA",
  "NgR3_siRNA"
)

meta_sub_100000 <- meta_clean_100000 %>%
  dplyr::filter(Condition %in% keep_conditions_100000)

counts_sub_100000 <- counts_mat_100000[, rownames(meta_sub_100000), drop = FALSE]

meta_sub_100000 <- meta_sub_100000[colnames(counts_sub_100000), , drop = FALSE]

stopifnot(identical(colnames(counts_sub_100000), rownames(meta_sub_100000)))

head(rownames(counts_sub_100000))

## 4) Create interaction model variables

meta_sub_100000$Perturbation_100000 <- dplyr::case_when(
  meta_sub_100000$Condition %in% c("IGG", "1:100000_EPR") ~ "EPR",
  meta_sub_100000$Condition %in% c("Scrambled_siRNA", "NgR3_siRNA") ~ "siRNA"
)

meta_sub_100000$Treatment_100000 <- dplyr::case_when(
  meta_sub_100000$Condition %in% c("IGG", "Scrambled_siRNA") ~ "Control",
  meta_sub_100000$Condition %in% c("1:100000_EPR", "NgR3_siRNA") ~ "Treated"
)

meta_sub_100000$Perturbation_100000 <- factor(
  meta_sub_100000$Perturbation_100000,
  levels = c("siRNA", "EPR")
)

meta_sub_100000$Treatment_100000 <- factor(
  meta_sub_100000$Treatment_100000,
  levels = c("Control", "Treated")
)

table(
  meta_sub_100000$Condition,
  meta_sub_100000$Perturbation_100000,
  meta_sub_100000$Treatment_100000
)

## 5) Build DESeq2 interaction model

dds_int_100000 <- DESeqDataSetFromMatrix(
  countData = counts_sub_100000,
  colData   = meta_sub_100000,
  design    = ~ Perturbation_100000 + Treatment_100000 + Perturbation_100000:Treatment_100000
)

head(rownames(dds_int_100000))

## 6) Filter low-count genes

keep_genes_100000 <- rowSums(counts(dds_int_100000) >= 10) >= 2
dds_int_100000 <- dds_int_100000[keep_genes_100000, ]

## 7) Run DESeq2

dds_int_100000 <- DESeq(dds_int_100000)

## 8) Inspect coefficient names

resultsNames_100000 <- resultsNames(dds_int_100000)
resultsNames_100000

saveRDS(
  dds_int_100000,
  "/home/yougr345/dds_interaction_model_100000EPR_vs_NgR3.rds"
)

## 9) Extract interaction term

interaction_name_100000 <- grep(
  "Perturbation.*Treatment|Treatment.*Perturbation",
  resultsNames_100000,
  value = TRUE
)

interaction_name_100000

stopifnot(length(interaction_name_100000) == 1)

res_interaction_100000 <- results(
  dds_int_100000,
  name = interaction_name_100000
)

summary(res_interaction_100000)

res_interaction_100000_df <- as.data.frame(res_interaction_100000)
res_interaction_100000_df$gene_id <- rownames(res_interaction_100000_df)

head(res_interaction_100000_df)

saveRDS(
  res_interaction_100000,
  "/home/yougr345/res_interaction_100000EPR_vs_NgR3.rds"
)

## 10) Rank genes by interaction Wald statistic

rank_df_100000 <- res_interaction_100000_df %>%
  dplyr::filter(!is.na(stat)) %>%
  dplyr::distinct(gene_id, .keep_all = TRUE)

ranked_genes_100000 <- rank_df_100000$stat
names(ranked_genes_100000) <- rank_df_100000$gene_id

# Positive values are relatively enriched in 1:100,000 EPR11334
# Negative values are relatively enriched in NgR3 siRNA
ranked_genes_100000 <- sort(ranked_genes_100000, decreasing = TRUE)

head(ranked_genes_100000)
tail(ranked_genes_100000)

saveRDS(
  ranked_genes_100000,
  "/home/yougr345/ranked_genes_INTERACTION_100000EPR_vs_NgR3_HALLMARK.rds"
)

## 11) Map Ensembl IDs to gene symbols

if (grepl("^ENSG", names(ranked_genes_100000)[1])) {
  
  ens_ids_100000 <- sub("\\..*$", "", names(ranked_genes_100000))
  
  gene_map_100000 <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys    = unique(ens_ids_100000),
    keytype = "ENSEMBL",
    columns = c("ENSEMBL", "SYMBOL")
  ) %>%
    dplyr::filter(!is.na(SYMBOL), SYMBOL != "") %>%
    dplyr::distinct(ENSEMBL, .keep_all = TRUE)
  
  ranked_df_symbols_100000 <- data.frame(
    ENSEMBL = ens_ids_100000,
    stat    = as.numeric(ranked_genes_100000),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::left_join(gene_map_100000, by = "ENSEMBL") %>%
    dplyr::filter(!is.na(SYMBOL), SYMBOL != "") %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::summarise(stat = stat[which.max(abs(stat))], .groups = "drop")
  
  ranked_genes_h_100000 <- ranked_df_symbols_100000$stat
  names(ranked_genes_h_100000) <- ranked_df_symbols_100000$SYMBOL
  ranked_genes_h_100000 <- sort(ranked_genes_h_100000, decreasing = TRUE)
  
} else {
  ranked_genes_h_100000 <- ranked_genes_100000
}


## Hallmark GSEA: 1:100,000 EPR11334 vs NgR3 siRNA interaction model
## Comparison: (1:100,000 EPR11334 vs IgG) vs (NgR3 siRNA vs scrambled siRNA)

## 1) Load Hallmark gene sets

msig_h_100000 <- msigdbr(
  species  = "Homo sapiens",
  category = "H"
)

pathways_h_100000 <- msig_h_100000 %>%
  split(x = .$gene_symbol, f = .$gs_name)

## 2) Run Hallmark fgsea

fgseaRes_interaction_H_100000 <- fgsea(
  pathways = pathways_h_100000,
  stats    = ranked_genes_h_100000,
  minSize  = 15,
  maxSize  = 500
)

fgseaRes_interaction_H_100000 <- fgseaRes_interaction_H_100000 %>%
  dplyr::arrange(padj)

print(head(fgseaRes_interaction_H_100000, 20))

sig_fgseaRes_interaction_H_100000 <- fgseaRes_interaction_H_100000 %>%
  dplyr::filter(padj < 0.05)

print(sig_fgseaRes_interaction_H_100000)

## 3) Save Hallmark GSEA results

fgseaRes_interaction_H_100000_clean <- fgseaRes_interaction_H_100000

fgseaRes_interaction_H_100000_clean$leadingEdge <- sapply(
  fgseaRes_interaction_H_100000_clean$leadingEdge,
  function(x) paste(x, collapse = ", ")
)

write.csv(
  fgseaRes_interaction_H_100000_clean,
  file = "~/fgseaRes_INTERACTION_100000EPR_vs_NgR3_HALLMARK_ALL.csv",
  row.names = FALSE
)

sig_fgseaRes_interaction_H_100000_clean <- fgseaRes_interaction_H_100000_clean %>%
  dplyr::filter(padj < 0.05)

write.csv(
  sig_fgseaRes_interaction_H_100000_clean,
  file = "~/fgseaRes_INTERACTION_100000EPR_vs_NgR3_HALLMARK_SIG.csv",
  row.names = FALSE
)


## C2 GSEA: 1:100,000 EPR11334 vs NgR3 siRNA interaction model
## Comparison: (1:100,000 EPR11334 vs IgG) vs (NgR3 siRNA vs scrambled siRNA)

## 1) Load C2 gene sets

msig_c2_100000 <- msigdbr(
  species  = "Homo sapiens",
  category = "C2"
)

pathways_c2_100000 <- msig_c2_100000 %>%
  split(x = .$gene_symbol, f = .$gs_name)

## 2) Run C2 fgsea

set.seed(123)

fgseaRes_interaction_C2_100000 <- fgsea(
  pathways = pathways_c2_100000,
  stats    = ranked_genes_h_100000,
  minSize  = 15,
  maxSize  = 500
)

fgseaRes_interaction_C2_100000 <- fgseaRes_interaction_C2_100000 %>%
  dplyr::arrange(padj)

print(head(fgseaRes_interaction_C2_100000, 20))

sig_fgseaRes_interaction_C2_100000 <- fgseaRes_interaction_C2_100000 %>%
  dplyr::filter(padj < 0.05)

print(sig_fgseaRes_interaction_C2_100000)

## 3) Save C2 GSEA results

fgseaRes_interaction_C2_100000_clean <- fgseaRes_interaction_C2_100000

fgseaRes_interaction_C2_100000_clean$leadingEdge <- sapply(
  fgseaRes_interaction_C2_100000_clean$leadingEdge,
  function(x) paste(x, collapse = ", ")
)

write.csv(
  fgseaRes_interaction_C2_100000_clean,
  file = "~/fgseaRes_INTERACTION_100000EPR_vs_NgR3_C2_ALL.csv",
  row.names = FALSE
)

sig_fgseaRes_interaction_C2_100000_clean <- fgseaRes_interaction_C2_100000_clean %>%
  dplyr::filter(padj < 0.05)

write.csv(
  sig_fgseaRes_interaction_C2_100000_clean,
  file = "~/fgseaRes_INTERACTION_100000EPR_vs_NgR3_C2_SIG.csv",
  row.names = FALSE
)



## Hallmark GSEA dot plot: 1:100,000 EPR11334 vs NgR3 siRNA interaction model
## Comparison: (1:100,000 EPR11334 vs IgG) vs (NgR3 siRNA vs scrambled siRNA)

## 1) Prepare significant Hallmark pathways

plot_df <- fgseaRes_interaction_H_100000_clean %>%
  filter(padj < 0.05) %>%
  mutate(
    pathway_label = pathway %>%
      gsub("^HALLMARK_", "", .) %>%
      gsub("_", " ", .) %>%
      tolower() %>%
      tools::toTitleCase(),
    
    pathway_label = case_when(
      pathway_label == "Tnfa Signaling Via Nfkb" ~ "TNFα/NFκB signalling",
      pathway_label == "Il6 Jak Stat3 Signaling" ~ "IL6–JAK–STAT3 signalling",
      pathway_label == "Il2 Stat5 Signaling" ~ "IL2–STAT5 signalling",
      pathway_label == "Tgf Beta Signaling" ~ "TGF-β signalling",
      pathway_label == "Myc Targets V1" ~ "Myc targets V1",
      pathway_label == "Myc Targets V2" ~ "Myc targets V2",
      pathway_label == "Uv Response Up" ~ "UV response up",
      pathway_label == "Uv Response Dn" ~ "UV response down",
      pathway_label == "Dna Repair" ~ "DNA repair",
      pathway_label == "E2f Targets" ~ "E2F targets",
      pathway_label == "Estrogen Response Early" ~ "Estrogen response early",
      pathway_label == "Estrogen Response Late" ~ "Estrogen response late",
      pathway_label == "Inflammatory Response" ~ "Inflammatory response",
      pathway_label == "Oxidative Phosphorylation" ~ "Oxidative phosphorylation",
      pathway_label == "Mitotic Spindle" ~ "Mitotic spindle",
      pathway_label == "Heme Metabolism" ~ "Heme metabolism",
      pathway_label == "Cholesterol Homeostasis" ~ "Cholesterol homeostasis",
      pathway_label == "Peroxisome" ~ "Peroxisome",
      TRUE ~ pathway_label
    ),
    
    category = case_when(
      pathway_label %in% c(
        "TNFα/NFκB signalling",
        "IL6–JAK–STAT3 signalling",
        "Inflammatory response",
        "TGF-β signalling"
      ) ~ "Inflammation / signalling",
      
      pathway_label %in% c(
        "Estrogen response early",
        "UV response up",
        "UV response down",
        "Heme metabolism"
      ) ~ "Stress / response",
      
      pathway_label %in% c(
        "Mitotic spindle"
      ) ~ "Cell division",
      
      pathway_label %in% c(
        "E2F targets",
        "Myc targets V1",
        "DNA repair"
      ) ~ "Proliferation / repair",
      
      pathway_label %in% c(
        "Oxidative phosphorylation",
        "Peroxisome",
        "Cholesterol homeostasis"
      ) ~ "Metabolism",
      
      TRUE ~ "Other"
    ),
    
    neg_log10_padj = -log10(padj)
  )

plot_df$pathway_label[plot_df$pathway == "HALLMARK_TNFA_SIGNALING_VIA_NFKB"] <- "TNFα/NFκB signalling"

## 2) Set pathway order

pathway_order <- c(
  "TNFα/NFκB signalling",
  "TGF-β signalling",
  "IL6–JAK–STAT3 signalling",
  "Mitotic spindle",
  "UV response down",
  "UV response up",
  "Inflammatory response",
  "Heme metabolism",
  "Estrogen response early",
  "Peroxisome",
  "Myc targets V1",
  "Cholesterol homeostasis",
  "Oxidative phosphorylation",
  "DNA repair",
  "E2F targets"
)

plot_df <- plot_df %>%
  filter(pathway_label %in% pathway_order)

plot_df$pathway_label <- factor(
  plot_df$pathway_label,
  levels = rev(pathway_order)
)

## 3) Create Hallmark dot plot

p <- ggplot(plot_df, aes(x = NES, y = pathway_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.7) +
  geom_point(aes(size = size, colour = neg_log10_padj), alpha = 0.95) +
  scale_colour_gradientn(
    colours = c("#2C7BB6", "#FEE08B", "#D7191C"),
    name = "-log10(padj)"
  ) +
  scale_size_continuous(
    name = "Gene set size",
    range = c(4, 13)
  ) +
  scale_x_continuous(
    limits = c(-2.4, 2.4),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_discrete(
    expand = expansion(mult = c(0.03, 0.02))
  ) +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = "Pathway"
  ) +
  annotate(
    "text",
    x = -1.55,
    y = 16,
    label = "Enriched in NgR3 siRNA",
    fontface = "bold",
    size = 4.8,
    hjust = 0.5
  ) +
  annotate(
    "text",
    x = 1.55,
    y = 16,
    label = "Enriched in 100,000 EPR11334",
    fontface = "bold",
    size = 4.8,
    hjust = 0.5
  ) +
  coord_cartesian(
    ylim = c(1, 15.2),
    clip = "off"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 20, hjust = 0.5, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 13, hjust = 0.5, margin = margin(b = 18)),
    axis.title.x = element_text(face = "bold", size = 16),
    axis.title.y = element_text(face = "bold", size = 16),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    panel.grid.major.y = element_line(colour = "grey85", linewidth = 0.8),
    panel.grid.major.x = element_line(colour = "grey85", linewidth = 0.8),
    panel.grid.minor = element_blank(),
    legend.title = element_text(face = "bold", size = 13),
    legend.text = element_text(size = 11),
    plot.margin = margin(t = 20, r = 30, b = 30, l = 20)
  )

p

## 4) Save Hallmark dot plot

png(
  "GSEA_Hallmark_normalised_dotplot_100000EPR_vs_NgR3siRNA_significant.png",
  width = 3000,
  height = 2400,
  res = 250
)

grid.newpage()

p_grob <- ggplotGrob(p)

pushViewport(
  viewport(
    x = 0.5,
    y = 0.5,
    width = 0.92,
    height = 0.88,
    just = c("center", "center")
  )
)

grid.draw(p_grob)
popViewport()

dev.off()


## C2 GSEA dot plot: 1:100,000 EPR11334 vs NgR3 siRNA interaction model
## Comparison: (1:100,000 EPR11334 vs IgG) vs (NgR3 siRNA vs scrambled siRNA)

## 1) Select significant C2 pathways for plotting

selected_pathways_c2 <- c(
  "CONCANNON_APOPTOSIS_BY_EPOXOMICIN_UP",
  "KEGG_APOPTOSIS",
  "BIOCARTA_IL6_PATHWAY",
  "WP_TGFBETA_RECEPTOR_SIGNALING",
  "ZHOU_INFLAMMATORY_RESPONSE_LIVE_UP",
  "KOHN_EMT_EPITHELIAL",
  "WP_REGULATION_OF_ACTIN_CYTOSKELETON",
  "WP_INTEGRINMEDIATED_CELL_ADHESION",
  "WANG_TUMOR_INVASIVENESS_UP",
  "REACTOME_CELL_CYCLE_CHECKPOINTS",
  "WP_DNA_REPLICATION",
  "DANG_MYC_TARGETS_UP",
  "REN_BOUND_BY_E2F",
  "KEGG_RIBOSOME",
  "REACTOME_CHOLESTEROL_BIOSYNTHESIS",
  "WP_ELECTRON_TRANSPORT_CHAIN_OXPHOS_SYSTEM_IN_MITOCHONDRIA"
)

plot_df <- fgseaRes_interaction_C2_100000_clean %>%
  filter(
    padj < 0.05,
    pathway %in% selected_pathways_c2
  ) %>%
  mutate(
    pathway_label = case_when(
      pathway == "CONCANNON_APOPTOSIS_BY_EPOXOMICIN_UP" ~ "Apoptosis (epoxomicin response)",
      pathway == "KEGG_APOPTOSIS" ~ "Apoptosis",
      pathway == "BIOCARTA_IL6_PATHWAY" ~ "IL6 signalling",
      pathway == "WP_TGFBETA_RECEPTOR_SIGNALING" ~ "TGFβ receptor signalling",
      pathway == "ZHOU_INFLAMMATORY_RESPONSE_LIVE_UP" ~ "Inflammatory response",
      pathway == "KOHN_EMT_EPITHELIAL" ~ "EMT / epithelial transition",
      pathway == "WP_REGULATION_OF_ACTIN_CYTOSKELETON" ~ "Actin cytoskeleton regulation",
      pathway == "WP_INTEGRINMEDIATED_CELL_ADHESION" ~ "Integrin-mediated cell adhesion",
      pathway == "WANG_TUMOR_INVASIVENESS_UP" ~ "Tumor invasiveness",
      pathway == "REACTOME_CELL_CYCLE_CHECKPOINTS" ~ "Cell cycle checkpoints",
      pathway == "WP_DNA_REPLICATION" ~ "DNA replication",
      pathway == "DANG_MYC_TARGETS_UP" ~ "MYC targets",
      pathway == "REN_BOUND_BY_E2F" ~ "E2F target genes",
      pathway == "KEGG_RIBOSOME" ~ "Ribosome",
      pathway == "REACTOME_CHOLESTEROL_BIOSYNTHESIS" ~ "Cholesterol biosynthesis",
      pathway == "WP_ELECTRON_TRANSPORT_CHAIN_OXPHOS_SYSTEM_IN_MITOCHONDRIA" ~ "Mitochondrial ETC / OXPHOS",
      TRUE ~ pathway
    ),
    neg_log10_padj = -log10(padj)
  )

## 2) Set pathway order

pathway_order <- c(
  "Apoptosis (epoxomicin response)",
  "Apoptosis",
  "IL6 signalling",
  "TGFβ receptor signalling",
  "Inflammatory response",
  "EMT / epithelial transition",
  "Actin cytoskeleton regulation",
  "Integrin-mediated cell adhesion",
  "Tumor invasiveness",
  "Cell cycle checkpoints",
  "DNA replication",
  "MYC targets",
  "E2F target genes",
  "Ribosome",
  "Cholesterol biosynthesis",
  "Mitochondrial ETC / OXPHOS"
)

plot_df <- plot_df %>%
  filter(pathway_label %in% pathway_order)

plot_df$pathway_label <- factor(
  plot_df$pathway_label,
  levels = rev(pathway_order)
)

n_pathways <- length(pathway_order)
label_y <- n_pathways + 0.85

## 3) Create C2 dot plot

p <- ggplot(plot_df, aes(x = NES, y = pathway_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.7) +
  geom_point(aes(size = size, colour = neg_log10_padj), alpha = 0.95) +
  scale_colour_gradientn(
    colours = c("#2C7BB6", "#FEE08B", "#D7191C"),
    name = "-log10(padj)"
  ) +
  scale_size_continuous(
    name = "Gene set size",
    range = c(4, 13)
  ) +
  scale_x_continuous(
    limits = c(-3.8, 3.8),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_discrete(
    expand = expansion(mult = c(0.03, 0.02))
  ) +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = "Pathway"
  ) +
  annotate(
    "text",
    x = -2.0,
    y = label_y,
    label = "Enriched in NgR3 siRNA",
    fontface = "bold",
    size = 4.8,
    hjust = 0.5
  ) +
  annotate(
    "text",
    x = 1.9,
    y = label_y,
    label = "Enriched in 100,000 EPR11334",
    fontface = "bold",
    size = 4.8,
    hjust = 0.5
  ) +
  coord_cartesian(
    ylim = c(1, n_pathways + 0.8),
    clip = "off"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 20, hjust = 0.5, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 13, hjust = 0.5, margin = margin(b = 18)),
    axis.title.x = element_text(face = "bold", size = 16),
    axis.title.y = element_text(face = "bold", size = 16),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 12),
    panel.grid.major.y = element_line(colour = "grey85", linewidth = 0.8),
    panel.grid.major.x = element_line(colour = "grey85", linewidth = 0.8),
    panel.grid.minor = element_blank(),
    legend.title = element_text(face = "bold", size = 13),
    legend.text = element_text(size = 11),
    plot.margin = margin(t = 35, r = 40, b = 35, l = 30)
  )

p

## 4) Save C2 dot plot

png(
  "C2_normalised_dotplot_100000EPR_vs_NgR3siRNA.png",
  width = 3200,
  height = 2500,
  res = 250
)

grid.newpage()

p_grob <- ggplotGrob(p)

pushViewport(
  viewport(
    x = 0.47,
    y = 0.5,
    width = 0.92,
    height = 0.90,
    just = c("center", "center")
  )
)

grid.draw(p_grob)
popViewport()

dev.off()