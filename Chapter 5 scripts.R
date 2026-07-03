## PCA for antibody-only samples
## IgG control vs 1:1,000 EPR11334 vs 1:100,000 EPR11334
## This section uses the existing DESeq2 object `dds`.
## Run this section after the main DESeq2 object has been created.

## 1) Subset DESeq2 object to antibody-only samples

dds_subset <- dds[, dds$Condition %in% c("IGG", "1:1000_EPR", "1:100000_EPR")]

dds_subset$Condition <- droplevels(dds_subset$Condition)

## 2) Check subset

colData(dds_subset)$Condition
table(dds_subset$Condition)

## 3) Run variance stabilising transformation on subset

vsd_subset <- vst(dds_subset, blind = TRUE)

## 4) Generate PCA data

pca_data <- plotPCA(
  vsd_subset,
  intgroup = "Condition",
  returnData = TRUE
)

## 5) Check PCA groups

unique(as.character(pca_data$Condition))
table(pca_data$Condition, useNA = "ifany")

## 6) Clean condition labels

pca_data$Condition <- as.character(pca_data$Condition)

pca_data$Condition[pca_data$Condition == "IGG"] <- "IgG control"
pca_data$Condition[pca_data$Condition == "1:1000_EPR"] <- "1:1,000 EPR11334"
pca_data$Condition[pca_data$Condition == "1:100000_EPR"] <- "1:100,000 EPR11334"

## 7) Set plotting order

pca_data$Condition <- factor(
  pca_data$Condition,
  levels = c(
    "IgG control",
    "1:1,000 EPR11334",
    "1:100,000 EPR11334"
  )
)

## 8) Check final labels

table(pca_data$Condition, useNA = "ifany")

## 9) Extract percentage variance explained

percentVar <- round(100 * attr(pca_data, "percentVar"))

## 10) Plot PCA

p <- ggplot(pca_data, aes(PC1, PC2, color = Condition)) +
  geom_point(size = 6) +
  xlab(paste0("PC1 (", percentVar[1], "%)")) +
  ylab(paste0("PC2 (", percentVar[2], "%)")) +
  theme_bw(base_size = 18) +
  coord_fixed()

print(p)

## 11) Save PCA plot

ggsave(
  filename = "/home/yougr345/PCA_EPR_only.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)



## Hallmark GSEA: IgG vs 1:1,000 EPR11334
## Comparison: IgG control vs 1:1,000 EPR11334
## This section uses the existing metadata object `meta` and count matrix `count_matrix`.
## Run this section after metadata cleaning and count matrix creation.

## 1) Create comparison factor

meta$Condition <- sub("_(\\d+)$", "", meta$Group)
meta$comparison_group <- NA_character_

meta$comparison_group[
  grepl("IGG", meta$Condition, ignore.case = TRUE)
] <- "IGG"

meta$comparison_group[
  grepl("EPR", meta$Condition, ignore.case = TRUE) &
    grepl("1000", meta$Condition) &
    !grepl("100000", meta$Condition)
] <- "EPR_1to1000"

meta$comparison_group <- factor(
  meta$comparison_group,
  levels = c("EPR_1to1000", "IGG")
)

## 2) Keep samples for this comparison

keep <- !is.na(meta$comparison_group)

meta_sub <- meta[keep, , drop = FALSE]
count_sub <- count_matrix[, rownames(meta_sub), drop = FALSE]

print(table(meta_sub$comparison_group))

## 3) Build DESeq2 dataset

dds_IGG_vs_1000EPR <- DESeqDataSetFromMatrix(
  countData = count_sub,
  colData   = meta_sub,
  design    = ~ comparison_group
)

dds_IGG_vs_1000EPR <- dds_IGG_vs_1000EPR[
  rowSums(counts(dds_IGG_vs_1000EPR)) > 10,
]

dds_IGG_vs_1000EPR <- DESeq(dds_IGG_vs_1000EPR)

resultsNames(dds_IGG_vs_1000EPR)

## 4) Get DESeq2 results

res_IGG_vs_1000EPR <- results(
  dds_IGG_vs_1000EPR,
  contrast = c("comparison_group", "IGG", "EPR_1to1000")
)

## 5) Create ranked gene list

gene_list_IGG_vs_1000EPR <- res_IGG_vs_1000EPR$log2FoldChange
names(gene_list_IGG_vs_1000EPR) <- rownames(res_IGG_vs_1000EPR)

gene_symbols_IGG_vs_1000EPR <- mapIds(
  org.Hs.eg.db,
  keys = names(gene_list_IGG_vs_1000EPR),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

gene_df_IGG_vs_1000EPR <- data.frame(
  symbol = gene_symbols_IGG_vs_1000EPR,
  log2FoldChange = gene_list_IGG_vs_1000EPR,
  stringsAsFactors = FALSE
)

gene_df_IGG_vs_1000EPR <- gene_df_IGG_vs_1000EPR[
  !is.na(gene_df_IGG_vs_1000EPR$symbol),
]

gene_df_IGG_vs_1000EPR <- gene_df_IGG_vs_1000EPR %>%
  group_by(symbol) %>%
  slice_max(order_by = abs(log2FoldChange), n = 1, with_ties = FALSE) %>%
  ungroup()

gene_list_IGG_vs_1000EPR <- gene_df_IGG_vs_1000EPR$log2FoldChange
names(gene_list_IGG_vs_1000EPR) <- gene_df_IGG_vs_1000EPR$symbol

gene_list_IGG_vs_1000EPR <- sort(
  gene_list_IGG_vs_1000EPR,
  decreasing = TRUE
)

## 6) Load Hallmark gene sets

msig_hallmark <- msigdbr(
  species = "Homo sapiens",
  category = "H"
)

pathways_hallmark <- msig_hallmark %>%
  split(x = .$gene_symbol, f = .$gs_name)

## 7) Run Hallmark fgsea

fgseaRes_IGG_vs_1000EPR <- fgsea(
  pathways = pathways_hallmark,
  stats = gene_list_IGG_vs_1000EPR,
  minSize = 15,
  maxSize = 500
)

fgseaRes_IGG_vs_1000EPR <- fgseaRes_IGG_vs_1000EPR %>%
  arrange(padj)

fgseaRes_IGG_vs_1000EPR[1:20, ]

## 8) Export Hallmark GSEA results

fgsea_export_IGG_vs_1000EPR <- fgseaRes_IGG_vs_1000EPR %>%
  mutate(
    leadingEdge = sapply(leadingEdge, paste, collapse = ", ")
  )

write.csv(
  fgsea_export_IGG_vs_1000EPR,
  "gsea_results_IGG_vs_1000EPR.csv",
  row.names = FALSE
)



## Hallmark GSEA heatmap: IgG vs 1:1,000 EPR11334
## Leading-edge genes from E2F targets and G2M checkpoint
## Pathways upregulated in 1:1,000 EPR11334 compared to IgG
## This section uses the existing DESeq2 object `dds_full`
## and the existing fgsea object `fgseaRes_IGG_vs_1000`.

## 1) Run VST from the full DESeq2 object

vsd <- vst(dds_full, blind = FALSE)
expr_mat <- assay(vsd)

## 2) Keep only IgG and 1:1000_EPR samples

samples_keep <- colData(dds_full)$Condition %in% c("IGG", "1:1000_EPR")

expr_sub <- expr_mat[, samples_keep, drop = FALSE]
meta_sub <- as.data.frame(colData(dds_full)[samples_keep, ])

table(meta_sub$Condition)

## 3) Select pathways of interest from this comparison

pathways_of_interest <- c(
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT"
)

fgsea_top <- fgseaRes_IGG_vs_1000 %>%
  as.data.frame() %>%
  filter(pathway %in% pathways_of_interest)

## 4) Extract leading-edge genes

leading_edge_list <- setNames(fgsea_top$leadingEdge, fgsea_top$pathway)

top_n <- 15
leading_edge_top_list <- lapply(leading_edge_list, function(x) head(x, top_n))
heatmap_genes <- unique(unlist(leading_edge_top_list))

## 5) Build pathway membership table

gene_pathway_df <- bind_rows(lapply(names(leading_edge_top_list), function(pw) {
  data.frame(
    Gene = leading_edge_top_list[[pw]],
    Pathway = pw,
    stringsAsFactors = FALSE
  )
}))

gene_pathway_df <- gene_pathway_df[!duplicated(gene_pathway_df$Gene), ]

gene_pathway_df$Pathway <- recode(
  gene_pathway_df$Pathway,
  "HALLMARK_E2F_TARGETS" = "E2F targets",
  "HALLMARK_G2M_CHECKPOINT" = "G2M checkpoint"
)

## 6) Map Ensembl IDs in the expression matrix to gene symbols

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

## If multiple Ensembl IDs map to the same gene symbol,
## keep the one with the highest mean expression

expr_df_annot$mean_expr <- rowMeans(expr_df_annot[, colnames(expr_sub), drop = FALSE])

expr_df_annot <- expr_df_annot %>%
  group_by(SYMBOL) %>%
  slice_max(order_by = mean_expr, n = 1, with_ties = FALSE) %>%
  ungroup()

expr_symbol_mat <- as.matrix(expr_df_annot[, colnames(expr_sub), drop = FALSE])
rownames(expr_symbol_mat) <- expr_df_annot$SYMBOL

## 7) Keep only leading-edge genes present after mapping

heatmap_genes_present <- intersect(gene_pathway_df$Gene, rownames(expr_symbol_mat))

expr_heatmap <- expr_symbol_mat[heatmap_genes_present, , drop = FALSE]

gene_pathway_df <- gene_pathway_df %>%
  filter(Gene %in% heatmap_genes_present)

## 8) Scale expression by gene

expr_heatmap_scaled <- t(scale(t(expr_heatmap)))
expr_heatmap_scaled <- expr_heatmap_scaled[complete.cases(expr_heatmap_scaled), ]

gene_pathway_df <- gene_pathway_df %>%
  filter(Gene %in% rownames(expr_heatmap_scaled))

## 9) Order samples: EPR first, IgG second

meta_sub$Replicate <- as.character(meta_sub$Replicate)

order_idx <- order(meta_sub$Condition == "IGG", as.numeric(meta_sub$Replicate))

expr_heatmap_scaled <- expr_heatmap_scaled[, order_idx, drop = FALSE]
meta_sub <- meta_sub[order_idx, , drop = FALSE]

## 10) Create cleaner sample labels

display_condition <- ifelse(
  meta_sub$Condition == "1:1000_EPR",
  "1:1,000 EPR11334",
  ifelse(
    meta_sub$Condition == "IGG",
    "IgG control",
    as.character(meta_sub$Condition)
  )
)

colnames(expr_heatmap_scaled) <- paste(display_condition, meta_sub$Replicate)

## 11) Cluster genes within each pathway block

pathway_levels <- c("E2F targets", "G2M checkpoint")

ordered_genes <- unlist(lapply(pathway_levels, function(pw) {
  genes_pw <- gene_pathway_df$Gene[gene_pathway_df$Pathway == pw]
  genes_pw <- intersect(genes_pw, rownames(expr_heatmap_scaled))
  
  if (length(genes_pw) > 1) {
    hc <- hclust(dist(expr_heatmap_scaled[genes_pw, , drop = FALSE]))
    genes_pw <- genes_pw[hc$order]
  }
  
  genes_pw
}))

expr_heatmap_scaled <- expr_heatmap_scaled[ordered_genes, , drop = FALSE]

gene_pathway_df <- gene_pathway_df %>%
  filter(Gene %in% ordered_genes) %>%
  slice(match(ordered_genes, Gene))

## 12) Create annotations

col_annot <- data.frame(
  Condition = factor(
    display_condition,
    levels = c("1:1,000 EPR11334", "IgG control")
  )
)

rownames(col_annot) <- colnames(expr_heatmap_scaled)

row_annot <- data.frame(
  Pathway = factor(
    gene_pathway_df$Pathway,
    levels = pathway_levels
  )
)

rownames(row_annot) <- gene_pathway_df$Gene

annotation_colors <- list(
  Condition = c(
    "1:1,000 EPR11334" = "black",
    "IgG control" = "grey70"
  ),
  Pathway = c(
    "E2F targets" = "#7B2CBF",
    "G2M checkpoint" = "#B8A1E3"
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
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize_col = 10,
  fontsize_row = 8,
  scale = "none",
  breaks = breaks,
  border_color = "grey60",
  color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
  silent = TRUE
)

## 15) Save heatmap

png(
  "1000EPR_vs_IGG_heatmap.png",
  width = 2800,
  height = 1600,
  res = 300
)

grid.newpage()
pushViewport(viewport(x = 0.52, y = 0.5, width = 0.96, height = 0.90))
grid.draw(p$gtable)
upViewport()

dev.off()



## C2 GSEA: IgG vs 1:100,000 EPR11334
## This section uses the existing metadata object `meta` and count matrix `count_matrix`.

## 1) Create comparison factor

meta$Condition <- sub("_(\\d+)$", "", meta$Group)
meta$comparison_group <- NA_character_

meta$comparison_group[
  grepl("IGG", meta$Condition, ignore.case = TRUE)
] <- "IGG"

meta$comparison_group[
  grepl("EPR", meta$Condition, ignore.case = TRUE) &
    grepl("100000", meta$Condition)
] <- "EPR_1to100000"

meta$comparison_group <- factor(
  meta$comparison_group,
  levels = c("EPR_1to100000", "IGG")
)

## 2) Keep samples for this comparison

keep <- !is.na(meta$comparison_group)

meta_sub <- meta[keep, , drop = FALSE]
count_sub <- count_matrix[, rownames(meta_sub), drop = FALSE]

print(table(meta_sub$comparison_group))

## 3) Build DESeq2 dataset

dds_IGG_vs_100000EPR <- DESeqDataSetFromMatrix(
  countData = count_sub,
  colData   = meta_sub,
  design    = ~ comparison_group
)

dds_IGG_vs_100000EPR <- dds_IGG_vs_100000EPR[
  rowSums(counts(dds_IGG_vs_100000EPR)) > 10,
]

dds_IGG_vs_100000EPR <- DESeq(dds_IGG_vs_100000EPR)

resultsNames(dds_IGG_vs_100000EPR)

## 4) Get DESeq2 results

res_IGG_vs_100000EPR <- results(
  dds_IGG_vs_100000EPR,
  contrast = c("comparison_group", "IGG", "EPR_1to100000")
)

## 5) Create ranked gene list

gene_list_IGG_vs_100000EPR <- res_IGG_vs_100000EPR$log2FoldChange
names(gene_list_IGG_vs_100000EPR) <- rownames(res_IGG_vs_100000EPR)

gene_symbols_IGG_vs_100000EPR <- mapIds(
  org.Hs.eg.db,
  keys = names(gene_list_IGG_vs_100000EPR),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

gene_df_IGG_vs_100000EPR <- data.frame(
  symbol = gene_symbols_IGG_vs_100000EPR,
  log2FoldChange = gene_list_IGG_vs_100000EPR,
  stringsAsFactors = FALSE
)

gene_df_IGG_vs_100000EPR <- gene_df_IGG_vs_100000EPR[
  !is.na(gene_df_IGG_vs_100000EPR$symbol),
]

gene_df_IGG_vs_100000EPR <- gene_df_IGG_vs_100000EPR %>%
  group_by(symbol) %>%
  slice_max(order_by = abs(log2FoldChange), n = 1, with_ties = FALSE) %>%
  ungroup()

gene_list_IGG_vs_100000EPR <- gene_df_IGG_vs_100000EPR$log2FoldChange
names(gene_list_IGG_vs_100000EPR) <- gene_df_IGG_vs_100000EPR$symbol

gene_list_IGG_vs_100000EPR <- sort(
  gene_list_IGG_vs_100000EPR,
  decreasing = TRUE
)

## 6) Load C2 gene sets

msig_c2 <- msigdbr(
  species = "Homo sapiens",
  collection = "C2"
)

pathways_c2 <- msig_c2 %>%
  split(x = .$gene_symbol, f = .$gs_name)

## 7) Run C2 fgsea

fgseaRes_c2_IGG_vs_100000EPR <- fgsea(
  pathways = pathways_c2,
  stats = gene_list_IGG_vs_100000EPR,
  minSize = 15,
  maxSize = 500
)

fgseaRes_c2_IGG_vs_100000EPR <- fgseaRes_c2_IGG_vs_100000EPR %>%
  arrange(padj)

fgseaRes_c2_IGG_vs_100000EPR[1:20, ]

## 8) Save C2 GSEA results

saveRDS(
  fgseaRes_c2_IGG_vs_100000EPR,
  "~/fgseaRes_c2_IGG_vs_100000EPR.rds"
)

fgsea_export_c2_IGG_vs_100000EPR <- fgseaRes_c2_IGG_vs_100000EPR %>%
  mutate(
    leadingEdge = sapply(leadingEdge, paste, collapse = ", ")
  )

saveRDS(
  fgsea_export_c2_IGG_vs_100000EPR,
  file = "~/fgsea_export_c2_IGG_vs_100000EPR.rds"
)

write.csv(
  fgsea_export_c2_IGG_vs_100000EPR,
  "gsea_results_C2_IGG_vs_100000EPR.csv",
  row.names = FALSE
)


## C2 GSEA dot plot: IgG vs 1:100,000 EPR11334
## Pathways decreased in 1:100,000 EPR11334 compared to IgG

## 1) Create plotting data

plot_df <- data.frame(
  pathway = c(
    "Ribosome",
    "Translation initiation",
    "Translation elongation",
    "SRP-dependent protein targeting",
    "Aerobic respiration / electron transport",
    "Oxidative phosphorylation",
    "GCN2 response to amino acid deficiency",
    "Cellular response to starvation",
    "ROBO signalling / migration",
    "SLIT/ROBO regulation"
  ),
  pathway_full = c(
    "KEGG_RIBOSOME",
    "REACTOME_EUKARYOTIC_TRANSLATION_INITIATION",
    "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION",
    "REACTOME_SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE",
    "REACTOME_AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT",
    "WP_ELECTRON_TRANSPORT_CHAIN_OXPHOS_SYSTEM_IN_MITOCHONDRIA",
    "REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY",
    "REACTOME_CELLULAR_RESPONSE_TO_STARVATION",
    "REACTOME_SIGNALING_BY_ROBO_RECEPTORS",
    "REACTOME_REGULATION_OF_EXPRESSION_OF_SLITS_AND_ROBOS"
  ),
  NES = c(
    2.078791823,
    1.991959054,
    2.076720178,
    2.050581578,
    1.659079928,
    1.829155223,
    1.886710997,
    1.676144325,
    1.702089048,
    1.935380958
  ),
  padj = c(
    0.000292856,
    0.000533318,
    0.000296780,
    0.000292856,
    0.012776841,
    0.024886433,
    0.009216133,
    0.040263763,
    0.012776841,
    0.000521661
  ),
  size = c(
    85,
    117,
    90,
    110,
    238,
    90,
    99,
    154,
    192,
    152
  )
)

## 2) Flip NES values for plotting direction

plot_df <- plot_df %>%
  mutate(
    NES_plot = -NES,
    neglog10_padj = -log10(padj)
  )

## 3) Set pathway order

plot_df$pathway <- factor(
  plot_df$pathway,
  levels = rev(c(
    "Ribosome",
    "Translation initiation",
    "Translation elongation",
    "SRP-dependent protein targeting",
    "Aerobic respiration / electron transport",
    "Oxidative phosphorylation",
    "GCN2 response to amino acid deficiency",
    "Cellular response to starvation",
    "ROBO signalling / migration",
    "SLIT/ROBO regulation"
  ))
)

## 4) Create dot plot

p <- ggplot(plot_df, aes(x = NES_plot, y = pathway)) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.7,
    colour = "grey35"
  ) +
  geom_point(
    aes(size = size, colour = neglog10_padj),
    alpha = 0.95
  ) +
  annotate(
    "text",
    x = -1.35,
    y = length(levels(plot_df$pathway)) + 0.55,
    label = "Decreased in 1:100,000 EPR11334",
    size = 6.2,
    fontface = "bold"
  ) +
  scale_size_continuous(
    name = "Gene set size",
    range = c(7, 24),
    breaks = c(100, 150, 200)
  ) +
  scale_colour_gradientn(
    name = expression(-log[10](padj)),
    colours = c("#2C7BB6", "#F6C85F", "#E85D3F", "#E41A1C")
  ) +
  scale_x_continuous(
    limits = c(-2.4, 0.2),
    breaks = c(-2, -1, 0),
    labels = c("-2", "-1", "0"),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = "Pathway"
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 13.5, hjust = 0.5, margin = margin(b = 24)),
    axis.title.x = element_text(size = 19, face = "bold"),
    axis.title.y = element_text(size = 19, face = "bold"),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 16),
    panel.grid.major.y = element_line(colour = "grey80", linewidth = 0.8),
    panel.grid.major.x = element_line(colour = "grey80", linewidth = 0.8),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14),
    legend.spacing.y = unit(0.4, "cm"),
    plot.margin = margin(t = 55, r = 55, b = 35, l = 140)
  ) +
  guides(
    size = guide_legend(order = 1, override.aes = list(alpha = 1)),
    colour = guide_colorbar(order = 2)
  )

print(p)

## 5) Save dot plot

ggsave(
  filename = "C2_Dot_plot_IgG_vs_100000_EPR.png",
  plot = p,
  width = 18,
  height = 10,
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)


## C2 heatmap: IgG vs 1:100,000 EPR11334
## This section uses the existing DESeq2 object `dds_IGG_vs_100000EPR`
## and the existing gene mapping object `gene_map`.

## 1) Create transformed object for this comparison

vsd_IGG_vs_100000EPR <- vst(dds_IGG_vs_100000EPR, blind = TRUE)

expr_mat_all <- assay(vsd_IGG_vs_100000EPR)
sample_info <- as.data.frame(colData(vsd_IGG_vs_100000EPR))

## 2) Identify sample columns

igg_cols <- rownames(sample_info)[
  grepl("^IGG$", sample_info$Condition, ignore.case = TRUE)
]

epr_cols <- rownames(sample_info)[
  grepl("100000|100,000|1:100000|1:100,000", sample_info$Condition, ignore.case = TRUE)
]

if (length(igg_cols) == 0) {
  stop("No IgG samples found in sample_info$Condition.")
}

if (length(epr_cols) == 0) {
  stop("No 1:100,000 EPR11334 samples found in sample_info$Condition.")
}

expr_mat <- expr_mat_all[, c(igg_cols, epr_cols), drop = FALSE]

## 3) Map Ensembl IDs to gene symbols using gene_map

gm <- as.data.frame(gene_map)

possible_ens_cols <- c(
  "ensembl_gene_id",
  "ensembl_id",
  "ensembl",
  "ENSEMBL",
  "gene_id",
  "Geneid"
)

possible_symbol_cols <- c(
  "gene_symbol",
  "symbol",
  "SYMBOL",
  "external_gene_name",
  "gene_name",
  "hgnc_symbol",
  "Gene",
  "gene"
)

ens_col <- possible_ens_cols[possible_ens_cols %in% colnames(gm)][1]
sym_col <- possible_symbol_cols[possible_symbol_cols %in% colnames(gm)][1]

if (is.na(ens_col) || length(ens_col) == 0) {
  stop("Could not find an Ensembl column in gene_map. Run colnames(gene_map) to inspect it.")
}

if (is.na(sym_col) || length(sym_col) == 0) {
  stop("Could not find a gene symbol column in gene_map. Run colnames(gene_map) to inspect it.")
}

map_df <- gm[, c(ens_col, sym_col)]

colnames(map_df) <- c("ensembl", "symbol")

map_df$ensembl <- sub("\\..*$", "", as.character(map_df$ensembl))
map_df$symbol <- as.character(map_df$symbol)

expr_ids <- sub("\\..*$", "", rownames(expr_mat))

expr_df <- data.frame(
  ensembl = expr_ids,
  expr_mat,
  check.names = FALSE
)

expr_df <- expr_df %>%
  left_join(map_df, by = "ensembl") %>%
  filter(!is.na(symbol), symbol != "")

expr_symbol_mat <- as.matrix(expr_df[, colnames(expr_mat), drop = FALSE])
rownames(expr_symbol_mat) <- expr_df$symbol

## If duplicated symbols exist, keep the row with the highest mean expression

expr_symbol_mat <- expr_symbol_mat[
  order(rowMeans(expr_symbol_mat), decreasing = TRUE),
  ,
  drop = FALSE
]

expr_symbol_mat <- expr_symbol_mat[
  !duplicated(rownames(expr_symbol_mat)),
  ,
  drop = FALSE
]

## 4) Define heatmap genes

genes_heatmap <- c(
  # Translation
  "EIF4A1", "EIF2S2", "EIF3K", "EEF1B2", "EEF1E1", "RPS27", "RPS29",
  
  # Mitochondrial respiration
  "COX17", "COX6A1", "COX5B", "NDUFA1", "NDUFB3", "UQCR10", "ATP5MF",
  
  # Stress response
  "CASTOR1", "LAMTOR2", "LAMTOR3", "LAMTOR4", "LAMTOR5", "RRAGB",
  
  # ROBO / motility
  "NTN1", "NELL2", "FLRT3", "COL4A5", "GPC1", "PFN1", "EVL", "PAK1"
)

## 5) Check genes present after mapping

genes_present <- genes_heatmap[
  genes_heatmap %in% rownames(expr_symbol_mat)
]

genes_missing <- setdiff(
  genes_heatmap,
  rownames(expr_symbol_mat)
)

message("Genes found: ", length(genes_present), " / ", length(genes_heatmap))

if (length(genes_missing) > 0) {
  message(
    "These genes were not found and will be omitted: ",
    paste(genes_missing, collapse = ", ")
  )
}

if (length(genes_present) < 20) {
  stop("Too many genes are missing after mapping. Check colnames(gene_map).")
}

expr_mat2 <- expr_symbol_mat[genes_present, , drop = FALSE]

expr_mat2 <- expr_mat2[
  match(genes_present, rownames(expr_mat2)),
  ,
  drop = FALSE
]

## 6) Scale expression by gene

expr_mat2_scaled <- t(scale(t(expr_mat2)))
expr_mat2_scaled[is.na(expr_mat2_scaled)] <- 0

## 7) Create column annotation

annotation_col <- data.frame(
  Condition = c(
    rep("IgG control", length(igg_cols)),
    rep("1:100,000 EPR11334", length(epr_cols))
  )
)

rownames(annotation_col) <- colnames(expr_mat2_scaled)

## 8) Create row annotation

annotation_row <- data.frame(
  Pathway = c(
    rep("Translation initiation", 7),
    rep("Mitochondrial respiration", 7),
    rep("Stress response", 6),
    rep("ROBO signalling / migration", 8)
  ),
  `Functional group` = c(
    rep("Translation", 7),
    rep("Metabolism", 7),
    rep("Stress response", 6),
    rep("Motility", 8)
  ),
  check.names = FALSE
)

rownames(annotation_row) <- rownames(expr_mat2_scaled)

## 9) Set annotation colours

annotation_colors <- list(
  Condition = c(
    "IgG control" = "black",
    "1:100,000 EPR11334" = "grey70"
  ),
  Pathway = c(
    "Translation initiation" = "#B7DFAE",
    "Mitochondrial respiration" = "#8FC97A",
    "Stress response" = "#5FA85C",
    "ROBO signalling / migration" = "#3F7F3F"
  ),
  `Functional group` = c(
    "Translation" = "#D7A9FF",
    "Metabolism" = "#BE7DF2",
    "Stress response" = "#9D4EDD",
    "Motility" = "#6F2DBD"
  )
)

## 10) Draw heatmap

ph <- pheatmap(
  expr_mat2_scaled,
  annotation_col = annotation_col,
  annotation_row = annotation_row,
  annotation_colors = annotation_colors,
  show_rownames = TRUE,
  labels_row = rownames(expr_mat2_scaled),
  show_colnames = FALSE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_row = 12,
  fontsize = 11,
  border_color = "grey65",
  color = colorRampPalette(c("#4575B4", "#F3E5AB", "#E85D3F"))(100),
  breaks = seq(-1.8, 1.8, length.out = 101),
  silent = TRUE
)

## 11) Preview heatmap

grid.newpage()

pushViewport(
  viewport(
    x = 0.54,
    y = 0.5,
    width = 0.88,
    height = 0.90
  )
)

grid.draw(ph$gtable)
popViewport()

## 12) Save heatmap

png(
  "Heatmap_IGG_vs_100000_EPR_C2_pathways.png",
  width = 2800,
  height = 3400,
  res = 220
)

grid.newpage()

pushViewport(
  viewport(
    x = 0.54,
    y = 0.5,
    width = 0.88,
    height = 0.90
  )
)

grid.draw(ph$gtable)
popViewport()

dev.off()



## C2 analysis: 1:1,000 EPR11334 vs 1:100,000 EPR11334
## This section uses the existing metadata object `meta` and count matrix `count_matrix`.

## 1) Create comparison factor

meta$Condition <- sub("_(\\d+)$", "", meta$Group)
meta$comparison_group <- NA_character_

meta$comparison_group[
  grepl("EPR", meta$Condition, ignore.case = TRUE) &
    grepl("1000", meta$Condition) &
    !grepl("100000", meta$Condition)
] <- "EPR_1to1000"

meta$comparison_group[
  grepl("EPR", meta$Condition, ignore.case = TRUE) &
    grepl("100000", meta$Condition)
] <- "EPR_1to100000"

meta$comparison_group <- factor(
  meta$comparison_group,
  levels = c("EPR_1to100000", "EPR_1to1000")
)

## 2) Keep samples for this comparison

keep <- !is.na(meta$comparison_group)

meta_sub <- meta[keep, , drop = FALSE]
count_sub <- count_matrix[, rownames(meta_sub), drop = FALSE]

print(table(meta_sub$comparison_group))

## 3) Build DESeq2 dataset

dds_1000EPR_vs_100000EPR <- DESeqDataSetFromMatrix(
  countData = count_sub,
  colData   = meta_sub,
  design    = ~ comparison_group
)

dds_1000EPR_vs_100000EPR <- dds_1000EPR_vs_100000EPR[
  rowSums(counts(dds_1000EPR_vs_100000EPR)) > 10,
]

dds_1000EPR_vs_100000EPR <- DESeq(dds_1000EPR_vs_100000EPR)

resultsNames(dds_1000EPR_vs_100000EPR)

## 4) Get DESeq2 results

res_1000EPR_vs_100000EPR <- results(
  dds_1000EPR_vs_100000EPR,
  contrast = c("comparison_group", "EPR_1to1000", "EPR_1to100000")
)

## 5) Create ranked gene list

gene_list_1000EPR_vs_100000EPR <- res_1000EPR_vs_100000EPR$log2FoldChange
names(gene_list_1000EPR_vs_100000EPR) <- rownames(res_1000EPR_vs_100000EPR)

gene_symbols_1000EPR_vs_100000EPR <- mapIds(
  org.Hs.eg.db,
  keys = names(gene_list_1000EPR_vs_100000EPR),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

gene_df_1000EPR_vs_100000EPR <- data.frame(
  symbol = gene_symbols_1000EPR_vs_100000EPR,
  log2FoldChange = gene_list_1000EPR_vs_100000EPR,
  stringsAsFactors = FALSE
)

gene_df_1000EPR_vs_100000EPR <- gene_df_1000EPR_vs_100000EPR[
  !is.na(gene_df_1000EPR_vs_100000EPR$symbol),
]

gene_df_1000EPR_vs_100000EPR <- gene_df_1000EPR_vs_100000EPR %>%
  group_by(symbol) %>%
  slice_max(order_by = abs(log2FoldChange), n = 1, with_ties = FALSE) %>%
  ungroup()

gene_list_1000EPR_vs_100000EPR <- gene_df_1000EPR_vs_100000EPR$log2FoldChange
names(gene_list_1000EPR_vs_100000EPR) <- gene_df_1000EPR_vs_100000EPR$symbol

gene_list_1000EPR_vs_100000EPR <- sort(
  gene_list_1000EPR_vs_100000EPR,
  decreasing = TRUE
)

## 6) Load C2 gene sets

msig_c2 <- msigdbr(
  species = "Homo sapiens",
  category = "C2"
)

pathways_c2 <- msig_c2 %>%
  split(x = .$gene_symbol, f = .$gs_name)

## 7) Run C2 fgsea

fgseaRes_c2_1000EPR_vs_100000EPR <- fgsea(
  pathways = pathways_c2,
  stats = gene_list_1000EPR_vs_100000EPR,
  minSize = 15,
  maxSize = 500
)

fgseaRes_c2_1000EPR_vs_100000EPR <- fgseaRes_c2_1000EPR_vs_100000EPR %>%
  arrange(padj)

fgseaRes_c2_1000EPR_vs_100000EPR[1:20, ]

## 8) Save C2 GSEA results

saveRDS(
  fgseaRes_c2_1000EPR_vs_100000EPR,
  file = "~/fgseaRes_c2_1000EPR_vs_100000EPR.rds"
)

fgsea_export_c2_1000EPR_vs_100000EPR <- fgseaRes_c2_1000EPR_vs_100000EPR %>%
  mutate(
    leadingEdge = sapply(leadingEdge, paste, collapse = ", ")
  )

saveRDS(
  fgsea_export_c2_1000EPR_vs_100000EPR,
  file = "~/fgsea_export_c2_1000EPR_vs_100000EPR.rds"
)

write.csv(
  fgsea_export_c2_1000EPR_vs_100000EPR,
  "gsea_results_C2_1000EPR_vs_100000EPR.csv",
  row.names = FALSE
)


## C2 dot plot: 1:1,000 EPR11334 vs 1:100,000 EPR11334

## 1) Start with C2 GSEA results

dot_df <- fgseaRes_c2_1000EPR_vs_100000EPR

## 2) Select pathways for dot plot

selected_pathways <- c(
  "ALONSO_METASTASIS_UP",
  "WINNEPENNINCKX_MELANOMA_METASTASIS_UP",
  "ZHOU_CELL_CYCLE_GENES_IN_IR_RESPONSE_24HR",
  "REACTOME_SYNTHESIS_OF_DNA",
  "REACTOME_CELLULAR_RESPONSE_TO_STARVATION",
  "REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY",
  "REACTOME_NONSENSE_MEDIATED_DECAY_NMD",
  "REACTOME_SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE",
  "WP_OXIDATIVE_PHOSPHORYLATION",
  "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",
  "KEGG_OXIDATIVE_PHOSPHORYLATION",
  "REACTOME_AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT",
  "REACTOME_TRANSLATION",
  "REACTOME_EUKARYOTIC_TRANSLATION_INITIATION",
  "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION",
  "KEGG_RIBOSOME"
)

## 3) Create pathway labels

pathway_labels <- c(
  "ALONSO_METASTASIS_UP" = "ALONSO metastasis up",
  "WINNEPENNINCKX_MELANOMA_METASTASIS_UP" = "Melanoma metastasis up",
  "ZHOU_CELL_CYCLE_GENES_IN_IR_RESPONSE_24HR" = "Cell cycle genes (IR, 24 h)",
  "REACTOME_SYNTHESIS_OF_DNA" = "Synthesis of DNA",
  "REACTOME_CELLULAR_RESPONSE_TO_STARVATION" = "Cellular response to starvation",
  "REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY" = "GCN2 response to amino acid deficiency",
  "REACTOME_NONSENSE_MEDIATED_DECAY_NMD" = "Nonsense-mediated decay (NMD)",
  "REACTOME_SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE" = "SRP-dependent protein targeting",
  "WP_OXIDATIVE_PHOSPHORYLATION" = "Oxidative phosphorylation (WP)",
  "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT" = "Respiratory electron transport",
  "KEGG_OXIDATIVE_PHOSPHORYLATION" = "Oxidative phosphorylation (KEGG)",
  "REACTOME_AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT" = "Aerobic respiration / electron transport",
  "REACTOME_TRANSLATION" = "Translation",
  "REACTOME_EUKARYOTIC_TRANSLATION_INITIATION" = "Eukaryotic translation initiation",
  "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION" = "Eukaryotic translation elongation",
  "KEGG_RIBOSOME" = "Ribosome"
)

## 4) Add pathway classes

pathway_class <- c(
  "ALONSO_METASTASIS_UP" = "Migration / invasion",
  "WINNEPENNINCKX_MELANOMA_METASTASIS_UP" = "Migration / invasion",
  "ZHOU_CELL_CYCLE_GENES_IN_IR_RESPONSE_24HR" = "Proliferation",
  "REACTOME_SYNTHESIS_OF_DNA" = "Proliferation",
  "REACTOME_CELLULAR_RESPONSE_TO_STARVATION" = "Stress adaptation",
  "REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY" = "Stress adaptation",
  "REACTOME_NONSENSE_MEDIATED_DECAY_NMD" = "Stress adaptation",
  "REACTOME_SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE" = "Stress adaptation",
  "WP_OXIDATIVE_PHOSPHORYLATION" = "Mitochondrial metabolism",
  "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT" = "Mitochondrial metabolism",
  "KEGG_OXIDATIVE_PHOSPHORYLATION" = "Mitochondrial metabolism",
  "REACTOME_AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT" = "Mitochondrial metabolism",
  "REACTOME_TRANSLATION" = "Translation",
  "REACTOME_EUKARYOTIC_TRANSLATION_INITIATION" = "Translation",
  "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION" = "Translation",
  "KEGG_RIBOSOME" = "Translation"
)

## 5) Create plotting data

plot_df <- dot_df %>%
  filter(pathway %in% selected_pathways) %>%
  mutate(
    pathway_label = pathway_labels[pathway],
    pathway_class = pathway_class[pathway],
    neglog10_padj = -log10(padj),
    pathway_label = factor(
      pathway_label,
      levels = pathway_labels[selected_pathways]
    )
  )

## 6) Create dot plot

p_dot <- ggplot(plot_df, aes(x = NES, y = pathway_label)) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey40",
    linewidth = 0.7
  ) +
  geom_point(
    aes(size = size, colour = neglog10_padj),
    alpha = 0.95
  ) +
  scale_colour_gradientn(
    colours = c("#2C7BB6", "#FEE08B", "#D7191C"),
    name = "-log10(padj)"
  ) +
  scale_size_continuous(
    name = "Gene set size",
    range = c(5, 16)
  ) +
  scale_x_continuous(
    limits = c(-0.5, 2.8)
  ) +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = "Pathway"
  ) +
  annotate(
    "text",
    x = 1.5,
    y = 16.35,
    label = "Enriched in 1:1,000 EPR11334",
    fontface = "bold",
    size = 4.8,
    hjust = 0.5
  ) +
  coord_cartesian(clip = "off") +
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
    plot.margin = margin(t = 30, r = 40, b = 30, l = 20)
  )

p_dot

## 7) Save dot plot

png(
  "GSEA_C2_1000_EPR_vs_100000_EPR_dotplot.png",
  width = 3800,
  height = 2400,
  res = 300
)

print(p_dot)

dev.off()


## C2 heatmap: 1:1,000 EPR11334 vs 1:100,000 EPR11334

## 1) Create transformed expression matrix for this comparison

vsd_1000EPR_vs_100000EPR <- vst(
  dds_1000EPR_vs_100000EPR,
  blind = TRUE
)

expr_mat_all <- assay(vsd_1000EPR_vs_100000EPR)
sample_info <- as.data.frame(colData(vsd_1000EPR_vs_100000EPR))

## 2) Identify sample columns

epr_1000_cols <- rownames(sample_info)[
  grepl("1000", sample_info$Condition) &
    !grepl("100000", sample_info$Condition)
]

epr_100000_cols <- rownames(sample_info)[
  grepl("100000", sample_info$Condition)
]

if (length(epr_1000_cols) == 0) {
  stop("No 1:1,000 EPR11334 samples found in sample_info$Condition.")
}

if (length(epr_100000_cols) == 0) {
  stop("No 1:100,000 EPR11334 samples found in sample_info$Condition.")
}

expr_mat <- expr_mat_all[, c(epr_1000_cols, epr_100000_cols), drop = FALSE]

## 3) Map Ensembl IDs to gene symbols

ensembl_ids <- rownames(expr_mat)
ensembl_ids_clean <- sub("\\..*$", "", ensembl_ids)

gene_map_heatmap <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(ensembl_ids_clean),
  keytype = "ENSEMBL",
  columns = c("ENSEMBL", "SYMBOL")
)

gene_map_heatmap <- gene_map_heatmap %>%
  filter(!is.na(SYMBOL), SYMBOL != "") %>%
  distinct(ENSEMBL, .keep_all = TRUE)

expr_df <- as.data.frame(expr_mat)
expr_df$ENSEMBL <- ensembl_ids_clean

expr_df_annot <- merge(expr_df, gene_map_heatmap, by = "ENSEMBL")

expr_df_annot$mean_expr <- rowMeans(
  expr_df_annot[, colnames(expr_mat), drop = FALSE],
  na.rm = TRUE
)

expr_df_annot <- expr_df_annot %>%
  group_by(SYMBOL) %>%
  slice_max(order_by = mean_expr, n = 1, with_ties = FALSE) %>%
  ungroup()

expr_symbol_mat <- as.matrix(expr_df_annot[, colnames(expr_mat), drop = FALSE])
rownames(expr_symbol_mat) <- expr_df_annot$SYMBOL

## 4) Define heatmap genes

genes_heatmap_28 <- c(
  "RPL21", "RPL34", "RPL23A", "RPS27A", "RPS3", "RPS6", "RPL18",
  "ATP5ME", "NDUFB3", "NDUFA1", "NDUFA9", "NDUFS3", "COX6A1", "COX5B",
  "PCNA", "PRIM1", "POLE2", "MCM6", "CDC6", "CCNA2", "RPA1",
  "CDH2", "LUM", "SPARC", "ITGAV", "PTPRZ1", "CA9", "ALDH1A2"
)

## 5) Check genes present

genes_present <- genes_heatmap_28[
  genes_heatmap_28 %in% rownames(expr_symbol_mat)
]

genes_missing <- setdiff(
  genes_heatmap_28,
  rownames(expr_symbol_mat)
)

cat("Number of requested genes found:", length(genes_present), "\n")
cat("Missing genes:\n")
print(genes_missing)

## 6) Subset expression matrix

expr_mat_heatmap <- expr_symbol_mat[genes_present, , drop = FALSE]

## 7) Scale expression by gene

expr_mat_heatmap_scaled <- t(scale(t(expr_mat_heatmap)))
expr_mat_heatmap_scaled[is.na(expr_mat_heatmap_scaled)] <- 0

## 8) Create column annotation

annotation_col <- data.frame(
  Condition = c(
    rep("1:1,000 EPR11334", length(epr_1000_cols)),
    rep("1:100,000 EPR11334", length(epr_100000_cols))
  )
)

rownames(annotation_col) <- colnames(expr_mat_heatmap_scaled)

## 9) Create row annotation

gene_group_df <- data.frame(
  Gene = genes_heatmap_28,
  `Functional group` = c(
    rep("Translation", 7),
    rep("Metabolism", 7),
    rep("Proliferation", 7),
    rep("Motility", 7)
  ),
  Pathway = c(
    rep("Ribosome / translation", 7),
    rep("Oxidative phosphorylation", 7),
    rep("DNA synthesis / cell cycle", 7),
    rep("Motility / metastasis", 7)
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

annotation_row <- gene_group_df[
  match(rownames(expr_mat_heatmap_scaled), gene_group_df$Gene),
  c("Functional group", "Pathway"),
  drop = FALSE
]

rownames(annotation_row) <- rownames(expr_mat_heatmap_scaled)

## 10) Set annotation colours

annotation_colors <- list(
  Condition = c(
    "1:1,000 EPR11334" = "black",
    "1:100,000 EPR11334" = "grey70"
  ),
  `Functional group` = c(
    "Translation" = "#88C999",
    "Metabolism" = "#4CAF50",
    "Proliferation" = "#7FBF7B",
    "Motility" = "#2E8B57"
  ),
  Pathway = c(
    "Ribosome / translation" = "#7B2CBF",
    "Oxidative phosphorylation" = "#9D4EDD",
    "DNA synthesis / cell cycle" = "#C77DFF",
    "Motility / metastasis" = "#5A189A"
  )
)

## 11) Preview heatmap

pheatmap(
  expr_mat_heatmap_scaled,
  annotation_col = annotation_col,
  annotation_row = annotation_row,
  annotation_colors = annotation_colors,
  show_rownames = TRUE,
  show_colnames = FALSE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_row = 12,
  border_color = "grey60"
)

## 12) Save heatmap

png(
  "C2_heatmap_1000_EPR_vs_100000_EPR.png",
  width = 2800,
  height = 3600,
  res = 250
)

ph <- pheatmap(
  expr_mat_heatmap_scaled,
  annotation_col = annotation_col,
  annotation_row = annotation_row,
  annotation_colors = annotation_colors,
  show_rownames = TRUE,
  show_colnames = FALSE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_row = 12,
  border_color = "grey60",
  silent = TRUE
)

grid.newpage()

pushViewport(
  viewport(
    x = 0.5,
    y = 0.5,
    width = 0.98,
    height = 0.90
  )
)

grid.draw(ph$gtable)
popViewport()

dev.off()