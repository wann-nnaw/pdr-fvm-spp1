## =========================================================
## 1. 加载R包
## =========================================================
library(DESeq2)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(grid)

## =========================================================
## 2. 设置输入文件和输出目录
## =========================================================
data_file <- "data/GSE102485_expressed_gene_reads.txt.gz"
out_dir <- "results/bulk_DEG"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

group_cols <- c(Control = "#3FA7D6", PDR = "#F25F5C")
volcano_cols <- c(Down = "#3A86FF", Not_sig = "#C9CDD6", Up = "#FF6B6B")
heat_cols <- colorRampPalette(c("#315C99", "#F7F7F2", "#F25F5C"))(100)

theme_nature <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_line(color = "black", linewidth = 0.35),
      axis.ticks.length = unit(2, "pt"),
      axis.text = element_text(color = "black", size = base_size - 1),
      axis.title = element_text(color = "black", size = base_size),
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0.5, margin = margin(b = 3)),
      plot.subtitle = element_text(size = base_size - 1, color = "grey45", hjust = 0.5, margin = margin(b = 6)),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.key.size = unit(9, "pt"),
      legend.title = element_text(size = base_size - 1, face = "bold"),
      legend.text = element_text(size = base_size - 1, color = "black"),
      legend.position = "right",
      legend.margin = margin(0, 0, 0, 4),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold", hjust = 0.5),
      plot.margin = margin(8, 8, 8, 8)
    )
}

## =========================================================
## 3. 读取GSE102485 raw counts矩阵
##    第1列是Gene，第2-31列是30个样本，后面是基因注释信息
## =========================================================
counts_raw <- read.delim(
  data_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

sample_id <- colnames(counts_raw)[2:31]
exclude_samples <- c("Bp2", "Bp3", "Bp4", "Rp1", "Rp2")
sample_id <- setdiff(sample_id, exclude_samples)

## =========================================================
## 4. 构建分组信息
##    IIr1、IIr2、Nr1、Nr2、Nr3作为Control，排除Bp2、Bp3、Bp4、Rp1、Rp2，其余样本作为PDR
## =========================================================
group_table <- data.frame(
  sample_id = sample_id,
  group = case_when(
    sample_id %in% c("IIr1", "IIr2", "Nr1", "Nr2", "Nr3") ~ "Control",
    TRUE ~ "PDR"
  ),
  stringsAsFactors = FALSE
)
group_table$group <- factor(group_table$group, levels = c("Control", "PDR"))
rownames(group_table) <- group_table$sample_id
write.csv(group_table, file.path(out_dir, "group_table.csv"), row.names = FALSE)

## =========================================================
## 5. 提取基因注释信息
##    counts文件后面存在重复Symbol/Chr/GeneType列，这里固定取第32-35列
## =========================================================
gene_anno <- counts_raw[, c(1, 32:35)]
colnames(gene_anno) <- c("Gene", "Symbol", "Chr", "GeneType", "Description")

## =========================================================
## 6. 提取counts矩阵并过滤低表达基因
##    DESeq2输入需要整数counts矩阵
## =========================================================
count_mat <- counts_raw[, sample_id]
count_mat <- as.matrix(count_mat)
mode(count_mat) <- "numeric"
count_mat <- round(count_mat)
rownames(count_mat) <- counts_raw$Gene

keep <- rowSums(count_mat) >= 10
count_mat <- count_mat[keep, ]
gene_anno <- gene_anno[keep, ]

write.csv(
  data.frame(Gene = rownames(count_mat), count_mat, check.names = FALSE),
  file.path(out_dir, "filtered_count_matrix.csv"),
  row.names = FALSE
)
write.csv(gene_anno, file.path(out_dir, "filtered_gene_annotation.csv"), row.names = FALSE)
saveRDS(
  list(
    count_mat = count_mat,
    gene_anno = gene_anno,
    group_table = group_table,
    sample_id = sample_id
  ),
  file.path(out_dir, "filtered_bulk_counts_inputs.rds")
)

## =========================================================
## 7. DESeq2差异分析：PDR vs Control
## =========================================================
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = group_table,
  design = ~ group
)

## =========================================================
## 8. 整理差异分析结果并筛选DEG
##    阈值：padj < 0.05 且 abs(log2FoldChange) > 0.5
## =========================================================
dds <- DESeq(dds)
res <- results(dds, contrast = c("group", "PDR", "Control"))
res_df <- as.data.frame(res)
res_df$Gene <- rownames(res_df)
res_df <- left_join(res_df, gene_anno, by = "Gene")
res_df <- res_df[, c("Gene", "Symbol", "Chr", "GeneType", "Description", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
res_df <- res_df[order(res_df$padj), ]

deg_df <- res_df %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 0.5)

write.csv(res_df, file.path(out_dir, "DESeq2_all_results.csv"), row.names = FALSE)
write.csv(deg_df, file.path(out_dir, "DESeq2_DEG_padj0.05_logFC0.5.csv"), row.names = FALSE)

## =========================================================
## 9. bulk差异基因与巨噬细胞4亚型交集基因取交集
##    bulk结果使用Symbol列，巨噬细胞文件使用gene列
## =========================================================
macrophage_file <- "results/venn/macrophage_4subtype_intersect_genes_direction_classified.csv"
macrophage_gene <- read.csv(macrophage_file, check.names = FALSE, stringsAsFactors = FALSE)
deg_macrophage_intersect <- inner_join(deg_df, macrophage_gene, by = c("Symbol" = "gene"))
write.csv(
  deg_macrophage_intersect,
  file.path(out_dir, "DESeq2_DEG_intersect_macrophage_4subtype_genes.csv"),
  row.names = FALSE
)

deg_gene_set <- unique(deg_df$Symbol[!is.na(deg_df$Symbol) & deg_df$Symbol != ""])
macrophage_gene_set <- unique(macrophage_gene$gene[!is.na(macrophage_gene$gene) & macrophage_gene$gene != ""])
intersect_gene_set <- intersect(deg_gene_set, macrophage_gene_set)
deg_only_n <- length(setdiff(deg_gene_set, macrophage_gene_set))
macrophage_only_n <- length(setdiff(macrophage_gene_set, deg_gene_set))
intersect_n <- length(intersect_gene_set)

theta <- seq(0, 2 * pi, length.out = 300)
venn_circle <- rbind(
  data.frame(x = -0.45 + 1.15 * cos(theta), y = 1.05 * sin(theta), set = "Bulk DEGs"),
  data.frame(x = 0.45 + 1.15 * cos(theta), y = 1.05 * sin(theta), set = "Macrophage DEGs")
)

p_venn <- ggplot(venn_circle, aes(x = x, y = y, fill = set, group = set)) +
  geom_polygon(
    alpha = 0.55,
    color = "white",
    linewidth = 1.2
  ) +
  annotate(
    "text",
    x = -0.95, y = 0,
    label = deg_only_n,
    size = 8,
    fontface = "bold",
    color = "black"
  ) +
  annotate(
    "text",
    x = 0, y = 0,
    label = intersect_n,
    size = 8,
    fontface = "bold",
    color = "black"
  ) +
  annotate(
    "text",
    x = 0.95, y = 0,
    label = macrophage_only_n,
    size = 8,
    fontface = "bold",
    color = "black"
  ) +
  annotate(
    "text",
    x = -0.85, y = 1.25,
    label = "Bulk DEGs",
    size = 6,
    fontface = "bold",
    color = "black"
  ) +
  annotate(
    "text",
    x = 0.85, y = 1.25,
    label = "Macrophage DEGs",
    size = 6,
    fontface = "bold",
    color = "black"
  ) +
  labs(
    title = "Venn plot"
  ) +
  scale_fill_manual(
    values = c(
      "Bulk DEGs" = "#F25F5C",
      "Macrophage DEGs" = "#3FA7D6"
    )
  ) +
  coord_equal(
    xlim = c(-1.9, 1.9),
    ylim = c(-1.25, 1.55),
    expand = FALSE
  ) +
  theme_void() +
  theme(
    plot.title = element_text(
      size = 20,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 14),
      color = "black"
    ),
    legend.position = "none",
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(10, 16, 10, 10)
  )

ggsave(
  file.path(out_dir, "venn_DEG_intersect_macrophage_genes.pdf"),
  p_venn,
  width = 6.2,
  height = 5.2,
  bg = "white"
)

ggsave(
  file.path(out_dir, "venn_DEG_intersect_macrophage_genes.png"),
  p_venn,
  width = 6.2,
  height = 5.2,
  dpi = 300,
  bg = "white"
)
## =========================================================
## 10. 火山图
## =========================================================
plot_df <- res_df %>% filter(!is.na(padj), !is.na(log2FoldChange))
plot_df$change <- "Not_sig"
plot_df$change[plot_df$padj < 0.05 & plot_df$log2FoldChange > 0.5] <- "Up"
plot_df$change[plot_df$padj < 0.05 & plot_df$log2FoldChange < -0.5] <- "Down"
plot_df$neg_log10_padj <- -log10(plot_df$padj)
plot_df$neg_log10_padj[is.infinite(plot_df$neg_log10_padj)] <- max(plot_df$neg_log10_padj[is.finite(plot_df$neg_log10_padj)], na.rm = TRUE)
label_df <- plot_df %>%
  filter(change != "Not_sig", !is.na(Symbol), Symbol != "") %>%
  arrange(padj) %>%
  slice_head(n = 12)

x_lim <- max(abs(plot_df$log2FoldChange), na.rm = TRUE)
x_lim <- ceiling(x_lim * 10) / 10

p_volcano <- ggplot(
  plot_df,
  aes(x = log2FoldChange, y = neg_log10_padj)
) +
  geom_point(
    aes(color = change),
    size = 4.8,
    alpha = 0.85,
    stroke = 0
  ) +
  geom_vline(
    xintercept = c(-0.5, 0.5),
    linetype = "dashed",
    color = "grey65",
    linewidth = 0.35
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "grey65",
    linewidth = 0.35
  ) +
  ggrepel::geom_text_repel(
    data = label_df,
    aes(label = Symbol),
    size = 3.8,
    color = "black",
    box.padding = 0.35,
    point.padding = 0.25,
    min.segment.length = 0,
    segment.color = "grey55",
    segment.linewidth = 0.25,
    max.overlaps = Inf,
    fontface = "bold"
  ) +
  scale_color_manual(
    values = volcano_cols,
    name = NULL,
    drop = FALSE
  ) +
  coord_cartesian(
    xlim = c(-x_lim, x_lim),
    clip = "off"
  ) +
  labs(
    title = "Volcano plot (GSE102485 PDR vs Control)",
    x = "log2 fold change",
    y = "-log10 adjusted P"
  ) +
  theme_nature(base_size = 11) +
  theme(
    legend.position = "top",
    legend.justification = "center",
    legend.title = element_blank(),
    legend.text = element_text(
      size = 12,
      face = "bold",
      color = "black"
    ),
    
    plot.title = element_text(
      size = 20,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 14),
      color = "black"
    ),
    
    axis.text.x = element_text(
      size = 14,
      face = "bold",
      color = "black"
    ),
    axis.text.y = element_text(
      size = 14,
      face = "bold",
      color = "black"
    ),
    axis.title.x = element_text(
      size = 16,
      face = "bold",
      color = "black"
    ),
    axis.title.y = element_text(
      size = 16,
      face = "bold",
      color = "black"
    ),
    
    panel.grid = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.5
    ),
    
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(10, 16, 10, 10)
  )

ggsave(
  file.path(out_dir, "volcano_DESeq2_PDR_vs_Control.pdf"),
  p_volcano,
  width = 6.2,
  height = 6.2,
  bg = "white"
)

ggsave(
  file.path(out_dir, "volcano_DESeq2_PDR_vs_Control.png"),
  p_volcano,
  width = 6.2,
  height = 6.2,
  dpi = 300,
  bg = "white"
)
## =========================================================
## 11. vst标准化，用于PCA和热图
## =========================================================
vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)

## =========================================================
## 12. PCA图
## =========================================================
pca <- prcomp(t(vsd_mat))
percent_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pca_df <- data.frame(
  name = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  group = group_table[rownames(pca$x), "group"],
  stringsAsFactors = FALSE
)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  stat_ellipse(
    aes(fill = group),
    type = "norm",
    geom = "polygon",
    alpha = 0.13,
    color = NA
  ) +
  geom_point(
    aes(fill = group),
    shape = 21,
    size = 4.8,
    stroke = 1.1,
    color = "white"
  ) +
  ggrepel::geom_text_repel(
    aes(label = name),
    size = 3.8,
    color = "black",
    fontface = "bold",
    box.padding = 0.35,
    point.padding = 0.25,
    min.segment.length = 0,
    segment.color = "grey55",
    segment.linewidth = 0.25,
    max.overlaps = Inf
  ) +
  scale_fill_manual(
    values = group_cols,
    name = NULL
  ) +
  labs(
    title = "PCA plot (GSE102485)",
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance")
  ) +
  theme_nature(base_size = 11) +
  theme(
    legend.position = "top",
    legend.justification = "center",
    legend.title = element_blank(),
    legend.text = element_text(
      size = 12,
      face = "bold",
      color = "black"
    ),
    
    plot.title = element_text(
      size = 20,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 14),
      color = "black"
    ),
    
    axis.text.x = element_text(
      size = 14,
      face = "bold",
      color = "black"
    ),
    axis.text.y = element_text(
      size = 14,
      face = "bold",
      color = "black"
    ),
    axis.title.x = element_text(
      size = 16,
      face = "bold",
      color = "black"
    ),
    axis.title.y = element_text(
      size = 16,
      face = "bold",
      color = "black"
    ),
    
    panel.grid = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.5
    ),
    
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(10, 16, 10, 10)
  )

ggsave(
  file.path(out_dir, "PCA_DESeq2_vst.pdf"),
  p_pca,
  width = 6.2,
  height = 6.2,
  bg = "white"
)

ggsave(
  file.path(out_dir, "PCA_DESeq2_vst.png"),
  p_pca,
  width = 6.2,
  height = 6.2,
  dpi = 300,
  bg = "white"
)

## =========================================================
## 13. Top 50 DEG热图
## =========================================================
top_gene <- deg_df %>%
  arrange(padj) %>%
  slice_head(n = 50) %>%
  pull(Gene)

heat_mat <- vsd_mat[top_gene, ]
heat_mat <- heat_mat - rowMeans(heat_mat)
rownames(heat_mat) <- deg_df$Symbol[match(top_gene, deg_df$Gene)]
rownames(heat_mat)[is.na(rownames(heat_mat)) | rownames(heat_mat) == ""] <- top_gene[is.na(rownames(heat_mat)) | rownames(heat_mat) == ""]
annotation_col <- data.frame(group = group_table$group)
rownames(annotation_col) <- group_table$sample_id
annotation_colors <- list(group = group_cols)
annotation_col <- annotation_col[colnames(heat_mat), , drop = FALSE]

annotation_col$group <- factor(
  annotation_col$group,
  levels = names(annotation_colors$group)
)

## 颜色图例刻度
legend_breaks <- c(-6, -4, -2, 0, 2, 4)

col_fun <- circlize::colorRamp2(
  seq(min(legend_breaks), max(legend_breaks), length.out = length(heat_cols)),
  heat_cols
)

## 控制标题和主图距离
ht_opt(RESET = TRUE)
ht_opt(TITLE_PADDING = unit(c(6, 14), "points"))

## 顶部分组注释
top_anno <- HeatmapAnnotation(
  group = annotation_col$group,
  col = list(group = annotation_colors$group),
  show_annotation_name = FALSE,
  simple_anno_size = unit(0.45, "cm"),
  annotation_legend_param = list(
    group = list(
      title = "group",
      title_gp = gpar(
        fontsize = 16,
        fontface = "bold",
        col = "black"
      ),
      labels_gp = gpar(
        fontsize = 16,
        fontface = "bold",
        col = "black"
      ),
      grid_width = unit(0.55, "cm"),
      grid_height = unit(0.55, "cm")
    )
  )
)

## 热图主体
p_heatmap <- Heatmap(
  heat_mat,
  name = "Expression",
  col = col_fun,
  
  top_annotation = top_anno,
  
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  
  show_column_names = TRUE,
  show_row_names = TRUE,
  
  row_names_side = "right",
  row_dend_side = "left",
  column_dend_side = "top",
  
  rect_gp = gpar(col = NA),
  border = FALSE,
  
  column_title = "Heatmap (PDR vs Control)",
  column_title_gp = gpar(
    fontsize = 28,
    fontface = "bold",
    col = "black"
  ),
  
  row_names_gp = gpar(
    fontsize = 12,
    fontface = "bold",
    col = "black"
  ),
  
  column_names_gp = gpar(
    fontsize = 14,
    fontface = "bold",
    col = "black"
  ),
  column_names_rot = 45,
  
  heatmap_legend_param = list(
    title = "",
    at = legend_breaks,
    labels = as.character(legend_breaks),
    labels_gp = gpar(
      fontsize = 16,
      fontface = "bold",
      col = "black"
    ),
    legend_height = unit(5.0, "cm"),
    grid_width = unit(0.45, "cm")
  )
)

pdf(
  file.path(out_dir, "heatmap_top50_DEG.pdf"),
  width = 10.5,
  height = 10
)

draw(
  p_heatmap,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legends = TRUE
)

dev.off()


png(
  file.path(out_dir, "heatmap_top50_DEG.png"),
  width = 3300,
  height = 3000,
  res = 300
)

draw(
  p_heatmap,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legends = TRUE
)

dev.off()
