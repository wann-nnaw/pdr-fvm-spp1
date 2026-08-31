## =========================================================
## 1. 加载R包
## =========================================================
library(readxl)
library(DESeq2)
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

## =========================================================
## 2. 设置输入文件、输出目录和目标基因
## =========================================================
out_dir <- "results/SPP1"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

disease_name <- "Proliferative diabetic retinopathy (PDR)"
target_genes <- c("C3", "FCGR2B", "CLEC5A", "GPNMB", "CD81")
spp1_signature_genes <- c("SPP1", target_genes)
pathway_required_gene <- "SPP1"
macrophage_types <- c("LAM", "IAM", "OSM", "TRM")

go_files <- c(
  IAM = "results/go_kegg/IAM_vs_others/GO_filtered_IAM_vs_others.csv",
  OSM = "results/go_kegg/OSM_vs_others/GO_filtered_OSM_vs_others.csv",
  LAM = "results/go_kegg/LAM_vs_others/GO_filtered_LAM_vs_others.csv",
  TRM = "results/go_kegg/TRM_vs_others/GO_filtered_TRM_vs_others.csv"
)

bulk_filtered_input_file <- "results/bulk_DEG/filtered_bulk_counts_inputs.rds"
bulk_deg_file <- "results/bulk_DEG/DESeq2_all_results.csv"
seurat_file <- "rda/UMAP_annotated.rda"

group_cols <- c(Control = "#3FA7D6", PDR = "#F25F5C")
subtype_cols <- c(LAM = "#2F6F73", IAM = "#B65C4A", OSM = "#6E5F8F", TRM = "#C49A3A")
gene_cols <- c(
  SPP1 = "#8E44AD",
  C3 = "#2E86AB",
  FCGR2B = "#F18F01",
  CLEC5A = "#C73E1D",
  GPNMB = "#3B8C5A",
  CD81 = "#5C677D"
)

white_plot_background <- theme(
  plot.background = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA),
  legend.background = element_rect(fill = "white", color = NA),
  legend.box.background = element_rect(fill = "white", color = NA),
  strip.background = element_rect(fill = "white", color = "grey80")
)

## =========================================================
## 3. 读取GO文件并筛选目标基因相关通路
##    部分文件后缀为csv但实际是xlsx，因此根据文件头PK判断格式
## =========================================================
read_go_file <- function(file) {
  file_head <- readBin(file, what = "raw", n = 2)
  if (identical(file_head, charToRaw("PK"))) {
    as.data.frame(read_excel(file))
  } else {
    read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
  }
}

go_hits_all <- data.frame()

for (subtype in names(go_files)) {
  go_df <- read_go_file(go_files[[subtype]])
  go_df$subtype <- subtype
  go_df$disease <- disease_name
  go_df$spp1_present <- sapply(go_df$geneID, function(x) {
    genes_in_term <- strsplit(as.character(x), "/")[[1]]
    pathway_required_gene %in% genes_in_term
  })
  go_df$other_target_genes <- sapply(go_df$geneID, function(x) {
    genes_in_term <- strsplit(as.character(x), "/")[[1]]
    paste(intersect(genes_in_term, target_genes), collapse = "/")
  })
  go_df$other_target_count <- sapply(go_df$other_target_genes, function(x) {
    ifelse(x == "", 0, length(strsplit(x, "/")[[1]]))
  })
  go_df$target_genes <- ifelse(
    go_df$other_target_genes == "",
    pathway_required_gene,
    paste(pathway_required_gene, go_df$other_target_genes, sep = "/")
  )
  go_df$target_count <- go_df$other_target_count + 1
  go_hits_all <- rbind(go_hits_all, go_df[go_df$spp1_present & go_df$other_target_count >= 1, ])
}

go_hits_all$GO_link <- paste0("https://www.ebi.ac.uk/QuickGO/term/", go_hits_all$ID)
go_hits_all$p.adjust <- as.numeric(go_hits_all$p.adjust)
go_hits_all$neg_log10_padj <- -log10(go_hits_all$p.adjust)

go_desc <- tolower(go_hits_all$Description)
go_hits_all$SPP1_related_category <- case_when(
  grepl("phagocytosis|endocytosis|engulfment", go_desc) ~ "Phagocytosis / endocytosis",
  grepl("immune response.*receptor|cell surface receptor", go_desc) ~ "Immune receptor signaling",
  grepl("inflammatory|humoral|complement|acute inflammatory", go_desc) ~ "Inflammation / complement",
  grepl("adhesion|migration|chemotaxis", go_desc) ~ "Adhesion / migration",
  grepl("cytokine", go_desc) ~ "Cytokine activity",
  grepl("lipid localization", go_desc) ~ "Lipid localization",
  grepl("osteoblast|ossification", go_desc) ~ "Tissue remodeling",
  grepl("response to external stimulus|cell development", go_desc) ~ "Cell development / stimulus response",
  grepl("integrin binding|extracellular matrix|tertiary granule|lysosome|vacuolar", go_desc) ~ "Integrin / lysosome / ECM",
  grepl("lymphocyte proliferation|leukocyte proliferation|mononuclear cell proliferation", go_desc) ~ "Leukocyte proliferation",
  TRUE ~ "Other"
)

go_hits_selected <- go_hits_all %>%
  arrange(p.adjust)

write.csv(go_hits_all, file.path(out_dir, "SPP1_related_GO_hits_all.csv"), row.names = FALSE)
write.csv(go_hits_selected, file.path(out_dir, "SPP1_related_GO_hits_selected.csv"), row.names = FALSE)

## =========================================================
## 4. SPP1相关GO通路气泡图和柱形图
## =========================================================
go_plot_df <- go_hits_selected %>%
  arrange(p.adjust, desc(target_count))

go_plot_df$Description_short <- ifelse(
  nchar(go_plot_df$Description) > 55,
  paste0(substr(go_plot_df$Description, 1, 52), "..."),
  go_plot_df$Description
)
go_plot_df$Description_plot <- ifelse(
  duplicated(go_plot_df$Description_short) | duplicated(go_plot_df$Description_short, fromLast = TRUE),
  paste0(go_plot_df$Description_short, " (", go_plot_df$subtype, ")"),
  go_plot_df$Description_short
)
go_plot_df$Description_plot <- factor(
  go_plot_df$Description_plot,
  levels = rev(go_plot_df$Description_plot)
)

p_go_bubble <- ggplot(
  go_plot_df,
  aes(x = subtype, y = Description_plot, size = target_count, color = neg_log10_padj)
) +
  geom_point(alpha = 0.9) +
  scale_color_gradient(low = "#74A9CF", high = "#B2182B") +
  scale_size(range = c(3, 9)) +
  labs(
    title = "SPP1-related GO pathways in macrophage subtypes",
    subtitle = disease_name,
    x = "Macrophage subtype",
    y = NULL,
    color = "-log10(adj.P)",
    size = "Target genes"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey35"),
    axis.text.y = element_text(size = 8),
    panel.grid.minor = element_blank()
  ) +
  white_plot_background

ggsave(file.path(out_dir, "SPP1_related_GO_bubble.pdf"), p_go_bubble, width = 8, height = 5.5, bg = "white")
ggsave(file.path(out_dir, "SPP1_related_GO_bubble.png"), p_go_bubble, width = 8, height = 5.5, dpi = 300, bg = "white")

p_go_bar <- ggplot(
  go_plot_df,
  aes(x = Description_plot, y = neg_log10_padj, fill = neg_log10_padj)
) +
  geom_col(width = 0.75, color = "black", linewidth = 0.2) +
  geom_text(
    aes(label = subtype),
    hjust = -0.15,
    size = 3.7,
    color = "#243447",
    fontface = "bold"
  ) +
  coord_flip() +
  scale_fill_gradient(low = "#D8ECF3", high = "#D95F4C") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(
    title = "Selected SPP1-related GO terms",
    subtitle = disease_name,
    x = NULL,
    y = "-log10 adjusted P value",
    fill = "-log10(adj.P)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey35"),
    axis.text.y = element_text(size = 8),
    legend.position = "right",
    panel.grid.minor = element_blank()
  ) +
  white_plot_background

ggsave(file.path(out_dir, "SPP1_related_GO_barplot.pdf"), p_go_bar, width = 8, height = 5.5, bg = "white")
ggsave(file.path(out_dir, "SPP1_related_GO_barplot.png"), p_go_bar, width = 8, height = 5.5, dpi = 300, bg = "white")

## =========================================================
## 5. 读取05_bulk_RNA_DESeq2.R保存的bulk counts和PDR / Control分组
## =========================================================
if (!file.exists(bulk_filtered_input_file)) {
  stop(
    "Cannot find ", bulk_filtered_input_file,
    ". Please run 05_bulk_RNA_DESeq2.R first."
  )
}

bulk_filtered_input <- readRDS(bulk_filtered_input_file)
count_mat <- bulk_filtered_input$count_mat
gene_anno <- bulk_filtered_input$gene_anno
group_table <- bulk_filtered_input$group_table
sample_id <- bulk_filtered_input$sample_id

group_table$group <- factor(group_table$group, levels = c("Control", "PDR"))
rownames(group_table) <- group_table$sample_id
count_mat <- count_mat[, sample_id, drop = FALSE]

## =========================================================
## 6. bulk vst标准化并计算SPP1 signature score
## =========================================================
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = group_table,
  design = ~ group
)
dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)

bulk_gene_expr <- as.data.frame(vsd_mat)
bulk_gene_expr$Gene <- rownames(bulk_gene_expr)
bulk_gene_expr <- left_join(gene_anno, bulk_gene_expr, by = "Gene")

bulk_spp1_expr <- bulk_gene_expr %>%
  filter(Symbol %in% spp1_signature_genes) %>%
  dplyr::select(Symbol, all_of(sample_id)) %>%
  group_by(Symbol) %>%
  summarise(across(everything(), mean), .groups = "drop")

bulk_expr_mat <- as.matrix(bulk_spp1_expr[, sample_id])
rownames(bulk_expr_mat) <- bulk_spp1_expr$Symbol

bulk_z_mat <- t(scale(t(bulk_expr_mat)))
bulk_z_mat[is.na(bulk_z_mat)] <- 0

bulk_score <- data.frame(
  sample_id = colnames(bulk_z_mat),
  SPP1_signature_score = colMeans(bulk_z_mat),
  stringsAsFactors = FALSE
) %>%
  left_join(group_table, by = "sample_id")

write.csv(bulk_score, file.path(out_dir, "bulk_SPP1_signature_score.csv"), row.names = FALSE)

signature_pvalue <- wilcox.test(
  SPP1_signature_score ~ group,
  data = bulk_score
)$p.value
signature_pvalue_df <- data.frame(
  comparison = "PDR_vs_Control",
  method = "Wilcoxon rank-sum test on VST z-score signature",
  p_value = signature_pvalue,
  p_label = paste0("Wilcoxon p = ", signif(signature_pvalue, 3))
)
write.csv(signature_pvalue_df, file.path(out_dir, "bulk_SPP1_signature_score_pvalue.csv"), row.names = FALSE)

signature_y <- max(bulk_score$SPP1_signature_score, na.rm = TRUE)
signature_range <- diff(range(bulk_score$SPP1_signature_score, na.rm = TRUE))
signature_bracket_y <- signature_y + signature_range * 0.18
signature_text_y <- signature_y + signature_range * 0.28

p_signature <- ggplot(bulk_score, aes(x = group, y = SPP1_signature_score, fill = group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.8, color = "black") +
  geom_jitter(width = 0.12, size = 2.6, shape = 21, color = "black") +
  annotate("segment", x = 1, xend = 2, y = signature_bracket_y, yend = signature_bracket_y) +
  annotate("segment", x = 1, xend = 1, y = signature_bracket_y, yend = signature_bracket_y - signature_range * 0.04) +
  annotate("segment", x = 2, xend = 2, y = signature_bracket_y, yend = signature_bracket_y - signature_range * 0.04) +
  annotate("text", x = 1.5, y = signature_text_y, label = signature_pvalue_df$p_label, size = 4.2, fontface = "bold") +
  scale_fill_manual(values = group_cols) +
  coord_cartesian(ylim = c(min(bulk_score$SPP1_signature_score), signature_text_y + signature_range * 0.08), clip = "off") +
  labs(
    title = "Bulk SPP1 signature score",
    subtitle = disease_name,
    x = NULL,
    y = "Mean z-score of SPP1 signature genes"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey35"),
    legend.position = "none",
    panel.grid.minor = element_blank()
  ) +
  white_plot_background

ggsave(file.path(out_dir, "bulk_SPP1_signature_score_boxplot.pdf"), p_signature, width = 4.6, height = 4.6, bg = "white")
ggsave(file.path(out_dir, "bulk_SPP1_signature_score_boxplot.png"), p_signature, width = 4.6, height = 4.6, dpi = 300, bg = "white")

## =========================================================
## 7. bulk中SPP1相关基因表达量箱线图
## =========================================================
bulk_expr_long <- bulk_spp1_expr %>%
  pivot_longer(cols = all_of(sample_id), names_to = "sample_id", values_to = "expression") %>%
  left_join(group_table, by = "sample_id")
bulk_expr_long$Symbol <- factor(bulk_expr_long$Symbol, levels = spp1_signature_genes)

write.csv(bulk_expr_long, file.path(out_dir, "bulk_PDR_Control_SPP1_genes_expression.csv"), row.names = FALSE)

bulk_deg_result <- read.csv(bulk_deg_file, check.names = FALSE, stringsAsFactors = FALSE)
bulk_gene_pvalue <- bulk_expr_long %>%
  group_by(Symbol) %>%
  summarise(
    y_max = max(expression, na.rm = TRUE),
    y_min = min(expression, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    bulk_deg_result %>%
      filter(Symbol %in% spp1_signature_genes) %>%
      dplyr::select(Symbol, Gene, log2FoldChange, pvalue, padj),
    by = "Symbol"
  ) %>%
  mutate(
    y_range = y_max - y_min,
    y_range = ifelse(y_range == 0, 1, y_range),
    bracket_y = y_max + y_range * 0.20,
    text_y = y_max + y_range * 0.33,
    p_label = paste0("padj = ", signif(padj, 3))
  )
write.csv(bulk_gene_pvalue, file.path(out_dir, "bulk_PDR_Control_SPP1_genes_DESeq2_pvalues.csv"), row.names = FALSE)

p_bulk_gene <- ggplot(bulk_expr_long, aes(x = group, y = expression, fill = group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.8, color = "black") +
  geom_jitter(width = 0.12, size = 1.7, shape = 21, color = "black") +
  geom_segment(
    data = bulk_gene_pvalue,
    aes(x = 1, xend = 2, y = bracket_y, yend = bracket_y),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = bulk_gene_pvalue,
    aes(x = 1, xend = 1, y = bracket_y, yend = bracket_y - y_range * 0.05),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = bulk_gene_pvalue,
    aes(x = 2, xend = 2, y = bracket_y, yend = bracket_y - y_range * 0.05),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = bulk_gene_pvalue,
    aes(x = 1.5, y = text_y, label = p_label),
    inherit.aes = FALSE,
    size = 3.6,
    fontface = "bold"
  ) +
  facet_wrap(~ Symbol, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = group_cols) +
  labs(
    title = "Bulk expression of SPP1-related genes",
    subtitle = disease_name,
    x = NULL,
    y = "VST-normalized expression"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey35"),
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  white_plot_background

ggsave(file.path(out_dir, "bulk_PDR_Control_SPP1_genes_boxplot.pdf"), p_bulk_gene, width = 8, height = 6, bg = "white")
ggsave(file.path(out_dir, "bulk_PDR_Control_SPP1_genes_boxplot.png"), p_bulk_gene, width = 8, height = 6, dpi = 300, bg = "white")

## =========================================================
## 8. 读取单细胞对象并提取巨噬细胞亚群
## =========================================================
load(seurat_file)
DefaultAssay(combined_anno) <- "SCT"

macrophage_cells <- subset(combined_anno, subset = cell_type %in% macrophage_types)
macrophage_cells$cell_type <- factor(macrophage_cells$cell_type, levels = macrophage_types)
Idents(macrophage_cells) <- macrophage_cells$cell_type

## =========================================================
## 9. 四类巨噬细胞背景中分别高亮单个亚群的SPP1表达
## =========================================================
make_spp1_palette <- function(high_col) {
  if (requireNamespace("colorspace", quietly = TRUE)) {
    colorRampPalette(c(
      "#F4F6F7",
      colorspace::lighten(high_col, 0.72),
      colorspace::lighten(high_col, 0.42),
      high_col,
      high_col,
      colorspace::darken(high_col, 0.22)
    ))(300)
  } else {
    colorRampPalette(c("#F4F6F7", high_col, high_col))(300)
  }
}

theme_spp1_tsne <- function(base_size = 14) {
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
      plot.title = element_text(
        size = base_size + 8, face = "bold",
        hjust = 0.5, color = "#1A1A2E", margin = margin(b = 3)
      ),
      plot.subtitle = element_text(
        size = base_size + 4, face = "bold",
        hjust = 0.5, margin = margin(b = 6)
      ),
      legend.position = "none",
      plot.margin = margin(8, 8, 8, 8)
    )
}

spp1_feature_df <- FetchData(macrophage_cells, vars = c("cell_type", "SPP1"))
spp1_feature_df$cell <- rownames(spp1_feature_df)
spp1_tsne <- as.data.frame(Embeddings(macrophage_cells, reduction = "tsne"))
spp1_tsne$cell <- rownames(spp1_tsne)
colnames(spp1_tsne)[1:2] <- c("tSNE_1", "tSNE_2")

spp1_feature_df <- left_join(spp1_feature_df, spp1_tsne, by = "cell")
spp1_feature_df$cell_type <- factor(spp1_feature_df$cell_type, levels = macrophage_types)
spp1_feature_cols <- c(
  LAM = "#1565C0",
  IAM = "#00695C",
  OSM = "#BF360C",
  TRM = "#6A1B9A"
)

# 四个面板必须共用同一数值范围、同一色带和同一图例，
# 否则不同色相/明度会让跨面板的表达量比较产生视觉偏差。
spp1_scale_limits <- c(
  0,
  unname(quantile(spp1_feature_df$SPP1, probs = 0.97, na.rm = TRUE))
)
spp1_scale_breaks <- pretty(spp1_scale_limits, n = 4)
spp1_scale_breaks <- spp1_scale_breaks[
  spp1_scale_breaks >= spp1_scale_limits[1] &
    spp1_scale_breaks <= spp1_scale_limits[2]
]
spp1_shared_palette <- c(
  "#F4F6F7", "#FEE8C8", "#FDBB84", "#E34A33", "#8C2D04"
)

spp1_subtype_summary <- spp1_feature_df %>%
  group_by(cell_type) %>%
  summarise(
    n_cells = n(),
    mean_SPP1 = mean(SPP1, na.rm = TRUE),
    median_SPP1 = median(SPP1, na.rm = TRUE),
    pct_SPP1_positive = mean(SPP1 > 0, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  arrange(desc(mean_SPP1))

write.csv(
  spp1_subtype_summary,
  file.path(out_dir, "SPP1_expression_summary_by_macrophage_subtype.csv"),
  row.names = FALSE
)

plot_spp1_one_subtype <- function(subtype, show_legend = TRUE) {
  background_df <- spp1_feature_df
  highlight_df <- spp1_feature_df %>%
    filter(cell_type == subtype) %>%
    arrange(SPP1)
  subtype_mean <- spp1_subtype_summary$mean_SPP1[
    spp1_subtype_summary$cell_type == subtype
  ]
  
  ggplot() +
    geom_point(
      data = background_df,
      aes(x = tSNE_1, y = tSNE_2),
      color = "#D5D8DC",
      size = 1.2,
      alpha = 0.65
    ) +
    geom_point(
      data = highlight_df,
      aes(x = tSNE_1, y = tSNE_2, color = SPP1),
      size = 2,
      alpha = 0.95
    ) +
    scale_color_gradientn(
      colors = spp1_shared_palette,
      limits = spp1_scale_limits,
      breaks = spp1_scale_breaks,
      oob = scales::squish,
      name = "SPP1 expression",
      guide = guide_colourbar(
        title.position = "top",
        barheight = unit(35, "mm"),
        barwidth = unit(4, "mm")
      )
    ) +
    labs(
      title = "SPP1 expression",
      subtitle = sprintf("%s  |  mean = %.2f", subtype, subtype_mean),
      x = "tSNE 1",
      y = "tSNE 2"
    ) +
    theme_spp1_tsne(base_size = 14) +
    theme(
      plot.subtitle = element_text(color = spp1_feature_cols[[subtype]]),
      legend.position = if (show_legend) "right" else "none",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10)
    )
}

for (subtype in macrophage_types) {
  p_spp1_subtype <- plot_spp1_one_subtype(subtype)
  ggsave(
    file.path(out_dir, paste0("SPP1_FeaturePlot_", subtype, "_only_in_4macrophages.pdf")),
    p_spp1_subtype,
    width = 5.2, height = 5.2,
    device = "pdf",
    bg = "white"
  )
  ggsave(
    file.path(out_dir, paste0("SPP1_FeaturePlot_", subtype, "_only_in_4macrophages.png")),
    p_spp1_subtype,
    width = 5.2, height = 5.2,
    dpi = 500,
    bg = "white"
  )
}

# 按各亚型的平均 SPP1 表达量从高到低排列面板。
# ggarrange(common.legend = TRUE) 保证四图只使用一个共享色标。
spp1_panel_order <- as.character(spp1_subtype_summary$cell_type)
spp1_panel_list <- lapply(
  spp1_panel_order,
  function(subtype) plot_spp1_one_subtype(subtype, show_legend = TRUE)
)

p_spp1_four_panel <- ggarrange(
  plotlist = spp1_panel_list,
  ncol = 4,
  nrow = 1,
  common.legend = TRUE,
  legend = "right",
  align = "hv"
)

ggsave(
  file.path(out_dir, "SPP1_FeaturePlot_4macrophages_shared_scale.pdf"),
  p_spp1_four_panel,
  width = 18.5, height = 5.2,
  device = "pdf",
  bg = "white"
)
ggsave(
  file.path(out_dir, "SPP1_FeaturePlot_4macrophages_shared_scale.png"),
  p_spp1_four_panel,
  width = 18.5, height = 5.2,
  dpi = 500,
  bg = "white"
)

## =========================================================
## 10. 巨噬细胞亚群中SPP1相关基因表达量图
## =========================================================
macrophage_expr <- FetchData(macrophage_cells, vars = c("cell_type", spp1_signature_genes))
macrophage_expr$cell <- rownames(macrophage_expr)

macrophage_expr_long <- macrophage_expr %>%
  pivot_longer(cols = all_of(spp1_signature_genes), names_to = "Symbol", values_to = "expression")
macrophage_expr_long$Symbol <- factor(macrophage_expr_long$Symbol, levels = spp1_signature_genes)

write.csv(macrophage_expr_long, file.path(out_dir, "macrophage_subtype_SPP1_genes_expression.csv"), row.names = FALSE)

p_macro_violin <- ggplot(macrophage_expr_long, aes(x = cell_type, y = expression, fill = cell_type)) +
  geom_violin(scale = "width", trim = TRUE, color = "black", linewidth = 0.25) +
  geom_boxplot(width = 0.13, outlier.shape = NA, color = "black", fill = "white", alpha = 0.75) +
  facet_wrap(~ Symbol, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = subtype_cols) +
  labs(
    title = "SPP1-related gene expression in macrophage subtypes",
    subtitle = disease_name,
    x = "Macrophage subtype",
    y = "SCT expression"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey35"),
    legend.position = "none",
    panel.grid.minor = element_blank()
  ) +
  white_plot_background

ggsave(file.path(out_dir, "macrophage_subtype_SPP1_genes_violin.pdf"), p_macro_violin, width = 9, height = 6.5, bg = "white")
ggsave(file.path(out_dir, "macrophage_subtype_SPP1_genes_violin.png"), p_macro_violin, width = 9, height = 6.5, dpi = 300, bg = "white")

p_macro_dot <- DotPlot(
  macrophage_cells,
  features = spp1_signature_genes,
  group.by = "cell_type",
  assay = "SCT"
) +
  scale_color_gradient(low = "#D9EAF7", high = "#B2182B") +
  labs(
    title = "SPP1-related genes across macrophage subtypes",
    subtitle = disease_name,
    x = NULL,
    y = "Macrophage subtype"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey35"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  white_plot_background

ggsave(file.path(out_dir, "macrophage_subtype_SPP1_genes_dotplot.pdf"), p_macro_dot, width = 7, height = 4.8, bg = "white")
ggsave(file.path(out_dir, "macrophage_subtype_SPP1_genes_dotplot.png"), p_macro_dot, width = 7, height = 4.8, dpi = 300, bg = "white")

## =========================================================
## 11. 基于 SPP1 signature score 构建诊断模型及 ROC 曲线
## =========================================================
library(pROC)
library(broom)

# 准备数据（确保分组为二分类因子）
diag_data <- bulk_score %>%
  mutate(group = factor(group, levels = c("Control", "PDR"))) %>%
  dplyr::select(group, score = SPP1_signature_score)

# 逻辑回归模型
diag_model <- glm(group ~ score, data = diag_data, family = binomial)

# 模型摘要（输出到控制台，同时保存至文件）
model_summary <- tidy(diag_model)
model_glance <- glance(diag_model)

write.csv(model_summary, file.path(out_dir, "diagnostic_model_coefficients.csv"), row.names = FALSE)
write.csv(model_glance, file.path(out_dir, "diagnostic_model_glance.csv"), row.names = FALSE)

# 预测概率
diag_data$pred_prob <- predict(diag_model, type = "response")

# ROC 曲线
roc_obj <- roc(diag_data$group, diag_data$pred_prob, levels = c("Control", "PDR"), direction = "<")
auc_val <- as.numeric(auc(roc_obj))
auc_ci <- ci.auc(roc_obj, method = "delong")  # 可选，用于置信区间
roc_pvalue <- model_summary$p.value[model_summary$term == "score"]
roc_p_label <- ifelse(
  roc_pvalue < 0.001,
  "P < 0.001",
  paste0("P = ", signif(roc_pvalue, 3))
)

# 提取 ROC 数据
roc_df <- data.frame(
  fpr = 1 - roc_obj$specificities,
  tpr = roc_obj$sensitivities
) %>%
  arrange(fpr, tpr)

# 绘制 ROC 曲线（风格与现有图形一致）
p_roc_diag <- ggplot(roc_df, aes(x = fpr, y = tpr)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "gray65", linewidth = 0.8) +
  geom_step(direction = "vh", color = "#2F9ED8", linewidth = 1.4) +
  annotate("label", x = 0.73, y = 0.22,
           label = paste0("AUC = ", sprintf("%.3f", auc_val), "\n", roc_p_label),
           fill = "white", color = "#23395B", size = 4.6) +
  coord_equal(xlim = c(-0.05, 1.05), ylim = c(-0.05, 1.05)) +
  scale_x_continuous(breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(
    title = "Diagnostic model based on SPP1 signature score",
    subtitle = paste0(disease_name, " (n = ", nrow(diag_data), ")"),
    x = "1 - Specificity",
    y = "Sensitivity"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, color = "black"),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey35"),
    axis.title = element_text(face = "bold", size = 13, color = "black"),
    axis.text = element_text(size = 11, color = "grey25"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "grey45", linewidth = 1.1)
  ) +
  white_plot_background

ggsave(file.path(out_dir, "diagnostic_ROC_SPP1_signature.pdf"), p_roc_diag, width = 5.5, height = 5.5, bg = "white")
ggsave(file.path(out_dir, "diagnostic_ROC_SPP1_signature.png"), p_roc_diag, width = 5.5, height = 5.5, dpi = 300, bg = "white")

# 保存预测结果
diag_data$sample_id <- rownames(diag_data)
write.csv(diag_data, file.path(out_dir, "diagnostic_predictions.csv"), row.names = FALSE)
