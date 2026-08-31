
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(CellChat)
library(patchwork)
library(scales)

## =========================================================
## 1. 路径、分组和配色
## =========================================================
out_dir      <- "results/SPP_positive"
plot_dir     <- file.path(out_dir, "plots")
deg_dir      <- file.path(out_dir, "DEG_GO_KEGG")
cellchat_dir <- file.path(out_dir, "CellChat_SPP1_focus")
dir.create(plot_dir,     recursive = TRUE, showWarnings = FALSE)
dir.create(deg_dir,      recursive = TRUE, showWarnings = FALSE)
dir.create(cellchat_dir, recursive = TRUE, showWarnings = FALSE)

seurat_file     <- "rda/UMAP_annotated.rda"
macrophage_types <- c("LAM", "IAM", "OSM", "TRM")
spp1_groups     <- c("SPP1- macrophage", "SPP1+ macrophage")
cellchat_idents <- c("SPP1- macrophage", "SPP1+ macrophage", "Endothelial", "Fibroblast")

## ── 配色方案说明 ──────────────────────────────────────────
##
##  tsne_group_cols  : tSNE 独立配色
##                     SPP1-  #7BAFC7  雾霭蓝（矿石蓝）
##                     SPP1+  #D4784A  铜锈橙（熔岩橙棕）
##
##  spp1_group_cols  : 组成条形图（subtype / sample）统一配色
##                     SPP1-  #8BB8B0  苔绿青（哑光水鸭）
##                     SPP1+  #C47960  赤陶砖红（暖褪色橙）
##
##  strength_cols    : CellChat 信号强度图独立配色
##                     SPP1-  #7A96B8  岩板蓝（低饱和钢蓝）
##                     SPP1+  #C4A045  古铜金（琥珀黄棕）
##
##  subtype_cols     : 巨噬细胞亚型（OSM 由原紫色改为林地绿）
##                     LAM  #3A7CA5  钢青蓝
##                     IAM  #B5502A  砖红棕
##                     OSM  #4D8B6F  林地绿（替换原紫色）
##                     TRM  #B08C3A  暖琥珀
##
## ─────────────────────────────────────────────────────────

# tSNE 独立配色
tsne_group_cols <- c(
  "SPP1- macrophage" = "#7BAFC7",
  "SPP1+ macrophage" = "#D4784A"
)

# 组成条形图统一配色（subtype + sample，内部一致）
spp1_group_cols <- c(
  "SPP1- macrophage" = "#8BB8B0",
  "SPP1+ macrophage" = "#C47960"
)

# 信号强度图独立配色
strength_cols <- c(
  "SPP1- macrophage" = "#7A96B8",
  "SPP1+ macrophage" = "#C4A045"
)

# 巨噬细胞亚型配色（OSM 去紫改绿）
subtype_cols <- c(
  LAM = "#3A7CA5",
  IAM = "#B5502A",
  OSM = "#4D8B6F",
  TRM = "#B08C3A"
)

pal_deg <- c(
  "Up"      = "#B2182B",
  "Down"    = "#2166AC",
  "Not_Sig" = "#D6CEC3"
)

go_ontology_colors <- c(
  "BP" = "#4E79A7",
  "MF" = "#59A14F",
  "CC" = "#E15759"
)

kegg_gradient_colors <- c("#6BAED6", "#74C69D", "#F2C14E", "#F28E2B", "#D95F59")

## =========================================================
## 2. 通用绘图函数
## =========================================================
theme_nature <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      panel.background  = element_rect(fill = "white", color = NA),
      plot.background   = element_rect(fill = "white", color = NA),
      panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid        = element_blank(),
      axis.line         = element_blank(),
      axis.ticks        = element_line(color = "black", linewidth = 0.35),
      axis.ticks.length = unit(2, "pt"),
      axis.text         = element_text(color = "black", size = base_size - 1),
      axis.title        = element_text(color = "black", size = base_size),
      plot.title        = element_text(size = base_size + 2, face = "bold",
                                       hjust = 0.5, margin = margin(b = 3)),
      plot.subtitle     = element_text(size = base_size - 1, color = "grey45",
                                       hjust = 0.5, margin = margin(b = 6)),
      legend.background = element_blank(),
      legend.key        = element_blank(),
      legend.key.size   = unit(9, "pt"),
      legend.title      = element_text(size = base_size - 1, face = "bold"),
      legend.text       = element_text(size = base_size - 1, color = "black"),
      legend.position   = "right",
      legend.margin     = margin(0, 0, 0, 4),
      strip.background  = element_blank(),
      strip.text        = element_text(size = base_size, face = "bold", hjust = 0.5),
      plot.margin       = margin(8, 8, 8, 8)
    )
}

save_plot <- function(plot, filename_base, width, height, dpi = 300) {
  ggsave(paste0(filename_base, ".pdf"), plot,
         width = width, height = height, device = "pdf", bg = "white")
  ggsave(paste0(filename_base, ".png"), plot,
         width = width, height = height, dpi = dpi, bg = "white")
}

empty_csv <- function(file) {
  write.csv(data.frame(), file, row.names = FALSE)
}

## =========================================================
## 3. 加载数据并定义 SPP1+/- 巨噬细胞
## =========================================================
load(seurat_file)
DefaultAssay(combined_anno) <- "SCT"

macrophage_cells <- subset(combined_anno, subset = cell_type %in% macrophage_types)
macrophage_cells$cell_type <- factor(macrophage_cells$cell_type, levels = macrophage_types)

spp1_expr   <- FetchData(macrophage_cells, vars = "SPP1")[, 1]
spp1_median <- median(spp1_expr, na.rm = TRUE)
macrophage_cells$SPP1_group <- ifelse(
  spp1_expr > spp1_median,
  "SPP1+ macrophage",
  "SPP1- macrophage"
)
macrophage_cells$SPP1_group <- factor(macrophage_cells$SPP1_group, levels = spp1_groups)
Idents(macrophage_cells) <- macrophage_cells$SPP1_group

group_summary <- data.frame(
  SPP1_median = spp1_median,
  SPP1_group  = names(table(macrophage_cells$SPP1_group)),
  cell_count  = as.integer(table(macrophage_cells$SPP1_group)),
  stringsAsFactors = FALSE
)
group_summary$proportion <- group_summary$cell_count / sum(group_summary$cell_count)
write.csv(group_summary, file.path(out_dir, "SPP1_positive_group_summary.csv"), row.names = FALSE)

subtype_group_counts <- as.data.frame(
  table(macrophage_cells$cell_type, macrophage_cells$SPP1_group),
  stringsAsFactors = FALSE
)
colnames(subtype_group_counts) <- c("cell_type", "SPP1_group", "cell_count")
subtype_group_counts <- subtype_group_counts %>%
  group_by(cell_type) %>%
  mutate(
    subtype_total = sum(cell_count),
    proportion    = cell_count / subtype_total,
    percent       = proportion * 100
  ) %>%
  ungroup()

sample_group_counts <- as.data.frame(
  table(macrophage_cells$sample, macrophage_cells$SPP1_group),
  stringsAsFactors = FALSE
)
colnames(sample_group_counts) <- c("sample", "SPP1_group", "cell_count")
sample_group_counts <- sample_group_counts %>%
  group_by(sample) %>%
  mutate(
    sample_total = sum(cell_count),
    proportion   = cell_count / sample_total,
    percent      = proportion * 100
  ) %>%
  ungroup()

subtype_sample_group_counts <- as.data.frame(
  table(macrophage_cells$sample, macrophage_cells$cell_type, macrophage_cells$SPP1_group),
  stringsAsFactors = FALSE
)
colnames(subtype_sample_group_counts) <- c("sample", "cell_type", "SPP1_group", "cell_count")
subtype_sample_group_counts <- subtype_sample_group_counts %>%
  group_by(sample, cell_type) %>%
  mutate(
    sample_subtype_total = sum(cell_count),
    proportion = ifelse(sample_subtype_total > 0,
                        cell_count / sample_subtype_total, NA_real_),
    percent    = proportion * 100
  ) %>%
  ungroup()

write.csv(subtype_group_counts,
          file.path(out_dir, "SPP1_positive_by_subtype.csv"),        row.names = FALSE)
write.csv(sample_group_counts,
          file.path(out_dir, "SPP1_positive_by_sample.csv"),          row.names = FALSE)
write.csv(subtype_sample_group_counts,
          file.path(out_dir, "SPP1_positive_by_subtype_sample.csv"), row.names = FALSE)

## =========================================================
## 4. 基础可视化
## =========================================================

## 4a. tSNE — 独立配色：雾霭蓝 vs 铜锈橙
p_tsne_group <- DimPlot(
  macrophage_cells,
  reduction = "tsne",
  group.by  = "SPP1_group",
  pt.size   = 2,
  alpha     = 1,
  cols      = tsne_group_cols          # 独立配色，不与条形图共享
) +
  labs(
    title    = "SPP1 macrophage states",
    x = "tSNE_1",
    y = "tSNE_2"
  ) +
  theme_nature(base_size = 12) +
  guides(color = guide_legend(override.aes = list(size = 4), title = NULL))+
  theme(
    legend.position = "top",
    # 主标题大小
    plot.title = element_text(size = 20, face = "bold",margin = margin(b = 14)),
    
    ## x/y 轴刻度文字大小
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    
    ## x/y 轴标题大小
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    
    ## 图例标题大小，比如 pval、Commun. Prob.
    legend.title = element_text(size = 12, face = "bold"),
    
    ## 图例内容大小，比如 3、min、max
    legend.text = element_text(size = 12)
  )

save_plot(p_tsne_group,
          file.path(plot_dir, "SPP1_positive_tsne"), width = 6, height = 6)

## 4b. 亚型组成：堆叠比例图（苔绿青 vs 赤陶砖红）
p_subtype_prop <- ggplot(
  subtype_group_counts,
  aes(x = cell_type, y = percent, fill = SPP1_group)
) +
  geom_col(width = 0.72, color = "black", linewidth = 0.25) +
  scale_fill_manual(values = spp1_group_cols) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.03)), limits = c(0, 100)) +
  labs(
    title = "Subtype composition",
    x     = "Macrophage subtype",
    y     = "Cell proportion (%)",
    fill  = NULL
  ) +
  theme_nature(base_size = 11) +
  theme(
    # 图例色块改成正方形
    legend.key.width  = unit(0.45, "cm"),
    legend.key.height = unit(0.45, "cm"),
    
    legend.position = "top",
    # 主标题大小
    plot.title = element_text(size = 26, face = "bold",margin = margin(b = 14)),
    
    ## x/y 轴刻度文字大小
    axis.text.x = element_text(size = 18, angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 18, face = "bold"),
    
    ## x/y 轴标题大小
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20, face = "bold"),

    ## 图例标题大小，比如 pval、Commun. Prob.
    legend.title = element_text(size = 16, face = "bold"),
    
    ## 图例内容大小，比如 3、min、max
    legend.text = element_text(size = 16),
    legend.margin = margin(b = 8)
  )


## 4c. 亚型组成：并排计数图（颜色与比例图一致）
p_subtype_count <- ggplot(
  subtype_group_counts,
  aes(x = cell_type, y = cell_count, fill = SPP1_group)
) +
  geom_col(position = position_dodge(width = 0.75),
           width = 0.68, color = "black", linewidth = 0.25) +
  scale_fill_manual(values = spp1_group_cols) +
  labs(
    title = "Subtype cell counts",
    x     = "Macrophage subtype",
    y     = "Cell count",
    fill  = NULL
  ) +
  theme_nature(base_size = 11) +
  theme(
    # 图例色块改成正方形
    legend.key.width  = unit(0.45, "cm"),
    legend.key.height = unit(0.45, "cm"),
    
    legend.position = "top",
    # 主标题大小
    plot.title = element_text(size = 26, face = "bold",margin = margin(b = 14)),
    
    ## x/y 轴刻度文字大小
    axis.text.x = element_text(size = 18, angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 18, face = "bold"),
    
    ## x/y 轴标题大小
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20, face = "bold"),
    
    ## 图例标题大小，比如 pval、Commun. Prob.
    legend.title = element_text(size = 16, face = "bold"),
    
    ## 图例内容大小，比如 3、min、max
    legend.text = element_text(size = 16),
    legend.margin = margin(b = 8)
  )


save_plot(p_subtype_prop,
          file.path(plot_dir, "SPP1_positive_subtype_prop"), width = 8, height = 7)

save_plot(p_subtype_count,
          file.path(plot_dir, "SPP1_positive_subtype_count"), width = 8, height = 7)

## 4d. 样本组成：堆叠比例图（颜色与亚型图一致）
p_sample_prop <- ggplot(
  sample_group_counts,
  aes(x = sample, y = percent, fill = SPP1_group)
) +
  geom_col(width = 0.72, color = "black", linewidth = 0.25) +
  scale_fill_manual(values = spp1_group_cols) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.03)), limits = c(0, 100)) +
  labs(
    title = "Sample composition",
    x     = "Sample",
    y     = "Cell proportion (%)",
    fill  = NULL
  ) +
  theme_nature(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"))+
  theme(
    # 图例色块改成正方形
    legend.key.width  = unit(0.45, "cm"),
    legend.key.height = unit(0.45, "cm"),
    
    legend.position = "top",
    # 主标题大小
    plot.title = element_text(size = 26, face = "bold",margin = margin(b = 14)),
    
    ## x/y 轴刻度文字大小
    axis.text.x = element_text(size = 18, angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 18, face = "bold"),
    
    ## x/y 轴标题大小
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20, face = "bold"),
    
    ## 图例标题大小，比如 pval、Commun. Prob.
    legend.title = element_text(size = 16, face = "bold"),
    
    ## 图例内容大小，比如 3、min、max
    legend.text = element_text(size = 16),
    legend.margin = margin(b = 8)
  )
save_plot(p_sample_prop,
          file.path(plot_dir, "SPP1_positive_sample_composition"), width = 8, height =7)

## =========================================================
## 5. 功能模块评分
## =========================================================
module_gene_sets <- list(

  # ── 脂质代谢 ─────────────────────────────────
  Lipid_Metabolism = c(
    "GPNMB", "PLA2G7", "LPL", "APOE", "FABP5",
    "TREM2", "CD36", "ABCA1", "PLIN2", "LIPA"
  ),

  # ── 氧化应激 ─────────────────────────────────
  Oxidative_Stress_Response = c(
    "TXNRD1", "GCLM", "GCLC", "SQSTM1",
    "HMOX1", "NQO1", "SOD2", "PRDX1",
    "SRXN1", "SESN2"
  ),

  # ── 糖酵解 ───────────────────────────────────
  Glycolysis = c(
    "SLC2A1", "SLC2A3", "HK2", "PFKP",
    "ALDOA", "PGK1", "ENO1", "PKM", "LDHA",
    "PFKFB3", "SLC16A3"
  ),

  # ── 氧化磷酸化 ───────────────────────────────
  OXPHOS = c(
    "NDUFB8", "UQCRC2", "COX5A",
    "ATP5F1B", "SDHA", "CS",
    "IDH2", "OGDH", "SUCLA2"
  ),

  # ── 组织驻留特征 ─────────────────────────────
  Tissue_Residency = c(
    "MRC1", "CD163", "FOLR2", "MAF",
    "F13A1", "LYVE1", "C1QA", "C1QB",
    "TIMD4", "VSIG4", "SIGLEC1"
  ),

  # ── 炎症激活 ─────────────────────────────────
  Inflammatory_Activation = c(
    "CCL3", "CCL4", "IL1B", "TNF", "NFKBIA",
    "IL6", "PTGS2", "S100A8", "CXCL2"
  ),

  # ── 内皮细胞互作 ─────────────────────────────
  Endothelial_Crosstalk = c(
    "VEGFA", "OSM", "MIF", "HIF1A",
    "CCL2", "ANGPT2", "CXCL1", "FGF2", "THBS1"
  ),

  # ── ECM重塑 ──────────────────────────────────
  ECM_Remodeling = c(
    "TGFB1", "PDGFB", "MMP9", "TIMP1",
    "FN1", "CTSB", "CTSD",
    "POSTN", "CTGF", "IL13RA2"
  )
)

present_module_genes <- lapply(module_gene_sets,
                               function(x) intersect(x, rownames(macrophage_cells)))
module_gene_summary <- data.frame(
  module             = names(module_gene_sets),
  requested_genes    = vapply(module_gene_sets, paste, character(1), collapse = "/"),
  present_genes      = vapply(present_module_genes, paste, character(1), collapse = "/"),
  present_gene_count = vapply(present_module_genes, length, integer(1)),
  stringsAsFactors   = FALSE
)
write.csv(module_gene_summary,
          file.path(out_dir, "SPP1_positive_module_gene_sets.csv"), row.names = FALSE)

score_one_module <- function(seurat_obj, genes) {
  genes <- intersect(genes, rownames(seurat_obj))
  if (length(genes) == 0) return(rep(NA_real_, ncol(seurat_obj)))
  mat        <- FetchData(seurat_obj, vars = genes)
  mat_scaled <- scale(mat)
  mat_scaled[is.na(mat_scaled)] <- 0
  rowMeans(mat_scaled)
}

module_score_df <- data.frame(
  cell       = colnames(macrophage_cells),
  SPP1_group = macrophage_cells$SPP1_group,
  cell_type  = macrophage_cells$cell_type,
  sample     = macrophage_cells$sample,
  stringsAsFactors = FALSE
)
for (module_name in names(module_gene_sets)) {
  module_score_df[[module_name]] <- score_one_module(
    macrophage_cells, module_gene_sets[[module_name]]
  )
}
write.csv(module_score_df,
          file.path(out_dir, "SPP1_positive_module_scores.csv"), row.names = FALSE)

module_score_long <- module_score_df %>%
  pivot_longer(
    cols      = all_of(names(module_gene_sets)),
    names_to  = "module",
    values_to = "score"
  ) %>%
  mutate(
    SPP1_group = factor(SPP1_group, levels = spp1_groups),
    cell_type  = factor(cell_type,  levels = macrophage_types),
    module     = factor(module,     levels = names(module_gene_sets))
  )

module_pvalues <- module_score_long %>%
  group_by(module) %>%
  summarise(
    p_value = wilcox.test(score ~ SPP1_group)$p.value,
    y       = max(score, na.rm = TRUE) + diff(range(score, na.rm = TRUE)) * 0.14,
    label   = paste0("p = ", signif(p_value, 3)),
    .groups = "drop"
  )
write.csv(module_pvalues,
          file.path(out_dir, "SPP1_positive_module_score_pvalues.csv"), row.names = FALSE)

## 5a. SPP1+/- 模块得分 — 纯箱型图（苔绿青 vs 赤陶砖红）
p_module_box <- ggplot(
  module_score_long,
  aes(x = SPP1_group, y = score, fill = SPP1_group)
) +
  geom_boxplot(
    width         = 0.55,
    outlier.size  = 0.8,
    outlier.alpha = 0.45,
    color         = "black",
    linewidth     = 0.3
  ) +
  geom_text(
    data        = module_pvalues,
    aes(x = 1.5, y = y, label = label),
    inherit.aes = FALSE,
    size        = 4.2,          # 显著性标注字体大小
    fontface    = "bold"
  ) +
  facet_wrap(~ module, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = spp1_group_cols) +
  labs(
    title = "Functional module scores in SPP1+/- macrophages",
    x     = NULL,
    y     = "Mean z-score"
  ) +
  theme_nature(base_size = 10) +
  theme(
    legend.position = "top",
    legend.title    = element_blank(),
    
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5,margin = margin(b = 14)),
    
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    
    axis.text.x = element_text(size = 11, angle = 30, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 13, face = "bold"),
    
    strip.text = element_text(size = 14, face = "bold"),

    legend.text = element_text(size = 13, face = "bold"),
    # 旋转后的横坐标标签需要额外的左右及底部留白，避免首个标签被裁切
    plot.margin = margin(t = 8, r = 8, b = 24, l = 30)
  ) +
  coord_cartesian(clip = "off")

save_plot(p_module_box,
          file.path(plot_dir, "SPP1_positive_module_score_boxplot"), width = 13, height = 8.5)

## 5b. 巨噬亚型模块得分 — 小提琴图（保留原有形式）
p_module_subtype <- ggplot(
  module_score_long,
  aes(x = cell_type, y = score, fill = cell_type)
) +
  geom_violin(
    scale = "width",
    trim = TRUE,
    color = "black",
    linewidth = 0.2,
    alpha = 0.82
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    color = "black",
    fill = "white",
    alpha = 0.72
  ) +
  facet_wrap(~ module, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = subtype_cols) +
  labs(
    title = "Functional module scores across macrophage subtypes",
    x     = NULL,
    y     = "Mean z-score"
  ) +
  theme_nature(base_size = 10) +
  theme(
    legend.position = "top",
    legend.title    = element_blank(),
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5,margin = margin(b = 14)),
    
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    
    axis.text.x = element_text(size = 12, angle = 30, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 13, face = "bold"),
    
    strip.text = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13, face = "bold")
  )
save_plot(p_module_subtype,
          file.path(plot_dir, "SPP1_positive_module_score_by_subtype"), width = 13, height = 7.5)

## =========================================================
## 6. 差异基因分析和火山图
## =========================================================
deg_raw <- FindMarkers(
  object          = macrophage_cells,
  ident.1         = "SPP1+ macrophage",
  ident.2         = "SPP1- macrophage",
  assay           = "SCT",
  recorrect_umi   = FALSE,
  min.pct         = 0.25,
  logfc.threshold = 0.25,
  test.use        = "wilcox"
)
deg_raw$gene <- rownames(deg_raw)
deg_annot <- deg_raw %>%
  mutate(
    change = case_when(
      p_val_adj < 0.05 & avg_log2FC >  0.25 ~ "Up",
      p_val_adj < 0.05 & avg_log2FC < -0.25 ~ "Down",
      TRUE                                    ~ "Not_Sig"
    ),
    neg_log10_padj = -log10(pmax(p_val_adj, .Machine$double.xmin)),
    change = factor(change, levels = c("Down", "Not_Sig", "Up"))
  ) %>%
  arrange(p_val_adj, desc(abs(avg_log2FC)))
write.csv(deg_annot, file.path(deg_dir, "SPP1pos_vs_SPP1neg_DEG.csv"), row.names = FALSE)

label_genes <- deg_annot %>%
  filter(change != "Not_Sig") %>%
  arrange(p_val_adj, desc(abs(avg_log2FC))) %>%
  slice_head(n = 18)

x_lim <- ceiling(max(abs(deg_annot$avg_log2FC), na.rm = TRUE) * 10) / 10

p_volcano <- ggplot(deg_annot, aes(x = avg_log2FC, y = neg_log10_padj)) +
  geom_point(aes(color = change), size = 5, alpha = 0.85, stroke = 0) +
  geom_vline(xintercept = c(-0.25, 0.25),
             linetype = "dashed", color = "grey65", linewidth = 0.35) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", color = "grey65", linewidth = 0.35) +
  ggrepel::geom_text_repel(
    data               = label_genes,
    aes(label = gene),
    size               = 4,
    color              = "black",
    box.padding        = 0.35,
    point.padding      = 0.25,
    max.overlaps       = Inf,
    min.segment.length = 0,
    fontface = "bold"
  ) +
  scale_color_manual(values = pal_deg, drop = FALSE) +
  coord_cartesian(xlim = c(-x_lim, x_lim)) +
  labs(
    title    = "Volcano plot: SPP1+ vs SPP1-",
    x     = "log2 fold change",
    y     = "-log10 adjusted P",
    color = NULL
  ) +
  theme_nature(base_size = 11) +
  theme(
    legend.position = "top",
    plot.title = element_text(size = 20, face = "bold",margin = margin(b = 14)),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    legend.margin = margin(b = 8)
  )

save_plot(p_volcano,
          file.path(plot_dir, "SPP1pos_vs_SPP1neg_volcano"), width = 6, height = 6)

top_heat_genes <- deg_annot %>%
  filter(change != "Not_Sig") %>%
  group_by(change) %>%
  slice_min(p_val_adj, n = 15, with_ties = FALSE) %>%
  ungroup() %>%
  pull(gene) %>%
  unique()

if (length(top_heat_genes) > 1) {
  set.seed(1234)
  heat_cells <- unlist(lapply(spp1_groups, function(group_name) {
    cells <- colnames(macrophage_cells)[macrophage_cells$SPP1_group == group_name]
    sample(cells, size = min(100, length(cells)))
  }))
  heat_df <- FetchData(macrophage_cells[, heat_cells], vars = top_heat_genes)
  heat_scaled <- scale(heat_df)
  heat_scaled[is.na(heat_scaled)] <- 0
  heat_plot_df <- as.data.frame(heat_scaled)
  heat_plot_df$cell <- rownames(heat_plot_df)
  heat_plot_df$SPP1_group <- macrophage_cells$SPP1_group[heat_plot_df$cell]
  heat_plot_long <- heat_plot_df %>%
    pivot_longer(cols = all_of(top_heat_genes), names_to = "gene", values_to = "scaled_expression") %>%
    mutate(
      SPP1_group = factor(SPP1_group, levels = spp1_groups),
      gene = factor(gene, levels = rev(top_heat_genes)),
      cell = factor(cell, levels = heat_cells)
    )
  p_heat <- ggplot(
    heat_plot_long,
    aes(x = cell, y = gene, fill = scaled_expression)
  ) +
    geom_raster() +
    facet_grid(
      . ~ SPP1_group,
      scales = "free_x",
      space = "free_x"
    ) +
    scale_fill_gradientn(
      colors = c(
        "#2166AC",
        "#67A9CF",
        "#F7F7F7",
        "#EF8A62",
        "#B2182B"
      ),
      limits = c(-2, 2),
      oob = squish,
      breaks = c(-2, 0, 2),
      labels = c("-2", "0", "2"),
      guide = guide_colorbar(
        barheight = unit(3.8, "cm"),
        barwidth  = unit(0.35, "cm")
      )
    ) +
    labs(
      title = "Heatmap: SPP1+ vs SPP1−",
      x = NULL,
      y = NULL,
      fill = "Scaled\nexpression"
    ) +
    theme_nature(base_size = 9) +
    theme(
      legend.position = "right",
      legend.justification = "top",
      legend.box.just = "top",
      
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.border = element_blank(),
      
      plot.title = element_text(
        size = 28,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 14)
      ),
      
      axis.text.y = element_text(
        size = 14,
        face = "bold"
      ),
      
      strip.text.x = element_text(
        size = 18,
        face = "bold",
        margin = margin(b = 8)
      ),
      
      ## 只修改标题和数字之间的距离
      legend.title = element_text(
        size = 16,
        face = "bold",
        lineheight = 0.9,
        margin = margin(b = 12)
      ),
      
      legend.text = element_text(
        size = 16,
        face = "bold"
      )
    )
  save_plot(p_heat, file.path(plot_dir, "SPP1pos_vs_SPP1neg_top_DEG_heatmap"), width = 10, height = 8)
}

## =========================================================
## 7. SPP1+ 上调基因的 GO 和 KEGG 富集分析
## =========================================================
sig_up_genes <- deg_annot %>%
  filter(change == "Up", gene != "SPP1") %>%
  arrange(p_val_adj, desc(avg_log2FC)) %>%
  pull(gene) %>%
  unique()

write.csv(data.frame(gene = sig_up_genes),
          file.path(deg_dir, "SPP1pos_up_genes_for_enrichment_excluding_SPP1.csv"),
          row.names = FALSE)

if (length(sig_up_genes) > 0) {
  entrez_df <- suppressMessages(tryCatch(
    bitr(sig_up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db),
    error = function(e) data.frame()
  ))
  write.csv(entrez_df,
            file.path(deg_dir, "SPP1pos_up_genes_entrez.csv"), row.names = FALSE)

  if (nrow(entrez_df) > 0) {
    cluster_gene <- unique(entrez_df$ENTREZID)

    ego <- tryCatch(
      enrichGO(
        gene          = cluster_gene,
        OrgDb         = org.Hs.eg.db,
        keyType       = "ENTREZID",
        ont           = "ALL",
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.2,
        readable      = TRUE,
        pool          = TRUE
      ),
      error = function(e) NULL
    )
    ego_results  <- if (is.null(ego)) data.frame() else as.data.frame(ego)
    ego_filtered <- if (nrow(ego_results) > 0 && "p.adjust" %in% colnames(ego_results)) {
      ego_results %>%
        filter(p.adjust < 0.05, Count >= 5, Count <= 500) %>%
        arrange(p.adjust)
    } else {
      data.frame()
    }
    write.csv(ego_results,  file.path(deg_dir, "SPP1pos_GO_results.csv"),   row.names = FALSE)
    write.csv(ego_filtered, file.path(deg_dir, "SPP1pos_GO_filtered.csv"),  row.names = FALSE)

    if (nrow(ego_filtered) > 0) {
      ego_df <- ego_filtered %>%
        group_by(ONTOLOGY) %>%
        slice_min(p.adjust, n = 6, with_ties = FALSE) %>%
        ungroup() %>%
        arrange(ONTOLOGY, p.adjust)
      p_go <- ggplot(
        ego_df,
        aes(
          x = -log10(p.adjust),
          y = reorder(Description, -log10(p.adjust))
        )
      ) +
        geom_col(aes(fill = ONTOLOGY), width = 0.72, show.legend = FALSE) +
        facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free") +
        scale_fill_manual(values = go_ontology_colors) +
        labs(
          title = "GO enrichment",
          x = "-log10 adjusted P",
          y = NULL
        ) +
        theme_nature(base_size = 10) +
        theme(
          # 主标题
          plot.title.position = "plot",
          plot.title = element_text(
            hjust = 0.5,
            size = 30,
            face = "bold",
            margin = margin(b = 26)
          ),

          # 坐标轴标题
          axis.title.x = element_text(
            size = 20,
            face = "bold"
          ),
          axis.title.y = element_text(
            size = 22,
            face = "bold"
          ),
          
          # 坐标轴刻度文字
          axis.text.x = element_text(
            size = 20,
            face = "bold"
          ),
          axis.text.y = element_text(
            size = 18,
            face = "bold"
          ),
          
          # 分面标题，也就是 BP / CC / MF 这类标题
          strip.text.y = element_text(
            size = 20,
            face = "bold"
          ),
          strip.text.x = element_text(
            size = 20,
            face = "bold"
          ),
          
          # 图例文字，虽然这里 show.legend = FALSE，保留也没问题
          legend.title = element_text(
            size = 13,
            face = "bold"
          ),
          legend.text = element_text(
            size = 12
          )
        )
      
      save_plot(
        p_go,
        file.path(plot_dir, "SPP1pos_GO_top_terms"),
        width = 10,
        height = 8
      )
    }

    kk <- tryCatch(
      enrichKEGG(
        gene          = cluster_gene,
        organism      = "hsa",
        pvalueCutoff  = 0.05,
        pAdjustMethod = "BH",
        qvalueCutoff  = 0.2
      ),
      error = function(e) NULL
    )
    kk_results  <- if (is.null(kk)) data.frame() else as.data.frame(kk)
    kk_filtered <- if (nrow(kk_results) > 0 && "p.adjust" %in% colnames(kk_results)) {
      kk_results %>% filter(p.adjust < 0.05) %>% arrange(p.adjust)
    } else {
      data.frame()
    }
    write.csv(kk_results,  file.path(deg_dir, "SPP1pos_KEGG_results.csv"),  row.names = FALSE)
    write.csv(kk_filtered, file.path(deg_dir, "SPP1pos_KEGG_filtered.csv"), row.names = FALSE)

    if (nrow(kk_filtered) > 0) {
      kegg_df <- kk_filtered %>%
        slice_min(p.adjust, n = 15, with_ties = FALSE) %>%
        mutate(
          Description_short = gsub(" - Homo sapiens \\(human\\)$", "", Description),
          Description_short = gsub(" pathway$", "", Description_short)
        )
      p_kegg <- ggplot(kegg_df,
                       aes(x = -log10(p.adjust),
                           y = reorder(Description_short, -log10(p.adjust)))) +
        geom_col(aes(fill = -log10(p.adjust)), width = 0.72) +
        scale_fill_gradientn(colors = kegg_gradient_colors, name = "-log10(adj.P)") +
        labs(
          title    = "Top KEGG pathways: SPP1+ macrophages",
          subtitle = "Up-regulated genes excluding SPP1",
          x = "-log10 adjusted P",
          y = NULL
        ) +
        theme_nature(base_size = 10) +
        theme(axis.text.y = element_text(size = 9))
      save_plot(p_kegg, file.path(plot_dir, "SPP1pos_KEGG_top_pathways"),
                width = 9, height = 6)
    }

    save(ego, kk, ego_results, ego_filtered, kk_results, kk_filtered,
         file = file.path(deg_dir, "SPP1pos_GO_KEGG_results.rda"))

  } else {
    empty_csv(file.path(deg_dir, "SPP1pos_GO_results.csv"))
    empty_csv(file.path(deg_dir, "SPP1pos_GO_filtered.csv"))
    empty_csv(file.path(deg_dir, "SPP1pos_KEGG_results.csv"))
    empty_csv(file.path(deg_dir, "SPP1pos_KEGG_filtered.csv"))
  }
} else {
  empty_csv(file.path(deg_dir, "SPP1pos_GO_results.csv"))
  empty_csv(file.path(deg_dir, "SPP1pos_GO_filtered.csv"))
  empty_csv(file.path(deg_dir, "SPP1pos_KEGG_results.csv"))
  empty_csv(file.path(deg_dir, "SPP1pos_KEGG_filtered.csv"))
}

## =========================================================
## 8. CellChat 分析：SPP1+/- 巨噬细胞与内皮细胞/成纤维细胞通讯
## =========================================================
run_spp1_cellchat <- function(combined_obj, macrophage_obj) {
  scRNA <- combined_obj
  scRNA$cellchat_ident <- as.character(scRNA$cell_type)

  macro_cells <- colnames(macrophage_obj)
  scRNA$cellchat_ident[macro_cells] <- as.character(macrophage_obj$SPP1_group)

  keep_cells  <- colnames(scRNA)[scRNA$cellchat_ident %in% cellchat_idents]
  scRNA_focus <- subset(scRNA, cells = keep_cells)
  scRNA_focus$cellchat_ident <- factor(scRNA_focus$cellchat_ident, levels = cellchat_idents)

  DefaultAssay(scRNA_focus) <- "RNA"
  scRNA_focus <- NormalizeData(scRNA_focus, verbose = FALSE)
  data_input  <- tryCatch(
    GetAssayData(scRNA_focus, assay = "RNA", layer = "data"),
    error = function(e) GetAssayData(scRNA_focus, assay = "RNA", slot = "data")
  )

  meta       <- scRNA_focus@meta.data
  meta$ident <- scRNA_focus$cellchat_ident
  cell.use   <- colnames(data_input)
  meta       <- meta[cell.use, ]

  focus_counts <- as.data.frame(table(meta$ident), stringsAsFactors = FALSE)
  colnames(focus_counts) <- c("ident", "cell_count")
  write.csv(focus_counts,
            file.path(cellchat_dir, "SPP1_focus_cell_counts.csv"), row.names = FALSE)

  cellchat <- createCellChat(object = data_input, meta = meta, group.by = "ident")
  cellchat <- addMeta(cellchat, meta = meta)
  cellchat <- setIdent(cellchat, ident.use = "ident")
  cellchat@idents <- droplevels(cellchat@idents)

  CellChatDB.use <- subsetDB(CellChatDB.human, search = "Secreted Signaling")
  cellchat@DB    <- CellChatDB.use

  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)

  saveRDS(cellchat, file.path(cellchat_dir, "SPP1_focus_cellchat.rds"))
  net_all <- subsetCommunication(cellchat)
  write.csv(net_all, file.path(cellchat_dir, "SPP1_focus_net_all.csv"), row.names = FALSE)

  spp1_net <- net_all %>%
    filter(
      source       %in% c("SPP1+ macrophage", "SPP1- macrophage"),
      target       %in% c("Endothelial", "Fibroblast"),
      pathway_name == "SPP1"
    )
  write.csv(spp1_net,
            file.path(cellchat_dir, "SPP1_focus_net_SPP1_pathway.csv"), row.names = FALSE)

  if (nrow(spp1_net) == 0) {
    writeLines(
      paste("No SPP1 pathway communication was detected by CellChat",
            "under the current filtering settings."),
      file.path(cellchat_dir, "SPP1_pathway_not_detected.txt")
    )
    return(list(cellchat = cellchat, spp1_net = spp1_net))
  }

  ## 配体-受体气泡图：放大点尺寸范围
  bubble_file_pdf <- file.path(cellchat_dir, "SPP1_pathway_bubble.pdf")
  bubble_file_png <- file.path(cellchat_dir, "SPP1_pathway_bubble.png")
  
  p_bubble <- netVisual_bubble(
    cellchat,
    sources.use    = c("SPP1+ macrophage", "SPP1- macrophage"),
    targets.use    = c("Endothelial", "Fibroblast"),
    signaling      = "SPP1",
    remove.isolate = FALSE
  ) +
    scale_size_continuous(range = c(4, 10)) +
    ggtitle("SPP1 signaling") +
    theme(
      ## 先统一压低整体文字
      text = element_text(size = 6, face = "bold"),
      
      ## 主标题
      plot.title = element_text(
        size = 14,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 12)
      ),
      
      ## 保持中间主图为正方形
      aspect.ratio = 1,
      
      ## x轴文字：一定要小，不然会爆炸
      axis.text.x = element_text(
        size = 6,
        angle = 45,
        hjust = 1,
        vjust = 1,
        face = "bold"
      ),
      
      ## y轴文字
      axis.text.y = element_text(
        size = 8,
        face = "bold"
      ),
      
      ## 去掉轴标题
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      
      ## 图例
      legend.position = "right",
      legend.justification = "top",
      legend.box.just = "top",
      legend.title = element_text(size = 6, face = "bold"),
      legend.text  = element_text(size = 6, face = "bold"),
      legend.key.width  = unit(0.35, "cm"),
      legend.key.height = unit(0.35, "cm"),
      
      ## 图例和主图距离
      legend.margin = margin(l = 1),
      legend.box.margin = margin(l = 2),
      
      ## 减少白色背景
      plot.margin = margin(t = 4, r = 4, b = 35, l = 4),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.2
      )
    ) +
    guides(
      size = guide_legend(
        title.theme = element_text(size = 10, face = "bold"),
        label.theme = element_text(size = 9, face = "bold")
      ),
      colour = guide_colorbar(
        title.theme = element_text(size = 10, face = "bold"),
        label.theme = element_text(size = 9, face = "bold")
      ),
      color = guide_colorbar(
        title.theme = element_text(size = 10, face = "bold"),
        label.theme = element_text(size = 9, face = "bold")
      )
    ) +
    coord_cartesian(clip = "off")
  
  ggsave(
    bubble_file_pdf,
    p_bubble,
    width = 5.2,
    height = 4.8,
    bg = "white"
  )
  
  ggsave(
    bubble_file_png,
    p_bubble,
    width = 5.2,
    height = 4.8,
    dpi = 300,
    bg = "white"
  )

  spp1_strength <- spp1_net %>%
    group_by(source, target) %>%
    summarise(
      communication_strength = sum(prob, na.rm = TRUE),
      n_interactions         = n(),
      ligand_receptor        = paste(interaction_name, collapse = ";"),
      .groups = "drop"
    )
  write.csv(spp1_strength,
            file.path(cellchat_dir,
                      "SPP1_signal_strength_to_Endothelial_Fibroblast.csv"),
            row.names = FALSE)

  ## 信号强度图：独立配色（岩板蓝 vs 古铜金）
  p_strength <- ggplot(
    spp1_strength,
    aes(x = target, y = communication_strength, fill = source)
  ) +
    geom_col(
      position  = position_dodge(width = 0.72),
      width     = 0.65,
      color     = "black",
      linewidth = 0.25
    ) +
    geom_text(
      aes(label = n_interactions),
      position = position_dodge(width = 0.72),
      vjust    = -0.35,
      size     = 5,
      fontface = "bold"
    ) +
    scale_fill_manual(values = strength_cols) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
    labs(
      title = "SPP1 signaling strength",
      x     = "Target cell type",
      y     = "Communication probability",
      fill  = NULL
    ) +
    theme_nature(base_size = 11) +
    theme(
      ## 图例色块正方形
      legend.key.width  = unit(0.45, "cm"),
      legend.key.height = unit(0.45, "cm"),
      
      ## 图例位置
      legend.position = "top",
      legend.text = element_text(size = 10, face = "bold"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.margin = margin(b = 6),
      legend.box.margin = margin(b = 4),
      
      ## 主标题
      plot.title = element_text(
        size = 20,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 10)
      ),
      
      ## 坐标轴刻度文字
      axis.text.x = element_text(size = 16, face = "bold"),
      axis.text.y = element_text(size = 16, face = "bold"),
      
      ## 坐标轴标题
      axis.title.x = element_text(size = 16, face = "bold", margin = margin(t = 8)),
      axis.title.y = element_text(size = 16, face = "bold", margin = margin(r = 8)),
      
      ## 减少外部空白背景
      plot.margin = margin(t = 8, r = 10, b = 8, l = 10)
    )
  
  save_plot(
    p_strength,
    file.path(cellchat_dir, "SPP1_signal_strength_to_Endothelial_Fibroblast"),
    width = 6,
    height = 5.5
  )

  list(cellchat = cellchat, spp1_net = spp1_net)
}


cellchat_result <- tryCatch(
  run_spp1_cellchat(combined_anno, macrophage_cells),
  error = function(e) {
    writeLines(paste("Focused CellChat failed:", conditionMessage(e)),
               file.path(cellchat_dir, "CellChat_error.txt"))
    empty_csv(file.path(cellchat_dir, "SPP1_focus_net_all.csv"))
    empty_csv(file.path(cellchat_dir, "SPP1_focus_net_SPP1_pathway.csv"))
    NULL
  }
)

## =========================================================
## 审稿补充分析：SPP1-high定义稳健性与样本层面统计
## =========================================================

review_spp1_dir <- "results/reviewer_supplement/SPP1_high_robustness"
dir.create(review_spp1_dir, recursive = TRUE, showWarnings = FALSE)

if (!exists("combined_anno")) {
  load(seurat_file)
  DefaultAssay(combined_anno) <- "SCT"
}
if (!exists("macrophage_cells")) {
  macrophage_cells <- subset(combined_anno, subset = cell_type %in% macrophage_types)
  macrophage_cells$cell_type <- factor(macrophage_cells$cell_type, levels = macrophage_types)
}

review_spp1_expr <- FetchData(macrophage_cells, vars = "SPP1")[, 1]
review_spp1_median <- median(review_spp1_expr, na.rm = TRUE)
review_spp1_q75 <- unname(quantile(review_spp1_expr, probs = 0.75, na.rm = TRUE))
review_spp1_fixed <- 1

review_group_df <- data.frame(
  cell = colnames(macrophage_cells),
  sample = as.character(macrophage_cells$sample),
  cell_type = as.character(macrophage_cells$cell_type),
  SPP1 = review_spp1_expr,
  median_split = ifelse(review_spp1_expr > review_spp1_median, "SPP1-high", "SPP1-low"),
  upper_quartile = ifelse(review_spp1_expr > review_spp1_q75, "SPP1-high", "SPP1-low"),
  expression_positive = ifelse(review_spp1_expr > 0, "SPP1-high", "SPP1-low"),
  fixed_cutoff_1 = ifelse(review_spp1_expr > review_spp1_fixed, "SPP1-high", "SPP1-low"),
  stringsAsFactors = FALSE
)

cohen_kappa_binary <- function(reference, alternative) {
  reference <- factor(reference, levels = c("SPP1-low", "SPP1-high"))
  alternative <- factor(alternative, levels = c("SPP1-low", "SPP1-high"))
  tab <- table(reference, alternative)
  n <- sum(tab)
  po <- sum(diag(tab)) / n
  pe <- sum(rowSums(tab) * colSums(tab)) / (n * n)
  if (isTRUE(all.equal(1, pe))) return(NA_real_)
  (po - pe) / (1 - pe)
}

review_classification_agreement <- lapply(
  c("upper_quartile", "expression_positive", "fixed_cutoff_1"),
  function(method) {
    tab <- table(review_group_df$median_split, review_group_df[[method]])
    data.frame(
      reference = "median_split",
      alternative = method,
      cells = nrow(review_group_df),
      median_high_cells = sum(review_group_df$median_split == "SPP1-high"),
      alternative_high_cells = sum(review_group_df[[method]] == "SPP1-high"),
      overall_agreement = mean(review_group_df$median_split == review_group_df[[method]]),
      cohens_kappa = cohen_kappa_binary(review_group_df$median_split, review_group_df[[method]]),
      stringsAsFactors = FALSE
    )
  }
) %>% bind_rows()

write.csv(review_group_df, file.path(review_spp1_dir, "SPP1_high_alternative_definitions.csv"), row.names = FALSE)
write.csv(review_classification_agreement, file.path(review_spp1_dir, "SPP1_high_definition_agreement.csv"), row.names = FALSE)

review_module_scores <- read.csv(file.path(out_dir, "SPP1_positive_module_scores.csv"), check.names = FALSE)
review_module_names <- setdiff(colnames(review_module_scores), c("cell", "SPP1_group", "cell_type", "sample"))
review_module_scores$SPP1_group <- factor(
  review_module_scores$SPP1_group,
  levels = c("SPP1- macrophage", "SPP1+ macrophage")
)

review_module_long <- review_module_scores %>%
  tidyr::pivot_longer(
    cols = all_of(review_module_names),
    names_to = "module",
    values_to = "score"
  )

review_sample_module <- review_module_long %>%
  group_by(sample, SPP1_group, module) %>%
  summarise(
    mean_score = mean(score, na.rm = TRUE),
    cells = sum(is.finite(score)),
    .groups = "drop"
  )

review_sample_module_wide <- review_sample_module %>%
  tidyr::pivot_wider(
    names_from = SPP1_group,
    values_from = c(mean_score, cells),
    names_sep = "__"
  ) %>%
  mutate(
    delta_SPP1pos_minus_SPP1neg =
      `mean_score__SPP1+ macrophage` - `mean_score__SPP1- macrophage`
  )

review_sample_module_tests <- review_sample_module_wide %>%
  group_by(module) %>%
  summarise(
    samples = sum(is.finite(delta_SPP1pos_minus_SPP1neg)),
    mean_delta = mean(delta_SPP1pos_minus_SPP1neg, na.rm = TRUE),
    median_delta = median(delta_SPP1pos_minus_SPP1neg, na.rm = TRUE),
    positive_direction_samples = sum(delta_SPP1pos_minus_SPP1neg > 0, na.rm = TRUE),
    paired_wilcox_p = ifelse(
      samples >= 3 && length(unique(delta_SPP1pos_minus_SPP1neg[is.finite(delta_SPP1pos_minus_SPP1neg)])) > 1,
      wilcox.test(delta_SPP1pos_minus_SPP1neg, mu = 0, exact = FALSE)$p.value,
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(paired_wilcox_p, method = "BH"))

write.csv(review_sample_module, file.path(review_spp1_dir, "module_scores_sample_group_means.csv"), row.names = FALSE)
write.csv(review_sample_module_wide, file.path(review_spp1_dir, "module_scores_sample_group_deltas.csv"), row.names = FALSE)
write.csv(review_sample_module_tests, file.path(review_spp1_dir, "module_scores_sample_level_tests.csv"), row.names = FALSE)

review_lmm_results <- lapply(review_module_names, function(module_name) {
  df <- review_module_long %>%
    filter(module == module_name, is.finite(score)) %>%
    mutate(SPP1_group = relevel(factor(SPP1_group), ref = "SPP1- macrophage"))
  if (!requireNamespace("lme4", quietly = TRUE) || length(unique(df$sample)) < 2) {
    return(data.frame(module = module_name, lrt_p = NA_real_, estimate = NA_real_))
  }
  full <- lme4::lmer(score ~ SPP1_group + (1 | sample), data = df, REML = FALSE)
  null <- lme4::lmer(score ~ 1 + (1 | sample), data = df, REML = FALSE)
  an <- anova(null, full)
  estimate <- unname(lme4::fixef(full)["SPP1_groupSPP1+ macrophage"])
  data.frame(module = module_name, lrt_p = an$`Pr(>Chisq)`[2], estimate = estimate)
}) %>% bind_rows() %>%
  mutate(p_adj = p.adjust(lrt_p, method = "BH"))

write.csv(review_lmm_results, file.path(review_spp1_dir, "module_scores_mixed_model_lme4.csv"), row.names = FALSE)

p_review_sample_module <- ggplot(
  review_sample_module,
  aes(SPP1_group, mean_score, group = sample, colour = sample)
) +
  geom_line(alpha = 0.75, linewidth = 0.55) +
  geom_point(size = 1.8) +
  facet_wrap(~module, scales = "free_y", ncol = 4) +
  labs(
    title = "Sample-level module-score directionality",
    x = NULL,
    y = "Sample mean module score"
  ) +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

ggsave(file.path(review_spp1_dir, "module_scores_sample_level_paired.pdf"), p_review_sample_module,
       width = 11, height = 7, bg = "white")
ggsave(file.path(review_spp1_dir, "module_scores_sample_level_paired.png"), p_review_sample_module,
       width = 11, height = 7, dpi = 300, bg = "white")

writeLines(
  c(
    "Reviewer supplement: SPP1-high robustness and sample-level statistics",
    paste0("cells=", nrow(review_group_df)),
    paste0("median_cutoff=", signif(review_spp1_median, 5)),
    paste0("upper_quartile_cutoff=", signif(review_spp1_q75, 5)),
    paste0("fixed_cutoff=", review_spp1_fixed),
    "Mixed models used score ~ SPP1_group + (1|sample)."
  ),
  file.path(review_spp1_dir, "SPP1_high_robustness_summary.txt")
)
