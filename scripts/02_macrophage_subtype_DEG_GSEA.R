load("rda/UMAP_annotated.rda")

# 加载依赖包
library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(gridExtra)
library(scales)

# ==============================================================================
# 1. 数据准备：提取 4 类巨噬细胞
# ==============================================================================

macrophage_types <- c("LAM", "IAM", "OSM", "TRM")
out_base_dir <- "results/macrophage_DEG_GSEA"
dir.create(out_base_dir, recursive = TRUE, showWarnings = FALSE)

required_annotation_columns <- c("cell_type", "cell_type_subtype")
missing_annotation_columns <- setdiff(required_annotation_columns, colnames(combined_anno@meta.data))
if (length(missing_annotation_columns) > 0) {
  stop(
    "Missing annotation columns in rda/UMAP_annotated.rda: ",
    paste(missing_annotation_columns, collapse = ", "),
    ". Re-run 01_DR_single_cell_processing_annotation.R first."
  )
}

macrophage_cells <- subset(
  combined_anno,
  subset = cell_type_subtype %in% macrophage_types
)
DefaultAssay(macrophage_cells) <- "SCT"
macrophage_cells$cell_type_subtype <- factor(
  macrophage_cells$cell_type_subtype,
  levels = macrophage_types
)
Idents(macrophage_cells) <- "cell_type_subtype"

subtype_cell_counts <- table(macrophage_cells$cell_type_subtype)
print(subtype_cell_counts)
if (any(subtype_cell_counts == 0)) {
  stop(
    "At least one macrophage subtype has no cells: ",
    paste(names(subtype_cell_counts)[subtype_cell_counts == 0], collapse = ", ")
  )
}

# ==============================================================================
# 2. Nature-style theme and palettes
# ==============================================================================

# Nature 风格主题
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

# 在 Nature 风格基础上放大标题、坐标轴和图例，生成更适合发表图的主题。
theme_publication <- function(base_size = 12) {
  theme_nature(base_size = base_size) +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5, margin = margin(b = 14)),
      plot.subtitle = element_text(size = 11, color = "grey35", hjust = 0.5, margin = margin(b = 10)),
      axis.text = element_text(size = 14, face = "bold", color = "black"),
      axis.title = element_text(size = 16, face = "bold", color = "black"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 12, face = "bold", color = "black"),
      legend.key.width = unit(0.45, "cm"),
      legend.key.height = unit(0.45, "cm"),
      strip.text = element_text(size = 14, face = "bold", color = "black")
    )
}

# 保存 
save_ggplot_white <- function(filename, plot, width, height, dpi = 300, device = NULL) {
  ggsave(
    filename,
    plot,
    width = width,
    height = height,
    dpi = dpi,
    device = device,
    bg = "white"
  )
}

pal_deg <- c(
  "Up"      = "#B2182B",
  "Down"    = "#2166AC",
  "Not_Sig" = "#D6CEC3"
)

pal_heat <- c(
  "#2166AC",  # 深蓝
  "#67A9CF",  # 亮蓝
  "#F7F7F7",  # 近白
  "#EF8A62",  # 亮红
  "#B2182B"   # 深红
)

# 为巨噬细胞亚型设置配色
subtype_pal <- c(
  "LAM" = "#2F6F73",
  "IAM" = "#B65C4A",
  "OSM" = "#6E5F8F",
  "TRM" = "#C49A3A"
)

# ==============================================================================
# 3. 巨噬细胞亚群注释前后 tSNE 与 marker 气泡图
# ==============================================================================

macrophage_subtype_marker_panel <- list(
  "LAM" = c("LGALS3", "ACP5", "PLA2G7", "GPNMB", "CTSD"),
  "IAM" = c("CCL4L2", "CCL3L1", "CCL3", "CCL4", "CD83"),
  "OSM" = c("TXNRD1", "SCD", "GCLM", "SQSTM1", "CTSB"),
  "TRM" = c("STAB1", "MAF", "DAB2", "F13A1", "MRC1", "CD163")
)
macrophage_subtype_marker_panel <- lapply(
  macrophage_subtype_marker_panel,
  function(x) unique(x[x %in% rownames(macrophage_cells)])
)
macrophage_subtype_marker_panel <- macrophage_subtype_marker_panel[
  lengths(macrophage_subtype_marker_panel) > 0
]

plot_macrophage_annotation_overview <- function(
    seurat_obj,
    marker_panel,
    subtype_pal,
    out_base_dir) {
  annotation_out_dir <- file.path(out_base_dir, "annotation_visualization")
  dir.create(annotation_out_dir, recursive = TRUE, showWarnings = FALSE)

  seurat_obj$cluster_before_annotation <- factor(
    as.character(seurat_obj$seurat_clusters),
    levels = c("0", "1", "2", "3")
  )
  seurat_obj$cell_type_subtype <- factor(
    seurat_obj$cell_type_subtype,
    levels = macrophage_types
  )

  cluster_pal <- setNames(unname(subtype_pal[macrophage_types]), c("0", "1", "2", "3"))
  dotplot_gradient <- c("#2F6F9F", "#8DB9D3", "#F6F0E8", "#D96B5A")

  tsne_theme <- theme_publication(base_size = 12) +
    theme(
      plot.title = element_text(size = 28, face = "bold", hjust = 0.5,
                                margin = margin(b = 14)),
      plot.subtitle = element_text(size = 16, face = "bold", hjust = 0.5,
                                   margin = margin(b = 10)),
      axis.text = element_text(size = 18, face = "bold", color = "black"),
      axis.title = element_text(size = 22, face = "bold", color = "black"),
      legend.title = element_text(size = 20, face = "bold"),
      legend.text = element_text(size = 18, face = "bold"),
      plot.margin = margin(16, 16, 16, 16)
    )

  p_tsne_before <- DimPlot(
    seurat_obj,
    reduction = "tsne",
    group.by = "cluster_before_annotation",
    cols = cluster_pal,
    label = TRUE,
    label.size = 7,
    repel = TRUE,
    pt.size = 2,
    alpha = 1
  ) +
    labs(
      title = "Macrophage clusters (pre-annotation)",
      subtitle = paste0(format(ncol(seurat_obj), big.mark = ","), " cells  ·  tSNE"),
      x = "tSNE 1", y = "tSNE 2", color = "Cluster"
    ) +
    tsne_theme

  p_tsne_after <- DimPlot(
    seurat_obj,
    reduction = "tsne",
    group.by = "cell_type_subtype",
    cols = subtype_pal,
    label = TRUE,
    label.size = 7,
    repel = TRUE,
    pt.size = 2,
    alpha = 1
  ) +
    labs(
      title = "Macrophage subtype annotation",
      subtitle = paste0(format(ncol(seurat_obj), big.mark = ","), " cells  ·  tSNE"),
      x = "tSNE 1", y = "tSNE 2", color = "Subtype"
    ) +
    tsne_theme

  dotplot_theme_before <- theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(size = 26, face = "bold", hjust = 0.5,
                                margin = margin(b = 14)),
      strip.background = element_rect(fill = "#F4F4F4", color = "grey55",
                                      linewidth = 0.35),
      strip.text.x = element_text(size = 15, face = "bold", color = "black"),
      axis.text.x = element_text(size = 13, face = "bold.italic", color = "black",
                                 angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 15, face = "bold", color = "black"),
      axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 10)),
      axis.title.x = element_blank(),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey35", fill = NA, linewidth = 0.4),
      panel.spacing.x = unit(2.5, "mm"),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(size = 15, face = "bold"),
      legend.text = element_text(size = 13, face = "bold"),
      plot.margin = margin(14, 14, 14, 14)
    )

  dotplot_theme_after <- dotplot_theme_before +
    theme(
      plot.title = element_text(size = 28, face = "bold", hjust = 0.5,
                                margin = margin(b = 16)),
      strip.text.x = element_text(size = 16, face = "bold", color = "black"),
      axis.text.y = element_text(size = 16, face = "bold", color = "black")
    )

  p_dot_before <- DotPlot(
    seurat_obj,
    features = marker_panel,
    group.by = "cluster_before_annotation",
    dot.scale = 7
  ) +
    scale_color_gradientn(
      colours = dotplot_gradient,
      name = "Scaled average\nexpression",
      limits = c(-2, 2),
      oob = squish
    ) +
    scale_size_continuous(range = c(0.4, 7), name = "Percent\nexpressed") +
    labs(title = "Macrophage marker expression (pre-annotation)", y = "Cluster") +
    dotplot_theme_before +
    guides(
      color = guide_colorbar(title.position = "top", title.hjust = 0.5,
                             barwidth = unit(42, "mm"), barheight = unit(4, "mm")),
      size = guide_legend(title.position = "top", title.hjust = 0.5,
                          nrow = 1, byrow = TRUE)
    )

  p_dot_after <- DotPlot(
    seurat_obj,
    features = marker_panel,
    group.by = "cell_type_subtype",
    dot.scale = 7
  ) +
    scale_color_gradientn(
      colours = dotplot_gradient,
      name = "Scaled average\nexpression",
      limits = c(-2, 2),
      oob = squish
    ) +
    scale_size_continuous(range = c(0.4, 7), name = "Percent\nexpressed") +
    labs(title = "Marker gene expression by macrophage subtype", y = "Subtype") +
    dotplot_theme_after +
    guides(
      color = guide_colorbar(title.position = "top", title.hjust = 0.5,
                             barwidth = unit(42, "mm"), barheight = unit(4, "mm")),
      size = guide_legend(title.position = "top", title.hjust = 0.5,
                          nrow = 1, byrow = TRUE)
    )

  plots <- list(
    Macrophage_clusters_before_annotation_tsne = p_tsne_before,
    Macrophage_subtypes_after_annotation_tsne = p_tsne_after,
    Macrophage_marker_before_annotation = p_dot_before,
    Macrophage_marker_after_annotation = p_dot_after
  )
  plot_sizes <- list(
    Macrophage_clusters_before_annotation_tsne = c(10, 8),
    Macrophage_subtypes_after_annotation_tsne = c(10, 8),
    Macrophage_marker_before_annotation = c(14, 9),
    Macrophage_marker_after_annotation = c(14, 9)
  )

  for (plot_name in names(plots)) {
    plot_size <- plot_sizes[[plot_name]]
    save_ggplot_white(
      file.path(annotation_out_dir, paste0(plot_name, ".pdf")),
      plots[[plot_name]], width = plot_size[1], height = plot_size[2],
      dpi = 500, device = "pdf"
    )
    save_ggplot_white(
      file.path(annotation_out_dir, paste0(plot_name, ".png")),
      plots[[plot_name]], width = plot_size[1], height = plot_size[2],
      dpi = 300
    )
  }

  invisible(plots)
}

macrophage_annotation_plots <- plot_macrophage_annotation_overview(
  seurat_obj = macrophage_cells,
  marker_panel = macrophage_subtype_marker_panel,
  subtype_pal = subtype_pal,
  out_base_dir = out_base_dir
)

# 根据某个亚型主色生成由浅到深的渐变色，用于区分不同样本
make_subtype_shades <- function(hex, n) {
  if (is.null(hex) || is.na(hex)) {
    hex <- "#4D4D4D"
  }
  if (n <= 1) {
    return(hex)
  }
  grDevices::colorRampPalette(c("#F3F0E8", hex))(n + 1)[-1]
}

# ==============================================================================
# 4. 函数：绘制各样本中的巨噬细胞亚型柱状图
# ==============================================================================

# 统计每个 sample 中四类巨噬细胞亚型的细胞数和占比，并按亚型分别输出柱状图。
plot_macrophage_subtype_abundance_by_sample <- function(seurat_obj, macrophage_types, out_base_dir) {
  sample_out_dir <- file.path(out_base_dir, "subtype_abundance")
  dir.create(sample_out_dir, recursive = TRUE, showWarnings = FALSE)

  sample_levels <- sort(unique(seurat_obj$sample))

  sample_total_counts <- as.data.frame(
    table(seurat_obj$sample),
    stringsAsFactors = FALSE
  )
  colnames(sample_total_counts) <- c("sample", "sample_macrophage_total")

  subtype_sample_counts <- as.data.frame(
    table(seurat_obj$sample, seurat_obj$cell_type_subtype),
    stringsAsFactors = FALSE
  )
  colnames(subtype_sample_counts) <- c("sample", "cell_type_subtype", "cell_count")

  subtype_sample_counts <- subtype_sample_counts %>%
    filter(cell_type_subtype %in% macrophage_types) %>%
    left_join(sample_total_counts, by = "sample") %>%
    mutate(
      sample = factor(sample, levels = sample_levels),
      cell_type_subtype = factor(cell_type_subtype, levels = macrophage_types),
      proportion = cell_count / sample_macrophage_total,
      percent = proportion * 100,
      label = paste0(cell_count, "\n", sprintf("%.1f%%", percent))
    ) %>%
    arrange(cell_type_subtype, sample)

  sample_plots <- lapply(macrophage_types, function(subtype) {
    plot_df <- subtype_sample_counts %>%
      filter(cell_type_subtype == subtype)

    sample_y_max <- max(plot_df$cell_count) * 1.22
    if (sample_y_max == 0) {
      sample_y_max <- 1
    }

    sample_pal <- setNames(
      make_subtype_shades(subtype_pal[[subtype]], nrow(plot_df)),
      as.character(plot_df$sample)
    )

    p_sample <- ggplot(plot_df, aes(x = sample, y = cell_count, fill = sample)) +
      geom_col(
        width = 0.68,
        color = "#2B2B2B",
        linewidth = 0.28
      ) +
      geom_text(
        aes(label = label),
        vjust = -0.35,
        size = 4.6,
        lineheight = 0.9,
        color = "#1F1F1F",
        fontface = "bold"
      ) +
      scale_fill_manual(values = sample_pal, guide = "none") +
      scale_y_continuous(
        limits = c(0, sample_y_max),
        expand = expansion(mult = c(0, 0.04))
      ) +
      labs(
        title = paste0(subtype, " macrophages across samples"),
        x = NULL,
        y = "Cell count"
      ) +
      theme_publication(base_size = 12) +
      theme(
        aspect.ratio = 1,
        axis.text.x = element_text(
          angle = 40,
          hjust = 1,
          face = "bold",
          size = 15,
          color = "black"
        ),
        axis.text.y = element_text(
          face = "bold",
          size = 15,
          color = "black"
        ),
        axis.title.y = element_text(
          size = 17,
          face = "bold"
        ),
        plot.title = element_text(
          size = 24,
          face = "bold",
          hjust = 0.5,
          margin = margin(b = 18)
        ),
        panel.border = element_rect(
          color = "black",
          fill = NA,
          linewidth = 0.45
        )
      )

    save_ggplot_white(
      file.path(sample_out_dir, paste0("Macrophage_subtype_abundance_by_sample_", subtype, ".pdf")),
      p_sample,
      width = 6.5,
      height = 6.5,
      device = "pdf",
      dpi = 500
    )

    save_ggplot_white(
      file.path(sample_out_dir, paste0("Macrophage_subtype_abundance_by_sample_", subtype, ".png")),
      p_sample,
      width = 6.5,
      height = 6.5,
      dpi = 300
    )

    p_sample
  })

  names(sample_plots) <- macrophage_types
  sample_plots
}

# ==============================================================================
# 5. 函数：差异分析、GSEA、火山图、热图、GSEA通路图
# ==============================================================================

# 指定巨噬细胞亚型与其他亚型：差异分析
run_deg <- function(seurat_obj, ident_name, assay = "SCT",
                    min.pct = 0.25, logfc.threshold = 0.25,
                    test.use = "wilcox") {
  degs <- FindMarkers(
    object = seurat_obj,
    ident.1 = ident_name,
    ident.2 = setdiff(macrophage_types, ident_name),
    assay = assay,
    recorrect_umi = FALSE,
    min.pct = min.pct,
    logfc.threshold = logfc.threshold,
    test.use = test.use
  )
  
  degs$gene <- rownames(degs)
  degs
}

# 标记 Up、Down 或 Not_Sig
annotate_deg <- function(deg_df, logfc_cutoff = 0.25, p_val_cutoff = 0.05) {
  deg_df %>%
    mutate(
      change = case_when(
        p_val_adj < p_val_cutoff & avg_log2FC > logfc_cutoff  ~ "Up",
        p_val_adj < p_val_cutoff & avg_log2FC < -logfc_cutoff ~ "Down",
        TRUE ~ "Not_Sig"
      ),
      neg_log10_padj = -log10(pmax(p_val_adj, .Machine$double.xmin)),
      change = factor(change, levels = c("Down", "Not_Sig", "Up"))
    )
}

# 按 avg_log2FC 构建排序基因列表，并运行 GSEA 富集分析
run_gsea <- function(deg_df) {
  gsea_df <- deg_df %>%
    arrange(desc(avg_log2FC)) %>%
    filter(!is.na(avg_log2FC)) %>%
    distinct(gene, .keep_all = TRUE)
  
  gene_list_gsea <- gsea_df$avg_log2FC
  names(gene_list_gsea) <- gsea_df$gene
  
  set.seed(1234)
  
  gsea_res <- gseGO(
    geneList       = gene_list_gsea,
    OrgDb          = org.Hs.eg.db,
    keyType        = "SYMBOL",
    ont            = "ALL",
    minGSSize      = 10,
    maxGSSize      = 500,
    pvalueCutoff   = 0.05,
    pAdjustMethod  = "BH",
    verbose        = FALSE
  )
  
  list(
    gsea_res = gsea_res,
    gsea_res_df = as.data.frame(gsea_res)
  )
}

# 火山图
plot_volcano <- function(deg_df, comparison_name, out_dir,
                         logfc_cutoff = 0.25, p_val_cutoff = 0.05,
                         extra_label_genes = character(0)) {
  degs_plot <- annotate_deg(deg_df, logfc_cutoff, p_val_cutoff) %>%
    arrange(change)
  
  top_label_genes <- degs_plot %>%
    filter(change != "Not_Sig") %>%
    arrange(p_val_adj, desc(abs(avg_log2FC))) %>%
    slice_head(n = 16)

  extra_label_genes <- unique(extra_label_genes)
  label_genes <- bind_rows(
    top_label_genes,
    degs_plot %>%
      filter(gene %in% extra_label_genes)
  ) %>%
    distinct(gene, .keep_all = TRUE)
  
  x_lim <- max(abs(degs_plot$avg_log2FC), na.rm = TRUE)
  x_lim <- ceiling(x_lim * 10) / 10
  
  p_volcano <-
    ggplot(degs_plot, aes(x = avg_log2FC, y = neg_log10_padj)) +
    geom_point(
      aes(color = change),
      size = 4.8,
      alpha = 0.85,
      stroke = 0
    ) +
    geom_vline(
      xintercept = c(-logfc_cutoff, logfc_cutoff),
      linetype = "dashed",
      linewidth = 0.35,
      color = "grey65"
    ) +
    geom_hline(
      yintercept = -log10(p_val_cutoff),
      linetype = "dashed",
      linewidth = 0.35,
      color = "grey65"
    ) +
    geom_text_repel(
      data = label_genes,
      aes(label = gene),
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
      values = pal_deg,
      breaks = c("Up", "Down", "Not_Sig"),
      labels = c("Up", "Down", "Not significant"),
      name = NULL,
      drop = FALSE
    ) +
    coord_cartesian(xlim = c(-x_lim, x_lim), clip = "off") +
    labs(
      title = paste0("Volcano plot (", gsub("_", " ", comparison_name), ")"),
      x = "log2 fold change",
      y = "-log10 adjusted P"
      ) +
    theme_publication(base_size = 11) +
    theme(
      legend.position = "top",
      legend.justification = "center",
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5, margin = margin(b = 14)),
      axis.text.x = element_text(size = 14, face = "bold", color = "black"),
      axis.text.y = element_text(size = 14, face = "bold", color = "black"),
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 12, face = "bold"),
      plot.margin = margin(10, 16, 10, 10)
    )
  
  save_ggplot_white(
    file.path(out_dir, paste0(comparison_name, "_volcano.pdf")),
    p_volcano,
    width = 6.2,
    height = 6.2,
    device = "pdf",
    dpi = 500
  )
  
  save_ggplot_white(
    file.path(out_dir, paste0(comparison_name, "_volcano.png")),
    p_volcano,
    width = 6.2,
    height = 6.2,
    dpi = 300
  )
  
  p_volcano
}

# 热图
plot_deg_heatmap <- function(seurat_obj, deg_df, ident_name, out_dir,
                             top_n_heatmap = 20,
                             max_cells_per_group = 100) {
  comparison_name <- paste0(ident_name, "_vs_others")
  
  degs_plot <- annotate_deg(deg_df)
  
  heatmap_genes <- degs_plot %>%
    filter(change %in% c("Up", "Down")) %>%
    group_by(change) %>%
    arrange(p_val_adj, desc(abs(avg_log2FC)), .by_group = TRUE) %>%
    slice_head(n = top_n_heatmap) %>%
    ungroup() %>%
    arrange(change, desc(avg_log2FC)) %>%
    pull(gene) %>%
    unique()
  
  heatmap_genes <- intersect(heatmap_genes, rownames(seurat_obj))
  
  if (length(heatmap_genes) <= 1) {
    return(NULL)
  }
  
  heatmap_cells <- subset(
    seurat_obj,
    subset = cell_type_subtype %in% macrophage_types
  )
  
  DefaultAssay(heatmap_cells) <- "SCT"
  
  heatmap_cells$comparison_group_raw <- ifelse(
    heatmap_cells$cell_type_subtype == ident_name,
    ident_name,
    "Other macrophages"
  )
  
  group_order_raw <- c(ident_name, "Other macrophages")
  group_order_plot <- c(ident_name, "Other macrophages")
  
  heatmap_cells$comparison_group_raw <- factor(
    heatmap_cells$comparison_group_raw,
    levels = group_order_raw
  )
  
  heatmap_cells$comparison_group <- factor(
    as.character(heatmap_cells$comparison_group_raw),
    levels = group_order_plot
  )
  
  set.seed(1234)
  heat_cells <- unlist(lapply(group_order_raw, function(group_name) {
    cells <- colnames(heatmap_cells)[as.character(heatmap_cells$comparison_group_raw) == group_name]
    if (length(cells) == 0) {
      return(character(0))
    }
    if (is.infinite(max_cells_per_group)) {
      return(cells)
    }
    sample(cells, size = min(max_cells_per_group, length(cells)))
  }), use.names = FALSE)
  
  heat_cells <- heat_cells[!is.na(heat_cells)]
  
  if (length(heat_cells) == 0) {
    return(NULL)
  }
  
  heatmap_cells <- heatmap_cells[, heat_cells]
  
  heat_df <- tryCatch(
    FetchData(heatmap_cells, vars = heatmap_genes, layer = "data"),
    error = function(e1) {
      tryCatch(
        FetchData(heatmap_cells, vars = heatmap_genes, slot = "data"),
        error = function(e2) FetchData(heatmap_cells, vars = heatmap_genes)
      )
    }
  )
  
  heat_df <- heat_df[, heatmap_genes, drop = FALSE]
  heat_scaled <- scale(as.matrix(heat_df))
  heat_scaled[is.na(heat_scaled)] <- 0
  heat_scaled <- pmax(pmin(heat_scaled, 2), -2)
  
  cell_info <- data.frame(
    cell = colnames(heatmap_cells),
    comparison_group = as.character(heatmap_cells$comparison_group),
    stringsAsFactors = FALSE
  ) %>%
    group_by(comparison_group) %>%
    mutate(cell_index = row_number()) %>%
    ungroup()
  
  n_genes <- length(heatmap_genes)
  
  heat_plot_long <- as.data.frame(as.table(heat_scaled), stringsAsFactors = FALSE)
  colnames(heat_plot_long) <- c("cell", "gene", "scaled_expression")
  heat_plot_long$scaled_expression <- as.numeric(heat_plot_long$scaled_expression)
  
  heat_plot_long <- heat_plot_long %>%
    left_join(cell_info, by = "cell") %>%
    mutate(
      comparison_group = factor(comparison_group, levels = group_order_plot),
      gene = factor(gene, levels = heatmap_genes),
      gene_index = n_genes - match(as.character(gene), heatmap_genes) + 1
    )
  
  legend_df <- cell_info %>%
    group_by(comparison_group) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    mutate(
      comparison_group = factor(comparison_group, levels = group_order_plot),
      gene_index = n_genes
    )
  
  group_pal <- c(
    "IAM" = "#ED716B",
    "LAM" = "#ED716B",
    "OSM" = "#ED716B",
    "TRM" = "#ED716B",
    "Other macrophages" = "#20B8C0"
  )
  group_pal <- group_pal[group_order_plot]
  
  bar_df <- cell_info %>%
    group_by(comparison_group) %>%
    summarise(
      xmin = min(cell_index) - 0.5,
      xmax = max(cell_index) + 0.5,
      .groups = "drop"
    ) %>%
    mutate(
      comparison_group = factor(comparison_group, levels = group_order_plot),
      ymin = n_genes + 0.62,
      ymax = n_genes + 1.08,
      fill_col = unname(group_pal[as.character(comparison_group)])
    )
  
  p_heatmap <- ggplot(
    heat_plot_long,
    aes(x = cell_index, y = gene_index, fill = scaled_expression)
  ) +
    geom_raster() +
    geom_rect(
      data = bar_df,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = bar_df$fill_col,
      color = NA
    ) +
    facet_grid(
      . ~ comparison_group,
      scales = "free_x",
      space = "free_x"
    ) +
    geom_point(
      data = legend_df,
      aes(x = cell_index, y = gene_index, color = comparison_group),
      inherit.aes = FALSE,
      alpha = 0,
      size = 4,
      show.legend = TRUE
    ) +
    scale_y_continuous(
      breaks = seq_len(n_genes),
      labels = rev(heatmap_genes),
      limits = c(0.5, n_genes + 1.12),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0))) +
    scale_fill_gradientn(
      colors = pal_heat,
      limits = c(-2, 2),
      oob = squish,
      breaks = c(-2, 0, 2),
      labels = c("-2", "0", "2"),
      name = "Scaled\nexpression",
      guide = guide_colorbar(
        barheight = unit(3.8, "cm"),
        barwidth = unit(0.35, "cm"),
        title.position = "top",
        order = 1
      )
    ) +
    scale_color_manual(
      values = group_pal,
      breaks = group_order_plot,
      labels = c(ident_name, "Other\nmacrophages"),
      name = "Identity",
      guide = guide_legend(
        override.aes = list(alpha = 1, size = 4),
        order = 2
      )
    ) +
    labs(
      title = paste0("Heatmap (", gsub("_", " ", comparison_name), ")"),
      x = NULL,
      y = NULL
    ) +
    theme_nature(base_size = 9) +
    theme(
      legend.position = "right",
      legend.justification = "top",
      legend.box.just = "top",
      legend.box = "vertical",
      legend.spacing.y = unit(0.18, "cm"),
      legend.box.margin = margin(t = 34, r = 0, b = 0, l = 8),
      
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.border = element_blank(),
      panel.spacing.x = unit(0, "pt"),
      panel.grid = element_blank(),
      
      plot.title = element_text(
        size = 28,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 14),
        color = "black"
      ),
      
      axis.text.y = element_text(
        size = 14,
        face = "bold",
        color = "black",
        hjust = 1,
        margin = margin(r = 4),
        lineheight = 0.82
      ),
      
      strip.text.x = element_text(
        size = 18,
        face = "bold",
        color = "black",
        hjust = 0.5,
        lineheight = 0.82,
        margin = margin(b = 8)
      ),
      strip.background = element_blank(),
      
      legend.title = element_text(
        size = 16,
        face = "bold",
        lineheight = 0.9,
        margin = margin(b = 12),
        color = "black"
      ),
      
      legend.text = element_text(
        size = 16,
        face = "bold",
        color = "black",
        lineheight = 0.85
      ),
      
      plot.margin = margin(8, 8, 8, 8)
    )
  
  save_ggplot_white(
    file.path(out_dir, paste0(comparison_name, "_top_DEG_heatmap.pdf")),
    p_heatmap,
    width = 10,
    height = 8,
    device = "pdf",
    dpi = 500
  )
  
  save_ggplot_white(
    file.path(out_dir, paste0(comparison_name, "_top_DEG_heatmap.png")),
    p_heatmap,
    width = 10,
    height = 8,
    dpi = 300
  )
  
  p_heatmap
}

# 把校正 P 值格式化为科学计数法字符串，方便在 GSEA 图中标注
format_p_adjust <- function(x) {
  ifelse(is.na(x), "NA", formatC(x, format = "e", digits = 3))
}

# 从 GSEA 结果对象中提取指定 GO 条目的基因集合。
get_gsea_gene_set <- function(gsea_res, go_id) {
  gene_sets <- tryCatch(gsea_res@geneSets, error = function(e) NULL)
  if (!is.null(gene_sets) && go_id %in% names(gene_sets)) {
    return(gene_sets[[go_id]])
  }

  result_df <- as.data.frame(gsea_res)
  if ("core_enrichment" %in% colnames(result_df) && go_id %in% result_df$ID) {
    return(strsplit(result_df$core_enrichment[result_df$ID == go_id][1], "/", fixed = TRUE)[[1]])
  }

  character(0)
}

# 为指定 GO 条目计算运行富集分数、命中位置和排序指标等绘图数据。
build_gsea_plot_data <- function(gsea_res, go_id, exponent = 1) {
  gene_list <- gsea_res@geneList
  gene_list <- sort(gene_list[!is.na(gene_list)], decreasing = TRUE)
  genes <- names(gene_list)
  gene_set <- intersect(get_gsea_gene_set(gsea_res, go_id), genes)

  if (length(gene_set) == 0) {
    stop("No overlap between GSEA gene set and ranked gene list: ", go_id)
  }

  hit_index <- genes %in% gene_set
  n_genes <- length(gene_list)
  n_hits <- sum(hit_index)
  hit_weights <- abs(gene_list[hit_index])^exponent
  hit_increase <- hit_weights / sum(hit_weights)
  miss_decrease <- ifelse(n_genes > n_hits, 1 / (n_genes - n_hits), 0)

  running_score <- numeric(n_genes)
  running_score[hit_index] <- hit_increase
  running_score[!hit_index] <- -miss_decrease
  running_score <- cumsum(running_score)

  list(
    es_df = data.frame(rank = seq_len(n_genes), running_score = running_score),
    hit_df = data.frame(rank = which(hit_index)),
    metric_df = data.frame(rank = seq_len(n_genes), metric = as.numeric(gene_list)),
    strip_df = data.frame(rank = seq_len(n_genes), y = 1, metric = as.numeric(gene_list))
  )
}

# 绘制单个 GO 条目的 GSEA 曲线、barcode 和 ranked metric 组合图。
plot_gsea_term <- function(gsea_res, go_id, go_name) {
  result_df <- as.data.frame(gsea_res)
  p_adjust <- result_df$p.adjust[match(go_id, result_df$ID)]
  plot_data <- build_gsea_plot_data(gsea_res, go_id)

  x_breaks <- pretty(plot_data$es_df$rank, n = 4)
  x_breaks <- x_breaks[x_breaks > 0 & x_breaks <= max(plot_data$es_df$rank)]

  x_scale <- scale_x_continuous(
    breaks = x_breaks,
    expand = expansion(mult = c(0.01, 0.01))
  )

  # GSEA 专用高级配色：深湖蓝曲线 + 鲜明红蓝 ranked metric 渐变。
  gsea_line_col <- "#006D77"
  gsea_high_col <- "#D94F45"
  gsea_mid_col  <- "#F7F3EF"
  gsea_low_col  <- "#3C5488"
  gsea_border_col <- "#111111"

  p_es <- ggplot(plot_data$es_df, aes(x = rank, y = running_score)) +
    geom_hline(yintercept = 0, color = "grey72", linewidth = 0.42) +
    geom_line(color = gsea_line_col, linewidth = 1.45, lineend = "round") +
    annotate(
      "label",
      x = Inf,
      y = Inf,
      label = paste0("p.adjust\n", format_p_adjust(p_adjust)),
      hjust = 1.02,
      vjust = 1.15,
      size = 4.2,
      fontface = "bold",
      label.size = 0.32,
      fill = "#FAF7F2",
      color = "black"
    ) +
    x_scale +
    labs(title = go_name, x = NULL, y = "Running Enrichment Score") +
    theme_nature(base_size = 12) +
    theme(
      panel.grid.major.x = element_line(color = "grey86", linewidth = 0.35),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.25),
      panel.border = element_rect(color = gsea_border_col, fill = NA, linewidth = 0.5),
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 13, face = "bold", color = "black"),
      axis.title.y = element_text(size = 15, face = "bold", margin = margin(r = 8)),
      plot.title = element_text(size = 28, face = "bold", hjust = 0.5, margin = margin(b = 22)),
      plot.margin = margin(8, 10, 2, 10)
    )

  p_hits <- ggplot() +
    geom_tile(
      data = plot_data$strip_df,
      aes(x = rank, y = 0, fill = metric),
      height = 0.30,
      width = 1
    ) +
    geom_segment(
      data = plot_data$hit_df,
      aes(x = rank, xend = rank, y = 0.18, yend = 0.98),
      linewidth = 0.45,
      color = "black"
    ) +
    scale_fill_gradient2(
      low = gsea_low_col,
      mid = gsea_mid_col,
      high = gsea_high_col,
      midpoint = 0,
      guide = "none"
    ) +
    x_scale +
    coord_cartesian(ylim = c(-0.18, 1.02), clip = "off") +
    theme_void(base_size = 12) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(color = gsea_border_col, fill = NA, linewidth = 0.45),
      plot.margin = margin(2, 10, 8, 10)
    )

  p_metric <- ggplot(plot_data$metric_df, aes(x = rank, y = metric, fill = metric)) +
    geom_col(width = 1, color = NA, alpha = 0.92) +
    scale_fill_gradient2(
      low = gsea_low_col,
      mid = "grey78",
      high = gsea_high_col,
      midpoint = 0,
      guide = "none"
    ) +
    x_scale +
    labs(x = "Rank in Ordered Dataset", y = "Ranked List Metric") +
    theme_nature(base_size = 12) +
    theme(
      panel.grid.major.x = element_line(color = "grey86", linewidth = 0.35),
      panel.grid.major.y = element_blank(),
      panel.border = element_rect(color = gsea_border_col, fill = NA, linewidth = 0.5),
      axis.line = element_blank(),
      axis.text = element_text(size = 13, face = "bold", color = "black"),
      axis.title.x = element_text(size = 15, face = "bold", margin = margin(t = 8)),
      axis.title.y = element_text(size = 15, face = "bold", margin = margin(r = 8)),
      plot.margin = margin(8, 10, 8, 10)
    )

  g_es <- ggplotGrob(p_es)
  g_hits <- ggplotGrob(p_hits)
  g_metric <- ggplotGrob(p_metric)
  max_widths <- grid::unit.pmax(g_es$widths, g_hits$widths, g_metric$widths)
  g_es$widths <- max_widths
  g_hits$widths <- max_widths
  g_metric$widths <- max_widths

  gridExtra::arrangeGrob(
    g_es,
    g_hits,
    g_metric,
    ncol = 1,
    heights = c(3.25, 0.72, 1.75)
  )
}

draw_gsea_plot <- function(p) {
  if (inherits(p, c("grob", "gtable", "gTree"))) {
    grid::grid.draw(p)
  } else {
    print(p)
  }
}

# 把单个 GSEA 图同时保存为 PDF 和 PNG。
save_gseaplot <- function(p, filename, width = 7.2, height = 7, dpi = 300) {
  pdf_file <- paste0(filename, ".pdf")
  png_file <- paste0(filename, ".png")
  
  pdf(pdf_file, width = width, height = height, useDingbats = FALSE)
  draw_gsea_plot(p)
  dev.off()
  
  png(png_file, width = width, height = height, units = "in", res = dpi, bg = "white")
  draw_gsea_plot(p)
  dev.off()
}

# 筛选显著正向富集的 Top GSEA 条目，逐个出图并生成合并图。
plot_top_gsea_terms <- function(gsea_res, gsea_res_df, comparison_name, out_dir,
                                top_n = 10) {
  top_gsea <- gsea_res_df %>%
    filter(NES > 0) %>%
    arrange(p.adjust) %>%
    slice_head(n = top_n)
  
  write.csv(
    top_gsea,
    file.path(out_dir, paste0("top_gsea_terms_", comparison_name, ".csv")),
    row.names = FALSE
  )
  
  if (nrow(top_gsea) == 0) {
    return(NULL)
  }
  
  gsea_plots <- lapply(seq_len(nrow(top_gsea)), function(i) {
    go_id <- top_gsea$ID[i]
    go_name <- as.character(top_gsea$Description[i])
    
    p <- plot_gsea_term(gsea_res, go_id, go_name)
    
    safe_name <- gsub("[^A-Za-z0-9]+", "_", go_name)
    safe_name <- gsub("^_|_$", "", safe_name)
    
    save_gseaplot(
      p = p,
      filename = file.path(out_dir, paste0("GSEA_", safe_name)),
      width = 7.2,
      height = 7,
      dpi = 300
    )
    
    p
  })
  
  names(gsea_plots) <- as.character(top_gsea$Description)
  
  gsea_grobs <- lapply(gsea_plots, function(p) {
    grid::grid.grabExpr(draw_gsea_plot(p))
  })
  
  title_grob <- grid::textGrob(
    paste0("Top GO-GSEA pathways (", gsub("_", " ", comparison_name), ")"),
    gp = grid::gpar(fontsize = 24, fontface = "bold")
  )
  
  p_gsea_combined <- gridExtra::arrangeGrob(
    grobs = gsea_grobs,
    ncol = 2,
    top = title_grob
  )
  
  pdf(
    file.path(out_dir, "GSEA_top_terms_combined.pdf"),
    width = 15,
    height = 32,
    useDingbats = FALSE
  )
  grid::grid.draw(p_gsea_combined)
  dev.off()
  
  png(
    file.path(out_dir, "GSEA_top_terms_combined.png"),
    width = 15,
    height = 32,
    units = "in",
    res = 300
  )
  grid::grid.draw(p_gsea_combined)
  dev.off()
  
  gsea_plots
}

# 完成一个巨噬细胞亚型 vs 其他亚型的 DEG、GSEA 和配套图表输出。
run_one_macrophage_comparison <- function(seurat_obj, ident_name, out_base_dir) {
  comparison_name <- paste0(ident_name, "_vs_others")
  out_dir <- file.path(out_base_dir, comparison_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  degs <- run_deg(seurat_obj, ident_name)
  write.csv(
    degs,
    file.path(out_dir, paste0("degs_", comparison_name, ".csv")),
    row.names = FALSE
  )
  
  degs_annotated <- annotate_deg(degs)
  write.csv(
    degs_annotated,
    file.path(out_dir, paste0("degs_", comparison_name, "_annotated.csv")),
    row.names = FALSE
  )
  
  gsea <- run_gsea(degs)
  write.csv(
    gsea$gsea_res_df,
    file.path(out_dir, paste0("gsea_", comparison_name, ".csv")),
    row.names = FALSE
  )
  
  volcano_extra_label_genes <- switch(
    ident_name,
    IAM = "SPP1",
    character(0)
  )
  plot_volcano(degs, comparison_name, out_dir,
               extra_label_genes = volcano_extra_label_genes)
  plot_deg_heatmap(seurat_obj, degs, ident_name, out_dir)
  plot_top_gsea_terms(gsea$gsea_res, gsea$gsea_res_df, comparison_name, out_dir)
  
  list(
    degs = degs,
    degs_annotated = degs_annotated,
    gsea = gsea
  )
}

# ==============================================================================
# 6. 运行四组比较：每一种巨噬细胞 vs 其他三种巨噬细胞
# ==============================================================================

macrophage_subtype_abundance_by_sample <- plot_macrophage_subtype_abundance_by_sample(
  seurat_obj = macrophage_cells,
  macrophage_types = macrophage_types,
  out_base_dir = out_base_dir
)

macrophage_deg_go_results <- lapply(
  macrophage_types,
  function(x) run_one_macrophage_comparison(macrophage_cells, x, out_base_dir)
)

names(macrophage_deg_go_results) <- macrophage_types

# ==============================================================================
# 7. 可视化：精选与巨噬细胞亚型功能相关的 GSEA 通路
# ==============================================================================

selected_gsea_out_dir <- file.path(out_base_dir, "selected_functional_gsea")
dir.create(selected_gsea_out_dir, recursive = TRUE, showWarnings = FALSE)

selected_gsea_patterns <- list(
  LAM = list(
    include = "lysosome|vacuolar|lytic vacuole|phagocyt|endocyt|fatty acid|lipid|sphingolipid|cholesterol|triglyceride|PPAR|foam cell|lipoprotein|steroid|atheroscler|neutral lipid|fatty liver|oxidative phosphorylation|electron transport|respiratory chain|reactive oxygen|ROS|MTORC1",
    preferred = c(
      "lysosome",
      "primary lysosome",
      "fatty acid metabolic process",
      "lipid metabolic process",
      "cholesterol homeostasis",
      "oxidative phosphorylation",
      "reactive oxygen species pathway"
    )
  ),
  IAM = list(
    include = "cytokine|chemokine|chemotaxis|taxis|leukocyte activation|leukocyte migration|inflammatory|immune|macrophage activation|NF-kB|NFKB|TNF|IL-1|IL-6|CXCL|CCL|monocyte|neutrophil|interferon|TLR|toll-like|innate immun|JAK|STAT",
    preferred = c(
      "tnfa signaling via nfkb",
      "inflammatory response",
      "cytokine activity",
      "cytokine receptor binding",
      "il6 jak stat3 signaling",
      "interferon gamma response",
      "leukocyte activation",
      "chemokine receptor binding",
      "chemotaxis"
    )
  ),
  OSM = list(
    include = "oncostatin|JAK|STAT|acute phase|oxidoreductase|oxidative|redox|glutathione|reactive oxygen|ROS|glycolytic|glycolysis|carbohydrate|carbon metabolism|HIF|hypoxia|pyridine|nicotinamide|OXPHOS|oxidative phosphorylation|electron transport|mitochondr",
    preferred = c(
      "il6 jak stat3 signaling",
      "hypoxia",
      "reactive oxygen species pathway",
      "oxidoreductase activity",
      "glycolytic process",
      "carbohydrate metabolic process",
      "oxidative phosphorylation"
    )
  ),
  TRM = list(
    include = "phagocytosis|phagosome|endocytosis|efferocytosis|Fc gamma|complement|lysosome|homeosta|ECM|extracellular matrix|tissue remodel|collagen|scavenger|MHC|antigen present|enzyme-linked receptor|TGF|transforming growth factor|receptor signaling|clathrin",
    preferred = c(
      "complement and coagulation cascades",
      "allograft rejection",
      "ecm receptor interaction",
      "endocytosis",
      "phagosome",
      "enzyme-linked receptor protein signaling pathway",
      "transforming growth factor beta receptor signaling pathway",
      "clathrin-dependent endocytosis"
    )
  )
)

selected_gsea_exclude_pattern <- paste(
  c(
    "ribosomal subunit",
    "structural constituent of ribosome",
    "cytosolic ribosome",
    "plasma membrane$",
    "cell periphery$",
    "developmental process$",
    "positive regulation of macromolecule metabolic process$",
    "positive regulation of gene expression$"
  ),
  collapse = "|"
)

# 清理 GSEA 通路名称中的特殊字符，生成适合作为文件名的字符串。
safe_gsea_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  gsub("^_|_$", "", x)
}

# 按亚型功能关键词筛选更有生物学解释意义的正向 GSEA 通路。
select_functional_gsea_terms <- function(subtype, top_n = 6) {
  comparison_name <- paste0(subtype, "_vs_others")
  gsea_file <- file.path(
    out_base_dir,
    comparison_name,
    paste0("gsea_", comparison_name, ".csv")
  )

  if (!file.exists(gsea_file)) {
    warning("GSEA result file not found: ", gsea_file)
    return(data.frame())
  }

  pattern_info <- selected_gsea_patterns[[subtype]]

  gsea_df <- read.csv(gsea_file, stringsAsFactors = FALSE, check.names = FALSE) %>%
    mutate(
      macro_subtype = subtype,
      p.adjust = as.numeric(p.adjust),
      NES = as.numeric(NES),
      Description_match = gsub("_", " ", Description),
      theme_related = grepl(pattern_info$include, Description_match, ignore.case = TRUE),
      preferred_term = grepl(
        paste(pattern_info$preferred, collapse = "|"),
        Description_match,
        ignore.case = TRUE
      ),
      broad_or_low_specificity = grepl(
        selected_gsea_exclude_pattern,
        Description_match,
        ignore.case = TRUE
      )
    ) %>%
    filter(
      !is.na(p.adjust),
      !is.na(NES),
      p.adjust < 0.05,
      NES > 0,
      theme_related,
      !broad_or_low_specificity
    ) %>%
    mutate(
      strong_signal = NES > 1.5,
      moderate_signal = NES > 1.2
    ) %>%
    arrange(desc(preferred_term), desc(strong_signal), desc(moderate_signal), desc(NES), p.adjust)

  if (nrow(gsea_df) == 0) {
    return(gsea_df)
  }

  strong_terms <- gsea_df %>%
    filter(strong_signal) %>%
    slice_head(n = top_n)

  if (nrow(strong_terms) >= top_n) {
    return(strong_terms)
  }

  bind_rows(
    strong_terms,
    gsea_df %>%
      filter(!ID %in% strong_terms$ID) %>%
      slice_head(n = top_n - nrow(strong_terms))
  ) %>%
    slice_head(n = top_n)
}

# 获取指定亚型的 GSEA 结果对象；若内存中没有，则尝试用 DEG 文件重建。
get_gsea_result_for_plotting <- function(subtype) {
  if (
    exists("macrophage_deg_go_results") &&
      !is.null(macrophage_deg_go_results[[subtype]]) &&
      !is.null(macrophage_deg_go_results[[subtype]]$gsea$gsea_res)
  ) {
    return(macrophage_deg_go_results[[subtype]]$gsea$gsea_res)
  }

  comparison_name <- paste0(subtype, "_vs_others")
  deg_file <- file.path(
    out_base_dir,
    comparison_name,
    paste0("degs_", comparison_name, ".csv")
  )

  if (!file.exists(deg_file)) {
    warning("DEG file not found for rebuilding GSEA object: ", deg_file)
    return(NULL)
  }

  run_gsea(read.csv(deg_file, stringsAsFactors = FALSE, check.names = FALSE))$gsea_res
}

selected_gsea_terms <- lapply(macrophage_types, select_functional_gsea_terms)
names(selected_gsea_terms) <- macrophage_types
selected_gsea_terms_df <- bind_rows(selected_gsea_terms)

write.csv(
  selected_gsea_terms_df,
  file.path(selected_gsea_out_dir, "selected_functional_gsea_terms.csv"),
  row.names = FALSE
)

selected_gsea_plots <- lapply(macrophage_types, function(subtype) {
  subtype_terms <- selected_gsea_terms[[subtype]]
  subtype_out_dir <- file.path(selected_gsea_out_dir, subtype)
  dir.create(subtype_out_dir, recursive = TRUE, showWarnings = FALSE)

  if (nrow(subtype_terms) == 0) {
    message("No selected functional GSEA terms for ", subtype)
    return(NULL)
  }

  gsea_res <- get_gsea_result_for_plotting(subtype)
  if (is.null(gsea_res)) {
    return(NULL)
  }

  subtype_plots <- lapply(seq_len(nrow(subtype_terms)), function(i) {
    go_id <- subtype_terms$ID[i]
    go_name <- subtype_terms$Description[i]

    p <- plot_gsea_term(gsea_res, go_id, go_name)
    safe_name <- safe_gsea_filename(paste(subtype, go_name, sep = "_"))

    save_gseaplot(
      p = p,
      filename = file.path(subtype_out_dir, paste0("Selected_GSEA_", safe_name)),
      width = 7.2,
      height = 7,
      dpi = 300
    )

    p
  })

  names(subtype_plots) <- subtype_terms$Description
  subtype_plots
})

names(selected_gsea_plots) <- macrophage_types
