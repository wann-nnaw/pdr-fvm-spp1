
library(dplyr)
library(readr)
library(ggplot2)
library(ggVennDiagram)
library(tidyr)

# ==============================================================================
# 1. 设置输入和输出目录
# ==============================================================================
input_dir <- "results/macrophage_DEG_GSEA"
out_dir <- "results/venn"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 2. 设置亚型和输入文件
# ==============================================================================
subtypes <- c("IAM", "LAM", "OSM", "TRM")
heatmap_subtypes <- c("OSM", "LAM", "IAM", "TRM")
input_files <- file.path(
  input_dir,
  paste0(subtypes, "_vs_others"),
  paste0("degs_", subtypes, "_vs_others_annotated.csv")
)
names(input_files) <- subtypes

# ==============================================================================
# 3. 读取并清洗 DEG 表
# ==============================================================================
deg_tables <- lapply(input_files, function(file) {
  dat <- read_csv(file, show_col_types = FALSE)

  dat %>%
    dplyr::mutate(gene = trimws(as.character(gene))) %>%
    dplyr::filter(change %in% c("Up", "Down")) %>%
    dplyr::filter(!is.na(gene), gene != "") %>%
    dplyr::distinct(gene, .keep_all = TRUE)
})

# ==============================================================================
# 4. 提取 4 个亚型共同 DEG
# ==============================================================================
gene_sets <- lapply(deg_tables, function(dat) dat$gene)
names(gene_sets) <- subtypes

intersect_genes <- Reduce(intersect, gene_sets) %>% sort()

prefixed_tables <- Map(function(dat, subtype) {
  dat %>%
    dplyr::filter(gene %in% intersect_genes) %>%
    dplyr::arrange(match(gene, intersect_genes)) %>%
    dplyr::rename_with(~ paste0(subtype, "_", .x), .cols = -gene)
}, deg_tables, subtypes)

intersect_table <- Reduce(function(x, y) dplyr::full_join(x, y, by = "gene"), prefixed_tables) %>%
  dplyr::select(gene, dplyr::everything()) %>%
  dplyr::arrange(gene)

intersect_csv <- file.path(out_dir, "macrophage_4subtype_intersect_genes_annotated.csv")
write_csv(intersect_table, intersect_csv)

# ==============================================================================
# 5. 汇总共同 DEG 的上下调方向
# ==============================================================================
change_cols <- paste0(subtypes, "_change")
logfc_cols <- paste0(subtypes, "_avg_log2FC")

direction_summary <- intersect_table %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    direction_class = dplyr::case_when(
      all(dplyr::c_across(dplyr::all_of(change_cols)) == "Up") ~ "All_Up",
      all(dplyr::c_across(dplyr::all_of(change_cols)) == "Down") ~ "All_Down",
      TRUE ~ "Mixed"
    ),
    up_in = paste(subtypes[dplyr::c_across(dplyr::all_of(change_cols)) == "Up"], collapse = ";"),
    down_in = paste(subtypes[dplyr::c_across(dplyr::all_of(change_cols)) == "Down"], collapse = ";"),
    direction_pattern = paste0(
      "Up: ", ifelse(up_in == "", "None", up_in),
      " | Down: ", ifelse(down_in == "", "None", down_in)
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(gene, direction_class, up_in, down_in, direction_pattern, dplyr::everything())

direction_csv <- file.path(out_dir, "macrophage_4subtype_intersect_genes_direction_classified.csv")
write_csv(direction_summary, direction_csv)

# ==============================================================================
# 6. 整理方向热图数据
# ==============================================================================
direction_long <- direction_summary %>%
  dplyr::select(gene, direction_class, direction_pattern, dplyr::all_of(change_cols), dplyr::all_of(logfc_cols)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(change_cols),
    names_to = "subtype",
    values_to = "change"
  ) %>%
  dplyr::mutate(subtype = sub("_change$", "", subtype)) %>%
  dplyr::left_join(
    direction_summary %>%
      dplyr::select(gene, dplyr::all_of(logfc_cols)) %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(logfc_cols),
        names_to = "subtype",
        values_to = "avg_log2FC"
      ) %>%
      dplyr::mutate(subtype = sub("_avg_log2FC$", "", subtype)),
    by = c("gene", "subtype")
  ) %>%
  dplyr::mutate(
    subtype = factor(subtype, levels = heatmap_subtypes),
    direction_symbol = ifelse(change == "Up", "Up", "Down")
  )

gene_order <- direction_summary %>%
  dplyr::arrange(dplyr::desc(OSM_avg_log2FC), gene) %>%
  dplyr::pull(gene)

direction_long <- direction_long %>%
  dplyr::mutate(gene = factor(gene, levels = rev(gene_order)))

heatmap_height <- max(7, length(gene_order) * 0.16 + 2)

# ==============================================================================
# 7. 绘制共同 DEG 方向热图
# ==============================================================================
direction_heatmap <- ggplot(
  direction_long,
  aes(x = subtype, y = gene, fill = avg_log2FC)
) +
  geom_tile(
    color = "#FFFFFF",
    linewidth = 0.35
  ) +
  scale_fill_gradient2(
    low = "#2F6FBB",
    mid = "#FFFFFF",
    high = "#C84C61",
    midpoint = 0,
    name = "avglog2FC",
    guide = guide_colorbar(
      barheight = unit(4.2, "cm"),
      barwidth = unit(0.45, "cm"),
      title.position = "top"
    )
  ) +
  labs(
    title = "Direction of Shared Significant DEGs",
    x = NULL,
    y = NULL
  ) +
  theme_nature(base_size = 11) +
  theme(
    plot.title = element_text(
      size = 26,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 14),
      color = "black"
    ),
    
    axis.text.x = element_text(
      size = 14,
      face = "bold",
      color = "black",
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    
    legend.position = "right",
    legend.justification = "top",
    legend.box.just = "top",
    legend.title = element_text(
      size = 16,
      face = "bold",
      color = "black",
      lineheight = 0.9,
      margin = margin(b = 10)
    ),
    legend.text = element_text(
      size = 14,
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

heatmap_png <- file.path(out_dir, "macrophage_4subtype_intersect_direction_heatmap.png")
heatmap_pdf <- file.path(out_dir, "macrophage_4subtype_intersect_direction_heatmap.pdf")

ggsave(
  heatmap_png,
  direction_heatmap,
  width = 7,
  height = 9,
  dpi = 300,
  bg = "white"
)

ggsave(
  heatmap_pdf,
  direction_heatmap,
  width = 7,
  height = 9,
  bg = "white"
)

# ==============================================================================
# 8. 绘制 4 个亚型 DEG 交集 Venn 图
# ==============================================================================
venn_plot <- ggVennDiagram(
  gene_sets,
  label = "count",
  label_alpha = 0,
  category.names = subtypes,
  set_color = c("#144B6E", "#0E8F87", "#B23A6F", "#D88A22"),
  set_size = 7,
  label_size = 6
) +
  scale_fill_gradientn(
    colours = c("#F7F4EE", "#BFE3DA", "#74B6C8", "#D481A3", "#D9A441"),
    values = scales::rescale(c(0, 0.15, 0.4, 0.7, 1))
  ) +
  labs(
    title = "Venn plot"
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.title = element_text(
      size = 28,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 16),
      color = "black"
    ),
    
    legend.position = "none",
    
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(10, 16, 10, 10)
  )

png_file <- file.path(out_dir, "macrophage_4subtype_Venn.png")
pdf_file <- file.path(out_dir, "macrophage_4subtype_Venn.pdf")

ggsave(
  png_file,
  venn_plot,
  width = 8,
  height = 7,
  dpi = 300,
  bg = "white"
)

ggsave(
  pdf_file,
  venn_plot,
  width = 8,
  height = 7,
  bg = "white"
)

