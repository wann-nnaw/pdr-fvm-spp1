load("rda/UMAP_annotated.rda")

library(CellChat)
library(Seurat)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

# 统一各子图的字体与字号。A–C 仅放大主标题，其余样式保持不变。
# Helvetica 是 PDF 设备原生支持的无衬线字体，在 macOS/Linux 上导出更稳定。
figure_font_family <- "Helvetica"
title_size_ac <- 28
title_size_lr <- 28

dir.create("results/cellchat", showWarnings = FALSE, recursive = TRUE)

# ==============================================
# 1. 读取 Seurat 对象
# ==============================================
scRNA <- combined_anno
DefaultAssay(scRNA) <- "RNA"

# 标准化
scRNA <- NormalizeData(scRNA)

# 使用标准化表达矩阵
data.input <- GetAssayData(scRNA, assay = "RNA", layer = "data")

# 提取元数据
meta = scRNA@meta.data 

# 确保细胞类型列名称正确
meta$ident <- meta$cell_type

# 保证顺序一致
cell.use <- colnames(data.input)
meta <- meta[cell.use, ]
all(rownames(meta) == colnames(data.input))

# ==============================================
# 3. 定义一个函数：自动跑全套 CellChat
# ==============================================
run_cellchat <- function(data_input, cell_meta, out_prefix){
  
  # 创建对象
  cellchat <- createCellChat(object = data_input,meta = cell_meta,group.by = "ident")
  cellchat <- addMeta(cellchat, meta = cell_meta)
  cellchat <- setIdent(cellchat, ident.use = "ident")
  cellchat@idents <- droplevels(cellchat@idents)
  
  # 数据库
  CellChatDB <- CellChatDB.human
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
  cellchat@DB <- CellChatDB.use
  
  ## 预处理
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  
  # 推断细胞通讯
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  # 保存结果
  saveRDS(cellchat, file = paste0(out_prefix, "_cellchat.rds"))
  write.csv(subsetCommunication(cellchat), paste0(out_prefix, "_net.csv"),  row.names = FALSE)
  
  # 可视化
  
  # 热图 
  mat <- cellchat@net$count
  col_fun <- colorRamp2(
    c(min(mat),
      quantile(mat, 0.2),
      quantile(mat, 0.4),
      quantile(mat, 0.6),
      quantile(mat, 0.8),
      max(mat)),
    c("#EAF4FF",  # 极浅蓝
      "#A7D8F0",  # 浅蓝
      "#5FA8D3",  # 中蓝
      "#F7E6A1",  # 柔和金黄
      "#F29A76",  # 珊瑚橙
      "#D95F59")  # 珊瑚红
  )
  
  pdf(paste0(out_prefix, "_heatmap.pdf"), width=8, height=7)
  draw(Heatmap(mat,
               name = "Interaction\nstrength",
               col = col_fun,
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               column_title = "Heatmap",
               column_title_gp = grid::gpar(
                 fontsize = title_size_ac,
                 fontface = "bold",
                 fontfamily = figure_font_family
               ),
               row_names_side = "left",
               column_names_rot = 45,
               row_names_gp = grid::gpar(
                 fontsize = 16, fontface = "bold", col = "black"
               ),
               column_names_gp = grid::gpar(
                 fontsize = 15, fontface = "bold", col = "black"
               ),
               rect_gp = grid::gpar(col = "white", lwd = 1),
               heatmap_legend_param = list(
                 title = "Strength",
                 title_gp = grid::gpar(fontsize = 16, fontface = "bold"),
                 labels_gp = grid::gpar(fontsize = 14, fontface = "bold"),
                 legend_direction = "vertical"
               )))
  dev.off()
  
  # 网络图
  groupSize <- as.numeric(table(cellchat@idents))
  
  pdf(paste0(out_prefix, "_circle_number.pdf"), width = 8, height = 8)
  par(font = 2)
  netVisual_circle(
    cellchat@net$count,
    vertex.weight = groupSize,
    weight.scale = TRUE,
    label.edge = FALSE,
    title.name = NULL,
    vertex.label.cex = 1.8,
    vertex.label.color = "black"
  )
  text(
    0, 1.5, "Number of interactions",
    cex = 2.6, font = 2, family = figure_font_family, xpd = NA
  )
  dev.off()
  
  pdf(paste0(out_prefix, "_circle_weight.pdf"),width = 8, height = 8)
  par(font = 2)
  netVisual_circle(
    cellchat@net$weight,
    vertex.weight = groupSize,
    weight.scale = TRUE,
    label.edge = FALSE,
    title.name = NULL,
    vertex.label.cex = 1.8,
    vertex.label.color = "black"
  )
  text(
    0, 1.5, "Strength of interactions",
    cex = 2.6, font = 2, family = figure_font_family, xpd = NA
  )
  dev.off()
  
  return(cellchat)
}

# ==============================================
# 4. 直接运行
# ==============================================
cellchat <- run_cellchat(data.input, meta, "results/cellchat/DR")

# ==============================================
# 5. 画分组气泡图
# ==============================================
plot_bubble <- function(cellchat_obj, title){
  p <- netVisual_bubble(
    cellchat_obj,
    sources.use = c("LAM", "IAM", "OSM", "TRM"),
    targets.use = c("Monocyte", "Microglia", "Fibroblast", "T cells", "Pericyte", "Endothelial", "DC"),
    remove.isolate = FALSE
  ) +
    ggtitle(title) +
    theme(
      # 横坐标文字设置
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        face = "bold",
        size = 13,
        family = figure_font_family,
        color = "black"
      ),
      
      axis.text.y = element_text(
        size = 15,
        family = figure_font_family,
        face = "bold",
        color = "black"
      ),
      # 删除 CellChat 自动生成的 source.target 和 interaction_name_2 轴标题，
      # 保留横轴的细胞组合与纵轴的配体–受体名称。
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(
        hjust = 0.5, size = title_size_lr, face = "bold",
        family = figure_font_family, margin = margin(b = 10)
      ),
      legend.title = element_text(size = 17, face = "bold", family = figure_font_family),
      legend.text = element_text(size = 15, face = "bold", family = figure_font_family),
      
      # 扩大底部边距，防止标签被切掉
      plot.margin = margin(t = 10, r = 10, b = 40, l = 10),
      
      # 去掉横坐标标签截断
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    
    # 强制不裁剪标签（最关键！）
    coord_cartesian(clip = "off")

  # CellChat 会从对象中携带默认标题（如 "DR"），在返回前强制覆盖。
  p$labels$title <- title
  return(p)
}

pdf("results/cellchat/DR_bubble.pdf", width = 11.5, height = 14)
p1 <- plot_bubble(cellchat, "Ligand–receptor communication")
# 打印前再次覆盖，避免 patchwork/CellChat 恢复对象默认标题。
p1$labels$title <- "Ligand–receptor communication"
print(patchwork::wrap_plots(p1))
dev.off()

# ==============================================
# 6. 配体-受体通讯桑基图
#    与上方气泡图使用相同的发送和接收细胞，并以通信概率表示流带宽度
# ==============================================
plot_lr_sankey <- function(cellchat_obj, title,
                           sources.use,
                           targets.use,
                           out_prefix,
                           top_n_lr = NULL,
                           force_ligands = character(0)) {
  if (!requireNamespace("ggalluvial", quietly = TRUE)) {
    stop("Package 'ggalluvial' is required. Install it with install.packages('ggalluvial').")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required. Install it with install.packages('dplyr').")
  }

  lr_flow <- subsetCommunication(
    cellchat_obj,
    sources.use = sources.use,
    targets.use = targets.use
  )

  required_cols <- c("source", "ligand", "receptor", "target", "prob")
  missing_cols <- setdiff(required_cols, colnames(lr_flow))
  if (length(missing_cols) > 0) {
    stop("Missing columns in CellChat communication table: ",
         paste(missing_cols, collapse = ", "))
  }

  lr_flow <- lr_flow |>
    dplyr::filter(
      source %in% sources.use,
      target %in% targets.use,
      is.finite(prob),
      prob > 0
    ) |>
    dplyr::group_by(source, ligand, receptor, target) |>
    dplyr::summarise(
      communication_prob = sum(prob, na.rm = TRUE),
      n_interactions = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      source = factor(source, levels = sources.use),
      target = factor(target, levels = targets.use)
    )

  if (nrow(lr_flow) == 0) {
    stop("No significant ligand-receptor communications were found for the selected cell groups.")
  }

  # 可选：每个 Sender 仅展示通信概率最高的配体-受体组合；
  # force_ligands 中指定的配体始终保留，不受 Top N 筛选影响。
  if (!is.null(top_n_lr)) {
    keep_lr <- lr_flow |>
      dplyr::group_by(source, ligand, receptor) |>
      dplyr::summarise(
        pair_prob = sum(communication_prob, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::group_by(source) |>
      dplyr::slice_max(order_by = pair_prob, n = top_n_lr, with_ties = FALSE) |>
      dplyr::ungroup() |>
      dplyr::select(source, ligand, receptor)

    forced_lr <- lr_flow |>
      dplyr::filter(as.character(ligand) %in% force_ligands) |>
      dplyr::distinct(source, ligand, receptor)

    keep_lr <- dplyr::bind_rows(keep_lr, forced_lr) |>
      dplyr::distinct(source, ligand, receptor)

    lr_flow <- lr_flow |>
      dplyr::semi_join(keep_lr, by = c("source", "ligand", "receptor"))
  }

  # 保存桑基图源数据，保证图中每条流均可追溯
  write.csv(
    lr_flow,
    paste0(out_prefix, "_sankey_source_data.csv"),
    row.names = FALSE
  )

  sender_colors <- c(
    "LAM" = "#4E79A7",
    "IAM" = "#59A14F",
    "OSM" = "#F28E2B",
    "TRM" = "#E15759"
  )
  flow_colors <- sender_colors[sources.use]
  missing_color <- is.na(flow_colors)
  if (any(missing_color)) {
    flow_colors[missing_color] <- grDevices::hcl.colors(sum(missing_color), "Dark 3")
  }

  p <- ggplot(
    lr_flow,
    aes(
      axis1 = source,
      axis2 = ligand,
      axis3 = receptor,
      axis4 = target,
      y = communication_prob
    )
  ) +
    ggalluvial::geom_alluvium(
      aes(fill = source),
      width = 0.12,
      alpha = 0.58,
      knot.pos = 0.35,
      discern = TRUE,
      color = NA
    ) +
    ggalluvial::geom_stratum(
      width = 0.16,
      fill = "white",
      color = "#4A4A4A",
      linewidth = 0.3,
      discern = TRUE
    ) +
    ggplot2::geom_text(
      stat = ggalluvial::StatStratum,
      aes(label = after_stat(stratum)),
      size = 5.2,
      family = figure_font_family,
      fontface = "bold",
      color = "black",
      check_overlap = TRUE,
      discern = TRUE
    ) +
    scale_x_discrete(
      limits = c("Sender", "Ligand", "Receptor", "Receiver"),
      expand = c(0.06, 0.06)
    ) +
    scale_fill_manual(values = flow_colors, drop = FALSE) +
    labs(
      title = title,
      subtitle = NULL,
      x = NULL,
      y = "Summed communication probability",
      fill = "Sender"
    ) +
    theme_classic(base_size = 15, base_family = figure_font_family) +
    theme(
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_text(face = "bold", size = 15, color = "black"),
      axis.text.y = element_text(face = "bold", size = 13, color = "black"),
      axis.title.y = element_text(face = "bold", size = 15),
      plot.title = element_text(
        hjust = 0.5, face = "bold", size = 30,
        margin = margin(b = 3)
      ),
      legend.title = element_text(size = 13, face = "bold"),
      legend.text = element_text(size = 12, face = "bold"),
      # 每幅图只有一个 Sender，底部图例信息重复，直接取消。
      legend.position = "none",
      plot.margin = margin(10, 20, 5, 20)
    ) +
    coord_cartesian(clip = "off")

  ggsave(
    paste0(out_prefix, "_sankey.pdf"),
    plot = p,
    # 中等画布尺寸：兼顾单图内容容量与插入组图后的可读性。
    width = 12,
    height = 8,
    units = "in",
    device = grDevices::pdf
  )

  if (requireNamespace("svglite", quietly = TRUE)) {
    ggsave(
      paste0(out_prefix, "_sankey.svg"),
      plot = p,
      width = 12,
      height = 8,
      units = "in",
      device = svglite::svglite
    )
  }

  return(p)
}

sources.use <- c("LAM", "IAM", "OSM", "TRM")
targets.use <- c("Monocyte", "Microglia", "Fibroblast", "T cells",
                 "Pericyte", "Endothelial", "DC")

p_sankey <- plot_lr_sankey(
  cellchat_obj = cellchat,
  title = "DR ligand-receptor communication flows",
  sources.use = sources.use,
  targets.use = targets.use,
  out_prefix = "results/cellchat/DR"
)

# ==============================================
# 7. 清晰版桑基图：按 Sender 拆分，Top 15，并保留 SPP1
# ==============================================
for (sender in sources.use) {
  plot_lr_sankey(
    cellchat_obj = cellchat,
    title = paste0("DR ligand-receptor flows: ", sender),
    sources.use = sender,
    targets.use = targets.use,
    out_prefix = paste0("results/cellchat/DR_top15_", sender),
    top_n_lr = 15,
    force_ligands = "SPP1"
  )
}
